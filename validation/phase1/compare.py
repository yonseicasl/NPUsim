#!/usr/bin/env python3
"""Phase-1 cross-simulator validation: NPUsim (systolic_array, 12x14 WS-mapped)
vs SCALE-Sim v2 (12x14 WS, CALC bandwidth) on identical AlexNet layer shapes.

Comparison semantics
--------------------
* SCALE-Sim "Total Cycles" = systolic compute-schedule cycles of ONE image and ONE
  group (fill/drain included, memory stalls excluded in CALC mode).
* NPUsim "Computation cycle" = serialized MAC-issue cycles of the busiest PE under
  its mapping (batch and groups already included by the mapping/network).
* To place both on the same axis, SCALE-Sim cycles are scaled by batch x groups.
* Both simulators also report a MAC-utilization-style number; we derive
  cycles-per-MAC to factor out shape volume.
"""
import csv, math, re, sys, os

ROOT = os.path.dirname(os.path.abspath(__file__))
# Optional argv[1] selects the NPUsim mapping/result set: "energy" (original
# eyeriss-RS mapping) or "matched" (WS-matched conv1/conv2 mapping).
MAPPING = sys.argv[1] if len(sys.argv) > 1 else "energy"
NPU_RESULT_DIR = os.path.join(ROOT, f"../../result/scalesim/alexnet/{MAPPING}")
SS_CONV_REPORT = os.path.join(ROOT, "out/npusim_val_convs/COMPUTE_REPORT.csv")
SS_FC_REPORT = os.path.join(ROOT, "out/npusim_val_fcs/COMPUTE_REPORT.csv")

BATCH = 4
# (npusim layer id, report tag, scale-sim row id, name, group factor, MACs/image all groups)
LAYERS = [
    (0,  "conv", 0, "conv1",   1, 96*55*55*11*11*3),
    (2,  "conv", 1, "conv2",   2, 256*27*27*5*5*48),
    (4,  "conv", 2, "conv3",   1, 384*13*13*3*3*256),
    (5,  "conv", 3, "conv4",   2, 384*13*13*3*3*192),
    (6,  "conv", 4, "conv5",   2, 256*13*13*3*3*192),
    (8,  "fc",   0, "fc1",     1, 4096*6*6*256),
    (9,  "fc",   1, "fc2",     1, 4096*4096),
    (10, "fc",   2, "fc3",     1, 1000*4096),
]

def parse_npusim_layer(path):
    txt = open(path).read()
    def grab(pattern):
        m = re.search(pattern, txt)
        return float(m.group(1).replace(",", "")) if m else None
    # PE-array operand-stream (interconnection) cycles: these carry the systolic
    # wavefront fill per stream (SY2), so the array-schedule proxy comparable to
    # SCALE-Sim's compute schedule is max(comp, IC_in, IC_w). The OUTPUT stream is
    # excluded: NPUsim's mapping writes partial sums through per reduction step,
    # whereas the WS array accumulates in place (documented modeling difference).
    sec = txt.split("PE array result")[1].split("Global buffer result")[0]
    ic = [float(x) for x in re.findall(r"\* (?:Input data|Weight|Output data)\s*:\s*([\d.]+) cycles", sec)[:3]]
    return {
        "comp_cycle": grab(r"Computation cycle\s*:\s*([\d.]+)"),
        "ic_in": ic[0], "ic_w": ic[1], "ic_out": ic[2],
        "latency": grab(r"Critical-path latency\s*:\s*([\d.]+)"),
    }

def parse_scalesim(path):
    rows = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            rows[int(row["LayerID"])] = {
                "total_cycles": float(row[" Total Cycles"]),
                "util": float(row[" Overall Util %"]),
                "mapping_eff": float(row[" Mapping Efficiency %"]),
            }
    return rows

def spearman(a, b):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0]*len(v)
        for pos, i in enumerate(order):
            r[i] = pos
        return r
    ra, rb = rank(a), rank(b)
    n = len(a)
    d2 = sum((ra[i]-rb[i])**2 for i in range(n))
    return 1 - 6*d2/(n*(n*n-1))

def main():
    reports = {"conv": parse_scalesim(SS_CONV_REPORT), "fc": parse_scalesim(SS_FC_REPORT)}
    print(f"{'layer':7s} {'NPU comp':>12s} {'NPU IC_in':>12s} {'NPU IC_w':>12s} "
          f"{'array-sched':>12s} {'SS x b x g':>12s} {'err%':>7s} {'comp-only%':>10s}")
    npu_v, ss_v, errs = [], [], []
    for npu_id, tag, ss_id, name, groups, macs in LAYERS:
        n = parse_npusim_layer(os.path.join(NPU_RESULT_DIR, f"layer_{npu_id}.txt"))
        s = reports[tag][ss_id]
        ss_scaled = s["total_cycles"]*BATCH*groups
        axis2 = max(n["comp_cycle"], n["ic_in"], n["ic_w"])
        err = (axis2 - ss_scaled)/ss_scaled*100.0
        err_comp = (n["comp_cycle"] - ss_scaled)/ss_scaled*100.0
        npu_v.append(axis2); ss_v.append(ss_scaled); errs.append(abs(err))
        print(f"{name:7s} {n['comp_cycle']:12,.0f} {n['ic_in']:12,.0f} {n['ic_w']:12,.0f} "
              f"{axis2:12,.0f} {ss_scaled:12,.0f} {err:+7.1f} {err_comp:+10.1f}")
    mape = sum(errs)/len(errs)
    rho = spearman(npu_v, ss_v)
    print(f"\narray-sched = max(comp, IC_in, IC_w): NPUsim's systolic array-schedule proxy")
    print(f"MAPE = {mape:.1f}%   Spearman rho = {rho:.3f}   (n={len(errs)} layers)")

if __name__ == "__main__":
    main()
