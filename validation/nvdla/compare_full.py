#!/usr/bin/env python3
"""NVDLA nv_small RTL vs NPUsim: cycle + DBB traffic, over every CRC-passing trace
that NPUsim can express (see gen_nvdla_workload.py for the exclusion list).

Golden side: validation/nvdla/golden_rtl_full.csv (extract_golden.py --write), from
the EMULATION-model Verilator binary (see PROVENANCE.md 2026-08-25 section).

Three comparisons per trace:
  cycle    Two comparables, chosen by what bounds the trace -- the same split the
           Eyeriss silicon comparison uses:
             compute-bound conv (MAC work >> DMA beats): RTL end-to-end vs NPUsim
               Compute-schedule latency. First-principles (ZERO fitted parameters):
               2026-08-25 measurement gives MAPE 6.8% / max -15.3% over the 7
               compute-bound traces, every error negative = the RTL's un-modeled
               per-stripe overhead (pipeline fills, stripe transitions), amortizing
               with size (-0.5% at K=270).
             memory-bound micro traces: Critical-path latency is the comparable.
               2026-08-26: the fixed pipeline constants are now declared in the
               config (launch 28 + DBB latency 32, both structural and exactly
               constant across every trace; conv-pipe drain 77 CALIBRATED on
               dc_1x1x8, n=1) -- so dc_1x1x8 is a CALIBRATION point, excluded
               from every aggregate and marked "cal". dc_4x1x8192's RTL time is
               bounded by per-stripe overhead (open item #3), not memory latency;
               it stays informational.
           No calibration/holdout split exists yet (plan Sec 3.2), so neither axis
           is a gate; numbers are printed for the record.
  logical  NPUsim DRAM serialized bytes per stream vs its own DECLARED encoding
           (mapping-padded C/K, pad-extended input extent) -- must be EXACT; a
           mismatch is machinery error, not modeling scope. Input became exact on
           2026-08-26 when the halo contract learned legacy-row R/S loops (union
           lower-bound correction; see mapping_table_t::input_halo_reuse()). The
           encoding's deviation from the true tensor volume (pad extension + atomic
           padding) is printed beside it as enc/vol so the distortion stays visible.
  physical RTL measured DBB bytes vs the closed-form atomic-padding contract:
             input  = W*H*ceil(C/8)*8
             weight = ceil(K*C*R*S/128)*128
             output = Wo*Ho*ceil(K/8)*8
           This is a LAW of the interface (verified exact on every passing trace);
           NPUsim does not model it yet, so it is reported as contract-vs-measured,
           not NPUsim-vs-measured.
"""
import csv, math, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_nvdla_workload import dims, short, EXCLUDED   # noqa: E402

SHORT_OVERRIDE = {"dc_1x1x8_1x1x8x1_int8_0": "nvdla_dc_1x1x8_1x1x8x1"}

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

def npusim_result(name):
    # The 1x1x8 point is a single [connected] layer (no shaping dummy), so its result
    # is layer_0; every generated conv network carries the dummy and compares layer_1.
    layer = "layer_0.txt" if name == "nvdla_dc_1x1x8_1x1x8x1" else "layer_1.txt"
    path = os.path.join(ROOT, "result/nvdla_small", name, "ws", layer)
    if not os.path.exists(path): return None
    text = open(path).read()
    crit = float(re.search(r"Critical-path latency :\s*([\d.]+)", text).group(1))
    sched = float(re.search(r"Compute-schedule latency :\s*([\d.]+)", text).group(1))
    section = text.split("Multi-chip <-> DRAM transactions", 1)[1]
    tx = re.findall(r"(Input data|Weight|Output data)\s*:\s*\d+/\d+/(\d+)", section)[:3]
    stream = {k: int(v)*8 for k, v in tx}
    return crit, sched, stream

def main():
    golden = {r["trace"]: r for r in csv.DictReader(open(
        os.path.join(ROOT, "validation/nvdla/golden_rtl_full.csv")))}
    print(f"{'trace':<30}{'RTL cyc':>9}{'NPU cyc':>10}{'err%':>8}{'crit%':>12}   "
          f"{'in':>7}{'wt':>7}{'out':>8}  {'enc/vol i/w/o':>15}  {'phys law':>8}")
    cyc_errs = []
    for trace, g in sorted(golden.items()):
        if not trace.startswith("dc_"): continue
        if trace in EXCLUDED or trace == "dc_24x33x55_5x5x55x25_int8_0": continue
        d = dims(trace)
        name = SHORT_OVERRIDE.get(trace, short(trace))
        res = npusim_result(name)
        if res is None:
            print(f"{trace:<30} (no NPUsim result yet)"); continue
        crit, sched, s = res
        rtl = int(g["cycles"])
        # Compute-bound iff the MAC schedule outweighs the total DBB beat count --
        # the same question "which resource bounds the end-to-end time" the RTL answers.
        dbb_beats = (int(g["in_phys"]) + int(g["wt_phys"]) + int(g["out_phys"])) // 8
        compute_bound = sched > dbb_beats
        CALIBRATION = {
            "dc_1x1x8_1x1x8x1_int8_0",        # drain constant (n=1)
            "dc_13x15x64_5x3x64x16_int8_0",   # stripe bubble B fit set --
            "dc_32x26x76_6x3x76x16_int8_0",   #   includes the K16/K270 twin that
            "dc_32x26x76_6x3x76x270_int8_0",  #   fixed the model's Kg-independence
            "dc_8x16x128_3x3x128x32_int8",
        }
        calibration_point = trace in CALIBRATION
        # Out of the stripe model's declared scope (calibrate_stripe_overhead.py): the
        # C-heavy 1x1-kernel stride-3 shape's overhead regime is neither the stripe bubble
        # nor memory latency. Informational, never aggregated.
        out_of_scope = trace == "dc_4x1x8192_1x1x8192x1_int8_0"
        npu_cyc = sched if compute_bound else crit
        err = (npu_cyc - rtl)/rtl*100
        if compute_bound and not calibration_point and not out_of_scope:
            cyc_errs.append(abs(err))
        W,H,C,R,S,K,Wo,Ho = d["W"],d["H"],d["C"],d["R"],d["S"],d["K"],d["Wo"],d["Ho"]
        st = d["sx"] if W > 1 else d["sy"]
        c_map = min(8,C)*math.ceil(C/min(8,C)); k_map = min(8,K)*math.ceil(K/min(8,K))
        enc = {"in": ((Ho-1)*st+R) * ((Wo-1)*st+S) * c_map,
               "wt": R*S*c_map*k_map, "out": Ho*Wo*k_map}
        vol = {"in": W*H*C, "wt": R*S*C*K, "out": Wo*Ho*K}
        law = {"in": W*H*math.ceil(C/8)*8, "wt": math.ceil(K*C*R*S/128)*128,
               "out": Wo*Ho*math.ceil(K/8)*8}
        meas = {"in": int(g["in_phys"]), "wt": int(g["wt_phys"]), "out": int(g["out_phys"])}
        npu = {"in": s.get("Input data",0), "wt": s.get("Weight",0), "out": s.get("Output data",0)}
        marks = {k: ("ok" if npu[k]==enc[k] else f"x{npu[k]/enc[k]:.2f}") for k in enc}
        drift = "/".join(f"{enc[k]/vol[k]:.2f}" for k in ("in","wt","out"))
        law_ok = all(law[k]==meas[k] for k in law)
        tag = ("cal" if calibration_point else
               "i" if out_of_scope else ("s" if compute_bound else "c"))
        crit_err = (crit - rtl)/rtl*100
        print(f"{trace:<30}{rtl:>9}{npu_cyc:>10.0f}{err:>+8.1f}{tag:>4}{crit_err:>+8.1f}   "
              f"{marks['in']:>7}{marks['wt']:>7}{marks['out']:>8}  {drift:>15}  {'exact' if law_ok else 'VIOLATED':>8}")
    if cyc_errs:
        print(f"\nHOLDOUT cycle (compute-schedule vs RTL): MAPE={sum(cyc_errs)/len(cyc_errs):.2f}%  "
              f"max={max(cyc_errs):.2f}%  (n={len(cyc_errs)}; calibrated model -- setup 137 + "
              f"stripe bubble 4.356 -- with parameters fitted ONLY on the rows tagged cal; "
              f"scope: stride-1/dilation-1 direct conv)")

if __name__ == "__main__":
    main()
