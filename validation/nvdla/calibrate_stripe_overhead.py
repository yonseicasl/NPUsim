#!/usr/bin/env python3
"""Per-stripe overhead: fit the one bubble parameter on a declared calibration set,
evaluate on held-out traces.

MODEL. After the pipeline constants (layer_setup_cycle = 137) the compute-bound
residual RTL - MACbound - setup clusters at 3.4-4.8 cycles per legacy
(C-group, r, s) reduction slice per output row on every stride-1/no-dilation conv
trace, and is INDEPENDENT of the kernel-group count: the 32x26x76 K16/K270 twin
changes weight blocks 17x with a near-identical residual (22,210 vs 22,898),
killing every per-weight-block form. So:

    overhead = B * ceil(C/8) * R * S * Ho

FORM PROVENANCE (stated, not hidden): the functional form was identified from the
full trace set's structure -- in particular the K16/K270 twin -- so BOTH twins are
placed in the calibration set. Only the parameter B is fitted (least squares
through the origin), and only on the calibration traces.

SCOPE. Valid for stride-1, dilation-1 direct convs. The non-unit-stride/dilation
traces put 9.4-29.4 cycles on the same axis (different stripe machinery) and the
degenerate shapes (1x1 kernel C-heavy, H=1) sit in other regimes entirely; all
stay OUT of both sets and OUT of the model's claimed scope.
"""
import csv, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_nvdla_workload import dims   # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SETUP = 137.0

CALIBRATION = ["dc_13x15x64_5x3x64x16_int8_0", "dc_32x26x76_6x3x76x16_int8_0",
               "dc_32x26x76_6x3x76x270_int8_0", "dc_8x16x128_3x3x128x32_int8"]
HOLDOUT = ["dc_14x7x49_3x4x49x32_int8_0", "dc_24x44x14_5x3x14x41_int8_0",
           "dc_35x22x54_6x8x54x29_int8_0"]

def slices_rows(trace):
    d = dims(trace)
    Cg = math.ceil(d["C"]/8); Kg = math.ceil(d["K"]/8)
    mac = d["Wo"]*d["Ho"]*d["R"]*d["S"]*Cg*Kg
    return Cg*d["R"]*d["S"]*d["Ho"], mac

def main():
    golden = {r["trace"]: int(r["cycles"]) for r in csv.DictReader(open(
        os.path.join(ROOT, "validation/nvdla/golden_rtl_full.csv")))}
    num = den = 0.0
    print("calibration set:")
    rows = []
    for t in CALIBRATION:
        n, mac = slices_rows(t)
        over = golden[t] - mac - SETUP
        rows.append((t, n, mac, over))
        num += over*n; den += n*n
        print(f"  {t:<32} slices*rows={n:>6}  overhead={over:>9.0f}  per={over/n:.2f}")
    B = num/den
    print(f"\nB (least squares through origin) = {B:.3f} cycles per (slice, output-row)")
    print("\ncalibration-set fit:")
    for t, n, mac, over in rows:
        pred = mac + SETUP + B*n
        print(f"  {t:<32} pred={pred:>9.0f}  RTL={golden[t]:>9}  err={(pred-golden[t])/golden[t]*100:>+6.2f}%")
    print("\nHOLDOUT evaluation (parameter never saw these):")
    errs = []
    for t in HOLDOUT:
        n, mac = slices_rows(t)
        pred = mac + SETUP + B*n
        e = (pred - golden[t])/golden[t]*100
        errs.append(abs(e))
        print(f"  {t:<32} pred={pred:>9.0f}  RTL={golden[t]:>9}  err={e:>+6.2f}%")
    print(f"\nholdout: MAPE={sum(errs)/len(errs):.2f}%  max={max(errs):.2f}%  (n={len(errs)})")
    return B

if __name__ == "__main__":
    main()
