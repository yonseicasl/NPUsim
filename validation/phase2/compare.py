#!/usr/bin/env python3
"""Phase-2 RTL validation: NPUsim (gemmini.cfg, 16x16 WS) vs Gemmini Verilator
RTL cycle counts on identical GEMM points.

NPUsim comparable quantity: the systolic compute schedule =
Computation cycle + Fold fill cycle (V2: RTL-calibrated weight-residency fold
fill + per-layer setup, produced directly by the simulator). RTL cycles come
from gemm_sweep.c's read_cycles() deltas (validation/phase2/rtl_cycles.txt).
"""
import os, re

ROOT = os.path.dirname(os.path.abspath(__file__))
POINTS = [(64,64,64),(128,128,128),(256,256,256),(16,512,512),(512,512,64),(512,512,512)]

RTL_CYCLES = {}
for line in open(os.path.join(ROOT, "rtl_cycles.txt")):
    m = re.search(r"M=(\d+) K=(\d+) N=(\d+) cycles=(\d+)", line)
    if m:
        RTL_CYCLES[(int(m.group(1)), int(m.group(2)), int(m.group(3)))] = int(m.group(4))

def parse_npusim(m, k, n):
    path = os.path.join(ROOT, f"../../result/gemmini/gemm_{m}x{k}x{n}/ws/layer_0.txt")
    txt = open(path).read()
    comp = float(re.search(r"Computation cycle\s*:\s*([\d.]+)", txt).group(1))
    ff = re.search(r"Fold fill cycle\s*:\s*([\d.]+)", txt)
    fold_fill = float(ff.group(1)) if ff else 0.0
    return comp, fold_fill

def main():
    print(f"{'M,K,N':>13s} {'NPU comp':>10s} {'fold fill':>10s} {'NPU sched':>10s} "
          f"{'RTL':>9s} {'err%':>7s}")
    errs = []
    for m, k, n in POINTS:
        comp, fold_fill = parse_npusim(m, k, n)
        rtl = RTL_CYCLES.get((m, k, n))
        if not rtl:
            continue
        sched = comp + fold_fill
        e = (sched - rtl)/rtl*100
        errs.append(abs(e))
        print(f"{m:4d},{k:4d},{n:4d} {comp:10,.0f} {fold_fill:10,.0f} {sched:10,.0f} "
              f"{rtl:9,d} {e:+7.1f}")
    print(f"\nMAPE = {sum(errs)/len(errs):.1f}%   max |err| = {max(errs):.1f}%  (n={len(errs)})")
    print("NPU sched = Computation cycle + Fold fill cycle "
          "(V2: 14 cyc/weight-residency fold + 2270 cyc/layer setup, RTL-calibrated)")

if __name__ == "__main__":
    main()
