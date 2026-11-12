#!/usr/bin/env python3
"""A `measure_tpu_v5e1.py --only balance` run: stream / compute balance for bf16 and int8 weights.

Loop rows (the stream bench's schedule) against three predictions whose constants were ALL fixed before the run
(int8_decompose.py: R = 752 GB/s, F = 2 us, w = 0.393 us from the 2026-11-04 int8 M = 128 cases):
  max+F     perfect overlap of weight streaming and compute
  pipe      n XLA windows (n from the run's own optimized HLO), each DMAed whole before its compute
  pipe+w    the same with a compute-side cost w per window
Isolated rows are a different schedule (operand copied first): kernel = device - copy, k-C its excess over the
compute roofline, eta = (copy + kernel - fused) / min(copy, kernel): the share of the shorter phase that streaming
hides. No simulator involved.

Usage: python3 balance_analyze.py RUN_DIR
"""
import collections
import csv
import os
import re
import statistics as st
import sys

from int8_decompose import F, PEAK, R, pipeline

W = 0.393


def windows(path):
    """Window count of the program's one convolution, from its window_config."""
    it = re.findall(r'"iteration_bounds":\[([^\]]*)\]', open(path).read())
    assert len(it) == 1, (path, it)
    n = 1
    for x in it[0].split(","):
        n *= int(x.strip('"'))
    return n


def main(run):
    rows = list(csv.DictReader(open(os.path.join(run, "balance.csv"))))
    iso = {(r["dtype"], r["M"], r["K"]): r for r in rows if r["schedule"] == "isolated"}
    err = collections.defaultdict(list)
    print(f"{'dtype':5} {'M':>4} {'K':>5} {'n':>3} {'C/S':>5} | {'slope':>7} {'fused':>7} {'/max':>5} |"
          f" {'max+F':>6} {'pipe':>6} {'pipe+w':>6} | {'kernel':>7} {'k-C':>6} {'eta':>5}")
    for r in rows:
        if r["schedule"] != "loop":
            continue
        dt, M, K = r["dtype"], int(r["M"]), int(r["K"])
        n = windows(os.path.join(run, "hlo", f"balance_{dt}_{M}x{K}_loop.txt"))
        S, C = K * K * (1 if dt == "int8" else 2) / R * 1e6, 2 * M * K * K / PEAK * 1e6
        f = float(r["fused_ops_sum_per_iter_us"])       # this run: one fused op per iteration (the ops column)
        pred = {"max+F": max(S, C) + F, "pipe": pipeline(n, S, C, 0.0), "pipe+w": pipeline(n, S, C, W)}
        regime = "C/S < 0.6" if C / S < 0.6 else "C/S 0.6-1.6" if C / S < 1.6 else "C/S > 1.6"
        for k, v in pred.items():
            for g in ("all", dt, regime):
                err[k, g].append(v / f - 1)
        i = iso[dt, r["M"], r["K"]]
        dev, cp = float(i["device_us"]), float(i["copy_us"])
        staged = f"{dev - cp:7.2f} {dev - cp - C:+6.2f} {(dev - f) / min(cp, dev - cp):5.2f}" if cp > 1 else "  (operand not staged)"
        print(f"{dt:5} {M:4d} {K:5d} {n:3d} {C / S:5.2f} | {float(r['loop_slope_us']):7.2f} {f:7.2f} {f / max(S, C):5.2f} | "
              + " ".join(f"{100 * (v / f - 1):+5.1f}%" for v in pred.values()) + f" | {staged}")
    print(f"\n{'error vs fused (mean / MAPE)':28s} | " + " | ".join(f"{k:>14s}" for k in ("max+F", "pipe", "pipe+w")))
    for g in ("all", "bf16", "int8", "C/S < 0.6", "C/S 0.6-1.6", "C/S > 1.6"):
        print(f"{g:24s} n={len(err['pipe', g]):2d} | " + " | ".join(
            f"{100 * st.mean(err[k, g]):+5.1f} / {100 * st.mean(map(abs, err[k, g])):4.1f}%" for k in ("max+F", "pipe", "pipe+w")))
    fit = os.path.join(run, "loopk_fit.csv")
    if os.path.exists(fit):
        d = [(r["M"], r["N"], float(r["per_kernel_us"])) for r in csv.DictReader(open(fit)) if r["kind"] == "loop_minus_unrolled"]
        print("\nloopk: what a fori_loop iteration adds to an unrolled kernel (us): " + ", ".join(f"{m}x{n}: {v:+.2f}" for m, n, v in d))


if __name__ == "__main__":
    main(sys.argv[1])
