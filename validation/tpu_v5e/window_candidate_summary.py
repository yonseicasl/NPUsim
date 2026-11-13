#!/usr/bin/env python3
"""Exact A / B / C error tables of a balance_npusim.py batch: errors recomputed, unrounded, from each run's own
report (critical-path cycles) and the unrounded measurement -- not from the CSV's one-decimal percentages.

Usage: python3 window_candidate_summary.py BATCH_DIR BALANCE_RUN_DIR [TAG_SUFFIX]     (TAG_SUFFIX e.g. _fusedout)
"""
import csv
import os
import re
import statistics as st
import sys

PEAK, R, CLOCK = 4 * 128 * 128 * 2 * 1.5e9, 750e9, 1.5e9


def load(batch, run, tag, suffix):
    meas = {(r["dtype"], r["M"], r["K"]): float(r["fused_ops_sum_per_iter_us"])
            for r in csv.DictReader(open(os.path.join(run, "balance.csv"))) if r["schedule"] == "loop"}
    rows = list(csv.DictReader(open(os.path.join(batch, f"balance_npusim_{tag}{suffix}.csv"))))
    for r in rows:
        k = (r["dtype"], r["M"], r["K"])
        for v in "ABC":
            a = f"bal_{r['dtype']}_{r['M']}x{r['K']}_{v}".lower()
            t = open(os.path.join(batch, f"generated_{tag}{suffix}", a, f"{a}_{a}_layer_0.txt")).read()
            us = float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", t).group(1)) / CLOCK * 1e6
            r[v] = 100 * (us / meas[k] - 1)
        M, K = int(r["M"]), int(r["K"])
        cs = (2 * M * K * K / PEAK) / (K * K * (1 if r["dtype"] == "int8" else 2) / R)
        r["regime"] = "memory-bound C/S<0.6" if cs < 0.6 else "balanced C/S 0.6-1.6" if cs < 1.6 else "compute-bound C/S>1.6"
        r["win"] = "full windows" if float(r["window_coverage"]) <= 1 else "partial edge windows"
    return rows


def main(batch, run, suffix=""):
    print("errors in %, unrounded, from the runs' reports; mean / MAPE; regime = (2MKN / 196.6 TFLOPS) / (weight bytes / 750 GB/s)\n")
    for tag in ("hbm2e_tpu_v5e_refsb", "analytic"):
        rows = load(batch, run, tag, suffix)
        print(f"== DRAM model: {tag}  ({len(rows)} cases)")
        print(f"{'group':32s}  n | {'A':>15s} | {'B':>15s} | {'C':>15s} | |err| C<A / C>A")
        for key in (lambda r: "all", lambda r: r["dtype"], lambda r: r["regime"], lambda r: r["dtype"] + ", " + r["regime"],
                    lambda r: r["win"]):
            groups = {}
            for r in rows:
                groups.setdefault(key(r), []).append(r)
            for g in sorted(groups):
                rs = groups[g]
                cell = lambda v: f"{st.mean(r[v] for r in rs):+6.2f} / {st.mean(abs(r[v]) for r in rs):5.2f}"
                b = sum(abs(r["C"]) < abs(r["A"]) for r in rs)
                print(f"{g:32s} {len(rs):2d} | {cell('A'):>15s} | {cell('B'):>15s} | {cell('C'):>15s} | {b} / {len(rs) - b}")
            print()
        print("regressions (|err| C > |err| A): signed A -> C, increase of |err| in points")
        for r in rows:
            if abs(r["C"]) > abs(r["A"]):
                print(f"   {r['dtype']} {r['M']}x{r['K']}: {r['A']:+.3f} -> {r['C']:+.3f}   +{abs(r['C']) - abs(r['A']):.3f}")
        print("worst C: " + "; ".join(f"{r['dtype']} {r['M']}x{r['K']} {r['C']:+.2f}" for r in sorted(rows, key=lambda r: -abs(r["C"]))[:6]) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:4])
