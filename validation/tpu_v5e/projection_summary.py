#!/usr/bin/env python3
"""A / B / C tables of a `balance_npusim.py --manifest` batch from its unrounded per-run values
(balance_npusim_<tag>_exact.json), in the groups fixed before the batch. Retrospective cohort: see the manifest.

Usage: python3 projection_summary.py BATCH_DIR MANIFEST.csv
"""
import csv
import json
import os
import statistics as st
import sys

PEAK, R = 4 * 128 * 128 * 2 * 1.5e9, 750e9


def main(batch, manifest):
    man = {(r["dtype"], int(r["M"]), int(r["K"]), int(r["N"])): r for r in csv.DictReader(open(manifest))}
    print("errors in % against the mean convolution-event duration of each case's ONE profiled invocation (8 correlated")
    print("events); unrounded; mean / MAPE. Retrospective, previously inspected cohort -- no generalization claim.")
    print("regime = (2MKN / 196.6 TFLOPS) / (weight bytes / 750 GB/s); direction = which dimensions XLA's windows split.\n")
    for tag in ("hbm2e_tpu_v5e_refsb_fusedout", "analytic_fusedout"):
        cases = json.load(open(os.path.join(batch, f"balance_npusim_{tag}_exact.json")))["cases"]
        for c in cases:
            m = man[c["dtype"], c["M"], c["K"], c["N"]]
            it = [int(x) for x in m["iteration_bounds"].split("x")]          # [N windows, M windows, K windows]
            cs = (2 * c["M"] * c["K"] * c["N"] / PEAK) / (c["K"] * c["N"] * (1 if c["dtype"] == "int8" else 2) / R)
            c.update(cal="w-calibration cases" if c["w_calibration_case"] else "other cases",
                     regime="memory-bound C/S<0.6" if cs < 0.6 else "balanced C/S 0.6-1.6" if cs < 1.6 else "compute-bound C/S>1.6",
                     win="full windows" if c["window_coverage"] <= 1 else "partial edge windows",
                     direction="windows over N" if it[2] == 1 else "windows over K" if it[0] == 1 else "windows over N and K",
                     shape=f"{c['K']}x{c['N']}")
        print(f"== DRAM model: {tag}  ({len(cases)} cases)")
        print(f"{'group':40s}  n | {'A':>15s} | {'B':>15s} | {'C':>15s} | |err| B<A  C<A  C<B")
        keys = [lambda c: "all 32 (calibration cases included)", lambda c: c["cal"], lambda c: "other cases, " + c["dtype"] if not c["w_calibration_case"] else None,
                lambda c: "other cases, " + c["regime"] if not c["w_calibration_case"] else None,
                lambda c: "other cases, " + c["win"] if not c["w_calibration_case"] else None,
                lambda c: "other cases, " + c["direction"] if not c["w_calibration_case"] else None,
                lambda c: "other cases, " + c["shape"] if not c["w_calibration_case"] else None]
        for key in keys:
            groups = {}
            for c in cases:
                if key(c):
                    groups.setdefault(key(c), []).append(c)
            for g in sorted(groups):
                cs_ = groups[g]
                cell = lambda v: f"{st.mean(c[v]['err_pct'] for c in cs_):+6.2f} / {st.mean(abs(c[v]['err_pct']) for c in cs_):5.2f}"
                closer = lambda x, y: sum(abs(c[x]['err_pct']) < abs(c[y]['err_pct']) for c in cs_)
                print(f"{g:40s} {len(cs_):2d} | {cell('A'):>15s} | {cell('B'):>15s} | {cell('C'):>15s} | {closer('B','A'):>6d} {closer('C','A'):>4d} {closer('C','B'):>4d}")
            print()
        print("cases where C is farther than A (signed A -> B -> C, increase of |err| in points):")
        for c in sorted(cases, key=lambda c: abs(c["A"]["err_pct"]) - abs(c["C"]["err_pct"])):
            d = abs(c["C"]["err_pct"]) - abs(c["A"]["err_pct"])
            if d > 0:
                print(f"   {c['dtype']} {c['M']:>3}x{c['K']}x{c['N']:<5} n={c['windows']:<3} {c['regime'][:13]:13s} {c['direction']:20s}"
                      f" {c['A']['err_pct']:+7.3f} -> {c['B']['err_pct']:+7.3f} -> {c['C']['err_pct']:+7.3f}   +{d:.3f}{'  (w-calibration)' if c['w_calibration_case'] else ''}")
        print("largest |err| of C: " + "; ".join(f"{c['dtype']} {c['M']}x{c['K']}x{c['N']} {c['C']['err_pct']:+.2f}"
                                               for c in sorted(cases, key=lambda c: -abs(c["C"]["err_pct"]))[:6]) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:3])
