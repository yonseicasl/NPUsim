#!/usr/bin/env python3
"""A `measure_tpu_v5e1.py --only tailprobe` run: the banked and the single-weight program of each case, side by side.
Reports event durations and compiler facts as recorded; subtracts nothing into a "tail".

Usage: python3 tailprobe_analyze.py RUN_DIR [BALANCE_RUN_DIR]     (the latter: banked vs that run's loop rows)
"""
import csv
import json
import os
import re
import sys


def main(run, balance=None):
    rows = {(r["dtype"], r["M"], r["K"], r["mode"]): r for r in csv.DictReader(open(os.path.join(run, "tailprobe.csv")))}
    bal = {}
    if balance:
        bal = {(r["dtype"], r["M"], r["K"]): float(r["fused_ops_sum_per_iter_us"])
               for r in csv.DictReader(open(os.path.join(balance, "balance.csv"))) if r["schedule"] == "loop"}
    meta = json.load(open(os.path.join(run, "tailprobe_meta.json")))
    print(f"{meta['device_kind']}, JAX {meta['jax']}, backend {meta['backend']}, quick {meta['quick']}")
    print(f"{'case':14s} | {'banked':>8s} {'single':>8s} {'single-banked':>14s} | {'balance run':>11s} | trace check (banked / single) | HLO facts equal | dynamic-slice banked / single")
    for p in csv.DictReader(open(os.path.join(run, "tailprobe_pairs.csv"))):
        k = (p["dtype"], p["M"], p["K"])
        b, s = rows[k + ("banked",)], rows[k + ("single",)]
        ub, us = (float(x["conv_fusion_us_per_event"]) if x["conv_fusion_us_per_event"] else float("nan") for x in (b, s))
        fb, fs = json.loads(b["hlo_facts"]), json.loads(s["hlo_facts"])
        ds = [len(re.findall(r" dynamic-slice\(", open(os.path.join(run, "hlo", f"tailprobe_{k[0]}_{k[1]}x{k[2]}_{m}.txt")).read()))
              for m in ("banked", "single")]
        ref = f"{bal[k]:8.3f} ({100 * (ub / bal[k] - 1):+.2f}%)" if k in bal else ""
        print(f"{k[0]} {k[1]:>3}x{k[2]:<5} | {ub:8.3f} {us:8.3f} {us - ub:+8.3f} ({100 * (us / ub - 1):+.2f}%) | {ref:>11s} | "
              f"{b['trace_check']} / {s['trace_check']} | {fb == fs} | {ds[0]} / {ds[1]}")
        print(f"   tagged tensors, copies, prefetch: banked {fb['tagged']} {fb['copy_start']} {fb['prefetch']}; "
              f"single {fs['tagged']} {fs['copy_start']} {fs['prefetch']}   windows {fb['windows'][0]['iteration_bounds']} x {fb['windows'][0]['kernel_window_bounds']}")
        print(f"   verdict: {p['verdict']}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
