#!/usr/bin/env python3
"""Memory-inclusive critical-path latency vs the Gemmini RTL golden (6 GEMM points).

check_timing.py validates the compute schedule (computation + fold fill). The RTL golden is
the end-to-end tiled_matmul_auto time, DMA included, so the critical path -- what a system
model such as MOSAIC consumes -- is compared against it here. Reads result/<accelerator>/
gemm_*/ws/network.txt (run `./npusim.sh run <accelerator> gemm_MxKxN ws` first).

Usage: check_critical_path.py [accelerator=gemmini_rtl] [--max-mape PCT]
"""
import argparse
import csv
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("accelerator", nargs="?", default="gemmini_rtl")
# Regression guard at the current level, not a fidelity claim: gemmini_rtl is at MAPE 47.6%
# (max +111%); the residual is the timeline's single-buffered multi-chip -> GLB boundary and
# the pipeline fill of large tiles (validation/phase2/README.md).
parser.add_argument("--max-mape", type=float, default=50.0)
args = parser.parse_args()

errors = []
for row in csv.DictReader(open(os.path.join(ROOT, "validation", "phase2", "golden_rtl_cycles.csv"))):
    net = f"gemm_{row['m']}x{row['k']}x{row['n']}"
    report = os.path.join(ROOT, "result", args.accelerator, net, "ws", "network.txt")
    text = open(report).read()
    crit = float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", text).group(1))
    rtl = float(row["rtl_cycles"])
    errors.append(crit / rtl - 1)
    print(f"  {net:18} RTL {rtl:9.0f}  critical path {crit:10.0f}  {100 * errors[-1]:+7.1f}%")
mape = 100 * sum(abs(e) for e in errors) / len(errors)
worst = 100 * max(errors, key=abs)
print(f"Gemmini critical path ({args.accelerator}): MAPE={mape:.2f}% max={worst:.2f}% (limit MAPE {args.max_mape:g}%)")
sys.exit(0 if mape <= args.max_mape else 1)
