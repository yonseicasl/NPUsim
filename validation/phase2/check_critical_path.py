#!/usr/bin/env python3
"""gemmini_rtl vs the Gemmini RTL golden (6 GEMM points): compute schedule and the memory-inclusive
critical path.

check_timing.py validates gemmini.cfg's compute schedule (computation + fold fill). The RTL golden is
the end-to-end tiled_matmul_auto time, DMA included, so the critical path -- what a system model such
as MOSAIC consumes -- is compared against it as well. gemmini_rtl uses the RTL's tiling and its own
fold-fill fit, so both metrics are reported for it. Reads result/<accelerator>/gemm_*/ws/network.txt
(run `./npusim.sh run <accelerator> gemm_MxKxN ws` first).

Usage: check_critical_path.py [accelerator=gemmini_rtl] [--max-sched PCT] [--max-crit PCT]
"""
import argparse
import csv
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("accelerator", nargs="?", default="gemmini_rtl")
# Regression guards at gemmini_rtl's current level (max |err| 11.8% / 18.9%), not fidelity claims:
# one fill per fold bounds the compute schedule near 12% with the RTL's fold counts (README.md).
parser.add_argument("--max-sched", type=float, default=12.0, help="max |error| of the compute schedule, %%")
parser.add_argument("--max-crit", type=float, default=20.0, help="max |error| of the critical path, %%")
args = parser.parse_args()

err = {"sched": [], "crit": []}
for row in csv.DictReader(open(os.path.join(ROOT, "validation", "phase2", "golden_rtl_cycles.csv"))):
    net = f"gemm_{row['m']}x{row['k']}x{row['n']}"
    text = open(os.path.join(ROOT, "result", args.accelerator, net, "ws", "network.txt")).read()
    sched = float(re.search(r"Compute-schedule latency\s*:\s*([\d.]+)", text).group(1))
    crit = float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", text).group(1))
    rtl = float(row["rtl_cycles"])
    err["sched"].append(sched / rtl - 1)
    err["crit"].append(crit / rtl - 1)
    print(f"  {net:18} RTL {rtl:9.0f}  schedule {sched:9.0f} {100 * err['sched'][-1]:+6.1f}%"
          f"  critical path {crit:10.0f} {100 * err['crit'][-1]:+7.1f}%")
ok = True
for key, label, limit in (("sched", "compute schedule", args.max_sched), ("crit", "critical path", args.max_crit)):
    mape = 100 * sum(map(abs, err[key])) / len(err[key])
    worst = 100 * max(err[key], key=abs)
    ok &= abs(worst) <= limit
    print(f"Gemmini {label} ({args.accelerator}): MAPE={mape:.2f}% max={worst:.2f}% (limit max {limit:g}%)")
sys.exit(0 if ok else 1)
