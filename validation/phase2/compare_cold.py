#!/usr/bin/env python3
"""NPUsim critical path vs the cold-L2 Gemmini RTL sweep (gemm_sweep_cold.c, golden_rtl_cycles_cold.csv).

Three NPUsim setups per point:
  A  gemmini_rtl                          memory path of the warm sweep (16 B per cycle, no DRAM timing)
  B  gemmini_rtl_mem                      500 MHz, 64-bit L2 <-> SimDRAM AXI port
  C  gemmini_rtl_mem + MOSAIC Ramulator2  B with SimDRAM's DDR3 (ddr3_1333_simdram.yaml) timing the
                                          DRAM streams (NPUsim built with MOSAIC_ROOT)
The warm golden (golden_rtl_cycles.csv: one A/B/C reused by all points) is listed beside the cold one.

Usage: python3 validation/phase2/compare_cold.py [--dram <MOSAIC DRAM yaml>]
"""
import argparse
import csv
import os
import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DRAM = os.path.join(ROOT, "..", "..", "..", "configs", "architecture", "dram", "ddr3_1333_simdram.yaml")


def run(accel, net, dram):
    env = dict(os.environ, NPUSIM_CONFIG_ROOT=os.path.join(ROOT, "configs"))
    env.pop("NPUSIM_MOSAIC_DRAM", None)
    if dram:
        env["NPUSIM_MOSAIC_DRAM"] = os.path.abspath(dram)
    with tempfile.TemporaryDirectory() as d:
        for x in ("datasets", "weights"):
            os.symlink(os.path.join(ROOT, "models", x), os.path.join(d, x))
        p = subprocess.run([os.path.join(ROOT, "models", "model"), "run", accel, net, "ws"], cwd=d, env=env,
                           capture_output=True, text=True)
        if p.returncode:
            raise RuntimeError(f"{accel} {net}: {(p.stdout + p.stderr)[-800:]}")
        text = open(os.path.join(d, f"{accel}_{net}.txt")).read()
    return float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", text).group(1))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dram", default=DEFAULT_DRAM)
    args = ap.parse_args()
    key = lambda r: (int(r["m"]), int(r["k"]), int(r["n"]))
    cold = {key(r): float(r["rtl_cycles"]) for r in csv.DictReader(open(os.path.join(HERE, "golden_rtl_cycles_cold.csv")))}
    warm = {key(r): float(r["rtl_cycles"]) for r in csv.DictReader(open(os.path.join(HERE, "golden_rtl_cycles.csv")))}
    setups = {"A": ("gemmini_rtl", None), "B": ("gemmini_rtl_mem", None), "C": ("gemmini_rtl_mem", args.dram)}
    jobs = [(s, p) for p in cold for s in setups]
    with ThreadPoolExecutor(8) as ex:
        res = dict(zip(jobs, ex.map(lambda j: run(*setups[j[0]][:1], "gemm_%dx%dx%d" % j[1], setups[j[0]][1]), jobs)))
    print(f"{'M x K x N':>16} {'RTL warm':>9} {'RTL cold':>9} | {'A':>9} {'B':>9} {'C':>9} | err vs cold A / B / C")
    errs = {s: [] for s in setups}
    for p in cold:
        for s in setups:
            errs[s].append(res[s, p] / cold[p] - 1)
        w = f"{warm[p]:9.0f}" if p in warm else f"{'-':>9}"
        print(f"{'%dx%dx%d' % p:>16} {w} {cold[p]:9.0f} | " + " ".join(f"{res[s, p]:9.0f}" for s in setups)
              + " | " + " / ".join(f"{100 * errs[s][-1]:+6.1f}%" for s in setups))
    for s in setups:
        e = errs[s]
        print(f"{s}: MAPE {100 * sum(map(abs, e)) / len(e):.1f}%  max {100 * max(e, key=abs):+.1f}%")


if __name__ == "__main__":
    main()
