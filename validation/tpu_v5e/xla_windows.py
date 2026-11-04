#!/usr/bin/env python3
"""XLA's tiling of each GEMM: the matmul fusion's window_config in the saved optimized HLO (RESULTS/hlo/<case>.txt).

iteration_bounds = windows per dimension (their product = how many tiles the kernel is pipelined over),
kernel_window_bounds = one weight window in (8, 128)-tile units, buffering_level = 2 for double buffering.
Writes RESULTS/xla_windows.csv. Usage: python3 xla_windows.py RESULTS
"""
import csv
import glob
import json
import math
import os
import re
import sys


def window(path):
    """window_config of the matmul fusion (the one with a convolution_algorithm_config) in an HLO text."""
    for line in open(path):
        if "convolution_algorithm_config" not in line or "window_config" not in line:
            continue
        cfg = json.loads(re.search(r"backend_config=(\{.*\})", line).group(1).split(", metadata=")[0])
        w = cfg["window_config"]
        vmem = sum(int(m["size"]) for m in cfg.get("used_scoped_memory_configs", []))
        return w, vmem
    return None, 0


def main():
    res = sys.argv[1]
    rows = []
    for path in sorted(glob.glob(os.path.join(res, "hlo", "gemm_*.txt"))):
        case = os.path.basename(path)[:-4]
        w, vmem = window(path)
        if w is None:
            continue
        M, K, N = (int(v) for v in case.split("_")[1].split("x"))
        it = [int(v) for v in w.get("iteration_bounds", [])]
        rows.append(dict(case=case, M=M, K=K, N=N, windows=math.prod(it) if it else 1,
                         iteration_bounds="x".join(map(str, it)),
                         kernel_window="x".join(w.get("kernel_window_bounds", [])),
                         input_window="x".join(w.get("input_window_bounds", [])),
                         buffering=w.get("buffering_level"), xla_estimated_cycles=w.get("estimated_cycles"),
                         scoped_vmem_bytes=vmem))
    out = os.path.join(res, "xla_windows.csv")
    with open(out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    for r in rows:
        if r["case"].endswith("_bf16_cold"):
            print(f"{r['M']:4d} {r['K']:6d} {r['N']:6d} windows {r['windows']:4d} ({r['iteration_bounds']:>8}) kernel "
                  f"{r['kernel_window']:>9} input {r['input_window']:>7} buf {r['buffering']} vmem {r['scoped_vmem_bytes'] / 2**20:.1f} MB")
    print(f"wrote {out} ({len(rows)} GEMM cases)")


if __name__ == "__main__":
    main()
