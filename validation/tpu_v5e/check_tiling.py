#!/usr/bin/env python3
"""compare_v5e.py's tiled estimate against a direct NPUsim run of the whole GEMM (no TPU data involved).

compare_v5e.py simulates one <= 512-row tile at 2 and 4 weight strips (the first tile with the layer setup,
the others without) and scales; here the same GEMM is simulated whole -- every strip as a DRAM iteration, row
tiles as the DRAM-level batch factor -- on the same config with program_setup_cycle 0.

Usage: python3 validation/tpu_v5e/check_tiling.py [M,K,N ...] [-j 4]   (needs models/model)
"""
import argparse
import math
import os
import re
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_v5e as c
palm = c.palm


def direct(shapes, jobs):
    tmp, cfgroot, rundir = palm.workdir((c.ACCEL,))
    p = os.path.join(cfgroot, "accelerators", f"{c.ACCEL}.cfg")
    text = open(p).read()
    open(p, "w").write(re.sub(r"(?m)^program_setup_cycle\s*=.*$", "program_setup_cycle = 0", text))
    todo = {}
    for M, K, N in shapes:
        mt = math.ceil(M / c.MTILE)
        assert M % mt == 0 and N % (palm.MXUS * c.DIM) == 0, "whole row tiles and strips only"
        name = f"whole_{M}x{K}x{N}"
        palm.gemm_files(cfgroot, rundir, name, c.ACCEL, M, K, N)
        if mt > 1:                            # rows: <= 512 per fold, the rest as DRAM-level batch tiles
            mf = os.path.join(cfgroot, "mappings", c.ACCEL, name, "ws.map")
            s = open(mf).read()
            s = s.replace(f"PE      = 1, {M},", f"PE      = 1, {M // mt},")
            s, n = re.subn(r"DRAM    = (\d+), 1,", lambda m: f"DRAM    = {m.group(1)}, {mt},", s)
            assert n == 1
            open(mf, "w").write(s)
        todo[M, K, N] = (cfgroot, rundir, name, c.ACCEL)
    with ThreadPoolExecutor(jobs) as ex:
        out = dict(zip(todo, ex.map(palm.run_gemm, todo.values())))
    shutil.rmtree(tmp)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("shapes", nargs="*", default=["512,2048,4096", "512,4096,4096", "1024,1024,1024", "128,4096,4096"])
    ap.add_argument("-j", "--jobs", type=int, default=4)
    a = ap.parse_args()
    shapes = [tuple(int(x) for x in s.split(",")) for s in a.shapes]
    est, _ = c.simulate(shapes, a.jobs)
    whole = direct(shapes, a.jobs)
    worst = 0.0
    for k in shapes:
        d = 100 * (est[k] / whole[k] - 1)
        worst = max(worst, abs(d))
        print(f"{k[0]:5d} x {k[1]:5d} x {k[2]:5d}  tiled {est[k] * 1e6:9.2f} us  whole {whole[k] * 1e6:9.2f} us  {d:+5.2f}%")
    print(f"max |tiled - whole| {worst:.2f}%")


if __name__ == "__main__":
    main()
