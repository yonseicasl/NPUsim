#!/usr/bin/env python3
"""TPU v5e GEMMs on NPUsim with the analytical [dram] (819 GB/s) vs MOSAIC's HBM2e (Ramulator2) vs measured.

Measured: tpu_v5e1_stream_bench.py's cold bf16 stream_gemm (weights streamed from HBM inside one program).
NPUsim: configs/accelerators/tpu_v5e.cfg with validation/palm's GEMM mapping (output columns over the 4 MXUs in
128-column weight strips); each shape is simulated at 1 and 2 strips and extended linearly in strips, as
compare_v5e.py does. With --mosaic-dram the dense tile streams' DRAM time comes from MOSAIC's Ramulator2
(NPUSIM_MOSAIC_DRAM), default configs/architecture/dram/hbm2e_tpu_v5e.yaml in the MOSAIC tree.

Usage: python3 validation/tpu_v5e/compare_hbm.py [RESULTS_DIR] [--mosaic-dram YAML] [-j 16]   (needs a MOSAIC build)
"""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "palm"))
import palm_decode as palm  # noqa: E402

ROOT, ACCEL, CLOCK_HZ, DIM = palm.ROOT, "tpu_v5e", 1.5e9, 128
palm.MXUS = 4
MOSAIC_HBM = os.path.join(ROOT, "..", "..", "..", "configs", "architecture", "dram", "hbm2e_tpu_v5e.yaml")


def run(job):
    """Critical-path seconds of an M x K x (strips x 4 x 128) GEMM, analytical or on a MOSAIC DRAM."""
    M, K, strips, dram = job
    tmp = tempfile.mkdtemp(prefix="hbm_")
    try:
        cfgroot, rundir = os.path.join(tmp, "configs"), os.path.join(tmp, "run")
        for d in ("accelerators", "networks", "mappings"):
            os.makedirs(os.path.join(cfgroot, d))
        os.makedirs(os.path.join(rundir, "weights"))
        os.symlink(os.path.join(ROOT, "models", "datasets"), os.path.join(rundir, "datasets"))
        shutil.copy(os.path.join(ROOT, "configs", "accelerators", f"{ACCEL}.cfg"), os.path.join(cfgroot, "accelerators"))
        name = f"v5e_{M}x{K}x{strips * palm.MXUS * DIM}"
        palm.gemm_files(cfgroot, rundir, name, ACCEL, M, K, strips * palm.MXUS * DIM)
        env = dict(os.environ, NPUSIM_CONFIG_ROOT=cfgroot)
        env.pop("NPUSIM_MOSAIC_DRAM", None)
        if dram:
            env["NPUSIM_MOSAIC_DRAM"] = os.path.abspath(dram)
        p = subprocess.run([os.path.join(ROOT, "models", "model"), "run", ACCEL, name, "ws"], cwd=rundir, env=env,
                           capture_output=True, text=True)
        if p.returncode:
            sys.exit(f"{name} failed:\n{(p.stdout + p.stderr)[-1500:]}")
        text = open(os.path.join(rundir, f"{ACCEL}_{name}.txt")).read()
        return float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", text).group(1)) / CLOCK_HZ
    finally:
        shutil.rmtree(tmp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results", nargs="?", default=os.path.join(HERE, "results", "2026-11-03_colab_v5e1_followup3"))
    ap.add_argument("--mosaic-dram", default=MOSAIC_HBM)
    ap.add_argument("-j", "--jobs", type=int, default=16)
    a = ap.parse_args()
    meas = [x for x in json.load(open(os.path.join(a.results, "results.json")))["stream_gemm"]
            if x["weights"] == "bf16" and x["cold"]]
    jobs = sorted({(x["M"], x["K"], s, d) for x in meas for s in (1, 2) for d in ("", a.mosaic_dram)})
    with ThreadPoolExecutor(a.jobs) as ex:
        sim = dict(zip(jobs, ex.map(run, jobs)))

    def at(M, K, N, d):
        strips = math.ceil(math.ceil(N / DIM) / palm.MXUS)
        t1, t2 = sim[M, K, 1, d], sim[M, K, 2, d]
        return (t1 + (strips - 1) * (t2 - t1)) * 1e6

    print(f"{'M':>4} {'K':>6} {'N':>6} | {'measured':>8} {'819 GB/s':>8} {'HBM2e':>8} us | {'err 819':>7} {'err HBM2e':>9}")
    errs = {"": [], a.mosaic_dram: []}
    for x in meas:
        M, K, N = x["M"], x["K"], x["N"]
        an, hb = at(M, K, N, ""), at(M, K, N, a.mosaic_dram)
        e1, e2 = 100 * (an / x["us"] - 1), 100 * (hb / x["us"] - 1)
        if M <= 128:
            errs[""].append(abs(e1)); errs[a.mosaic_dram].append(abs(e2))
        print(f"{M:4d} {K:6d} {N:6d} | {x['us']:8.1f} {an:8.1f} {hb:8.1f}    | {e1:+6.1f}% {e2:+8.1f}%")
    print(f"MAPE M <= 128 (n={len(errs[''])}): analytical 819 GB/s {sum(errs['']) / len(errs['']):.1f}%,"
          f" MOSAIC HBM2e {sum(errs[a.mosaic_dram]) / len(errs[a.mosaic_dram]):.1f}%")


if __name__ == "__main__":
    main()
