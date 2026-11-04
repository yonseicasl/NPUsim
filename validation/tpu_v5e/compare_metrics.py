#!/usr/bin/env python3
"""TPU v5e decode GEMMs, NPUsim vs measured on more than latency: MAC count, HBM bytes, effective HBM bandwidth,
MXU utilization.

Measured (tpu_v5e1_stream_bench.py stream_gemm, cold bf16, weights streamed inside one program): time; the work is
the GEMM's -- 2*M*K*N FLOPs, K*N + M*K + M*N bf16 elements of HBM traffic (each operand read once, the output
written once); XProf's hlo_stats confirms these for the traced cases (xprof_stats.py). Derived: effective HBM
bandwidth = bytes / time, MXU utilization = FLOP rate / 197 TFLOP/s.
NPUsim (tpu_v5e.cfg, validation/palm's mapping, 1 and 2 strips extended linearly as compare_v5e.py does): critical
path, busy scalar-MAC cycles (MACs issued, padding included), DRAM transactions x line size (bytes it moves), and the
same derived quantities.

Usage: python3 compare_metrics.py [RESULTS_DIR] [--mosaic-dram YAML] [-j 16] [--dump M,K,STRIPS]
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
PEAK_FLOPS = 4 * DIM * DIM * 2 * CLOCK_HZ          # 197 TFLOP/s
LINE_BYTES = 8192 // 8                             # tpu_v5e.cfg [dram] line_size (bits) / 8


def report(M, K, strips, dram=""):
    """NPUsim report text of an M x K x (strips x 4 x 128) GEMM."""
    tmp = tempfile.mkdtemp(prefix="met_")
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
        return open(os.path.join(rundir, f"{ACCEL}_{name}.txt")).read()
    finally:
        shutil.rmtree(tmp)


def metrics(job):
    """(cycles, busy MACs, DRAM bytes by type) of one NPUsim run."""
    M, K, s, dram = job
    text = report(M, K, s, dram)
    num = lambda pat: float(re.search(pat, text).group(1))
    dram_bytes = {}
    sect = text[text.index("Multi-chip <-> DRAM transactions"):]
    for t in ("Input data", "Weight", "Output data"):
        dram_bytes[t] = float(re.search(rf"\* {t}\s*:\s*(\d+)/", sect).group(1)) * LINE_BYTES
    return {"cycles": num(r"Critical-path latency\s*:\s*([\d.]+)"),
            "macs": num(r"Busy scalar-MAC cycles\s*:\s*([\d.]+)"), **dram_bytes}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results", nargs="?", default=os.path.join(HERE, "results", "2026-11-03_colab_v5e1_followup3"))
    ap.add_argument("--mosaic-dram", default="")
    ap.add_argument("-j", "--jobs", type=int, default=16)
    ap.add_argument("--dump")
    a = ap.parse_args()
    if a.dump:
        M, K, s = (int(v) for v in a.dump.split(","))
        print(report(M, K, s, a.mosaic_dram))
        return
    r = json.load(open(os.path.join(a.results, "results.json")))
    meas = [x for x in r["stream_gemm"] + r["odd_gemm"] if x["weights"] == "bf16" and x["cold"]]
    jobs = sorted({(x["M"], x["K"], s, a.mosaic_dram) for x in meas for s in (1, 2)})
    with ThreadPoolExecutor(a.jobs) as ex:
        sim = dict(zip(jobs, ex.map(metrics, jobs)))

    def at(M, K, N):
        """NPUsim quantities at N, linear in weight strips between the 1- and 2-strip runs."""
        strips = math.ceil(math.ceil(N / DIM) / palm.MXUS)
        m1, m2 = sim[M, K, 1, a.mosaic_dram], sim[M, K, 2, a.mosaic_dram]
        return {k: m1[k] + (strips - 1) * (m2[k] - m1[k]) for k in m1}

    rows = []
    print(f"{'M':>4} {'K':>6} {'N':>6} | {'time us':>15} | {'MACs N/real':>11} {'bytes N/real':>12} |"
          f" {'eff. HBM GB/s':>13} | {'MXU util %':>11}")
    print(f"{'':>18} | {'meas':>7} {'NPUsim':>7} | {'':>11} {'':>12} | {'meas':>6} {'NPUsim':>6} | {'meas':>5} {'NPUsim':>5}")
    for x in meas:
        M, K, N = x["M"], x["K"], x["N"]
        n = at(M, K, N)
        t_meas, t_sim = x["us"] * 1e-6, n["cycles"] / CLOCK_HZ
        real_bytes = 2 * (K * N + M * K + M * N)
        sim_bytes = n["Input data"] + n["Weight"] + n["Output data"]
        row = dict(M=M, K=K, N=N, meas_us=x["us"], npusim_us=t_sim * 1e6, time_err_pct=100 * (t_sim / t_meas - 1),
                   mac_ratio=n["macs"] / (M * K * N), byte_ratio=sim_bytes / real_bytes,
                   meas_hbm_gbps=real_bytes / t_meas / 1e9, npusim_hbm_gbps=sim_bytes / t_sim / 1e9,
                   meas_mxu_pct=100 * 2 * M * K * N / t_meas / PEAK_FLOPS,
                   npusim_mxu_pct=100 * 2 * n["macs"] / t_sim / PEAK_FLOPS)
        rows.append(row)
        print(f"{M:4d} {K:6d} {N:6d} | {row['meas_us']:7.1f} {row['npusim_us']:7.1f} | {row['mac_ratio']:11.3f}"
              f" {row['byte_ratio']:12.3f} | {row['meas_hbm_gbps']:6.0f} {row['npusim_hbm_gbps']:6.0f} |"
              f" {row['meas_mxu_pct']:5.1f} {row['npusim_mxu_pct']:5.1f}")
    out = os.path.join(a.results, "compare_metrics" + ("_mosaic" if a.mosaic_dram else "") + ".csv")
    with open(out, "w", newline="") as fh:
        import csv
        w = csv.DictWriter(fh, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
