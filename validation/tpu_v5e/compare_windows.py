#!/usr/bin/env python3
"""TPU v5e GEMMs on NPUsim with XLA's own tiling: latency and the multi-metric view, three setups.

  A  validation/palm's mapping (one weight strip of 4 MXUs x 128 columns per DRAM tile, 1- and 2-strip runs
     extended linearly), layer_setup_cycle = 0 -- what compare_metrics.py ran
  B  A with layer_setup_cycle = 330 (0.22 us per kernel, measured: tpu_v5e.cfg)
  C  XLA's tiling (xla_windows.py: each GEMM's matmul window_config in the run's optimized HLO). NPUsim simulates
     ONE window -- its K rows and 128-column blocks, the blocks on the 4 MXUs in XLA's number of rounds with the
     fewest padded blocks (128-column granularity: 9 blocks = 3 MXUs x 3 rounds) -- giving its DRAM time D and
     on-chip time C (the slowest of NoP, GLB, PE array, PE). The kernel runs XLA's way: a window is DMAed whole
     before its compute, the next one is fetched meanwhile (buffering_level 2):
         t = 330 cycles + D + (n - 1) max(D, C) + C,   n = windows along N x windows along K
     MACs and weight bytes are per window x n; the input is read once per K window (XLA: HBM bytes = weights +
     input). Composed from NPUsim's per-stage costs, not fitted. Finding: putting a window's column blocks alone
     on the MXUs leaves MXUs idle for narrow windows (M 512: 2-block windows -> 2 of 4 MXUs, +84%); XLA keeps all
     four busy over another dimension.
  D  the two concerns apart. Work (1, 128-column granularity): validation/palm's mapping on all 4 MXUs, the 1- and
     2-strip runs extended to ceil(N / 128) / 4 strips (fractional, so only 128-column padding) at setup 0.
     Pipeline (3): that GEMM's per-stage totals split over XLA's n windows, run XLA's way (above):
         t = 330 + D/n + (n - 1) max(D/n, C/n) + C/n,   D = DRAM busy, C = the slowest on-chip stage's busy.

Measured: the same run's stream_gemm / odd_gemm (cold bf16); its work is XLA's (HBM bytes = weights + input,
xprof_stats.py; the bench reduces the output inside the fusion). NPUsim bytes are its input + weight DRAM traffic.

Usage: python3 compare_windows.py [RESULTS_DIR] [-j 16]   (RESULTS_DIR holds results.json and xla_windows.csv)
"""
import argparse
import csv
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "palm"))
import palm_decode as palm  # noqa: E402
from compare_metrics import ACCEL, CLOCK_HZ, DIM, LINE_BYTES, PEAK_FLOPS, ROOT  # noqa: E402

palm.MXUS = MXUS = 4


def window_shape(M, K, N, nb_win, kw_win):
    """One XLA window of an M x K x N GEMM: (K rows, column blocks, windows along N, windows along K)."""
    blocks = math.ceil(N / DIM)
    nb, kw = min(nb_win, blocks), min(kw_win, K)
    return kw, nb, math.ceil(blocks / nb), math.ceil(K / kw)


def xla_mapping(path, M, kw, nb):
    """Rewrite the ws.map palm.gemm_files wrote for one window (M x kw x nb*128): its blocks on the MXUs in XLA's
    number of rounds with the fewest padded blocks."""
    rounds = math.ceil(nb / MXUS)
    mxus = math.ceil(nb / rounds)
    rows = {"MAC": [1, 1], "PE": [1, M], "PE_X": [DIM, 1], "PE_Y": [1, 1], "GLB": [rounds, 1],
            "CHIPS_X": [1, 1], "CHIPS_Y": [mxus, 1], "DRAM": [1, 1]}
    c = {"PE_Y": DIM, "GLB": math.ceil(kw / DIM)}
    N, K = nb * DIM, kw
    with open(path, "w") as f:
        f.write(f"# K, B, P, Q, C, R, S, H, W, GROUP, STRIDE\n\n[connected]\n# {N}, {M}, 1, 1, {K}, 1, 1, 0, 0, 1, 1,\n")
        for lv, (k, b) in rows.items():
            f.write(f"{lv:7s} = {k}, {b}, 1, 1, {c.get(lv, 1)}, 1, 1, 0, 0, 1, 1,\n")


STAGES = ("DRAM", r"Multi-chip \(NoP\)", "Global buffer", "PE array", r"PE \(compute\+LB\)")


def run(job):
    """(cycles, busy MACs, input bytes, weight bytes, stage busy cycles) of one NPUsim run; job = (M, K, N, setup,
    window or None). With a window, K and N are the window's (kw, nb x 128)."""
    M, K, N, setup, win = job
    tmp = tempfile.mkdtemp(prefix="win_")
    try:
        cfgroot, rundir = os.path.join(tmp, "configs"), os.path.join(tmp, "run")
        for d in ("accelerators", "networks", "mappings"):
            os.makedirs(os.path.join(cfgroot, d))
        os.makedirs(os.path.join(rundir, "weights"))
        os.symlink(os.path.join(ROOT, "models", "datasets"), os.path.join(rundir, "datasets"))
        cfg = open(os.path.join(ROOT, "configs", "accelerators", f"{ACCEL}.cfg")).read()
        cfg = re.sub(r"(?m)^layer_setup_cycle\s*=.*$", f"layer_setup_cycle = {setup}", cfg)
        open(os.path.join(cfgroot, "accelerators", f"{ACCEL}.cfg"), "w").write(cfg)
        name = f"g{M}x{K}x{N}"
        palm.gemm_files(cfgroot, rundir, name, ACCEL, M, K, N)
        if win:
            xla_mapping(os.path.join(cfgroot, "mappings", ACCEL, name, "ws.map"), M, K, N // DIM)
        p = subprocess.run([os.path.join(ROOT, "models", "model"), "run", ACCEL, name, "ws"], cwd=rundir,
                           env=dict(os.environ, NPUSIM_CONFIG_ROOT=cfgroot), capture_output=True, text=True)
        if p.returncode:
            sys.exit(f"{name} {win} failed:\n{(p.stdout + p.stderr)[-1500:]}")
        text = open(os.path.join(rundir, f"{ACCEL}_{name}.txt")).read()
        sect = text[text.index("Multi-chip <-> DRAM transactions"):]
        nin, nw = (float(re.search(rf"\* {t}\s*:\s*(\d+)/", sect).group(1)) * LINE_BYTES for t in ("Input data", "Weight"))
        busy = text[text.index("Busy cycles (ratio"):]
        stages = [float(re.search(rf"\* {st}\s*:\s*([\d.]+)", busy).group(1)) for st in STAGES]
        return (float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", text).group(1)),
                float(re.search(r"Busy scalar-MAC cycles\s*:\s*([\d.]+)", text).group(1)), nin, nw, stages)
    finally:
        shutil.rmtree(tmp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results", nargs="?", default=os.path.join(HERE, "results", "2026-11-04_colab_v5e1_profile"))
    ap.add_argument("-j", "--jobs", type=int, default=16)
    a = ap.parse_args()
    r = json.load(open(os.path.join(a.results, "results.json")))
    meas = [x for x in r["stream_gemm"] + r["odd_gemm"] if x["weights"] == "bf16" and x["cold"]]
    wins = {}
    for w in csv.DictReader(open(os.path.join(a.results, "xla_windows.csv"))):
        if w["case"].endswith("_bf16_cold"):
            kr, kc = (int(v) for v in w["kernel_window"].split("x"))
            wins[int(w["M"]), int(w["K"]), int(w["N"])] = (kc, kr * 8)   # (128-column blocks, K rows) per window
    shapes = sorted({(x["M"], x["K"], x["N"]) for x in meas if (x["M"], x["K"], x["N"]) in wins})

    strip = lambda N: math.ceil(math.ceil(N / DIM) / MXUS)
    jobs = set()
    for M, K, N in shapes:
        for setup in (0, 330):
            jobs |= {(M, K, s * MXUS * DIM, setup, None) for s in (1, 2)}
        kw, nb, _, _ = window_shape(M, K, N, *wins[M, K, N])
        jobs.add((M, kw, nb * DIM, 0, wins[M, K, N]))
    jobs = sorted(jobs, key=lambda j: -j[0] * j[1] * j[2])
    with ThreadPoolExecutor(a.jobs) as ex:
        sim = dict(zip(jobs, ex.map(run, jobs)))

    def setup_ab(M, K, N, setup):
        v1, v2 = sim[M, K, MXUS * DIM, setup, None], sim[M, K, 2 * MXUS * DIM, setup, None]
        s = strip(N)
        return tuple(x1 + (s - 1) * (x2 - x1) for x1, x2 in zip(v1[:4], v2[:4]))

    def setup_d(M, K, N):
        v1, v2 = sim[M, K, MXUS * DIM, 0, None], sim[M, K, 2 * MXUS * DIM, 0, None]
        s = math.ceil(N / DIM) / MXUS
        lin = lambda a, b: a + (s - 1) * (b - a)
        stages = [lin(a, b) for a, b in zip(v1[4], v2[4])]
        _, _, n_n, n_k = window_shape(M, K, N, *wins[M, K, N])
        n = n_n * n_k
        d, c = stages[0] / n, max(stages[1:]) / n
        return (330 + d + (n - 1) * max(d, c) + c,) + tuple(lin(a, b) for a, b in zip(v1[1:4], v2[1:4]))

    def setup_c(M, K, N):
        kw, nb, n_n, n_k = window_shape(M, K, N, *wins[M, K, N])
        _, macs, nin, nw, st = sim[M, kw, nb * DIM, 0, wins[M, K, N]]
        n = n_n * n_k
        d, c = st[0], max(st[1:])
        return 330 + d + (n - 1) * max(d, c) + c, n * macs, n_k * nin, n * nw

    by = {(x["M"], x["K"], x["N"]): x["us"] for x in meas}
    rows, errs = [], {k: {"decode": [], "compute": [], "small": []} for k in "ABCD"}
    print(f"{'M':>4} {'K':>6} {'N':>6} | {'meas':>7} | {'A us':>7} {'B us':>7} {'C us':>7} {'D us':>7} | {'err A':>6} {'B':>6}"
          f" {'C':>6} {'D':>6} | {'MAC/real A':>10} {'D':>5} | {'bytes/XLA A':>11} {'D':>5} | {'HBM GB/s meas':>13} {'D':>4} |"
          f" {'MXU% meas':>9} {'D':>5}")
    for M, K, N in shapes:
        t = by[M, K, N]
        res = {"A": setup_ab(M, K, N, 0), "B": setup_ab(M, K, N, 330), "C": setup_c(M, K, N), "D": setup_d(M, K, N)}
        xbytes = 2 * (K * N + M * K)
        grp = "decode" if M <= 128 and K * N >= 4096 * 4096 else ("compute" if M >= 512 and K >= 4096 else "small")
        row = dict(M=M, K=K, N=N, window_blocks=wins[M, K, N][0], window_k=wins[M, K, N][1], meas_us=t)
        for k, (cyc, macs, nin, nw) in res.items():
            nb = nin + nw
            us = cyc / CLOCK_HZ * 1e6
            row.update({f"{k}_us": round(us, 3), f"{k}_err_pct": round(100 * (us / t - 1), 2),
                        f"{k}_mac_ratio": round(macs / (M * K * N), 4), f"{k}_byte_ratio": round(nb / xbytes, 4),
                        f"{k}_hbm_gbps": round(nb / us / 1e3, 1), f"{k}_mxu_pct": round(100 * 2 * macs / (us * 1e-6) / PEAK_FLOPS, 2)})
            errs[k][grp].append(abs(100 * (us / t - 1)))
        row.update(meas_hbm_gbps=round(xbytes / t / 1e3, 1), meas_mxu_pct=round(100 * 2 * M * K * N / (t * 1e-6) / PEAK_FLOPS, 2))
        rows.append(row)
        print(f"{M:4d} {K:6d} {N:6d} | {t:7.2f} | {row['A_us']:7.2f} {row['B_us']:7.2f} {row['C_us']:7.2f} {row['D_us']:7.2f} |"
              f" {row['A_err_pct']:+6.1f} {row['B_err_pct']:+6.1f} {row['C_err_pct']:+6.1f} {row['D_err_pct']:+6.1f} |"
              f" {row['A_mac_ratio']:10.3f} {row['D_mac_ratio']:5.3f} | {row['A_byte_ratio']:11.3f} {row['D_byte_ratio']:5.3f} |"
              f" {row['meas_hbm_gbps']:13.0f} {row['D_hbm_gbps']:4.0f} | {row['meas_mxu_pct']:9.1f} {row['D_mxu_pct']:5.1f}")
    for k in "ABCD":
        print(f"MAPE {k}: " + ", ".join(f"{g} {statistics.mean(v):.1f}% (n={len(v)})" for g, v in errs[k].items() if v)
              + f", all {statistics.mean([e for v in errs[k].values() for e in v]):.1f}%")
    out = os.path.join(a.results, "compare_windows.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
