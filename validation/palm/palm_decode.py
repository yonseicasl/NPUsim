#!/usr/bin/env python3
"""PaLM decode latency on TPU v4: NPUsim (per-chip GEMMs) + KV-cache read + ICI communication,
against Pope et al., "Efficiently Scaling Transformer Inference", MLSys 2023 (golden_pope2023_decode.csv).

Partitioning follows the paper: feedforward and attention projections in the 2D weight-stationary
layout (weights E_x F_yz, X*YZ = chips, X and YZ minimizing its communication time 2BL(E/X + F/YZ)),
PaLM's parallel attention/feedforward block fused into one input and one output matmul
(W_in x2 for SwiGLU + W_Q + W_K/W_V, and W_out + W_O), multiquery attention sharded over batch.

Per chip and decode step:
  GEMMs   NPUsim critical path on configs/accelerators/tpu_v4[_int8w].cfg (8 MXUs of 128x128, HBM2
          1200 GB/s): fused input [B, E/X] x [E/X, (2F + H*dh)/YZ + 2*dh], fused output
          [B, (F + H*dh)/YZ] x [(F + H*dh)/YZ, E/X], and once per step the unembedding
          [B, E/X] x [E/X, V/YZ].
  KV      (B/chips) sequences x context x K,V x dh x 2 B / HBM bandwidth. Memory-bound: QK and AV do
          ~24 MAC per KV byte against the chip's ~115 MAC per HBM byte.
  ICI     the paper's 2D weight-stationary time 2B(E/X + (F + H*dh)/YZ) x 2 B / 270 GB/s plus the
          all-to-all of Q into and out of the batch sharding (2 x B*H*dh/YZ x 2 B / 270 GB/s).
          Bandwidth terms only: no per-collective latency.
The paper overlaps communication with compute (Looped CollectiveEinsum, "partially or fully"), so
both bounds are reported: serial (GEMM + KV + ICI) and overlapped (max(GEMM + KV, ICI)).
Not modeled: layernorm, softmax and other vector work, int8 dequantization, sampling, host.

Usage: python3 validation/palm/palm_decode.py   (needs models/model built)
"""
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
HERE = os.path.dirname(os.path.abspath(__file__))

CLOCK_HZ = 1.05e9          # TPU v4
HBM_BPS = 1200e9
ICI_BPS = 270e9            # per-chip interconnect bandwidth quoted by the paper
ACT_BYTES = 2              # bf16 activations and KV cache
DIM, MXUS = 128, 8
MOSAIC_DRAM = None         # a MOSAIC DRAM yaml: NPUsim's DRAM streams on Ramulator2 (needs a MOSAIC build)
CONTEXT = 2048 + 32        # mean context while generating 64 tokens after a 2048-token prompt

MODELS = {  # PaLM (Chowdhery et al. 2022; 540B also in Pope et al. Table D.1)
    "palm540b": dict(layers=118, E=18432, F=73728, H=48, dh=256, V=256000),
    "palm62b": dict(layers=64, E=8192, F=32768, H=32, dh=256, V=256000),
}


def split(chips, E, F):
    """(X, YZ) with X * YZ = chips minimizing E/X + F/YZ (ties: the larger X)."""
    return min(((x, chips // x) for x in range(1, chips + 1) if chips % x == 0),
               key=lambda p: (E / p[0] + F / p[1], -p[0]))


def gemm_files(cfgroot, rundir, name, accel, M, K, N):
    """NPUsim net + mapping for C[M x N] = A[M x K] B[K x N] on one TPU v4 chip: output columns split
    over the 8 MXUs, one 128-column weight strip per DRAM iteration, M rows streamed per fold."""
    h = next((c for c in (64, 32, 16, 8, 4, 2) if K % c == 0 and K // c <= 4096), 1)
    with open(os.path.join(cfgroot, "networks", f"{name}.cfg"), "w") as f:
        f.write(f"[net]\nheight={h}\nwidth={K // h}\nchannels=1\nbatch={M}\n\n[data]\ndataset=imagenet\n"
                f"test=datasets/gemm/test.lst\nlabels=datasets/gemm/labels.lst\nweight=weights/{name}.wgh\n"
                f"top=1\n\n[connected]\noutput={N}\nactivation=linear\n\n[softmax]\ngroups=1\n\n[cost]\ntype=l2\n")
    with open(os.path.join(rundir, "weights", f"{name}.wgh"), "wb") as f:
        f.truncate(4 * (N + K * N))                     # sparse zero weights (timing build)
    strips = math.ceil(math.ceil(N / DIM) / MXUS)
    rows = {"MAC": [1, 1], "PE": [1, M], "PE_X": [DIM, 1], "PE_Y": [1, 1], "GLB": [1, 1],
            "CHIPS_X": [1, 1], "CHIPS_Y": [MXUS, 1], "DRAM": [strips, 1]}
    c = {"PE_Y": DIM, "GLB": math.ceil(K / DIM)}
    mdir = os.path.join(cfgroot, "mappings", accel, name)
    os.makedirs(mdir, exist_ok=True)
    with open(os.path.join(mdir, "ws.map"), "w") as f:
        f.write(f"# K, B, P, Q, C, R, S, H, W, GROUP, STRIDE\n\n[connected]\n# {N}, {M}, 1, 1, {K}, 1, 1, 0, 0, 1, 1,\n")
        for lv, (k, b) in rows.items():
            f.write(f"{lv:7s} = {k}, {b}, 1, 1, {c.get(lv, 1)}, 1, 1, 0, 0, 1, 1,\n")


def run_gemm(args):
    cfgroot, rundir, name, accel = args
    env = dict(os.environ, NPUSIM_CONFIG_ROOT=cfgroot)
    env.pop("NPUSIM_MOSAIC_DRAM", None)
    if MOSAIC_DRAM:                                   # time the DRAM tile streams on MOSAIC's Ramulator2
        env["NPUSIM_MOSAIC_DRAM"] = MOSAIC_DRAM
    p = subprocess.run([os.path.join(ROOT, "models", "model"), "run", accel, name, "ws"], cwd=rundir, env=env,
                       capture_output=True, text=True)
    if p.returncode:
        sys.exit(f"{name} failed:\n{(p.stdout + p.stderr)[-1500:]}")
    text = open(os.path.join(rundir, f"{accel}_{name}.txt")).read()
    return float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", text).group(1)) / CLOCK_HZ


def workdir(accels):
    """(tmp, config root, run dir) holding copies of configs/accelerators/<accels>.cfg."""
    tmp = tempfile.mkdtemp(prefix="palm_")
    cfgroot, rundir = os.path.join(tmp, "configs"), os.path.join(tmp, "run")
    for d in ("accelerators", "networks", "mappings"):
        os.makedirs(os.path.join(cfgroot, d))
    os.makedirs(os.path.join(rundir, "weights"))
    os.symlink(os.path.join(ROOT, "models", "datasets"), os.path.join(rundir, "datasets"))
    for a in accels:
        shutil.copy(os.path.join(ROOT, "configs", "accelerators", f"{a}.cfg"), os.path.join(cfgroot, "accelerators"))
    return tmp, cfgroot, rundir


def main():
    golden = list(csv.DictReader(open(os.path.join(HERE, "golden_pope2023_decode.csv"))))
    tmp, cfgroot, rundir = workdir(("tpu_v4", "tpu_v4_int8w"))

    cases, jobs = [], []
    for g in golden:
        m, n, B = MODELS[g["model"]], int(g["chips"]), int(g["batch"])
        accel = "tpu_v4_int8w" if g["weights"] == "int8" else "tpu_v4"
        Fq = m["F"] + m["H"] * m["dh"]                    # YZ-partitioned width after the fused input matmul
        X, YZ = split(n, m["E"], Fq)
        shapes = {"in": (B, m["E"] // X, (2 * m["F"] + m["H"] * m["dh"]) // YZ + 2 * m["dh"]),
                  "out": (B, Fq // YZ, m["E"] // X),
                  "unembed": (B, m["E"] // X, m["V"] // YZ)}
        names = {}
        for part, (M, K, N) in shapes.items():
            # Every MXU repeats the same 128-column weight strip per DRAM iteration and the simulated
            # time is linear in the strip count (checked to 0.2% at 1-3 strips), while NPUsim's wall
            # time grows with M x N: simulate 1 and 2 strips and extrapolate to the real count.
            names[part] = [f"palm_{accel}_{M}x{K}x{s * MXUS * DIM}" for s in (1, 2)]
            for s, name in zip((1, 2), names[part]):
                gemm_files(cfgroot, rundir, name, accel, M, K, s * MXUS * DIM)
                jobs.append((cfgroot, rundir, name, accel))
            names[part].append(math.ceil(math.ceil(N / DIM) / MXUS))
        cases.append((g, m, n, B, X, YZ, Fq, shapes, names))
    jobs = list({j[2]: j for j in jobs}.values())
    with ThreadPoolExecutor(8) as ex:
        sim = dict(zip((j[2] for j in jobs), ex.map(run_gemm, jobs)))
    times = {}
    for case in cases:
        for one, two, strips in case[-1].values():
            times[one, strips] = sim[one] + (strips - 1) * (sim[two] - sim[one])

    print(f"{'case':28} {'X x YZ':7} {'GEMM':>8} {'KV':>7} {'ICI':>7} {'unemb':>6} | {'serial':>7} {'overlap':>7}"
          f" {'roofline':>8} | {'measured':>8}  serial/overlap err")
    for g, m, n, B, X, YZ, Fq, shapes, names in cases:
        L = m["layers"]
        t = {part: times[one, strips] for part, (one, _, strips) in names.items()}
        gemm = t["in"] + t["out"]
        kv = max(B / n, 1) * CONTEXT * 2 * m["dh"] * ACT_BYTES / HBM_BPS
        ici = (2 * B * (m["E"] / X + Fq / YZ) + 2 * B * m["H"] * m["dh"] / YZ) * ACT_BYTES / ICI_BPS
        unemb = t["unembed"]
        serial = L * (gemm + kv + ici) + unemb
        overlap = L * max(gemm + kv, ici) + unemb
        wbytes = 1 if g["weights"] == "int8" else 2
        params = L * (3 * m["E"] * m["F"] + 2 * m["E"] * m["H"] * m["dh"] + 2 * m["E"] * m["dh"]) + m["E"] * m["V"]
        roof = max(params * wbytes / n / HBM_BPS, 2 * params * B / n / (MXUS * DIM * DIM * 2 * CLOCK_HZ)) \
            + L * kv
        meas = float(g["latency_s_64_tokens"]) / 64
        ms = lambda t: f"{t * 1e3:7.2f}"
        print(f"{g['model']} {n:2d}ch B{B:<3d} {g['weights']:5} {X}x{YZ:<5} {ms(L * gemm):>8} {ms(L * kv)} {ms(L * ici)}"
              f" {ms(unemb):>6} | {ms(serial)} {ms(overlap)} {ms(roof):>8} | {ms(meas):>8}  "
              f"{100 * (serial / meas - 1):+6.1f}% / {100 * (overlap / meas - 1):+6.1f}%")
    print("(ms per decoded token; GEMM/KV/ICI summed over layers)")
    shutil.rmtree(tmp)


if __name__ == "__main__":
    main()
