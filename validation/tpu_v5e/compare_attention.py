#!/usr/bin/env python3
"""NPUsim's attention consumer ([kvcache] attention = 1) vs decode attention measured on TPU v5e-1.

Measured: tpu_v5e1_stream_bench.py's `attention` section (Llama-3-8B: 32 query heads, 8 KV heads, head_dim
128, bf16 KV in (B, kv, L, d); distinct caches each iteration, so every step reads its KV from HBM).
NPUsim: configs/accelerators/tpu_v5e.cfg plus a [kvcache] section with the same geometry, batch B and context L,
appended to a tiny projection layer; the report's "Attention step" (the appended decode-attention step: K and V
supply from the DRAM model, QK / softmax / AV) is compared. The step is KV-read bound on v5e, so its rates are
structural, not fitted: QK/AV at the 4 MXUs' 4 x 128^2 MAC/cycle, softmax at ~5 VPU ops per score over 8 x 128
lanes (0.005 cycles/score). two_pass: the bench's einsum -> softmax -> einsum scores all cached tokens before AV.

Usage: python3 validation/tpu_v5e/compare_attention.py [RESULTS_DIR] [-j 16]   (needs models/model built)
"""
import argparse
import json
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

ROOT, ACCEL, CLOCK_HZ, HBM = palm.ROOT, "tpu_v5e", 1.5e9, 819e9
palm.MXUS = 4
QH, KVH, DH, KV_BYTES = 32, 8, 128, 2
KVCACHE = """
[kvcache]
attention = 1
n_q_heads = {qh}
n_kv_heads = {kvh}
head_dim = {dh}
kv_precision_bytes = {kvb}
context_length = {L}
batch = {B}
attention_algorithm = two_pass
kv_schedule = streaming
kv_tile_bytes = 1048576
attention_macs_per_cycle = 65536
softmax_cycles_per_element = 0.005
profile_reference = validation/tpu_v5e/compare_attention.py (structural rates; KV supply from [dram])
"""


def run(case):
    """Attention step seconds for one (B, L), and the whole report's critical path."""
    B, L = case
    tmp = tempfile.mkdtemp(prefix="attn_")
    try:
        cfgroot, rundir = os.path.join(tmp, "configs"), os.path.join(tmp, "run")
        for d in ("accelerators", "networks", "mappings"):
            os.makedirs(os.path.join(cfgroot, d))
        os.makedirs(os.path.join(rundir, "weights"))
        os.symlink(os.path.join(ROOT, "models", "datasets"), os.path.join(rundir, "datasets"))
        with open(os.path.join(ROOT, "configs", "accelerators", f"{ACCEL}.cfg")) as f:
            cfg = f.read()
        with open(os.path.join(cfgroot, "accelerators", f"{ACCEL}.cfg"), "w") as f:
            f.write(cfg + KVCACHE.format(qh=QH, kvh=KVH, dh=DH, kvb=KV_BYTES, L=L, B=B))
        name = f"attn_b{B}_l{L}"
        palm.gemm_files(cfgroot, rundir, name, ACCEL, B, 128, 512)     # tiny projection the step follows
        env = dict(os.environ, NPUSIM_CONFIG_ROOT=cfgroot)
        env.pop("NPUSIM_MOSAIC_DRAM", None)
        p = subprocess.run([os.path.join(ROOT, "models", "model"), "run", ACCEL, name, "ws"], cwd=rundir, env=env,
                           capture_output=True, text=True)
        if p.returncode:
            sys.exit(f"{name} failed:\n{(p.stdout + p.stderr)[-1500:]}")
        text = open(os.path.join(rundir, f"{ACCEL}_{name}.txt")).read()
        step = float(re.search(r"Attention step\s*:\s*([\d.]+)", text).group(1))
        return step / CLOCK_HZ
    finally:
        shutil.rmtree(tmp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results", nargs="?", default=os.path.join(HERE, "results", "2026-11-03_colab_v5e1_followup3"))
    ap.add_argument("-j", "--jobs", type=int, default=16)
    a = ap.parse_args()
    meas = json.load(open(os.path.join(a.results, "results.json")))["attention"]
    cases = [(x["B"], x["L"]) for x in meas]
    with ThreadPoolExecutor(a.jobs) as ex:
        sim = dict(zip(cases, ex.map(run, cases)))
    print(f"{'B':>4} {'L':>6} {'KV MB':>8} | {'measured':>9} {'NPUsim':>9} {'roofline':>9} us | {'err':>7}")
    errs = {"B >= 8": [], "B = 1": []}
    for x in meas:
        B, L = x["B"], x["L"]
        kv = 2 * B * L * KVH * DH * KV_BYTES
        s = sim[B, L] * 1e6
        err = 100 * (s / x["us"] - 1)
        errs["B >= 8" if B >= 8 else "B = 1"].append(abs(err))
        print(f"{B:4d} {L:6d} {kv / 2 ** 20:8.1f} | {x['us']:9.1f} {s:9.1f} {kv / HBM * 1e6:9.1f}    | {err:+6.1f}%")
    print("MAPE " + ", ".join(f"{k} {sum(v) / len(v):.1f}% (n={len(v)})" for k, v in errs.items() if v))


if __name__ == "__main__":
    main()
