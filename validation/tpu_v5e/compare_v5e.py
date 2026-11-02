#!/usr/bin/env python3
"""NPUsim (configs/accelerators/tpu_v5e.cfg) vs Cloud TPU v5e measurements from measure_tpu.py.

  GEMM    every gemm.csv shape on NPUsim with validation/palm's mapping (output columns over the 4 MXUs
          in 128-column weight strips). Rows go in tiles of <= 512; one tile is simulated at 1 and 2
          strips and scaled, since the time is linear in strips (palm_decode.py) and in row tiles
          (1024 rows = 2 x 512 to 1.5%). Against the measured device time (wall time without a trace).
  decode  decode.csv per-scope times: the four projection GEMMs vs the rest (attention, norms, SwiGLU,
          psum) -- the time NPUsim's PaLM model leaves out -- and, at tp=1, NPUsim's projection GEMMs
          against the measured GEMM scopes.

Usage: python3 validation/tpu_v5e/compare_v5e.py RESULTS_DIR [-j 16]   (needs models/model built)
Writes RESULTS_DIR/compare_gemm.csv.
"""
import argparse
import csv
import math
import os
import shutil
import statistics
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "palm"))
import palm_decode as palm  # noqa: E402

palm.MXUS, palm.CLOCK_HZ = 4, 1.5e9        # TPU v5e; palm.gemm_files / run_gemm read these
ACCEL, DIM, MTILE = "tpu_v5e", 128, 512
PEAK, HBM = 4 * DIM * DIM * 2 * 1.5e9, 819e9
D, HKV, DH = 4096, 8, 128                  # measure_tpu.py's Llama-3-8B decode layer
PROJ = {"qkv": (4096, 6144), "o_proj": (4096, 4096), "ffn_up": (4096, 28672), "ffn_down": (14336, 4096)}


def num(r, k):
    return float(r[k]) if r.get(k) not in (None, "") else None


def simulate(shapes, jobs):
    """{(M, K, N): seconds} on NPUsim."""
    tmp, cfgroot, rundir = palm.workdir((ACCEL,))
    plan, todo = {}, {}
    for M, K, N in shapes:
        mt = math.ceil(M / MTILE)
        mtile = math.ceil(M / mt)
        names = [f"v5e_{mtile}x{K}x{s * palm.MXUS * DIM}" for s in (1, 2)]
        for s, name in zip((1, 2), names):
            if name not in todo:
                palm.gemm_files(cfgroot, rundir, name, ACCEL, mtile, K, s * palm.MXUS * DIM)
                todo[name] = (cfgroot, rundir, name, ACCEL)
        plan[M, K, N] = (mt, names, math.ceil(math.ceil(N / DIM) / palm.MXUS))
    print(f"NPUsim: {len(todo)} runs on {jobs} threads ...", flush=True)
    with ThreadPoolExecutor(jobs) as ex:
        sim = dict(zip(todo, ex.map(palm.run_gemm, todo.values())))
    shutil.rmtree(tmp)
    return {k: mt * (sim[a] + (strips - 1) * (sim[b] - sim[a])) for k, (mt, (a, b), strips) in plan.items()}


def roofline(M, K, N):
    return max(2 * M * K * N / PEAK, 2 * (M * K + K * N + M * N) / HBM)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("-j", "--jobs", type=int, default=16)
    a = ap.parse_args()
    gemm = list(csv.DictReader(open(os.path.join(a.results, "gemm.csv"))))
    dpath = os.path.join(a.results, "decode.csv")
    decode = list(csv.DictReader(open(dpath))) if os.path.exists(dpath) else []
    shapes = {(int(r["M"]), int(r["K"]), int(r["N"])) for r in gemm}
    shapes |= {(int(r["batch"]), *kn) for r in decode if r["tp"] == "1" for kn in PROJ.values()}
    sim = simulate(sorted(shapes), a.jobs)

    out, errs = [], {}
    print(f"\n{'M':>5} {'K':>6} {'N':>6} | {'measured':>9} {'NPUsim':>9} {'roofline':>9} us | {'err':>7}  src")
    for r in gemm:
        M, K, N = int(r["M"]), int(r["K"]), int(r["N"])
        src = "device" if num(r, "device_us") is not None else "wall"
        meas = num(r, "device_us") if src == "device" else num(r, "wall_us")
        s, roof = sim[M, K, N] * 1e6, roofline(M, K, N) * 1e6
        err = 100 * (s / meas - 1)
        errs.setdefault("decode (M<=128)" if M <= 128 else "other", []).append(abs(err))
        out.append(dict(M=M, K=K, N=N, measured_us=meas, src=src, npusim_us=round(s, 2), roofline_us=round(roof, 2),
                        err_pct=round(err, 1)))
        print(f"{M:5d} {K:6d} {N:6d} | {meas:9.2f} {s:9.2f} {roof:9.2f}    | {err:+6.1f}%  {src}")
    allerr = [e for v in errs.values() for e in v]
    print("MAPE " + ", ".join(f"{k} {statistics.mean(v):.1f}% (n={len(v)})" for k, v in errs.items())
          + f", all {statistics.mean(allerr):.1f}%")
    with open(os.path.join(a.results, "compare_gemm.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, list(out[0]))
        w.writeheader()
        w.writerows(out)

    if not decode:
        return
    if all(num(r, "device_us") is None for r in decode):
        print("\ndecode: no device trace (wall time only), per-scope breakdown unavailable")
        return
    print(f"\n{'tp':>2} {'B':>4} {'L':>5} | {'total':>8} {'GEMM':>8} {'attn':>7} {'KVroof':>7} {'psum':>7} {'other':>7}"
          f" {'nonGEMM':>7} | {'NPUsim GEMM':>11} {'err':>7}   (us per layer)")
    for r in decode:
        tp, B, L, tot = int(r["tp"]), int(r["batch"]), int(r["context"]), num(r, "device_us")
        sc = {k[:-3]: float(v) for k, v in r.items() if k.endswith("_us") and k not in ("wall_us", "device_us") and v}
        g = sum(sc.get(s, 0) for s in PROJ)
        ps = sc.get("psum_attn", 0) + sc.get("psum_ffn", 0)
        att = sc.get("attn", 0)
        kvroof = B * L * HKV * DH * 2 * 2 / tp / HBM * 1e6
        line = (f"{tp:2d} {B:4d} {L:5d} | {tot:8.1f} {g:8.1f} {att:7.1f} {kvroof:7.1f} {ps:7.1f} {tot - g - att - ps:7.1f}"
                f" {100 * (1 - g / tot):6.1f}% |")
        if tp == 1:
            sg = sum(sim[(B, *kn)] for kn in PROJ.values()) * 1e6
            line += f" {sg:11.1f} {100 * (sg / g - 1) if g else float('nan'):+6.1f}%"
        print(line)
    print("GEMM = qkv + o_proj + ffn_up + ffn_down scopes; other = norms, SwiGLU, residual, unattributed ops")


if __name__ == "__main__":
    main()
