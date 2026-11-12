#!/usr/bin/env python3
"""Single-chip Cloud TPU measurements for NPUsim fidelity (validation/tpu_v5e/README.md). JAX only, no npusim
needed: copy this one file to a Colab (v5e-1 TPU) or Kaggle notebook, or a TRC/Cloud TPU VM, and run it. The
multi-chip experiments are in measure_tpu_v5e8.py, which imports this file for its measurement machinery.

  gemm.csv         C = A @ B in bf16, square / decode / prefill / odd shapes
  vector.csv       rmsnorm, softmax, SwiGLU, residual add
  decode.csv       one Llama-3-8B-shaped decode layer on one chip, time per named scope (qkv, attn, o_proj,
                   ffn_*, norms)
  smallk.csv       (--only smallk, not in the default set) fixed per-kernel time: K sweeps at small M, N,
                   the smallest elementwise kernel, chains of n matmuls in one program; fits in smallk_fit.csv
  loopk.csv        (--only loopk, not in the default set) per-iteration cost of a while loop: n dependent
                   matmuls as a fori_loop over stacked weights vs unrolled; fits in loopk_fit.csv
  env.json         devices, versions

Device time comes from the profiler's perfetto trace (XLA module and op events on TPU:0); wall time
(async-pipelined calls, block at the end) is recorded too and is the fallback when no trace is found.

Usage: python3 measure_tpu_v5e1.py [--out DIR] [--quick] [--keep-traces] [--only gemm,vector,decode,smallk,loopk]
"""
import argparse
import csv
import glob
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import sys
import time

try:
    import jax
except ImportError:  # bare TPU VM: Kaggle/Colab images ship JAX, a fresh Cloud TPU VM may not
    import subprocess
    import sys
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "jax[tpu]"])
    import jax
import jax.numpy as jnp
from jax.sharding import Mesh, NamedSharding, PartitionSpec as P

try:
    from jax import shard_map
except ImportError:
    from jax.experimental.shard_map import shard_map

BF16, F32 = jnp.bfloat16, jnp.float32
ARGS = None
SCOPES = ("norm1", "qkv", "attn", "o_proj", "psum_attn", "norm2", "ffn_up", "ffn_act", "ffn_down", "psum_ffn")


# ---------------------------------------------------------------- timing

def parse_trace(d):
    """(module_us list, {op label: total us}) on the first matching device, or None without a trace."""
    paths = glob.glob(os.path.join(d, "**", "perfetto_trace.json.gz"), recursive=True)
    if not paths:
        return None
    ev = json.load(gzip.open(paths[0]))
    ev = ev.get("traceEvents", []) if isinstance(ev, dict) else ev
    meta = lambda kind: {(e["pid"], e.get("tid")) if kind == "thread_name" else e["pid"]: e["args"]["name"]
                         for e in ev if e.get("ph") == "M" and e.get("name") == kind}
    procs, threads = meta("process_name"), meta("thread_name")
    pids = [p for p, n in procs.items() if re.search(ARGS.device_match, n)]
    if not pids:
        return None
    pid = sorted(pids, key=lambda p: procs[p])[0]
    mods, ops = [], {}
    for e in ev:
        if e.get("ph") != "X" or e.get("pid") != pid:
            continue
        t = threads.get((pid, e.get("tid")), "")
        if "Module" in t:
            mods.append(float(e["dur"]))
        elif "Ops" in t:
            a = json.dumps(e.get("args", {}))
            scope = next((s for s in SCOPES if f"/{s}/" in a or f"/{s}\"" in a), None)
            key = f"{scope}|{e['name']}" if scope else f"-|{e['name']}"
            ops[key] = ops.get(key, 0.0) + float(e["dur"])
    return mods, ops


def bench(f, sets, tag):
    """dict(wall_us, device_us, ops) for one call of jitted f. sets: argument tuples used in turn, so that no
    call reuses the previous call's buffers (XLA keeps an unchanged input on chip across executions: decode
    GEMMs with <= ~118 MB of weights ran 3-4x faster than HBM allows when called on the same buffers)."""
    sets = sets if isinstance(sets, list) else [sets]
    k = iter(range(10 ** 9))
    call = lambda: f(*sets[next(k) % len(sets)])
    out = call()
    jax.block_until_ready(out)
    t = time.perf_counter()
    jax.block_until_ready(call())
    n = max(5, min(ARGS.max_iters, int(0.3 / max(time.perf_counter() - t, 1e-6))))
    for _ in range(3):
        out = call()
    jax.block_until_ready(out)
    t = time.perf_counter()
    for _ in range(n):
        out = call()
    jax.block_until_ready(out)
    res = {"wall_us": (time.perf_counter() - t) / n * 1e6, "device_us": None, "ops": None, "iters": n,
           "sets": len(sets)}
    d = os.path.join(ARGS.out, "traces", tag)
    try:
        with jax.profiler.trace(d, create_perfetto_trace=True):
            for _ in range(n):
                out = call()
            jax.block_until_ready(out)
        parsed = parse_trace(d)
        if parsed and parsed[0]:
            mods, ops = parsed
            res["device_us"] = statistics.median(mods[-n:])
            res["ops"] = {k: v / len(mods) for k, v in sorted(ops.items(), key=lambda kv: -kv[1])}
    except Exception as e:  # profiler unavailable: wall time only
        print(f"  [{tag}] profiler: {type(e).__name__}: {e}")
    if not ARGS.keep_traces:
        shutil.rmtree(d, ignore_errors=True)
    return res


def write(name, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(ARGS.out, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {name} ({len(rows)} rows)")


def rnd(key, shape, sharding=None):
    x = jax.random.normal(jax.random.PRNGKey(key), shape, F32).astype(BF16)
    return jax.device_put(x, sharding) if sharding else x


def nsets(set_bytes):
    """Distinct argument sets to rotate: >= --rotate-mb in total, 2..64 (1 with --rotate-mb 0: warm reuse)."""
    if ARGS.rotate_mb <= 0:
        return 1
    return max(2, min(64, math.ceil(ARGS.rotate_mb * 2 ** 20 / set_bytes)))


def fmt(r):
    dev = f"{r['device_us']:9.1f}" if r["device_us"] is not None else "        -"
    return f"wall {r['wall_us']:9.1f} us  device {dev} us"


# ---------------------------------------------------------------- experiments

def gemm():
    sq = [128, 256, 512, 1024, 2048, 4096, 8192]
    proj = [(4096, 6144), (4096, 4096), (4096, 28672), (14336, 4096)]   # Llama-3-8B qkv, o, gate|up, down
    shapes = [(s, s, s) for s in sq] + [(m, k, n) for m in (1, 8, 32, 128, 512, 2048) for k, n in proj] \
        + [(100, 1000, 3000), (777, 1234, 555), (1000, 1000, 1000), (384, 768, 3072)]
    if ARGS.quick:
        shapes = [(256, 256, 256), (8, 4096, 4096), (777, 1234, 555)]
    dev = jax.devices()[0]
    f = jax.jit(lambda a, b: a @ b)
    rows = []
    for M, K, N in shapes:
        sets = [(jax.device_put(rnd(2 * i + 1, (M, K)), dev), jax.device_put(rnd(2 * i + 2, (K, N)), dev))
                for i in range(nsets(2 * (M * K + K * N)))]
        r = bench(f, sets, f"gemm_{M}x{K}x{N}")
        t = r["device_us"] or r["wall_us"]
        rows.append(dict(M=M, K=K, N=N, dtype="bf16", sets=r["sets"], wall_us=round(r["wall_us"], 2),
                         device_us=r["device_us"] and round(r["device_us"], 2),
                         tflops=round(2 * M * K * N / t / 1e6, 2),
                         ops=r["ops"] and ";".join(f"{k}={v:.2f}" for k, v in list(r["ops"].items())[:4])))
        print(f"gemm {M:5d} x {K:5d} x {N:5d}  {fmt(r)}  {rows[-1]['tflops']:7.2f} TFLOPS")
    write("gemm.csv", rows)


def rmsnorm(x, g):
    x32 = x.astype(F32)
    return (x32 * jax.lax.rsqrt(jnp.mean(x32 * x32, -1, keepdims=True) + 1e-5) * g.astype(F32)).astype(BF16)


def vector():
    fns = {"rmsnorm": lambda x, y: rmsnorm(x, y[0]),
           "softmax": lambda x, y: jax.nn.softmax(x.astype(F32), -1).astype(BF16),
           "swiglu": lambda x, y: jax.nn.silu(x) * y,
           "add": lambda x, y: x + y}
    rows = []
    for M in ((8,) if ARGS.quick else (8, 128, 2048)):
        for D in ((4096,) if ARGS.quick else (4096, 14336)):
            sets = [(rnd(2 * i + 3, (M, D)), rnd(2 * i + 4, (M, D))) for i in range(nsets(4 * M * D))]
            for name, fn in fns.items():
                r = bench(jax.jit(fn), sets, f"vec_{name}_{M}x{D}")
                rows.append(dict(op=name, M=M, D=D, wall_us=round(r["wall_us"], 2),
                                 device_us=r["device_us"] and round(r["device_us"], 2)))
                print(f"vector {name:8s} {M:5d} x {D:5d}  {fmt(r)}")
    write("vector.csv", rows)


def decode(tps=(1,)):
    """The decode layer at each tensor-parallel degree in tps (measure_tpu_v5e8.py passes 1, 2, 4, 8; Megatron
    split: column-split qkv / gate-up, row-split o / down + psum)."""
    D, H, HKV, DH, F = 4096, 32, 8, 128, 14336            # Llama-3-8B layer
    devs = jax.devices()
    tps = [t for t in tps if t <= len(devs) and HKV % t == 0]
    grid = [(1, 2048), (32, 2048)] if ARGS.quick else [(b, l) for b in (1, 8, 32, 128) for l in (512, 2048, 8192)]
    rows, ops_all = [], {}
    for tp in tps:
        mesh = Mesh(devs[:tp], ("x",))
        h, hkv = H // tp, HKV // tp

        def layer(x, g1, wqkv, wo, g2, wgu, wd, kc, vc):
            B = x.shape[0]
            with jax.named_scope("norm1"):
                y = rmsnorm(x, g1)
            with jax.named_scope("qkv"):
                qkv = y @ wqkv
            with jax.named_scope("attn"):                 # KV-cache append not modeled (small)
                q = qkv[:, :h * DH].reshape(B, hkv, h // hkv, DH)
                s = jnp.einsum("bgqd,blgd->bgql", q, kc, preferred_element_type=F32) / math.sqrt(DH)
                p = jax.nn.softmax(s, -1).astype(BF16)
                o = jnp.einsum("bgql,blgd->bgqd", p, vc).reshape(B, h * DH)
            with jax.named_scope("o_proj"):
                o = o @ wo
            if tp > 1:
                with jax.named_scope("psum_attn"):
                    o = jax.lax.psum(o, "x")
            x = x + o
            with jax.named_scope("norm2"):
                y = rmsnorm(x, g2)
            with jax.named_scope("ffn_up"):
                gu = y @ wgu
            with jax.named_scope("ffn_act"):
                a = jax.nn.silu(gu[:, :F // tp]) * gu[:, F // tp:]
            with jax.named_scope("ffn_down"):
                o = a @ wd
            if tp > 1:
                with jax.named_scope("psum_ffn"):
                    o = jax.lax.psum(o, "x")
            return x + o

        # per-chip shard layout: wqkv columns [q heads | k head(s) | v head(s)] of this chip,
        # wgu columns [gate | up] of this chip's F slice, kc/vc this chip's KV heads.
        spec = dict(x=P(), g1=P(), wqkv=P(None, "x"), wo=P("x", None), g2=P(), wgu=P(None, "x"), wd=P("x", None),
                    kc=P(None, None, "x", None), vc=P(None, None, "x", None))
        f = jax.jit(shard_map(layer, mesh=mesh, in_specs=tuple(spec.values()), out_specs=P()) if tp > 1 else layer)
        put = lambda k, a, name: rnd(k, a, NamedSharding(mesh, spec[name]))
        # consecutive calls use different layers' weights and KV cache, as a real decode step does
        nl = nsets(2 * (D * (H + 2 * HKV) * DH + H * DH * D + 3 * D * F) // tp)
        ws = [[put(20 * i + 11, (1, D), "g1"), put(20 * i + 12, (D, tp * (h + 2 * hkv) * DH), "wqkv"),
               put(20 * i + 13, (H * DH, D), "wo"), put(20 * i + 14, (1, D), "g2"),
               put(20 * i + 15, (D, 2 * F), "wgu"), put(20 * i + 16, (F, D), "wd")] for i in range(nl)]
        for B, L in grid:
            kv = [(put(20 * i + 17, (B, L, HKV, DH), "kc"), put(20 * i + 18, (B, L, HKV, DH), "vc")) for i in range(nl)]
            x = put(19, (B, D), "x")
            tag = f"decode_tp{tp}_b{B}_l{L}"
            r = bench(f, [(x, *ws[i], *kv[i]) for i in range(nl)], tag)
            scoped = {}
            for k, v in (r["ops"] or {}).items():
                scoped[k.split("|")[0]] = scoped.get(k.split("|")[0], 0.0) + v
            rows.append(dict(tp=tp, batch=B, context=L, sets=r["sets"], wall_us=round(r["wall_us"], 2),
                             device_us=r["device_us"] and round(r["device_us"], 2),
                             **{f"{s}_us": round(scoped[s], 2) for s in (*SCOPES, "-") if s in scoped}))
            ops_all[tag] = r["ops"]
            print(f"decode tp{tp} B{B:<4d} L{L:<5d} {fmt(r)}  " + " ".join(f"{s}={v:.1f}" for s, v in scoped.items()))
            del kv
    write("decode.csv", rows)
    json.dump(ops_all, open(os.path.join(ARGS.out, "decode_ops.json"), "w"), indent=1)


def linfit(xs, ys):
    """(intercept, slope) of the least-squares line."""
    mx, my = statistics.mean(xs), statistics.mean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx if sxx else 0.0
    return my - b * mx, b


def smallk():
    """Fixed per-kernel time (not in the default --only): what a kernel costs beyond its bytes and FLOPs.

    ksweep  C[M x N] = A[M x K] @ B[K x N] at fixed small M, N and K 128..16384: time = fixed + per-byte x
            bytes, so the K = 0 intercept is the fixed cost. device_us is the XLA module (program launch and
            finish included), ops_us the sum of its op events (the kernels alone).
    floor   x + 1 on one vreg ([8, 128]) up to 4 MB: the smallest kernel.
    chain   n dependent matmuls y = y @ W_i (W_i distinct, K x K) in ONE program, n = 1..32: the slope over n
            is what one more kernel costs inside a program (as in a model's layer), the n = 1 intercept the
            program's own overhead.
    oddchain the odd small GEMMs (M 100-512, K = N 1000-2048) as chains of 1 and 8 in one program: the per-kernel
            slope is their in-program cost (no one-kernel program fixed part, no operand copies).
    copy_us is the copy-start/done time inside each module (operands staged into on-chip memory).
    Writes smallk.csv (all points) and smallk_fit.csv (intercepts and slopes).
    """
    dev = jax.devices()[0]
    rows, fits = [], []
    put = lambda k, s, scale=1.0: jax.device_put((rnd(k, s).astype(F32) * scale).astype(BF16), dev)
    ops_sum = lambda r: r["ops"] and round(sum(r["ops"].values()), 3)
    copy_sum = lambda r: r["ops"] and round(sum(v for k, v in r["ops"].items() if "copy-" in k), 3)

    ks = (128, 512, 2048) if ARGS.quick else (128, 256, 512, 1024, 2048, 4096, 8192, 16384)
    mn = [(8, 128)] if ARGS.quick else [(8, 128), (8, 1024), (128, 128), (128, 1024)]
    f = jax.jit(lambda a, b: a @ b)
    for M, N in mn:
        pts = []
        for K in ks:
            sets = [(put(2 * i + 1, (M, K)), put(2 * i + 2, (K, N))) for i in range(nsets(2 * (M * K + K * N)))]
            r = bench(f, sets, f"smallk_{M}x{K}x{N}")
            rows.append(dict(kind="ksweep", M=M, K=K, N=N, n=1, bytes=2 * (M * K + K * N + M * N), sets=r["sets"],
                             wall_us=round(r["wall_us"], 3), device_us=r["device_us"] and round(r["device_us"], 3),
                             ops_us=ops_sum(r), copy_us=copy_sum(r)))
            pts.append(rows[-1])
            print(f"smallk ksweep {M:4d} x {K:5d} x {N:5d}  {fmt(r)}  ops {rows[-1]['ops_us']}")
        for col in ("device_us", "ops_us"):
            xy = [(p["bytes"], p[col]) for p in pts if p[col] is not None]
            if len(xy) >= 2:
                a, b = linfit(*zip(*xy))
                fits.append(dict(kind="ksweep", M=M, N=N, time=col, fixed_us=round(a, 3),
                                 per_byte_GBps=round(1e-3 / b, 1) if b > 0 else None))
                print(f"  fit {col}: fixed {a:.3f} us, {1e-3 / b if b > 0 else float('nan'):.1f} GB/s")

    shapes = [(8, 128), (8, 1024)] if ARGS.quick else [(8, 128), (8, 1024), (128, 1024), (1024, 1024)]
    g = jax.jit(lambda x: x + 1)
    for M, D in shapes:
        sets = [(put(i + 100, (M, D)),) for i in range(nsets(2 * M * D))]
        r = bench(g, sets, f"smallk_floor_{M}x{D}")
        rows.append(dict(kind="floor", M=M, K=D, N=0, n=1, bytes=4 * M * D, sets=r["sets"],
                         wall_us=round(r["wall_us"], 3), device_us=r["device_us"] and round(r["device_us"], 3),
                         ops_us=ops_sum(r), copy_us=copy_sum(r)))
        print(f"smallk floor  {M:4d} x {D:5d}  {fmt(r)}  ops {rows[-1]['ops_us']}")

    ns = (1, 2, 4, 8) if ARGS.quick else (1, 2, 4, 8, 16, 32)
    for M, K in ([(8, 128)] if ARGS.quick else [(8, 128), (8, 1024), (128, 1024)]):
        pts = []
        for n in ns:
            def h(x, *ws):
                for w in ws:
                    x = x @ w
                return x
            # weights scaled by 1/sqrt(K) keep the chained values in range
            sets = [(put(1000 * i, (M, K)), *[put(1000 * i + j + 1, (K, K), K ** -0.5) for j in range(n)])
                    for i in range(nsets(2 * (M * K + n * K * K)))]
            r = bench(jax.jit(h), sets, f"smallk_chain_{M}x{K}_n{n}")
            rows.append(dict(kind="chain", M=M, K=K, N=K, n=n, bytes=2 * (M * K + n * K * K + M * K), sets=r["sets"],
                             wall_us=round(r["wall_us"], 3), device_us=r["device_us"] and round(r["device_us"], 3),
                             ops_us=ops_sum(r), copy_us=copy_sum(r)))
            pts.append(rows[-1])
            print(f"smallk chain  {M:4d} x {K:5d} n={n:2d}  {fmt(r)}  ops {rows[-1]['ops_us']}")
        xy = [(p["n"], p["device_us"]) for p in pts if p["device_us"] is not None]
        if len(xy) >= 2:
            a, b = linfit(*zip(*xy))
            fits.append(dict(kind="chain", M=M, N=K, time="device_us", fixed_us=round(a, 3), per_kernel_us=round(b, 3)))
            print(f"  fit chain: program {a:.3f} us + {b:.3f} us per kernel")

    # The odd small GEMMs (K = N around multiples of 128) inside one program: 1 and 8 dependent y = y @ W_i; the
    # per-kernel slope (t8 - t1) / 7 is their cost without the one-kernel program's fixed part and copies.
    odd = [(128, 1024), (128, 1100)] if ARGS.quick else \
        [(m, k) for m in (100, 128, 200, 256, 512) for k in (1000, 1024, 1100, 1152, 1280, 2000, 2048)]
    for M, K in odd:
        t = {}
        for n in (1, 8):
            def h(x, *ws):
                for w in ws:
                    x = x @ w
                return x
            sets = [(put(1000 * i + 7, (M, K)), *[put(1000 * i + j + 8, (K, K), K ** -0.5) for j in range(n)])
                    for i in range(nsets(2 * (M * K + n * K * K)))]
            r = bench(jax.jit(h), sets, f"smallk_oddchain_{M}x{K}_n{n}")
            rows.append(dict(kind="oddchain", M=M, K=K, N=K, n=n, bytes=2 * (M * K + n * K * K + M * K),
                             sets=r["sets"], wall_us=round(r["wall_us"], 3),
                             device_us=r["device_us"] and round(r["device_us"], 3), ops_us=ops_sum(r),
                             copy_us=copy_sum(r)))
            t[n] = r["device_us"]
            print(f"smallk oddchain {M:4d} x {K:5d} n={n}  {fmt(r)}  ops {rows[-1]['ops_us']} copy {rows[-1]['copy_us']}")
        if t[1] is not None and t[8] is not None:
            fits.append(dict(kind="oddchain", M=M, N=K, time="device_us", fixed_us=round(t[1] - (t[8] - t[1]) / 7, 3),
                             per_kernel_us=round((t[8] - t[1]) / 7, 3)))
    write("smallk.csv", rows)
    write("smallk_fit.csv", fits)


def loopk():
    """Per-iteration cost of a while loop (not in the default --only): n dependent y = y @ W_i (distinct K x K
    weights) once as a jax.lax.fori_loop over a stacked [n, K, K] array (one dynamic slice per iteration, the
    structure of a scanned layer stack) and once unrolled in the same program, n = 1, 4, 16. The slope over n
    of each is its per-kernel cost; loop - unrolled is what a loop iteration adds. Writes loopk.csv and
    loopk_fit.csv (copy_us: the module's copy-start/done time)."""
    dev = jax.devices()[0]
    put = lambda k, s, scale=1.0: jax.device_put((rnd(k, s).astype(F32) * scale).astype(BF16), dev)
    ops_sum = lambda r: r["ops"] and round(sum(r["ops"].values()), 3)
    copy_sum = lambda r: r["ops"] and round(sum(v for k, v in r["ops"].items() if "copy-" in k), 3)
    shapes = [(8, 128), (128, 1024)] if ARGS.quick else [(8, 128), (100, 1000), (128, 1024), (256, 2048), (8, 4096)]
    ns = (1, 4) if ARGS.quick else (1, 4, 16)
    rows, fits = [], []

    def looped(x, ws):
        return jax.lax.fori_loop(0, ws.shape[0], lambda i, y: y @ ws[i], x)

    def unrolled(x, ws):
        for i in range(ws.shape[0]):
            x = x @ ws[i]
        return x

    for M, K in shapes:
        pts = {"loop": [], "unrolled": []}
        for n in ns:
            # weights scaled by 1/sqrt(K) keep the chained values in range; distinct sets rotate as elsewhere
            sets = [(put(5000 * i + 1, (M, K)), put(5000 * i + 2, (n, K, K), K ** -0.5))
                    for i in range(nsets(2 * (M * K + n * K * K)))]
            for kind, fn in (("loop", looped), ("unrolled", unrolled)):
                r = bench(jax.jit(fn), sets, f"loopk_{kind}_{M}x{K}_n{n}")
                rows.append(dict(kind=kind, M=M, K=K, N=K, n=n, bytes=2 * (M * K + n * K * K + M * K), sets=r["sets"],
                                 wall_us=round(r["wall_us"], 3), device_us=r["device_us"] and round(r["device_us"], 3),
                                 ops_us=ops_sum(r), copy_us=copy_sum(r)))
                pts[kind].append(rows[-1])
                print(f"loopk {kind:8s} {M:4d} x {K:5d} n={n:2d}  {fmt(r)}  copy {rows[-1]['copy_us']}")
        slope = {}
        for kind, p in pts.items():
            xy = [(q["n"], q["device_us"]) for q in p if q["device_us"] is not None]
            if len(xy) >= 2:
                a, b = linfit(*zip(*xy))
                slope[kind] = b
                fits.append(dict(kind=kind, M=M, N=K, fixed_us=round(a, 3), per_kernel_us=round(b, 3)))
        if len(slope) == 2:
            fits.append(dict(kind="loop_minus_unrolled", M=M, N=K, fixed_us="",
                             per_kernel_us=round(slope["loop"] - slope["unrolled"], 3)))
            print(f"  {M} x {K}: loop {slope['loop']:.3f} us, unrolled {slope['unrolled']:.3f} us per kernel "
                  f"-> {slope['loop'] - slope['unrolled']:+.3f} us per loop iteration")
    write("loopk.csv", rows)
    write("loopk_fit.csv", fits)


def balance():
    """bf16 vs int8 weights where HBM streaming and MXU compute balance (not in the default --only).

    Stored data: with int8 weights at M = 128 (compute ~ stream) the GEMM fusion of the stream bench takes
    1.2-1.36 x max(stream, compute); bf16 reaches compute ~ stream at M ~ 256, never measured with more than 3
    windows. Both dtypes run here on the SAME schedule as that data, tpu_v5e1_stream_bench.py's gemm_loop: one
    fori_loop iteration is (x + carry) @ W[i % R].astype(bf16), a float32 sum of the result into the carry, W a
    stack of R distinct K x K copies (>= --rotate-mb in total). Per (dtype, K, M), three scopes, never mixed:
      loop_slope_us  (median wall(2n) - median wall(n)) / n: one loop iteration, wall clock, as in the stored csv
      fused_ops_sum_per_iter_us  profiled device time of ALL fusion / convolution ops of the 8-iteration program,
                     / 8. NOT the GEMM alone: the loop's input add, output sum and scalar fusions are in it
                     unless the HLO shows them fused into the convolution. module_us_n8 is the whole program;
                     ops is its complete op breakdown, and the traces are kept (traces/), so each op can be
                     attributed with the saved HLO afterwards
      schedule=isolated rows: the same GEMM as a one-op program (operands staged first, a different schedule) --
                     a control only; compare loop rows with loop rows.
    Reading, on loop rows (indications to follow up, none a proof by itself -- a conversion pipelined under the
    compute would not show as a ratio): bf16 at M ~ 256 also ~1.2 x max -> not int8, a pipeline effect of any dtype; only int8
    -> a dtype-specific stage; excess / windows constant over K -> per window; excess / weight elements constant
    -> per element. Writes balance.csv and the optimized HLO (window_config) of every program under hlo/."""
    if jax.default_backend() != "tpu" and not ARGS.quick:
        raise SystemExit("balance: no TPU backend (full sizes are TPU-only; --quick is the CPU smoke test)")
    ARGS.keep_traces = True                       # needed to attribute ops to HLO instructions later
    ms = (32, 128) if ARGS.quick else (32, 64, 128, 192, 256, 384, 512)
    ks = (512,) if ARGS.quick else (4096, 8192)

    def loop(n, x, ws):
        def body(i, c):
            w = jax.lax.dynamic_index_in_dim(ws, i % ws.shape[0], keepdims=False)
            y = (x + c.astype(x.dtype)) @ w.astype(BF16)
            return c + jnp.sum(y.astype(F32)) * 1e-9
        return jax.lax.fori_loop(0, n, body, jnp.float32(0))
    run_loop = jax.jit(loop, static_argnums=0)
    isolated = jax.jit(lambda a, b: a @ b.astype(BF16))

    def stack(key, shape, dtype):                 # generated on the device, as the stream bench does
        if dtype == "int8":
            gen = lambda k: jax.random.randint(k, shape, -8, 8, jnp.int8)
        else:
            gen = lambda k: (jax.random.normal(k, shape, F32) * 0.02).astype(BF16)
        return jax.jit(gen)(jax.random.PRNGKey(key))

    def hlo(tag, f, *args):
        try:
            open(os.path.join(ARGS.out, "hlo", tag + ".txt"), "w").write(f.lower(*args).compile().as_text())
        except Exception as e:
            print(f"  [{tag}] hlo: {type(e).__name__}: {e}")

    def slope(x, ws, n, repeats=5):
        def wall(k):
            t = time.perf_counter()
            jax.block_until_ready(run_loop(k, x, ws))
            return time.perf_counter() - t
        for k in (n, 2 * n):
            wall(k)
        t1 = statistics.median(wall(n) for _ in range(repeats))
        t2 = statistics.median(wall(2 * n) for _ in range(repeats))
        return max(t2 - t1, 0.0) / n * 1e6

    rows = []
    os.makedirs(os.path.join(ARGS.out, "hlo"), exist_ok=True)
    opstr = lambda r: r["ops"] and ";".join(f"{k}={v:.3f}" for k, v in r["ops"].items())   # every op
    for dtype, wbytes in (("bf16", 2), ("int8", 1)):
        for K in ks:
            R = 1 if ARGS.rotate_mb <= 0 else max(2, math.ceil(ARGS.rotate_mb * 2 ** 20 / (K * K * wbytes)))
            ws = stack(K + wbytes, (R, K, K), dtype)
            for M in ms:
                x = stack(M + 7, (M, K), "bf16")
                tag = f"balance_{dtype}_{M}x{K}"
                base = dict(dtype=dtype, M=M, K=K, N=K, copies=R)
                hlo(tag + "_loop", run_loop, 8, x, ws)
                us = slope(x, ws, 10 if M * K * K > 2 ** 33 else 40)
                r = bench(lambda a, b: run_loop(8, a, b), (x, ws), tag + "_loop")
                fus = r["ops"] and sum(v for k, v in r["ops"].items() if "fusion" in k or "convolution" in k) / 8
                rows.append(dict(base, schedule="loop", loop_slope_us=round(us, 3),
                                 module_us_n8=r["device_us"] and round(r["device_us"], 3),
                                 fused_ops_sum_per_iter_us=fus and round(fus, 3), ops=opstr(r)))
                print(f"balance loop     {dtype} {M:4d} x {K:5d}  slope {us:9.2f} us  fused ops / iter {fus and round(fus, 2)}")
                hlo(tag + "_isolated", isolated, x, ws[0])
                r = bench(isolated, [(x, ws[i]) for i in range(R)], tag + "_isolated")
                rows.append(dict(base, schedule="isolated", wall_us=round(r["wall_us"], 3),
                                 device_us=r["device_us"] and round(r["device_us"], 3),
                                 copy_us=r["ops"] and round(sum(v for k, v in r["ops"].items() if "copy-" in k), 3),
                                 ops=opstr(r)))
                print(f"balance isolated {dtype} {M:4d} x {K:5d}  {fmt(r)}")
            del ws
    write("balance.csv", rows)


def hlo_facts(text):
    """The compiler's choices in an optimized HLO text, as written -- nothing inferred about hardware.
    windows: the window_config fields of each convolution; convolutions: (result type, operand types);
    tagged: non-scalar tensors carrying a memory-space tag S(n) (untagged = the default space); copies and
    cross-program prefetch; fusion / while counts."""
    types = dict(re.findall(r"(?m)^\s*(?:ROOT\s+)?%([\w.\-]+) = ([a-z]\w*\[[\d,]*\](?:\{[^}]*\})?) ", text))
    convs = []
    for res, args in re.findall(r"= ([a-z]\w*\[[\d,]*\](?:\{[^}]*\})?) convolution\(([^)]*)\)", text):
        ops = []
        for arg in args.split(","):
            arg = arg.strip()
            m = re.match(r"([a-z]\w*\[[\d,]*\])", arg)                # operand printed with its type
            ops.append(m.group(1) if m else re.sub(r"\{.*", "", types.get(arg.lstrip("%"), "?")))
        convs.append((re.sub(r"\{.*", "", res), tuple(ops)))
    keys = ("iteration_bounds", "kernel_window_bounds", "input_window_bounds", "output_window_bounds")
    found = {k: re.findall(rf'"{k}":\[([^\]]*)\]', text) for k in keys}
    n = len(found["iteration_bounds"])
    windows = [{k: found[k][i].replace('"', "") if i < len(found[k]) else None for k in keys} for i in range(n)]
    tagged = sorted({f"{shape}:S({space})" for shape, space in
                     re.findall(r"([a-z]\w*\[\d+,[\d,]+\])\{[^}]*S\((\d+)\)[^}]*\}", text)})
    return dict(convolutions=convs, windows=windows, tagged=tagged,
                copy_start=len(re.findall(r" copy-start\(", text)), copy_done=len(re.findall(r" copy-done\(", text)),
                prefetch="cross_program_prefetch_index" in text,
                fusions=len(re.findall(r" fusion\(", text)), whiles=len(re.findall(r" while\(", text)))


def op_events(d, per_module):
    """Raw device events of a kept trace: module count and, per op name, (count, total us); the events whose
    hlo_category says convolution; sha256 of the trace files. `verified` only when exactly `per_module`
    convolution-fusion events were recorded per module -- the trace-level evidence that the GEMM ran once per
    requested loop iteration. None without a readable device trace."""
    files = sorted(glob.glob(os.path.join(d, "**", "*"), recursive=True))
    hashes = {os.path.basename(f): hashlib.sha256(open(f, "rb").read()).hexdigest() for f in files if os.path.isfile(f)}
    paths = [f for f in files if f.endswith("perfetto_trace.json.gz")]
    if not paths:
        return None
    ev = json.load(gzip.open(paths[0]))
    ev = ev.get("traceEvents", []) if isinstance(ev, dict) else ev
    name = lambda kind, key: {key(e): e["args"]["name"] for e in ev if e.get("ph") == "M" and e.get("name") == kind}
    procs, threads = name("process_name", lambda e: e["pid"]), name("thread_name", lambda e: (e["pid"], e.get("tid")))
    pids = sorted((p for p, n in procs.items() if re.search(ARGS.device_match, n)), key=lambda p: procs[p])
    if not pids:
        return None
    modules, ops, conv = 0, {}, [0, 0.0]
    for e in ev:
        if e.get("ph") != "X" or e.get("pid") != pids[0]:
            continue
        t = threads.get((pids[0], e.get("tid")), "")
        if "Module" in t:
            modules += 1
        elif "Ops" in t:
            c = ops.setdefault(e["name"], [0, 0.0])
            c[0], c[1] = c[0] + 1, c[1] + float(e["dur"])
            if "convolution" in str(e.get("args", {}).get("hlo_category", "")):
                conv[0], conv[1] = conv[0] + 1, conv[1] + float(e["dur"])
    return dict(modules=modules, ops=ops, conv_events=conv[0], conv_us=conv[1], files=hashes,
                verified=modules > 0 and conv[0] == per_module * modules)


def pair_verdict(banked, single, K, banked_ev=None, single_ev=None):
    """What the two programs' artifacts allow to be said -- (verdict, notes). Deliberately conservative:
      * differing HLO signatures (convolution types, window geometry, fusion / loop counts) -> inconclusive;
      * equal signatures are only "the same observed signatures": global counts do not show that the GEMM is
        inside the loop. That needs the trace: exactly one convolution-fusion event per requested iteration in
        BOTH programs (op_events().verified); without it -> inconclusive;
      * weight staging is never concluded. A memory-space tag on a weight-shaped tensor, or extra copy /
        prefetch ops, in the single-weight program are listed as RAW evidence (a copy may be of the input, a tag
        is not a traced dataflow), and their absence is "not established", not "not staged"."""
    diff = [f"{k}: banked {banked[k]} vs single {single[k]}"
            for k in ("convolutions", "windows", "fusions", "whiles") if banked[k] != single[k]]
    if not banked["convolutions"] or not banked["windows"]:
        diff.append("no convolution / window_config in the HLO (not a TPU compile?)")
    if diff:
        return "inconclusive: the two HLO texts do not show the same signatures", diff
    if not (banked_ev and single_ev and banked_ev["verified"] and single_ev["verified"]):
        got = lambda e: "no device trace" if not e else f"{e['conv_events']} convolution-fusion events in {e['modules']} modules"
        return ("inconclusive: same observed HLO signatures, but one convolution-fusion event per loop iteration "
                "is not verified in both traces"), [f"banked: {got(banked_ev)}", f"single: {got(single_ev)}"]
    weightish = lambda f: {t for t in f["tagged"] if re.search(rf"\[(?:\d+,)?{K},{K}\]", t)}
    raw = [f"memory-space tag on a weight-shaped tensor only in the single program: {sorted(weightish(single) - weightish(banked))}"] \
        if weightish(single) - weightish(banked) else []
    if single["copy_start"] != banked["copy_start"] or single["prefetch"] != banked["prefetch"]:
        raw.append(f"copy-start {banked['copy_start']} -> {single['copy_start']}, prefetch {banked['prefetch']} -> "
                   f"{single['prefetch']} (which operand is copied is not identified)")
    return ("same observed HLO signatures, one convolution-fusion event per iteration in both traces; weight "
            "staging NOT ESTABLISHED from these artifacts" + (" -- raw evidence listed, not a conclusion" if raw else "")
            + "; any time difference is a two-program contrast, not a measured DMA tail"), raw


def tailprobe():
    """Optional probe for the memory-bound bf16 over-prediction of the window-schedule candidate (not in the
    default --only). The four regressing cases (bf16 32 / 64 x 4096 / 8192) and one int8 control, each as TWO
    programs with the SAME loop body as `balance` ((x + carry) @ W[i % R].astype(bf16), f32 sum into the carry):
      banked  W = a bank of R distinct weights (>= --rotate-mb): every iteration streams a new weight
      single  W = one weight, reused by every iteration
    Per program (tailprobe.csv), each number with its scope:
      wall_n_s / wall_2n_s, loop_slope_us   raw host wall samples of the n- and 2n-iteration loop and the SIGNED
                                            slope (median(2n) - median(n)) / n; slope_flag marks a nonpositive one
      modules_traced, op_event_counts,      raw profiler events of the 8-iteration program
      op_event_us
      conv_fusion_events_per_module,        the events whose hlo_category is a convolution; a per-event time is
      conv_fusion_us_per_event              given only when there are exactly 8 per module (else blank)
      all_fusion_event_sum_per_requested_iter_us   every fusion / convolution event summed, / 8 -- may contain
                                            fusions outside the loop; NOT a loop-kernel time
      hlo_facts, hlo_sha256, trace_sha256   compiler output (static) and file hashes
    Per case (tailprobe_pairs.csv): pair_verdict() -- inconclusive unless the HLO signatures match AND both
    traces show one convolution-fusion event per iteration; staging is never concluded, only raw evidence listed.
    What this is NOT: a measurement of the time between the last DMA and the end of the fused op (there are no
    per-DMA events); nothing is subtracted between programs. Full sizes need a TPU; --quick on CPU checks only
    that the code runs."""
    if jax.default_backend() != "tpu" and not ARGS.quick:
        raise SystemExit("tailprobe: no TPU backend (full sizes are TPU-only; --quick is the CPU functional check)")
    ARGS.keep_traces = True
    ITERS = 8                                     # loop iterations of the profiled program
    cases = [("bf16", 32, 512), ("int8", 32, 512)] if ARGS.quick else \
        [("bf16", 32, 4096), ("bf16", 64, 4096), ("bf16", 32, 8192), ("bf16", 64, 8192), ("int8", 32, 4096)]

    def loop(n, x, ws):                           # identical to balance(): same body, dependence and reduction
        def body(i, c):
            w = jax.lax.dynamic_index_in_dim(ws, i % ws.shape[0], keepdims=False)
            y = (x + c.astype(x.dtype)) @ w.astype(BF16)
            return c + jnp.sum(y.astype(F32)) * 1e-9
        return jax.lax.fori_loop(0, n, body, jnp.float32(0))
    run_loop = jax.jit(loop, static_argnums=0)

    def stack(key, shape, dtype):
        if dtype == "int8":
            gen = lambda k: jax.random.randint(k, shape, -8, 8, jnp.int8)
        else:
            gen = lambda k: (jax.random.normal(k, shape, F32) * 0.02).astype(BF16)
        return jax.jit(gen)(jax.random.PRNGKey(key))

    def walls(x, ws, n, repeats=5):               # raw samples, nothing clipped
        def wall(k):
            t = time.perf_counter()
            jax.block_until_ready(run_loop(k, x, ws))
            return time.perf_counter() - t
        for k in (n, 2 * n):
            wall(k)
        return [wall(n) for _ in range(repeats)], [wall(2 * n) for _ in range(repeats)]

    os.makedirs(os.path.join(ARGS.out, "hlo"), exist_ok=True)
    rows, pairs = [], []
    for dtype, M, K in cases:
        wbytes = 1 if dtype == "int8" else 2
        R = max(2, math.ceil(max(ARGS.rotate_mb, 1) * 2 ** 20 / (K * K * wbytes)))
        x = stack(M + 7, (M, K), "bf16")
        bank = stack(K + wbytes, (R, K, K), dtype)
        facts, evs = {}, {}
        for mode, ws in (("banked", bank), ("single", stack(K + wbytes + 100, (1, K, K), dtype))):
            tag = f"tailprobe_{dtype}_{M}x{K}_{mode}"
            try:
                text = run_loop.lower(ITERS, x, ws).compile().as_text()
            except Exception as e:
                text = ""
                print(f"  [{tag}] hlo: {type(e).__name__}: {e}")
            open(os.path.join(ARGS.out, "hlo", tag + ".txt"), "w").write(text)
            facts[mode] = hlo_facts(text)
            n = 10 if M * K * K > 2 ** 33 else 40
            wn, w2n = walls(x, ws, n)
            us = (statistics.median(w2n) - statistics.median(wn)) / n * 1e6
            r = bench(lambda a, b: run_loop(ITERS, a, b), (x, ws), tag)
            ev = evs[mode] = op_events(os.path.join(ARGS.out, "traces", tag), ITERS)
            allfus = ev and sum(t for k, (c, t) in ev["ops"].items() if "fusion" in k or "convolution" in k)
            rows.append(dict(
                dtype=dtype, M=M, K=K, N=K, mode=mode, weights_in_loop=ws.shape[0], loop_iterations_n=n,
                wall_n_s=json.dumps(wn), wall_2n_s=json.dumps(w2n), loop_slope_us=round(us, 3),
                slope_flag="ok" if us > 0 else "nonpositive: wall slope unusable",
                requested_modules=r["iters"], modules_traced=ev and ev["modules"],
                module_us_n8=r["device_us"] and round(r["device_us"], 3),
                op_event_counts=ev and ";".join(f"{k}={c}" for k, (c, t) in ev["ops"].items()),
                op_event_us=ev and ";".join(f"{k}={t:.3f}" for k, (c, t) in ev["ops"].items()),
                conv_fusion_events_per_module=ev and ev["modules"] and round(ev["conv_events"] / ev["modules"], 3),
                conv_fusion_us_per_event=round(ev["conv_us"] / ev["conv_events"], 3) if ev and ev["verified"] else None,
                all_fusion_event_sum_per_requested_iter_us=ev and ev["modules"] and round(allfus / ev["modules"] / ITERS, 3),
                hlo_facts=json.dumps(facts[mode]), hlo_sha256=hashlib.sha256(text.encode()).hexdigest(),
                trace_sha256=ev and json.dumps(ev["files"])))
            print(f"tailprobe {mode:6s} {dtype} {M:3d} x {K:5d}  slope {us:9.2f} us ({rows[-1]['slope_flag']})  "
                  f"conv fusion / event {rows[-1]['conv_fusion_us_per_event']}")
        verdict, notes = pair_verdict(facts["banked"], facts["single"], K, evs["banked"], evs["single"])
        pairs.append(dict(dtype=dtype, M=M, K=K, verdict=verdict, notes=" | ".join(notes),
                          windows_banked=json.dumps(facts["banked"]["windows"]),
                          windows_single=json.dumps(facts["single"]["windows"])))
        print(f"  -> {verdict}" + (f" ({'; '.join(notes)})" if notes else ""))
        del bank
    write("tailprobe.csv", rows)
    write("tailprobe_pairs.csv", pairs)
    dev = jax.devices()[0]
    json.dump(dict(script_sha256=hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
                   argv=sys.argv, quick=ARGS.quick, backend=jax.default_backend(), jax=jax.__version__,
                   device_kind=dev.device_kind, device=str(dev), rotate_mb=ARGS.rotate_mb, cases=cases,
                   profiled_loop_iterations=ITERS,
                   kinds=dict(wall="host wall clock samples and their signed slope",
                              events="profiler event counts / durations of the device (one event per op execution)",
                              hlo_facts="compiler output, static", verdict="from the two HLO texts and traces only"),
                   not_measured="time between the last DMA and the end of the fused op; measured HBM bytes; "
                                "where the reused weight physically resides",
                   cpu_note="a non-TPU run is a functional smoke check, not hardware validation"),
              open(os.path.join(ARGS.out, "tailprobe_meta.json"), "w"), indent=1)


def run(experiments, default_only, out_name, argv=None):
    """Parse the common options, write env.json and run the selected experiments (name -> function)."""
    global ARGS
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"/kaggle/working/{out_name}" if os.path.isdir("/kaggle/working")
                    else f"{out_name}_out")
    ap.add_argument("--quick", action="store_true", help="a few shapes per experiment (smoke test)")
    ap.add_argument("--keep-traces", action="store_true")
    ap.add_argument("--only", default=default_only, help="comma list of " + ",".join(experiments))
    ap.add_argument("--max-iters", type=int, default=200)
    ap.add_argument("--rotate-mb", type=float, default=512,
                    help="rotate distinct input buffers totalling this much (0: reuse one set, warm on-chip)")
    ap.add_argument("--device-match", default=r"TPU:0\b", help="regex for the trace process to read")
    ARGS = ap.parse_args(argv)
    os.makedirs(ARGS.out, exist_ok=True)
    devs = jax.devices()
    env = dict(jax=jax.__version__, backend=jax.default_backend(), n_devices=len(devs),
               device_kind=devs[0].device_kind, devices=[str(d) for d in devs],
               coords=[list(getattr(d, "coords", []) or []) for d in devs],   # TPU chip position in the torus
               date=time.strftime("%Y-%m-%d %H:%M:%S %Z"))
    json.dump(env, open(os.path.join(ARGS.out, "env.json"), "w"), indent=1)
    print(json.dumps(env))
    for name in ARGS.only.split(","):
        if name not in experiments:
            raise SystemExit(f"unknown experiment {name!r}; choose from {', '.join(experiments)}")
        t = time.time()
        experiments[name]()
        print(f"== {name}: {time.time() - t:.0f} s")
    print(f"results in {ARGS.out}")


def main(argv=None):
    run({"gemm": gemm, "vector": vector, "decode": decode, "smallk": smallk, "loopk": loopk, "balance": balance,
         "tailprobe": tailprobe},
        "gemm,vector,decode", "tpu_v5e1", argv)


if __name__ == "__main__":
    main()
