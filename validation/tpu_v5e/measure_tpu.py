#!/usr/bin/env python3
"""Cloud TPU measurements for NPUsim fidelity (validation/tpu_v5e/README.md). JAX only, no npusim needed:
copy this one file to a Kaggle TPU notebook (TpuV5E8) or a TRC/Cloud TPU VM and run it.

  gemm.csv         one chip: C = A @ B in bf16, square / decode / prefill / odd shapes
  vector.csv       one chip: rmsnorm, softmax, SwiGLU, residual add
  collectives.csv  all chips: psum, all_gather, psum_scatter, all_to_all vs bytes per chip
  decode.csv       one Llama-3-8B-shaped decode layer, on one chip and tensor-parallel over all chips,
                   time per named scope (qkv, attn, o_proj, ffn_*, norms, psum)
  env.json         devices, versions

Device time comes from the profiler's perfetto trace (XLA module and op events on TPU:0); wall time
(async-pipelined calls, block at the end) is recorded too and is the fallback when no trace is found.

Usage: python3 measure_tpu.py [--out DIR] [--quick] [--keep-traces] [--only gemm,vector,...]
"""
import argparse
import csv
import glob
import gzip
import json
import math
import os
import re
import shutil
import statistics
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


def collectives():
    devs = jax.devices()
    n = len(devs)
    if n < 2:
        print("collectives: one device, skipped")
        return
    mesh = Mesh(devs, ("x",))
    ops = {"psum": lambda x: jax.lax.psum(x, "x"),
           "all_gather": lambda x: jax.lax.all_gather(x, "x", tiled=True),
           "psum_scatter": lambda x: jax.lax.psum_scatter(x, "x", tiled=True),
           "all_to_all": lambda x: jax.lax.all_to_all(x, "x", 0, 0, tiled=True)}
    rows = []
    for lg in ((10, 20) if ARGS.quick else range(10, 27, 2)):        # 1 KiB .. 64 MiB per chip
        elems = max(n, (2 ** lg // 2) // n * n)
        x = rnd(5, (n * elems,), NamedSharding(mesh, P("x")))
        for name, op in ops.items():
            f = jax.jit(shard_map(op, mesh=mesh, in_specs=P("x"), out_specs=P("x")))
            r = bench(f, (x,), f"coll_{name}_{lg}")
            rows.append(dict(op=name, chips=n, bytes_per_chip=elems * 2, wall_us=round(r["wall_us"], 2),
                             device_us=r["device_us"] and round(r["device_us"], 2)))
            print(f"collective {name:12s} {elems * 2:>10d} B/chip  {fmt(r)}")
    write("collectives.csv", rows)


def decode():
    D, H, HKV, DH, F = 4096, 32, 8, 128, 14336            # Llama-3-8B layer
    devs = jax.devices()
    tps = [1] + ([len(devs)] if len(devs) > 1 and HKV % len(devs) == 0 else [])
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


def main():
    global ARGS
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/kaggle/working/tpu_v5e" if os.path.isdir("/kaggle/working") else "tpu_v5e_out")
    ap.add_argument("--quick", action="store_true", help="a few shapes per experiment (smoke test)")
    ap.add_argument("--keep-traces", action="store_true")
    ap.add_argument("--only", default="gemm,vector,collectives,decode")
    ap.add_argument("--max-iters", type=int, default=200)
    ap.add_argument("--rotate-mb", type=float, default=512,
                    help="rotate distinct input buffers totalling this much (0: reuse one set, warm on-chip)")
    ap.add_argument("--device-match", default=r"TPU:0\b", help="regex for the trace process to read")
    ARGS = ap.parse_args()
    os.makedirs(ARGS.out, exist_ok=True)
    devs = jax.devices()
    env = dict(jax=jax.__version__, backend=jax.default_backend(), n_devices=len(devs),
               device_kind=devs[0].device_kind, devices=[str(d) for d in devs],
               date=time.strftime("%Y-%m-%d %H:%M:%S %Z"))
    json.dump(env, open(os.path.join(ARGS.out, "env.json"), "w"), indent=1)
    print(json.dumps(env))
    for name in ARGS.only.split(","):
        t = time.time()
        globals()[name]()
        print(f"== {name}: {time.time() - t:.0f} s")
    print(f"results in {ARGS.out}")


if __name__ == "__main__":
    main()
