# %% [markdown]
# # TPU v5e-1 (Colab) follow-up measurements for NPUsim
#
# Runtime -> Change runtime type -> **v5e-1 TPU**, then Run All (~20-30 min). Results go to
# /content/tpu_v5e1_followup/ and are zipped for download at the end.
#
# Follows up the isolated Colab runs of `measure_tpu.py` (results/2026-11-03_colab_v5e1_*), which time one
# GEMM or one decode layer per XLA program. With fresh buffers such a program copies an operand HBM -> on-chip
# memory and only then runs the fusion, serially, so it overstates an HBM-bound GEMM; a real decode step
# streams many different weights inside one program. Everything here runs inside one program:
#
# 1. stream_gemm: GEMMs over R distinct weights in one loop (cold, the decode case), and over one weight
#    (warm, on-chip reuse) for contrast; bf16 and int8 weights (dequantized in the kernel)
# 2. decode_stack: R distinct Llama-3-8B decode layers (weights + KV cache) in one loop, full and GEMM-only
# 3. attention: decode attention alone (GQA 32/8 heads, d 128) over a rotated KV cache, vs its KV-read roofline
# 4. hbm_read: read-only streams 64 MB .. 4 GB (sustained HBM read bandwidth; spec 819 GB/s)
# 5. odd_gemm: dimensions around multiples of 128 (padding/tiling granularity)
#
# TPU_BENCH_SECTIONS (comma list, default all) runs only some sections; in Colab put
#   %env TPU_BENCH_SECTIONS=decode_stack,attention,hbm_read
# in a cell above this one.
#
# Also, for every case (TPU_BENCH_TRACE=all, the default; `some` = a few, `none`): the optimized HLO, XLA's cost
# analysis and a profiler trace (hlo/, traces/, cost_analysis). xprof_stats.py turns the traces into per-op FLOPs,
# bytes, FLOP rate, HBM / VMEM bandwidth and bound-by (XProf hlo_stats / roofline / op_profile). The `counters`
# section traces one GEMM under candidate profiler options to find one that fills XProf's perf_counters
# (hardware counters; empty with the default options). The HBM memory stats of the run (env.memory_stats).
#
# Timing: each case runs n and 2n dependent iterations inside one jitted fori_loop; time per iteration =
# (t(2n) - t(n)) / n, median of 5 (no host dispatch). Each iteration consumes the whole previous result (a full
# reduction or the output itself), so XLA cannot shrink the work to the elements a scalar carry would use.
# Iteration i reads copy i % R of R distinct buffers totalling >= TPU_BENCH_ROTATE_MB (512 MB), so nothing is
# served from on-chip memory across iterations; `copies` records R.
#
# `gemm.csv` holds the cold bf16 GEMMs in measure_tpu.py's format, so `compare_v5e.py <dir>` compares them with
# NPUsim directly; `colab_summary.py <dir>` prints the rest.

# %%
import csv, json, math, os, platform, shutil, statistics, sys, time, traceback
from functools import partial

import jax
import jax.numpy as jnp
import numpy as np
from jax import lax

OUT = os.environ.get("TPU_BENCH_OUT", "/content/tpu_v5e1_followup")
QUICK = os.environ.get("TPU_BENCH_QUICK") == "1"    # tiny sizes (CPU smoke test)
ROTATE_MB = float(os.environ.get("TPU_BENCH_ROTATE_MB", "512"))
SECTIONS = os.environ.get("TPU_BENCH_SECTIONS",
                          "stream_gemm,decode_stack,attention,hbm_read,odd_gemm,counters").split(",")
TRACE = os.environ.get("TPU_BENCH_TRACE", "all")      # all | some | none
os.makedirs(OUT, exist_ok=True)
dev = jax.devices()[0]
ONE = jax.sharding.SingleDeviceSharding(dev)
RESULTS = {"env": {"jax": jax.__version__, "python": sys.version.split()[0], "platform": platform.platform(),
                   "device_kind": dev.device_kind, "n_devices": len(jax.devices()), "quick": QUICK,
                   "rotate_mb": ROTATE_MB, "sections": SECTIONS, "date": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())},
           "stream_gemm": [], "decode_stack": [], "attention": [], "hbm_read": [], "odd_gemm": [], "counters": [],
           "cost_analysis": {},
           "errors": []}
print(RESULTS["env"])
if jax.devices()[0].platform != "tpu" and not QUICK and os.environ.get("TPU_BENCH_ALLOW_CPU") != "1":
    # JAX silently falls back to the CPU when it finds no TPU; full sizes there take hours and measure nothing.
    raise SystemExit(f"No TPU: JAX sees {jax.devices()}. Select a TPU runtime (Colab: Runtime -> Change runtime "
                     "type -> v5e-1 TPU; Kaggle: Accelerator TPU v5e-8). If one is selected, this JAX has no TPU "
                     "support: run  !pip install -U 'jax[tpu]'  then restart the session. "
                     "(TPU_BENCH_QUICK=1 or TPU_BENCH_ALLOW_CPU=1 to run on the CPU anyway.)")


def save():
    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(RESULTS, f, indent=1)


def per_iter(run, *args, n=20, repeats=5):
    """Seconds per loop iteration of run(n, *args) (n static), from the slope between n and 2n."""
    for k in (n, 2 * n):
        jax.block_until_ready(run(k, *args))
    def wall(k):
        t = time.perf_counter()
        jax.block_until_ready(run(k, *args))
        return time.perf_counter() - t
    t1 = statistics.median(wall(n) for _ in range(repeats))
    t2 = statistics.median(wall(2 * n) for _ in range(repeats))
    return max(t2 - t1, 0.0) / n


def copies(nbytes, at_least=2):
    """Distinct copies to cycle through: >= ROTATE_MB in total, at least `at_least` (1 with ROTATE_MB 0)."""
    if ROTATE_MB <= 0:
        return 1
    return max(at_least, math.ceil(ROTATE_MB * 2 ** 20 / nbytes))


def make(shape, dtype, seed, fill=None):
    """A random (or constant) array generated on the TPU, so large stacks never pass through the host."""
    def gen(k):
        if fill is not None:
            return jnp.full(shape, fill, dtype)
        if dtype == jnp.int8:
            return jax.random.randint(k, shape, -8, 8, jnp.int8)
        return (jax.random.normal(k, shape, jnp.float32) * 0.02).astype(dtype)
    return jax.jit(gen, out_shardings=ONE)(jax.random.PRNGKey(seed))


def pick(stack, i):
    """Copy i % R of a stack of R copies (dynamic, so XLA cannot hoist it out of the loop). Fine for weights: the
    slice fuses into the dot. Not for large KV caches or streamed arrays: there XLA materializes the slice (an extra
    HBM round trip); those use UNROLL separate buffers in turn inside the loop body instead."""
    return lax.dynamic_index_in_dim(stack, i % stack.shape[0], keepdims=False)


UNROLL = 8   # at most this many distinct buffers are unrolled in a loop body; more (small ones) use pick()


def guarded(section, params, fn):
    try:
        out = fn()
        RESULTS[section].append(dict(params, **out))
        print(section, params, {k: (f"{v:.3e}" if isinstance(v, float) else v) for k, v in out.items()})
    except Exception as e:                          # keep going; record what failed
        RESULTS["errors"].append({"section": section, "params": params, "error": repr(e)[:500]})
        print("ERROR", section, params, repr(e)[:200])
        traceback.print_exc(limit=1)
    save()


def profiler_options(advanced):
    """jax.profiler.ProfileOptions with an advanced_configuration dict, or None if this JAX has no such option."""
    if not advanced or not hasattr(jax.profiler, "ProfileOptions"):
        return None
    o = jax.profiler.ProfileOptions()
    o.advanced_configuration = advanced
    return o


def trace(tag, run, *args, perfetto=False, advanced=None, wanted=True):
    """For one case: the optimized HLO (OUT/hlo/tag.txt: fusions, async copies, on-chip memory-space placement),
    XLA's cost analysis (flops, bytes accessed) and a profiler trace (OUT/traces/tag: xplane.pb, plus a perfetto
    copy if asked; XProf derives hlo_stats / op_profile / roofline / memory_profile from it). Skipped with
    TPU_BENCH_TRACE=none, and for cases not `wanted` with TPU_BENCH_TRACE=some."""
    if TRACE == "none" or (TRACE == "some" and not wanted):
        return
    try:
        c = run.lower(*args).compile()
        os.makedirs(os.path.join(OUT, "hlo"), exist_ok=True)
        with open(os.path.join(OUT, "hlo", tag + ".txt"), "w") as fh:
            fh.write(c.as_text())
        ca = c.cost_analysis()
        ca = ca[0] if isinstance(ca, (list, tuple)) and ca else ca
        RESULTS["cost_analysis"][tag] = {k: v for k, v in (ca or {}).items() if isinstance(v, (int, float))}
        kw = {"create_perfetto_trace": perfetto}
        opts = profiler_options(advanced)
        if opts is not None:
            kw["profiler_options"] = opts
        with jax.profiler.trace(os.path.join(OUT, "traces", tag), **kw):
            jax.block_until_ready(run(*args))
    except Exception as e:
        RESULTS["errors"].append({"section": "trace", "params": {"tag": tag}, "error": repr(e)[:300]})


# %% [markdown]
# ## 1. GEMMs streaming distinct weights (cold) vs one weight (warm)

# %%
@partial(jax.jit, static_argnums=0)
def gemm_loop(n, x, ws):
    def body(i, c):
        y = (x + c.astype(x.dtype)) @ pick(ws, i).astype(jnp.bfloat16)
        return c + jnp.sum(y.astype(jnp.float32)) * 1e-9
    return lax.fori_loop(0, n, body, jnp.float32(0))


def gemm(M, K, N, wdtype, cold=True, tag=None):
    wt = jnp.int8 if wdtype == "int8" else jnp.bfloat16
    wbytes = K * N * jnp.dtype(wt).itemsize
    R = copies(wbytes) if cold else 1
    x = make((M, K), jnp.bfloat16, M + 7)
    ws = make((R, K, N), wt, K + N)
    s = per_iter(gemm_loop, x, ws, n=10 if M * K * N > 2 ** 33 else 40)
    if tag:                                         # streamed loop vs the same GEMM as a one-op program
        trace(tag, gemm_loop, 8, x, ws, perfetto=True)
        trace(tag + "_isolated", jax.jit(lambda a, b: a @ b.astype(jnp.bfloat16)), x, ws[1 % R], perfetto=True)
    else:
        trace(f"gemm_{M}x{K}x{N}_{wdtype}_{'cold' if cold else 'warm'}", gemm_loop, 8, x, ws, wanted=False)
    del ws
    return {"us": s * 1e6, "tflops": 2 * M * K * N / s / 1e12 if s else None,
            "weight_gbps": wbytes / s / 1e9 if s else None, "copies": R}


D, H, HKV, DH, F = (256, 4, 2, 32, 512) if QUICK else (4096, 32, 8, 128, 14336)    # Llama-3-8B
PROJ = {"qkv": (D, (H + 2 * HKV) * DH), "o_proj": (H * DH, D), "ffn_up": (D, 2 * F), "ffn_down": (F, D)}
MS = [1, 8] if QUICK else [1, 8, 32, 128, 512]
for M in (MS if "stream_gemm" in SECTIONS else []):
    for name, (K, N) in PROJ.items():
        for wd in ("bf16", "int8"):
            for cold in (True, False):
                tag = f"stream_gemm_{M}x{K}x{N}_{wd}" if (M == 8 and name == "o_proj" and wd == "bf16" and cold) else None
                guarded("stream_gemm", {"M": M, "K": K, "N": N, "proj": name, "weights": wd, "cold": cold},
                        partial(gemm, M, K, N, wd, cold, tag))


# %% [markdown]
# ## 2. Decode layers, R distinct layers in one program (Llama-3-8B, one chip)

# %%
def decode_stack(B, L, gemm_only=False, n_layers=4, tag=None):
    shapes = {"g1": ((D,), jnp.bfloat16), "wqkv": (PROJ["qkv"], jnp.bfloat16), "wo": (PROJ["o_proj"], jnp.bfloat16),
              "g2": ((D,), jnp.bfloat16), "wgu": (PROJ["ffn_up"], jnp.bfloat16), "wd": (PROJ["ffn_down"], jnp.bfloat16),
              "k": ((B, HKV, L, DH), jnp.bfloat16), "v": ((B, HKV, L, DH), jnp.bfloat16)}
    per_layer = sum(math.prod(s) * jnp.dtype(t).itemsize for s, t in shapes.values())
    R = min(UNROLL, copies(per_layer))              # 2 at full size: a layer is >= 436 MB
    p = [{name: make(s, t, 100 * j + seed, fill=1 if name in ("g1", "g2", "k", "v") else None)
          for seed, (name, (s, t)) in enumerate(shapes.items())} for j in range(R)]
    rms = lambda v, g: (v * lax.rsqrt(jnp.mean(v.astype(jnp.float32) ** 2, -1, keepdims=True) + 1e-6)).astype(v.dtype) * g

    def layer(x, w):
        h = x if gemm_only else rms(x, w["g1"])
        q, kn, vn = jnp.split(h @ w["wqkv"], [H * DH, (H + HKV) * DH], axis=-1)
        if gemm_only:
            o = q
        else:
            qg = q.reshape(B, HKV, H // HKV, DH)
            s = jnp.einsum("bkgd,bkld->bkgl", qg, w["k"]) / float(np.sqrt(DH))
            pr = jax.nn.softmax(s.astype(jnp.float32), -1).astype(jnp.bfloat16)
            o = jnp.einsum("bkgl,bkld->bkgd", pr, w["v"]).reshape(B, H * DH)
        x = x + o @ w["wo"]
        h = x if gemm_only else rms(x, w["g2"])
        gu = h @ w["wgu"]
        a = gu[:, :F] * gu[:, F:] if gemm_only else jax.nn.silu(gu[:, :F]) * gu[:, F:]
        return (x + a @ w["wd"]).astype(jnp.bfloat16)

    @partial(jax.jit, static_argnums=0)
    def run(n, x, p):
        def body(i, xx):                            # the R distinct layers, one after the other
            for w in p:
                xx = layer(xx, w)
            return xx
        return lax.fori_loop(0, n, body, x)

    x0 = make((B, D), jnp.bfloat16, 99, fill=0.01)
    s = per_iter(run, x0, p, n=max(1, n_layers // R)) / R
    trace(tag or f"decode_b{B}_l{L}_{'gemm_only' if gemm_only else 'full'}", run, max(1, n_layers // R), x0, p,
          perfetto=bool(tag), wanted=bool(tag))
    del p
    return {"us_per_layer": s * 1e6, "copies": R, "layer_mb": per_layer / 2 ** 20}


for B in (([1, 8] if QUICK else [1, 8, 32, 128]) if "decode_stack" in SECTIONS else []):
    for L in ([128] if QUICK else [512, 2048, 8192]):
        tag = "decode_stack_b8_l2048" if (B == 8 and L == 2048) else None
        guarded("decode_stack", {"B": B, "L": L, "variant": "full"}, partial(decode_stack, B, L, False, tag=tag))
        guarded("decode_stack", {"B": B, "L": L, "variant": "gemm_only"}, partial(decode_stack, B, L, True))


# %% [markdown]
# ## 3. Decode attention alone over a rotated KV cache

# %%
def attend(q, k, v):
    s = jnp.einsum("bkgd,bkld->bkgl", q, k) / float(np.sqrt(DH))
    pr = jax.nn.softmax(s.astype(jnp.float32), -1).astype(jnp.bfloat16)
    return jnp.einsum("bkgl,bkld->bkgd", pr, v)     # the output is the next query


def attention(B, L):
    """KV in (B, KV heads, L, d), the layout serving stacks keep. Up to UNROLL distinct caches unrolled in the loop
    body; smaller caches (more copies) are picked from a stack, where the slice stays on chip."""
    kv_bytes = 2 * B * L * HKV * DH * 2
    R = copies(kv_bytes)
    q0 = make((B, HKV, H // HKV, DH), jnp.bfloat16, 3)
    if R <= UNROLL:
        kv = [(make((B, HKV, L, DH), jnp.bfloat16, 0, fill=1), make((B, HKV, L, DH), jnp.bfloat16, 0, fill=1))
              for _ in range(R)]

        @partial(jax.jit, static_argnums=0)
        def run(n, q, kv):
            def body(i, q):
                for k, v in kv:
                    q = attend(q, k, v)
                return q
            return lax.fori_loop(0, n, body, q)
        s = per_iter(run, q0, kv, n=max(2, 20 // R)) / R
        trace(f"attn_b{B}_l{L}", run, 2, q0, kv, wanted=False)
    else:
        kv = (make((R, B, HKV, L, DH), jnp.bfloat16, 1, fill=1), make((R, B, HKV, L, DH), jnp.bfloat16, 2, fill=1))

        @partial(jax.jit, static_argnums=0)
        def run(n, q, kv):
            return lax.fori_loop(0, n, lambda i, q: attend(q, pick(kv[0], i), pick(kv[1], i)), q)
        s = per_iter(run, q0, kv, n=20)
        trace(f"attn_b{B}_l{L}", run, 4, q0, kv, wanted=False)
    del kv
    return {"us": s * 1e6, "kv_gbps": kv_bytes / s / 1e9 if s else None, "copies": R,
            "mode": "unrolled" if R <= UNROLL else "stacked",
            "vs_roofline_819": s / (kv_bytes / 819e9) if s else None}


for B in (([1, 8] if QUICK else [1, 8, 32, 128]) if "attention" in SECTIONS else []):
    for L in ([128] if QUICK else [512, 2048, 8192, 32768]):
        if B * L <= 2 ** 20:                        # <= 4 GB of KV per layer
            guarded("attention", {"B": B, "L": L}, partial(attention, B, L))


# %% [markdown]
# ## 4. HBM read bandwidth

# %%
@partial(jax.jit, static_argnums=0)
def read_loop(n, arrays):
    def body(i, c):
        for a in arrays:                            # each distinct (rows, 4096) array read in full; the max with the
            # carry makes every read depend on the previous one, else XLA hoists the loop-invariant sum out of the loop
            c = c + jnp.sum(jnp.maximum(a, c.astype(a.dtype)).astype(jnp.float32)) * 1e-12
        return c
    return lax.fori_loop(0, n, body, jnp.float32(0))


def hbm_read(mbytes):
    """Read-only stream: R >= 2 distinct 2-D arrays (>= 512 MB) read in turn, unrolled (no slicing)."""
    nbytes = mbytes * 2 ** 20
    R = min(UNROLL, copies(nbytes))
    arrays = [make((nbytes // 2 // 4096, 4096), jnp.bfloat16, 0, fill=1) for _ in range(R)]
    s = per_iter(read_loop, arrays, n=max(2, 20 // R)) / R
    trace(f"hbm_{mbytes}mb", read_loop, 2, arrays, wanted=False)
    del arrays
    return {"us": s * 1e6, "read_gbps": nbytes / s / 1e9 if s else None, "copies": R}


for mb in (([1, 8] if QUICK else [64, 256, 1024, 4096]) if "hbm_read" in SECTIONS else []):
    guarded("hbm_read", {"mbytes": mb}, partial(hbm_read, mb))


# %% [markdown]
# ## 5. Dimensions around multiples of 128

# %%
for M in (([100, 128] if QUICK else [100, 128, 200, 256, 512]) if "odd_gemm" in SECTIONS else []):
    for KN in ([1000, 1024] if QUICK else [1000, 1024, 1100, 1152, 1280, 2000, 2048]):
        guarded("odd_gemm", {"M": M, "K": KN, "N": KN, "weights": "bf16", "cold": True},
                partial(gemm, M, KN, KN, "bf16", True))


# %% [markdown]
# ## Hardware counters: one GEMM traced under candidate profiler options
#
# XProf's perf_counters / utilization_viewer tools come back empty with the default options. These candidates (not
# documented; names seen in XProf / Cloud TPU profiler options) are tried on one case, each into traces/counters_<i>;
# xprof_stats.py --list-tools then shows which trace carries counters. A rejected option is recorded and skipped.

# %%
COUNTER_OPTIONS = [{}, {"tpu_trace_mode": "TRACE_COMPUTE_AND_SYNC"}, {"performance_counter_sets": "all"},
                   {"tpu_performance_counter_sets": "all"}, {"enable_perf_counters": True}]


def counters():
    M, K, N = (8, 256, 256) if QUICK else (8, 4096, 4096)
    x = make((M, K), jnp.bfloat16, 7)
    ws = make((copies(K * N * 2), K, N), jnp.bfloat16, 8)
    for i, adv in enumerate(COUNTER_OPTIONS):
        n_err = len(RESULTS["errors"])
        trace(f"counters_{i}", gemm_loop, 8, x, ws, advanced=adv)
        RESULTS["counters"].append({"trace": f"counters_{i}", "advanced_configuration": adv,
                                    "ok": len(RESULTS["errors"]) == n_err})
    del ws


if "counters" in SECTIONS and TRACE != "none":
    try:
        counters()
    except Exception as e:
        RESULTS["errors"].append({"section": "counters", "params": {}, "error": repr(e)[:500]})
    save()


# %% [markdown]
# ## Save, CSV for compare_v5e.py, zip

# %%
try:
    RESULTS["env"]["memory_stats"] = dev.memory_stats()        # HBM bytes_limit / peak_bytes_in_use
except Exception:
    pass
save()
with open(os.path.join(OUT, "gemm.csv"), "w", newline="") as fh:   # measure_tpu.py format (device_us = slope)
    w = csv.DictWriter(fh, ["M", "K", "N", "dtype", "sets", "wall_us", "device_us", "tflops", "ops"])
    w.writeheader()
    for r in RESULTS["stream_gemm"] + RESULTS["odd_gemm"]:
        if r["weights"] == "bf16" and r["cold"]:
            w.writerow(dict(M=r["M"], K=r["K"], N=r["N"], dtype="bf16", sets=r["copies"], wall_us=round(r["us"], 2),
                            device_us=round(r["us"], 2), tflops=round(r["tflops"] or 0, 2), ops=""))
with open(os.path.join(OUT, "env.json"), "w") as fh:
    json.dump(RESULTS["env"], fh, indent=1)
print("errors:", len(RESULTS["errors"]))
archive = shutil.make_archive(OUT, "zip", OUT)
print("wrote", archive)
try:
    from google.colab import files
    files.download(archive)
except Exception:
    pass
