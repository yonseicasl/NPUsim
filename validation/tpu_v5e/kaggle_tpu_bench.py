# %% [markdown]
# # TPU v5e-8 measurements for NPUsim fidelity
#
# Run on a Kaggle notebook with the "TPU v5e-8" accelerator (Run All, ~30-45 min). Everything is written
# to /kaggle/working/tpu_v5e_results.json (plus a profiler trace in /kaggle/working/trace); download
# both from the notebook's Output tab.
#
# 1. single-chip GEMMs: decode-style GEMVs (M = 1..512) to compute-bound GEMMs, bf16 weights and int8
#    weights dequantized in the kernel
# 2. HBM bandwidth: elementwise streams over large arrays
# 3. collectives over the 8 chips: all-reduce, all-gather, reduce-scatter, 4 KB .. 128 MB per chip
# 4. one Llama-2-70B-shaped decode layer, tensor-parallel over 8 chips, B = 1..128, context 2048, and
#    variants without collectives / attention / vector ops that split its time into components
#
# Timing: each op runs n and 2n times inside one jitted fori_loop whose iterations depend on each
# other (no hoisting, no host dispatch in the loop); time per iteration = (t(2n) - t(n)) / n, the
# median of several runs.

# %%
import json, os, platform, statistics, sys, time, traceback
from functools import partial

import jax
import jax.numpy as jnp
import numpy as np
from jax import lax
from jax.sharding import Mesh, PartitionSpec as P

try:
    from jax import shard_map as _shard_map         # newer JAX
except ImportError:
    from jax.experimental.shard_map import shard_map as _shard_map


def shard_map(f, **kw):
    """shard_map without the replication check (its flag is check_vma in newer JAX, check_rep before)."""
    for flag in ("check_vma", "check_rep"):
        try:
            return _shard_map(f, **kw, **{flag: False})
        except TypeError:
            pass
    return _shard_map(f, **kw)

OUT = os.environ.get("TPU_BENCH_OUT", "/kaggle/working")
QUICK = os.environ.get("TPU_BENCH_QUICK") == "1"    # tiny sizes (CPU smoke test)
os.makedirs(OUT, exist_ok=True)
devs = jax.devices()
RESULTS = {"env": {"jax": jax.__version__, "python": sys.version.split()[0], "platform": platform.platform(),
                   "devices": [str(d) for d in devs], "device_kind": devs[0].device_kind,
                   "n_devices": len(devs), "quick": QUICK},
           "gemm": [], "hbm": [], "collective": [], "decode_layer": [], "errors": []}
print(RESULTS["env"])


def save():
    with open(os.path.join(OUT, "tpu_v5e_results.json"), "w") as f:
        json.dump(RESULTS, f, indent=1)


def per_iter(run, *args, n=20, repeats=5):
    """Seconds per loop iteration of run(n, *args) (n static), from the slope between n and 2n."""
    for k in (n, 2 * n):                            # compile both and warm up
        jax.block_until_ready(run(k, *args))
    def wall(k):
        t = time.perf_counter()
        jax.block_until_ready(run(k, *args))
        return time.perf_counter() - t
    t1 = statistics.median(wall(n) for _ in range(repeats))
    t2 = statistics.median(wall(2 * n) for _ in range(repeats))
    return max(t2 - t1, 0.0) / n


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


# %% [markdown]
# ## 1. Single-chip GEMMs

# %%
dev0 = devs[0]
rng = np.random.default_rng(0)


@partial(jax.jit, static_argnums=0)
def gemm_loop(n, x, w):
    def body(i, c):
        y = (x + c.astype(x.dtype)) @ w.astype(jnp.bfloat16)
        return c + y[0, 0].astype(jnp.float32) * 1e-6
    return lax.fori_loop(0, n, body, jnp.float32(0))


def gemm(M, K, N, wdtype):
    x = jax.device_put(jnp.asarray(rng.standard_normal((M, K), dtype=np.float32), jnp.bfloat16), dev0)
    w = rng.integers(-8, 8, (K, N), dtype=np.int8) if wdtype == "int8" else rng.standard_normal((K, N), dtype=np.float32)
    w = jax.device_put(jnp.asarray(w, jnp.int8 if wdtype == "int8" else jnp.bfloat16), dev0)
    s = per_iter(gemm_loop, x, w, n=10 if M * K * N > 2 ** 33 else 40)
    wbytes = K * N * (1 if wdtype == "int8" else 2)
    return {"seconds": s, "tflops": 2 * M * K * N / s / 1e12 if s else None,
            "weight_gbps": wbytes / s / 1e9 if s else None}


KN = [(1024, 1024)] if QUICK else [(4096, 4096), (8192, 8192), (8192, 28672), (28672, 8192)]
MS = [1, 16] if QUICK else [1, 8, 32, 64, 128, 256, 512]
for K, N in KN:
    for M in MS:
        for wd in ("bf16", "int8"):
            guarded("gemm", {"M": M, "K": K, "N": N, "weights": wd}, partial(gemm, M, K, N, wd))
for S in ([256] if QUICK else [1024, 2048, 4096, 8192]):
    guarded("gemm", {"M": S, "K": S, "N": S, "weights": "bf16"}, partial(gemm, S, S, S, "bf16"))


# %% [markdown]
# ## 2. HBM bandwidth

# %%
@partial(jax.jit, static_argnums=0)
def stream_loop(n, a):
    def body(i, c):
        b = a * (1 + c).astype(a.dtype)             # read a, write b
        return c + b[0].astype(jnp.float32) * 1e-9
    return lax.fori_loop(0, n, body, jnp.float32(0))


def hbm(mbytes):
    a = jax.device_put(jnp.ones((mbytes * 2 ** 20 // 2,), jnp.bfloat16), dev0)
    s = per_iter(stream_loop, a, n=20)
    return {"seconds": s, "gbps_read_plus_write": 2 * mbytes * 2 ** 20 / s / 1e9 if s else None}


for mb in ([1, 8] if QUICK else [1, 16, 128, 1024, 4096]):
    guarded("hbm", {"mbytes": mb}, partial(hbm, mb))


# %% [markdown]
# ## 3. Collectives over all chips

# %%
mesh = Mesh(np.array(devs), ("i",))
ND = len(devs)


def collective(kind, per_chip_bytes):
    n_el = max(per_chip_bytes // 2, ND)
    n_el -= n_el % ND
    x = jax.device_put(jnp.ones((ND * n_el,), jnp.bfloat16), jax.sharding.NamedSharding(mesh, P("i")))
    op = {"all_reduce": lambda v: lax.psum(v, "i"),
          "all_gather": lambda v: lax.all_gather(v, "i", tiled=True),
          "reduce_scatter": lambda v: lax.psum_scatter(v, "i", tiled=True)}[kind]

    @partial(jax.jit, static_argnums=0)
    def run(n, x):
        def local(xl):
            def body(i, c):
                y = op(xl + c.astype(xl.dtype))
                return c + y[0].astype(jnp.float32) * 1e-9
            return lax.fori_loop(0, n, body, jnp.float32(0))
        return shard_map(local, mesh=mesh, in_specs=P("i"), out_specs=P())(x)
    s = per_iter(run, x, n=20)
    return {"seconds": s, "algbw_gbps": per_chip_bytes / s / 1e9 if s else None}


for kind in ("all_reduce", "all_gather", "reduce_scatter"):
    for b in ([4096, 2 ** 20] if QUICK else [4096, 65536, 2 ** 20, 2 ** 24, 2 ** 27]):
        guarded("collective", {"kind": kind, "per_chip_bytes": b, "chips": ND}, partial(collective, kind, b))


# %% [markdown]
# ## 4. Decode layer, tensor-parallel over the chips (Llama-2-70B shape)

# %%
CFG = dict(E=256, F=512, H=8, KVH=8, DH=32, CTX=128) if QUICK else \
      dict(E=8192, F=28672, H=64, KVH=8, DH=128, CTX=2048)


def decode_layer(B, wdtype, comm=True, attn=True, vector=True, n_layers=4, trace=False):
    E, F, H, KVH, DH, CTX = (CFG[k] for k in ("E", "F", "H", "KVH", "DH", "CTX"))
    h, kvh = H // ND, max(KVH // ND, 1)             # heads per chip
    wt = jnp.int8 if wdtype == "int8" else jnp.bfloat16
    def W(*shape):                                   # global weight, sharded by the spec given at use
        return np.asarray(rng.integers(-8, 8, shape, dtype=np.int8) if wdtype == "int8" else rng.standard_normal(shape, dtype=np.float32) * 0.02)
    sh = lambda spec: jax.sharding.NamedSharding(mesh, spec)
    params = {
        "wqkv": jax.device_put(jnp.asarray(W(E, ND * (h + 2 * kvh) * DH), wt), sh(P(None, "i"))),
        "wo": jax.device_put(jnp.asarray(W(ND * h * DH, E), wt), sh(P("i", None))),
        "wg": jax.device_put(jnp.asarray(W(E, F), wt), sh(P(None, "i"))),
        "wu": jax.device_put(jnp.asarray(W(E, F), wt), sh(P(None, "i"))),
        "wd": jax.device_put(jnp.asarray(W(F, E), wt), sh(P("i", None))),
        "k": jax.device_put(jnp.ones((B, ND * kvh, CTX, DH), jnp.bfloat16), sh(P(None, "i"))),
        "v": jax.device_put(jnp.ones((B, ND * kvh, CTX, DH), jnp.bfloat16), sh(P(None, "i"))),
    }
    x0 = jax.device_put(jnp.ones((B, E), jnp.bfloat16) * 0.01, sh(P()))
    specs = {"wqkv": P(None, "i"), "wo": P("i", None), "wg": P(None, "i"), "wu": P(None, "i"),
             "wd": P("i", None), "k": P(None, "i"), "v": P(None, "i")}
    red = (lambda v: lax.psum(v, "i")) if comm else (lambda v: v)
    norm = (lambda v: v * lax.rsqrt(jnp.mean(v.astype(jnp.float32) ** 2, -1, keepdims=True) + 1e-6).astype(v.dtype)) \
        if vector else (lambda v: v)
    act = (lambda g, u: jax.nn.silu(g) * u) if vector else (lambda g, u: g * u)
    bf = lambda w: w.astype(jnp.bfloat16)

    def layer(x, p):
        y = norm(x) @ bf(p["wqkv"])
        q, k_new, v_new = jnp.split(y, [h * DH, (h + kvh) * DH], axis=-1)
        if attn:
            q = q.reshape(B, kvh, h // kvh, DH)
            s = jnp.einsum("bkgd,bktd->bkgt", q, p["k"]) / float(np.sqrt(DH))   # python float keeps bf16
            pr = jax.nn.softmax(s.astype(jnp.float32), -1).astype(jnp.bfloat16) if vector else s
            o = jnp.einsum("bkgt,bktd->bkgd", pr, p["v"]).reshape(B, h * DH)
        else:
            o = q
        x = x + red(o @ bf(p["wo"]))
        z = norm(x)
        return (x + red(act(z @ bf(p["wg"]), z @ bf(p["wu"])) @ bf(p["wd"]))).astype(jnp.bfloat16)

    @partial(jax.jit, static_argnums=0)
    def run(n, x, params):
        def local(x, p):
            return lax.fori_loop(0, n, lambda i, xx: layer(xx, p), x)
        return shard_map(local, mesh=mesh, in_specs=(P(), specs), out_specs=P())(x, params)

    s = per_iter(run, x0, params, n=n_layers)
    if trace:
        jax.profiler.start_trace(os.path.join(OUT, "trace"))
        jax.block_until_ready(run(n_layers, x0, params))
        jax.profiler.stop_trace()
    return {"seconds_per_layer": s}


for wd in ("bf16", "int8"):
    for B in ([1, 8] if QUICK else [1, 8, 32, 64, 128]):
        base = {"B": B, "weights": wd, "chips": ND, **CFG}
        guarded("decode_layer", dict(base, variant="full"),
                partial(decode_layer, B, wd, trace=(B == 8 and wd == "bf16")))
        guarded("decode_layer", dict(base, variant="no_comm"), partial(decode_layer, B, wd, comm=False))
        guarded("decode_layer", dict(base, variant="no_attn"), partial(decode_layer, B, wd, attn=False))
        guarded("decode_layer", dict(base, variant="no_vector"), partial(decode_layer, B, wd, vector=False))
        guarded("decode_layer", dict(base, variant="gemm_only"),
                partial(decode_layer, B, wd, comm=False, attn=False, vector=False))

# %%
save()
print("errors:", len(RESULTS["errors"]))
print("wrote", os.path.join(OUT, "tpu_v5e_results.json"))
