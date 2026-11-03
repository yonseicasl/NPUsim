# %% [markdown]
# # TPU v5e-8 measurements for NPUsim fidelity
#
# Run on a Kaggle notebook with the "TPU v5e-8" accelerator (Run All, ~30-45 min). Everything is written
# to /kaggle/working/tpu_v5e_results.json (plus a profiler trace in /kaggle/working/trace); download
# both from the notebook's Output tab.
#
# 1. single-chip GEMMs: decode-style GEMVs (M = 1..512) to compute-bound GEMMs, bf16 weights and int8
#    weights dequantized in the kernel
# 2. HBM bandwidth: read-only streams over large arrays
# 3. collectives over the 8 chips: all-reduce, all-gather, reduce-scatter, 4 KB .. 128 MB per chip
# 4. one Llama-2-70B-shaped decode layer, tensor-parallel over 8 chips, B = 1..128, context 2048, and
#    variants without collectives / attention / vector ops that split its time into components
#
# Timing: each op runs n and 2n times inside one jitted fori_loop whose iterations depend on each
# other (no hoisting, no host dispatch in the loop); time per iteration = (t(2n) - t(n)) / n, the
# median of several runs. Each iteration consumes all of the previous result (a full reduction): a scalar
# carry taken from one element would let XLA sink that slice into the op and compute one element.
#
# Cold operands: TPU v5e keeps an operand of up to ~128 MiB on chip when every iteration reads the same
# buffer, so a loop over one weight measures on-chip reuse, not HBM. Iteration i reads copy i % R of a
# stack of R distinct buffers (weights, KV cache, streamed arrays) totalling >= TPU_BENCH_ROTATE_MB
# (default 512 MB), as consecutive layers of a real model do. TPU_BENCH_ROTATE_MB=0 reuses one buffer
# (warm). Each result records its `copies`.

# %%
import json, math, os, platform, statistics, sys, time, traceback
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
ROTATE_MB = float(os.environ.get("TPU_BENCH_ROTATE_MB", "512"))
os.makedirs(OUT, exist_ok=True)
devs = jax.devices()
RESULTS = {"env": {"jax": jax.__version__, "python": sys.version.split()[0], "platform": platform.platform(),
                   "devices": [str(d) for d in devs], "device_kind": devs[0].device_kind,
                   "n_devices": len(devs), "quick": QUICK, "rotate_mb": ROTATE_MB},
           "gemm": [], "hbm": [], "collective": [], "decode_layer": [], "errors": []}
print(RESULTS["env"])
if jax.devices()[0].platform != "tpu" and not QUICK and os.environ.get("TPU_BENCH_ALLOW_CPU") != "1":
    # JAX silently falls back to the CPU when it finds no TPU; full sizes there take hours and measure nothing.
    raise SystemExit(f"No TPU: JAX sees {jax.devices()}. Select a TPU runtime (Colab: Runtime -> Change runtime "
                     "type -> v5e-1 TPU; Kaggle: Accelerator TPU v5e-8). If one is selected, this JAX has no TPU "
                     "support: run  !pip install -U 'jax[tpu]'  then restart the session. "
                     "(TPU_BENCH_QUICK=1 or TPU_BENCH_ALLOW_CPU=1 to run on the CPU anyway.)")


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


def copies(nbytes):
    """How many distinct copies of an operand of nbytes (per chip) to cycle through: >= ROTATE_MB in total, so
    no copy is still on chip when it is read again; 1 if one copy is already that large or ROTATE_MB is 0."""
    return 1 if ROTATE_MB <= 0 or nbytes >= ROTATE_MB * 2 ** 20 else math.ceil(ROTATE_MB * 2 ** 20 / nbytes)


def make(shape, dtype, sharding, seed, fill=None):
    """A random (or constant) array generated on the device(s), so large stacks never pass through the host."""
    def gen(k):
        if fill is not None:
            return jnp.full(shape, fill, dtype)
        if dtype == jnp.int8:
            return jax.random.randint(k, shape, -8, 8, jnp.int8)
        return (jax.random.normal(k, shape, jnp.float32) * 0.02).astype(dtype)
    return jax.jit(gen, out_shardings=sharding)(jax.random.PRNGKey(seed))


def pick(stack, i):
    """Copy i % R of a stack of R copies (dynamic, so XLA cannot hoist it out of the loop)."""
    return lax.dynamic_index_in_dim(stack, i % stack.shape[0], keepdims=False)


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
def gemm_loop(n, x, ws):
    def body(i, c):
        y = (x + c.astype(x.dtype)) @ pick(ws, i).astype(jnp.bfloat16)
        return c + jnp.sum(y.astype(jnp.float32)) * 1e-9          # all of y, so no slice sinks into the dot
    return lax.fori_loop(0, n, body, jnp.float32(0))


def gemm(M, K, N, wdtype):
    x = jax.device_put(jnp.asarray(rng.standard_normal((M, K), dtype=np.float32), jnp.bfloat16), dev0)
    wt = jnp.int8 if wdtype == "int8" else jnp.bfloat16
    wbytes = K * N * (1 if wdtype == "int8" else 2)
    R = copies(wbytes)
    ws = make((R, K, N), wt, jax.sharding.SingleDeviceSharding(dev0), seed=K + N + M)
    s = per_iter(gemm_loop, x, ws, n=10 if M * K * N > 2 ** 33 else 40)
    del ws
    return {"seconds": s, "tflops": 2 * M * K * N / s / 1e12 if s else None,
            "weight_gbps": wbytes / s / 1e9 if s else None, "copies": R}


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
def stream_loop(n, stack):
    def body(i, c):
        return c + jnp.sum(pick(stack, i).astype(jnp.float32)) * 1e-12   # reads the whole copy
    return lax.fori_loop(0, n, body, jnp.float32(0))


def hbm(mbytes):
    R = copies(mbytes * 2 ** 20)
    stack = make((R, mbytes * 2 ** 20 // 2), jnp.bfloat16, jax.sharding.SingleDeviceSharding(dev0), 0, fill=1)
    s = per_iter(stream_loop, stack, n=20)
    del stack
    return {"seconds": s, "read_gbps": mbytes * 2 ** 20 / s / 1e9 if s else None, "copies": R}


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
                return c + jnp.sum(y.astype(jnp.float32)) * 1e-12
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
    sh = lambda spec: jax.sharding.NamedSharding(mesh, spec)
    # one layer's global shapes and shardings; R layers of them stacked on a leading axis, iteration i
    # running layer i % R (distinct weights and KV cache, as consecutive layers of the model)
    shapes = {"wqkv": ((E, ND * (h + 2 * kvh) * DH), wt, P(None, "i")), "wo": ((ND * h * DH, E), wt, P("i", None)),
              "wg": ((E, F), wt, P(None, "i")), "wu": ((E, F), wt, P(None, "i")), "wd": ((F, E), wt, P("i", None)),
              "k": ((B, ND * kvh, CTX, DH), jnp.bfloat16, P(None, "i")),
              "v": ((B, ND * kvh, CTX, DH), jnp.bfloat16, P(None, "i"))}
    per_chip = sum(math.prod(s) * jnp.dtype(t).itemsize for s, t, _ in shapes.values()) // ND
    R = copies(per_chip)
    specs = {name: P(None, *spec) for name, (_, _, spec) in shapes.items()}
    params = {name: make((R, *s), t, sh(specs[name]), seed, fill=1 if name in ("k", "v") else None)
              for seed, (name, (s, t, _)) in enumerate(shapes.items())}
    x0 = jax.device_put(jnp.ones((B, E), jnp.bfloat16) * 0.01, sh(P()))
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
            return lax.fori_loop(0, n, lambda i, xx: layer(xx, jax.tree_util.tree_map(lambda a: pick(a, i), p)), x)
        return shard_map(local, mesh=mesh, in_specs=(P(), specs), out_specs=P())(x, params)

    s = per_iter(run, x0, params, n=n_layers)
    if trace:
        jax.profiler.start_trace(os.path.join(OUT, "trace"))
        jax.block_until_ready(run(n_layers, x0, params))
        jax.profiler.stop_trace()
    del params
    return {"seconds_per_layer": s, "copies": R}


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
