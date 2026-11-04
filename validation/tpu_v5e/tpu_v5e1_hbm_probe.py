# %% [markdown]
# # TPU v5e-1 HBM probe: address mapping and refresh (Colab)
#
# Runtime -> Change runtime type -> **v5e-1 TPU**, then Run All (~15-25 min). Results go to
# /content/tpu_v5e1_hbm_probe/ and are zipped for download at the end.
#
# How does the TPU keep 89-92% of its HBM peak on long streams (validation/tpu_v5e/REPORT.md)? A refresh postponed
# and caught up later can only be hidden for ~31 us (8 x tREFI), yet a 5.7 ms stream kept 92%; refresh is either
# overlapped with reads from other banks (per-bank REFSB plus a mapping that spreads streams) or simply cheap. No
# DRAM counter is exposed, so this measures bandwidth under controlled access patterns, for comparison with MOSAIC's
# Ramulator2 HBM2e under candidate mappings and refresh schemes:
#
# 1. seq_duration: one sequential stream kept up for ~1 ms, ~10 ms, ~200 ms (does efficiency decay?)
# 2. stride: 4 KB chunks every S tiles, S = 1 .. 16384 (bandwidth dips = strides that alias channels / banks)
# 3. chunk: chunks of 4 .. 128 KB at a fixed 1 MB stride (channel interleave granularity)
# 4. burst: 4 / 16 / 64 MB bursts separated by ~20 us of HBM idle (VPU-only work) vs back-to-back. Refresh postponed
#    into idle time makes bursts faster than the stream; refresh hidden under other banks makes no difference.
# 5. xla_stride: section 2 with plain XLA slicing (a fallback if Pallas fails)
#
# Sections 1-4 use Pallas with explicit DMAs (HBM -> VMEM, `inflight` in flight), so the access pattern is exactly
# the one written: an array (G, S, 8, 128) f32 is read as chunks [g, 0:c] -- one (8, 128) f32 tile is 4 KB, the
# HBM layout unit. Timing: n and 2n passes inside one kernel, time per pass = (t(2n) - t(n)) / n, median of 5.
# TPU_BENCH_SECTIONS (comma list) runs only some sections.

# %%
import json, math, os, platform, shutil, statistics, sys, time, traceback
from functools import partial

import jax
import jax.numpy as jnp
from jax import lax
from jax.experimental import pallas as pl
from jax.experimental.pallas import tpu as pltpu

OUT = os.environ.get("TPU_BENCH_OUT", "/content/tpu_v5e1_hbm_probe")
QUICK = os.environ.get("TPU_BENCH_QUICK") == "1"   # tiny sizes, Pallas interpret mode (CPU smoke test)
SECTIONS = os.environ.get("TPU_BENCH_SECTIONS", "seq_duration,stride,chunk,burst,xla_stride").split(",")
os.makedirs(OUT, exist_ok=True)
dev = jax.devices()[0]
RESULTS = {"env": {"jax": jax.__version__, "python": sys.version.split()[0], "platform": platform.platform(),
                   "device_kind": dev.device_kind, "quick": QUICK, "sections": SECTIONS,
                   "date": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())},
           "seq_duration": [], "stride": [], "chunk": [], "burst": [], "xla_stride": [], "errors": []}
print(RESULTS["env"])
if dev.platform != "tpu" and not QUICK and os.environ.get("TPU_BENCH_ALLOW_CPU") != "1":
    raise SystemExit(f"No TPU: JAX sees {jax.devices()}. Colab: Runtime -> Change runtime type -> v5e-1 TPU; if one "
                     "is selected, run  !pip install -U 'jax[tpu]'  and restart the session.")

TILE = 8 * 128 * 4                                  # bytes of one (8, 128) f32 tile
ARRAY_BYTES = (8 << 20) if QUICK else (2 << 30)     # the probed array
INTERPRET = pltpu.InterpretParams() if QUICK else False


def save():
    with open(os.path.join(OUT, "results.json"), "w") as f:
        json.dump(RESULTS, f, indent=1)


def guarded(section, params, fn):
    try:
        out = fn()
        RESULTS[section].append(dict(params, **out))
        print(section, params, {k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in out.items()})
    except Exception as e:
        RESULTS["errors"].append({"section": section, "params": params, "error": repr(e)[:600]})
        print("ERROR", section, params, repr(e)[:300])
        traceback.print_exc(limit=2)
    save()


def per_pass(run, n, *args, repeats=5):
    """Seconds per pass of run(n, *args) (n static), from the slope between n and 2n passes. The arrays go in as
    arguments: an array the function closes over becomes a compile-time constant (XLA may fold the reads away)."""
    for k in (n, 2 * n):
        jax.block_until_ready(run(k, *args))
    def wall(k):
        t = time.perf_counter()
        jax.block_until_ready(run(k, *args))
        return time.perf_counter() - t
    t1 = statistics.median(wall(n) for _ in range(repeats))
    t2 = statistics.median(wall(2 * n) for _ in range(repeats))
    return max(t2 - t1, 1e-12) / n


def passes_for(seconds_per_pass_guess, target_s):
    return max(1, int(target_s / max(seconds_per_pass_guess, 1e-9)))


# %% [markdown]
# ## The DMA stream kernel

# %%
def make_stream(G, S, c, k, inflight, dmas_per_burst, spin_iters, with_dma=True):
    """run(n, x) -> (8, 128) over x = array(G, S, 8, 128): n bursts; a burst issues dmas_per_burst DMAs of k chunks of c tiles (chunk g of the
    (G, S, 8, 128) array = x[g, 0:c]), `inflight` in flight, cycling through the array, then spins spin_iters VPU
    iterations with HBM idle. with_dma=False times the spin alone."""
    groups = G // k

    def kernel(n_bursts, x_hbm, o_ref, buf, sem):
        def burst(b, acc):
            if with_dma:
                def copy(j, slot):
                    g = ((b * dmas_per_burst + j) % groups) * k
                    return pltpu.make_async_copy(x_hbm.at[pl.ds(g, k), pl.ds(0, c)], buf.at[slot], sem.at[slot])
                for s in range(min(inflight, dmas_per_burst)):
                    copy(s, s).start()

                def body(j, acc):
                    slot = j % inflight
                    copy(j, slot).wait()
                    acc = acc + buf[slot, 0, 0]
                    @pl.when(j + inflight < dmas_per_burst)
                    def _():
                        copy(j + inflight, slot).start()
                    return acc
                acc = lax.fori_loop(0, dmas_per_burst, body, acc)
            if spin_iters:
                acc = lax.fori_loop(0, spin_iters, lambda i, a: a * 1.0000001 + 1e-7, acc)
            return acc
        o_ref[...] = lax.fori_loop(0, n_bursts, burst, jnp.zeros((8, 128), jnp.float32))

    @partial(jax.jit, static_argnums=0)
    def run(n, x):
        return pl.pallas_call(
            partial(kernel, n), out_shape=jax.ShapeDtypeStruct((8, 128), jnp.float32),
            in_specs=[pl.BlockSpec(memory_space=pltpu.ANY)],
            scratch_shapes=[pltpu.VMEM((inflight, k, c, 8, 128), jnp.float32), pltpu.SemaphoreType.DMA((inflight,))],
            interpret=INTERPRET)(x)
    return run


def array(G, S):
    return jax.jit(lambda: jnp.ones((G, S, 8, 128), jnp.float32))()


def stream_bw(S, c, target_s=0.03, inflight=4, dma_bytes=512 * 1024):
    """Useful GB/s of reading c-tile chunks every S tiles (c <= S; c = S is a sequential stream) over the array,
    kept up for ~target_s per measurement."""
    G = ARRAY_BYTES // (S * TILE)
    k = max(1, min(G, dma_bytes // (c * TILE)))
    G -= G % k
    per_dma = k * c * TILE
    x = array(G, S)
    run = make_stream(G, S, c, k, inflight, dmas_per_burst=G // k, spin_iters=0)
    guess = per_pass(run, 1, x, repeats=1)           # one pass over the array
    n = passes_for(guess, target_s)
    s = per_pass(run, n, x)
    del x
    nbytes = (G // k) * per_dma
    return {"gbps": nbytes / s / 1e9, "bytes_per_pass": nbytes, "passes": n, "dma_bytes": per_dma,
            "chunks": G, "pass_us": s * 1e6}


# %% [markdown]
# ## 1. Sequential stream kept up for ~1 ms .. ~200 ms

# %%
for target in (([0.001, 0.01]) if QUICK else [0.001, 0.01, 0.2]) if "seq_duration" in SECTIONS else []:
    guarded("seq_duration", {"target_s": target}, partial(stream_bw, 128, 128, target))


# %% [markdown]
# ## 2. Stride sweep: 4 KB chunks every S tiles

# %%
STRIDES = [1, 2, 4, 8] if QUICK else [1, 2, 3, 4, 5, 8, 9, 16, 17, 32, 33, 64, 65, 127, 128, 129, 256, 257, 512,
                                      1024, 1025, 2048, 4096, 8192, 16384]
for S in (STRIDES if "stride" in SECTIONS else []):
    guarded("stride", {"stride_tiles": S, "stride_bytes": S * TILE, "chunk_tiles": 1}, partial(stream_bw, S, 1))


# %% [markdown]
# ## 3. Chunk size at a fixed 1 MB stride

# %%
for c in (([1, 2]) if QUICK else [1, 2, 4, 8, 16, 32]) if "chunk" in SECTIONS else []:
    S = 4 if QUICK else 256                          # 1 MB
    guarded("chunk", {"chunk_tiles": c, "chunk_bytes": c * TILE, "stride_tiles": S}, partial(stream_bw, S, c))


# %% [markdown]
# ## 4. Bursts with HBM idle in between vs back-to-back

# %%
def burst_bw(burst_mb, idle_us=20.0, inflight=4):
    """Read burst_mb sequentially, then ~idle_us of VPU-only work; burst time = (burst + idle) - idle alone."""
    c, dma = 128, 512 * 1024
    G = ARRAY_BYTES // (c * TILE)                    # sequential: chunk = stride = 128 tiles
    k = max(1, min(G, dma // (c * TILE)))
    G -= G % k
    dmas = max(1, min(G // k, burst_mb * 2 ** 20 // (k * c * TILE)))
    x = array(G, c)
    calib_iters = 100 if QUICK else 10000
    spin_only = make_stream(G, c, c, k, inflight, dmas, calib_iters, with_dma=False)
    it_s = per_pass(spin_only, 4, x, repeats=3) / calib_iters   # seconds per spin iteration
    spin = int(min(1e7, max(1, idle_us * 1e-6 / max(it_s, 1e-12))))
    both = make_stream(G, c, c, k, inflight, dmas, spin)
    idle = make_stream(G, c, c, k, inflight, dmas, spin, with_dma=False)
    nbytes = dmas * k * c * TILE
    target = 0.05
    n = passes_for(per_pass(both, 1, x, repeats=1), target)
    t_both, t_idle = per_pass(both, n, x), per_pass(idle, n, x)
    stream = make_stream(G, c, c, k, inflight, dmas, 0)
    t_stream = per_pass(stream, passes_for(per_pass(stream, 1, x, repeats=1), target), x)
    del x
    return {"burst_gbps": nbytes / max(t_both - t_idle, 1e-12) / 1e9, "stream_gbps": nbytes / t_stream / 1e9,
            "idle_us": t_idle * 1e6, "burst_us": (t_both - t_idle) * 1e6, "bytes": nbytes, "spin_iters": spin}


for mb in (([4]) if QUICK else [4, 16, 64]) if "burst" in SECTIONS else []:
    guarded("burst", {"burst_mb": mb}, partial(burst_bw, mb))


# %% [markdown]
# ## 5. Stride sweep with plain XLA (fallback)

# %%
def xla_stride(S, c=1):
    G = ARRAY_BYTES // (S * TILE)
    x = array(G, S)

    @partial(jax.jit, static_argnums=0)
    def run(n, x):
        def body(i, acc):     # the max with the carry keeps the strided read inside the loop
            return acc + jnp.sum(jnp.maximum(x[:, :c], acc)) * 1e-12
        return lax.fori_loop(0, n, body, jnp.float32(0))
    n = passes_for(per_pass(run, 1, x, repeats=1), 0.03)
    s = per_pass(run, n, x)
    del x
    return {"gbps": G * c * TILE / s / 1e9, "pass_us": s * 1e6}


for S in (([1, 4]) if QUICK else [1, 2, 4, 8, 16, 32, 64, 128, 256, 1024, 4096, 16384]) if "xla_stride" in SECTIONS else []:
    guarded("xla_stride", {"stride_tiles": S}, partial(xla_stride, S))


# %% [markdown]
# ## Save and zip

# %%
save()
print("errors:", len(RESULTS["errors"]))
archive = shutil.make_archive(OUT, "zip", OUT)
print("wrote", archive)
try:
    from google.colab import files
    files.download(archive)
except Exception:
    pass
