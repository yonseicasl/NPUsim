# TPU v5e-8 measurements (Kaggle)

Measured ground truth for NPUsim's TPU modeling (refs/tpu/README.md): single-chip GEMMs, HBM bandwidth,
8-chip collectives, and a tensor-parallel Llama-2-70B-shaped decode layer split into components. The
PaLM decode comparison (validation/palm/) leaves 41-60% of the measured latency unattributed; the decode
layer variants here measure where such time goes on real TPUs.

## Run on Kaggle

1. kaggle.com -> Create -> New Notebook. File -> Import Notebook -> upload `kaggle_tpu_bench.ipynb`
   (or paste `kaggle_tpu_bench.py` into one cell).
2. Settings (right panel) -> Accelerator -> **TPU v5e-8**. Internet is not needed.
3. Run All (about 30-45 minutes; Kaggle allows 20 TPU hours per week and 9 hours per session).
4. Output tab -> download `tpu_v5e_results.json` and the `trace/` directory.
5. Put them on this server under `validation/tpu_v5e/results/` (e.g. `results/2026-11-xx/`).

The first cell prints the JAX version and the devices; check that eight TPU v5 lite devices are listed.
A failing measurement is recorded under `errors` in the JSON and the run continues.

## What is measured

| section | content | timing |
|---|---|---|
| gemm | M in 1..512 x (K, N) in {4096^2, 8192^2, 8192x28672, 28672x8192}, bf16 and int8 weights (dequantized in the kernel); square bf16 GEMMs 1024..8192 | slope of n vs 2n dependent iterations in one jitted `fori_loop` |
| hbm | elementwise stream over 1 MB .. 4 GB | same |
| collective | all-reduce, all-gather, reduce-scatter over 8 chips, 4 KB .. 128 MB per chip | same, inside `shard_map` |
| decode_layer | Llama-2-70B layer (E 8192, F 28672, 64 heads, 8 KV heads, d_head 128), TP 8, B 1..128, context 2048, bf16 and int8 weights; variants full / no_comm / no_attn / no_vector / gemm_only; one profiler trace (B 8, bf16) | per layer, n vs 2n layers |

**Cold operands.** A loop that reads the same buffer every iteration measures on-chip reuse, not HBM:
TPU v5e keeps an operand of up to ~128 MiB on chip across iterations (and across executions), so decode
GEMMs with small weights ran 3-4x faster than HBM allows (Colab v5e-1, 2026-11-03). Every section above
except the collectives therefore cycles iteration i through copy i % R of R distinct buffers (weights,
KV cache, streamed array; a decode iteration runs layer i % R), R chosen so the copies total >= 512 MB
(`copies` in each result; at most ~1 GB per chip). `TPU_BENCH_ROTATE_MB=<MB>` changes the total;
`TPU_BENCH_ROTATE_MB=0` reuses one buffer (warm, the old behaviour) for comparison.

`make_notebook.py` regenerates the notebook from the script. `TPU_BENCH_QUICK=1 TPU_BENCH_OUT=<dir>` with
`JAX_PLATFORMS=cpu XLA_FLAGS=--xla_force_host_platform_device_count=8` runs a tiny CPU smoke test.
