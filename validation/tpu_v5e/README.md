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
| hbm | read-only stream over 1 MB .. 4 GB | same |
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

# TPU v5e-1 follow-up (Colab)

`colab_v5e1_bench.py` / `.ipynb` follows up the isolated Colab runs of `measure_tpu.py` (one GEMM or one decode
layer per XLA program; `results/2026-11-03_colab_v5e1_*` on branch `tpu-measure`). With fresh buffers such a
program copies an operand HBM -> on-chip memory and then runs the fusion, serially, so it overstates an HBM-bound
GEMM; a decode step streams many different weights inside one program. Here everything runs inside one program
(n vs 2n dependent iterations of a `fori_loop`, R distinct operand copies, >= 512 MB):

| section | measures | answers |
|---|---|---|
| stream_gemm | Llama-3-8B projections (M 1..512) streaming distinct weights (cold) / one weight (warm), bf16 and int8 | HBM-bound GEMM time with copy/compute overlap; on-chip reuse factor; int8 dequant cost |
| decode_stack | R distinct Llama-3-8B layers (weights + KV) in one loop, B 1..128, L 512..8192, full / GEMM-only | layer time without the one-op-program copy; non-GEMM share |
| attention | decode attention alone over a rotated KV cache, B 1..128, L 512..32768 | KV-read efficiency vs roofline (open in validation/palm) |
| hbm_read | read-only streams 64 MB .. 4 GB | sustained HBM read bandwidth (to check a MOSAIC HBM2e config) |
| odd_gemm | M 100..512 x K=N 1000..2048 around multiples of 128 | padding / tiling granularity |

For two cases (8x4096x4096 streamed and as a one-op program; decode B 8 L 2048) it also keeps the optimized HLO
(`hlo/`), XLA's cost analysis and a profiler trace (`traces/`, for TensorBoard's profile plugin: trace viewer,
op_profile, hlo_stats, roofline, memory_profile), and the run's HBM memory stats.

Colab: Runtime -> Change runtime type -> v5e-1 TPU; upload the notebook (File -> Upload notebook) and Run All
(~20-30 min). It downloads `tpu_v5e1_followup.zip`; unzip it under `results/<date>_colab_v5e1_followup/`, then

```
LD_LIBRARY_PATH=<npusim libs> python3 validation/tpu_v5e/compare_v5e.py results/<dir>     # tpu-measure branch
python3 validation/tpu_v5e/colab_summary.py results/<dir> [--isolated results/2026-11-03_colab_v5e1_rot]
```

`gemm.csv` there is in measure_tpu.py's format (device_us = loop slope), so compare_v5e.py takes it as is and
writes `compare_gemm.csv`, which colab_summary.py then uses for the NPUsim columns.
