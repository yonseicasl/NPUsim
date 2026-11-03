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

`tpu_v5e1_stream_bench.py` / `.ipynb` follows up the isolated Colab runs of `measure_tpu.py` (one GEMM or one decode
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

## Raw Colab archives

`results/2026-11-03_colab_v5e1_warm.zip` (one input buffer reused) and `results/2026-11-03_colab_v5e1_rot.zip`
(rotated buffers) are the archives the two Colab v5e-1 runs of `measure_tpu.py` produced, unchanged
(gemm/vector/decode CSVs, decode_ops.json, env.json; no traces kept). Branch `tpu-measure` has them unpacked
under `results/2026-11-03_colab_v5e1_{warm,rot}/` with their comparisons.

## Results: 2026-11-03 Colab v5e-1 follow-up (`results/2026-11-03_colab_v5e1_followup/`)

TPU v5 lite, JAX 0.11.2, no errors (raw archive `results/2026-11-03_colab_v5e1_followup.zip`).
`compare.log` / `compare_gemm.csv`: compare_v5e.py (branch tpu-measure) on this run's gemm.csv.

- **Streamed decode GEMMs** (Llama-3-8B projections, M 1..128, weights read from HBM inside one program):
  NPUsim -0.3..-8% (16 points, all under), e.g. 8x4096x4096 46.8 us measured vs 43.4. Weights stream at
  690-750 GB/s (85-92% of 819), which accounts for the gap. The isolated runs' -22..-35% (rotated: serial
  copy + fusion) and +74..+312% (warm: weight on chip) were measurement artifacts: streamed 8x4096x4096 is 46.8
  us vs 56.4 rotated (copy 44.8 + fusion 11.6) and 11.7 warm. Inside one program warm = cold (ratio ~1.0).
- **int8 weights**: 0.50-0.56x the bf16 time at M <= 32 (HBM-bound, half the bytes), 0.6x at M 128, ~1.0x at
  M 512 (compute-bound); dequantization is free at decode sizes.
- **Small / unaligned GEMMs** (odd_gemm, M 100..512, K = N 1000..2048): NPUsim -23..-36% at K=N 1000-1024,
  -10..-20% at 2000-2048, -12..+7% at 1100-1280. 1000 and 1024 measure the same (padding to 128 is right, and
  NPUsim agrees). The pattern points to a fixed per-GEMM cost (~1-2 us) NPUsim lacks and to its partition of
  output columns over the 4 MXUs in 128-column strips being coarser than XLA's (1024 -> 1100: NPUsim 1.6x,
  measured 1.27x). Not yet modeled.
- **decode_stack GEMM-only** per layer: 580 us at B 1 (NPUsim's 4 GEMMs -3.9%), 637-676 us at B 8..128, ~50 us
  above the sum of the four streamed GEMMs (588-606); NPUsim -10..-14% there. Remeasure with the fixed bench.
- **Not valid in this run** (bench artifacts, found in the saved HLO; fixed, rerun pending): attention and the
  full decode layer above ~128 MB of KV picked the cache out of a stacked array with a dynamic slice and
  relaid (B, L, kv, d) for the einsum -- two extra HBM copies, 3.6x the KV roofline instead of ~1.3-1.5x; and
  hbm_read summed a 1-D array (slow on TPU) and at 1-4 GB used one buffer, so XLA hoisted the read out of the
  loop. The bench now keeps KV as (B, kv, L, d), unrolls up to 8 distinct buffers in the loop body instead of
  slicing large ones, and reads >= 2 distinct 2-D arrays. Rerun only those sections with
  `%env TPU_BENCH_SECTIONS=decode_stack,attention,hbm_read` in a cell above the first code cell.

## Results: 2026-11-03 Colab v5e-1 follow-up, rerun (`results/2026-11-03_colab_v5e1_followup2/`)

Same chip, the fixed bench (KV as (B, kv, L, d), distinct KV/layer buffers unrolled in the loop body), all
sections; no errors (raw archive `results/2026-11-03_colab_v5e1_followup2.zip`).

- **Reproducible**: the 80 streamed GEMMs match the first run within 1.3% (median; max 6.8%), NPUsim MAPE 5.0% ->
  4.7%; the small odd GEMMs within 2.8% (median; a few us-scale cases up to 35%), NPUsim MAPE 16.1% -> 15.5%.
- **Attention is a KV read at the sustained HBM rate**: 1.03-1.12x the KV roofline at 819 GB/s (733-797 GB/s) for
  B >= 8, 1.16-1.7x at B 1 (short reads). The earlier 3.6x (stacked KV) and the isolated runs' 1.25-1.35x carried
  the cache relayout copy that a (B, L, kv, d) layout forces.
- **A decode layer is its GEMMs plus its attention**: full layer = GEMM-only layer + attention within 2.4% over
  all 12 (B, L). NPUsim's four projection GEMMs are -3..-13% of the GEMM-only layer (-10.6% at B 1, -5..-6% at
  B 8, -3..-6% at B 32, -11..-13% at B 128); the GEMM-only layer runs 1-12% over the four streamed GEMMs alone
  (residual adds, the gate/up split).
- **hbm_read is still not valid**: with the copies unrolled, each read no longer depended on the loop carry and
  XLA hoisted the loop-invariant sums out of the loop (0 or absurd GB/s). Fixed: each read takes a max with the
  carry; the compiled HLO now keeps every reduce inside the while body (it had none before). The Kaggle bench's
  HBM section had the same flaws (1-D, stacked, one buffer at >= 512 MB) and is fixed the same way. Rerun with
  `%env TPU_BENCH_SECTIONS=hbm_read` (a few minutes).

## Results: 2026-11-03 Colab v5e-1, third run (`results/2026-11-03_colab_v5e1_followup3/`)

`tpu_v5e1_stream_bench` with the hbm_read fix, all sections, no errors (raw archive
`results/2026-11-03_colab_v5e1_followup3.zip`).

- **Sustained HBM read: 732-757 GB/s** (64 MB 732, 256 MB 749, 1 GB 755, 4 GB 757), 89-92% of the 819 GB/s spec.
  The same rate shows up in every HBM-bound measurement: streamed decode weights 690-750 GB/s, attention's KV
  read 733-797 GB/s. This is the target for a MOSAIC HBM2e config (2 stacks x 409.6 GB/s) to reproduce.
- **Reproducible against the second run** (median / max |difference|): streamed GEMMs 0.9% / 5.8%, decode
  layers 0.7% / 5.1%, attention 0.3% / 28.8%, small odd GEMMs 6.3% / 70.9% -- the large maxima are us-scale
  cases (B 1 L 512 attention, M 100 K=N 1000-1024 GEMMs).

### What the three Colab runs establish for NPUsim (TPU v5e, one chip)

| quantity | measured | NPUsim |
|---|---|---|
| decode-shaped GEMMs, weights streamed (M 1..128) | HBM-bound at 690-750 GB/s | -0.3..-8% (MAPE ~5%), assumes 819 GB/s |
| compute-bound GEMMs | 187-192 of 197 TFLOPS | -1..+12% |
| int8 weights at M <= 32 | half the bf16 time | not compared yet (tpu_v5e int8 config) |
| decode attention, B >= 8 | KV bytes / ~745 GB/s (1.03-1.12x roofline) | not modeled in the PaLM flow |
| decode layer | GEMM-only layer + attention (within 2.4%) | 4 projection GEMMs -3..-13% of the GEMM-only layer |
| small / unaligned GEMMs | | -10..-36% at K=N 1000-1024 / 2000-2048: fixed per-GEMM cost, MXU strip partition |
