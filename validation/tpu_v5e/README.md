# NPUsim vs Cloud TPU v5e (measured)

The PaLM comparison (`validation/palm/`) leaves NPUsim 41-60% under the published TPU v4 decode latency,
and public data cannot say where the rest goes. This directory measures it on a real TPU: per-GEMM
device times to check `configs/accelerators/tpu_v5e.cfg`, and a decode layer broken down by named scope
to size what NPUsim does not model (attention, vector work, collectives).

| file | runs on | does |
|---|---|---|
| `measure_tpu.py` | TPU (JAX only, one file) | GEMM sweep, vector ops, collectives, Llama-3-8B decode layer at tp=1 and tp=8 -> CSVs |
| `compare_v5e.py` | here (needs `models/model`) | NPUsim `tpu_v5e` on the measured GEMM shapes; decode GEMM vs non-GEMM share |
| `test_parse_trace.py` | anywhere with JAX | checks the profiler-trace parser on a synthetic TPU trace |
| `kernel-metadata.json` | Kaggle CLI | private script kernel on `TpuV5E8` |

`tpu_v5e.cfg` is `tpu_v4.cfg` with 4 MXUs (one TensorCore), 1500 MHz (197 TFLOPS / (4 x 128^2 x 2)) and
HBM2e 819 GB/s. Nothing in it is fitted to TPU data yet.

## Measuring

Kaggle (free, 20 TPU h/week, TPU v5e-8, no v6e; refs/tpu/README.md section 4). The account needs phone
verification for TPUs and internet.

1. Web: New Notebook -> Settings -> Accelerator `TPU v5e-8`, Internet on. Upload `measure_tpu.py` (Add
   Input -> upload, or paste it into a cell) and run `!python measure_tpu.py` (quick check:
   `--quick`). Results land in `/kaggle/working/tpu_v5e/`; download them from the Output tab.
2. CLI: `pip install kaggle`, put the API token in `~/.kaggle/kaggle.json`, replace `KAGGLE_USERNAME` in
   `kernel-metadata.json`, then
   ```
   kaggle kernels push -p validation/tpu_v5e
   kaggle kernels status KAGGLE_USERNAME/npusim-tpu-v5e-measure
   kaggle kernels output KAGGLE_USERNAME/npusim-tpu-v5e-measure -p validation/tpu_v5e/results/<date>
   ```
   `machine_shape` is what selects the v5e-8 per a 2026-10 report; check the run's `env.json`
   (`device_kind`) in case the CLI ignored it.

TRC or a paid Cloud TPU VM (`v5litepod-1` for gemm/vector/decode tp=1, `v5litepod-8` for all):
`python3 measure_tpu.py --out results/<date>` on the VM (installs `jax[tpu]` if missing).

A full run takes about 10-20 minutes. Device times come from the profiler (`jax.profiler.trace`,
perfetto export, XLA module/op events on `TPU:0`); if the trace has no TPU process the CSVs keep wall
time only and `compare_v5e.py` says so. `--keep-traces` keeps the raw traces for TensorBoard.

## Comparing

```
LD_LIBRARY_PATH=<npusim libs> python3 validation/tpu_v5e/compare_v5e.py validation/tpu_v5e/results/<date>
```

~15 minutes on 16 threads. Prints measured vs NPUsim vs roofline per GEMM (MAPE split into decode-shaped
M <= 128 and the rest) and writes `compare_gemm.csv`; then per decode case the measured layer time, its
projection GEMMs, attention against its KV-read roofline, psum, and the non-GEMM share -- at tp=1 also
NPUsim's GEMM sum against the measured GEMM scopes. Commit the result CSVs under `results/<date>/`.

Not covered: int8 weights (PaLM's low-latency setup), v6e (256x256 MXU; Kaggle has none), multi-chip
NPUsim. XLA may fuse ops across named scopes (e.g. a norm into the next matmul); `decode_ops.json`
keeps the per-op times for checking the attribution.

## Results

### 2026-11-03, Colab TPU v5e-1, warm (`results/2026-11-03_colab_v5e1_warm/`)

First run, before buffer rotation: every call reused the same input buffers. JAX 0.11.2, `TPU v5 lite`.

- GEMMs that are compute-bound or stream more than the on-chip memory holds match NPUsim: M >= 512
  -6..+12% (MAPE 9%, 8192^3 0.0%), HBM-bound M = 1 and the 235 MB 4096x28672 weight -6..-2%. Large
  GEMMs reach 187-192 of 197 TFLOPS; HBM-bound ones ~750 GB/s (92% of 819).
- Decode-shaped GEMMs at M = 8..128 with <= 118 MB of weights ran at 2.7-3.3 TB/s, 3-4x what HBM allows
  (wall time too), so NPUsim reads +74..+312%. An unchanged input apparently stays on chip across executions
  (presumably XLA cross-program prefetch, not confirmed: `copy-start`/`copy-done` in the trace; the 235 MB weight and the non-dot
  M = 1 path do not benefit). A real decode step never reuses a layer's weights, so these points are not
  usable; `measure_tpu.py` now rotates distinct buffers (`--rotate-mb`, default 512; 0 reproduces this).
- Decode layer (Llama-3-8B, tp=1): at B = 1 (not affected) NPUsim's four GEMMs are -4..+2% of the
  measured scopes. Norms, SwiGLU and residual take 3-5 us per layer (~1%). Attention runs at 1.25-1.35x
  its KV-read roofline and dominates as B x L grows (49% of the layer at B 8 / L 8192, 93% at B 128 /
  L 8192). Single-chip vector work therefore does not explain the PaLM gap; collectives (tp=8, v5e-8)
  are the next measurement.

### 2026-11-03, Colab TPU v5e-1, rotated buffers (`results/2026-11-03_colab_v5e1_rot/`)

Same chip and JAX, `measure_tpu.py` with buffer rotation (>= 512 MB of distinct inputs; decode rotates
two whole layers incl. KV cache). `compare_v5e.py` now also reports the module time minus its
`copy-start`/`copy-done` ops ("kernel").

- With fresh buffers a GEMM whose operand fits on chip runs as a serial copy HBM -> on-chip memory
  followed by the fusion (33.5 MB in 44.8 us, 117 MB in 155.6 us: ~750 GB/s); the fusion time is the
  warm time. Operands that do not fit (235 MB) and the M = 1 path stream from HBM, also at ~750 GB/s.
  The sustained HBM read rate is therefore ~750 GB/s (92% of 819) across 14 independent points.
- Compute-bound GEMMs (M >= 512 or >= 1024^3): kernel -1..+12% (MAPE 9.6%). HBM-streaming GEMMs: -6..-2%
  (NPUsim assumes 819 GB/s). Decode-shaped GEMMs that fit on chip, isolated: module -22..-35%, because
  NPUsim overlaps the HBM stream with compute while XLA serializes copy and fusion in a one-op program.
- Decode layer (Llama-3-8B, tp=1), where the copies run inside the layer program: NPUsim's four GEMMs
  are -3..-9% of the measured GEMM scopes plus the unattributed ops (the copies), -4..-12% at B = 1.
  Attention 1.25-1.35x its KV-read roofline; norms, SwiGLU, residual 4-6 us (~1%).
- `tpu_v5e.cfg` keeps the specified 819 GB/s: the 8% gap to the sustained ~750 GB/s is DRAM behavior
  (refresh, bank and row effects) that the ideal `[dram]` model does not have, not a parameter to fit.
  It is the expected source of the -2..-9% HBM-bound errors; MOSAIC's Ramulator2 DRAM co-simulation is
  the way to model it. Open: attention's KV-read efficiency (~0.75-0.8 of roofline) in validation/palm.

### NPUsim on MOSAIC's Ramulator 2.1 (same rotated results)

`compare_v5e.py --dram validation/tpu_v5e/hbm2e_v5e.yaml` times NPUsim's DRAM tile streams on MOSAIC's
Ramulator2 (`NPUSIM_MOSAIC_DRAM`; the binary must be built with `MOSAIC_ROOT`). `hbm2e_v5e.yaml`: v5e's
HBM2e as 32 64-bit pseudo channels at 3.2 Gbps (819.2 GB/s), HBM2 timings in ns, nBL 1 as in MOSAIC's
HBM3 configs. MOSAIC's `dram_bw_probe` (sequential reads, 16 MB) sustains 674.8 GB/s on it, 82% of
nominal -- the same 82% as MOSAIC's calibrated HBM3 -- against the TPU's measured ~750 GB/s (92%).

Fixed on the way: with Ramulator the first weight strip's DRAM fill is exposed and partly hidden from
the second strip on, so extrapolating from 1 and 2 strips put 8192^3 under the compute roofline.
`compare_v5e.py` now extrapolates from 2 and 4 strips (analytic results move by < 3%).

| | analytic `[dram]` (819 GB/s, ideal) | Ramulator2 HBM2e (82% sustained) |
|---|---|---|
| compute-bound GEMMs, kernel | -0..+7% (2048^3 +15%) | same (DRAM hidden): -0..+7% |
| HBM-streaming GEMMs (M = 1, 235 MB weight; 10 points) | -6.4..-1.9% | +5.2..+10.3% |
| decode layer GEMMs (tp=1, 12 cases), mean abs | 5.8% (-12.4..-2.9%) | 6.3% (-1.3..+9.3%) |

The two DRAM models bracket the TPU: ideal streaming is ~8% fast, Ramulator's HBM2e ~8-10% slow. The
gap is the DRAM controller's sustained efficiency (82% vs 92%: refresh policy, address mapping, row
handling), which is DRAM modeling with the TPU's 750 GB/s as its evidence, not an NPUsim parameter.

Refresh. `hbm2e_v5e_norefresh.yaml` (NoRefresh) sustains 726.9 GB/s in `dram_bw_probe` (88.7%, vs 82.4%
with AllBank), yet NPUsim's results are identical to the cycle with either config: `ramulator_ns()`
(components/dram.cc) times each tile stream on a fresh Ramulator2 device, ~1 MiB long (one MXU's 4096 x
128 bf16 weight strip; smaller tiles are replicated to 1 MiB), i.e. ~1.5 us, shorter than the first
refresh at tREFI 3.9 us. Refresh never fires in NPUsim's co-simulation, so the Ramulator column above is
already the refresh-free bound. NPUsim streams at ~715 GB/s there (1 x 4096 x 8192: 67 MB in 93.8 us),
the probe's refresh-free 727 GB/s minus the per-stream start; the remaining gap to the TPU's ~750 GB/s
is Ramulator's row/command model and the fresh device per stream, not refresh. (HBM2_R20 has no
per-bank refresh manager: PerBank looks up nREFIpb, the 2.0 port names it nREFISB.)

### 2026-11-03, Colab TPU v5e-1, follow-up (`results/2026-11-03_colab_v5e1_followup/`)

Measured by a separate script (not `measure_tpu.py`; `results.json` keeps all of it, `hlo/` three HLO
dumps; traces not committed): GEMMs streamed in a `while` loop over distinct weights (bf16 and int8
weights), a stacked decode layer (`while` over layers, `full` and `gemm_only`), attention alone, odd-shaped
GEMMs, an HBM read test.

- Streamed decode-shaped GEMMs (M 1-128, the four Llama-3-8B projections; 20): no serial copy as in the
  isolated runs, 690-755 GB/s. NPUsim analytic mean -4.6% (MAPE 5.0%), Ramulator2 HBM2e +5.0% (6.0%):
  the same bracket as before, now on the cases that matter for decode. M 512: -5..+4% on both.
- Odd-shaped small GEMMs (M 100-512, K = N 1000-2048; 35): NPUsim short, analytic MAPE 15%, Ramulator 11%.
  Padding agrees (1100 and 1152 take the same time on both); a fixed ~1.5-3 us per kernel is missing
  (`layer_setup_cycle` 0, no first-DMA latency in `tpu_v5e.cfg`).
- Stacked decode `gemm_only`: B = 1 745-751 GB/s (analytic -4..-5%, Ramulator +7..+8%); B >= 8 645-685
  GB/s (analytic -10..-14%, Ramulator -2..+1%). The B 8 trace shows why: qkv/o_proj weights are first
  copied out of the stacked array (`dynamic-slice` fusions and `slice-done` waits, ~180 us per layer)
  and the dots then read the copies -- the cost of slicing stacked weights in a `while` loop, not of
  streaming them.
- Attention at >= 268 MB of KV runs at 226 GB/s (3.6x the KV roofline) while <= 67 MB runs at 550-610
  GB/s. GQA is not materialized (the einsum is grouped); the HLO slices the layer's KV out of the stack
  and relayouts it (`copy_bitcast_fusion`, layout {3,1,2,0}), presumably through HBM once it exceeds the
  on-chip memory. Implementation-dependent; not used for NPUsim.
- `hbm_read` is invalid (64 MB at 94 GB/s, 4 GB in 1 us: the read was elided or not timed on device).

### 2026-11-03, Colab TPU v5e-1, follow-up 2 (`results/2026-11-03_colab_v5e1_followup2/`)

The follow-up script re-run with the stacked artifacts removed: the decode stack now passes each layer
its own weight and KV buffers (no `while`/`dynamic_slice`) with head-major KV `[B, Hkv, L, D]`, and
attention runs unrolled. `compare_stack.py` compares the layers; `compare_v5e.py --accel tpu_v5e_int8w
--gemm gemm_int8.csv` the int8-weight GEMMs (`tpu_v5e_int8w.cfg`: int8 weights, bf16 activations).

| | analytic `[dram]` | Ramulator2 HBM2e |
|---|---|---|
| decode layer, full (GEMMs + KV at 819 GB/s; 12 B x L) | mean -8.4%, MAPE 8.4% | mean -0.5%, MAPE 3.8% |
| decode layer, gemm_only (12) | mean -8.1%, MAPE 8.1% | mean +3.5%, MAPE 3.9% |
| streamed decode GEMMs, bf16 (20; mean of both runs) | mean -4.3%, MAPE 4.8% | mean +5.2%, MAPE 6.1% |
| streamed GEMMs, int8 weights, M <= 128 (16) | MAPE 8.3% | MAPE 3.5% |
| streamed GEMMs, int8 weights, M = 512 (4) | MAPE 1.7% | MAPE 1.7% |
| odd small GEMMs (35; mean of both runs) | mean -15.1%, MAPE 15.3% | mean -7.4%, MAPE 9.9% |

- Attention alone (head-major KV) runs at 1.10x the KV read at 819 GB/s (~745 GB/s) from 67 MB of KV
  up, the HBM rate the weight streams sustain; small cases are latency-bound (x1.7 at 1 MB). The 1.25-1.35x
  in the rotated `measure_tpu.py` run and the 3.6x of the first follow-up came from the KV layout
  (`[B, L, Hkv, D]`) and the stacked slicing, not from the hardware.
- Run to run (the 55 GEMMs measured twice): streamed 1.4% mean (max 6.8%), odd small kernels 4.8% (max
  35%), which bounds what the small-kernel errors can say.
- int8 weights stream at 650-740 GB/s of weight bytes; at M = 512 the int8 GEMMs take the bf16 time
  (~180 TFLOPS), i.e. the MXU runs bf16 after dequantization, as `tpu_v5e_int8w.cfg` assumes.
- `hbm_read` is still invalid in this run (64 MB in 0.69 us); fixed in follow-up 3.

### 2026-11-03, Colab TPU v5e-1, follow-up 3 (`results/2026-11-03_colab_v5e1_followup3/`)

Third run of the same script with `hbm_read` fixed. `compare_gemm*.csv` here are follow-up 2's (the same
shapes: only their `npusim_us` is used, by `compare_stack.py`).

- HBM read: 64 MB 731.6 GB/s, 256 MB 749.4, 1 GB 754.9, 4 GB 756.9 GB/s -- 89-92% of the specified 819,
  the ~750 GB/s the GEMM and attention streams sustain, now measured directly. Same yardstick for the DRAM
  models: Ramulator2 HBM2e probe 674.8 GB/s (AllBank refresh) / 726.9 (no refresh), NPUsim's ~1 MiB
  Ramulator streams ~715; the analytical `[dram]` 819.
- Decode layers: full NPUsim+Ramulator2 mean -0.7% / MAPE 3.5% (analytic -8.6% / 8.6%), gemm_only +3.6% /
  3.6% (-8.0% / 8.0%).
- Over the three runs (55 GEMMs): streamed decode GEMMs analytic mean -4.3% / MAPE 4.9%, Ramulator2 +5.2%
  / 6.2%, run-to-run CV 0.8% (max 2.8%); odd small GEMMs analytic -14.6% / 15.1%, Ramulator2 -6.9% / 9.7%,
  CV 4.3% (max 21%). Decode layers differ by 1.0% (max 5.1%) between runs 2 and 3, streamed GEMMs 1.2%,
  attention 4.6% (large cases at 1.10-1.12x the KV roofline again).
