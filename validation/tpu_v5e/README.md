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
