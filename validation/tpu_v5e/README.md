# NPUsim vs Cloud TPU v5e (measured)

The PaLM comparison (`validation/palm/`) leaves NPUsim 41-60% under the published TPU v4 decode latency,
and public data cannot say where the rest goes. This directory measures it on a real TPU: per-GEMM
device times to check `configs/accelerators/tpu_v5e.cfg`, and a decode layer broken down by named scope
to size what NPUsim does not model (attention, vector work, collectives).

| file | runs on | does |
|---|---|---|
| `measure_tpu_v5e1.py` | one TPU chip (JAX only, one file) | GEMM sweep, vector ops, Llama-3-8B decode layer, `smallk`, `loopk` -> CSVs; the shared measurement, trace and output code |
| `measure_tpu_v5e8.py` | 2+ chips (JAX only, one self-contained file) | collectives over 2/4/8 chips, the decode layer at tp 1/2/4/8; also takes the single-chip experiments with `--only` |
| `compare_v5e.py` | here (needs `models/model`) | NPUsim `tpu_v5e` on the measured GEMM shapes; decode GEMM vs non-GEMM share |
| `compare_stack.py` | here | decode layers and attention of a follow-up `results.json` vs NPUsim |
| `test_parse_trace.py` | anywhere with JAX | checks the profiler-trace parser and `linfit` on synthetic data, and that the two scripts' shared functions are identical |

The single-chip script was `measure_tpu.py` until 2026-11-06 (renamed when the multi-chip experiments moved
out); the results below name it as it was at the time.

`tpu_v5e.cfg` is `tpu_v4.cfg` with 4 MXUs (one TensorCore), 1500 MHz (197 TFLOPS / (4 x 128^2 x 2)) and
HBM2e 819 GB/s; the TPU-measured values in it are listed in its header.

## Measuring

One chip (Colab `v5e-1 TPU`, Kaggle, or a `v5litepod-1` VM): paste `measure_tpu_v5e1.py` into a
`%%writefile /content/measure_tpu_v5e1.py` cell and run `%run /content/measure_tpu_v5e1.py --out <dir>`
(default sections gemm, vector, decode; `--only smallk`, `--only loopk` for the fixed-cost experiments).
Use `%run`, not `!python`: the notebook kernel holds the TPU once JAX is imported.

Eight chips (Kaggle `TPU v5e-8`, free 20 TPU h/week with queues of hours; refs/tpu/README.md section 4; or a
`v5litepod-8` VM): `measure_tpu_v5e8.py` alone, which carries the same measurement code as the single-chip
file (a copy, kept identical by `test_parse_trace.py`):

```
%%writefile /kaggle/working/measure_tpu_v5e8.py      (then the file's contents)
%run /kaggle/working/measure_tpu_v5e8.py              (collectives,decode -> /kaggle/working/tpu_v5e8/)
%run /kaggle/working/measure_tpu_v5e8.py --only loopk # a single-chip experiment not measured yet
```

On a VM: `python3 measure_tpu_v5e8.py --out results/<date>` (installs `jax[tpu]` if missing). Change the shared
code in both files.

Fixed per-kernel time: `measure_tpu_v5e1.py --only smallk` (a few minutes on one chip) writes `smallk.csv` --
K sweeps at M 8/128, N 128/1024 (the K = 0 intercept is the fixed cost, `device_us` with the program launch,
`ops_us` the kernels alone), `x + 1` from one vreg to 4 MB, and chains of 1-32 dependent matmuls in one
program (slope = cost of one more kernel inside a program, intercept = the program's own) -- and the fits
in `smallk_fit.csv`.

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

### 2026-11-04, Colab TPU v5e-1, fixed per-kernel time (`results/2026-11-04_colab_v5e1_smallk/`)

`measure_tpu_v5e1.py --only smallk`; `gemm_smallk.csv` is its K sweep for `compare_v5e.py`.

| | TPU | NPUsim analytic | NPUsim Ramulator2 |
|---|---|---|---|
| K sweep intercept, isolated kernel (M 8/128, N 128/1024) | 0.95-1.19 us | 0.11-0.76 us | 0.11-0.78 us |
| K sweep slope | 543-653 GB/s | 408 (N 128 at M 8), 771-808 GB/s | 408, 715-717 GB/s |
| smallest kernel (`x + 1`, one vreg) | 0.57 us | - | - |
| chain of n matmuls in one program | 0.69 us + 0.22 us per kernel | - | - |
| K sweep MAPE (32) | | 38.1% | 33.6% |

- An isolated kernel pays ~1 us that does not scale with bytes; module time equals the op sum (launch is
  outside the device time), so it is inside the kernel: the first operand's HBM latency and DMA start,
  the array fill/drain and the result write-back, not overlapped with anything. In one program the next
  kernel's operands are fetched during the current one: the chain costs 0.69 us once plus 0.22 us per
  kernel (8 x 128 x 128), and 2 MB-weight kernels chain at 757 GB/s (2.77 us each). The decomposition is
  inferred from the fits, not measured part by part. In a decode layer (~600 us, a few dozen kernels) the
  per-kernel part is ~0.1%; the large isolated-kernel MAPE is the one-kernel-program case.
- NPUsim models every GEMM as an isolated layer with its own fill but no per-kernel control cost
  (`layer_setup_cycle` 0 at the time). Applied since: 0.22 us = 330 cycles at 1.5 GHz (below).
- Narrow N: at M 8, N 128 NPUsim runs one MXU and its 128-element-per-cycle port (384 GB/s) limits it to
  408 GB/s, while the TPU streams at 653 GB/s; how the TPU spreads such a GEMM is not known.
- `compare_v5e.py` fix: N with fewer than 4 128-column blocks ran padded to 4 MXUs (N 128 simulated as
  512, 4x the weights); now only nb MXUs work (8 x 4096 x 128: 2.85 us vs 2.75 measured). In the rotated
  run 256^3 moves from +30% to -4% (kernel), 128^3 from +4% to -49% (its missing fixed time now shows). A
  partial last round of strips still counts whole: charging its column share was tried and put the K = N
  1100-1280 GEMMs at -20..-32% instead of -2..-20%.

### `layer_setup_cycle = 330` in `tpu_v5e.cfg` / `tpu_v5e_int8w.cfg` (2026-11-04)

The per-kernel cost inside one program from the smallk chain (0.22 us at 1500 MHz), the one TPU-measured
value in the v5e configs; everything else stays specified. It is added once per layer to the compute
schedule and overlaps the DRAM stream, so it shows where compute is the critical path (8 x 128 x 128:
+330 cycles) and hides behind HBM-bound streams (8 x 4096 x 512: +10 cycles). All comparisons above were
re-run with it (mean / MAPE, analytic | Ramulator2):

| | before | after |
|---|---|---|
| smallk K sweep (32) | -33.7 / 38.1% \| -29.2 / 33.6% | -28.9 / 33.2% \| -24.9 / 29.3% |
| rotated small (128^3, 256^3, odd; kernel) | +21.1 / 39.1% \| +31.4 / 47.0% | +28.4 / 35.2% \| +38.0 / 44.2% |
| odd small GEMMs, follow-up 2 (35) | -14.4 / 15.5% \| -6.6 / 10.2% | -13.9 / 15.2% \| -6.3 / 10.1% |
| streamed decode GEMMs bf16 (20) | -4.1 / 4.7% \| +5.5 / 6.2% | unchanged |
| int8 streamed, M <= 128 (16) | -8.3 / 8.3% \| +1.9 / 3.5% | unchanged |
| compute-bound GEMMs, rotated (kernel) | +5.9 / 7.5% \| +8.3 / 8.7% | +6.3 / 7.5% \| +8.5 / 8.6% |
| decode layer full, follow-up 3 (12) | -8.6 / 8.6% \| -0.7 / 3.5% | unchanged |

### Small GEMMs: program setup, finer column split, K split (2026-11-04)

Three changes, each from its own evidence (mean / MAPE, analytic | Ramulator2; "before" = the previous
commit):

1. `program_setup_cycle = 1035` (new `[accelerator]` knob in NPUsim, `scheduler/npu.cc`: added once per run
   to the network totals, default 0): the smallk chain intercept, 0.69 us at 1500 MHz. One NPUsim run is one
   program; `compare_v5e.py --in-program` zeroes it for kernels measured inside a loop or chain.
2. A partial last round of 128-column strips is charged its column share (XLA tiles finer: K = N
   1024 -> 1100 measured 1.27x, column share 1.22x, whole rounds 1.6x).
3. Fewer than 4 column blocks: the idle MXUs split K (partial sums), so all 4 stream weights (one MXU's
   128-element port capped N = 128 at 384 GB/s; the TPU streams 621 GB/s there).

| (measurement context) | before | after |
|---|---|---|
| smallk K sweep, isolated (32) | -28.9 / 33.2% \\| -24.9 / 29.3% | -8.8 / 16.1% \\| -2.9 / 12.2% |
| smallk2 K sweep, isolated, rerun (32) | | -8.9 / 16.2% \\| -3.1 / 12.2% |
| rotated small GEMMs, isolated (6, module) | -20.2 / 23.8% \\| -14.6 / 22.0% | -1.5 / 16.6% \\| +3.9 / 12.9% |
| compute-bound GEMMs, isolated (kernel) | +6.3 / 7.5% \\| +8.5 / 8.6% | +7.2 / 7.2% \\| +9.0 / 9.0% |
| odd small GEMMs in a `while` loop, follow-up 2 (35) | -13.9 / 15.2% \\| -6.3 / 10.1% | -22.2 / 23.3% \\| -15.2 / 17.0% |
| odd small GEMMs chained in one program, per kernel (35) | | +33.6 / 33.6% \\| +45.7 / 45.7% |
| streamed decode GEMMs, int8 GEMMs, decode layers | | unchanged |

The smallk2 run (`results/2026-11-04_colab_v5e1_smallk2/`, raw `.zip` beside it) repeats smallk within
1% and adds `copy_us` and the odd-shape chains (`gemm_oddchain.csv` = their per-kernel slopes):

- Isolated kernels are a serial copy plus a kernel: XLA stages the weights into on-chip memory (0.4 us +
  bytes / ~750 GB/s, the sustained HBM read rate) and only then runs the kernel. That is the 543-653 GB/s
  of the isolated K sweep: no lower HBM bandwidth, the copy does not overlap the compute.
- Chained in one program, an odd small GEMM costs its weight bytes at ~712 GB/s: its start overlaps the
  previous kernel. NPUsim charges every GEMM its own fill, drain and first DMA, ~1-1.3 us too much each
  (+34%). Inside a `while` loop (follow-up 2) the same shapes take 1-3 us more per iteration than in the
  chain, so NPUsim is short there. Open: an inter-layer overlap model (layer i+1's prologue under layer
  i's tail) for chained kernels, the per-iteration cost of loops, XLA's serial staging for isolated ones.

### Inter-layer overlap, copy-aware comparison, loop measurement (2026-11-04)

- NPUsim `[accelerator] interlayer_overlap` (default 0; 1 in `tpu_v5e*.cfg`): each layer records its prologue
  (until the last stage starts its first tile) and tail (after the first stage, DRAM, finishes) from the
  per-tile pipeline (`pipeline_timeline_cycles` now returns both); the network rollup hides
  min(tail_i, prologue_{i+1}) at every boundary -- the next layer's weights stream while the previous one
  drains, its compute still starts after it. Layers with SFU / decompression / KV / attention extensions keep
  0. Checked on chained 128 x 1024 x 1024 networks (1131 cycles hidden per boundary) and in
  `unittest/validation_test.cc`; `run_validation.sh` and `run_stats_timeline_validation.sh` pass.
- `compare_v5e.py --chain` compares per-kernel slopes of a chain with NPUsim's GEMM minus that overlap
  (implies `--in-program`). Odd small GEMMs chained in one program (35, smallk2):

  | | no overlap | with overlap |
  |---|---|---|
  | analytic, all | +33.6 / 33.6% | +11.3 / 11.9% |
  | analytic, M 512 | +26.4 / 26.4% | -0.3 / 2.7% |
  | Ramulator2, all | +45.7 / 45.7% | +23.2 / 23.4% |

  Ramulator2 stays higher: each NPUsim tile stream starts a fresh device and pays its first access.
- XLA's serial staging is kept out of NPUsim: where a module's copy time is recorded (`ops` copy entries, or a
  `copy_us` column as in smallk2's `gemm_smallk.csv`), `compare_v5e.py` also reports kernel = measured - copy.
  That kernel time is a VMEM-resident run, so it is meaningful against NPUsim for compute-bound isolated GEMMs
  only; memory-bound isolated GEMMs are dominated by the copy and are not used as fidelity points.
- `compare_v5e.py --in-program` read and rewrote the config copy in one expression, truncating it first: every
  in-program comparison in `c4a9b06` failed and kept the previous CSVs (same conditions: program setup 0).
  Fixed; all of them re-run.
- `measure_tpu_v5e1.py --only loopk` (not yet run on a TPU): n = 1, 4, 16 dependent matmuls as a `fori_loop` over
  stacked weights and unrolled, M x K = 8 x 128, 100 x 1000, 128 x 1024, 256 x 2048, 8 x 4096; the
  loop - unrolled slope is the per-iteration cost of a while loop (`loopk_fit.csv`).

### Small and odd GEMMs: activation residency, staging what-if, tile granularity (2026-11-05)

- **Activation residency** (NPUsim `[accelerator] activation_residency`, default 0; 1 in `tpu_v5e*.cfg`): in the
  network path, the activation between mapped layers stays in the GLB when it fits its tensor residency
  capacity (`apply_graph_residency`, as the executable-IR path already did); 2 also keeps the run's own input and
  output, for a kernel in the middle of a chain (`compare_v5e.py --chain` sets it). Evidence: a chained odd GEMM
  costs its weight bytes only (100 x 1000: 2 MB in 2.81 us, 712 GB/s), so its input and output do not touch HBM.
  Odd GEMMs chained in one program (35):

  | | interlayer overlap only | + activation residency |
  |---|---|---|
  | analytic, all | +11.3 / 11.9% | -5.5 / 6.4% |
  | analytic, M <= 256 | +14.2 / 14.2% | -6.7 / 7.3% |
  | Ramulator2, all | +23.2 / 23.4% | +2.0 / 4.1% |
  | M 512 (both) | 2.1-2.7% | 2.5-2.6% |

  Chained kernels now bracket the TPU like the large decode GEMMs do (analytic ~5% short, Ramulator2 ~2% long).
- **Isolated kernels, staging not overlapped** (what-if `configs/accelerators/tpu_v5e_serial.cfg`: `[multi_chip]
  double_buffer = 0`, no outstanding-request override, so DRAM and on-chip stages add per tile). XLA copies an
  on-chip-sized operand before an isolated kernel (smallk2 `copy_us`: 0.4 us + bytes / ~755 GB/s):

  | (mean / MAPE, analytic) | overlapped (`tpu_v5e`) | not overlapped (`tpu_v5e_serial`) |
  |---|---|---|
  | smallk2 K sweep, copy >= 50% of the module (18) | -16.0 / 18.8% | +10.4 / 10.4% |
  | smallk2 K sweep, copy < 50% (14) | +0.1 / 13.0% | +13.4 / 13.4% |
  | smallk2 K sweep, all (Ramulator2) | -1.9 / 12.2% | +17.6 / 17.6% |
  | rotated, copied decode shapes M 8-128 (10) | -25.0 / 25.0% | +11.3 / 11.9% |
  | rotated, HBM-streaming, no copy (8) | -1.9 / 4.7% | +39.8 / 39.8% |
  | rotated, compute-bound (15) | -9.3 / 11.6% | -6.2 / 17.3% |

  Serial staging fits only where XLA stages (isolated program, operand on chip) and then overshoots by ~10%:
  the TPU's VMEM-resident kernel still overlaps its on-chip stages, which this config serializes too. Not
  applied; a conditional staging mode (isolated program and operand <= on-chip memory) would be XLA behavior
  in NPUsim, left for a decision.
- **Why small N = 1024 GEMMs read long in isolation** (K sweep K <= 256: +5..+26%). Not a per-strip fill: NPUsim
  declares none (`Stripe-transition: not modeled`; fold fill = the 330 setup cycles). The cause is tile
  granularity: one DRAM tile per strip (8 x 128 x 1024: two 128 KB tiles), so the five stages barely overlap
  (busy DRAM 522 / NoP 274 / GLB 128 / PE array 256 / PE 346 cycles, critical path 1526, prologue 1180).
  Moving part of K to the DRAM loop makes it worse, not finer (8 x 2048 x 1024: prologue 0.76 / 1.53 / 3.06 /
  6.11 us at 1 / 2 / 4 / 8 DRAM K splits; the split K serializes through the partial sums), and K <= 256 has no
  K to split. A finer DMA granularity inside a DRAM tile (the next stage starting on the first chunk) is a
  pipeline-model change, not a mapping one; not done.
- **v5e-8 code** (now `measure_tpu_v5e8.py`): collectives now sweep 2, 4, 8 chips (rotated buffers) and decode runs tensor
  parallel 1, 2, 4, 8; `env.json` records each chip's torus coordinates. For a Kaggle TpuV5E8 session:
  `%run measure_tpu_v5e8.py` (collectives and decode; the single-chip sections are measured already).

## Fidelity improvement plan (2026-11-11)

Pooled MAPE over the validated sets (102 points: decode layers, streamed decode GEMMs, int8 GEMMs, chained
odd GEMMs, compute-bound GEMMs): analytic 6.5%, Ramulator2 `hbm2e_v5e.yaml` 5.0%. With the out-of-scope sets
(isolated small GEMMs, GEMMs in a `while` loop; 169 points): 11.8% / 8.9%. What is left, by size:

| # | error source | where it shows | fix | needs |
|---|---|---|---|---|
| 1 | sustained HBM rate: TPU 92% of 819 GB/s, analytic 100%, Ramulator2 (`hbm2e_v5e.yaml`, fresh device per 1 MiB stream) ~87% | decode layers, streamed and int8 GEMMs, attention: analytic -5..-8%, Ramulator2 ~+5% | a DRAM model that sustains the measured rate for the right reason | no new measurement |
| 2 | compute-bound GEMMs read +7..9% (2048^3 +15%); the TPU reaches 95-97% of peak | large M | find the term in NPUsim's stage breakdown (row-tile fill, padding) | no new measurement |
| 3 | coarse DRAM tiles: stages cannot overlap inside a tile | isolated small N = 1024 GEMMs +5..+26% | pipeline model: the next stage starts on a tile's first chunk | no new measurement |
| 4 | per-iteration cost of a `while` loop | GEMMs in a loop -22% | model from `loopk` | `measure_tpu_v5e1.py --only loopk` |
| 5 | inter-chip collectives | PaLM multi-chip -41..-60% | ICI latency/bandwidth model | `measure_tpu_v5e8.py` on a v5e-8 |

Item 1, DRAM side, exists on branch `mosaic` (parallel session; MOSAIC `npusim-ddr3`, Ramulator2 `378827a`):
`configs/architecture/dram/hbm2e_tpu_v5e.yaml` (Ramulator 2.1 native HBM2, public specs: 698 GB/s with
all-bank refresh, 748 without) and `hbm2e_tpu_v5e_refsb.yaml` (`rw_beside_bank_refresh` + `IdleBankRefresh`:
per-bank refresh hidden behind reads of other banks, 751 GB/s = 91.7%, inside the measured 732-757); decode
GEMM MAPE 5.4% analytic -> 3.2% there. This branch does not re-implement it: item 1 here is to run this
branch's comparison sets on that config (below). The two branches have diverged (25 / 23 commits) and both
edit `tpu_v5e.cfg` and this README; merging is left to the owner.

### Item 1: the validated sets on `hbm2e_tpu_v5e_refsb` (2026-11-11)

`compare_v5e.py --dram <MOSAIC>/configs/architecture/dram/hbm2e_tpu_v5e_refsb.yaml` (per-bank refresh hidden
behind reads of other banks; MOSAIC `npusim-ddr3`, Ramulator2 `378827a`) on this branch's sets; outputs
`compare_*hbm2e_tpu_v5e_refsb*`. One 8 x 4096 x 2048 GEMM streams at 775 GB/s analytic, 690 on `hbm2e_v5e.yaml`,
706 on `hbm2e_tpu_v5e.yaml` (all-bank refresh), 715 on refsb. Mean / MAPE:

| set | analytic `[dram]` | `hbm2e_v5e.yaml` | `hbm2e_tpu_v5e_refsb.yaml` |
|---|---|---|---|
| decode layer, full (12) | -8.6 / 8.6% | -0.7 / 3.5% | -3.2 / 4.0% |
| decode layer, gemm_only (12) | -8.0 / 8.0% | +3.6 / 3.6% | -0.0 / 2.4% |
| streamed decode GEMMs bf16 (20) | -4.1 / 4.7% | +5.5 / 6.2% | +2.5 / 3.2% |
| int8 streamed GEMMs (20) | -6.9 / 7.0% | +1.3 / 3.2% | -1.6 / 3.1% |
| chained odd GEMMs (35) | -5.5 / 6.4% | +2.0 / 4.1% | -0.5 / 3.1% |
| compute-bound isolated GEMMs, kernel (15) | +7.2 / 7.2% | +9.0 / 9.0% | +8.2 / 8.2% |
| **pooled (102)** | **-4.0 / 6.5%** | **+3.3 / 5.0%** | **+0.8 / 4.0%** |

- The HBM-bound sets move to 3.1-3.2% with a mean within 2.5% of zero; without the compute-bound set the
  pooled MAPE of the other 87 points is 3.3%.
- The full decode layer is -3.2% while its GEMMs are -0.0%: `compare_stack.py` adds the KV read at the
  specified 819 GB/s (NPUsim keeps KV streams analytical on MOSAIC DRAM too), ~8% short of the measured rate.
  `hbm2e_v5e.yaml`'s 3.5% was its +3.6% GEMMs cancelling that. Putting the KV stream on the DRAM model is the
  fix (the attention consumer on branch `mosaic` does it).
- Compute-bound GEMMs do not depend on the DRAM model (+7..9%): item 2, now the largest term.

### Item 2: the layer timeline's tile clock and the DRAM-level loop (2026-11-11, simulator change)

Cause of the compute-bound +7..9%, found in NPUsim's own stage report (no TPU number involved). A 512 x 2048
tile at 1 / 2 / 4 weight strips (512 columns each, one DRAM iteration per strip), cycles:

| strips | compute schedule | roofline | critical path | prologue |
|---|---|---|---|---|
| 1 | 8522 | 8192 | 10581 | 1862 |
| 2 | 16714 | 16384 | 20067 | 3353 |
| 4 | 33098 | 32768 | 39432 | 6334 |

The compute schedule is the roofline plus the 330 setup cycles; the excess is the prologue, which is exactly
(DRAM + NoP + GLB + PE-array stage totals) / 16. The layer timeline's tile clock was the GLB repetition count
alone (16 = K / 128 folds), while the DRAM-level loop is simulated iteration by iteration and its stage totals are
summed -- so with D strips one timeline tile spanned D strips and the pipeline fill was charged D times.

`[accelerator] timeline_dram_tiles` (default 0: the previous clock, results unchanged; 1 in `tpu_v5e.cfg`,
`tpu_v5e_int8w.cfg`, `tpu_v5e_serial.cfg`): tile clock = GLB repetitions x the DRAM mapping row's iterations
(temporal factors only; chips are other rows); the off-chip and output-commit tiles scale with it; stage totals
are untouched. Same tile with the knob on: 10581 / 18390 / 34682 (prologue 1862 / 1676 / 1584). Unit check
`run_outer_iterations` in `unittest/stats_timeline_test.cc` (a 100-cycle DRAM stage before a 200-cycle compute
stage: 250 -> 212.5 cycles, derived before running); `run_validation.sh` (617 configs) and
`run_stats_timeline_validation.sh` pass with the knob off.

`compare_v5e.py` with it: constants of a run are charged once -- `program_setup_cycle` once per program (zeroed
in every simulated piece, added once unless `--in-program`), `layer_setup_cycle` by the first row tile only (the
other row tiles run on `<accel>_inner` with it 0); before, a 2048-row GEMM paid both four times. `err_kernel`
keeps the program setup in the prediction (a row without copies has one prediction for module and kernel);
`npusim_noprog_us` / `err_kernel_noprog_pct` are a sensitivity only. `scale_time()` is checked by
`test_compare_scale.py` on mocked times; `check_tiling.py` compares the tiled estimate with a direct run of the
whole GEMM.

The TPU comparison numbers of every set with this change are a separate entry below once the re-runs finish;
`compare_*serial*` and `compare_*norefresh*` (what-ifs) predate it and were not re-run.

int8 weights with the knob (Codex's representative regression, `results/2026-11-11_codex_audit/core_representative/`,
refsb, M 1 / 8, K = N = 4096): bf16 +2.9 / +3.9% -> -0.5 / +0.5%, int8 -5.8 / -4.0% -> -8.9 / -7.2%. The knob stays on
in `tpu_v5e_int8w.cfg`: it corrects the timeline, which does not depend on the weight format, and the int8 deficit
is in the measurement, not in the knob. (Superseded 2026-11-12: this paragraph went on to say that int8 weights
"stream 3-8% slower per byte"; that was a ratio of total times. The fit below finds the same per-byte rate for both
dtypes and a fixed part per loop iteration -- see "int8 weights: what the error is".) Before the knob the per-strip
fill covered part of it.

**Unresolved: row tiling in `compare_v5e.py` (found by Codex, 2026-11-11).** For M > 512 the script simulates
one <= 512-row tile as its own run and sums the tiles. A direct run of the whole GEMM (`check_tiling.py`) gives
1024 x 128 x 1024: tiled 6.63 us vs whole 5.45 us, **+21.7%** (`results/2026-11-11_codex_audit/whole_gemm_small.log`).
The earlier "1024 rows = 2 x 512 to 1.5%" was one point (1024 x 4096 x 512, before `timeline_dram_tiles`) and does
not hold in general: each simulated tile pays its own pipeline fill and drain, so the error grows as the compute
per tile shrinks. This is the comparison script's approximation, not the simulator's result -- the same kind of
per-iteration fill the knob removed inside one run, reintroduced by splitting the run. Not fixed here (no overlap
subtracted, nothing fitted); the cause at the row boundary (fill, input reuse) needs its own analysis, and the
clean route is to simulate M > 512 GEMMs whole (DRAM-level batch factor), which costs simulation time and failed
once at 4 row tiles. Affected: the 10 isolated GEMMs with M > 512 in `2026-11-03_colab_v5e1_rot` (10 of the 15
compute-bound points; the script prints a warning listing them). Not affected (M <= 512): streamed decode GEMMs,
int8 GEMMs, chained odd GEMMs, the K sweeps, the decode layers.

`check_tiling.py`, default shapes, on commit `7044920`'s binary (`results/2026-11-11_check_tiling/`):

| shape | tiled estimate | whole GEMM | difference | what it tests |
|---|---|---|---|---|
| 512 x 2048 x 4096 | 44.84 us | 45.06 us | -0.47% | strip extrapolation (8 strips from 2 and 4) |
| 512 x 4096 x 4096 | 88.86 us | 89.06 us | -0.22% | strip extrapolation |
| 128 x 4096 x 4096 | 44.41 us | 44.42 us | -0.02% | strip extrapolation |
| 1024 x 1024 x 1024 | 13.54 us | 14.17 us | -4.44% | two row tiles |

The strip extrapolation (every M <= 512 set) agrees with the whole run within 0.5%, slightly low. Row tiling is
off in both directions (-4.4% here, +21.7% at 1024 x 128 x 1024): unvalidated, as above.

### Every set with `timeline_dram_tiles` on (2026-11-11 re-run; binary of commit `7044920`)

Mean / MAPE, before (commit `cb08bea`) -> after. The decode layers are `compare_stack.py` on follow-up 3's
measurements with follow-up 2's re-run GEMM CSVs and the KV read at 819 GB/s (run in a scratch directory;
`results/2026-11-03_colab_v5e1_followup3/compare_*` are the pre-change copies and were not rewritten).

| set | analytic `[dram]` | `hbm2e_v5e.yaml` | `hbm2e_tpu_v5e_refsb.yaml` |
|---|---|---|---|
| streamed decode GEMMs bf16 (20) | -4.1 / 4.7% -> -8.0 / 8.0% | +5.5 / 6.2% -> +1.5 / 4.6% | +2.5 / 3.2% -> -1.5 / 2.5% |
| int8 streamed GEMMs (20) | -6.9 / 7.0% -> -11.1 / 11.1% | +1.3 / 3.2% -> -2.9 / 4.8% | -1.6 / 3.1% -> -5.8 / 6.0% |
| chained odd GEMMs (35) | -5.5 / 6.4% -> -9.1 / 9.1% | +2.0 / 4.1% -> -1.7 / 3.8% | -0.5 / 3.1% -> -4.2 / 4.7% |
| compute-bound isolated, kernel (15) | +7.2 / 7.2% -> -1.8 / 3.8% | +9.0 / 9.0% -> -0.0 / 4.4% | +8.2 / 8.2% -> -0.8 / 4.4% |
| -- M <= 512 (5) | +5.3 / 5.3% -> +0.0 / 4.3% | +7.0 / 7.0% -> +1.6 / 5.7% | +6.4 / 6.4% -> +1.0 / 5.2% |
| -- M > 512, row-tiled estimate (10) | +8.1 / 8.1% -> -2.8 / 3.6% | +10.1 / 10.1% -> -0.8 / 3.8% | +9.1 / 9.1% -> -1.8 / 4.1% |
| decode layer, full (12) | -8.6 / 8.6% -> -10.7 / 10.7% | -0.7 / 3.5% -> -2.8 / 3.9% | -3.2 / 4.0% -> -5.3 / 5.3% |
| decode layer, gemm_only (12) | -8.0 / 8.0% -> -11.2 / 11.2% | +3.6 / 3.6% -> +0.4 / 2.7% | -0.0 / 2.4% -> -3.3 / 3.8% |
| **pooled (102)** | **-4.0 / 6.5% -> -8.4 / 8.7%** | **+3.3 / 5.0% -> -1.2 / 4.3%** | **+0.8 / 4.0% -> -3.6 / 4.6%** |

- The change does what it was for: compute-bound GEMMs +8.2% -> -0.8% on refsb (the 5 points without row
  tiling +6.4% -> +1.0%; the 10 row-tiled ones carry the unvalidated row estimate, see above).
- It also takes ~4% off every other set, because their GEMMs have 2-56 strips too. So the pooled figure on
  refsb does NOT improve: 4.0% -> 4.6% (4.6% on the 92 points without row tiling), and NPUsim is now short on
  nearly everything (pooled mean -3.6%). The per-strip fill had been standing in for costs NPUsim does not
  model; with it gone they show: int8 weights (-5.8%; the TPU streams int8 3-8% slower per byte, cause not
  established), chained small GEMMs (-4.2%), the KV read at 819 GB/s in the decode layer (-5.3% full vs -3.3%
  gemm_only).
- Separate effect, not part of the table: with Codex's KV stream model (`kv_streams.py`, refsb probe) the
  decode layer on refsb is -2.8 / 2.9% and the pooled 102 points -3.3 / 4.3% (before the core change, KV-only:
  3.76%).
- `hbm2e_v5e.yaml` (82% sustained, too slow a DRAM) now pools best at 4.3% (mean -1.2%): a slower DRAM
  cancelling the same missing costs, not a better model.

### int8 weights: what the error is (2026-11-12; no simulator change)

With Codex's audit (`INT8_AUDIT.md`, `audit_int8.py`, `results/2026-11-12_int8_audit/`). `int8_decompose.py`
reproduces the numbers here from stored files (`results/2026-11-12_int8_core/decompose.log`); every constant in it
is fitted to TPU data and none is in the simulator.

**NPUsim's int8 path.** `tpu_v5e_int8w.cfg` differs from `tpu_v5e.cfg` in one line (`weight_format = int8`): same
mapping, tiles and schedule, half the weight bytes. No accounting error was found in it. What it cannot express is
that XLA compiles a different schedule per dtype (128 x 4096 x 4096: 8 windows, `kernel_window` 512x4, with bf16
weights and 11, 512x3, with int8 -- raw compiler metadata; the cast is folded into a mixed-type convolution, there is no separate conversion op or bf16 copy in HBM).

**Scope of the measured times.** `gemm.csv` / `gemm_int8.csv` of follow-up 2 are *loop slopes*: the stream bench's
`(median wall(2n) - median wall(n)) / n`, written into `device_us` too. One iteration is `(x + carry) @
W.astype(bf16)`, a float32 sum of the result and the carry update -- not a profiled kernel. `compare_v5e.py` now
records `timing_kind` (`--timing device_trace | device_slope | loop_slope | wall`; unstated, it is guessed and
marked `(inferred)`, or `unknown` when a file has only `device_us`). The stored files: follow-up 2 `gemm.csv`,
`gemm_int8.csv` = loop_slope; smallk2 `gemm_oddchain.csv` = device_slope; `rot/gemm.csv` = device_trace.

**HBM-bound (M <= 32), loop-slope scope:** `t = F + weight_bytes / r` over the four projections.

| | bf16 F (us) | bf16 r (GB/s) | int8 F (us) | int8 r (GB/s) |
|---|---|---|---|---|
| TPU, follow-up 2 and 3, M 1 / 8 / 32 | 1.3-2.9 | 751-757 | 1.4-3.3 | 733-765 |
| NPUsim, refsb | 0.25-0.60 | 737-742 | 0.14-0.50 | 731-742 |

The per-byte rate is the same for both dtypes and NPUsim has it within 2%. The difference is the fixed part per loop
iteration, ~2 us on the TPU and ~0.3 us in NPUsim, the same for bf16 and int8; with half the bytes it weighs twice
as much (4096 x 4096: 2 us of 46 us vs of 25 us). In bf16 NPUsim's slightly slower rate covers it. What the 2 us
is (loop body, fusion start-up) is not attributed; in the profiled fusion scope the int8 error of the same
predictions is smaller (Codex: MAPE 6.0% loop vs 4.9% fusion) but not gone.

**Compute ~ stream (int8, M = 128), device-trace scope (2026-11-04 profile run, fusion self time):** the fusion
takes 1.20-1.36 x max(stream, compute) on 11-75 windows, where a perfectly overlapped window pipeline gives
1.01-1.09; NPUsim, which overlaps perfectly, is 10-15% short on these four shapes (follow-up 2, loop slope). bf16 at compute ~ stream exists only with 2-3
windows (256 x 2000^2, 256 x 2048^2), where first DMA + last compute explain it. Whether this is int8 or "any
dtype at compute ~ stream with many windows" cannot be told from stored data.

**A hypothesis that fits and is not adopted:** a compute-side cost per XLA window, hidden under DMA slack when
HBM-bound. Estimated on the four int8 M = 128 cases (0.36-0.41 us) and applied to the other 42 multi-window cases
of that run: mean -3.3% / MAPE 4.1% without it, -0.0% / 2.7% with it (bf16 compute-bound -9.6% -> -3.1%). Against
it: the value is not one number (0.24-0.32 us on compute-bound cases with >= 56 windows, 0.17-0.19 on isolated
kernels), nothing shows what happens at a window boundary, and M = 512 int8 ~ bf16 only excludes a large *serial*
per-weight conversion, not a pipelined one. A per-window knob set from this would be a fitted constant.

**Discriminating measurement (one chip; not run; indications, no single ratio is a proof):** `measure_tpu_v5e1.py --only balance` -- the stream bench's loop for
bf16 and int8, K = N = 4096 and 8192, M 32-512, saving the loop slope, the profiled time of the loop's fused ops per iteration (all fusion ops, not
the GEMM alone), the complete op list, the traces and the optimized HLO of each program, with the one-op program
as a separately labelled control. bf16 at M ~ 256 also at ~1.2 x max
-> a pipeline effect to model for every dtype; only int8 -> a dtype-specific stage; excess per window constant
over K -> per window, per weight element constant -> per element.

### Balance run (2026-11-12, Colab v5e-1): the stream / compute excess is not an int8 effect

`results/2026-11-12_colab_v5e1_balance/` (`balance.csv`, `loopk*.csv`, the optimized HLO of all 56 programs;
`balance_analyze.py` -> `balance_analyze.log`; raw archive with traces: `SOURCE.sha256`). 28 loop cases: bf16 and
int8 weights, K = N = 4096 and 8192, M 32-512, on the stream bench's loop.

- **Scope check.** Each loop program has two profiled ops, the `while` and one fused op per iteration; the wall
  loop slope and that fused op's device time agree to 1.3% on average (5.2% at most); the loop itself costs
  0.4 us per iteration. The ~2 us fixed part seen before is inside the fused op.
- **bf16 shows the same excess at compute ~ stream.** fused / max(stream, compute) at C/S ~ 1: bf16 1.27 (K 4096,
  M 256), 1.18 (K 8192); int8 1.36, 1.26. It rises from ~1.0-1.1 at C/S 0.1-0.5 and falls back to 1.06-1.14 at
  C/S 2-4 for both dtypes. So the int8 M = 128 error is the balanced regime, which int8 reaches at half the M.
- **XLA's windows follow M and dtype** (HLO `iteration_bounds`): K 4096 bf16 has 8 windows up to M 192, 11 at
  M 256, 16 at M 512; int8 11 up to M 192, 16 above; K 8192: 32-52 (bf16), 39-64 (int8).
- **Predictions fixed before the run** (R 752 GB/s, F 2 us, w 0.393 us per window, all from earlier data), error
  against the fused op, mean / MAPE:

| | perfect overlap (max + F) | window pipeline | window pipeline + w |
|---|---|---|---|
| all 28 | -8.0 / 8.0% | -5.4 / 5.8% | -0.2 / 2.6% |
| bf16 (14) | -7.5 / 7.6% | -4.5 / 5.1% | -1.8 / 2.7% |
| int8 (14) | -8.4 / 8.4% | -6.2 / 6.6% | +1.5 / 2.6% |
| C/S < 0.6 (10) | -2.1 / 2.2% | -0.0 / 1.3% | +0.7 / 1.7% |
| C/S 0.6-1.6 (10) | -13.8 / 13.8% | -9.7 / 9.7% | -2.4 / 4.2% |
| C/S > 1.6 (8) | -8.1 / 8.1% | -6.6 / 6.6% | +1.5 / 2.0% |

  NPUsim today behaves like the first column. The per-window cost predicted a new run, new shapes and bf16 in the
  balanced regime without refitting; worst cases bf16 384 x 4096 (-8.3%) and 192 x 8192 (-7.3%).
- **Control, a different schedule** (one-op program, operand copied first when <= ~67 MB): the copy runs at 750
  GB/s and the kernel after it is the compute roofline + 2.5-7 us (K 4096) or 7-15 us (int8 K 8192), with no
  bump at M 256. Streaming hides 0.71-0.94 of the shorter of the two phases (`eta`). So the excess belongs to
  streaming and computing at the same time; whether it sits at window boundaries (what `w` assumes) or is
  interference between DMA and compute is still not shown, and the staged kernels suggest less than w per window
  (0.2-0.3 us).
- **loopk.** A `fori_loop` iteration adds +0.19 / +0.39 / +0.36 / -0.13 / +0.74 us to an unrolled kernel (five
  shapes): small. The -22% on odd GEMMs inside a `while` loop is not the loop. Unexplained: the unrolled
  per-kernel cost of 100 x 1000 x 1000 here is 4.1 us (single window: stream + compute in series), where the
  smallk2 `oddchain` slope gave 2.81 us for the same GEMM (that program staged 11 us of copies); the
  "chained odd GEMM" set rests on the latter schedule.

Not changed: simulator and configs. What this run supports modelling next, for every dtype: weights arriving in
XLA-sized windows with a cost per window on the compute side -- its value measured (0.2-0.4 us), its mechanism open.

### Window-schedule candidate in NPUsim (2026-11-12; opt-in, not adopted)

`[accelerator] offchip_window_count` and `window_switch_cycle` (both default 0; no config uses them) and the driver
`balance_npusim.py`. Implementation, tests, the 28-case A / B / C results on both DRAM models, the three
regressions, the residency fallback of 12 cases and what is not established:
[`results/2026-11-12_window_candidate/REVIEW.md`](results/2026-11-12_window_candidate/REVIEW.md).
Primary comparison (DRAM traffic as the compiler accounts the measured loop, `fused_output_consumer`; refsb, same
cases and runs): A (config as is) MAPE 8.07%, B (window schedule) 5.84%, C (+ 590 cycles per window, fixed before
the measurement) 2.73%; bf16 6.89 -> 3.13%, int8 9.25 -> 2.32%; 24 of 28 closer, four memory-bound bf16 cases
worse (+0.5 .. +3.3 points) through the on-chip weight-load tail, whose existence on the TPU is the open
measurement. An earlier batch on another activation assumption (2.93%) is kept as history; the first batches of
that driver ran on wrong inputs and are withdrawn (`INVALID_SHARED_NETWORK/`).

### Tail probe (optional; prepared 2026-11-12, NOT run -- needs one TPU v5e chip)

`measure_tpu_v5e1.py --only tailprobe` (TPU only at full size; `--quick` on CPU is a functional check). For the
four memory-bound bf16 cases the window-schedule candidate over-predicts (32 / 64 x 4096 / 8192) and one int8
control (32 x 4096), the same loop body as `balance` runs as two programs: **banked** (a bank of distinct
weights, a new one streamed every iteration) and **single** (one weight reused by every iteration).

Three kinds of numbers come out, and they are kept apart:

| kind | columns | what it is |
|---|---|---|
| wall clock | `wall_n_s`, `wall_2n_s`, `loop_slope_us`, `slope_flag` | raw host samples of the n- and 2n-iteration loop and their SIGNED slope; a nonpositive slope is flagged, not clipped |
| profiler events | `modules_traced`, `op_event_counts`, `op_event_us`, `conv_fusion_events_per_module`, `conv_fusion_us_per_event`, `all_fusion_event_sum_per_requested_iter_us` | device events, one per op execution; no per-DMA events. The per-event convolution time is given only when the trace has exactly 8 convolution-fusion events per 8-iteration module; the all-fusion sum may include fusions outside the loop and is not a loop-kernel time |
| compiler facts | `hlo_facts`, `hlo/*.txt` | read from the optimized HLO: convolution operand types, window geometry, memory-space tags, copies / prefetch, fusion and loop counts -- static, not measured |

`tailprobe_pairs.csv` gives one conservative verdict per case:

- *inconclusive* when the two HLO texts do not show the same signatures (convolution types, window geometry,
  fusion / loop counts), or when one convolution-fusion event per loop iteration is not verified in BOTH traces
  (equal global counts in the HLO do not show that the GEMM stayed inside the loop);
- otherwise *same observed signatures ...; weight staging NOT ESTABLISHED*. Staging is never concluded: a
  memory-space tag on a weight-shaped tensor or extra copy / prefetch ops in the single-weight program are listed
  as raw evidence (a copy may be of the input; a tag is not a traced dataflow), and their absence does not show
  the weight was not staged.

What it cannot give: the time between the end of the last DMA and the end of the fused op -- there is no event
for it -- nor where the reused weight physically resides. A banked-vs-single difference is a difference between
two programs; it is not to be reported as a measured tail, and no model parameter is to be fitted to it. One
run may well end inconclusive. `test_tailprobe.py` (CPU, functional only): `hlo_facts` on HLO saved from the
balance run, `op_events` on synthetic traces (8 per module verified; 9, 1 and none are not), `pair_verdict`
cases (an extra copy alone or a tag alone never yields a staging conclusion), the `--quick` run and the TPU
guard. `tailprobe_meta.json` records script hash, argv, backend, device kind, JAX version; each row the HLO and
trace file hashes and requested vs traced module counts.

Run (Colab, TPU v5e-1): `%run measure_tpu_v5e1.py --only tailprobe`, then zip the output directory (it keeps
`traces/`).
