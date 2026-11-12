# Window-schedule candidate: implementation and result note (2026-11-12, for joint review)

Status: **candidate, opt-in, not adopted.** Both knobs default to 0; no existing config, baseline or published
number changes. 590 cycles is an empirically fitted effective cost, not a hardware constant. Nothing here is a
gain over an earlier cohort: A, B and C are compared on the same 28 cases, same runs, same DRAM model.

## Final status (2026-11-12)

- **Primary comparison: the matched-traffic batch, `full_fused_output/`** (section "Full batch with
  `fused_output_consumer`"): the DRAM traffic of every run is what the compiler accounts for the measured
  program. refsb, A / B / C MAPE 8.07 / 5.84 / 2.73% (C bias -0.39%); bf16 6.89 -> 3.13%, int8 9.25 -> 2.32%;
  24 of 28 closer (14 of 14 int8); **four regressions, all memory-bound bf16**. Exact figures:
  `codex_independent_review/npusim_codex_full_fused_review.json` (Codex's independent check of both batches).
- The residency2 batch (section "Results") is a valid historical run on a different activation assumption; it is
  no longer the primary comparison.
- 590 cycles is unchanged since before the measurement; nothing was refitted. The traffic is the compiler's
  static MATRIX-payload accounting (weights and the input matrix read; no output matrix), not measured traffic:
  the fused op's 4-byte scalar write, its scalar reads and its extra work (input add, output conversion and
  reduction) are not modelled. All knobs opt-in, default 0; no config, no push, no merge.
- Open: the memory-bound bf16 over-prediction (section "Remaining evidence gap").

## What was implemented

Core (`scheduler/stats.cc/.h`, `npu.cc/.h`), `[accelerator]`:

- `offchip_window_count` (0 = off). The off-chip stages deliver the layer's data in that many windows instead of
  once per tile of the timeline clock: a window is a run of consecutive tiles, its whole off-chip cost sits on its
  first tile, DRAM and the multi-chip link carry it at once (cut-through: max of the two, not the sum), and the
  multi-chip -> GLB boundary's one tile of look-ahead becomes one window (windows double-buffered as wholes).
  The on-chip stages keep the fine clock. Refused with an error: more windows than tiles; a single-buffered or
  GLB-bypassed multi-chip -> GLB boundary. The DRAM <-> multi-chip depth is left as configured.
- `window_switch_cycle` (0): one compute-side cost per window, on the window's last tile, joined to the compute
  axis of that tile inside the PE stage's own per-axis combination (hidden when another PE axis is longer);
  `busy_cycle_pe` then reports the scheduled sum.
- **It is a count-based, equal-share schedule**: equal off-chip share per window, near-equal runs of tiles. The
  compiler's window geometry and partial edge windows are not represented (14 of the 28 cases have them).

Tests: `unittest/stats_timeline_test.cc` `run_window_cases`, 25 hand-derived cases (off, DMA-hidden,
compute-exposed, crossover, 1 / 3 windows, uneven 3 + 2 with and without w, count vs total-cost counterexample,
cut-through incl. link-bound, explicit outstanding limits, LB-bound / partly exposed / single-buffered LB, PE busy
56 / 50 / 50). `run_validation.sh` passes; with the knobs off, 6 GEMM runs are cycle-identical to commit
`0252029`'s binary. Driver: `balance_npusim.py`; `test_balance_inputs.py` (inputs of every run, incl. the
shared-network bug reproduced and rejected).

## The experiment

28 loop cases of `results/2026-11-12_colab_v5e1_balance` (bf16 / int8 weights, K = N 4096 / 8192, M 32-512).
Measured = the loop's one fused op per iteration (device trace). NPUsim: each GEMM simulated whole (no strip
extrapolation), `program_setup_cycle` 0, `activation_residency` 2 requested, window count = `iteration_bounds` of
the case's own HLO (an input, not predicted). A = config as is, B = + `offchip_window_count`, C = B +
`window_switch_cycle` 590 (0.393 us, fixed before the measurement, one value for both dtypes).
Provenance: `FROZEN.sha256`, `balance_npusim_<dram>.json` (binary, library, configs, git, HLO and csv hashes),
`generated_<dram>/<run>/` (cfg, network, mapping, both NPUsim reports of all 168 runs). Codex checked both
batches independently (unique inputs, MAC counts, report = CSV, hashes). `INVALID_SHARED_NETWORK/` holds the
earlier batches, which ran every case on the last shape's network: **withdrawn, not to be cited.**

## Results of the first batch, `--activations residency2` (historical; `summary.txt`)

Mean / MAPE in %, error vs the measured fused op.

Regime = compute / stream with compute = 2 M K N / 196.6 TFLOPS and stream = weight bytes / 750 GB/s.

| refsb DRAM | n | A | B | C | cases where C is closer than A |
|---|---|---|---|---|---|
| memory-bound, C/S < 0.6 | 10 | -4.09 / 4.23 | -1.11 / 2.66 | -0.39 / 2.37 | 7 of 10 |
| balanced, 0.6-1.6 | 10 | -13.96 / 13.96 | -9.99 / 9.99 | -3.32 / 4.45 | 10 of 10 |
| compute-bound, > 1.6 | 8 | -8.45 / 8.45 | -6.84 / 6.84 | +1.29 / 1.73 | 8 of 8 |
| bf16 | 14 | -7.49 / 7.60 | -3.90 / 5.01 | -1.62 / 3.07 | 11 of 14 |
| int8 | 14 | -10.23 / 10.23 | -7.93 / 7.93 | -0.30 / 2.79 | 14 of 14 |
| full windows | 14 | -6.20 / 6.30 | -3.28 / 4.39 | +0.52 / 2.14 | 11 of 14 |
| partial edge windows | 14 | -11.53 / 11.53 | -8.56 / 8.56 | -2.43 / 3.71 | 14 of 14 |
| activations actually on chip | 16 | -9.07 / 9.16 | -5.37 / 6.34 | -1.53 / 2.87 | 13 of 16 |
| activations to / from DRAM (fallback) | 12 | -8.58 / 8.58 | -6.64 / 6.64 | -0.20 / 3.00 | 12 of 12 |
| (all 28, mixed residency) | 28 | -8.86 / 8.91 | -5.92 / 6.47 | -0.96 / 2.93 | 25 of 28 |

All aggregates in this note are computed from unrounded per-run errors (`window_candidate_summary.py` ->
`summary_exact.txt`; the earlier `summary.txt` used the CSV's one-decimal percentages).

Analytic DRAM (819 GB/s, ~9% faster than the measured copy rate): A -12.57 / 12.57, B -9.72 / 9.72,
C -4.32 / 5.53; its memory-bound group stays at -7.91% (the DRAM rate), compute-bound +1.16 / 1.73, balanced
-5.11 / 6.20.

Worst C (refsb): bf16 384x4096 -8.42, bf16 256x4096 -6.96, bf16 256x8192 -6.29, int8 128x4096 -5.73,
int8 192x8192 +5.37, int8 128x8192 -5.26: the balanced regime is improved, not solved.

## The three regressions (refsb; all memory-bound bf16)

| case | A | B | C |
|---|---|---|---|
| bf16 32 x 4096 | -2.26% | +3.24% | +4.09% |
| bf16 32 x 8192 | +0.72% | +2.14% | +2.36% |
| bf16 64 x 8192 | -0.07% | +1.34% | +1.56% |

From the saved reports, in cycles (1500 = 1 us), identical in all three: **B - A = 3840, C - B = 590**.

- C - B is exactly one `window_switch_cycle`: when memory-bound only the last window's cost is exposed. It is
  not zero: it adds 0.85, 0.22 and 0.22 points to the signed error of these three cases.
- B - A is one window of the **"PE array" stage** -- the on-chip weight load into the MXUs, which this config
  models at 128 weights / cycle / MXU: 4096 cycles per bf16 window here, independent of M -- minus the two fine
  tiles (256 cycles) that A already exposed. Per window the reports give DRAM 8475, link 4096, GLB 2048,
  PE array 4096, PE 1034-2058; tail after the first stage: A 352-384, B 4192-4224, C 4782-4814. int8 32 x 4096
  shows the same structure (B - A 1344 = 1489 - 145).
- **Demonstrated** (simulator): fetching windows whole moves the last window's on-chip weight load, a stage that
  existed before, from under the DMA to after it. **Not demonstrated** (TPU): that the chip spends that time
  after the last DMA. The only measured hint is cross-program -- loop fused op minus the one-op program's copy
  (46.49 - 44.76 = 1.73 us at 32 x 4096) against NPUsim's 3.2 us tail -- and a cross-program difference is not an
  attributed tail. No change to that stage is proposed from this.

## Activation residency: requested, not always granted

`activation_residency = 2` asks for the run's input and output to stay on chip. NPUsim grants it only when the
whole tensor fits the GLB's residency capacity (2 MiB here). The reports show: 16 cases "input GLB-resident,
output retained in GLB"; **12 cases "input DRAM-backed, output materialized to DRAM"** (M K >= 1.5 M elements:
bf16 / int8 384, 512 x 4096 and 192-512 x 8192). So those 12 carry 2 x M x K x 2 bytes of DRAM traffic per GEMM
that the other 16 do not, and neither matches the benchmark, whose fused op reduces the output on chip and reads
its input every iteration. The 12 are mostly compute-bound, where DRAM is hidden, but two balanced ones are
affected: bf16 256 x 8192 is DRAM-bound in NPUsim on refsb (C - B = 0.39 us only, error -6.3%) while the same
case on the analytic DRAM exposes the window cost (-6.7% from -14.0%). Their good average (C 3.0%) does not
validate either residency assumption. Buffers were not enlarged to make them fit.
What the fused op actually touches is in the next section.

## Byte accounting of the fused op (added 2026-11-12, from the kept traces; no simulation)

`trace_bytes.py` reads each program's XProf `xplane.pb`; output
`results/2026-11-12_colab_v5e1_balance/trace_bytes.csv` (302 op records of 86 programs).

**Kind of evidence: static compiler accounting, not measured traffic.** `bytes_accessed` and
`memory_access_breakdown` are stats of the op's event *metadata* (one record per HLO op, shared by its 1600
executions), written by the compiler's cost analysis. The traces contain no per-DMA events (the other device
lines are empty) and no hardware byte counters. Memory spaces are raw ids: 1 = operands with the default layout
(no `S(n)` tag in the HLO), 3 = the destination of a copy tagged `S(1)`; XProf names them HBM and VMEM, which the
artifacts cannot confirm.

- **Loop (all 28 cases):** the one fused op per iteration reads `weight bytes + M x K x 2` from space 1 and
  writes 4 bytes; it touches no space 3. Operands in the HLO: the weight stack and the `M x K` bf16 input, both
  default layout; two `S(6)` scalars; result `f32[]`. So per iteration the compiler accounts **weights + the
  input activation read, and no M x N output anywhere** -- the input add and the output reduction are inside
  the op.
- **One-op program (control):** the staged operand is the weights when they are <= 67 MB (copy reads space 1,
  writes space 3; the fusion then reads the weights from space 3, the input from space 1 and writes the
  `M x N` output to space 1); with 134 MB of bf16 weights the *input* is staged instead and the weights are read
  from space 1. A different schedule and different traffic from the loop, as assumed.

Against that, NPUsim's two outcomes for these runs:

| | weights | input read (M K 2 B) | output write (M N 2 B) |
|---|---|---|---|
| compiler accounting of the loop | yes | yes | no |
| 16 cases, residency granted | yes | **no** | no |
| 12 cases, capacity fallback | yes | yes | **yes** |

The missing input read is 0.8% (M 32) to 6.25% (M 256) of the bf16 weight bytes at K 4096 and twice that for
int8; the extra output write is the same amount for the 12 fallback cases (up to 12.5% of the int8 weights).

**Smallest correction this supports (proposed, not implemented):** a residency setting that drops only the
run's *output* materialization -- consumed on chip as it is produced, so without the whole-tensor capacity
condition, which is what a fused reduction does -- while the input is read from DRAM. One mode, default
unchanged, unit-tested on the DRAM traffic it produces; then re-run the 28 cases (A / B / C unchanged
otherwise). Expected direction, from the byte fractions only: slightly more DRAM time for the 16 (the three
regressions move further up by <= 1%), less for the 12; bf16 256 x 8192 on refsb becomes compute-bound and
exposes the window cost. It is not expected to remove the balanced-regime residual.

**Still missing, and not obtainable from these artifacts:** when inside the fused op the DMA ends (no per-DMA
or per-window events), hence no same-program attribution of the on-chip weight-load tail; and measured (not
compiler-accounted) HBM bytes.

## What can and cannot be said

- In NPUsim itself, on a new measurement, the window schedule plus one fixed per-window cost removes most of the
  balanced-regime and compute-bound deficit for both dtypes with one value: refsb 8.9% -> 2.9%, 25 of 28 closer.
- It makes three memory-bound bf16 cases worse: their absolute error rises by 1.83, 1.64 and 1.49 percentage
  points (signed shifts +6.35, +1.64, +1.63). Primarily the whole-window dependency (B - A = 3840 cycles), with
  an additional exposed final-window cost (C - B = 590 cycles: +0.85, +0.22, +0.22 points of the signed shift).
- Overfitting is NOT ruled out. A value fixed before this measurement transferring to it supports that one
  scalar; but it is the same hardware and compiler, the window counts are measured inputs, and the schedule was
  designed after seeing these data. That is not validation in an independent environment.
- Not established: the physical mechanism of the per-window cost; exact window geometry; the residency of the
  loop's activations; transfer to isolated / staged programs or to shapes without an HLO. 28 cases from one
  session and one compiler version are not 28 independent environments.

## Proposed follow-up (nothing changed yet)

1. Keep the knobs opt-in; no config promotion on these numbers.
2. Done: byte accounting from the kept traces (section above). Next: the output-only residency mode, tested,
   and the 28 cases re-run on it.
3. For the on-chip weight-load tail: the kept traces cannot attribute it (one event per fused op). Without
   per-DMA events nothing here measures "DMA end to op end" directly; a loop on one reused weight is only an
   additional contrast (and the compiler decides whether that weight is staged at all). See the last section.
4. Geometry: represent full and edge windows (unequal shares) only if 2-3 show the equal-share approximation
   matters; the 14 partial-window cases are 3.7% against 2.1%.
5. Only after 2-4: the older bf16 / int8 stream sets (follow-up 2, with HLO windows where an HLO exists), and a
   decision on promoting a named calibration profile.

## Full batch with `fused_output_consumer` (added 2026-11-12; `full_fused_output/`)

Same 28 cases, same 590 cycles, same window counts; only the activation traffic changed to what the compiler
accounts for the loop (weights + M x K x 2 read, no output) -- checked in every one of the 168 runs from the
`[dram] line_size` of its saved config. A traffic-consistency correction: it is not offered as a timing gain,
and the comparison below is A / B / C within this batch (`full_fused_output/summary.txt` has every split and the
per-case change against the residency2 batch). Mean / MAPE in %:

| refsb DRAM | n | A | B | C | C closer than A |
|---|---|---|---|---|---|
| memory-bound | 10 | -2.66 / 3.02 | +0.32 / 1.92 | +1.03 / 1.78 | 6 of 10 |
| balanced | 10 | -12.81 / 12.81 | -8.85 / 8.88 | -3.07 / 4.41 | 10 of 10 |
| compute-bound | 8 | -8.46 / 8.46 | -6.96 / 6.96 | +1.18 / 1.80 | 8 of 8 |
| bf16 | 14 | -6.63 / 6.89 | -3.10 / 4.72 | -0.97 / 3.13 | 10 of 14 |
| int8 | 14 | -9.25 / 9.25 | -6.96 / 6.96 | +0.20 / 2.32 | 14 of 14 |
| full windows | 14 | -5.50 / 5.76 | -2.63 / 4.25 | +1.16 / 2.20 | 10 of 14 |
| partial edge windows | 14 | -10.38 / 10.38 | -7.43 / 7.43 | -1.94 / 3.25 | 14 of 14 |
| all | 28 | -7.94 / 8.07 | -5.03 / 5.84 | -0.39 / 2.73 | 24 of 28 |

Analytic DRAM: A -11.96 / 11.96, B -9.13 / 9.13, C -3.75 / 4.94 (compute-bound +1.05 / 1.80%, balanced
-4.85 / 5.91%; memory-bound stays short by the 819 GB/s). Exact tables: `full_fused_output/summary_exact.txt`.

- **Four regressions now, all memory-bound bf16** (absolute error, points): 32 x 4096 +3.35 (-1.50 -> +4.85),
  32 x 8192 +1.64, 64 x 8192 +1.63, 64 x 4096 +0.49. With the input read restored the A baseline of these cases
  is closer (bf16 memory-bound A: -2.91 -> -1.61%), so the whole-window tail costs more against it: for
  memory-bound bf16 as a group the candidate is slightly WORSE on MAPE (A 2.22%, C 2.37%; signed -1.61 -> +2.29).
- The balanced residual is unchanged in size (-3.07 / 4.41%; worst bf16 384 x 4096 -8.92, 256 x 8192 -6.62,
  256 x 4096 -6.51, int8 192 x 8192 +5.31, int8 128 x 4096 -5.31); compute-bound is within 2%.
- int8: 14 of 14 closer, C +0.20 / 2.32%.

Reading: with consistent traffic the picture of the residency2 batch holds -- the window schedule with one
fixed per-window cost removes most of the balanced and compute-bound deficit (refsb 8.07% -> 2.73%), and it
over-predicts memory-bound bf16 through the on-chip weight-load tail. Whether that tail is physical is still
the open measurement; until then the candidate should not be enabled for memory-bound streaming.

### Effect of the traffic correction alone (refsb, variant C, same cases)

Reported separately because it is a semantic correction, not a model improvement: overall MAPE 2.93 -> 2.73%;
int8 2.79 -> 2.32% (better); bf16 3.07 -> 3.13% (slightly worse); 16 cases closer, 12 farther. The 16 cases whose
activations had been on chip: 2.87 -> 2.33%; the 12 that had fallen back to DRAM: 3.00 -> 3.25%. Matching the
program's traffic does not guarantee a better time in each case.

The four regressions, exactly (signed A -> C, increase of absolute error): bf16 32 x 4096 -1.50 -> +4.85%
(+3.35 points); 32 x 8192 +1.11 -> +2.76% (+1.64); 64 x 8192 +0.71 -> +2.34% (+1.63); 64 x 4096 -2.86 -> +3.35%
(+0.49). Nominal regime MAPE of C: memory-bound 1.78%, balanced 4.41%, compute-bound 1.80%; full windows 2.20%,
partial 3.25%; worst bf16 384 x 4096 -8.92%.

## Remaining evidence gap, and the smallest measurement that closes it

**Gap.** In NPUsim the window schedule exposes, after the last DMA, one window of the on-chip weight-load stage
plus one window cost: 3.2 us for a bf16 window (4783 cycles), 1.4 us for an int8 one (2127). Whether the TPU
spends that is not known: the traces have one event per fused op.

A conditional diagnostic only, NOT a bound and NOT a measured tail: writing the measured fused op of the
memory-bound rows as `bytes / r + c0` through the two points K = 4096 and 8192 (bytes = weights + input) gives
c0 = 1.84 us (bf16 M 32), 2.51 (bf16 M 64), 1.42 (int8 M 32), 1.32 (int8 M 64). Those two points are different
compiled programs (4x the weights, 2x the activations, another window count and possibly another effective
rate), the line has no residual degree of freedom, and c0 is whatever an assumed straight line leaves over. It
neither bounds the last-window time nor tests the weight-load stage. No conclusion is drawn from it.

**Next observation needed, and what it can and cannot give.** No available trace has per-DMA events, so the
"DMA end to fused-op end" time cannot be measured directly with these tools; every option below is a contrast
between programs, and one Colab run may be inconclusive because the compiler chooses the memory space, the
window geometry and what it hoists.

Prepared, not run: `measure_tpu_v5e1.py --only tailprobe` (see README "Tail probe"). For the four regressing
bf16 cases and one int8 control it runs the same loop body twice -- on a bank of distinct weights (streaming)
and on ONE reused weight -- and records, per program, the fused-op event duration, the wall loop slope, the
full HLO and trace, and the compiler's actual choices (weight layout and memory space, copies, convolution
operand types, fusions per iteration, window geometry). The pair is reported as comparable only if those
choices match apart from the weight's staging, and as "staged" only if the HLO shows it; otherwise the result
is flagged inconclusive. It reports durations and compiler facts; it does not subtract programs into a "tail".

Until such an observation exists the candidate is not to be enabled for memory-bound streaming, and no
parameter is refitted.
