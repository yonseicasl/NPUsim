# TPU v5e int8 path audit — 2026-11-12

The current evidence does **not** establish a missing standalone dequantization
stage. It establishes differences in compiler scheduling and a mismatch between
the measured loop and a standalone GEMM comparison. A fitted int8 delay would
conceal these distinctions. This audit changes no simulator latency constants.

## Reproduction and evidence

Run `python3 audit_int8.py --source <original validation/tpu_v5e> --output
<new results directory>`. The source must contain followup2, followup3 and the
2026-11-04 profile run. On washington the source is the npusim-engine NPUsim
checkout, and the output is the tpu-measure worktree's
`results/2026-11-12_int8_audit/`.

The script verifies unique shape/dtype/cold keys and 80 rows per run (20 shapes,
2 dtypes, cold/warm), joins 120 bf16/int8 pairs, extracts 293 optimized HLO window
records and 118 GEMM fusion profile records, and records SHA-256 for 137 input
files (139 when both cached prediction CSVs are supplied). `paired.csv` contains the unrounded measured times. `audit.json` preserves
HLO window metadata, profile durations, source paths and hashes. Multiple HLO
window records can belong to one module; counts are not distinct kernels.

## What the experiment actually measures

`tpu_v5e1_stream_bench.py:per_iter` measures
`(median(wall(2*n)) - median(wall(n))) / n`, after compilation and warmup, with
`block_until_ready`. Its exporter writes this value into **both** `wall_us` and
`device_us`. The latter column is therefore not a directly profiled device
kernel time. A comparison row labelled `src=device` is ambiguous for this data.

One dependent loop iteration computes `(x + carry) @ weights.astype(bfloat16)`,
converts and sums the output, and updates the carry. Program startup largely
cancels in the slope, while per-iteration input arithmetic and output reduction
remain. Isolated GEMM latency, loop-slope time and fusion profile time must be
reported as different scopes. The existing `--in-program` option removes program
startup but does not model every fused loop operation.

This is int8 **weight storage** with bf16 activations/result semantics, not a
benchmark of native int8 activation-by-weight arithmetic. The source does a cast,
with no scale multiplication or zero-point subtraction. Optimized HLO for
`gemm_8x4096x6144_int8_cold` passes a bf16 activation and s8 weight directly to
`convolution.5` inside the output-reduction fusion. The source weight conversion
is folded into that operation. This does not prove conversion is free; it means
that a separate full bf16 weight tensor in HBM is not supported by this artifact.

The current followup2 predictions already use `--in-program`:
`npusim_us == npusim_noprog_us`. No further 0.69 µs subtraction is warranted.

With `--prediction-bf16 <CSV> --prediction-int8 <CSV>`, the audit compares
unchanged cached predictions with each measured scope. It rejects predictions
that still include program startup and writes `scope_comparison.csv`.
On the Nov4 run (20 shapes per dtype and scope):

| Dtype | Loop-slope MAPE | Profile-fusion MAPE |
|---|---:|---:|
| bf16 | 3.132% | 2.254% |
| int8 | 6.033% | 4.877% |

This is a change of comparison scope, **not a simulator improvement**. A fusion
also includes compiler-fused arithmetic around GEMM, so it is not pure matrix
multiplication. Int8 underprediction remains: mean signed error −4.554% against
the fusion. Reusing these predictions against the original followup2 loop
measurements reproduces bf16 MAPE 2.528% and int8 MAPE 5.970%.

## Dtype changes the compiled schedule

For the Nov4 `8×4096×6144` cold case:

| HLO field | bf16 weights | int8 weights |
|---|---|---|
| `iteration_bounds` | 10,1,1 | 2,1,7 |
| `kernel_window_bounds` | 512,5 | 80,24 |
| `input_window_bounds` | 1,32 | 1,5 |
| `output_window_bounds` | 1,5 | 1,24 |
| Profile fusion average | 68.78578125 µs | 35.30578125 µs |

These are raw compiler metadata; their units are not assumed to be literal
matrix elements. The difference shows that simply halving weight bytes while
keeping an identical schedule misses compiler choices. It does not alone prove
which schedule change accounts for the timing error. `estimated_cycles` is a
compiler estimate, not a hardware counter.

## Why a constant int8 penalty is not justified

Define the diagnostic `int8_us - bf16_us/2`. It is **not conversion time**:
halving bf16 also halves costs shared by both dtypes, and compute-bound GEMMs
should not halve their execution time at all.

For M=8, shapes ordered as (4096,6144), (4096,4096), (4096,28672),
(14336,4096), cold loop residuals are:

| Run | Residuals, µs |
|---|---|
| followup2 | 1.297, 1.968, −0.188, 0.346 |
| followup3 | 1.481, 0.925, 1.115, 1.171 |
| Nov4 profile run | 0.751, 0.825, 1.148, 0.793 |

The Nov4 M=1, K=N=4096 loop residual is −0.505 µs; the corresponding profiled
fusion residual is +1.007 µs. These scopes cannot be substituted silently.
For M=512 the diagnostic rises sharply because computation dominates; this is
not evidence for hundreds of microseconds of conversion overhead.

`warm` in these data means one weight copy, not demonstrated on-chip residency.
Individual int8 weight matrices contain 16–112 MiB, so the warm/cold contrast
cannot be used as a clean conversion-only experiment. Logical activation/output
bytes are recorded for sensitivity analysis; they are not asserted to equal
physical HBM traffic under compiler tiling and fusion.

## Follow-up and acceptance criteria

1. Make comparison timing provenance explicit. Preserve published baselines and
   report a correctly scoped in-program comparison separately.
2. Compare dtype-dependent XLA partitioning against NPUsim's mapping, transfer
   chunking and fill/overlap model. Any core change needs an independently
   derivable accounting or scheduling argument, not only lower fitted error.
3. Keep bf16 controls and memory/compute regimes in regression tests. Publish
   improvements and regressions separately rather than hiding them in a mean.
4. If existing artifacts cannot isolate the cause, a future one-chip experiment
   should save loop slope, GEMM fusion duration, full optimized HLO and identical
   dependency controls for both dtypes. It does not require eight chips. No new
   TPU measurements were performed for this audit; the user's v5e-8 allocation
   remains unavailable.

Claude owns core/frontend changes and their validation; Codex owns this audit.
The experiment artifacts alone do not establish a completed core accuracy gain.

## Completed collaboration follow-up

Local commits on washington's `tpu-measure` branch:

- `5f59cef`: independent audit and preserved evidence.
- `45122a3`: timing provenance, reproducible decomposition, balance experiment,
  and correction of the earlier claim that int8 has a lower per-byte stream rate.
- `29df8f0`: full balance run requires TPU; causal readings explicitly qualified.
- `314db94`: complete op/trace retention; aggregate fusion time explicitly named
  `fused_ops_sum_per_iter_us`, not mislabelled as GEMM-only time.

The comparator now preserves ambiguous input as `unknown`, marks inferred
provenance, and supports explicit `--timing`. A mixed missing-device-time crash
found during review was fixed and added to tests. Independent comparison,
trace-parser and KV tests passed (KV: five cases). A temporary CLI integration
check verified explicit overrides, ambiguous defaults and unchanged numeric
predictions. Claude also completed a CPU quick smoke run of the new experiment;
this is not TPU validation.

The prepared command on an available single-chip TPU runtime is:

```sh
python3 measure_tpu_v5e1.py --only balance --keep-traces --out balance_run
```

This uses the original dependent streaming loop for both dtypes, across
M=32,64,128,192,256,384,512 and K=N=4096,8192. It saves loop slope, module time,
aggregate fused-op time, all op records, traces and HLO. Isolated controls are
labelled separately. Full execution refuses a CPU backend.

No core latency/configuration change was justified by current evidence. The
original int8 loop MAPE remains 5.970%; this phase improved validation provenance
and prepared the experiment needed to choose a structural model correction.
Nothing was pushed in this collaboration phase. Hardware confirmation remains
pending TPU availability; no background monitoring of allocation is configured.
