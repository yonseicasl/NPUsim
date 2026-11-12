# TPU v5e comparison audit — 2026-11-11

This audit was performed with the existing Claude session `cfc77b0e-68da-4b85-8105-5c2791807a97`
(`tpu 정합성 평가 계획`) on washington. No new TPU measurements were taken. Kaggle v5e-8 allocation is
still pending; multi-chip collective and tensor-parallel fidelity remains unvalidated.

The working tree is
`/home/bogil/MOSAIC/.claude/worktrees/npusim-engine/simulator/3rd/npusim/.claude/worktrees/tpu-measure`,
branch `tpu-measure`, baseline `cb08bea`. Original measurement CSVs and JSONs are preserved.

## Consistent DRAM supply for KV reads

The old decode comparison used Ramulator for projection GEMMs but divided KV bytes by the ideal
819 GB/s for every DRAM model. `compare_stack.py --kv-times TAG=CSV` can now use a modeled KV read
for the corresponding DRAM tag. Default behavior remains unchanged. `kv_streams.py` produces this
CSV from the existing MOSAIC `dram_bw_probe` with the same YAML used for the GEMMs. It never consumes
measured TPU timings as modeling parameters.

Geometry: bf16, 8 KV heads, head dimension 128; total bytes = `B * L * 8 * 128 * 2 * 2`.
The probe verifies actual read bytes and zero write bytes. It uses clock ticks rather than the rounded
printed microsecond value. Up to 32 MiB is directly simulated; larger contiguous streams use the
16-to-32 MiB simulated slope with startup counted once. Every row labels `exact` or `extrapolated`.
A separate 64 MiB simulation checks that extrapolation: 89.254375 us predicted vs 89.251875 us simulated,
relative difference **0.002801%**. This supports the local DRAM extrapolation, not full attention fidelity.
Set `--max-probe-mib` at least as large as the largest KV stream to disable extrapolation for all cases.

This is a contiguous memory supply approximation. It omits QK/softmax/AV, vector work, layout costs,
competing streams and tensor-parallel communication. Small attention cases remain latency-sensitive.
The NPUsim core's attention implementation is not changed by this audit.
Claude's independent review identified an additional distinction: the projection GEMMs use NPUsim's
fresh approximately 1 MiB Ramulator streams (about 715 GB/s for this setting), while KV uses the probe's
full-pressure contiguous stream (about 751 GB/s). Thus the YAML is shared but the request frontends differ.
The KV-only gain partly reflects this frontend difference; it is a comparison-model improvement, not a
validated unified attention/DRAM implementation.

| Same existing measured cohort | Baseline MAPE | KV supply model MAPE |
|---|---:|---:|
| Full decode layer, 12 points | 3.969% | 2.247% |
| Full decode mean signed error | -3.206% | -0.671% |
| Fixed pooled cohort, 102 points | 3.959% | 3.756% |
| Compute-bound isolated GEMM, 15 points | 8.225% | 8.225% |

These numbers isolate the KV comparison change. They do not include subsequent Claude GEMM fixes.
The 102 points comprise full decode (12), streamed bf16 (20), streamed int8 (20), chained odd GEMMs (35)
and compute-bound isolated kernels (15). `gemm_only` is a diagnostic, not an additional independent cohort.
This is retrospective validation, not a held-out general-accuracy claim. Warm buffer reuse, while-loop
artifacts and isolated copy-dominated small kernels remain outside this cohort.

## Reproduce on washington

From the worktree above:

```bash
M=/home/bogil/MOSAIC/.claude/worktrees/npusim-engine
V=validation/tpu_v5e
A=$V/results/2026-11-11_codex_audit
B=$A/baseline_results
python3 "$V/kv_streams.py" "$B/2026-11-03_colab_v5e1_followup3" \
  --dram "$M/configs/architecture/dram/hbm2e_tpu_v5e_refsb.yaml" \
  --probe "$M/simulator/build/bin/dram_bw_probe" --out "$A/kv_refsb.csv"
python3 "$V/compare_stack.py" "$B/2026-11-03_colab_v5e1_followup3" \
  --dram-tag hbm2e_v5e hbm2e_tpu_v5e_refsb \
  --kv-times "hbm2e_tpu_v5e_refsb=$A/kv_refsb.csv" --output "$A/stack_refsb.json"
python3 "$V/audit_fidelity.py" "$B" --output "$A/baseline.json"
python3 "$V/audit_fidelity.py" "$B" \
  --stack-json "$A/stack_refsb.json" --output "$A/kv_updated.json"
python3 -m unittest discover -s "$V" -p test_kv_streams.py -v
/home/bogil/.claude/jobs/cfc77b0e/tmp/jv/bin/python "$V/test_parse_trace.py"
```

The KV CSV sidecar records probe/YAML/input SHA-256, exact commands and raw probe outputs. Comparison and
cohort JSONs record input hashes. `dram_environment.json` additionally records the actual shared Ramulator library hash. The library must match the existing environment
(MOSAIC worktree `5c16d92`, Ramulator `378827a`); a probe executable hash alone does not freeze its libraries.

Validation completed: five KV tests; original six comparison metric groups unchanged without the new option;
measurement parser and single-/multi-chip script shared-function parity test; real 12-case comparison;
102-point cohort reproduction; separate 64 MiB DRAM extrapolation check.

## Review shared with Claude

The GEMM comparison multiplies small row-tile simulations for large M. Per-program startup must be charged
once, and per-layer setup must not be duplicated for every row tile. Claude owns that correction and the
compute pipeline investigation. A review caught an initial proposal to remove all program setup from the
kernel metric: subtracting copy time from the measurement does not establish that all program setup has
vanished. In particular, copy=0 must preserve equal module/kernel comparisons unless another component is
independently established. Keep the constant in the primary metric; treat any alternative as sensitivity.

The direct stage breakdown shows prologue increasing with strip count, explaining much of the compute
excess. The fixed-cost comparison bug and pipeline modeling issue must be evaluated separately.


## First run after TPU allocation

Use the existing `measure_tpu_v5e8.py` with `--quick --keep-traces` first. Inspect `env.json`: the backend
must be TPU and the expected eight devices must be present before interpreting any performance numbers.
The script intentionally supports CPU smoke tests and reduced device counts, so a successful run alone
is not proof of eight-chip TPU validation. Keep profiler traces; require device timings in collectives
and decode output, and then run the full 2/4/8-chip sweep. Run `--only loopk --keep-traces` separately
for the still-unmeasured loop cost. Do not infer collective latency from single-chip measurements.

Baseline comparison inputs are also preserved under
`results/2026-11-11_codex_audit/baseline_results/`; use that directory as the audit results root if subsequent
core comparisons replace generated CSVs in the original directories.

## Core correction and representative validation

Claude's local commits `7044920` and `7f39fbc` add the `timeline_dram_tiles` option and document the int8
residual. Default is off for other accelerators; TPU v5e profiles enable it. The timeline previously summed
work over DRAM iterations but divided it into only one iteration's tile count. The new option includes those
temporal DRAM iterations in the tile count without changing stage busy totals. For a 512x2048 tile and
1/2/4 column strips, disabled cycles remain 10581/20067/39432 exactly. Enabled cycles are
10581/18390/34682; the four-strip prologue decreases from 6334 to 1584 cycles, while compute stays 33098.
These are structural simulator checks, not errors against a TPU measurement.

The comparison also charges program startup once per program and layer setup only on the first row tile.
Primary kernel error retains the program startup. Removing it is now a separately labeled sensitivity.

The independent 7-case refsb regression uses identical existing device measurements before and after:

| Representative subset | Count | Baseline MAPE | Core-corrected MAPE |
|---|---:|---:|---:|
| bf16, M=1/8, K=N=4096 | 2 | 3.352% | 0.461% |
| int8 weights, same shapes | 2 | 4.893% | **8.018% (regression)** |
| Chained odd GEMMs | 3 | 2.634% | 1.871% |

`core_representative_summary.json` contains all individual predictions and input hashes. These seven cases
are not the 102-point aggregate and cannot establish a new pooled MAPE.

The int8 regression is retained and explicitly reported, not hidden by restoring the known timeline error
only for int8. The existing TPU measurements show a lower byte rate for int8 than bf16; the previous excess
fill cost partly cancelled this residual. Its cause is unresolved. Conversion/dequantization is a hypothesis,
not a measured decomposition, and no fitted penalty was added. The full updated cohort remains pending.

Core validation passed: 617 configurations, 21 Python frontend tests, existing stats/timeline regressions,
the new hand-derived outer-iteration test (250 vs 212.5 cycles with identical busy totals), and
`test_compare_scale.py`. Claude owns the long full-cohort reevaluation and larger whole-GEMM checks, with
separate reporting after completion. No push or branch merge was performed.


### Direct whole-GEMM counterexample — unresolved

The final independent `check_tiling.py 1024,128,1024 -j 2` completed with a tiled estimate of **6.63 us**
and direct whole-GEMM time **5.45 us**, a **+21.74%** difference (see `whole_gemm_small.log`). This is a
validation finding, not a passing accuracy check. It is outside the fixed 102-point cohort, but proves
that the row-scaling shortcut has no general 1.5% accuracy bound. Repeated fill/drain and different reuse
across row boundaries are candidates, not established causes. The core DRAM tile-count correction and
once-only setup accounting do not resolve this approximation. No speculative overlap subtraction was added.
Claude was asked to warn on row-tiled estimates and retain this case as a follow-up reproduction.

The current result is therefore a tested structural correction plus a more explicit comparison pipeline,
with quantified bf16/chain and KV-model improvements, an int8 regression, and a known row-scaling limitation.
It is not a certification of general NPUsim accuracy. Full-cohort core reevaluation and multi-chip TPU
validation remain separate pending work.
