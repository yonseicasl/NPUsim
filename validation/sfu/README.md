# SFU validation

Regression coverage for the Special Function Unit (`components/sfu.{h,cc}`), implementing
the test strategy of [plan/plan_sfu.md](../../plan/plan_sfu.md).

## What the SFU models

A per-chip post-processing component between the final accumulator value and the
output-format cast (`placement = post_accumulator`, `fusion = fused`, `queue_depth = 1`
serial contract). It charges cycles, internal ingress/egress traffic and dynamic/static
energy for the activation applied to each network-valid output element exactly once —
after all reduction, after repetition scaling, with mapping padding excluded. Nebula's
`forward()` remains the only functional owner of activation values; the SFU's pure
evaluators exist for these tests only.

- Element-wise operations: `linear` (bypass), `relu`, `leaky`, `elu`, `sigmoid`
  (= Nebula `logistic`), `tanh`, `hsigmoid`, `hswish`, plus the LLM operations `gelu` and
  `silu`/`swish` (SFU-level; Nebula's frontend has no enum for them yet).
- Standalone softmax (Phase 7): the multi-pass microprogram
  `max -> subtract -> exp -> sum -> reciprocal -> normalize` built from the micro-op
  profiles `vmax`, `vadd`, `exp`, `recip`, `vmul`; rows distribute across every chip's
  SFU, and the operand tensor's streaming is charged per `softmax_operand_residency`
  (`dram`: DRAM device + off-chip link + GLB staging/feed ports, the default matching the
  simulator's materialize-every-layer flow; `glb`: on-chip retained, GLB feed/result
  ports only, with a capacity fail-fast).
- Streaming overlap (Phase 5): `queue_depth = 1` keeps the serial contract (the SFU
  window is appended to the critical path -- bit-identical to the pre-streaming model);
  `queue_depth >= 2` adds the SFU as a sixth per-tile stage of the layer's pipeline
  timeline behind a bounded output-tile queue. Final outputs release on the LAST tile of
  each reduction group (`datatype_repetitions[OUTPUT]` of the tile clock), a fast SFU
  hides behind the producer except fill/drain, and a slow one back-pressures the PE stage
  (reported as producer queue stall). Hand-calculated cases live in
  `unittest/stats_timeline_test.cc` (`run_sfu_streaming_cases`).
- final_output_tile events (Phase 2): each multi-chip -> DRAM output commit (the
  reduction-complete, once-per-element boundary that also charges the final cast) is one
  `final_output_tile` event; the layer's SFU work arrives in
  `commits x output-repetitions` serial invocations, each with its own setup and
  pipeline fill. The identity gate -- committed elements must reproduce the mapped
  output volume exactly -- guards the hookup: on a mismatch the model falls back to the
  layer-granular single invocation and the report says so (`Output tile commits` line).
  Spill/reload independence is structural (psum traffic never crosses that boundary).
- RTL calibration (Phase 8): `nvdla_small_sfu.cfg` carries an RTL-measured nv_small SDP
  timing profile (core II = 1.0000, X1-ALU fill = 151 cycles) with declared
  `profile_reference` provenance; see [SDP_CALIBRATION.md](SDP_CALIBRATION.md) for the
  measurements, the runner (`run_sdp_rtl_traces.sh`, `NVDLA_PROJECT`-parameterized),
  the extractor (`calibrate_sdp.py <project>`), and the root-cause record of the
  `ew_*` trace failures on nv_small (NO LUT/EW engine — spec
  `SDP_LUT_DISABLE`/`SDP_EW_DISABLE`). A profile without declared provenance is
  reported `Timing provenance : NOT DECLARED -- ... not calibration-grade`.
- **LUT-class calibration (2026-08-28)**: the smallest LUT-enabled spec,
  `nv_small_256_full`, was brought up on the same harness
  (`nvdla_tmake_driver.py` + `gen_verilator_filelist.py`; local patches recorded in
  SDP_CALIBRATION.md) and all official `ew_*` traces **PASS their goldens** there —
  measuring the SDP LUT engine at **II = 1 elem/cycle exactly and fill = 156 cycles,
  LUT-table-content independent** (exp ≡ linear-LE ≡ linear-LO at 453 cycles / 297
  elements), with the X1 path spec-exact at 2 elem/cycle. A ready-to-declare
  `[sfu]` profile for sigmoid/tanh/exp/recip is in SDP_CALIBRATION.md.
- Architecture allowlist: `[sfu] supported_ops = linear,relu,...` makes out-of-list
  operations fail fast (up-front, before any layer simulates) — a missing execution
  unit is an architecture fact, not an uncalibrated cost. Standalone softmax requires
  its whole micro-op set to be listed. Linear is always implied.
- Softmax groups: Nebula now honors `[softmax] groups` (one independent normalization
  per batch x group span; fail-fast on non-dividing groups), and the SFU model derives
  rows = batch x groups, length = size/groups from the same fields
  (`gemm256_g4` fixture: 1024 x 64 -> 28,812-cycle multi-pass window, hand-checked).
- Softmax operand streaming now also charges DRAM open-page row activations when the
  [dram] row model is declared (same resolution as dram_t: tRC when calibrated, flat
  row_miss otherwise; bank parallelism helps latency, never energy), and the event
  count feeds the unpriced-active check for `row_miss_energy`.

## Running

```sh
bash validation/sfu/run_sfu_validation.sh
```

Builds `sfu_test.cc` against `components/sfu.cc` (no simulator build required) and runs:

| group | checks |
|---|---|
| bypass/boundary | linear and 0-element invocations charge no element, cycle, energy or event |
| lane tail | `N < lanes`, `N = lanes`, `N = lanes + 1` chunk/cycle/tail-utilization exactness |
| formula | `setup + latency + (chunks-1) x II`; relu/leaky profile independence |
| energy linearity | doubling `sfu_op_energy_relu` doubles only the op-energy axis |
| parallelism | intra-chip units: latency = busiest unit, work/energy conserved |
| event identity | chip-partition element shares sum to `B x K x P x Q`; one op per element |
| softmax | full hand calculation (window, ops, scratchpad reads/writes, per-pass energy); empty/length-1 boundaries |
| unpriced/uncalibrated | active ops with missing energy keys surface as unpriced events; defaulted timing profiles are flagged UNCALIBRATED |
| mapping | Nebula activation names resolve (`logistic`→sigmoid, `swish`→silu); unsupported names refuse — no silent ReLU fallback |
| functional | evaluators vs reference math (leaky slope 0.1, relu6-based hard ops, GELU tanh-approx vs erf, softmax sums to 1) |
| fail-fast | 0 lanes/units, `queue_depth != 1`, unknown placement/fusion, 0 II, profiles on the linear bypass, `strict_profiles` with a missing active profile — all abort (subprocess modes) |

## End-to-end fixtures

- `configs/accelerators/gemmini_sfu.cfg` — gemmini baseline + a fully profiled `[sfu]`
  section (PLACEHOLDER, uncalibrated costs; Phase-8 targets NVDLA SDP for calibration).
  Kept separate from `gemmini.cfg` so SFU-disabled baselines stay numerically identical.
- `configs/accelerators/gemmini_sfu_unpriced.cfg` — same, minus `sfu_op_energy_relu`; a
  relu run must report the unpriced active SFU event.
- `configs/accelerators/gemmini_sfu_stream.cfg` — same as `gemmini_sfu.cfg` with
  `queue_depth = 2`: the relu layer's 4,096 SFU cycles collapse to a 16-cycle drain tail
  on the critical path (4,080 hidden), while the depth-1 twin still exposes all 4,096.

Useful runs (`cd models; NPUSIM_CONFIG_ROOT=../configs ./model run <acc> <net> ws`):

- `gemmini_sfu` × `gemm256` — linear bypass (zero-cost SFU event) + standalone softmax
  layer on the SFU multi-pass model.
- `gemmini_sfu` × `gemm256_relu` — fused relu: +4096 serial SFU cycles on the critical
  path, DRAM traffic identical to the `gemmini` run (fused invariant).
- `gemmini` × `gemm256_relu` — no `[sfu]`: legacy numbers unchanged, report carries
  `Activation scope : NOT MODELED`.

## Contract notes

- Runs without `[sfu]` keep their numeric baselines; the energy breakdown keeps six rows
  and the undercount denominator stays 6 (`energy_cost_schema_t::components_in_scope()`).
- With `[sfu]`, the SFU is a seventh energy component; missing `sfu_read/write_energy`
  mark it PARTIAL and any active-but-unpriced SFU event blocks absolute energy/power
  qualification like every other unpriced event.
- `pe_t::activation()` (the old hardcoded ReLU) is deprecated and aborts if called.
- Standalone softmax scope: SFU-internal cycles/traffic/energy only; operand streaming
  between the memory hierarchy and the SFU is stated as not-included in the layer report
  (materialized-tensor accounting is a later phase).
