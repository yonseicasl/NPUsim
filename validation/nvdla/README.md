# NVDLA RTL fidelity — status (2026-08-25)

Supersedes the 2026-08-20 first-result state (kept in git history and in
`.claude/worktrees/nvdla-fidelity`). Executes the plan
`assessment/NVDLA_RTL_FIDELITY_PLAN_2026-08-20.md` through Phase 1-2 data
collection: the multi-Atomic-C CRC blocker is ROOT-CAUSED AND FIXED, all 13
official nv_small `dc_*` traces are CRC-verified golden, and the first
full-suite NPUsim comparison exists on the current HEAD.

## Result (2026-08-25, `python3 validation/nvdla/compare_full.py`)

```
cycle HOLDOUT (3 traces):   compute-schedule MAPE 0.55%, max 0.75%
                            (calibrated model; parameters never saw these traces)
critical path (all 9):      -11.7%..+22.2%, MAPE 11.7% (informational)
calibration set (5):        1x1x8 (drain constant, -15.7%) + 4 stripe-fit traces
                            (+0.4/-0.1/-0.0/+3.0%)
out of scope (4x1x8192):    -64.2% informational -- 1x1-kernel C-heavy stride-3
                            regime, outside the stripe model's declared scope
DRAM traffic: input/weight/output ALL EXACT vs the declared encoding on every
              comparable trace (sub-beat outputs round up to one 8B DBB beat).
              RTL physical DBB traffic obeys the closed-form atomic-padding law
              exactly on all 14 traces.
```

Files:
- `golden_rtl_full.csv` — 14 CRC-passing traces: cycles + per-stream DBB bytes
  (regenerate: `extract_golden.py --write`; refuses non-PASS runs)
- `gen_nvdla_workload.py` — trace registers -> NPUsim network/mapping/weights
- `compare_full.py` — the comparison; prints per-trace verdicts
- `PROVENANCE.md` — RTL commit, toolchain, patches, root-cause record

## What was wrong before (and its status)

The 2026-08-20 state had every C>8 trace failing CRC with all-zero output; the
root cause (PROVENANCE.md 2026-08-25 section) is the TSMC `RAMDP_*_E2`
behavioral models' `assign #1` write-pulse choreography, which never lands a
write under Verilator; the MCIF read-response FIFOs sit on exactly those
models, so DBB read DATA (weights first — they fill the FIFO past its flop
bypass) arrived as zeros while every handshake and every cycle count stayed
correct. `+define+EMULATION` swaps in NVIDIA's own simple behavioral RAM
models and fixes it with no RTL change; tick counts are bit-identical.

Consequences worth keeping in mind:
- Cycle counts extracted before the fix were TIMING-correct (proven by
  identical tick counts), but golden status still requires the CRC PASS, so
  only EMULATION-binary runs feed `golden_rtl_full.csv`.
- The failure boundary was never "C > 8"; it was "enough in-flight beats".

## Energy (2026-08-26)

`nvdla_small_cacti22.cfg` is the energy-calibrated variant (timing identical to
nvdla_small.cfg). Component provenance and every derivation live in its header and
validation/phase4/PROVENANCE.md; gate P4-10 pins the declared costs to the CACTI/
DRAMsim3 reference rows and pins the REFUSAL: with MAC/registers/tree/control
undeclared, the fixture is denied a wattage, its total says ESTIMATED, and the
undercount is visible per component -- the same honest direction as gemmini_cacti22.
Calibrated: CBUF SRAM (CACTI 7 @22nm, RTL macro geometry), CACC assembly RAM
(per-byte, via the accumulator key), DBBIF (DRAMsim3 DDR3-1600). The relative
breakdown behaves physically: the DMA-heavy 1x1-kernel trace is 72.6%% DRAM while
the K=270 compute-heavy trace is 8.4%% DRAM.

## Known-open items (each names its owner)

1. **Input DRAM traffic axis — CLOSED 2026-08-26.** The halo contract now folds
   legacy-row R/S loops into the union extent (full_h = (P_tot-1)*stride+R_tot)
   and stats applies it in BOTH directions: coalescing down stays gated on the
   ring working set fitting the GLB (a retention claim), while raising the
   repetition-scaled undercount UP to the union is unconditional (the union is
   a lower bound under any loop order and buffer size). All 9 comparable traces
   now read exact against the declared encoding; the Eyeriss/Gemmini baselines
   are bit-identical (their legacy rows carry no filter loops). Hand identity
   pinned in unittest (`validate_input_halo_contract`, legacy-R/S case:
   replicated 10,752 < unique 16,384, matching the dc_13x15x64 measurement).
2. **Cycle model for memory-bound traces — MODELED 2026-08-26.** The fixed
   pipeline constants are declared as `[spatial_arch] layer_setup_cycle = 137`:
   launch 28 (enable->first read request; exactly constant on all 13 conv
   traces) + DBB latency 32 (exactly constant on all 14 AND independently
   derivable from the harness's AXI_R_LATENCY=32) + conv-pipe drain 77
   (calibrated on dc_1x1x8, the only trace whose compute tail is zero -- that
   trace moves from the holdout to a labeled calibration point). Measurement:
   `measure_pipeline_constants.py`. dc_1x1x8 critical path: -76.9% -> -18.3%;
   the residual is composition, not magnitude -- NPUsim overlaps the setup with
   the DRAM stream where the RTL runs them serially, and shrinking the constant
   to absorb that would be curve-fitting, so it stays as measured.
   dc_4x1x8192 remains informational: its RTL time is per-stripe-overhead
   bound (item 3), not latency bound.
3. **Per-stripe overhead — CALIBRATED 2026-08-26.** Model: one fixed bubble per
   legacy (C-group, r, s) reduction slice per output row,
   `stripe_transition_cycle = 4.356` (`[spatial_arch]`), count =
   `mapping_table_t::stripe_transition_count()` = legacy C*R*S*P, added once per
   layer beside layer setup. Kg-independence is MEASURED (the 32x26x76 K16/K270
   twin: 17x the weight blocks, near-identical residual 22,210 vs 22,898 --
   which kills every per-weight-block form; both twins therefore sit in the
   calibration set, since they selected the form). B fitted by least squares on
   the 4-trace calibration set only; the 3-trace holdout reads MAPE 0.55% / max
   0.75% (calibrate_stripe_overhead.py, reproduced by the simulator itself via
   compare_full.py). Scope: stride-1/dilation-1 direct conv -- the non-unit
   stride/dilation traces put 9.4-29.4 cycles on the same axis and stay outside
   the claim, as does the degenerate C-heavy 1x1-kernel shape (-64.2%
   informational).
4. **Excluded shapes.** dilation (2 traces), asymmetric stride (1), R>H (1):
   NPUsim's mapping table cannot express them; listed by
   `gen_nvdla_workload.py` with reasons.
5. **CACC psum residency — MODELED 2026-08-26.** Two corrections together:
   (a) `psum_must_leave_array()` is order-aware for SAME-LEVEL reduction/output
   coexistence -- the within-level order is not unknown, it is the declared
   array<->GLB stationary type the scheduler itself sequences the legacy-GLB
   offsets by; OUTPUT_STATIONARY completes each output tile's reduction before
   advancing, so the psum stays at the array edge (nv_small: in CACC) and never
   crosses the GLB. The ACROSS-level case (Eyeriss conv3) stays unconditional --
   nesting, not ordering. Hand identities pinned in unittest.
   (b) nvdla_small.cfg's `pe_stationary` corrected weight->output stationary:
   nv_small's CSC completes all (r,s,c-group) slices per output stripe before
   moving on; the per-lane weight hold WITHIN a stripe is `mac_stationary`, a
   different level. Result: psum GLB round trips are GONE (8x16x128: psum legs
   0 loads / 476 stores == exactly one final write-out per output tile, the T12
   identity), with compute-schedule, DRAM traffic, and both calibrations
   bit-invariant.

   **Array streaming granularity — MODELED 2026-08-26.** The critical-path
   inflation was five stacked artifacts of charging streams one op at a time on
   serialized shared media the RTL does not have. Each got its own structural
   declaration (every knob defaults to the historical behavior; Eyeriss/Gemmini
   are bit-identical):
     - intra-array wt bus width 512 b + [spatial_arch] line_size 64:512:64 (the
       CSC weight bus moves a 64 B atomic-kernel slice per cycle);
     - writeback_injection_parallel: per-PE write-backs are dedicated wires into
       the CMAC tree/mac2accu bus -- time overlaps, energy/transactions per wire
       unchanged;
     - array_fabric_separate / [separate] fabric_separate: dat, wt and accu ride
       separate buses; the link/overlap axes combine by MAX (T10's rule family);
     - operand_streams_pipelined: MAC operand registers load in lockstep with
       compute (RTL sustains 1 op/lane/cycle) -- one suppression point zeroes
       only the exposed-cycle views of the MAC<->LB streams;
     - GLB converted shared -> separate (D_BANK partitions dat/wt into disjoint
       banks read concurrently; static 68/52/8 KB split approximates the
       per-layer D_BANK programming) and the inherited-from-Eyeriss weight
       bypass corrected to CBUF-resident (the measured fetch-once + per-stripe
       re-read behavior); streams_pipelined amortizes the per-transfer packet
       fill across the stream.
   Ladder measured on 8x16x128: 5,003,743 -> 822,548 -> 221,282 -> 152,738 ->
   93,947 cycles (53x). Critical path across all 9 traces: -11.7%..+22.2%,
   MAPE 11.7% (informational; systematically positive residual = remaining
   un-modeled overlaps). Compute-schedule, DRAM traffic, and both calibrations
   are bit-invariant throughout.

## Reproduction

```
# goldens (needs the EMULATION binary, PROVENANCE.md):
cd ext/nvdla/outdir/nv_small/verilator_emu/run_<trace> && ../VNV_nvdla trace.bin
python3 validation/nvdla/extract_golden.py --write

# NPUsim side:
python3 validation/nvdla/gen_nvdla_workload.py
./npusim.sh run nvdla_small nvdla_<trace> ws
python3 validation/nvdla/compare_full.py
```
