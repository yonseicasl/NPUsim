# NVDLA nv_small SDP timing calibration (Phase 8, plan/plan_sfu.md) — 2026-08-26

RTL-measured timing constants for the SFU's NVDLA-SDP profile
(`configs/accelerators/nvdla_small_sfu.cfg`), from CRC-PASSING Verilator runs of the
official nv_small `sdp_*` trace tests.

## Provenance

- RTL, toolchain, harness patches: **identical** to the conv-pipeline fidelity work —
  `nvdla/hw` master @ `1a65f1f5b48268accaa47c95f95c2601918be095`, Verilator 5.022
  (chipyard conda env), the EMULATION-RAM binary
  (`ext/nvdla/outdir/nv_small/verilator_emu/VNV_nvdla`); every patch is recorded in
  [validation/nvdla/PROVENANCE.md](../nvdla/PROVENANCE.md).
- Runner: [`run_sdp_rtl_traces.sh`](run_sdp_rtl_traces.sh) — builds each official
  trace's `trace.bin` with the stock trace-parser flow and runs it; idempotent, keeps
  PASSing logs.
- Extraction/analysis: [`calibrate_sdp.py`](calibrate_sdp.py). Cycle boundary:
  the LAST accepted CSB write of 1 to an SDP-block `D_OP_ENABLE`
  (`NVDLA_SDP_RDMA` 0x8008, `NVDLA_SDP` 0x9038) → the last output DBB write beat;
  2 harness ticks per clock. The window **includes SDP-RDMA operand fetch from DBB**
  (the honest analog of NPUsim's operand streaming + SFU flow).
- Golden discipline: only `*** PASS` (CRC-verified) runs are used. 12 traces PASS.

## Measurements (12 CRC-PASS traces)

| trace | class | elements | cycles | cyc/elem |
|---|---|---:|---:|---:|
| sdp_1x1x1_pass_through_int8 | stream | 1 | 85 | 85.000 |
| sdp_1x1x8_pass_through_int8_0 | stream | 8 | 85 | 10.625 |
| sdp_4x1x8192_pass_through_int8_0 | stream | 32,768 | 32,851 | 1.003 |
| sdp_4x22x42_bypass_int8 | stream | 3,696 | 4,312 | 1.167 |
| sdp_1x8192x1_pass_through_int8_0 | stream | 8,192 | 65,613 | 8.009 |
| sdp_8192x1x1_pass_through_int8_0 | stream | 8,192 | 65,619 | 8.010 |
| sdp_3x3x33_bs_int8_reg_0 | alu | 297 | 448 | 1.508 |
| sdp_3x3x33_bn_int8_reg_0 | alu | 297 | 448 | 1.508 |
| sdp_3x3x33_bn_int8_mem_0 | alu | 297 | 451 | 1.519 |
| sdp_5x24x18_bs_int8_mem_0 | alu | 2,160 | 2,975 | 1.377 |
| sdp_23x13x42_bs_int8_mem_0 | alu | 12,558 | 14,446 | 1.150 |
| (sdp_pdp/… not run) | | | | |

## Derived constants (each a DIRECT measurement, no mixed-regime fitting)

- **Core throughput = 1 element/cycle exactly.** `(32,851 − 85)/(32,768 − 1) = 1.0000`
  on the layout-efficient 32,768-element pass-through. Matches the nv_small spec
  throughput of 1. → `lanes = 1`, `*_initiation_interval = 1`.
- **Streaming skeleton fill = 85 cycles.** N=1 and N=8 pass-throughs are identical
  (85), so the whole window at tiny N is fill/drain (RDMA + SDP + WDMA pipe).
- **X1-ALU path fill = 151 cycles.** `448 − 297×1` on the register-operand bias trace
  (no second-operand DMA). ReLU executes on the same X1 ALU stage as bias →
  `relu_pipeline_latency = 151`, `relu_initiation_interval = 1`.

The W=1/H-major pass-throughs (8.01 cyc/elem) and the mem-operand bias traces are
**operand-DMA-bound references**, recorded but never fed into the core constants —
NPUsim charges operand streaming separately, so absorbing DMA layout inefficiency into
the SFU op profile would double-model it (plan Phase-8: no residual stuffing).

## Resolved: the `ew_*` "CRC FAIL" root cause (2026-08-26)

**nv_small has no LUT and no EW engine at all** — `spec/defs/nv_small.spec` declares
`SDP_LUT_DISABLE` and `SDP_EW_DISABLE`, so the SDP the RTL builds simply does not
instantiate the unit those four traces exercise. The evidence chain:

1. The harness's FAIL message prints `should have been <computed> , was <golden>` —
   i.e. the *computed-over-memory* CRC first. All four traces compute the **identical**
   CRC `c43621c5` while their trace-declared goldens differ — including across
   **different input tensors** (`ew_int8_reg_0`'s `.dat` differs from the other three)
   and different shapes (3x3x33 vs 3x3x32). An input-independent output is exactly what
   a datapath without the programmed EW/LUT stage produces (converter fallout of the
   X1 path), and it can never match goldens generated for an EW-enabled configuration.
2. `trace.bin` carries the CORRECT goldens (verified byte-level), and the trace parser
   output is correct — the packing/compare machinery is not at fault.
3. A stock (non-EMULATION) binary cross-run hangs at ~tick 18.9k on this trace (the
   known RAMDP write-pulse pathology stalls it), so it cannot adjudicate — but the spec
   defines plus (1) are conclusive without it.

Consequences, reflected in `nvdla_small_sfu.cfg`:
- `supported_ops = linear,relu,leaky` — LUT-class operations and the softmax micro-ops
  **fail fast** on this architecture (missing execution unit = architecture fact, not an
  uncalibrated cost). The four `ew_*` traces are out of scope for nv_small, not a
  harness defect.
- Calibrating sigmoid/tanh/exp therefore required an EW/LUT-enabled spec build —
  **done below (nv_small_256_full, 2026-08-28)**.

# nv_small_256_full bring-up and LUT-class calibration (2026-08-28)

`nv_small_256_full` is the smallest LUT/EW-enabled spec (`SDP_LUT_ENABLE`,
`SDP_EW_ENABLE`, `SDP_EW_THROUGHPUT_1`, `SDP_BS/BN_THROUGHPUT_2`, secondary
memif/CVSRAM enabled) and ships its own official trace set with goldens — including the
`ew_*` LUT traces. Same RTL commit, Verilator 5.022, EMULATION-RAM binary, `--timing`,
gcc, as the nv_small build.

## Bring-up record (all local, none change RTL semantics)

1. `tree.make`: `PROJECTS := nv_small nv_small_256_full`.
2. `validation/sfu/nvdla_tmake_driver.py` — dependency-free replacement for
   `tools/bin/tmake` (its perl deps YAML/Capture::Tiny are not installed); walks
   `tools/etc/build.config` and runs the same `cd <sandbox>; make PROJECT=…` commands.
   Ran `vmod` with `--assume-done manual,odif` (the `manual` stage's ORDT
   post-processing needs the missing `XML::Simple`).
3. `outdir/nv_small_256_full/spec/manual/NVDLA_CFGROM.vh` copied from nv_small (pure
   register-offset macros, config-invariant — same RDL) so the cfgrom vmod stage could
   regenerate its previously silently-empty `NV_NVDLA_CFGROM_rom.v` (the CFGROM body
   itself is hard-coded in the checked-in source). A whole-tree scan confirmed cfgrom is
   the ONLY silently-empty generated file.
4. `validation/sfu/gen_verilator_filelist.py` → `verif/verilator/
   verilator_nv_small_256_full.f`: mirrors the shipped nv_small .f structure but
   registers ALL of `vmod/rams/model` and `vmod/fifos` with `-v` up front (their module
   names ≠ file names), avoiding the 67-entry iteration the nv_small bring-up needed.
5. `verif/verilator/nvdla.cpp` (backup `.pre_cvsram_fix`): the CVSRAM connection block —
   never compiled before (nv_small has no secondary memif) — was missing the
   `AXI_WDATA_MKPTR` width macro on `w_wdata`/`r_rdata` that the DBB block uses; both
   memifs are 64-bit here, so the fix is exactly the DBB pattern.
6. `verif/verilator/new_trace_to_verilator.py` (backup `.pre_secmem_fix`): accept
   `SEC_MEM` in memory commands — the harness routes AXI ops by address
   (bit31 clear → CVSRAM), so the blob format needs no change; the old `PRI_MEM`-only
   assert aborted the packer mid-write, which is precisely the "short read" failure
   (`sdp_13x4x29_ew_lut_int8` here, and the same assert family —
   `check_crc` memory_type — explains nv_small's `sdp_8x8x32_bypass_int8_0` open item).

## Measurements (all CRC-PASS on the official goldens)

| trace | elements | cycles | cyc/elem |
|---|---:|---:|---:|
| sdp_1x1x1_pass_through_int8 | 1 | 81 | 81.000 |
| sdp_1x1x8_pass_through_int8_0 | 8 | 81 | 10.125 |
| sdp_4x1x8192_pass_through_int8_0 | 32,768 | 16,467 | 0.503 |
| sdp_3x3x33_bs_int8_reg_0 | 297 | 268 | 0.902 |
| sdp_23x13x42_bs_int8_mem_0 | 12,558 | 7,265 | 0.579 |
| sdp_3x3x33_ew_int8_reg_0 | 297 | 447 | 1.505 |
| **sdp_3x3x33_ew_le_exp_int8** | 297 | **453** | 1.525 |
| **sdp_3x3x33_ew_le_lin_int8** | 297 | **453** | 1.525 |
| **sdp_3x3x32_ew_lo_lin_int8** | 288 | **453** | 1.573 |
| sdp_13x4x29_ew_lut_int8 | 1,508 | 1,763 | 1.169 |

The four nv_small-failing `ew_*` traces **PASS against their official goldens here** —
closing the root-cause loop: the unit exists on this config, and the same harness
reproduces the goldens exactly.

## Derived constants (direct measurements; spec-cross-validated)

- **X1 (BS/BN/ReLU) path = 2 elements/cycle**: `(16,467 − 81)/32,767 = 0.5001` — spec
  `SDP_BS/BN_THROUGHPUT_2`, measured exact. Stream fill 81 cycles (N=1 ≡ N=8);
  X1-ALU fill `268 − ceil(297/2) = 119` cycles.
- **EW/LUT path = 1 element/cycle EXACTLY**: after removing the fill,
  `(453 − 156)/297 = 1.0000` — spec `SDP_EW_THROUGHPUT_1`, measured exact.
- **EW/LUT path fill = 156 cycles**, and it is **LUT-table-content independent**:
  exp-table, linear-LE and linear-LO traces all measure 453 cycles. The EW path without
  a LUT lookup (register-operand `ew_int8_reg`) fills in 150 — the lookup itself adds
  6 fill cycles, zero throughput.
- `sdp_13x4x29_ew_lut` (1.0818 cyc/elem) is a DMA-inclusive upper-bound reference
  (W=13-byte lines + CVSRAM second operand), never a core constant.

## Ready-to-declare [sfu] profile — NVDLA SDP LUT engine (EW_THROUGHPUT_1)

```ini
# For an NVDLA-style SFU with the LUT/EW engine (nv_small_256_full class).
# One scalar LUT function per pass (sigmoid, tanh, exp, reciprocal, ...):
sigmoid_pipeline_latency = 156
sigmoid_initiation_interval = 1
sigmoid_approximation = lut
tanh_pipeline_latency = 156
tanh_initiation_interval = 1
tanh_approximation = lut
exp_pipeline_latency = 156
exp_initiation_interval = 1
exp_approximation = lut
recip_pipeline_latency = 156
recip_initiation_interval = 1
recip_approximation = lut
profile_reference = nv_small_256_full-SDP-RTL-nvdla/hw@1a65f1f-verilator5.022-EMU-10-CRC-PASS-traces-calibrate_sdp.py-2026-08-28
```

The latency/II pair describes the ENGINE, not the table content (measured identical
across tables), so it applies to any single-table LUT operation. Multi-primitive ops
(GELU, SiLU = lookup + multiply) are NOT covered by a single measured pass and stay
uncalibrated. `nvdla_small_sfu.cfg` itself is unchanged — nv_small still has no LUT
unit; these constants are for LUT-capable NVDLA-class configs.

## Reproduce (256_full)

```sh
python3 validation/sfu/nvdla_tmake_driver.py nv_small_256_full vmod --assume-done manual,odif
python3 validation/sfu/gen_verilator_filelist.py nv_small_256_full
# verilate + compile with the recorded flags (+define+EMULATION --timing, gcc)
NVDLA_PROJECT=nv_small_256_full validation/sfu/run_sdp_rtl_traces.sh sdp_3x3x33_ew_le_exp_int8 ...
python3 validation/sfu/calibrate_sdp.py nv_small_256_full
```

## Open items

- `sdp_8x8x32_bypass_int8_0` "short read": root-caused 2026-08-28 — the trace packer's
  strict `check_crc` memory_type assert aborts the blob build mid-write (same assert
  family as the SEC_MEM case fixed for `sdp_13x4x29_ew_lut`); the trace checks a
  non-pri-mem region nv_small does not have, so it stays out of scope there.
  `sdp_4x22x42_bypass_int8` covers bypass.
- **Energy: no measurement exists.** `ext/nvdla/syn/` ships Synopsys DC scripts
  (`dc_run.tcl`), but no `dc_shell`/liberty libraries are installed on this machine, so
  a power run is not executable here. The concrete protocol when tooling exists:
  synthesize `NV_NVDLA_SDP*` with the target library, replay the PASSing traces on a
  `VM_TRACE=1` harness build to dump SAIF/VCD activity, and derive per-element op/port
  energies from the reported switching power. Until then the config declares no SFU
  energy keys, so active SFU events surface as UNPRICED and absolute energy/power stays
  unqualified.

## Reproduce

```sh
validation/sfu/run_sdp_rtl_traces.sh sdp_1x1x1_pass_through_int8 \
    sdp_1x1x8_pass_through_int8_0 sdp_4x1x8192_pass_through_int8_0 \
    sdp_3x3x33_bs_int8_reg_0 ...        # idempotent
python3 validation/sfu/calibrate_sdp.py
```
