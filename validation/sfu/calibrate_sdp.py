#!/usr/bin/env python3
"""Phase-8 (plan_sfu.md): SDP timing calibration from nv_small RTL trace runs.

Extracts the SDP-active cycle window from CRC-PASSING Verilator runs of the official
nv_small `sdp_*` trace tests (run via run_sdp_rtl_traces.sh on the EMULATION binary --
see validation/nvdla/PROVENANCE.md for the RTL commit, toolchain and patches), then fits

    cycles = fixed + slope x elements

per SDP operating class. The mapping onto the NPUsim [sfu] model is exact:
with lanes = L and setup_cycle = 0, the model's window for one invocation is
pipeline_latency + (ceil(N/L) - 1) x II, so `fixed` calibrates <op>_pipeline_latency and
`slope x L` calibrates <op>_initiation_interval.

Cycle boundary (SDP-only traces; same tick convention as validation/nvdla):
  start = the LAST accepted CSB write of 1 to an SDP-block D_OP_ENABLE
          (NVDLA_SDP_RDMA 0x8008 -> CSB word 0x2002, NVDLA_SDP 0x9038 -> word 0x240e)
  end   = the last output DBB write beat ("DBB: write, last tick")
  cycles = (end_tick - start_tick) / 2      (the harness runs 2 ticks per clock)

This window INCLUDES the SDP-RDMA operand fetch from DBB, which is the honest analog of
the NPUsim standalone-activation flow (operand streaming + SFU); the report labels the
fixed term accordingly. Only PASSing runs are used; FAILing runs (the four nv_small
`ew_*` LUT traces as of 2026-08-26) are listed as excluded, never fitted.
"""
import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
# Optional argv[1] selects the NVDLA project outdir (default nv_small). The LUT-class
# analysis only produces numbers on an EW/LUT-enabled config (nv_small_256_full):
# nv_small's SDP has no such unit (spec SDP_LUT_DISABLE/SDP_EW_DISABLE).
PROJECT = sys.argv[1] if len(sys.argv) > 1 else "nv_small"
EMU = os.path.join(ROOT, f"ext/nvdla/outdir/{PROJECT}/verilator_emu")

# Per-config SDP datapath widths, from spec/defs/<project>.spec -- the fill extraction
# must subtract the STREAMING term at the config's own throughput, or a 2-wide pipe
# yields a nonsense negative fill.
PROJECT_INFO = {
    "nv_small":           {"alu_throughput": 1, "ew_throughput": None},  # no EW/LUT unit
    "nv_small_256_full":  {"alu_throughput": 2, "ew_throughput": 1},
}

OP_ENABLE_WORDS = {0x8008//4, 0x9038//4}   # SDP_RDMA, SDP

def classify(name):
    if "pass_through" in name or "bypass" in name:
        return "stream"          # streaming skeleton: no arithmetic engaged
    if "_bs_" in name or "_bn_" in name:
        return "alu"             # X1 ALU path (bias/batch-norm; same path ReLU uses)
    if "_ew_" in name:
        return "lut"             # EW/LUT engine (exp/sigmoid/tanh class)
    return "other"

def elements(name):
    m = re.search(r"sdp_(\d+)x(\d+)x(\d+)_", name + "_")
    if not m:
        return None
    w, h, c = (int(g) for g in m.groups())
    return w*h*c

def extract(run_dir):
    log = os.path.join(run_dir, "run.log")
    if not os.path.exists(log):
        return None
    text = open(log, errors="replace").read()
    if "*** PASS" not in text:
        return ("NOT_PASS", None)
    start = None
    for m in re.finditer(r"\((\d+)\) write to nvdla: addr ([0-9a-f]+), data 0*1\b", text):
        if int(m.group(2), 16) in OP_ENABLE_WORDS:
            start = int(m.group(1))
    end = None
    for m in re.finditer(r"\((\d+)\) DBB: write, last tick", text):
        end = int(m.group(1))
    if start is None or end is None:
        return ("NO_BOUNDARY", None)
    return ("PASS", (end - start)//2)

def width(name):
    m = re.search(r"sdp_(\d+)x(\d+)x(\d+)_", name + "_")
    return int(m.group(1)) if m else None

def main():
    rows = {}
    excluded = []
    for entry in sorted(os.listdir(EMU)):
        if not entry.startswith("run_sdp_"):
            continue
        name = entry[len("run_"):]
        result = extract(os.path.join(EMU, entry))
        if result is None:
            continue
        status, cycles = result
        count = elements(name)
        if status != "PASS" or count is None:
            excluded.append((name, status))
            continue
        rows[name] = (classify(name), count, cycles)

    print(f"{'trace':44s} {'class':7s} {'elements':>9s} {'cycles':>9s} {'cyc/elem':>9s}")
    for name, (kind, count, cycles) in rows.items():
        print(f"{name:44s} {kind:7s} {count:9d} {cycles:9d} {cycles/count:9.3f}")
    for name, status in excluded:
        print(f"{name:44s} EXCLUDED ({status})")
    print()

    # Regime-aware reading -- a single line does NOT fit the stream class, because the
    # W=1 / H-major surfaces are DBB-layout-degenerate (one element per line, so the
    # window is operand-DMA-bound, not SDP-core-bound). Each calibrated constant below
    # is one direct measurement, never a residual of a mixed fit (plan_sfu.md Phase-8:
    # no residual stuffing).
    def get(name):
        return rows.get(name, (None, None, None))

    _, n1, c1 = get("sdp_1x1x1_pass_through_int8")
    _, n8, c8 = get("sdp_1x1x8_pass_through_int8_0")
    _, nbig, cbig = get("sdp_4x1x8192_pass_through_int8_0")
    _, nalu, calu = get("sdp_3x3x33_bs_int8_reg_0")

    info = PROJECT_INFO.get(PROJECT, {"alu_throughput": 1, "ew_throughput": None})
    alu_tp = info["alu_throughput"]
    alu_fill = None
    if c1 is not None:
        print(f"stream fill/drain constant (N=1)          : {c1} cycles"
              + (f"  (N=8 identical: {c8})" if c8 == c1 else f"  (N=8: {c8})"))
    if cbig is not None and c1 is not None:
        per_elem = (cbig - c1)/(nbig - 1)
        print(f"SDP core throughput (4x1x8192, layout-ok) : "
              f"{cbig}/{nbig} = {cbig/nbig:.3f} cyc/elem; "
              f"(cycles - fill)/(N-1) = {per_elem:.4f} -> {alu_tp} elem/cycle"
              f" (spec SDP_BS/BN_THROUGHPUT_{alu_tp})")
    if calu is not None and nalu is not None:
        streaming = (nalu + alu_tp - 1)//alu_tp
        alu_fill = calu - streaming
        print(f"X1-ALU path fill (3x3x33 bs, reg operand) : {calu} - ceil({nalu}/{alu_tp})"
              f" = {alu_fill} cycles  -> relu_pipeline_latency = {alu_fill}")
    print()
    print("recorded as DMA-inclusive references, NOT core constants:")
    for name in ("sdp_1x8192x1_pass_through_int8_0", "sdp_8192x1x1_pass_through_int8_0",
                 "sdp_4x22x42_bypass_int8", "sdp_5x24x18_bs_int8_mem_0",
                 "sdp_23x13x42_bs_int8_mem_0", "sdp_3x3x33_bn_int8_mem_0",
                 "sdp_3x3x33_bn_int8_reg_0"):
        kind, count, cycles = get(name)
        if cycles is not None:
            print(f"    {name:42s} {cycles:9d} cycles / {count} elements"
                  f" = {cycles/count:.3f} cyc/elem")
    print()
    print(f"-> [sfu] calibrated mapping ({PROJECT} SDP): num_units_per_chip = 1,"
          f" lanes = {alu_tp},")
    print("   setup_cycle = 0, relu_initiation_interval = 1, relu_pipeline_latency = "
          + (str(alu_fill) if alu_fill is not None else "n/a"))

    # LUT/EW class (X2/Y path): only measurable on an EW/LUT-enabled config. The EW
    # engine's own throughput is the bottleneck of these traces (spec EW_THROUGHPUT);
    # the fill comes from the small LE-exp trace after removing the streaming term at
    # that throughput, and LUT-content independence is cross-checked across exp vs
    # linear LE vs LO tables. The large trace is a DMA-inclusive upper-bound reference
    # (its W=13 lines and CVSRAM operand pay beat inefficiency), never the core II.
    lut_rows = {n: v for n, v in rows.items() if v[0] == "lut"}
    if lut_rows and info["ew_throughput"]:
        ew_tp = info["ew_throughput"]
        print()
        _, n_exp, c_exp = get("sdp_3x3x33_ew_le_exp_int8")
        _, n_big, c_big = get("sdp_13x4x29_ew_lut_int8")
        if c_exp is not None:
            ew_fill = c_exp - (n_exp + ew_tp - 1)//ew_tp
            print(f"EW/LUT path fill (3x3x33 le_exp)          : {c_exp} - "
                  f"ceil({n_exp}/{ew_tp}) = {ew_fill} cycles")
            print(f"EW/LUT path II                            : {ew_tp} elem/cycle"
                  f" (spec SDP_EW_THROUGHPUT_{ew_tp}; per-element after fill:"
                  f" {(c_exp - ew_fill)/n_exp:.4f})")
        if c_big is not None and c_exp is not None and n_big != n_exp:
            ii = (c_big - c_exp)/(n_big - n_exp)
            print(f"    DMA-inclusive reference (13x4x29 spread): "
                  f"({c_big} - {c_exp})/({n_big} - {n_exp}) = {ii:.4f} cyc/elem")
        for name in ("sdp_3x3x33_ew_le_lin_int8", "sdp_3x3x32_ew_lo_lin_int8",
                     "sdp_3x3x33_ew_int8_reg_0"):
            kind, count, cycles = get(name)
            if cycles is not None:
                print(f"    cross-check {name:36s}: {cycles} cycles / {count} elements")
        print("    -> LUT-class (sigmoid/tanh/exp/recip) profile: lanes = "
              f"{ew_tp}, initiation_interval = 1, pipeline_latency = EW/LUT fill above"
              " (LUT table content independent)")
    elif PROJECT == "nv_small":
        print("   LUT-class ops: ARCHITECTURALLY ABSENT on nv_small (spec"
              " SDP_LUT_DISABLE/SDP_EW_DISABLE); calibrated on nv_small_256_full.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
