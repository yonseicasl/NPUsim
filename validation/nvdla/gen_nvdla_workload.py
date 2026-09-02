#!/usr/bin/env python3
"""Generate NPUsim workloads (network cfg + mapping + weights) for NVDLA nv_small traces.

Dimensions are parsed from each trace's REGISTER PROGRAMMING (the trace_parser sequence
log), not from the trace's file name -- the name's WxHxC / RxSxCxK convention is not
authoritative (dc_13x15x64_5x3x64x16 actually programs R=3,S=5) and several traces
program don't-care fields with junk that must not be read as configuration (the
CRC-passing dc_1x1x8 point programs stride=1x2, dilation=4x14 -- both meaningless for
its 1x1 kernel).

Selection rules (each exclusion is a statement about NPUsim expressiveness, recorded
here rather than silently skipped):
  - dilation > 1 with R,S > 1: NPUsim's mapping table has no dilation dimension.
    Excludes dc_6x8x192 (dil 4x4), dc_8x8x36_dilation (dil 2x2).
  - asymmetric stride: the mapping's STRIDE column is scalar. Excludes dc_24x33x55 (4x3).
  - kernel taller than the input with implicit padding rows (H=1, R=3, pad_top=15):
    excludes dc_8192x1x1 -- NPUsim's derived input extent cannot describe it.

Padding: nebula supports symmetric padding only (output = (H + 2*pad_h - R)/stride + 1),
while the traces pad asymmetrically. Only the TOTAL pad enters the size arithmetic, so
pad_total = (Ho-1)*stride + R - H is split as pad_h = total//2; an ODD total leaves one
row/column that is expressed by enlarging the input instead (H+1) -- the smallest
encoding that reproduces the exact programmed output dims. The enlargement (and only
it) overstates the input volume; compare_full.py reports NPUsim against this DECLARED
encoding (must be exact -- that validates the machinery) and the encoding against the
true tensor volume (the documented distortion) as separate columns.

Layer-0 dummy: nebula's OpenCV loader accepts only 1- or 3-channel network inputs
(ext/nebula/models/networks/convolutional.cc), so each network carries a 1x1x3->C
convolution as layer 0 purely to shape the input; the NVDLA layer is layer 1 and
result/…/layer_1.txt is what the comparison reads. Layer 0's timing is not compared.

Mapping policy v1 (documented limitation): Atomic-C=8 on PE_Y, Atomic-K=8 on PE_X,
all remaining loops (K/8, C/8, R, S, P, Q) at the GLB level -- data physically lives
in CBUF. NPUsim's psum-retention rule then sees a reduction split above the array with
output loops beside it and charges GLB psum round trips; the real nv_small accumulates
those partials inside CACC without touching CBUF. Cycle comparisons made with this
mapping therefore overstate GLB output traffic by construction; the DBB (DRAM-level)
traffic is unaffected because no DRAM-level factor exceeds 1 for CBUF-resident traces.
Modeling CACC as the PE-array temporal buffer with a capacity-gated retention rule is
the known follow-up (same class of change as E20-3), not a config knob.
"""
import re, os, sys, math

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BASE = os.path.join(ROOT, "ext/nvdla/outdir/nv_small/verilator/test")

EXCLUDED = {
    "dc_6x8x192_3x3x192x32_int8_0": "effective dilation 4x4 (no dilation dim in NPUsim mapping)",
    "dc_8x8x36_4x4x36x16_dilation_int8_0": "effective dilation 2x2",
    "dc_24x33x55_5x5x55x25_int8_0": "asymmetric stride 4x3 (mapping STRIDE is scalar)",
    "dc_8192x1x1_2x3x1x41_int8_0": "R=3 > H=1 with pad_top=15 (input extent not expressible)",
}

def parse_regs(trace):
    log = os.path.join(BASE, trace, "trace_parser_cmd_sequence_command.log")
    regs = {}
    for line in open(log):
        m = re.match(r"WRITE (NVDLA_\w+) (\S+) NA ([0-9A-Fa-f]+) NA", line)
        if m:
            regs[(m.group(1), m.group(2))] = int(m.group(3), 16)
    return regs

def dims(trace):
    g = parse_regs(trace)
    din0 = g[("NVDLA_CSC","D_DATAIN_SIZE_EXT_0_0")]; din1 = g[("NVDLA_CSC","D_DATAIN_SIZE_EXT_1_0")]
    w0 = g[("NVDLA_CSC","D_WEIGHT_SIZE_EXT_0_0")];   w1 = g[("NVDLA_CSC","D_WEIGHT_SIZE_EXT_1_0")]
    do0 = g[("NVDLA_CSC","D_DATAOUT_SIZE_0_0")];     do1 = g[("NVDLA_CSC","D_DATAOUT_SIZE_1_0")]
    st = g[("NVDLA_CSC","D_CONV_STRIDE_EXT_0")]
    return {
        "W": (din0 & 0x1fff) + 1, "H": ((din0 >> 16) & 0x1fff) + 1, "C": (din1 & 0x1fff) + 1,
        "S": (w0 & 0x1f) + 1, "R": ((w0 >> 16) & 0x1f) + 1,
        "K": ((w1 >> 16) & 0x1fff) + 1,
        "Wo": (do0 & 0x1fff) + 1, "Ho": ((do0 >> 16) & 0x1fff) + 1,
        "sx": (st & 7) + 1, "sy": ((st >> 16) & 7) + 1,
    }

def short(trace):
    return "nvdla_" + trace.replace("_int8_0", "").replace("_int8", "")

def emit(trace):
    d = dims(trace)
    if trace in EXCLUDED:
        print(f"SKIP {trace}: {EXCLUDED[trace]}"); return None
    if d["sx"] != d["sy"] and d["H"] > 1 and d["W"] > 1:
        print(f"SKIP {trace}: asymmetric stride {d['sx']}x{d['sy']}"); return None
    stride = d["sx"] if d["W"] > 1 else d["sy"]
    name = short(trace)
    W,H,C,R,S,K,Wo,Ho = d["W"],d["H"],d["C"],d["R"],d["S"],d["K"],d["Wo"],d["Ho"]
    # Symmetric-pad encoding of the (possibly asymmetric) programmed padding: only the
    # total enters nebula's size arithmetic; an odd total becomes one extra input row/col.
    tot_h = (Ho - 1)*stride + R - H
    tot_w = (Wo - 1)*stride + S - W
    if tot_h < 0 or tot_w < 0:
        print(f"SKIP {trace}: negative pad total ({tot_h},{tot_w})"); return None
    H_enc, pad_h = H + (tot_h % 2), tot_h // 2
    W_enc, pad_w = W + (tot_w % 2), tot_w // 2

    # ---- network: dummy 1x1 3->C conv (layer 0) + the NVDLA conv (layer 1) ----
    os.makedirs(f"{ROOT}/configs/networks", exist_ok=True)
    with open(f"{ROOT}/configs/networks/{name}.cfg", "w") as f:
        f.write(f"""# NVDLA nv_small official trace {trace}: input {W}x{H}x{C},
# kernel (R x S) {R}x{S}, K={K}, stride {stride} -> output {Wo}x{Ho}x{K}.
# Encoded as {W_enc}x{H_enc} input with symmetric pad ({pad_h},{pad_w}) -- see the
# padding note in gen_nvdla_workload.py; input-volume delta is reported, not hidden.
# Dims decoded from the trace's CSC register programming (see gen_nvdla_workload.py).
# Layer 0 is a shaping dummy (nebula's loader takes only 1/3-channel inputs);
# the comparison reads layer_1.txt only.
[net]
height={H_enc}
width={W_enc}
channels=3
batch=1

[data]
dataset=imagenet
test=datasets/imagenet/test.lst
labels=datasets/imagenet/labels.lst
weight=weights/{name}.wgh
top=1

[convolutional]
filters={C}
size=1
stride=1
padding=0
activation=linear

[convolutional]
filters={K}
filter_height={R}
filter_width={S}
stride={stride}
padding_height={pad_h}
padding_width={pad_w}
activation=linear

[cost]
type=l2
""")
    # ---- mapping ----
    os.makedirs(f"{ROOT}/configs/mappings/nvdla_small/{name}", exist_ok=True)
    k_pe = min(8, K); c_pe = min(8, C)
    k_glb = math.ceil(K / k_pe); c_glb = math.ceil(C / c_pe)
    def row(vals): return ", ".join(str(v) for v in vals) + ","
    with open(f"{ROOT}/configs/mappings/nvdla_small/{name}/ws.map", "w") as f:
        f.write("# The order of mapping table is\n# K, B, P, Q, C, R, S, H, W, GROUP, STRIDE\n\n")
        f.write(f"# layer 0: shaping dummy 1x1x3->{C} conv on {W_enc}x{H_enc}; everything at GLB.\n")
        f.write("[convolutional]\n")
        f.write(f"MAC     = {row([1,1,1,1,1,1,1,0,0,1,1])}\n")
        f.write(f"PE      = {row([1,1,1,1,1,1,1,0,0,1,1])}\n")
        f.write(f"PE_X    = {row([min(8,C),1,1,1,1,1,1,0,0,1,1])}\n")
        f.write(f"PE_Y    = {row([1,1,1,1,3,1,1,0,0,1,1])}\n")
        f.write(f"GLB     = {row([math.ceil(C/min(8,C)),1,H_enc,W_enc,1,1,1,0,0,1,1])}\n")
        f.write(f"CHIPS_X = {row([1,1,1,1,1,1,1,0,0,1,1])}\n")
        f.write(f"CHIPS_Y = {row([1,1,1,1,1,1,1,0,0,1,1])}\n")
        f.write(f"DRAM    = {row([1,1,1,1,1,1,1,0,0,1,1])}\n\n")
        f.write(f"# layer 1: {trace}. Atomic-K={k_pe} on PE_X, Atomic-C={c_pe} on PE_Y,\n")
        f.write(f"# remaining K/{k_pe}={k_glb}, C/{c_pe}={c_glb}, R={R}, S={S}, P={Ho}, Q={Wo} at GLB (CBUF-resident).\n")
        f.write("[convolutional]\n")
        f.write(f"MAC     = {row([1,1,1,1,1,1,1,0,0,1,stride])}\n")
        f.write(f"PE      = {row([1,1,1,1,1,1,1,0,0,1,stride])}\n")
        f.write(f"PE_X    = {row([k_pe,1,1,1,1,1,1,0,0,1,stride])}\n")
        f.write(f"PE_Y    = {row([1,1,1,1,c_pe,1,1,0,0,1,stride])}\n")
        f.write(f"GLB     = {row([k_glb,1,Ho,Wo,c_glb,R,S,0,0,1,stride])}\n")
        f.write(f"CHIPS_X = {row([1,1,1,1,1,1,1,0,0,1,stride])}\n")
        f.write(f"CHIPS_Y = {row([1,1,1,1,1,1,1,0,0,1,stride])}\n")
        f.write(f"DRAM    = {row([1,1,1,1,1,1,1,0,0,1,stride])}\n")
    # ---- weights: zero-filled floats, layer0 (C + 3*C) + layer1 (K + C*R*S*K) ----
    n_floats = (C + 3*C) + (K + C*R*S*K)
    with open(f"{ROOT}/models/weights/{name}.wgh", "wb") as f:
        f.write(b"\x00" * (4 * n_floats))
    print(f"{name}: {W}x{H}x{C} k{R}x{S} K{K} s{stride} -> {Wo}x{Ho}  (PE {k_pe}x{c_pe}, GLB K{k_glb} C{c_glb} P{Ho} Q{Wo})")
    return name

if __name__ == "__main__":
    traces = sys.argv[1:] or sorted(
        t for t in os.listdir(BASE) if t.startswith("dc_") and t != "dc_1x1x8_1x1x8x1_int8_0")
    for t in traces:
        if os.path.exists(os.path.join(BASE, t, "trace_parser_cmd_sequence_command.log")):
            emit(t)
