#!/usr/bin/env python3
"""Eyeriss external-traffic comparison methodology, single-sourced.

Both the phase-3 report (compare.py) and the acceptance gate (validation/check_timing.py)
compare NPUsim's traffic counters against the Eyeriss chip's JSSC'17 Table V figures. Doing it
twice invited exactly the drift it caused: the gate grew a RAW point comparison that is not
apples-to-apples, while the report already had the better-founded treatment. The methodology
therefore lives here and both import it.

TWO AXES, TWO DIFFERENT PROBLEMS -- and neither is a plain point comparison:

DRAM. Table V's DRAM figures are measured AFTER the chip's run-length compression of
ifmaps/ofmaps. NPUsim moves dense data, so comparing its bytes to the chip's compressed bytes
overstates the error by the compression ratio. The comparable reference adds the compression
savings back using the chip's own reported ifmap-zero fractions (RLC ~ 21.33 bits/nonzero; a
layer's ofmap zeros are the next layer's ifmap zeros).

GLB. Table V's GLB figure is SRAM access traffic, and it depends on per-layer mapper
parameters the paper does not publish -- so a point target does not exist. What can be stated
is a BAND:
    lower = distinct dense tensor volumes            (perfect on-chip retention)
    upper = NPUsim's literal per-repetition streaming (retains nothing, re-streams every pass)
Only SRAM-MEDIATED streams belong on this axis: a GLB-bypassed datatype crosses the on-chip
fabric (and NPUsim charges it since the P1-B bypass fix) but never touches the SRAM.
A chip figure INSIDE the band means the model brackets reality.

CORRECTED 2026-08-25. This paragraph used to say that a chip figure ABOVE the upper bound proved
a traffic SOURCE missing from the model that no mapping could reach. Two measurements refuted it,
in order:
  * Moving conv3's output factors up the hierarchy -- same work, same loop bounds -- takes its GLB
    traffic 9.2 -> 32.2 -> 101.2 MB, and the chip's 50.2 MB sits inside that range. The psum GLB
    path exists and fires; how much of it fires is a tiling decision (E20-3).
  * The remaining gap was not a source at all but a HALF-COUNTED one. The psum load was suppressed
    by a latch that never reset within a layer, so the array wrote its psum tile out once per
    revisit and read it back once per output offset (conv3: 19656 write-backs, 312 loads). Pairing
    the two took the upper bound to 73.1 MB and put every layer inside the band (E20-3b).
So read a chip figure above the upper bound as a question about THIS model's accounting -- start
with whether each modeled movement is charged in both directions -- not as proof of an absent
mechanism.
"""
import re
from pathlib import Path

MHZ = 200.0
# JSSC'17 Table V: processing latency (ms), total latency (ms), active PEs,
# GLB accesses (MB), DRAM accesses (MB), ifmap zeros (%)
CHIP = {
    "conv1": (16.5, 20.9, 154, 18.5, 5.0, 0.01),
    "conv2": (39.2, 41.9, 135, 77.6, 4.0, 38.7),
    "conv3": (21.8, 23.6, 156, 50.2, 3.0, 72.5),
    "conv4": (16.0, 18.4, 156, 37.4, 2.1, 79.3),
    "conv5": (10.0, 10.5, 156, 24.9, 1.3, 77.6),
}
LAYER_FILES = {"conv1": 0, "conv2": 2, "conv3": 4, "conv4": 5, "conv5": 6}
RLC_BITS_PER_NONZERO = 21.33
# AlexNet conv shapes, batch 4: (K, P, Q, C, R, S, stride, groups)
SHAPES = {
    "conv1": (96, 55, 55, 3, 11, 11, 4, 1),
    "conv2": (256, 27, 27, 96, 5, 5, 1, 2),
    "conv3": (384, 13, 13, 256, 3, 3, 1, 1),
    "conv4": (384, 13, 13, 384, 3, 3, 1, 2),
    "conv5": (256, 13, 13, 384, 3, 3, 1, 2),
}
LABEL = {"Input data": "input", "Weight": "weight", "Output data": "output"}


def volumes() -> dict:
    """Distinct (fetch-once, dense, halo-free) tensor volumes in MB at 16 bit, batch 4.
    Their sum is the rigorous LOWER bound of GLB<->array traffic (perfect retention) and the
    dense reference the RLC adjustment adds back onto."""
    out = {}
    for name, (K, P, Q, C, R, S, stride, groups) in SHAPES.items():
        ih, iw = (P - 1)*stride + R, (Q - 1)*stride + S
        out[name] = (4*C*ih*iw*2/1e6, K*C*R*S//groups*2/1e6, 4*K*P*Q*2/1e6)
    return out


def rlc_ratio(zero_pct: float) -> float:
    return min(1.0, (1.0 - zero_pct/100.0)*RLC_BITS_PER_NONZERO/16.0)


def dense_equivalent_dram(name: str, volume: dict) -> float:
    """Chip DRAM traffic converted to a dense-equivalent, so it is comparable to NPUsim's
    dense bytes: add back what run-length compression saved on the ifmap and the ofmap."""
    names = list(CHIP)
    index = names.index(name)
    zeros_in = CHIP[name][5]
    # A layer's ofmap zeros are the next layer's ifmap zeros; the last layer reuses its own.
    zeros_out = CHIP[names[index + 1]][5] if index + 1 < len(names) else zeros_in
    input_mb, _, output_mb = volume[name]
    return CHIP[name][4] + input_mb*(1 - rlc_ratio(zeros_in)) + output_mb*(1 - rlc_ratio(zeros_out))


def parse_traffic(layer_path: Path) -> dict:
    """NPUsim traffic counters for one layer, in MB at the configured link widths
    (eyeriss.cfg: GLB<->array link 256 b -> 32 B/txn, multi-chip<->DRAM 64 b -> 8 B/txn).
    `glb_sram_mb` excludes GLB-bypassed datatypes; `glb_fabric_mb` includes them."""
    text = layer_path.read_text(encoding="utf-8")
    glb_section = text.split("Global buffer result", 1)[1]
    table = glb_section.split("serialized", 1)[1].split("# of request", 1)[0]
    bypassed = set(re.search(r"GLB-bypassed \(direct stream\)\s*:(.*)",
                             glb_section).group(1).split())
    rows = re.findall(r"(Input data|Weight|Output data)\s*:\s*\d+/\d+/(\d+)", table)
    dram_section = text.split("Multi-chip <-> DRAM transactions", 1)[1]
    dram_tx = sum(int(v) for v in re.findall(r"/(\d+)\n", dram_section)[:3])
    return {
        "glb_sram_mb": sum(int(v) for n, v in rows if LABEL[n] not in bypassed)*32/1e6,
        "glb_fabric_mb": sum(int(v) for _, v in rows)*32/1e6,
        "dram_mb": dram_tx*8/1e6,
        "bypassed": bypassed,
    }


def evaluate(result_dir: Path) -> dict:
    """Per-layer traffic verdicts on both axes for one Eyeriss result directory."""
    volume = volumes()
    verdicts = {}
    for name, layer_index in LAYER_FILES.items():
        measured = parse_traffic(result_dir / f"layer_{layer_index}.txt")
        chip_dense = dense_equivalent_dram(name, volume)
        lower = sum(volume[name])
        chip_glb = CHIP[name][3]
        upper = measured["glb_sram_mb"]
        verdicts[name] = {
            "dram_mb": measured["dram_mb"],
            "dram_dense_reference": chip_dense,
            "dram_error_pct": (measured["dram_mb"]/chip_dense - 1)*100.0,
            "glb_lower": lower,
            "glb_chip": chip_glb,
            "glb_upper": upper,
            "glb_in_band": lower <= chip_glb <= upper,
            "glb_above_upper": chip_glb > upper,
        }
    return verdicts
