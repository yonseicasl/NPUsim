#!/usr/bin/env python3
"""Phase-3 silicon validation: NPUsim (eyeriss.cfg, 12x14 RS, silicon.map) vs the
measured Eyeriss chip (Chen et al., JSSC Jan. 2017, Table V).

Chip conditions: AlexNet, batch 4, 1.0 V, core 200 MHz / link 60 MHz, 16-bit.
Comparable latency = paper's "Processing Latency" (DRAM traffic fully overlapped)
vs NPUsim's Computation cycle + Fold fill cycle (V2, silicon-calibrated u=0.11).
Secondary axes: DRAM / GLB<->array traffic in MB (16-bit values).
"""
import os, re

# The traffic-comparison methodology (RLC dense-equivalent DRAM, GLB band, SRAM-vs-fabric
# distinction) is shared with validation/check_timing.py so the two cannot drift apart.
from traffic_reference import CHIP, LAYER_FILES, MHZ, rlc_ratio, RLC_BITS_PER_NONZERO

ROOT = os.path.dirname(os.path.abspath(__file__))

# Distinct (fetch-once, dense, unpadded-output/halo-free-input) tensor volumes in
# MB at 16 bit -- the rigorous LOWER bound of GLB<->array traffic (perfect on-chip
# retention) and the dense reference for the RLC adjustment. b=4.
def volumes():
    out = {}
    for name, (K, P, Q, C, R, S, st, g) in {
        "conv1": (96, 55, 55, 3, 11, 11, 4, 1),
        "conv2": (256, 27, 27, 96, 5, 5, 1, 2),
        "conv3": (384, 13, 13, 256, 3, 3, 1, 1),
        "conv4": (384, 13, 13, 384, 3, 3, 1, 2),
        "conv5": (256, 13, 13, 384, 3, 3, 1, 2)}.items():
        ih, iw = (P - 1)*st + R, (Q - 1)*st + S
        out[name] = (4*C*ih*iw*2/1e6, K*C*R*S//g*2/1e6, 4*K*P*Q*2/1e6)
    return out

def parse(layer_index):
    path = os.path.join(ROOT, f"../../result/eyeriss/alexnet/silicon/layer_{layer_index}.txt")
    txt = open(path).read()
    comp = float(re.search(r"Computation cycle\s*:\s*([\d.]+)", txt).group(1))
    ff = re.search(r"Fold fill cycle\s*:\s*([\d.]+)", txt)
    fold = float(ff.group(1)) if ff else 0.0
    # serialized link transactions are link-bitwidth units (eyeriss.cfg):
    # GLB<->array link = 256 b -> 32 B/txn, multi-chip<->DRAM link = 64 b -> 8 B/txn.
    glb_sec = txt.split("Global buffer result")[1]
    # The chip's Table V "GLB (MB)" is SRAM access traffic, so only SRAM-MEDIATED streams
    # belong on this axis. A GLB-bypassed datatype (Eyeriss streams filter weights from the
    # chip ingress straight into the PE spads) traverses the on-chip fabric -- and NPUsim
    # charges it for that since the P1-B bypass fix -- but never touches the SRAM. Counting it
    # here inflated the band's upper bound about 7x and made the band verdict pass vacuously.
    table = glb_sec.split("serialized")[1].split("# of request")[0]
    bypassed = set(re.search(r"GLB-bypassed \(direct stream\)\s*:(.*)", glb_sec).group(1).split())
    label = {"Input data": "input", "Weight": "weight", "Output data": "output"}
    glb_tx = sum(int(serialized) for name, serialized in re.findall(
                 r"(Input data|Weight|Output data)\s*:\s*\d+/\d+/(\d+)", table)
                 if label[name] not in bypassed)
    dram_sec = txt.split("Multi-chip <-> DRAM transactions")[1]
    dram_tx = sum(int(x) for x in re.findall(r"/(\d+)\n", dram_sec)[:3])
    return comp, fold, glb_tx*32/1e6, dram_tx*8/1e6

def main():
    vols = volumes()
    names = list(LAYER_FILES)
    print(f"{'layer':>6s} {'NPU sched(ms)':>13s} {'chip proc(ms)':>13s} {'err%':>6s}  "
          f"{'NPU DRAM':>8s} {'chip-dense':>10s} {'x':>5s}  "
          f"{'GLB lower':>9s} {'chip':>6s} {'upper':>7s} {'band':>5s}")
    errs, dram_errs, in_band_all = [], [], True
    for i, (name, li) in enumerate(LAYER_FILES.items()):
        comp, fold, glb_mb, dram_mb = parse(li)
        ms = (comp + fold)/(MHZ*1000.0)
        chip_ms, _, _, chip_glb, chip_dram, z_in = CHIP[name]
        e = (ms - chip_ms)/chip_ms*100
        errs.append(abs(e))
        # RLC dense-equivalent chip DRAM: add back the ifmap/ofmap compression
        # savings (a layer's ofmap zeros = the next layer's ifmap zeros; the last
        # layer's ofmap uses its own ifmap-zero fraction as the estimate).
        in_v, wt_v, out_v = vols[name]
        z_out = CHIP[names[i+1]][5] if i + 1 < len(names) else z_in
        chip_dense = chip_dram + in_v*(1 - rlc_ratio(z_in)) + out_v*(1 - rlc_ratio(z_out))
        dram_errs.append(abs(dram_mb/chip_dense - 1)*100)
        # GLB band: lower = distinct volumes (perfect on-chip retention), upper = what THIS mapping
        # streams through the GLB SRAM.
        #
        # CORRECTION (2026-08-20). This comment used to claim that a chip figure above the upper
        # bound meant the model was "missing a traffic SOURCE" that no mapping could reach. That is
        # WRONG, and measuring refuted it. conv3's GLB traffic is mapping-controlled: moving output
        # factors from the array levels up to the GLB -- same total work, same loop bounds -- gives
        #     baseline 9.2 MB  ->  K moved up 32.2 MB  ->  K and P moved up 101.2 MB
        # and the chip's 50.2 MB sits inside that range. The psum GLB path exists and fires; how
        # much of it fires is a tiling decision.
        #
        # What IS true is narrower and worth stating precisely: the model cannot reach the chip's
        # GLB traffic AND its measured latency at the same time. Those same probes take conv3's
        # computation from 3.83M cycles (matching the chip's 4.36M within 2.4%) to 15.3M and 199M,
        # because moving spatial factors up de-parallelizes the array. So the upper bound here is
        # not a hard ceiling on the model -- it is the traffic of the one mapping that also matches
        # the latency. Closing the gap needs psum spilling driven by GLB CAPACITY rather than by the
        # tiling hierarchy alone (E20-3), not a new traffic source bolted on.
        glb_lower = in_v + wt_v + out_v
        in_band = glb_lower <= chip_glb <= glb_mb
        in_band_all &= in_band
        print(f"{name:>6s} {ms:13.2f} {chip_ms:13.1f} {e:+6.1f}  "
              f"{dram_mb:8.1f} {chip_dense:10.1f} {dram_mb/chip_dense:5.2f}  "
              f"{glb_lower:9.1f} {chip_glb:6.1f} {glb_mb:7.1f} {'ok' if in_band else 'FAIL':>5s}")
    print(f"\nlatency MAPE = {sum(errs)/len(errs):.1f}%   max |err| = {max(errs):.1f}%  (n={len(errs)})")
    print(f"DRAM vs RLC-dense-equivalent chip: MAPE = {sum(dram_errs)/len(dram_errs):.1f}%   "
          f"max = {max(dram_errs):.1f}%")
    print(f"GLB band verdict: chip within [perfect-retention, literal-streaming] on "
          f"{'ALL' if in_band_all else 'NOT all'} layers")
    if not in_band_all:
        print("  Layers marked FAIL have chip GLB traffic above what THIS mapping streams through "
              "the GLB. Measured, not assumed: moving output factors up the hierarchy takes conv3 "
              "from 9.2 MB to 32.2 MB to 101.2 MB, so the chip's 50.2 MB IS reachable by a mapping "
              "in this model -- but only at 4x to 50x the computation cycles the chip measures. "
              "The gap is that psum spilling follows the tiling hierarchy here rather than GLB "
              "CAPACITY, so latency and GLB traffic cannot both be matched (E20-3).")
    print("  (An exact point match would additionally need the chip's unpublished per-layer "
          "mapper parameters.)")
    print("NPU sched = (Computation + Fold fill) / 200MHz;  traffic in MB (16-bit data)")

if __name__ == "__main__":
    main()
