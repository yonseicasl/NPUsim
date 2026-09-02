#!/usr/bin/env python3
"""Measure the nv_small pipeline's fixed launch/latency/drain constants from RTL logs.

The end-to-end window (last engine D_OP_ENABLE accept -> last output DBB write beat)
decomposes into
    launch      = first DBB read request  - enable accept
    mem_latency = first DBB read data     - first DBB read request
    stream      = last DBB read data      - first DBB read data
    drain       = last output DBB write   - last DBB read data
launch, mem_latency and (for a trace whose compute finishes under the read stream)
drain are PIPELINE DEPTHS -- structural constants of the RTL and its testbench, not
workload properties. mem_latency's expected value is knowable without any golden at
all: the harness's AXIResponder pre-fills its response pipe with AXI_R_LATENCY=32
bubbles advancing one per cycle (verif/verilator/nvdla.cpp), i.e. ~33 cycles.

CALIBRATION DISCIPLINE: constants are averaged ONLY over the traces NPUsim cannot
express (gen_nvdla_workload.py's exclusion list + sdp) -- disjoint from every trace
the cycle comparison evaluates, so the evaluated traces stay a true holdout for
these constants.  All ticks are halved (2 ticks/cycle).
"""
import re, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_nvdla_workload import EXCLUDED   # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
EMU = os.path.join(ROOT, "ext/nvdla/outdir/nv_small/verilator_emu")
OP_ENABLE_WORDS = {0x3010//4, 0x4008//4, 0x5008//4, 0x6008//4, 0x7008//4, 0x9038//4}

def gaps(trace):
    log = os.path.join(EMU, f"run_{trace}", "run.log")
    if not os.path.exists(log): return None
    text = open(log, errors="replace").read()
    if "*** PASS" not in text: return None
    enable = None
    for m in re.finditer(r"\((\d+)\) write to nvdla: addr ([0-9a-f]+), data 0*1\b", text):
        if int(m.group(2), 16) in OP_ENABLE_WORDS:
            enable = int(m.group(1))
    ar = [int(m.group(1)) for m in re.finditer(r"\((\d+)\) DBB: read request", text)]
    rd = [int(m.group(1)) for m in re.finditer(r"\((\d+)\) DBB: read push", text)]
    wr = [int(m.group(1)) for m in re.finditer(r"\((\d+)\) DBB: write, last tick", text)]
    if enable is None or not wr: return None
    t = lambda a, b: (b - a) // 2
    if not ar:      # SDP pass-through style: no conv reads... still has reads (MRDMA)
        return None
    return {
        "launch": t(enable, ar[0]),
        "mem_latency": t(ar[0], rd[0]),
        "stream": t(rd[0], rd[-1]),
        "drain": t(rd[-1], wr[-1]),
        "end_to_end": t(enable, wr[-1]),
    }

def main():
    calib = sorted(set(EXCLUDED) | {"sdp_1x1x1_pass_through_int8"})
    rows = []
    print(f"{'trace':<38}{'set':>6}{'launch':>8}{'memlat':>8}{'stream':>9}{'drain':>8}{'e2e':>9}")
    for d in sorted(os.listdir(EMU)):
        if not d.startswith("run_"): continue
        trace = d[4:]
        g = gaps(trace)
        if g is None: continue
        which = "CAL" if trace in calib else "hold"
        rows.append((trace, which, g))
        print(f"{trace:<38}{which:>6}{g['launch']:>8}{g['mem_latency']:>8}"
              f"{g['stream']:>9}{g['drain']:>8}{g['end_to_end']:>9}")
    cal = [g for _, w, g in rows if w == "CAL"]
    if cal:
        for key in ("launch", "mem_latency", "drain"):
            vals = [g[key] for g in cal]
            print(f"\nCAL {key}: mean={sum(vals)/len(vals):.1f} min={min(vals)} max={max(vals)}")
        setup = sum(g["launch"] + g["mem_latency"] + g["drain"] for g in cal)/len(cal)
        print(f"CAL setup constant (launch + mem_latency + drain): {setup:.1f} cycles")

if __name__ == "__main__":
    main()
