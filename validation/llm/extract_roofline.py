#!/usr/bin/env python3
"""Extract roofline data from NPUsim LLM projection-GEMM result files.

Case Study 1 (evaluation discussion 2026-09-03 Sec 2) compares THREE lines per workload:

  1. machine roofline      -- min(P_peak, BW * AI_machine), with AI_machine from the
                              ALGORITHMIC minimum DRAM traffic (each GEMM tensor moved
                              once: M*K + K*N + M*N bytes at int8);
  2. mapping-aware bound   -- min(P_peak, BW * AI_mapped), with AI_mapped from the DRAM
                              traffic the MAPPING actually generated (measured payload
                              transactions x line bytes) -- tiling/reuse included, temporal
                              stalls not;
  3. effective performance -- N_MAC / (cycles / freq) from the simulated schedule.

The machine-vs-mapping gap is traffic the mapping added (refetch/padding); the
mapping-vs-effective gap is the architecture-dependent execution overhead (utilization,
fill/drain, buffering, imperfect overlap) that only the timeline sees.

For each gemmini_llm_*_layer_0.txt, pull:
  N_MAC              = "# of computations"
  cycles             = "Compute-schedule latency" (the validated metric)
  DRAM bytes         = (input+weight+output payload transactions) x DRAM line bytes
GEMM shapes (M, K, N) come from gen_llm_configs.py (same directory), so the algorithmic
traffic is derived from the same source that generated the workloads.

Gemmini roofline constants (gemmini.cfg): 16x16 INT8 array @200MHz, DRAM 8192 MB/s,
DRAM line 8 bytes. Emits a CSV to stdout.

Run from models/ (where the result files land):
  python3 ../validation/llm/extract_roofline.py > roofline.csv
"""
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_llm_configs import MODELS, M_SETS   # noqa: E402  (single source for GEMM shapes)

# Gemmini config constants (configs/accelerators/gemmini.cfg)
ARRAY_MACS = 16 * 16
FREQ_HZ = 200e6
PEAK_MAC_PER_S = ARRAY_MACS * FREQ_HZ          # 5.12e10 MAC/s
BW_BYTES_PER_S = 8192 * 1e6                     # 8.192e9 B/s (bandwidth = 8192 MB/s)
DRAM_LINE_BYTES = 8                             # dram line_size = 8, int8

RESULT_GLOB = "gemmini_llm_*_layer_0.txt"


def parse(path):
    text = open(path, errors="replace").read()
    def find(pat):
        m = re.search(pat, text)
        return m.group(1) if m else None
    n_mac = find(r"# of computations\s*:\s*(\d+)")
    cycles = find(r"Compute-schedule latency\s*:\s*([\d.]+)")
    crit = find(r"Critical-path latency\s*:\s*([\d.]+)")
    # DRAM payload transactions: "* <name> : payload/metadata/serialized"
    block = re.search(r"Multi-chip <-> DRAM transactions.*?"
                      r"Input data\s*:(\d+)/\d+/\d+.*?"
                      r"Weight\s*:(\d+)/\d+/\d+.*?"
                      r"Output data\s*:(\d+)/\d+/\d+", text, re.S)
    if not (n_mac and cycles and block):
        return None
    in_txn, wt_txn, out_txn = (int(block.group(i)) for i in (1, 2, 3))
    dram_bytes = (in_txn + wt_txn + out_txn) * DRAM_LINE_BYTES
    n_mac = int(n_mac)
    cycles = float(cycles)
    ai = n_mac / dram_bytes if dram_bytes else float("inf")
    p_eff = n_mac / (cycles / FREQ_HZ)
    p_bound = min(PEAK_MAC_PER_S, BW_BYTES_PER_S * ai)
    eff = p_eff / p_bound if p_bound else 0.0
    return {
        "n_mac": n_mac, "cycles": cycles, "crit": crit,
        "dram_bytes": dram_bytes, "ai": ai,
        "p_eff": p_eff, "p_bound": p_bound, "eff": eff,
        "in_txn": in_txn, "wt_txn": wt_txn, "out_txn": out_txn,
    }


def machine_bound(model, label, phase):
    """Machine-roofline line: algorithmic minimum DRAM traffic (each tensor once, int8:
    input M*K + weight K*N + output M*N bytes) -> AI_machine and its roofline bound."""
    shapes = {lbl: (k, n) for lbl, k, n in MODELS.get(model, [])}
    if label not in shapes or phase not in M_SETS:
        return None
    k, n = shapes[label]
    m = M_SETS[phase]
    algo_bytes = m*k + k*n + m*n
    ai_machine = (m*k*n) / algo_bytes
    p_machine = min(PEAK_MAC_PER_S, BW_BYTES_PER_S * ai_machine)
    return {"algo_bytes": algo_bytes, "ai_machine": ai_machine, "p_machine": p_machine}


def name_fields(fname):
    # gemmini_llm_<model>_<label>_<phase>_layer_0.txt
    stem = fname[len("gemmini_llm_"):-len("_layer_0.txt")]
    phase = "decode" if stem.endswith("_decode") else "prefill"
    stem = stem[: -len("_" + phase)]
    parts = stem.split("_", 1)
    model = parts[0]
    label = parts[1] if len(parts) > 1 else ""
    return model, label, phase


def result_paths():
    """Yield (workload_stem, layer_0 path). Supports both layouts: the current
    result/gemmini/llm_*/<mapping>/layer_0.txt tree (repo root or cwd) and the legacy flat
    gemmini_llm_*_layer_0.txt files in the cwd."""
    repo = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    seen = set()
    for base in (os.getcwd(), repo):
        for path in sorted(glob.glob(os.path.join(base, "result", "gemmini",
                                                  "llm_*", "*", "layer_0.txt"))):
            stem = os.path.basename(os.path.dirname(os.path.dirname(path)))
            if stem not in seen:
                seen.add(stem)
                yield "gemmini_" + stem + "_layer_0.txt", path
    for path in sorted(glob.glob(os.path.join(os.getcwd(), RESULT_GLOB))):
        stem = os.path.basename(path)[len("gemmini_"):-len("_layer_0.txt")]
        if stem not in seen:
            seen.add(stem)
            yield os.path.basename(path), path


def main():
    rows = []
    for fname, path in result_paths():
        d = parse(path)
        if not d:
            continue
        model, label, phase = name_fields(fname)
        rows.append((model, label, phase, d))

    # Three lines per workload (Case Study 1): machine roofline >= mapping-aware bound >=
    # effective. gap_mapping = mapping-added traffic; gap_exec = execution overhead the
    # timeline exposes below the mapping-aware bound.
    print("model,label,phase,N_MAC,cycles,DRAM_bytes,algo_bytes,"
          "AI_machine,AI_mapped,P_machine_GMAC_s,P_mapbound_GMAC_s,P_eff_GMAC_s,"
          "gap_mapping,gap_exec,efficiency")
    for model, label, phase, d in rows:
        mb = machine_bound(model, label, phase)
        algo_bytes = mb["algo_bytes"] if mb else ""
        ai_machine = f"{mb['ai_machine']:.4f}" if mb else ""
        p_machine = mb["p_machine"] if mb else None
        gap_mapping = f"{d['p_bound']/p_machine:.4f}" if p_machine else ""
        p_machine_s = f"{p_machine/1e9:.3f}" if p_machine else ""
        print(f"{model},{label},{phase},{d['n_mac']},{d['cycles']:.0f},"
              f"{d['dram_bytes']},{algo_bytes},"
              f"{ai_machine},{d['ai']:.4f},"
              f"{p_machine_s},{d['p_bound']/1e9:.3f},{d['p_eff']/1e9:.3f},"
              f"{gap_mapping},{d['eff']:.4f},{d['eff']:.4f}")
    print(f"# {len(rows)} configs; PEAK={PEAK_MAC_PER_S/1e9:.2f} GMAC/s, "
          f"BW={BW_BYTES_PER_S/1e9:.2f} GB/s", file=sys.stderr)
    print("# columns: P_machine >= P_mapbound >= P_eff; gap_mapping = P_mapbound/P_machine"
          " (mapping-added traffic), gap_exec = P_eff/P_mapbound (execution overhead)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
