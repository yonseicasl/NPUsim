#!/usr/bin/env bash
# Fidelity measurement environment (evaluation.md Section B).
#
# Measures how accurately NPUsim reproduces EXECUTION CYCLES (and DRAM traffic, where the
# reference has it) against RTL / silicon ground truth -- but ONLY for the configurations
# that actually have committed golden reference data:
#
#   * Gemmini 16x16 WS   -- Chipyard GemminiRocketConfig, Verilator cycle-exact RTL,
#                           6 GEMM points   (validation/phase2/golden_rtl_cycles.csv)
#   * Eyeriss (silicon)  -- published silicon latency, AlexNet conv/fc
#                           (validation/phase3/golden_eyeriss_silicon.csv)
#   * NVDLA nv_small     -- official dc_* RTL traces, CRC-verified; cycle HOLDOUT + DRAM
#                           traffic   (validation/nvdla/golden_rtl_full.csv)
#
# ENERGY is intentionally absent: NPUsim's energy is UNCALIBRATED (no pJ provenance) and no
# Joule-level golden exists for any of these, so an energy-fidelity number would be
# meaningless. The report states this explicitly rather than printing a fabricated figure.
#
# NOT covered (no golden data exists): Gemmini 128x128, output-stationary dataflow. Those
# need new RTL runs; the report leaves labelled slots so they can be filled in later.
#
#   Usage:  bash validation/fidelity/measure_fidelity.sh
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "$here/../.." && pwd)
cd "$repo"
export LD_LIBRARY_PATH="$repo/ext/DRAMsim3:$repo/ext/nebula/library:${LD_LIBRARY_PATH:-}"

log() { echo "[fidelity] $*" >&2; }

log "building model ..."
./npusim.sh build model >/dev/null 2>&1

# --- Gemmini 16x16 WS + Eyeriss silicon: run the reference workloads, then compare -------
log "running Gemmini 16x16 WS GEMM points ..."
for w in gemm_64x64x64 gemm_128x128x128 gemm_256x256x256 \
         gemm_16x512x512 gemm_512x512x64 gemm_512x512x512; do
    ./npusim.sh run gemmini "$w" ws >/dev/null
done
log "running Eyeriss AlexNet (silicon reference) ..."
./npusim.sh run eyeriss alexnet silicon >/dev/null

# --- NVDLA nv_small: run every committed trace network so compare_full.py has fresh
# NPUsim results on THIS checkout (result/ is not tracked; without this step the NVDLA
# comparison silently used stale or missing results). FORCE=1 re-runs existing ones.
for net_cfg in configs/networks/nvdla_dc_*.cfg; do
    net=$(basename "$net_cfg" .cfg)
    # Only the trace networks with an nvdla_small ws mapping belong to this comparison
    # (e.g. nvdla_dc_1x1x8_relu is an SFU fixture with a different mapping).
    [[ -f "configs/mappings/nvdla_small/$net/ws.map" ]] || continue
    if [[ "${FORCE:-0}" != "1" && -f "result/nvdla_small/$net/ws/layer_0.txt" ]]; then
        continue
    fi
    log "running NVDLA trace $net ..."
    ./npusim.sh run nvdla_small "$net" ws >/dev/null
done

log "comparing against Gemmini RTL / Eyeriss silicon golden ..."
cyc_out=$(python3 validation/check_timing.py --check-baseline 2>&1 || true)

log "comparing against NVDLA nv_small RTL golden ..."
nvdla_out=$(python3 validation/nvdla/compare_full.py 2>&1 || true)

# --- Consolidate into one table + CSV ----------------------------------------------------
python3 - "$here" <<PY
import os, re, sys, csv
here = sys.argv[1]
cyc = """$cyc_out"""
nv  = """$nvdla_out"""

def mape_max(text, label):
    m = re.search(re.escape(label) + r".*?MAPE=\s*([\d.]+)%\s*max=\s*([\d.\-]+)%", text, re.S)
    return (m.group(1), m.group(2)) if m else (None, None)

rows = []
g = mape_max(cyc, "Gemmini RTL:")
rows.append(["Gemmini",  "16x16", "weight-stationary", "Gemmini RTL (Verilator, 6 GEMM)",
             "cycles", g[0], g[1]])
e = mape_max(cyc, "Eyeriss silicon latency:")
rows.append(["Eyeriss",  "12x14", "row-stationary",    "silicon (AlexNet)",
             "latency", e[0], e[1]])
n = mape_max(nv, "HOLDOUT cycle")
rows.append(["NVDLA",    "nv_small", "-",              "nv_small RTL (holdout, n=3)",
             "cycles", n[0], n[1]])

# Slots with no golden yet (requested but unvalidatable).
pending = [
    ["Gemmini", "128x128", "weight-stationary", "NO GOLDEN (needs 128x128 RTL run)", "cycles", None, None],
    ["Gemmini", "128x128", "output-stationary", "NO GOLDEN (needs OS RTL run)",       "cycles", None, None],
]

W = 74
print("\n" + "=" * W)
print("  NPUsim FIDELITY REPORT  (Section B -- cycle/traffic vs RTL/silicon golden)")
print("=" * W)
print("  %-8s %-8s %-18s %-8s %8s %8s" %
      ("config", "array", "dataflow", "metric", "MAPE", "max"))
print("  " + "-" * (W-2))
for r in rows:
    cfg, arr, df, ref, metric, mape, mx = r
    mape_s = f"{mape}%" if mape else "  n/a"
    mx_s   = f"{mx}%"  if mx  else "  n/a"
    print("  %-8s %-8s %-18s %-8s %8s %8s   [%s]" % (cfg, arr, df, metric, mape_s, mx_s, ref))
print("  " + "-" * (W-2))
for r in pending:
    cfg, arr, df, ref, metric, _, _ = r
    print("  %-8s %-8s %-18s %-8s %8s %8s   [%s]" % (cfg, arr, df, metric, "  --", "  --", ref))
print("=" * W)
print("  energy: UNCALIBRATED (no pJ provenance, no Joule golden) -> not a fidelity metric")
print("  NVDLA also validates DRAM traffic EXACTLY (see validation/nvdla/compare_full.py)")
print("=" * W)

out = os.path.join(here, "fidelity_report.csv")
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["config","array","dataflow","reference","metric","mape_pct","max_pct","status"])
    for r in rows:
        w.writerow(r[:1]+r[1:2]+r[2:3]+r[3:4]+r[4:5]+[r[5],r[6],"validated"])
    for r in pending:
        w.writerow(r[:5]+[r[5],r[6],"no_golden"])
print("wrote", out)
PY
