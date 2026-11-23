#!/usr/bin/env python3
"""Decode-latency-first LOGIC screening under the fixed target memory -- GENERATE ONLY: hardware points, workload
scenarios, accelerator configs and first-order supply bounds. Nothing is simulated.

Fixed: the target DRAM of target_profile.json (four dies, 2 x 2, one tier; 1024 channels, 24 GB, 32 TB/s nominal
peak in aggregate). It is NOT a search axis. MOSAIC's DRAM model cannot represent it as given (consistency
report: 4 banks per channel, 6 GB per die, timings and the four interface domains are missing), so every point
names the memory SCENARIO it would run on, and none of them is "the target":
  hbm2e_refsb          the validated HBM2e model -- reference only
  ideal_aggregate      ONE endpoint of 1024 channels x 512 DQ x 500 Mbps (stand-in yaml; 32.8 TB/s nominal, capacity
                       wrong by construction): an idealized aggregate-supply screen, no four-domain locality, no
                       contention between dies
Decision variables here (compute baseline tpu_v5e, 1.5 GHz, 128 x 128 arrays; tpu_v7.cfg is separate, unmeasured):
  mxus     number of MXUs (weight ingest scales with it: one 128-element bf16 row = 256 B per cycle per MXU)
  feed_kb  DRAM line = HBM-to-VMEM staging link width per cycle (1 KB = 1.536 TB/s as in tpu_v5e.cfg)
MXUs, SRAM, feed and interface logic share one logic-area budget and NO AREA MODEL exists: every point other than
the first is an exploratory resource trade-off, with no feasibility or optimality claim.
Workload scenarios (no model is selected; shapes of an 8B-class decoder are a labelled scenario): batch-1 matrix-
vector products, one small-batch sensitivity, one batch-1 KV attention step. Batch-1 on NPUsim: the input vector is
padded to whole DRAM lines, and the TPU's own M = 1 schedule was outside the hardware-compared cohort.
The bound columns are ARITHMETIC (bytes over the slowest of three declared rates), written down so a run can be
compared with what the configuration implies; they are not results.

Usage: python3 logic_manifest.py OUT_DIR
"""
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import dse_pilot as dp

CLOCK = dp.ca.CLOCK_HZ
LOGIC = [(4, 1, "tpu_v5e.cfg as is (baseline)"), (4, 32, "feed widened only"), (8, 32, "2 x MXUs, wide feed"),
         (16, 32, "4 x MXUs, wide feed"), (16, 4, "4 x MXUs, moderate feed")]
MEMORY = {"hbm2e_refsb": dict(nominal=819.2e9, note="validated HBM2e model; reference only"),
          "ideal_aggregate": dict(nominal=1024 * 512 * 500e6 / 8, note="stand-in single endpoint, 1024 ch x 512 DQ x 500 Mbps; "
                                  "NOT the target organisation (capacity, banks, timings, four domains)")}
WORK = {"gemv_4096x4096": dict(kind="gemm", M=1, K=4096, N=4096), "gemv_4096x14336": dict(kind="gemm", M=1, K=4096, N=14336),
        "gemm_m8_4096x4096": dict(kind="gemm", M=8, K=4096, N=4096), "attn_b1_l2048": dict(kind="attention", B=1, L=2048)}


def logic_cfg(base, mxus, feed_kb):
    """tpu_v5e.cfg with `mxus` MXUs and a staging link / DRAM line of feed_kb KB per cycle."""
    cfg = base
    cfg, a = re.subn(r"(?m)^num_chips=4$", f"num_chips={mxus}", cfg)
    lo, hi = cfg.index("[multi_chip]"), cfg.index("[dram]")
    mc, dr = cfg[lo:hi], cfg[hi:]
    bits = feed_kb * 8192
    mc, b = re.subn(r"(?m)^height = 4$", f"height = {mxus}", mc)
    mc, c = re.subn(r"(?m)^bandwidth = 1536000$", f"bandwidth = {int(feed_kb * 1024 * CLOCK / 1e6)}", mc)
    mc, d = re.subn(r"(?m)^bitwidth = 8192$", f"bitwidth = {bits}", mc)
    mc, e = re.subn(r"(?m)^line_size = 8192:8192:8192$", f"line_size = {bits}:{bits}:{bits}", mc)
    dr, f = re.subn(r"(?m)^bitwidth= 8192$", f"bitwidth= {bits}", dr)
    dr, g = re.subn(r"(?m)^line_size=8192:8192:8192$", f"line_size={bits}:{bits}:{bits}", dr)
    assert (a, b, c, d, e, f, g) == (1,) * 7, (a, b, c, d, e, f, g)
    return cfg[:lo] + mc + dr


def main():
    out = os.path.abspath(sys.argv[1])
    if os.path.exists(os.path.join(out, "manifest.csv")):
        raise SystemExit(f"{out} already holds a manifest; nothing was changed")
    os.makedirs(os.path.join(out, "accelerators"), exist_ok=True)
    base = open(os.path.join(dp.ca.ROOT, "configs", "accelerators", f"{dp.ca.ACCEL}.cfg")).read()
    rows = []
    for mxus, feed, what in LOGIC:
        name = f"mxu{mxus}_feed{feed}k"
        cfg = logic_cfg(base, mxus, feed)
        assert (mxus, feed) != (4, 1) or cfg == base
        open(os.path.join(out, "accelerators", f"tpu_v5e_{name}.cfg"), "w").write(cfg)
        ingest = mxus * dp.palm.DIM * dp.ACT * CLOCK            # weights shifted into the arrays, bytes / s
        link = feed * 1024 * CLOCK
        for mem, mf in MEMORY.items():
            rate = min(mf["nominal"], link, ingest)
            limiter = min((("memory nominal peak", mf["nominal"]), ("staging link", link), ("MXU weight ingest", ingest)), key=lambda x: x[1])[0]
            for wn, w in WORK.items():
                r = dict(point=f"{wn}__{mem}__{name}", workload=wn, memory_scenario=mem, memory_note=mf["note"], logic_point=name,
                         logic_what=what, mxus=mxus, array="128 x 128", feed_kb_per_cycle=feed, compute_baseline="tpu_v5e",
                         accelerator_cfg=f"accelerators/tpu_v5e_{name}.cfg",
                         area="none: exploratory resource trade-off" if (mxus, feed) != (4, 1) else "baseline",
                         memory_nominal_TB_per_s=round(mf["nominal"] / 1e12, 3), staging_link_TB_per_s=round(link / 1e12, 3),
                         mxu_weight_ingest_TB_per_s=round(ingest / 1e12, 3), first_order_limiter=limiter)
                if w["kind"] == "gemm":
                    wbytes = w["K"] * w["N"] * dp.ACT
                    r.update(weight_bytes=wbytes, whole_strips=w["N"] % (mxus * dp.palm.DIM) == 0,
                             weight_supply_bound_us=round(wbytes / rate * 1e6, 3),
                             compute_bound_us=round(w["M"] * w["K"] * w["N"] / (mxus * dp.palm.DIM ** 2) / CLOCK * 1e6, 3))
                else:
                    kv = w["B"] * w["L"] * 2 * dp.ca.KVH * dp.ca.DH * dp.ca.KV_BYTES
                    r.update(kv_read_bytes=kv, kv_supply_bound_us=round(kv / min(mf["nominal"], link) * 1e6, 3),
                             note="KV read does not pass the MXU weight path in the model; QK / AV / softmax at configured rates")
                r["status"] = "generated; not simulated"
                rows.append(r)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(out, "manifest.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, keys)
        wr.writeheader()
        wr.writerows(rows)
    json.dump(dict(target_profile="../target_profile.json", fixed="the target DRAM is not a search axis; no point runs on the target organisation itself",
                   logic_points=[dict(mxus=m, feed_kb=f, what=w) for m, f, w in LOGIC], memory_scenarios=MEMORY, workloads=WORK,
                   objective="batch-1 decode latency first; prefill and throughput are not optimized",
                   model_gaps=["target DRAM organisation not representable (banks, capacity, timings)", "four memory-interface domains / locality",
                               "no area, power or thermal model", "standalone vector / SFU operations outside the attention step",
                               "batch-1 input padded to DRAM lines", "256 x 256 arrays not generated (tpu_v7.cfg shows the scaling)",
                               "long-stream rule unchecked for the stand-in yaml"],
                   rows=rows), open(os.path.join(out, "manifest.json"), "w"), indent=1)
    print(f"{len(LOGIC)} logic points x {len(MEMORY)} memory scenarios x {len(WORK)} workloads = {len(rows)} points; nothing simulated")
    seen = set()
    for r in rows:
        if r["workload"] == "gemv_4096x4096" and (r["memory_scenario"], r["logic_point"]) not in seen:
            seen.add((r["memory_scenario"], r["logic_point"]))
            print(f"  {r['memory_scenario']:16s} {r['logic_point']:16s} memory {r['memory_nominal_TB_per_s']:6.2f} | link {r['staging_link_TB_per_s']:6.2f} | "
                  f"MXU ingest {r['mxu_weight_ingest_TB_per_s']:5.2f} TB/s -> {r['first_order_limiter']:18s} | batch-1 4096x4096: supply bound "
                  f"{r['weight_supply_bound_us']:7.2f} us, compute bound {r['compute_bound_us']:.3f} us")


if __name__ == "__main__":
    main()
