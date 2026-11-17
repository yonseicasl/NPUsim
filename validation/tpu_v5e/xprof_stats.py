#!/usr/bin/env python3
"""Per-op TPU profile (FLOPs, bytes, utilization, bound-by) from saved jax.profiler traces, via XProf's converters.

Every trace directory under RESULTS/traces/<case>/ (as tpu_v5e1_stream_bench.py writes them) is converted with
XProf's hlo_stats and roofline_model tools; their rows go to RESULTS/xprof/{hlo_stats,roofline}.csv with a `case`
column, and op_profile's per-program totals to RESULTS/xprof/programs.csv. Prints the top ops of each case.

Usage: python3 xprof_stats.py RESULTS [--case SUBSTRING] [--top 8] [--list-tools]
Needs xprof (`pip install xprof`; on the server PYTHONPATH=<pip --target dir>).
"""
import argparse
import csv
import glob
import json
import os

from xprof.convert import raw_to_tool_data as r2t

SHOW = ("category", "hlo_op_name", "occurrences", "total_self_time", "measured_flop_rate", "model_flop_rate",
        "measured_memory_bw", "hbm_bw", "operational_intensity", "bound_by")


def table(data):
    """XProf's {cols, rows} data table -> list of dicts (column id -> value)."""
    d = json.loads(data) if isinstance(data, (str, bytes)) else data
    if isinstance(d, list):          # roofline_model returns [table, device info, ...]
        d = d[0]
    cols = [c["id"] for c in d["cols"]]
    return [{c: (v.get("v") if isinstance(v, dict) else v) for c, v in zip(cols, r["c"])} for r in d.get("rows", [])]


def convert(xplanes, tool):
    data, _ = r2t.xspace_to_tool_data(xplanes, tool, {})
    return data.decode() if isinstance(data, bytes) else data


def programs(op_profile):
    """Per-program totals of op_profile: FLOPs utilization, bandwidth utilizations, raw time / FLOPs / bytes."""
    root = json.loads(op_profile)["byProgram"]
    out = []
    for p in root.get("children", []):
        m = p.get("metrics", {})
        bw = m.get("bandwidthUtils", [])
        by = m.get("rawBytesAccessedArray", [])
        out.append({"program": p.get("name"), "time_us": m.get("rawTime", 0) / 1e6, "flops": m.get("rawFlops"),
                    "flops_util": m.get("flops"), **{f"bw_util_{i}": v for i, v in enumerate(bw)},
                    **{f"bytes_{i}": v for i, v in enumerate(by)}})
    return out


def write(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)


def fmt(v):
    if isinstance(v, float):
        return f"{v:.4g}"
    return v[:48] if isinstance(v, str) else v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--case", default="")
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--list-tools", action="store_true")
    ap.add_argument("--counters", action="store_true", help="only report the hardware-counter tables' row counts")
    a = ap.parse_args()
    out = os.path.join(a.results, "xprof")
    os.makedirs(out, exist_ok=True)
    hlo, roof, progs = [], [], []
    for case_dir in sorted(glob.glob(os.path.join(a.results, "traces", "*"))):
        case = os.path.basename(case_dir)
        if a.case not in case:
            continue
        xplanes = sorted(glob.glob(os.path.join(case_dir, "**", "*.xplane.pb"), recursive=True))
        if not xplanes:
            continue
        if a.list_tools:
            print(case, r2t.xspace_to_tool_names(xplanes))
        if a.counters:   # hardware counters: are XProf's perf_counters / utilization_viewer tables filled?
            for tool in ("perf_counters", "utilization_viewer"):
                try:
                    n = len(table(convert(xplanes, tool)))
                except Exception as e:
                    n = f"error {e!r}"[:80]
                print(f"{case}: {tool} rows = {n}")
            continue
        rows = table(convert(xplanes, "hlo_stats"))
        hlo += [dict(case=case, **r) for r in rows]
        roof += [dict(case=case, **r) for r in table(convert(xplanes, "roofline_model"))]
        progs += [dict(case=case, **p) for p in programs(convert(xplanes, "op_profile"))]
        print(f"\n== {case}: {len(rows)} ops")
        for r in sorted(rows, key=lambda r: -(r.get("total_self_time") or 0))[:a.top]:
            print("  " + " | ".join(f"{k}={fmt(r.get(k))}" for k in SHOW if k in r))
    if a.counters:   # counter check only: keep the per-op tables already written
        return
    write(os.path.join(out, "hlo_stats.csv"), hlo)
    write(os.path.join(out, "roofline.csv"), roof)
    write(os.path.join(out, "programs.csv"), progs)
    print(f"\nwrote {out}/hlo_stats.csv ({len(hlo)} rows), roofline.csv ({len(roof)}), programs.csv ({len(progs)})")


if __name__ == "__main__":
    main()
