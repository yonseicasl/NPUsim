#!/usr/bin/env python3
"""Decode layers and attention of a follow-up run (results.json: decode_stack, attention) vs NPUsim.

NPUsim layer = its four Llama-3-8B projection GEMMs at batch B (summed from compare_gemm[_<yaml>].csv, which
compare_v5e.py writes for the same results dir: the streamed decode shapes are in gemm.csv) + the KV-cache read
at the specified 819 GB/s, as validation/palm models a layer (NPUsim keeps KV streams analytical even on MOSAIC
DRAM). gemm_only is compared with the GEMMs alone, full with GEMMs + KV.
Optionally --kv-times TAG=CSV replaces that tag's KV roofline with DRAM supply
predictions from kv_streams.py. This does not model full attention computation.
--output writes predictions, metrics and input hashes without modifying measurements.

Usage: python3 validation/tpu_v5e/compare_stack.py RESULTS_DIR [--dram-tag hbm2e_v5e ...]
"""
import argparse
import csv
import json
import math
from pathlib import Path

from kv_streams import digest, kv_bytes
import os
import statistics

PROJ = [(4096, 6144), (4096, 4096), (4096, 28672), (14336, 4096)]
HKV, DH, HBM = 8, 128, 819e9


def load_kv_times(path):
    """Fail closed on wrong geometry, duplicate cases or invalid model predictions."""
    values = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            case = int(row["B"]), int(row["L"])
            if case in values:
                raise ValueError(f"duplicate KV case {case} in {path}")
            us = float(row["kv_us"])
            if int(row["kv_bytes"]) != kv_bytes(*case):
                raise ValueError(f"wrong KV geometry for {case} in {path}")
            if not math.isfinite(us) or us <= 0:
                raise ValueError(f"invalid KV time for {case} in {path}")
            if row["method"] not in ("exact", "extrapolated"):
                raise ValueError(f"unknown KV prediction method in {path}")
            values[case] = (us, row["method"])
    return values


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--dram-tag", nargs="*", default=["hbm2e_v5e"], help="compare_gemm_<tag>.csv besides analytic")
    ap.add_argument("--kv-times", action="append", default=[], metavar="TAG=CSV",
                    help="use modeled KV stream times for one DRAM tag; generate with kv_streams.py")
    ap.add_argument("--output", type=Path, help="optional JSON predictions, metrics and input hashes")
    a = ap.parse_args()
    d = json.load(open(os.path.join(a.results, "results.json")))
    models = {"analytic": "compare_gemm.csv", **{t: f"compare_gemm_{t}.csv" for t in a.dram_tag}}
    kv_models, kv_paths = {}, {}
    for spec in a.kv_times:
        tag, sep, path = spec.partition("=")
        if not sep or tag not in models or tag in kv_models:
            ap.error("--kv-times requires a unique known TAG=CSV")
        kv_models[tag], kv_paths[tag] = load_kv_times(path), path
    # Validate coverage before reporting any metrics.
    for tag, values in kv_models.items():
        needed = {(r["B"], r["L"]) for r in d["decode_stack"] if r["variant"] == "full"}
        if needed - values.keys():
            ap.error(f"missing KV cases for {tag}: {sorted(needed - values.keys())}")
    gemm = {}
    for name, f in models.items():
        sim = {(int(r["M"]), int(r["K"]), int(r["N"])): float(r["npusim_us"])
               for r in csv.DictReader(open(os.path.join(a.results, f)))}
        gemm[name] = lambda B, sim=sim: sum(sim[(B, K, N)] for K, N in PROJ)
    stack = {(r["B"], r["L"], r["variant"]): r["us_per_layer"] for r in d["decode_stack"]}
    errs, predictions = {}, []
    if kv_models:
        print("KV model: contiguous DRAM supply only; no QK/softmax/AV or vector costs.")
        print("Per-row method records exact vs extrapolated stream times; this is not full attention validation.")
    print(f"{'B':>4} {'L':>5} {'variant':9} {'measured':>9} {'KVroof':>7} | "
          + " | ".join(f"{n:>10} {'err':>6}" for n in models) + "   (us per layer)")
    for (B, L, var), meas in sorted(stack.items()):
        kv = B * L * HKV * DH * 2 * 2 / HBM * 1e6 if var == "full" else 0.0
        cells = []
        for n in models:
            modeled_kv, method = kv_models[n][B, L] if n in kv_models and var == "full" else (kv, "analytic" if var == "full" else "none")
            pred = gemm[n](B) + modeled_kv
            e = pred / meas - 1
            errs.setdefault((n, var), []).append(e)
            predictions.append(dict(B=B, L=L, variant=var, model=n, measured_us=meas,
                                    gemm_us=gemm[n](B), kv_us=modeled_kv, kv_method=method,
                                    predicted_us=pred, err_pct=100 * e))
            cells.append(f"{pred:10.1f} {100 * e:+5.1f}%")
        print(f"{B:4d} {L:5d} {var:9} {meas:9.1f} {kv:7.1f} | " + " | ".join(cells))
    for (n, var), v in errs.items():
        print(f"{var:9} {n:10}: mean {100 * statistics.mean(v):+.1f}%  MAPE {100 * statistics.mean(map(abs, v)):.1f}%"
              f"  (n={len(v)})")
    if a.output:
        paths = [os.path.join(a.results, "results.json")]
        paths += [os.path.join(a.results, f) for f in models.values()] + list(kv_paths.values())
        metrics = [dict(model=n, variant=var, n=len(v), mean_pct=100 * statistics.mean(v),
                        mape_pct=100 * statistics.mean(map(abs, v))) for (n, var), v in errs.items()]
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(json.dumps(dict(predictions=predictions, metrics=metrics,
            inputs={str(Path(p).resolve()): digest(p) for p in paths},
            kv_models=kv_paths, limitation="GEMM sum + KV supply; not a full attention/decoder simulation"), indent=2) + "\n")
    print("\nattention alone: measured / KV read at 819 GB/s")
    for r in d["attention"]:
        print(f"  B {r['B']:4d} L {r['L']:6d} {r.get('mode', ''):9} {r['us']:9.1f} us  {r['kv_gbps']:6.0f} GB/s"
              f"  x{r['vs_roofline_819']:.2f}")


if __name__ == "__main__":
    main()
