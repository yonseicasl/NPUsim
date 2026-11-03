#!/usr/bin/env python3
"""Decode layers and attention of a follow-up run (results.json: decode_stack, attention) vs NPUsim.

NPUsim layer = its four Llama-3-8B projection GEMMs at batch B (summed from compare_gemm[_<yaml>].csv, which
compare_v5e.py writes for the same results dir: the streamed decode shapes are in gemm.csv) + the KV-cache read
at the specified 819 GB/s, as validation/palm models a layer (NPUsim keeps KV streams analytical even on MOSAIC
DRAM). gemm_only is compared with the GEMMs alone, full with GEMMs + KV.

Usage: python3 validation/tpu_v5e/compare_stack.py RESULTS_DIR [--dram-tag hbm2e_v5e ...]
"""
import argparse
import csv
import json
import os
import statistics

PROJ = [(4096, 6144), (4096, 4096), (4096, 28672), (14336, 4096)]
HKV, DH, HBM = 8, 128, 819e9


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--dram-tag", nargs="*", default=["hbm2e_v5e"], help="compare_gemm_<tag>.csv besides analytic")
    a = ap.parse_args()
    d = json.load(open(os.path.join(a.results, "results.json")))
    models = {"analytic": "compare_gemm.csv", **{t: f"compare_gemm_{t}.csv" for t in a.dram_tag}}
    gemm = {}
    for name, f in models.items():
        sim = {(int(r["M"]), int(r["K"]), int(r["N"])): float(r["npusim_us"])
               for r in csv.DictReader(open(os.path.join(a.results, f)))}
        gemm[name] = lambda B, sim=sim: sum(sim[(B, K, N)] for K, N in PROJ)
    stack = {(r["B"], r["L"], r["variant"]): r["us_per_layer"] for r in d["decode_stack"]}
    errs = {}
    print(f"{'B':>4} {'L':>5} {'variant':9} {'measured':>9} {'KV':>7} | "
          + " | ".join(f"{n:>10} {'err':>6}" for n in models) + "   (us per layer)")
    for (B, L, var), meas in sorted(stack.items()):
        kv = B * L * HKV * DH * 2 * 2 / HBM * 1e6 if var == "full" else 0.0
        cells = []
        for n in models:
            pred = gemm[n](B) + kv
            e = pred / meas - 1
            errs.setdefault((n, var), []).append(e)
            cells.append(f"{pred:10.1f} {100 * e:+5.1f}%")
        print(f"{B:4d} {L:5d} {var:9} {meas:9.1f} {kv:7.1f} | " + " | ".join(cells))
    for (n, var), v in errs.items():
        print(f"{var:9} {n:10}: mean {100 * statistics.mean(v):+.1f}%  MAPE {100 * statistics.mean(map(abs, v)):.1f}%"
              f"  (n={len(v)})")
    print("\nattention alone: measured / KV read at 819 GB/s")
    for r in d["attention"]:
        print(f"  B {r['B']:4d} L {r['L']:6d} {r.get('mode', ''):9} {r['us']:9.1f} us  {r['kv_gbps']:6.0f} GB/s"
              f"  x{r['vs_roofline_819']:.2f}")


if __name__ == "__main__":
    main()
