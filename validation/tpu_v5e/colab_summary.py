#!/usr/bin/env python3
"""Summary of tpu_v5e1_stream_bench.py results: cold vs warm and int8 GEMMs, decode layers vs their GEMMs (and NPUsim's,
if compare_v5e.py has written compare_gemm.csv in the same directory), attention and HBM read bandwidth.

Usage: python3 validation/tpu_v5e/colab_summary.py RESULTS_DIR [--isolated OTHER_RESULTS_DIR]
--isolated: a measure_tpu.py run (one GEMM per program) to compare the streamed GEMMs with.
"""
import argparse
import csv
import json
import os

HBM = 819e9


def ratio(a, b):
    return a / b if a and b else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--isolated")
    a = ap.parse_args()
    r = json.load(open(os.path.join(a.results, "results.json")))
    print(r["env"], "errors:", len(r["errors"]))
    npu_path = os.path.join(a.results, "compare_gemm.csv")
    npu = {(int(x["M"]), int(x["K"]), int(x["N"])): float(x["npusim_us"]) for x in csv.DictReader(open(npu_path))} \
        if os.path.exists(npu_path) else {}
    iso = {}
    if a.isolated:
        iso = {(int(x["M"]), int(x["K"]), int(x["N"])): float(x["device_us"])
               for x in csv.DictReader(open(os.path.join(a.isolated, "gemm.csv"))) if x["device_us"]}

    g = {(x["M"], x["K"], x["N"], x["weights"], x["cold"]): x for x in r["stream_gemm"]}
    print(f"\n{'M':>4} {'K':>6} {'N':>6} | {'cold':>8} {'warm':>8} {'int8':>8} us | {'GB/s':>6} {'warm/cold':>9}"
          f" {'int8/bf16':>9} | {'isolated':>8} {'NPUsim':>8} {'err':>7}")
    for (M, K, N, wd, cold), x in g.items():
        if wd != "bf16" or not cold:
            continue
        warm, i8 = g.get((M, K, N, "bf16", False)), g.get((M, K, N, "int8", True))
        s = npu.get((M, K, N))
        print(f"{M:4d} {K:6d} {N:6d} | {x['us']:8.1f} {warm['us'] if warm else 0:8.1f} {i8['us'] if i8 else 0:8.1f}    | "
              f"{x['weight_gbps'] or 0:6.0f} {ratio(warm and warm['us'], x['us']):9.2f} {ratio(i8 and i8['us'], x['us']):9.2f} | "
              f"{iso.get((M, K, N), 0):8.1f} {s or 0:8.1f} {100 * (ratio(s, x['us']) - 1) if s else 0:+6.1f}%")

    print(f"\n{'B':>4} {'L':>5} | {'layer':>8} {'GEMM-only':>9} {'4 GEMMs':>8} {'NPUsim':>8} {'err':>7} | {'non-GEMM':>8}  (us/layer)")
    d = {(x["B"], x["L"], x["variant"]): x for x in r["decode_stack"]}
    for (B, L, var), x in d.items():
        if var != "full":
            continue
        go = d.get((B, L, "gemm_only"))
        projs = [k for k in g if k[0] == B and k[3] == "bf16" and k[4]]
        gsum = sum(g[k]["us"] for k in projs) if len(projs) == 4 else 0
        ssum = sum(npu.get(k[:3], 0) for k in projs) if npu and len(projs) == 4 else 0
        ref = go["us_per_layer"] if go else 0
        print(f"{B:4d} {L:5d} | {x['us_per_layer']:8.1f} {ref:9.1f} {gsum:8.1f} {ssum:8.1f} "
              f"{100 * (ratio(ssum, ref) - 1) if ssum and ref else 0:+6.1f}% | {x['us_per_layer'] - ref:8.1f}")
    print("GEMM-only: the layer without norms, attention and SwiGLU (its four projections + residual adds);"
          " NPUsim err is against it")

    print(f"\n{'B':>4} {'L':>6} | {'us':>8} {'KV GB/s':>8} {'x roofline(819)':>15}")
    for x in r["attention"]:
        print(f"{x['B']:4d} {x['L']:6d} | {x['us']:8.1f} {x['kv_gbps'] or 0:8.0f} {x['vs_roofline_819'] or 0:15.2f}")
    print(f"\n{'MB':>6} | {'read GB/s':>9} {'of 819':>7}")
    for x in r["hbm_read"]:
        print(f"{x['mbytes']:6d} | {x['read_gbps'] or 0:9.0f} {100 * (x['read_gbps'] or 0) * 1e9 / HBM:6.1f}%")
    if r["odd_gemm"]:
        print(f"\n{'M':>4} {'K=N':>5} | {'us':>8} {'TFLOPS':>7} {'NPUsim':>8} {'err':>7}")
        for x in r["odd_gemm"]:
            s = npu.get((x["M"], x["K"], x["N"]))
            print(f"{x['M']:4d} {x['K']:5d} | {x['us']:8.1f} {x['tflops'] or 0:7.1f} {s or 0:8.1f} "
                  f"{100 * (ratio(s, x['us']) - 1) if s else 0:+6.1f}%")


if __name__ == "__main__":
    main()
