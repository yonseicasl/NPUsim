#!/usr/bin/env python3
"""Parse the context-length compression sweep (validation/decomp/out/) into results.csv.

Speedup is critical-path(dense) / critical-path(this), matched to the no-compression run at
the SAME experiment / model / context. Prints per-experiment summaries:
  A crossover (weight vs KV compression across context), B generality, C weight x KV
  heatmap, D KV-decoder ablation.
"""
import csv
import glob
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_ROOT = os.path.join(HERE, "out")
MANIFEST = os.path.join(HERE, "runs.csv")
RESULTS = os.path.join(HERE, "results.csv")


def num(s):
    m = re.search(r"[-+]?\d*\.?\d+", s)
    return float(m.group()) if m else None


def parse_layer(path):
    text = open(path).read()
    def grab(pat):
        m = re.search(pat, text)
        return float(m.group(1)) if m else None
    d = {
        "crit": grab(r"Critical-path latency\s*:\s*([\d.]+)"),
        "kv_supply": grab(r"KV supply DRAM cycles\s*:\s*([\d.]+)"),
        "kv_decoder_cyc": grab(r"decoder reconstitute\s*:\s*([\d.]+)"),
        "kv_compression": grab(r"KV compression\s*:\s*([\d.]+)x"),
        "weight_compression": grab(r"Compression \(dense/comp\):\s*([\d.]+)x"),
    }
    m = re.search(r"Bottleneck level\s*:\s*(.+)", text)
    d["bottleneck"] = m.group(1).strip() if m else "?"
    return d


def find_layer_file(run_dir):
    c = glob.glob(os.path.join(run_dir, "*_layer_0.txt"))
    return c[0] if c else None


def main():
    runs = list(csv.DictReader(open(MANIFEST)))
    rows = []
    for r in runs:
        row = dict(r)
        lf = find_layer_file(os.path.join(OUT_ROOT, r["run_id"]))
        row["parsed"] = 0
        if lf:
            row.update(parse_layer(lf))
            row["parsed"] = 1
        rows.append(row)

    # Dense baseline per (experiment, model, context).
    dense = {}
    for row in rows:
        if int(row.get("is_dense", 0)) and row.get("crit"):
            dense[(row["experiment"], row["model"], row["context"])] = row["crit"]
    for row in rows:
        base = dense.get((row["experiment"], row["model"], row["context"]))
        cp = row.get("crit")
        row["dense_crit"] = base
        row["speedup"] = (base / cp) if (base and cp) else None

    cols = ["run_id", "experiment", "variant", "model", "context", "weight_cr", "kv_cr",
            "kv_decoder", "is_dense", "parsed", "crit", "dense_crit", "speedup",
            "kv_supply", "kv_decoder_cyc", "kv_compression", "weight_compression", "bottleneck"]
    with open(RESULTS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print("parsed %d/%d runs -> %s\n" % (sum(r["parsed"] for r in rows), len(rows), RESULTS))

    def sp(row):
        return row.get("speedup")

    # ---- A: crossover ------------------------------------------------------
    print("=== A. Crossover (Qwen3 FFN-up): speedup vs dense, per context ===")
    print("  %-8s %-10s %8s %8s %8s" % ("context", "bottleneck", "weight", "kv", "both"))
    a = [r for r in rows if r["experiment"] == "A"]
    ctxs = sorted({int(r["context"]) for r in a})
    for ctx in ctxs:
        cell = {r["variant"]: r for r in a if int(r["context"]) == ctx}
        bn = cell.get("dense", {}).get("bottleneck", "?")
        def g(v):
            s = sp(cell.get(v, {})); return ("%.3f" % s) if s else "-"
        print("  %-8d %-10s %8s %8s %8s" % (ctx, bn, g("weight"), g("kv"), g("both")))

    # ---- B: generality -----------------------------------------------------
    print("\n=== B. Generality: KV-compression (CR4) speedup by model x context ===")
    b = [r for r in rows if r["experiment"] == "B" and r["variant"] == "kv"]
    models = sorted({r["model"] for r in b})
    bctx = sorted({int(r["context"]) for r in b})
    print("  %-12s " % "model" + "".join("%9d" % c for c in bctx))
    for m in models:
        cells = []
        for c in bctx:
            hit = [r for r in b if r["model"] == m and int(r["context"]) == c]
            s = sp(hit[0]) if hit else None
            cells.append("%9s" % (("%.3f" % s) if s else "-"))
        print("  %-12s " % m + "".join(cells))

    # ---- C: weight x KV heatmap -------------------------------------------
    print("\n=== C. Weight x KV compression heatmap (Qwen3, ctx=32768): speedup ===")
    c = [r for r in rows if r["experiment"] == "C"]
    wcrs = sorted({float(r["weight_cr"]) for r in c})
    kcrs = sorted({float(r["kv_cr"]) for r in c})
    print("  wCR\\kvCR " + "".join("%8.0f" % k for k in kcrs))
    for wc in wcrs:
        cells = []
        for kc in kcrs:
            hit = [r for r in c if float(r["weight_cr"]) == wc and float(r["kv_cr"]) == kc]
            s = sp(hit[0]) if hit else None
            cells.append("%8s" % (("%.3f" % s) if s else "-"))
        print("  %6.0f   " % wc + "".join(cells))

    # ---- D: KV decoder ablation -------------------------------------------
    print("\n=== D. KV-decoder ablation (Qwen3, ctx=32768, kvCR4): decoder throughput ===")
    d = [r for r in rows if r["experiment"] == "D" and r["variant"] == "kv"]
    print("  %-16s %10s %10s %10s" % ("decoder B/cyc", "decoderCyc", "speedup", "bottleneck"))
    for r in sorted(d, key=lambda x: -float(x["kv_decoder"] if x["kv_decoder"] else 0)):
        dec = float(r["kv_decoder"]) if r["kv_decoder"] else 0
        label = "ideal (0)" if dec == 0 else "%g" % dec
        s = sp(r)
        print("  %-16s %10s %10s %10s" % (
            label, r.get("kv_decoder_cyc", "-"), ("%.3f" % s) if s else "-", r.get("bottleneck", "?")))


if __name__ == "__main__":
    main()
