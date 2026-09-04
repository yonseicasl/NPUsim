#!/usr/bin/env python3
"""Section D (compression) sweep configs, parameterized by KV-cache CONTEXT LENGTH.

Redesign rationale (see validation/decomp/kv_crossover.py): at the RTL-calibrated Gemmini
DRAM, a decode projection is PE-array-fill-bound; what actually makes decode memory-bound
is the KV-cache DRAM read, which grows with context length. So the sweep's memory-bound
axis is CONTEXT LENGTH (realistic), not an artificial DRAM-bandwidth knob. Two independent
compression components are studied together:
  - weight decompression ([decomp]): shrinks WEIGHT DRAM traffic.
  - KV compression ([kvcache] kv_compression_ratio): shrinks the KV-cache read.

Every point reuses the EXISTING decode network + mapping files and is launched with
`model run-ir <accel.cfg> <net.cfg> <map.map>`. The KV footprint per token is a model
property (2 * n_kv_heads * head_dim * precision), so it lives in each accel config.

Experiments:
  A. Crossover           : Qwen3 FFN-up, context x {dense, weightCR4, kvCR4, bothCR4}
  B. Generality          : {Qwen3,TinyLlama,Gemma} x context x {dense, kvCR4}
  C. Weight x KV heatmap : Qwen3 at a KV-dominant context, weight_CR x kv_CR grid
  D. KV-decoder ablation : Qwen3 at a KV-dominant context, kvCR4, decoder throughput sweep

Writes configs/accelerators/decomp_sweep/*.cfg and validation/decomp/runs.csv.
"""
import csv
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BASE_ACCEL = os.path.join(ROOT, "configs", "accelerators", "gemmini.cfg")
OUT_DIR = os.path.join(ROOT, "configs", "accelerators", "decomp_sweep")
NET_DIR = os.path.join(ROOT, "configs", "networks")
MAP_DIR = os.path.join(ROOT, "configs", "mappings", "gemmini")
MANIFEST = os.path.join(os.path.dirname(__file__), "runs.csv")

# Weight-decompression placeholders (uncalibrated). decoder_ratio 1 = decoder keeps pace.
W_DECODER_RATIO = 1.0
W_STARTUP = 32
W_QUEUE = 4
W_OVERLAP = 1
W_DECODER_ENERGY = 0.05
W_STATIC_ENERGY = 0.001
KV_READ_ENERGY = 0.02

# int8 per-token KV footprint (K+V = kv_out_combined) and a representative decode projection.
MODELS = {
    "qwen3":     ("llm_qwen3_ffn_up_decode",      1024),   # GQA kv_out=1024
    "tinyllama": ("llm_tinyllama_ffn_up_decode",   256),   # aggressive GQA
    "gemma":     ("llm_gemma_slide_ffn_up_decode", 256),
}
MAIN = "qwen3"
CONTEXTS = [512, 2048, 8192, 16384, 32768, 65536, 131072]
MEM_BOUND_CTX = 32768                                       # KV-dominant operating point


def num_tag(x):
    return ("%g" % x).replace(".", "p")


def decomp_section(weight_cr):
    if weight_cr <= 1.0:
        return ""
    return ("\n[decomp]\ncompression_ratio = %g\ndecoder_ratio = %g\nstartup_cycles = %d\n"
            "input_queue_depth = %d\noverlap = %d\ngranularity = weight_tile\n"
            "decomp_decoder_energy = %g\ndecomp_static_energy = %g\n" % (
                weight_cr, W_DECODER_RATIO, W_STARTUP, W_QUEUE, W_OVERLAP,
                W_DECODER_ENERGY, W_STATIC_ENERGY))


def kvcache_section(kv_bpt, context, kv_cr=1.0, kv_decoder=0.0):
    text = ("\n[kvcache]\nkv_bytes_per_token = %d\ncontext_length = %d\n"
            "kv_compression_ratio = %g\nkvcache_read_energy = %g\n" % (
                kv_bpt, context, kv_cr, KV_READ_ENERGY))
    if kv_decoder > 0.0:
        text += "kv_decoder_bytes_per_cycle = %g\n" % kv_decoder
    return text


def write_accel(name, kv_bpt, context, weight_cr=1.0, kv_cr=1.0, kv_decoder=0.0):
    base = open(BASE_ACCEL).read()
    text = base + decomp_section(weight_cr) + kvcache_section(kv_bpt, context, kv_cr, kv_decoder)
    with open(os.path.join(OUT_DIR, name + ".cfg"), "w") as f:
        f.write(text)


def net_path(wl):  return os.path.join(NET_DIR, wl + ".cfg")
def map_path(wl):  return os.path.join(MAP_DIR, wl, "ws.map")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for f in os.listdir(OUT_DIR):
        if f.endswith(".cfg"):
            os.remove(os.path.join(OUT_DIR, f))
    runs = []

    def add(experiment, name, workload, model, context, weight_cr, kv_cr, kv_decoder, variant):
        runs.append({
            "run_id": "%s__%s" % (experiment, name),
            "experiment": experiment, "variant": variant,
            "accel_cfg": os.path.join(OUT_DIR, name + ".cfg"),
            "network_cfg": net_path(workload), "mapping": map_path(workload),
            "workload": workload, "model": model, "context": context,
            "weight_cr": weight_cr, "kv_cr": kv_cr, "kv_decoder": kv_decoder,
            # "dense" here = the no-compression baseline AT this context (KV read present).
            "is_dense": int(variant == "dense"),
        })

    # ---- A. Crossover: Qwen3 FFN-up, context x compressor -------------------
    wl, bpt = MODELS[MAIN]
    for ctx in CONTEXTS:
        for variant, wcr, kcr in [("dense", 1.0, 1.0), ("weight", 4.0, 1.0),
                                  ("kv", 1.0, 4.0), ("both", 4.0, 4.0)]:
            name = "decomp_A_%s_ctx%d_%s" % (MAIN, ctx, variant)
            write_accel(name, bpt, ctx, weight_cr=wcr, kv_cr=kcr)
            add("A", name, wl, MAIN, ctx, wcr, kcr, 0.0, variant)

    # ---- B. Generality: models x context x {dense, kvCR4} -------------------
    for model, (wl, bpt) in MODELS.items():
        for ctx in [2048, 8192, 32768, 131072]:
            for variant, kcr in [("dense", 1.0), ("kv", 4.0)]:
                name = "decomp_B_%s_ctx%d_%s" % (model, ctx, variant)
                write_accel(name, bpt, ctx, kv_cr=kcr)
                add("B", name, wl, model, ctx, 1.0, kcr, 0.0, variant)

    # ---- C. Weight x KV compression heatmap (KV-dominant context) -----------
    wl, bpt = MODELS[MAIN]
    for wcr in [1.0, 2.0, 4.0]:
        for kcr in [1.0, 2.0, 4.0]:
            variant = "dense" if (wcr == 1.0 and kcr == 1.0) else "grid"
            name = "decomp_C_w%s_kv%s" % (num_tag(wcr), num_tag(kcr))
            write_accel(name, bpt, MEM_BOUND_CTX, weight_cr=wcr, kv_cr=kcr)
            add("C", name, wl, MAIN, MEM_BOUND_CTX, wcr, kcr, 0.0, variant)

    # ---- D. KV-decoder ablation (KV-dominant context, kvCR4) ----------------
    # Decompression is not free: sweep the KV decoder throughput (dense B/cycle). 0 = ideal.
    wl, bpt = MODELS[MAIN]
    # dense baseline (no compression) for the speedup reference at this context.
    write_accel("decomp_D_dense", bpt, MEM_BOUND_CTX)
    add("D", "decomp_D_dense", wl, MAIN, MEM_BOUND_CTX, 1.0, 1.0, 0.0, "dense")
    for dec in [0.0, 64.0, 16.0, 4.0, 1.0]:
        name = "decomp_D_dec%s" % (num_tag(dec) if dec > 0 else "ideal")
        write_accel(name, bpt, MEM_BOUND_CTX, kv_cr=4.0, kv_decoder=dec)
        add("D", name, wl, MAIN, MEM_BOUND_CTX, 1.0, 4.0, dec, "kv")

    cols = ["run_id", "experiment", "variant", "accel_cfg", "network_cfg", "mapping",
            "workload", "model", "context", "weight_cr", "kv_cr", "kv_decoder", "is_dense"]
    with open(MANIFEST, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(runs)

    by_exp = {}
    for r in runs:
        by_exp[r["experiment"]] = by_exp.get(r["experiment"], 0) + 1
    print("wrote %d accelerator configs to %s" %
          (len([f for f in os.listdir(OUT_DIR) if f.endswith('.cfg')]), OUT_DIR))
    print("wrote %d runs to %s" % (len(runs), MANIFEST))
    for e in sorted(by_exp):
        print("  experiment %s: %d runs" % (e, by_exp[e]))


if __name__ == "__main__":
    main()
