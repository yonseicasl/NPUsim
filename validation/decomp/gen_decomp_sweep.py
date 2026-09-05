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

Experiments (context axis -- KV interaction):
  A. Crossover           : Qwen3 FFN-up, context x {dense, weightCR4, kvCR4, bothCR4}
  B. Generality          : {Qwen3,TinyLlama,Gemma} x context x {dense, kvCR4}
  C. Weight x KV heatmap : Qwen3 at a KV-dominant context, weight_CR x kv_CR grid
  D. KV-decoder ablation : Qwen3 at a KV-dominant context, kvCR4, decoder throughput sweep

Experiments (weight-decompression axis -- CS2, evaluation discussion 2026-09-03 Sec 1;
no [kvcache], DRAM bandwidth modeled by scaling read/write/transfer_cycle, at a
bandwidth-constrained point where the weight supply is the bottleneck):
  W1. CR x decoder heatmap : weight CR x weight decoder_ratio (break-even + decoder wall)
  W2. Bandwidth transition : weight CR x DRAM bandwidth (DRAM -> decoder -> compute)
  W3. Queue x overlap      : input_queue_depth x overlap on/off
  W4. Startup sensitivity  : decoder startup_cycles sweep
  W5. Bypass break-even    : near-1 CR x metadata_scale_bytes_per_tile (dense bypass edge)

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


def _scale_cycle_value(rhs, factor):
    """Scale an 'a:b:c' or single-value cycle spec by factor (rounded, >= 1)."""
    return ":".join(str(max(1, int(round(float(p)*factor)))) for p in rhs.split(":"))


def scale_dram_bandwidth(base_text, bw_scale):
    """Model DRAM bandwidth by the per-access cost: scale the [dram] section's read_cycle /
    write_cycle / transfer_cycle by 1/bw_scale (the `bandwidth=` line does not drive the
    critical path). bw_scale = 1.0 keeps the RTL-calibrated cost."""
    if bw_scale == 1.0:
        return base_text
    factor = 1.0/bw_scale
    out, section = [], None
    for line in base_text.splitlines(keepends=True):
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1]
        if section == "dram":
            key = stripped.replace(" ", "")
            for name in ("read_cycle=", "write_cycle=", "transfer_cycle="):
                if key.startswith(name):
                    rhs = stripped.split("=", 1)[1].strip()
                    line = "%s = %s\n" % (name[:-1], _scale_cycle_value(rhs, factor))
                    break
        out.append(line)
    return "".join(out)


def decomp_section(weight_cr, dr=None, queue=None, overlap=None, startup=None,
                   metadata=0.0):
    if weight_cr <= 1.0 and metadata <= 0.0:
        return ""
    return ("\n[decomp]\ncompression_ratio = %g\ndecoder_ratio = %g\nstartup_cycles = %d\n"
            "input_queue_depth = %d\noutput_buffer_tiles = 2\noverlap = %d\n"
            "granularity = weight_tile\nmetadata_scale_bytes_per_tile = %g\n"
            "decomp_decoder_energy = %g\ndecomp_static_energy = %g\n" % (
                weight_cr,
                W_DECODER_RATIO if dr is None else dr,
                W_STARTUP if startup is None else startup,
                W_QUEUE if queue is None else queue,
                W_OVERLAP if overlap is None else overlap,
                metadata, W_DECODER_ENERGY, W_STATIC_ENERGY))


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


def write_weight_accel(name, weight_cr, bw=1.0, dr=None, queue=None, overlap=None,
                       startup=None, metadata=0.0):
    """Weight-decompression-only config (no [kvcache]): base Gemmini with the DRAM cost
    scaled to model bandwidth, plus the [decomp] section with the given knobs."""
    base = scale_dram_bandwidth(open(BASE_ACCEL).read(), bw)
    text = base + decomp_section(weight_cr, dr=dr, queue=queue, overlap=overlap,
                                 startup=startup, metadata=metadata)
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

    def add(experiment, name, workload, model, context, weight_cr, kv_cr, kv_decoder,
            variant, bw=1.0, dr="", queue="", overlap="", startup="", metadata=""):
        runs.append({
            "run_id": "%s__%s" % (experiment, name),
            "experiment": experiment, "variant": variant,
            "accel_cfg": os.path.join(OUT_DIR, name + ".cfg"),
            "network_cfg": net_path(workload), "mapping": map_path(workload),
            "workload": workload, "model": model, "context": context,
            "weight_cr": weight_cr, "kv_cr": kv_cr, "kv_decoder": kv_decoder,
            # "dense" here = the no-compression baseline AT this operating point.
            "is_dense": int(variant == "dense"),
            "bw": bw, "dr": dr, "queue": queue, "overlap": overlap,
            "startup": startup, "metadata": metadata,
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

    # ======== CS2: weight-decompression axis (no [kvcache]) ====================
    # The weight supply must be the bottleneck for these axes to bite, so W1/W3/W4/W5 run
    # at a bandwidth-constrained point (W_BW); W2 sweeps the bandwidth itself across the
    # PE-fill -> DRAM transition.
    W_BW = 0.05
    wl, _ = MODELS[MAIN]

    def w_dense(bw):
        name = "decomp_W_dense_bw%s" % num_tag(bw)
        if not os.path.exists(os.path.join(OUT_DIR, name + ".cfg")):
            write_weight_accel(name, 1.0, bw=bw)
            add("W", name, wl, MAIN, 0, 1.0, 1.0, 0.0, "dense", bw=bw)
        return name

    # ---- W1. CR x weight-decoder heatmap (@ W_BW) ----------------------------
    w_dense(W_BW)
    for cr in [1.25, 1.5, 2.0, 3.0, 4.0]:
        for dr in [0.25, 0.5, 1.0, 2.0, 4.0]:
            name = "decomp_W1_cr%s_dr%s" % (num_tag(cr), num_tag(dr))
            write_weight_accel(name, cr, bw=W_BW, dr=dr)
            add("W1", name, wl, MAIN, 0, cr, 1.0, 0.0, "weight", bw=W_BW, dr=dr)

    # ---- W2. Bandwidth transition: CR x DRAM bandwidth -----------------------
    for bw in [1.0, 0.5, 0.25, 0.1, 0.05, 0.025]:
        w_dense(bw)
        for cr in [2.0, 4.0]:
            name = "decomp_W2_cr%s_bw%s" % (num_tag(cr), num_tag(bw))
            write_weight_accel(name, cr, bw=bw)
            add("W2", name, wl, MAIN, 0, cr, 1.0, 0.0, "weight", bw=bw, dr=W_DECODER_RATIO)

    # ---- W3. Queue depth x overlap (@ W_BW, CR2, decoder just keeping pace) ---
    for queue in [1, 2, 4, 8]:
        for overlap in [0, 1]:
            name = "decomp_W3_q%d_ov%d" % (queue, overlap)
            write_weight_accel(name, 2.0, bw=W_BW, queue=queue, overlap=overlap)
            add("W3", name, wl, MAIN, 0, 2.0, 1.0, 0.0, "weight", bw=W_BW,
                dr=W_DECODER_RATIO, queue=queue, overlap=overlap)

    # ---- W4. Startup-latency sensitivity (@ W_BW, CR2) ------------------------
    for startup in [0, 32, 128, 512, 2048]:
        name = "decomp_W4_st%d" % startup
        write_weight_accel(name, 2.0, bw=W_BW, startup=startup)
        add("W4", name, wl, MAIN, 0, 2.0, 1.0, 0.0, "weight", bw=W_BW,
            dr=W_DECODER_RATIO, startup=startup)

    # ---- W5. Dense-bypass break-even: near-1 CR x per-tile metadata (@ W_BW) --
    # With metadata charged PER TILE, a low CR can compress below dense on values yet
    # exceed it in total -> the engine bypasses to dense. This sweep walks that edge.
    for cr in [1.02, 1.05, 1.1, 1.25, 1.5]:
        for meta in [0.0, 8.0, 32.0]:
            name = "decomp_W5_cr%s_m%s" % (num_tag(cr), num_tag(meta))
            write_weight_accel(name, cr, bw=W_BW, metadata=meta)
            add("W5", name, wl, MAIN, 0, cr, 1.0, 0.0, "weight", bw=W_BW,
                dr=W_DECODER_RATIO, metadata=meta)

    cols = ["run_id", "experiment", "variant", "accel_cfg", "network_cfg", "mapping",
            "workload", "model", "context", "weight_cr", "kv_cr", "kv_decoder", "is_dense",
            "bw", "dr", "queue", "overlap", "startup", "metadata"]
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
