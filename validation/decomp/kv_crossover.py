#!/usr/bin/env python3
"""Weight-vs-KV DRAM crossover + compression effectiveness (evaluation.md Sec 4).

SCOPE: the [kvcache] component is a KV-cache TRAFFIC PROXY (components/kvcache.h) -- it
injects the aggregate KV read volume, with no attention consumer, no KV address stream,
and no KV write/cache state. Every result below is therefore a traffic-sensitivity
analysis, not an architecture-level KV-cache execution study.

Two views of the decode memory-bound picture:

  1. ANALYTICAL -- per decode step, per transformer layer, how the two DRAM traffics scale
     with context length: weight (all projections + FFN, loaded once at M=1) is fixed;
     KV-cache read (2 * n_kv_heads * head_dim * context) grows linearly. Reports the
     crossover context where KV overtakes weight for each model.

  2. SIMULATION -- Qwen3 FFN-up projection with the [kvcache] component across context
     lengths, showing the PE-fill -> KV-DRAM bottleneck transition, and how much a WEIGHT
     compressor (Amdahl-limited by KV) vs a KV compressor (attacks the dominant traffic)
     actually speeds the decode step up.

Run:  python3 validation/decomp/kv_crossover.py
"""
import os
import re
import subprocess
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MODEL = os.path.join(ROOT, "models", "model")
NET = os.path.join(ROOT, "configs", "networks", "llm_qwen3_ffn_up_decode.cfg")
MAP = os.path.join(ROOT, "configs", "mappings", "gemmini", "llm_qwen3_ffn_up_decode", "ws.map")
BASE_ACCEL = os.path.join(ROOT, "configs", "accelerators", "gemmini.cfg")
ENV = dict(os.environ, LD_LIBRARY_PATH=":".join([
    os.path.join(ROOT, "ext", "DRAMsim3"),
    os.path.join(ROOT, "ext", "nebula", "library"),
    os.environ.get("LD_LIBRARY_PATH", ""),
]))

# int8 weight/KV. Per-layer weight = Q,K,V,O + gate,up,down; KV/token = K+V (kv_out combined).
MODELS = {  # name: (d_model, q_out, kv_out, d_ff)
    "Qwen3-0.6B":      (1024, 2048, 1024, 3072),
    "TinyLlama-1.1B":  (2048, 2048,  256, 5632),
    "Gemma-E2B(slide)":(1536, 2048,  256, 6144),
}


def layer_weight_bytes(d, qo, kv, dff):
    return d*qo + d*(kv//2) + d*(kv//2) + qo*d + d*dff + d*dff + dff*d


def analytical_crossover():
    print("=" * 74)
    print("1. ANALYTICAL: weight vs KV-cache DRAM per decode step, per layer (int8)")
    print("=" * 74)
    print("%-18s %11s %9s | KV/weight @ ctx  | crossover ctx" % ("model", "weight(MB)", "KV/tok"))
    print("%-18s %11s %9s | 2K   8K   32K    | (KV > weight)" % ("", "", ""))
    for name, (d, qo, kv, dff) in MODELS.items():
        w = layer_weight_bytes(d, qo, kv, dff)
        kvt = kv
        ratios = [kvt*L/w for L in (2048, 8192, 32768)]
        xover = w/kvt
        print("%-18s %11.2f %9d | %.2f %.2f %.2f | L=%d tokens" % (
            name, w/1e6, kvt, ratios[0], ratios[1], ratios[2], round(xover)))
    print("  -> below the crossover WEIGHT dominates DRAM; above it the KV cache does.")
    print("  -> aggressive GQA (TinyLlama/Gemma, kv=256) keeps KV small: weight dominates")
    print("     to very long context. Qwen3 (kv=1024) crosses in a realistic range.\n")


def run(accel_text):
    with tempfile.NamedTemporaryFile("w", suffix=".cfg", delete=False, dir="/tmp") as f:
        f.write(accel_text)
        path = f.name
    outdir = tempfile.mkdtemp(dir="/tmp")
    os.symlink(os.path.join(ROOT, "models", "datasets"), os.path.join(outdir, "datasets"))
    subprocess.run([MODEL, "run-ir", path, NET, MAP], cwd=outdir, env=ENV,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    layer = None
    for fn in os.listdir(outdir):
        if fn.endswith("_layer_0.txt"):
            layer = os.path.join(outdir, fn)
    text = open(layer).read() if layer else ""
    def grab(pat):
        m = re.search(pat, text)
        return float(m.group(1)) if m else None
    return {
        "crit": grab(r"Critical-path latency\s*:\s*([\d.]+)"),
        "bottleneck": (re.search(r"Bottleneck level\s*:\s*(.+)", text).group(1).strip()
                       if "Bottleneck level" in text else "?"),
        "weight_dram": grab(r"Access cycle\s*\n(?:.*\n)*?\s*\* Weight\s*:\s*([\d.]+) cycles"),
        "kv_supply": grab(r"KV supply DRAM cycles\s*:\s*([\d.]+)"),
    }


def kvsec(ctx, cr=1.0):
    return ("\n[kvcache]\nkv_bytes_per_token = 1024\ncontext_length = %d\n"
            "kv_compression_ratio = %g\nkvcache_read_energy = 0.02\n" % (ctx, cr))


def decompsec(cr):
    return ("\n[decomp]\ncompression_ratio = %g\ndecoder_ratio = 1.0\nstartup_cycles = 32\n"
            "input_queue_depth = 4\noverlap = 1\ngranularity = weight_tile\n"
            "decomp_decoder_energy = 0.05\ndecomp_static_energy = 0.001\n" % cr)


def simulation():
    base = open(BASE_ACCEL).read()
    print("=" * 74)
    print("2. SIMULATION: Qwen3 FFN-up decode, [kvcache] context sweep (RTL-calibrated bw)")
    print("=" * 74)
    print("%-8s %13s %-10s | %8s %8s | %s" % (
        "context", "KVsupply(cyc)", "bottleneck", "wCR4", "kvCR4", "which compressor wins"))
    for ctx in (512, 2048, 8192, 16384, 32768, 65536, 131072):
        dense = run(base + kvsec(ctx))
        wcr = run(base + kvsec(ctx) + decompsec(4.0))
        kcr = run(base + kvsec(ctx, cr=4.0))
        w_sp = dense["crit"]/wcr["crit"] if wcr["crit"] else 0
        k_sp = dense["crit"]/kcr["crit"] if kcr["crit"] else 0
        kv = dense["kv_supply"] or 0
        winner = "weight" if w_sp > k_sp else "KV"
        print("%-8d %13.0f %-10s | %7.3fx %7.3fx | %s" % (
            ctx, kv, dense["bottleneck"], w_sp, k_sp, winner))
    print("\n  weight compression (wCR4): best at SHORT context (weight is the DRAM traffic);")
    print("     Amdahl-limited once KV dominates -> its speedup DECAYS as context grows.")
    print("  KV compression (kvCR4): negligible at short context, but attacks the dominant")
    print("     traffic at long context -> its speedup GROWS with context. The two curves")
    print("     cross where weight vs KV stop/start dominating the DRAM traffic.")


if __name__ == "__main__":
    analytical_crossover()
    simulation()
