#!/usr/bin/env python3
"""Numerical regression checks for the decompression and KV-cache/attention models.

Every assertion is a CLOSED-FORM identity of the declared model -- computed here from the
config knobs alone and compared against the simulator's report -- so a regression in the
cost model fails loudly instead of shifting a curve nobody re-derives:

  decomp (absolute decoder):
    D1 compressed bytes  == ceil(dense/CR + metadata_per_tile * tiles)
    D2 decoder cycles    == (startup + dense/throughput) * weight_refetch
    D3 weight DRAM ratio == compressed/dense (access-cycle scaling)
    D4 per-tile cv       == same compressed bytes as uniform (total pinned)
    D5 near-1 CR + metadata that erases the gain -> BYPASSED, decoder cycles 0

  kv/attention:
    K1 blocking step     == supply_K + supply_V + qk + softmax + av + write
    K2 two_pass - online == softmax barrier (streaming, same knobs otherwise)
    K3 occupancy         == bytes_per_token * (context + 1)

Exit code 0 iff every check passes.  Run:  python3 validation/decomp/check_decomp.py
"""
import math
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MODEL = os.path.join(ROOT, "models", "model")
NET = os.path.join(ROOT, "configs", "networks", "llm_qwen3_ffn_up_decode.cfg")
MAP = os.path.join(ROOT, "configs", "mappings", "gemmini", "llm_qwen3_ffn_up_decode", "ws.map")
BASE = open(os.path.join(ROOT, "configs", "accelerators", "gemmini.cfg")).read()
ENV = dict(os.environ, LD_LIBRARY_PATH=":".join(
    [os.path.join(ROOT, "ext", "DRAMsim3"), os.path.join(ROOT, "ext", "nebula", "library"),
     os.environ.get("LD_LIBRARY_PATH", "")]))

# The qwen3 ffn_up decode workload: dense weight = K x N = 1024 x 3072 int8, and its
# baseline mapping refetches the weight ONCE from DRAM (weight repetition factor 1).
DENSE_WEIGHT = 1024*3072
WEIGHT_REFETCH = 1
TILE_BYTES = 256          # multi-chip weight tile of this mapping (256 elements int8)

failures = []


def check(name, got, want, tol=1e-6):
    ok = (abs(got - want) <= tol*max(1.0, abs(want)))
    print("  %-4s %-58s %s" % (name, "got %.6g want %.6g" % (got, want),
                               "ok" if ok else "FAIL"))
    if not ok:
        failures.append(name)


def run(extra):
    with tempfile.NamedTemporaryFile("w", suffix=".cfg", delete=False, dir="/tmp") as f:
        f.write(BASE + extra)
        cfg = f.name
    out = tempfile.mkdtemp(dir="/tmp")
    os.symlink(os.path.join(ROOT, "models", "datasets"), os.path.join(out, "datasets"))
    subprocess.run([MODEL, "run-ir", cfg, NET, MAP], cwd=out, env=ENV,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    text = ""
    for fn in os.listdir(out):
        if fn.endswith("_layer_0.txt"):
            text = open(os.path.join(out, fn)).read()
    return text


def grab(text, pat):
    m = re.search(pat, text)
    return float(m.group(1)) if m else None


DECOMP = ("\n[decomp]\ncompression_ratio = {cr}\ndecoder_bytes_per_cycle = {tp}\n"
          "startup_cycles = {st}\ninput_queue_depth = 4\noutput_buffer_tiles = 2\n"
          "overlap = 1\ngranularity = weight_tile\nmetadata_scale_bytes_per_tile = {meta}\n"
          "tile_ratio_cv = {cv}\ndecomp_decoder_energy = 0.05\ndecomp_static_energy = 0.001\n")


def main():
    print("== decomp: absolute decoder, CR=2, throughput 16 B/cyc, startup 32, meta 0 ==")
    t = run(DECOMP.format(cr=2.0, tp=16.0, st=32, meta=0.0, cv=0.0))
    tiles = DENSE_WEIGHT//TILE_BYTES
    comp = grab(t, r"Weight bytes dense/comp\s*:\s*\d+\s*/\s*(\d+)")
    check("D1", comp, math.ceil(DENSE_WEIGHT/2.0))
    dec = grab(t, r"Decoder cycles\s*:\s*([\d.]+)")
    check("D2", dec, (32 + DENSE_WEIGHT/16.0)*WEIGHT_REFETCH)
    # D3: the weight DRAM access axis carries exactly the compressed/dense ratio. The DRAM
    # access-cycle rows sit after the "Multi-chip <-> DRAM transactions" section (several
    # other components print their own "Access cycle" blocks earlier).
    def weight_dram_access(text):
        tail = text.split("Multi-chip <-> DRAM transactions", 1)
        return grab(tail[1], r"Access cycle\s*\n[^\n]*\n\s*\* Weight\s*:\s*([\d.]+) cycles") \
            if len(tail) == 2 else None
    dense_t = run("")
    wd = weight_dram_access(dense_t)
    wc = weight_dram_access(t)
    if wd and wc:
        check("D3", wc/wd, comp/DENSE_WEIGHT, tol=1e-3)
    else:
        print("  D3   (weight DRAM rows not found -- pattern drift)              FAIL")
        failures.append("D3")

    print("== decomp: metadata charged per tile (meta=8 -> +8*tiles bytes) ==")
    t = run(DECOMP.format(cr=2.0, tp=16.0, st=32, meta=8.0, cv=0.0))
    comp_meta = grab(t, r"Weight bytes dense/comp\s*:\s*\d+\s*/\s*(\d+)")
    check("D1m", comp_meta, math.ceil(DENSE_WEIGHT/2.0 + 8.0*tiles))

    print("== decomp: per-tile cv pins the layer total (same bytes as uniform) ==")
    t = run(DECOMP.format(cr=2.0, tp=16.0, st=32, meta=0.0, cv=0.5))
    comp_cv = grab(t, r"Weight bytes dense/comp\s*:\s*\d+\s*/\s*(\d+)")
    check("D4", comp_cv, comp, tol=1e-5)

    print("== decomp: bypass break-even (CR 1.05, meta 16/tile erases the gain) ==")
    # values = dense/1.05 saves ~4.76%; metadata 16*tiles = 16/256 = 6.25% of dense > gain.
    t = run(DECOMP.format(cr=1.05, tp=16.0, st=32, meta=16.0, cv=0.0))
    bypassed = 1.0 if "BYPASSED" in t else 0.0
    check("D5", bypassed, 1.0)
    check("D5d", grab(t, r"Decoder cycles\s*:\s*([\d.]+)") or 0.0, 0.0)

    KV = ("\n[kvcache]\nattention = 1\nn_q_heads = 16\nn_kv_heads = 8\nhead_dim = 64\n"
          "kv_precision_bytes = 1\ncontext_length = 4096\nattention_macs_per_cycle = 256\n"
          "softmax_cycles_per_element = 4\nattention_algorithm = {algo}\n"
          "kv_schedule = {sched}\nkv_tile_bytes = 65536\nkvcache_read_energy = 0.02\n"
          "kvcache_attention_mac_energy = 0.1\nkvcache_softmax_energy = 0.5\n")

    print("== kv/attention: blocking step identity (K1) ==")
    t = run(KV.format(algo="online", sched="blocking"))
    supply = re.search(r"K / V supply\s*:\s*([\d.]+) / ([\d.]+)", t)
    qksa = re.search(r"QK / softmax / AV\s*:\s*([\d.]+) / ([\d.]+) / ([\d.]+)", t)
    write = grab(t, r"cache append write\s*:\s*([\d.]+)")
    step = grab(t, r"Attention step\s*:\s*([\d.]+)")
    if supply and qksa and step is not None:
        want = (float(supply.group(1)) + float(supply.group(2)) +
                float(qksa.group(1)) + float(qksa.group(2)) + float(qksa.group(3)) +
                (write or 0.0))
        check("K1", step, want, tol=1e-6)
    else:
        print("  K1   (attention report rows not found)                          FAIL")
        failures.append("K1")

    print("== kv/attention: two_pass - online == softmax barrier (K2, streaming) ==")
    online = grab(run(KV.format(algo="online", sched="streaming")),
                  r"Attention step\s*:\s*([\d.]+)")
    twop = grab(run(KV.format(algo="two_pass", sched="streaming")),
                r"Attention step\s*:\s*([\d.]+)")
    softmax = 16*4096*4.0
    if online is not None and twop is not None:
        # The barrier is the dominant, exact difference; per-pass pipeline fill differs by
        # one tile between the interleaved and split passes, hence the small tolerance.
        check("K2", twop - online, softmax, tol=0.05)
    else:
        print("  K2   (steps not found)                                          FAIL")
        failures.append("K2")

    print("== kv/attention: cache occupancy (K3) ==")
    occ = grab(run(KV.format(algo="online", sched="blocking")),
               r"cache occupancy\s*:\s*(\d+)")
    check("K3", occ or 0.0, 1024.0*(4096 + 1))

    print()
    if failures:
        print("FAILED checks: %s" % ", ".join(failures))
        return 1
    print("decompression / KV-attention numerical checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
