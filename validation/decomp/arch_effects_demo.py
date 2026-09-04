#!/usr/bin/env python3
"""The two architecture-level effects an analytical layer-mean model cannot represent
(evaluation discussion 2026-09-03, Secs 5 and 7). Both experiments hold the ANALYTICAL
inputs fixed -- same average compression ratio and DRAM traffic, same KV traffic -- and
vary only what the timeline sees.

  1. Per-tile compression variation ([decomp] tile_ratio_cv): identical layer-average CR
     and identical DRAM bytes; only the per-tile compressed sizes (bursty arrival into the
     bounded decoder queue) differ -> the decoder exposure and critical path move.

  2. KV supply scheduling ([kvcache] kv_schedule): identical KV traffic; only the schedule
     (blocking / streaming / double_buffered, bounded KV buffer) differs -> the exposed vs
     hidden KV latency and critical path move. SCOPE: this is a supply-scheduling PROXY --
     the consumer is the mapped layer's window, not attention execution (QK/softmax/AV are
     not modeled; see components/kvcache.h for the full non-modeled list).

Run:  python3 validation/decomp/arch_effects_demo.py
"""
import os
import re
import subprocess
import tempfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MODEL = os.path.join(ROOT, "models", "model")
NET = os.path.join(ROOT, "configs", "networks", "llm_qwen3_ffn_up_decode.cfg")
MAP = os.path.join(ROOT, "configs", "mappings", "gemmini", "llm_qwen3_ffn_up_decode", "ws.map")
BASE = open(os.path.join(ROOT, "configs", "accelerators", "gemmini.cfg")).read()
ENV = dict(os.environ, LD_LIBRARY_PATH=":".join(
    [os.path.join(ROOT, "ext", "DRAMsim3"), os.path.join(ROOT, "ext", "nebula", "library"),
     os.environ.get("LD_LIBRARY_PATH", "")]))


def run(accel_text):
    with tempfile.NamedTemporaryFile("w", suffix=".cfg", delete=False, dir="/tmp") as f:
        f.write(accel_text)
        cfg = f.name
    out = tempfile.mkdtemp(dir="/tmp")
    os.symlink(os.path.join(ROOT, "models", "datasets"), os.path.join(out, "datasets"))
    subprocess.run([MODEL, "run-ir", cfg, NET, MAP], cwd=out, env=ENV,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    text = ""
    for fn in os.listdir(out):
        if fn.endswith("_layer_0.txt"):
            text = open(os.path.join(out, fn)).read()
    def g(pat):
        m = re.search(pat, text)
        return float(m.group(1)) if m else None
    return {
        "crit": g(r"Critical-path latency\s*:\s*([\d.]+)"),
        "comp_bytes": g(r"Weight bytes dense/comp\s*:\s*\d+\s*/\s*(\d+)"),
        "dec_exposed": g(r"\* on critical path\s*:\s*([\d.]+)"),
        "dec_stall": None,
        "kv_exposed": g(r"exposed on crit path\s*:\s*([\d.]+)"),
        "kv_hidden": g(r"hidden behind layer\s*:\s*([\d.]+)"),
    }


DECOMP = ("\n[decomp]\ncompression_ratio = 2.0\ndecoder_ratio = 1.0\nstartup_cycles = 32\n"
          "input_queue_depth = 2\noutput_buffer_tiles = 2\noverlap = 1\n"
          "granularity = weight_tile\ntile_ratio_cv = %g\n"
          "decomp_decoder_energy = 0.05\ndecomp_static_energy = 0.001\n")

KV = ("\n[kvcache]\nkv_bytes_per_token = 1024\ncontext_length = 32768\n"
      "kvcache_read_energy = 0.02\nkv_schedule = %s\n%s")


def main():
    print("=" * 76)
    print("1. Per-tile compression variation (same avg CR=2, same DRAM bytes)")
    print("   [Qwen3 FFN-up decode, decoder_ratio=1, queue=2]")
    print("=" * 76)
    print("  %-8s %14s %14s %14s" % ("cv", "comp_bytes", "crit_path", "dec_exposed"))
    base_crit = None
    for cv in (0.0, 0.3, 0.6):
        r = run(BASE + DECOMP % cv)
        if base_crit is None:
            base_crit = r["crit"]
        delta = (r["crit"]/base_crit - 1.0)*100 if base_crit else 0.0
        print("  %-8g %14.0f %14.1f %14.1f   (%+.2f%% vs uniform)" % (
            cv, r["comp_bytes"] or 0, r["crit"] or 0, r["dec_exposed"] or 0, delta))
    print("  -> identical analytical inputs; only the arrival pattern differs.\n")

    print("=" * 76)
    print("2. KV supply scheduling (same KV traffic, ctx=32768)")
    print("   [tile 256KB where scheduled; buffer default: streaming/db=2]")
    print("=" * 76)
    print("  %-22s %14s %14s %14s" % ("schedule", "crit_path", "kv_exposed", "kv_hidden"))
    tile = "kv_tile_bytes = 262144\n"
    for sched, extra in (("aggregate", ""), ("blocking", ""), ("streaming", tile),
                         ("streaming", tile + "kv_buffer_tiles = 1\n"),
                         ("double_buffered", tile)):
        r = run(BASE + KV % (sched, extra))
        label = sched + (" (buf=1)" if "kv_buffer_tiles = 1" in extra else "")
        def cell(v):
            return "%14.1f" % v if v is not None else "%14s" % "(dram axis)"
        print("  %-22s %14.1f %s %s" % (
            label, r["crit"] or 0, cell(r["kv_exposed"]), cell(r["kv_hidden"])))
    print("  -> same traffic, same roofline bound; the schedule decides what is exposed.")


if __name__ == "__main__":
    main()
