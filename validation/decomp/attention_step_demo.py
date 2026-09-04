#!/usr/bin/env python3
"""Architecture-Level Impact of KV-Cache Streaming and Buffering (evaluation discussion
2026-09-03 Sec 7 -- the REAL consumer, not the proxy).

One Qwen3-0.6B decode step (16 q-heads / 8 kv-heads / head_dim 64, int8 KV) with the
[kvcache] attention consumer: QK/softmax/AV execute at declared rates, KV tiles carry a
real dependency through the bounded buffer, the KV read/write is a dedicated stream on the
live DRAM model, and the step appends the current token's K/V to the cache.

Every run below moves IDENTICAL KV traffic and computes IDENTICAL attention work -- only
the supply schedule, the buffer, and the softmax structure (online vs two_pass barrier)
differ. The spread between rows is therefore pure architecture-level scheduling effect,
which is exactly what an analytical KV-bytes/bandwidth model cannot show.

Remaining scope limits (stated, not hidden): the step models ONE decode step at the
declared context (no multi-step cache growth), assumes sequential K/V layout for the
row-activation model, and appends after the mapped projection layer (Q dependency).

Run:  python3 validation/decomp/attention_step_demo.py
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

ATTN = """
[kvcache]
attention = 1
n_q_heads = 16
n_kv_heads = 8
head_dim = 64
kv_precision_bytes = 1
context_length = {ctx}
attention_macs_per_cycle = 256
softmax_cycles_per_element = 4
attention_algorithm = {algo}
kv_schedule = {sched}
kv_tile_bytes = 262144
{extra}kvcache_attention_mac_energy = 0.1
kvcache_softmax_energy = 0.5
"""


def run(ctx, algo, sched, extra=""):
    with tempfile.NamedTemporaryFile("w", suffix=".cfg", delete=False, dir="/tmp") as f:
        f.write(BASE + ATTN.format(ctx=ctx, algo=algo, sched=sched, extra=extra))
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
        return float(m.group(1)) if m else 0.0
    return {
        "step": g(r"Attention step\s*:\s*([\d.]+)"),
        "exposed": g(r"KV exposed / hidden\s*:\s*([\d.]+)"),
        "hidden": g(r"KV exposed / hidden\s*:\s*[\d.]+ / ([\d.]+)"),
        "stall": g(r"buffer-full stall\s*:\s*([\d.]+)"),
    }


def main():
    ctx = 32768
    print("=" * 78)
    print("Attention step (REAL consumer), Qwen3 decode, ctx=%d -- identical KV traffic" % ctx)
    print("and attention work in every row; only the schedule/buffer/algorithm differ.")
    print("=" * 78)
    print("  %-34s %13s %13s %13s" % ("schedule / algorithm", "step(cyc)", "KV exposed", "KV hidden"))
    rows = [
        ("blocking / online", "online", "blocking", ""),
        ("streaming(buf2) / online", "online", "streaming", ""),
        ("streaming(buf1) / online", "online", "streaming", "kv_buffer_tiles = 1\n"),
        ("double_buffered / online", "online", "double_buffered", ""),
        ("streaming(buf2) / two_pass", "two_pass", "streaming", ""),
        ("double_buffered / two_pass", "two_pass", "double_buffered", ""),
    ]
    base = None
    for label, algo, sched, extra in rows:
        r = run(ctx, algo, sched, extra)
        if base is None:
            base = r["step"]
        print("  %-34s %13.0f %13.0f %13.0f   (%+.2f%% vs blocking)" % (
            label, r["step"], r["exposed"], r["hidden"],
            (r["step"]/base - 1.0)*100 if base else 0.0))
    print("\n  online hides the softmax under the KV stream; two_pass pays it as a BARRIER")
    print("  between the K pass (QK) and the V pass (AV). The supply schedule then decides")
    print("  how much of the KV stream the attention compute can hide.")


if __name__ == "__main__":
    main()
