#!/usr/bin/env python3
"""measure_tpu_v5e1.py tailprobe: FUNCTIONAL checks only (CPU; no TPU, no timing claim, no hardware validation).
hlo_facts() on HLO texts saved from a real TPU run; op_events() on synthetic traces; pair_verdict() never
concluding more than its inputs show; the --quick run on the CPU backend; the TPU guard for full sizes.
Run with a Python that has jax."""
import argparse
import copy
import csv
import gzip
import json
import math
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.argv = [sys.argv[0]]
import measure_tpu_v5e1 as m

H = os.path.join(HERE, "results", "2026-11-12_colab_v5e1_balance", "hlo")
loop = m.hlo_facts(open(os.path.join(H, "balance_bf16_256x4096_loop.txt")).read())
assert loop["windows"] == [dict(iteration_bounds="11,1,1", kernel_window_bounds="512,3",
                                input_window_bounds="32,32", output_window_bounds="32,3")], loop["windows"]
assert loop["convolutions"] == [("bf16[256,4096]", ("bf16[256,4096]", "bf16[4096,4096]"))], loop["convolutions"]
assert loop["whiles"] >= 1 and loop["fusions"] >= 1 and loop["copy_start"] == 0 and not loop["prefetch"], loop
iso = m.hlo_facts(open(os.path.join(H, "balance_bf16_256x4096_isolated.txt")).read())
assert iso["copy_start"] == 1 and iso["copy_done"] == 1 and iso["prefetch"] and iso["whiles"] == 0, iso
assert "bf16[4096,4096]:S(1)" in iso["tagged"], iso["tagged"]
i8 = m.hlo_facts(open(os.path.join(H, "balance_int8_128x4096_loop.txt")).read())
assert i8["windows"][0]["kernel_window_bounds"] == "512,3" and i8["convolutions"][0][1][1] == "s8[4096,4096]", i8


# op_events(): raw counts from a trace; "verified" only with exactly 8 convolution-fusion events per module
def trace(conv_per_module, modules=2, extra_fusion=0):
    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, "plugins", "profile", "run"))
    ev = [dict(ph="M", name="process_name", pid=3, args=dict(name="/device:TPU:0")),
          dict(ph="M", name="process_name", pid=9, args=dict(name="/host:CPU")),
          dict(ph="M", name="thread_name", pid=3, tid=1, args=dict(name="XLA Modules")),
          dict(ph="M", name="thread_name", pid=3, tid=2, args=dict(name="XLA Ops"))]
    for i in range(modules):
        ev.append(dict(ph="X", pid=3, tid=1, name="jit_loop", ts=1000 * i, dur=500.0))
        ev.append(dict(ph="X", pid=3, tid=2, name="while.4", ts=1000 * i, dur=450.0, args=dict(hlo_category="while")))
        for j in range(conv_per_module):
            ev.append(dict(ph="X", pid=3, tid=2, name="convert_reduce_fusion.2", ts=1000 * i + 50 * j, dur=50.0,
                           args=dict(hlo_category="convolution fusion")))
        for j in range(extra_fusion):                       # e.g. a one-time copy fusion outside the loop
            ev.append(dict(ph="X", pid=3, tid=2, name="fusion.7", ts=1000 * i + 900, dur=30.0,
                           args=dict(hlo_category="data formatting")))
        ev.append(dict(ph="X", pid=9, tid=5, name="convolution on the host", ts=1000 * i, dur=7.0))
    with gzip.open(os.path.join(d, "plugins", "profile", "run", "perfetto_trace.json.gz"), "wt") as f:
        json.dump(dict(traceEvents=ev), f)
    return d


m.ARGS = argparse.Namespace(device_match=r"TPU:0\b")
good = m.op_events(trace(8), 8)
assert good["verified"] and good["modules"] == 2 and good["conv_events"] == 16 and good["conv_us"] == 800.0, good
assert good["ops"]["while.4"] == [2, 900.0] and len(good["files"]["perfetto_trace.json.gz"]) == 64
nine = m.op_events(trace(9), 8)
assert not nine["verified"] and nine["conv_events"] == 18
hoisted = m.op_events(trace(1), 8)                          # the GEMM ran once per module, not once per iteration
assert not hoisted["verified"]
mixed = m.op_events(trace(8, extra_fusion=1), 8)            # another fusion does not count as the convolution
assert mixed["verified"] and mixed["conv_events"] == 16 and mixed["ops"]["fusion.7"] == [2, 60.0]
assert m.op_events(tempfile.mkdtemp(), 8) is None           # no trace

# pair_verdict(): never more than the artifacts show
same = copy.deepcopy(loop)
v, notes = m.pair_verdict(loop, same, 4096, good, good)
assert "weight staging NOT ESTABLISHED" in v and not notes and "not a measured DMA tail" in v, v
copies = copy.deepcopy(loop); copies["copy_start"] += 1    # an extra copy alone is not weight staging
v, notes = m.pair_verdict(loop, copies, 4096, good, good)
assert "NOT ESTABLISHED" in v and "raw evidence" in v and "which operand is copied is not identified" in notes[0], (v, notes)
tagged = copy.deepcopy(loop); tagged["tagged"] = loop["tagged"] + ["bf16[1,4096,4096]:S(1)"]
v, notes = m.pair_verdict(loop, tagged, 4096, good, good)  # a tag is raw evidence, not a traced dataflow
assert "NOT ESTABLISHED" in v and "raw evidence" in v and "S(1)" in notes[0], (v, notes)
for bad in (None, nine, hoisted):                           # loop membership not verified in the trace
    v, notes = m.pair_verdict(loop, same, 4096, good, bad)
    assert v.startswith("inconclusive") and "not verified in both traces" in v, v
other = copy.deepcopy(loop); other["windows"][0]["iteration_bounds"] = "8,1,1"
assert m.pair_verdict(loop, other, 4096, good, good)[0].startswith("inconclusive")
more = copy.deepcopy(loop); more["fusions"] += 1
assert m.pair_verdict(loop, more, 4096, good, good)[0].startswith("inconclusive")
assert m.pair_verdict(m.hlo_facts(""), m.hlo_facts(""), 4096, good, good)[0].startswith("inconclusive")
assert m.pair_verdict(loop, iso, 4096, good, good)[0].startswith("inconclusive")
for v in (m.pair_verdict(loop, x, 4096, good, good)[0] for x in (same, copies, tagged)):
    assert "staged per" not in v and "reused weight staged" not in v and "NOT staged" not in v, v

# the experiment on the CPU backend: runs, keeps raw wall samples and a signed slope, concludes nothing
out = tempfile.mkdtemp()
script = os.path.join(HERE, "measure_tpu_v5e1.py")
p = subprocess.run([sys.executable, script, "--only", "tailprobe", "--quick", "--max-iters", "5", "--rotate-mb", "2",
                    "--out", out], capture_output=True, text=True)
assert p.returncode == 0, p.stderr[-800:]
rows = list(csv.DictReader(open(os.path.join(out, "tailprobe.csv"))))
assert [(r["dtype"], r["mode"]) for r in rows] == [("bf16", "banked"), ("bf16", "single"), ("int8", "banked"), ("int8", "single")]
for r in rows:
    wn, w2n = json.loads(r["wall_n_s"]), json.loads(r["wall_2n_s"])
    assert len(wn) == len(w2n) == 5 and all(math.isfinite(x) and x > 0 for x in wn + w2n), r
    assert math.isfinite(float(r["loop_slope_us"])) and r["slope_flag"] in ("ok", "nonpositive: wall slope unusable"), r
    assert (float(r["loop_slope_us"]) > 0) == (r["slope_flag"] == "ok")
    assert r["conv_fusion_us_per_event"] == "" and (r["weights_in_loop"] == "1") == (r["mode"] == "single"), r
    assert len(r["hlo_sha256"]) == 64
pairs = list(csv.DictReader(open(os.path.join(out, "tailprobe_pairs.csv"))))
assert len(pairs) == 2 and all(q["verdict"].startswith("inconclusive") for q in pairs), pairs
meta = json.load(open(os.path.join(out, "tailprobe_meta.json")))
assert meta["backend"] == "cpu" and meta["quick"] and meta["device_kind"] and meta["profiled_loop_iterations"] == 8
assert "not_measured" in meta and "smoke check" in meta["cpu_note"] and len(meta["script_sha256"]) == 64
assert len(os.listdir(os.path.join(out, "hlo"))) == 4
p = subprocess.run([sys.executable, script, "--only", "tailprobe", "--out", out], capture_output=True, text=True)
assert p.returncode != 0 and "no TPU backend" in (p.stderr + p.stdout), p.stderr[-300:]
print("ok (functional checks on CPU only; no hardware validation)")
