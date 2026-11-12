#!/usr/bin/env python3
"""measure_tpu_v5e1.py tailprobe: FUNCTIONAL checks only (CPU; no TPU, no timing claim).
hlo_facts() on HLO texts saved from a real TPU run, pair_verdict() on constructed facts, the --quick run of the
experiment on the CPU backend (files written, verdict present), and the TPU guard for full sizes.
Run with a Python that has jax."""
import copy
import csv
import json
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
assert not [t for t in loop["tagged"] if "4096,4096" in t], loop["tagged"]     # the weights: default memory space

iso = m.hlo_facts(open(os.path.join(H, "balance_bf16_256x4096_isolated.txt")).read())
assert iso["copy_start"] == 1 and iso["copy_done"] == 1 and iso["prefetch"] and iso["whiles"] == 0, iso
assert "bf16[4096,4096]:S(1)" in iso["tagged"], iso["tagged"]                  # the staged copy of the weights

i8 = m.hlo_facts(open(os.path.join(H, "balance_int8_128x4096_loop.txt")).read())
assert i8["windows"][0]["kernel_window_bounds"] == "512,3" and i8["convolutions"][0][1][1] == "s8[4096,4096]", i8

# pair_verdict: only what the two HLO texts show
v, why = m.pair_verdict(loop, copy.deepcopy(loop), 4096)
assert v.startswith("comparable schedule, reused weight NOT staged") and not why, (v, why)
staged = copy.deepcopy(loop); staged["tagged"] = loop["tagged"] + ["bf16[1,4096,4096]:S(1)"]
v, why = m.pair_verdict(loop, staged, 4096)
assert v.startswith("comparable schedule, reused weight staged per the HLO") and "not a measured DMA tail" in v, v
other = copy.deepcopy(loop); other["windows"][0]["iteration_bounds"] = "8,1,1"
v, why = m.pair_verdict(loop, other, 4096)
assert v.startswith("inconclusive") and "windows" in why[0], (v, why)
hoisted = copy.deepcopy(loop); hoisted["fusions"] += 1
assert m.pair_verdict(loop, hoisted, 4096)[0].startswith("inconclusive")
empty = m.hlo_facts("")
assert m.pair_verdict(empty, empty, 4096)[0].startswith("inconclusive")        # e.g. a CPU compile: nothing to compare
v, why = m.pair_verdict(loop, iso, 4096)                                       # a one-op program is another schedule
assert v.startswith("inconclusive"), v

# the experiment itself on CPU: runs, writes its files, labels the pair inconclusive (no TPU window_config)
out = tempfile.mkdtemp()
script = os.path.join(HERE, "measure_tpu_v5e1.py")
p = subprocess.run([sys.executable, script, "--only", "tailprobe", "--quick", "--max-iters", "5", "--rotate-mb", "2",
                    "--out", out], capture_output=True, text=True)
assert p.returncode == 0, p.stderr[-800:]
rows = list(csv.DictReader(open(os.path.join(out, "tailprobe.csv"))))
assert [(r["dtype"], r["mode"]) for r in rows] == [("bf16", "banked"), ("bf16", "single"), ("int8", "banked"), ("int8", "single")]
assert all(r["weights_in_loop"] == ("1" if r["mode"] == "single" else r["weights_in_loop"]) and float(r["loop_slope_us"]) >= 0 for r in rows)
assert all(int(r["weights_in_loop"]) >= 2 for r in rows if r["mode"] == "banked")
pairs = list(csv.DictReader(open(os.path.join(out, "tailprobe_pairs.csv"))))
assert len(pairs) == 2 and all(q["verdict"].startswith("inconclusive") for q in pairs), pairs
meta = json.load(open(os.path.join(out, "tailprobe_meta.json")))
assert meta["backend"] == "cpu" and meta["quick"] and "not_measured" in meta and len(meta["script_sha256"]) == 64
assert len(os.listdir(os.path.join(out, "hlo"))) == 4

p = subprocess.run([sys.executable, script, "--only", "tailprobe", "--out", out], capture_output=True, text=True)
assert p.returncode != 0 and "no TPU backend" in (p.stderr + p.stdout), p.stderr[-300:]
print("ok (functional checks on CPU only; no hardware validation)")
