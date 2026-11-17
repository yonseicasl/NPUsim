#!/usr/bin/env python3
"""chainprobe (measure_tpu_v5e1.py) on CPU: functional only -- rows, kept HLO, refusal to touch earlier output.
Usage: <python with jax> test_chainprobe.py"""
import csv
import glob
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
out = tempfile.mkdtemp()
cmd = [sys.executable, os.path.join(HERE, "measure_tpu_v5e1.py"), "--only", "chainprobe", "--quick", "--out", out]
p = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, JAX_PLATFORMS="cpu"))
assert p.returncode == 0, p.stdout[-2000:] + p.stderr[-2000:]
rows = list(csv.DictReader(open(os.path.join(out, "chainprobe.csv"))))
assert [(r["kind"], r["n"]) for r in rows] == [("args", "1"), ("stacked", "1"), ("args", "2"), ("stacked", "2")], rows
for r in rows:
    text = open(os.path.join(out, "hlo", f"chainprobe_{r['kind']}_{r['M']}x{r['K']}_n{r['n']}.txt")).read()
    assert text and json.loads(r["hlo_facts"]) is not None and len(r["hlo_sha256"]) == 64 and float(r["wall_us"]) > 0
    assert r["weight_bytes"] == str(2 * int(r["n"]) * int(r["K"]) ** 2) and r["backend"] == "cpu"
before = {f: os.path.getmtime(f) for f in glob.glob(os.path.join(out, "**", "*"), recursive=True)}
p = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, JAX_PLATFORMS="cpu"))
assert p.returncode != 0 and "already holds chainprobe output" in p.stderr + p.stdout, p.stderr[-500:]
assert before == {f: os.path.getmtime(f) for f in glob.glob(os.path.join(out, "**", "*"), recursive=True)}
# without --quick a non-TPU backend is refused before anything is measured
out2 = tempfile.mkdtemp()
p = subprocess.run(cmd[:4] + ["--out", out2], capture_output=True, text=True, env=dict(os.environ, JAX_PLATFORMS="cpu"))
assert "chainprobe: no TPU" in p.stdout + p.stderr and not os.path.exists(os.path.join(out2, "chainprobe.csv")), p.stdout[-800:]
print("chainprobe checks passed")
