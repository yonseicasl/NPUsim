#!/usr/bin/env python3
"""measure_tpu_v5e8.py as a notebook would run it (CPU, 8 virtual devices; functional checks only):
the kernel's "-f <connection file>" argument is ignored, other unknown arguments are refused, a session without a
TPU gets an instruction instead of a JAX error, and a failing experiment or case does not lose the others.
Run with a Python that has jax."""
import csv
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ENV = dict(os.environ, XLA_FLAGS="--xla_force_host_platform_device_count=8")


def cell(code, devices=8):
    """Run code as a pasted cell would: the script's namespace with __name__ != "__main__", kernel-style argv."""
    src = f"""
import sys
sys.argv = ["ipykernel_launcher.py", "-f", "/root/.local/share/jupyter/runtime/kernel-abc.json"]
sys.path.insert(0, {HERE!r})
import measure_tpu_v5e8 as m
{code}
"""
    env = dict(ENV, XLA_FLAGS=f"--xla_force_host_platform_device_count={devices}")
    return subprocess.run([sys.executable, "-c", src], capture_output=True, text=True, env=env)


out = tempfile.mkdtemp()
# pasted cell, default options, no TPU: an instruction, not a traceback from inside JAX
p = cell(f"m.main(['--out', {out!r}])")
assert "not on a TPU" in (p.stdout + p.stderr) and "Accelerator" in (p.stdout + p.stderr), p.stderr[-600:]
assert "unrecognized arguments" not in p.stderr, p.stderr[-300:]
err = json.load(open(os.path.join(out, "errors.json"))) if os.path.exists(os.path.join(out, "errors.json")) else {}
assert not err, err                                           # refused up front, not recorded as a crash

# the kernel's -f is ignored when main() takes sys.argv itself; anything else unknown is refused
p = cell("import sys; sys.argv += ['--quick', '--max-iters', '3', '--rotate-mb', '1', '--only', 'collectives', '--out', %r]; m.main()" % out)
assert p.returncode == 0 and "wrote collectives.csv" in p.stdout, p.stderr[-800:]
rows = list(csv.DictReader(open(os.path.join(out, "collectives.csv"))))
assert len(rows) == 8 and {r["op"] for r in rows} == {"psum", "all_gather", "psum_scatter", "all_to_all"}
assert all(r["chips"] == "8" and float(r["wall_us"]) > 0 for r in rows)
p = cell("import sys; sys.argv += ['--quick', '--bogus']; m.main()")
assert p.returncode != 0 and "unrecognized arguments: --bogus" in p.stderr, p.stderr[-300:]

# one chip only: said so
p = cell(f"m.main(['--quick', '--only', 'collectives', '--out', {out!r}])", devices=1)
assert "needs >= 2 chips" in (p.stdout + p.stderr), (p.stdout + p.stderr)[-400:]

# a failing case is recorded and the sweep continues; a failing experiment does not stop the next one
out2 = tempfile.mkdtemp()
p = cell(f"""
real = m.jax.lax.all_gather
def broken(*a, **k): raise RuntimeError("injected failure")
m.jax.lax.all_gather = broken
def boom(): raise ValueError("injected experiment failure")
m.run({{"boom": boom, "collectives": m.collectives}}, "boom,collectives", "tpu_v5e8",
      ['--quick', '--max-iters', '3', '--rotate-mb', '1', '--out', {out2!r}])
""")
assert p.returncode == 0, p.stderr[-800:]
assert "== boom: FAILED" in p.stdout and "FAILED: boom (errors.json)" in p.stdout, p.stdout[-600:]
assert "injected experiment failure" in json.load(open(os.path.join(out2, "errors.json")))["boom"]
rows = list(csv.DictReader(open(os.path.join(out2, "collectives.csv"))))
bad = [r for r in rows if r["op"] == "all_gather"]
assert len(rows) == 8 and len(bad) == 2 and all("injected failure" in r["error"] for r in bad), rows
assert all(float(r["wall_us"]) > 0 for r in rows if r["op"] != "all_gather")
print("ok (functional checks on CPU only)")
