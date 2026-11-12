#!/usr/bin/env python3
"""parse_trace on a synthetic TPU-shaped perfetto trace, linfit, and that measure_tpu_v5e8.py's shared code
matches measure_tpu_v5e1.py's (no TPU needed; needs jax importable)."""
import gzip, json, os, sys, tempfile, types
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import measure_tpu_v5e1 as m

m.ARGS = types.SimpleNamespace(device_match=r"TPU:0\b")
M = lambda pid, tid, kind, name: {"ph": "M", "pid": pid, "tid": tid, "name": kind, "args": {"name": name}}
X = lambda pid, tid, name, dur, **a: {"ph": "X", "pid": pid, "tid": tid, "name": name, "dur": dur, "ts": 0, "args": a}
ev = [M(1, None, "process_name", "/host:CPU"), M(3, None, "process_name", "/device:TPU:0"),
      M(4, None, "process_name", "/device:TPU:1"), M(3, 1, "thread_name", "XLA Modules"),
      M(3, 2, "thread_name", "XLA Ops"), M(4, 1, "thread_name", "XLA Modules"),
      X(3, 1, "jit_layer", 100.0), X(3, 1, "jit_layer", 120.0), X(4, 1, "jit_layer", 999.0),
      X(3, 2, "fusion.1", 30.0, tf_op="jit(layer)/qkv/dot_general"),
      X(3, 2, "fusion.1", 30.0, tf_op="jit(layer)/qkv/dot_general"),
      X(3, 2, "all-reduce", 10.0, tf_op="jit(layer)/psum_attn/psum"), X(3, 2, "copy.2", 4.0)]
d = tempfile.mkdtemp()
os.makedirs(os.path.join(d, "plugins", "profile", "run"))
with gzip.open(os.path.join(d, "plugins", "profile", "run", "perfetto_trace.json.gz"), "wt") as f:
    json.dump({"traceEvents": ev}, f)
mods, ops = m.parse_trace(d)
assert mods == [100.0, 120.0], mods
assert ops == {"qkv|fusion.1": 60.0, "psum_attn|all-reduce": 10.0, "-|copy.2": 4.0}, ops
assert m.parse_trace(tempfile.mkdtemp()) is None
a, b = m.linfit([1, 2, 3], [3, 5, 7])
assert abs(a - 1) < 1e-9 and abs(b - 2) < 1e-9, (a, b)

# measure_tpu_v5e8.py is self-contained: its measurement machinery must stay identical to measure_tpu_v5e1.py's
import inspect
import measure_tpu_v5e8 as m8
for name in ("parse_trace", "bench", "write", "rnd", "nsets", "fmt", "rmsnorm", "linfit", "run",
             "gemm", "vector", "decode", "smallk", "loopk"):
    assert inspect.getsource(getattr(m, name)) == inspect.getsource(getattr(m8, name)), f"{name} differs"
assert m.SCOPES == m8.SCOPES
print("ok")
