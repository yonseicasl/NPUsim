#!/usr/bin/env python3
"""parse_trace on a synthetic TPU-shaped perfetto trace (no TPU needed; needs jax importable)."""
import gzip, json, os, sys, tempfile, types
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import measure_tpu as m

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
print("ok")
