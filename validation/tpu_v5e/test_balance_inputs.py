#!/usr/bin/env python3
"""balance_npusim.py gives every run its own, correct inputs (no simulator needed): several shapes x A/B/C are
generated, then every saved network and mapping must match the shape in its run's name; and check_inputs()
rejects a run that reads another shape's network -- the shared-name bug of the first batches."""
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import balance_npusim as b

run = os.path.join(HERE, "results", "2026-11-12_colab_v5e1_balance")
out = tempfile.mkdtemp()
only = "bf16:32:4096,bf16:512:8192,int8:128:4096,int8:64:8192"
subprocess.run([sys.executable, os.path.join(HERE, "balance_npusim.py"), run, out, "--only", only, "--generate-only"],
               check=True, capture_output=True)
gen = os.path.join(out, "generated_analytic")
runs = sorted(os.listdir(gen))
assert len(runs) == 12, runs                                   # 4 shapes x A, B, C
for accel in runs:
    dt, M, K, v = re.fullmatch(r"bal_(bf16|int8)_(\d+)x(\d+)_([abc])", accel).groups()
    M, K = int(M), int(K)
    net = open(os.path.join(gen, accel, "network.cfg")).read()
    g = lambda key: int(re.search(rf"(?m)^{key}=(\d+)", net).group(1))
    assert (g("batch"), g("height") * g("width"), g("output")) == (M, K, K), (accel, net)
    assert f"# {K}, {M}, 1, 1, {K}," in open(os.path.join(gen, accel, "mappings", accel, "ws.map")).read(), accel
    cfg = open(os.path.join(gen, accel, accel + ".cfg")).read()
    assert ("weight_format = int8" in cfg) == (dt == "int8"), accel
    assert ("offchip_window_count" in cfg) == (v != "a") and ("window_switch_cycle" in cfg) == (v == "c"), accel
    assert "activation_residency = 2" in cfg and "program_setup_cycle = 0" in cfg, accel
    assert "fused_output_consumer" not in cfg, accel
shutil.rmtree(out)

# --activations fused_output: its own artifact names and config lines; the default mode is untouched above
out = tempfile.mkdtemp()
subprocess.run([sys.executable, os.path.join(HERE, "balance_npusim.py"), run, out, "--only", "int8:128:4096",
                "--generate-only", "--activations", "fused_output"], check=True, capture_output=True)
gen = os.path.join(out, "generated_analytic_fusedout")
for accel in sorted(os.listdir(gen)):
    cfg = open(os.path.join(gen, accel, accel + ".cfg")).read()
    assert "fused_output_consumer = 1" in cfg and "activation_residency = 0" in cfg, accel
    assert "activation_residency = 2" not in cfg, accel
shutil.rmtree(out)

# the old bug, reproduced: two shapes written under one name -> the first run's network is the second shape
tmp, cfgroot, rundir = b.palm.workdir(("tpu_v5e",))
for M, K in ((32, 4096), (512, 8192)):
    b.palm.gemm_files(cfgroot, rundir, "g", f"acc{M}", M, K, K)
b.check_inputs(cfgroot, "g", "acc512", 512, 8192, 8192)
try:
    b.check_inputs(cfgroot, "g", "acc32", 32, 4096, 4096)
    raise SystemExit("check_inputs accepted a run whose network is another shape")
except AssertionError as e:
    assert "network is (512, 8192, 8192)" in str(e), e
shutil.rmtree(tmp)
print("ok")
