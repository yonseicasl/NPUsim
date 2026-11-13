#!/usr/bin/env python3
"""[accelerator] fused_output_consumer is refused, before any layer runs, for what it cannot describe -- and
accepted for a one-layer GEMM, where it removes exactly the output's DRAM traffic (needs models/model; ~10 s)."""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "palm"))
import palm_decode as palm

palm.MXUS, palm.CLOCK_HZ = 1, 1.5e9
tmp, cfgroot, rundir = palm.workdir(("tpu_v5e",))
base = open(os.path.join(cfgroot, "accelerators", "tpu_v5e.cfg")).read()
assert "fused_output_consumer" not in base                     # no config uses it by default


def run(accel, cfg_text, name, two_layers=False):
    open(os.path.join(cfgroot, "accelerators", accel + ".cfg"), "w").write(cfg_text)
    palm.gemm_files(cfgroot, rundir, name, accel, 8, 128, 128)
    if two_layers:                                             # a second mapped layer in network and mapping
        net = os.path.join(cfgroot, "networks", name + ".cfg")
        s = open(net).read()
        open(net, "w").write(s.replace("[softmax]", "[connected]\noutput=128\nactivation=linear\n\n[softmax]"))
        mp = os.path.join(cfgroot, "mappings", accel, name, "ws.map")
        body = open(mp).read()
        open(mp, "w").write(body + "\n" + body[body.index("[connected]"):])
        with open(os.path.join(rundir, "weights", name + ".wgh"), "wb") as f:
            f.truncate(2 * 4 * (128 + 128 * 128))
    p = subprocess.run([os.path.join(palm.ROOT, "models", "model"), "run", accel, name, "ws"], cwd=rundir,
                       env=dict(os.environ, NPUSIM_CONFIG_ROOT=cfgroot), capture_output=True, text=True)
    report = os.path.join(rundir, f"{accel}_{name}_layer_0.txt")
    return p, (open(report).read() if os.path.exists(report) else "")


def dram_tx(report):
    blk = report[report.index("Multi-chip <-> DRAM transactions"):][:400]
    return tuple(int(re.search(rf"\* {k}\s*:\s*(\d+)", blk).group(1)) for k in ("Input data", "Weight", "Output data"))


off = re.sub(r"(?m)^activation_residency\s*=.*$", "activation_residency = 0", base)
on = off.replace("activation_residency = 0", "activation_residency = 0\nfused_output_consumer = 1")

p0, r0 = run("g_off", off, "one_off")
p1, r1 = run("g_on", on, "one_on")
assert p0.returncode == 0 and p1.returncode == 0, (p0.stderr[-300:], p1.stderr[-300:])
t0, t1 = dram_tx(r0), dram_tx(r1)
assert t0[2] > 0 and t1 == (t0[0], t0[1], 0), (t0, t1)         # input and weight kept, output removed
assert "fused_output_consumer: input DRAM-backed, output consumed on chip" in r1 and "fused_output_consumer" not in r0

p, r = run("g_conflict", on.replace("activation_residency = 0", "activation_residency = 1"), "conflict")
assert p.returncode != 0 and "cannot be combined with activation_residency" in p.stderr and not r, p.stderr[-300:]

p, r = run("g_two", on, "two", two_layers=True)
assert p.returncode != 0 and "exactly one mapped layer" in p.stderr and "2 mapped layers" in p.stderr, p.stderr[-400:]
assert not r, "a layer ran before the refusal"

p, r = run("g_two_off", off, "two_off", two_layers=True)       # the same two-layer run is fine without the flag
assert p.returncode == 0, p.stderr[-400:]
shutil.rmtree(tmp)
print("ok")
