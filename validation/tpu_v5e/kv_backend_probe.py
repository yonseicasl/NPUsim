#!/usr/bin/env python3
"""Which timing backend does NPUsim's attention consumer ([kvcache] attention = 1) use for its KV read?

Each (B, L) case runs twice on identical inputs: without and with $NPUSIM_MOSAIC_DRAM (a MOSAIC Ramulator2 yaml).
Everything is kept under OUT/<case>_<backend>/ (accelerator cfg with the [kvcache] section, network, mapping,
both reports, stdout). From the reports: the DRAM timing-model line, the tile-stream DRAM busy of the small
projection the step follows (this MUST change with the yaml, or the yaml was not in effect), and the attention
step with its K / V supply. If the K / V supply is identical in both runs, the KV read does not go through the
yaml's DRAM model. No measurement is involved.

Usage: python3 kv_backend_probe.py DRAM_YAML OUT_DIR [B:L ...]      (needs models/model)
"""
import csv
import hashlib
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_attention as ca
palm = ca.palm
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()


def run(B, L, yaml, out):
    tag = f"attn_b{B}_l{L}_{'mosaic' if yaml else 'analytic'}"
    d = os.path.join(out, tag)
    if os.path.exists(d):
        raise SystemExit(f"{d} exists; nothing was changed -- use a fresh OUT_DIR")
    cfgroot, rundir = os.path.join(d, "configs"), os.path.join(d, "run")
    for sub in ("accelerators", "networks", "mappings"):
        os.makedirs(os.path.join(cfgroot, sub))
    os.makedirs(os.path.join(rundir, "weights"))
    os.symlink(os.path.join(ca.ROOT, "models", "datasets"), os.path.join(rundir, "datasets"))
    cfg = open(os.path.join(ca.ROOT, "configs", "accelerators", f"{ca.ACCEL}.cfg")).read()
    open(os.path.join(cfgroot, "accelerators", f"{ca.ACCEL}.cfg"), "w").write(
        cfg + ca.KVCACHE.format(qh=ca.QH, kvh=ca.KVH, dh=ca.DH, kvb=ca.KV_BYTES, L=L, B=B))
    name = f"attn_b{B}_l{L}"
    palm.gemm_files(cfgroot, rundir, name, ca.ACCEL, B, 128, 512)
    env = dict(os.environ, NPUSIM_CONFIG_ROOT=cfgroot)
    env.pop("NPUSIM_MOSAIC_DRAM", None)
    if yaml:
        env["NPUSIM_MOSAIC_DRAM"] = yaml
    p = subprocess.run([os.path.join(ca.ROOT, "models", "model"), "run", ca.ACCEL, name, "ws"], cwd=rundir, env=env,
                       capture_output=True, text=True)
    open(os.path.join(d, "stdout.txt"), "w").write(p.stdout + p.stderr)
    shutil.rmtree(os.path.join(rundir, "weights"))            # zero-filled placeholder, not evidence
    if p.returncode:
        raise SystemExit(f"{tag} failed: {(p.stdout + p.stderr)[-800:]}")
    net, layer = (open(os.path.join(rundir, f"{ca.ACCEL}_{name}{s}.txt")).read() for s in ("", "_layer_0"))
    both = net + layer
    g = lambda rx: (re.search(rx, both) or [None, None])[1]
    return dict(B=B, L=L, backend_env=os.path.basename(yaml) if yaml else "(unset)",
                kv_read_bytes=B * L * ca.KVH * ca.DH * ca.KV_BYTES * 2,
                dram_timing_model=g(r"DRAM timing model\s*:\s*(.*)"),
                projection_dram_busy_cycles=g(r"\* DRAM\s*:\s*([\d.]+) cycles \("),
                attention_step_cycles=g(r"Attention step\s*:\s*([\d.]+)"),
                k_supply_cycles=g(r"K / V supply\s*:\s*([\d.]+) /"), v_supply_cycles=g(r"K / V supply\s*:\s*[\d.]+ / ([\d.]+)"),
                kv_supply_dram_cycles=g(r"KV supply DRAM cycles\s*:\s*([\d.]+)"),
                kv_schedule=g(r"KV schedule\s*:\s*(.*)"),
                report_sha256=sha(os.path.join(rundir, f"{ca.ACCEL}_{name}_layer_0.txt")), run_dir=tag)


def main():
    yaml, out = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    cases = [tuple(int(x) for x in c.split(":")) for c in sys.argv[3:]] or [(1, 2048), (8, 2048)]
    os.makedirs(out, exist_ok=True)
    rows = [run(B, L, y, out) for B, L in cases for y in (None, yaml)]
    with open(os.path.join(out, "kv_backend_probe.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for a, b in zip(rows[::2], rows[1::2]):
        same = lambda k: a[k] == b[k]
        print(f"B {a['B']} L {a['L']} ({a['kv_read_bytes']} KV bytes read):")
        print(f"   DRAM timing model   unset: {a['dram_timing_model']}")
        print(f"                       yaml : {b['dram_timing_model']}")
        print(f"   projection DRAM busy {a['projection_dram_busy_cycles']} -> {b['projection_dram_busy_cycles']} cycles"
              f"  ({'UNCHANGED: yaml not in effect?' if same('projection_dram_busy_cycles') else 'changes with the yaml'})")
        for k in ("k_supply_cycles", "v_supply_cycles", "kv_supply_dram_cycles", "attention_step_cycles"):
            print(f"   {k:22s} {a[k]} -> {b[k]}  ({'IDENTICAL' if same(k) else 'differs'})")
        rate = a["kv_read_bytes"] / (float(a["kv_supply_dram_cycles"]) / ca.CLOCK_HZ) / 1e9
        print(f"   KV bytes / KV supply DRAM time = {rate:.1f} GB/s in both runs")


if __name__ == "__main__":
    main()
