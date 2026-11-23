#!/usr/bin/env python3
"""Which backend times the KV read of NPUsim's attention consumer ([kvcache] attention = 1)? Small gate, no TPU data.

Each (B, L) case runs twice on identical inputs: without and with $NPUSIM_MOSAIC_DRAM (a MOSAIC Ramulator2 yaml);
the first case also with kv_schedule = blocking. OUT/<run>/ keeps the accelerator cfg with its [kvcache] section,
network, mapping, both reports and stdout; kv_backend_probe.csv / .json hold what the reports say and the checks:
  inputs      binary, libraries, yaml, base config and scripts hashed before and after the runs
  bytes       reported KV read and append bytes against the geometry (B x L x 2 x kv heads x head dim x bytes)
  no backend  with --baseline (an earlier probe directory): the runs without a yaml reproduce it -- K / V supply,
              attention step, projection DRAM busy, and the layer report byte for byte
  backend     K / V supply differ from the no-yaml run only by the DRAM read time; with --diag (the csv of
              unittest/run_mosaic_kv_stream_test.sh, same binary) the reported K / V DRAM read equals the stream
              helper's time for the same bytes, exact up to 32 MiB and by the stated long-stream rule beyond;
              supply = DRAM read + decoder share (nothing added twice); QK / softmax / AV unchanged
Scope of every number: the attention step appended after a small dummy projection; the K / V stream time is a
fresh-stream TOTAL per pass that the streaming schedule spreads uniformly over its tiles -- not per-tile arrivals.
QK / AV and softmax rates are configured assumptions (compare_attention.py), not hardware-validated.

Usage: python3 kv_backend_probe.py DRAM_YAML OUT_DIR [B:L ...] [--baseline DIR] [--diag CSV]   (needs models/model)
"""
import argparse
import csv
import hashlib
import json
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
MIB = 1 << 20


def frozen(yaml):
    """sha256 of everything a run executes or reads."""
    model = os.path.join(ca.ROOT, "models", "model")
    # Every shared object the binary loads, resolved under THIS environment (the runs inherit it): `name => path`
    # lines and the ones listed by absolute path (libnpusim, libnebula). An unresolved dependency is an error.
    ldd = subprocess.run(["ldd", model], capture_output=True, text=True).stdout
    if "not found" in ldd:
        raise SystemExit("unresolved libraries of models/model under this LD_LIBRARY_PATH:\n" + ldd)
    libs = sorted(set(re.findall(r"(?m)^\s*(?:\S+ => )?(/\S+) \(0x", ldd)))
    assert any("libnpusim" in l for l in libs) and any("libramulator" in l for l in libs), libs
    files = [model, yaml, os.path.join(ca.ROOT, "configs", "accelerators", f"{ca.ACCEL}.cfg"), os.path.abspath(__file__),
             os.path.join(HERE, "compare_attention.py"), os.path.join(HERE, "..", "palm", "palm_decode.py")] + libs
    return {os.path.relpath(f, ca.ROOT) if f.startswith(ca.ROOT) else f: sha(f) for f in files}


def run(B, L, yaml, out, schedule="streaming"):
    tag = f"attn_b{B}_l{L}_{'mosaic' if yaml else 'analytic'}" + ("" if schedule == "streaming" else f"_{schedule}")
    d = os.path.join(out, tag)
    if os.path.exists(d):
        raise SystemExit(f"{d} exists; nothing was changed -- use a fresh OUT_DIR")
    cfgroot, rundir = os.path.join(d, "configs"), os.path.join(d, "run")
    for sub in ("accelerators", "networks", "mappings"):
        os.makedirs(os.path.join(cfgroot, sub))
    os.makedirs(os.path.join(rundir, "weights"))
    os.symlink(os.path.join(ca.ROOT, "models", "datasets"), os.path.join(rundir, "datasets"))
    cfg = open(os.path.join(ca.ROOT, "configs", "accelerators", f"{ca.ACCEL}.cfg")).read()
    kv = ca.KVCACHE.format(qh=ca.QH, kvh=ca.KVH, dh=ca.DH, kvb=ca.KV_BYTES, L=L, B=B)
    assert kv.count("kv_schedule = streaming") == 1
    kv = kv.replace("kv_schedule = streaming", f"kv_schedule = {schedule}")
    open(os.path.join(cfgroot, "accelerators", f"{ca.ACCEL}.cfg"), "w").write(cfg + kv)
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
    f = lambda rx: (lambda v: None if v is None else float(v))(g(rx))
    read_dram = re.search(r"K / V DRAM read\s*:\s*([\d.]+) / ([\d.]+) cycles, decoder excluded \((\d+) B per pass; "
                          r"K (\w+), V (\w+)\)", layer)
    dense, comp = (int(x) for x in re.search(r"KV bytes dense/comp\s*:\s*(\d+) / (\d+) B", layer).groups())
    occupancy = int(g(r"cache occupancy\s*:\s*(\d+) B after append"))
    qk = re.search(r"QK / softmax / AV\s*:\s*([\d.]+) / ([\d.]+) / ([\d.]+)", layer).groups()
    return dict(B=B, L=L, schedule=schedule, backend_env=os.path.basename(yaml) if yaml else "(unset)",
                dram_timing_model=g(r"DRAM timing model\s*:\s*(.*)"),
                kv_read_bytes_reported=comp, kv_dense_bytes_reported=dense, kv_append_bytes_reported=occupancy - dense,
                kv_read_bytes_geometry=B * L * 2 * ca.KVH * ca.DH * ca.KV_BYTES,
                kv_append_bytes_geometry=B * 2 * ca.KVH * ca.DH * ca.KV_BYTES,
                projection_dram_busy_cycles=f(r"\* DRAM\s*:\s*([\d.]+) cycles \("),
                attention_step_cycles=f(r"Attention step\s*:\s*([\d.]+)"),
                k_supply_cycles=f(r"K / V supply\s*:\s*([\d.]+) /"), v_supply_cycles=f(r"K / V supply\s*:\s*[\d.]+ / ([\d.]+)"),
                kv_supply_dram_cycles=f(r"KV supply DRAM cycles\s*:\s*([\d.]+)"),
                decoder_cycles=f(r"decoder reconstitute\s*:\s*([\d.]+)"),
                k_dram_read_cycles=read_dram and float(read_dram.group(1)),
                v_dram_read_cycles=read_dram and float(read_dram.group(2)),
                read_stream_bytes=read_dram and int(read_dram.group(3)),
                k_stream=read_dram and read_dram.group(4), v_stream=read_dram and read_dram.group(5),
                append_cycles=f(r"cache append write\s*:\s*([\d.]+)"),
                append_label="analytical" if re.search(r"cache append write\s*:.*\(analytical\)", layer) else "",
                qk_softmax_av_cycles="/".join(qk), report_sha256=sha(os.path.join(rundir, f"{ca.ACCEL}_{name}_layer_0.txt")),
                run_dir=tag)


def helper_cycles(diag, nbytes, freq_mhz):
    """([dram]-clock cycles, 'exact' | 'extrapolated') of one read stream by the stream helper's own samples."""
    exact, t16, t32 = {}, None, None
    for r in csv.reader(open(diag)):
        if r[0] == "two_streams":
            exact[int(r[1])], exact[2 * int(r[1])] = float(r[2]), float(r[3])
        if r[0] == "extrapolation":
            t16, t32 = float(r[1]), float(r[2])
    if nbytes <= 32 * MIB:
        return (exact[nbytes] * freq_mhz / 1e3, "exact") if nbytes in exact else (None, "exact")
    return (t32 + (nbytes - 32 * MIB) * (t32 - t16) / (16 * MIB)) * freq_mhz / 1e3, "extrapolated"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("yaml")
    ap.add_argument("out")
    ap.add_argument("cases", nargs="*", default=["1:2048", "8:2048", "32:8192"])
    ap.add_argument("--baseline", help="an earlier probe directory (kv_backend_probe.csv) the no-yaml runs must reproduce")
    ap.add_argument("--diag", help="mosaic_kv_stream_diag.csv of unittest/run_mosaic_kv_stream_test.sh (same binary)")
    a = ap.parse_args()
    yaml, out = os.path.abspath(a.yaml), os.path.abspath(a.out)
    cases = [tuple(int(x) for x in c.split(":")) for c in a.cases]
    os.makedirs(out, exist_ok=True)
    if os.path.exists(os.path.join(out, "kv_backend_probe.csv")):
        raise SystemExit(f"{out} already holds a probe; nothing was changed -- use a fresh OUT_DIR")
    before = frozen(yaml)
    rows = [run(B, L, y, out) for B, L in cases for y in (None, yaml)]
    rows += [run(*cases[0], y, out, "blocking") for y in (None, yaml)]
    after = frozen(yaml)
    cfg = open(os.path.join(ca.ROOT, "configs", "accelerators", f"{ca.ACCEL}.cfg")).read()
    freq = float(re.search(r"(?m)^frequency\s*=\s*([\d.]+)", cfg[cfg.index("[dram]"):]).group(1))
    base = {}
    if a.baseline:
        base = {(int(r["B"]), int(r["L"]), r["backend_env"]): r for r in csv.DictReader(open(os.path.join(a.baseline, "kv_backend_probe.csv")))}
    checks = [("inputs unchanged during the runs (binary, libraries, yaml, config, scripts)", before == after)]
    near = lambda x, y: x is not None and y is not None and abs(x - y) <= 0.051      # reports print 0.1 cycle
    for r in rows:
        tag = r["run_dir"]
        checks.append((f"{tag}: reported KV read / append bytes = geometry",
                       r["kv_read_bytes_reported"] == r["kv_dense_bytes_reported"] == r["kv_read_bytes_geometry"]
                       and r["kv_append_bytes_reported"] == r["kv_append_bytes_geometry"]))
        checks.append((f"{tag}: KV supply = K + V supply, each = its DRAM read + half the decoder",
                       abs(r["kv_supply_dram_cycles"] - r["k_supply_cycles"] - r["v_supply_cycles"]) <= 0.151 and  # three printed values
                       (r["k_dram_read_cycles"] is None or
                        (near(r["k_supply_cycles"], r["k_dram_read_cycles"] + r["decoder_cycles"] / 2) and
                         near(r["v_supply_cycles"], r["v_dram_read_cycles"] + r["decoder_cycles"] / 2)))))
        b = base.get((r["B"], r["L"], r["backend_env"])) if r["schedule"] == "streaming" else None
        if b and r["backend_env"] == "(unset)":
            checks.append((f"{tag}: reproduces the baseline run without a yaml (supply, step, projection DRAM busy; "
                           "layer report byte for byte)",
                           all(near(r[k], float(b[k])) for k in ("k_supply_cycles", "v_supply_cycles", "attention_step_cycles",
                                                                  "projection_dram_busy_cycles"))
                           and r["report_sha256"] == b["report_sha256"]))
        if b and r["backend_env"] != "(unset)":
            checks.append((f"{tag}: projection DRAM busy as in the baseline run with the yaml",
                           near(r["projection_dram_busy_cycles"], float(b["projection_dram_busy_cycles"]))))
        if r["backend_env"] != "(unset)":
            checks.append((f"{tag}: K / V read time reported as backend-derived, append labelled analytical",
                           r["k_dram_read_cycles"] is not None and r["append_label"] == "analytical"
                           and r["read_stream_bytes"] == r["kv_read_bytes_reported"] // 2))
            if a.diag and r["k_dram_read_cycles"] is not None:
                want, kind = helper_cycles(a.diag, r["read_stream_bytes"], freq)
                r["helper_stream_cycles"], r["helper_stream_kind"] = want, kind
                checks.append((f"{tag}: K and V DRAM read = the stream helper for {r['read_stream_bytes']} B ({kind})",
                               near(r["k_dram_read_cycles"], want) and near(r["v_dram_read_cycles"], want)
                               and r["k_stream"] == r["v_stream"] == kind))
    for u, m in zip(rows[::2], rows[1::2]):                   # pairs: without / with the yaml
        tag = m["run_dir"]
        checks.append((f"{tag}: QK / softmax / AV, append and bytes equal to the run without the yaml",
                       all(u[k] == m[k] for k in ("qk_softmax_av_cycles", "append_cycles", "kv_read_bytes_reported",
                                                  "kv_append_bytes_reported"))))
        checks.append((f"{tag}: K / V supply now differ from the analytical ones, and are not their sum with them",
                       m["k_supply_cycles"] != u["k_supply_cycles"] and m["k_dram_read_cycles"] is not None
                       and near(m["k_supply_cycles"], m["k_dram_read_cycles"] + m["decoder_cycles"] / 2)))
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(out, "kv_backend_probe.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    json.dump(dict(frozen_before=before, frozen_after=after, dram_frequency_mhz=freq, baseline=a.baseline, diag=a.diag,
                   git_head=subprocess.run(["git", "-C", ca.ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
                   git_dirty=subprocess.run(["git", "-C", ca.ROOT, "status", "--porcelain", "-uno"], capture_output=True, text=True).stdout.split("\n"),
                   checks=[dict(check=c, ok=bool(ok)) for c, ok in checks], rows=rows),
              open(os.path.join(out, "kv_backend_probe.json"), "w"), indent=1)
    rate = lambda r: r["kv_read_bytes_reported"] / (r["kv_supply_dram_cycles"] / ca.CLOCK_HZ) / 1e9
    for u, m in zip(rows[::2], rows[1::2]):
        print(f"B {u['B']} L {u['L']} {u['schedule']} ({u['kv_read_bytes_reported']} KV bytes read, {u['kv_append_bytes_reported']} appended):")
        print(f"   projection DRAM busy   {u['projection_dram_busy_cycles']} -> {m['projection_dram_busy_cycles']} cycles")
        print(f"   K / V supply           {u['k_supply_cycles']} / {u['v_supply_cycles']} -> {m['k_supply_cycles']} / {m['v_supply_cycles']} cycles"
              f"   (yaml: {m['k_stream']} / {m['v_stream']}, helper {m.get('helper_stream_cycles')})")
        print(f"   KV bytes / supply time {rate(u):.1f} GB/s without the yaml, {rate(m):.1f} GB/s with it")
        print(f"   QK / softmax / AV      {u['qk_softmax_av_cycles']} -> {m['qk_softmax_av_cycles']}; append {u['append_cycles']} -> {m['append_cycles']} cycles")
        print(f"   attention step         {u['attention_step_cycles']} -> {m['attention_step_cycles']} cycles")
    bad = [c for c, ok in checks if not ok]
    print(f"{len(checks) - len(bad)} of {len(checks)} checks passed" + "".join(f"\n   FAILED: {c}" for c in bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
