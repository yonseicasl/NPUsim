#!/usr/bin/env python3
"""TPU x 3D-stacked-DRAM design-space exploration, first pilot: one NPUsim TPU configuration, two memory models,
three layer-level workloads -- six simulations, every input and output kept. No TPU data, nothing fitted.

Chip (fixed): configs/accelerators/tpu_v5e.cfg -- four 128 x 128 MXUs, 1.5 GHz, its declared buffers and mappings;
the [dram] clock equals the compute clock (the MOSAIC stream times are converted with [dram] frequency and
consumed in that one cycle domain). Window knobs stay at their default 0.
Memory (MOSAIC Ramulator2 yaml through NPUSIM_MOSAIC_DRAM; tile streams and the attention K / V read):
  hbm2e_refsb     HBM2e beside the die -- the model compared with TPU v5e hardware; the paired baseline
  hbdram_example  hybrid-bonded DRAM on logic -- the repository's EXISTING example profile
                  (thermal_sweep/hb3d_example.yaml). An existing model scenario, not a validated device.
Workloads (synthetic, layer level -- NOT end-to-end model throughput, and no claim that a model fits the memory):
  gemm_m32     bf16 GEMM 32 x 4096 x 4096, fused_output scope (input and weights read from DRAM once, the output
  gemm_m512    bf16 GEMM 512 x 4096 x 4096  matrix never written back); program_setup_cycle 0
  attn_b8      bf16 KV attention step, B 8, L 2048 (8 KV heads x 128, two_pass, streaming; the geometry of the
               reviewed KV gate). Metric: the attention step appended after a small dummy projection.
OUT/: manifest.json (memory organisation, capacity, nominal interface rate and transaction size -- derived, see
memory_facts), preflight.json, runs/<point>/ (config root, reports, stdout, /usr/bin/time), pilot.csv / .json.

Usage: python3 dse_pilot.py MOSAIC_ROOT OUT_DIR [--generate-only] [-j 2]          (needs models/model)
"""
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tpu_v5e"))
import compare_attention as ca                      # sets the v5e geometry: 4 MXUs, 1.5 GHz
import kv_backend_probe as kvp
palm = ca.palm
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
MIB, GIB = 1 << 20, 1 << 30
ACT = 2                                              # bf16 bytes
ACC = 4                                              # fp32 accumulator bytes (the GLB output partition)
MEMORIES = {"hbm2e_refsb": ("configs/architecture/dram/hbm2e_tpu_v5e_refsb.yaml",
                            "HBM2e stacks beside the die; the model compared with TPU v5e hardware (baseline)"),
            "hbdram_example": ("configs/architecture/dram/thermal_sweep/hb3d_example.yaml",
                               "hybrid-bonded DRAM on logic; the repository's existing example profile -- a model "
                               "scenario, not a validated device")}
WORKLOADS = {"gemm_m32": dict(kind="gemm", M=32, K=4096, N=4096), "gemm_m512": dict(kind="gemm", M=512, K=4096, N=4096),
             "attn_b8": dict(kind="attention", B=8, L=2048)}


def memory_facts(root, yaml):
    """Organisation, capacity and nominal interface rate DERIVED from the yaml (and, for HBDRAM, from the preset
    tables of MOSAIC's Ramulator2 python spec); transaction size and tCK from the runtime (dram_bw_probe)."""
    text = open(yaml).read()
    probe = subprocess.run([os.path.join(root, "simulator", "build", "bin", "dram_bw_probe"), yaml, "1", "0"],
                           capture_output=True, text=True).stdout
    m = re.search(r"channels=(\d+) tx=(\d+)B tCK=([\d.]+)ns", probe)
    if not m:
        raise SystemExit(f"dram_bw_probe does not run {yaml}:\n{probe[-400:]}")
    f = dict(yaml=os.path.relpath(yaml, root), yaml_sha256=sha(yaml), runtime_channels=int(m.group(1)),
             transaction_bytes=int(m.group(2)), tck_ns=float(m.group(3)))
    hb = re.search(r"impl:\s*HBDRAM\s+org:\s+preset:\s*(\w+)\s+channel:\s*(\d+)\s+timing:\s+preset:\s*(\w+)", text)
    if hb:
        spec = open(os.path.join(root, "simulator/src/dram/ramulator2/python/ramulator/dram/hbdram.py")).read()
        org = json.loads(re.search(rf'"{hb.group(1)}":\s*(\{{[^}}]*\}})', spec).group(1))
        rate = int(re.search(rf'"{hb.group(3)}":\s*\{{"rate":\s*(\d+)', spec).group(1))
        ch = int(hb.group(2))
        bits = org["bankgroup"] * org["bank"] * org["row"] * org["column"] * org["dq"]
        f.update(impl="HBDRAM", org_preset=hb.group(1), timing_preset=hb.group(3), channels=ch, dq_per_channel=org["dq"],
                 total_dq=ch * org["dq"], rate_mbps_per_dq=rate, capacity_bytes=ch * bits // 8,
                 nominal_interface_gbps=ch * org["dq"] * rate * 1e6 / 8 / 1e9,
                 derivation="capacity = channels x bankgroup x bank x row x column x dq bits (preset table); "
                            "nominal rate = channels x dq x rate")
    else:
        counts = [json.loads(c) for c in re.findall(r"count:\s*(\[[^\]]*\])", text)]
        dq = [int(d) for d in re.findall(r"(?m)^\s*dq:\s*(\d+)", text)]
        rates = [int(t) for t in re.findall(r"timing:\s*\[(\d+),", text)]
        assert counts and len(counts) == len(dq) == len(rates) and len(set(map(tuple, counts))) == 1, "unexpected HBM yaml"
        c = counts[0]                               # channel, pseudo-channel, bank group, bank, row, column
        bits = math.prod(c) * dq[0]
        f.update(impl=re.search(r"dram:\s+impl:\s*(\w+)", text).group(1), channels=len(counts), dq_per_channel=c[1] * dq[0],
                 total_dq=len(counts) * c[1] * dq[0], rate_mbps_per_dq=rates[0], capacity_bytes=len(counts) * bits // 8,
                 nominal_interface_gbps=len(counts) * c[1] * dq[0] * rates[0] * 1e6 / 8 / 1e9,
                 derivation="capacity = controllers x prod(org count) x dq bits; nominal rate = controllers x "
                            "pseudo-channels x dq x rate")
    assert f["runtime_channels"] == f["channels"], f
    return f


def frozen(yamls):
    f = kvp.frozen(yamls[0])
    for y in yamls[1:]:
        f[y] = sha(y)
    for s in (os.path.abspath(__file__), os.path.join(HERE, "..", "tpu_v5e", "kv_backend_probe.py")):
        f[os.path.relpath(s, ca.ROOT)] = sha(s)
    return f


def live_set(w):
    """Bytes one MXU's VMEM share ([separate], declared PER MXU in the accelerator cfg) must hold for a GEMM: the
    whole activation (each MXU reads all K inputs), one 128-column weight tile, one tile of fp32 accumulators.
    Independent of the MXU count. Arithmetic on declared sizes -- the simulator does not enforce capacity."""
    return dict(Input=w["M"] * w["K"] * ACT, Weight=w["K"] * palm.DIM * ACT, Output=w["M"] * palm.DIM * ACC)


def generate(out, point, w, yaml, cfg=None, mxus=None):
    """The run directory of one point; returns (dir, network name). cfg / mxus: another TPU point's accelerator
    text and MXU count -- the mapping is REGENERATED for it (strips of mxus x 128 columns), not reused."""
    d = os.path.join(out, "runs", point)
    if os.path.exists(d):
        raise SystemExit(f"{d} exists; nothing was changed -- use a fresh OUT_DIR")
    cfgroot, rundir = os.path.join(d, "configs"), os.path.join(d, "run")
    for sub in ("accelerators", "networks", "mappings"):
        os.makedirs(os.path.join(cfgroot, sub))
    os.makedirs(os.path.join(rundir, "weights"))
    os.symlink(os.path.join(ca.ROOT, "models", "datasets"), os.path.join(rundir, "datasets"))
    cfg = cfg or open(os.path.join(ca.ROOT, "configs", "accelerators", f"{ca.ACCEL}.cfg")).read()
    assert "offchip_window_count" not in cfg and "window_switch_cycle" not in cfg     # default 0: not a v5e fit
    if w["kind"] == "gemm":
        name, dims = f"gemm_{w['M']}x{w['K']}x{w['N']}", (w["M"], w["K"], w["N"])
        cfg = re.sub(r"(?m)^program_setup_cycle\s*=.*$", "program_setup_cycle = 0", cfg)
        cfg, k = re.subn(r"(?m)^activation_residency\s*=.*$", "activation_residency = 0\nfused_output_consumer = 1", cfg)
        assert k == 1
    else:
        name, dims = f"attn_b{w['B']}_l{w['L']}", (w["B"], 128, 512)           # the dummy projection of the gate
        cfg += ca.KVCACHE.format(qh=ca.QH, kvh=ca.KVH, dh=ca.DH, kvb=ca.KV_BYTES, L=w["L"], B=w["B"])
    open(os.path.join(cfgroot, "accelerators", f"{ca.ACCEL}.cfg"), "w").write(cfg)
    base_mxus = palm.MXUS
    palm.MXUS = mxus or base_mxus                 # gemm_files reads it: CHIPS_Y and the strip width
    try:
        palm.gemm_files(cfgroot, rundir, name, ca.ACCEL, *dims)
    finally:
        palm.MXUS = base_mxus
    json.dump(dict(point=point, workload=w, dram_yaml=yaml, network=name, mxus=mxus or base_mxus),
              open(os.path.join(d, "point.json"), "w"), indent=1)
    return d, name


def simulate(d, name, yaml, timeout):
    rundir = os.path.join(d, "run")
    env = dict(os.environ, NPUSIM_CONFIG_ROOT=os.path.join(d, "configs"))
    env.pop("NPUSIM_MOSAIC_DRAM", None)
    if yaml:
        env["NPUSIM_MOSAIC_DRAM"] = yaml
    t = time.time()
    p = subprocess.run(["/usr/bin/time", "-v", "-o", os.path.join(d, "time.txt"), os.path.join(ca.ROOT, "models", "model"),
                        "run", ca.ACCEL, name, "ws"], cwd=rundir, env=env, capture_output=True, text=True, timeout=timeout)
    open(os.path.join(d, "stdout.txt"), "w").write(p.stdout + p.stderr)
    shutil.rmtree(os.path.join(rundir, "weights"))            # zero-filled placeholder, not evidence
    if p.returncode:
        raise SystemExit(f"{d} failed ({p.returncode}): {(p.stdout + p.stderr)[-600:]}")
    return round(time.time() - t, 1), int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", open(os.path.join(d, "time.txt")).read()).group(1))


def parse(d, name, w):
    rundir = os.path.join(d, "run")
    layer = open(os.path.join(rundir, f"{ca.ACCEL}_{name}_layer_0.txt")).read()
    stdout = open(os.path.join(d, "stdout.txt")).read()
    cfg = open(os.path.join(d, "configs", "accelerators", f"{ca.ACCEL}.cfg")).read()
    f = lambda rx, t=layer: (lambda m: None if m is None else float(m.group(1)))(re.search(rx, t))
    r = dict(dram_timing_model=re.search(r"DRAM timing model\s*:\s*(.*)", layer).group(1).strip(),
             residency=re.search(r"Tensor residency\s*:\s*(.*)", layer).group(1).strip(),
             layer_critical_path_cycles=f(r"Critical-path latency\s*:\s*([\d.]+)"),
             bottleneck=re.search(r"Bottleneck level\s*:\s*(.*)", layer).group(1).strip())
    busy = layer[layer.index("Busy cycles (ratio of critical path)"):]
    for key, label in (("dram", r"DRAM"), ("nop", r"Multi-chip \(NoP\)"), ("glb", r"Global buffer"), ("pe_array", r"PE array"),
                       ("pe_compute", r"PE \(compute\+LB\)")):
        m = re.search(rf"\* {label}\s*:\s*([\d.]+) cycles \(\s*([\d.]+) %\)", busy)
        r[f"busy_{key}_cycles"], r[f"busy_{key}_pct"] = float(m.group(1)), float(m.group(2))
    blk = layer[layer.index("Multi-chip <-> DRAM transactions"):][:400]
    bits = [int(b) for b in re.search(r"(?m)^line_size\s*=\s*([\d:]+)", cfg[cfg.index("[dram]"):]).group(1).split(":")]
    tx = [int(re.search(rf"\* {k}\s*:\s*(\d+)", blk).group(1)) for k in ("Input data", "Weight", "Output data")]
    r["dram_payload_bytes"] = dict(zip(("input", "weight", "output"), (x * b // 8 for x, b in zip(tx, bits))))
    r["scalar_macs"] = f(r"Busy scalar-MAC cycles\s*:\s*([\d.]+)")
    glb = {k: int(float(re.search(rf"{k} size\s*:\s*([\d.]+) KB", stdout[stdout.index("Global buffer"):]).group(1)) * 1024)
           for k in ("Input", "Weight", "Output")}
    r["glb_declared_bytes"] = glb
    mp = [[int(x) for x in l.split("=")[1].split(",")[:5]] for l in open(os.path.join(d, "configs", "mappings", ca.ACCEL, name, "ws.map"))
          if re.match(r"\s*\w+\s*=", l)]
    r["mapping_products_KBC"] = [math.prod(row[i] for row in mp) for i in (0, 1, 4)]
    if w["kind"] == "gemm":
        M, K, N = w["M"], w["K"], w["N"]
        r.update(metric="GEMM critical path", latency_cycles=r["layer_critical_path_cycles"],
                 expected_payload_bytes=dict(input=M * K * ACT, weight=K * N * ACT, output=0), expected_macs=M * K * N,
                 live_set_bytes=live_set(w),
                 approximation="tile streams: a tile is charged its share of a >= 1 MiB sequential stream on a fresh "
                               "device (amortized start); contiguous, uncontended")
        r["mapping_ok"] = r["mapping_products_KBC"] == [N, M, K]
    else:
        B, L = w["B"], w["L"]
        a = kvp_row(layer)
        r.update(a, metric="attention step (appended after a dummy projection)", latency_cycles=a["attention_step_cycles"],
                 expected_kv_read_bytes=B * L * 2 * ca.KVH * ca.DH * ca.KV_BYTES,
                 expected_kv_append_bytes=B * 2 * ca.KVH * ca.DH * ca.KV_BYTES,
                 approximation=f"K and V: two fresh independent read streams ({a['k_stream']} / {a['v_stream']}), each "
                               "total spread uniformly over the supply tiles; QK / AV / softmax rates are configured "
                               "assumptions; append, energy and link / row counts analytical; uncompressed bf16")
    return r


def kvp_row(layer):
    g = lambda rx: re.search(rx, layer)
    rd = g(r"K / V DRAM read\s*:\s*([\d.]+) / ([\d.]+) cycles, decoder excluded \((\d+) B per pass; K (\w+), V (\w+)\)")
    dense, comp = (int(x) for x in g(r"KV bytes dense/comp\s*:\s*(\d+) / (\d+) B").groups())
    qk = g(r"QK / softmax / AV\s*:\s*([\d.]+) / ([\d.]+) / ([\d.]+)").groups()
    ex = g(r"KV exposed / hidden\s*:\s*([\d.]+) / ([\d.]+)").groups()
    return dict(attention_step_cycles=float(g(r"Attention step\s*:\s*([\d.]+)").group(1)),
                k_dram_read_cycles=rd and float(rd.group(1)), v_dram_read_cycles=rd and float(rd.group(2)),
                read_stream_bytes=rd and int(rd.group(3)), k_stream=rd and rd.group(4), v_stream=rd and rd.group(5),
                kv_read_bytes_reported=comp, kv_append_bytes_reported=int(g(r"cache occupancy\s*:\s*(\d+) B").group(1)) - dense,
                qk_cycles=float(qk[0]), softmax_cycles=float(qk[1]), av_cycles=float(qk[2]),
                kv_exposed_cycles=float(ex[0]), kv_hidden_cycles=float(ex[1]),
                kv_stall_cycles=float(g(r"buffer-full stall\s*:\s*([\d.]+)").group(1)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mosaic_root")
    ap.add_argument("out")
    ap.add_argument("--generate-only", action="store_true")
    ap.add_argument("-j", "--jobs", type=int, default=2)
    ap.add_argument("--timeout", type=float, default=7200)
    a = ap.parse_args()
    root, out = os.path.abspath(a.mosaic_root), os.path.abspath(a.out)
    assert a.jobs <= 2, "the pilot is capped at 2 concurrent simulations"
    if os.path.exists(os.path.join(out, "manifest.json")):
        raise SystemExit(f"{out} already holds a pilot; nothing was changed -- use a fresh OUT_DIR")
    os.makedirs(out, exist_ok=True)
    yamls = {m: os.path.join(root, rel) for m, (rel, _) in MEMORIES.items()}
    mem = {m: dict(memory_facts(root, y), architecture=MEMORIES[m][1]) for m, y in yamls.items()}
    cfg = open(os.path.join(ca.ROOT, "configs", "accelerators", f"{ca.ACCEL}.cfg")).read()
    dram_mhz = float(re.search(r"(?m)^frequency\s*=\s*([\d.]+)", cfg[cfg.index("[dram]"):]).group(1))
    assert dram_mhz * 1e6 == ca.CLOCK_HZ, "the pilot keeps the [dram] clock equal to the compute clock"
    before = frozen(list(yamls.values()))
    points = {f"{wn}__{m}": (w, m) for wn, w in WORKLOADS.items() for m in MEMORIES}
    json.dump(dict(chip=dict(config=f"configs/accelerators/{ca.ACCEL}.cfg", mxus=palm.MXUS, array=f"{palm.DIM} x {palm.DIM}",
                             clock_hz=ca.CLOCK_HZ, dram_clock_mhz=dram_mhz, window_knobs="default 0",
                             dma_knobs="none set: 64-transaction requests, unlimited in flight, no path latency"),
                   memories=mem, workloads=WORKLOADS, points={p: dict(workload=w, memory=m) for p, (w, m) in points.items()},
                   frozen_before=before,
                   git_head=subprocess.run(["git", "-C", ca.ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
                   git_dirty=subprocess.run(["git", "-C", ca.ROOT, "status", "--porcelain", "-uno"], capture_output=True, text=True).stdout.split("\n")),
              open(os.path.join(out, "manifest.json"), "w"), indent=1)
    dirs = {p: generate(out, p, w, yamls[m]) for p, (w, m) in points.items()}
    pre_d, pre_name = generate(out, "preflight_gemm_m32__analytic", WORKLOADS["gemm_m32"], None)
    kv = WORKLOADS["attn_b8"]
    need = max(WORKLOADS["gemm_m512"]["K"] * WORKLOADS["gemm_m512"]["N"] * ACT,
               kv["B"] * kv["L"] * 2 * ca.KVH * ca.DH * ca.KV_BYTES)
    pre = dict(largest_live_tensor_bytes=need, capacity_ok={m: f["capacity_bytes"] >= need for m, f in mem.items()},
               note="capacity is checked for the pilot's own tensors only; whether a whole model fits is NOT examined")
    if a.generate_only:
        json.dump(pre, open(os.path.join(out, "preflight.json"), "w"), indent=1)
        print(f"generated {len(dirs)} points under {out}/runs; nothing simulated")
        return
    # preflight: the memory-heavy GEMM on the analytical DRAM -- run-time estimate and the reference that the
    # MOSAIC runs must differ from
    wall, rss = simulate(pre_d, pre_name, None, a.timeout)
    ref = parse(pre_d, pre_name, WORKLOADS["gemm_m32"])
    pre.update(analytic_run=dict(wall_s=wall, peak_rss_kb=rss, busy_dram_cycles=ref["busy_dram_cycles"],
                                 latency_cycles=ref["latency_cycles"], dram_timing_model=ref["dram_timing_model"]),
               runtime_estimate="about the analytic run per GEMM M 32 point; front-end load grows with the batch (M 512, B 8)")
    json.dump(pre, open(os.path.join(out, "preflight.json"), "w"), indent=1)
    print(f"preflight: analytic GEMM M 32 in {wall} s, {rss / 1e6:.1f} GB peak; capacity {pre['capacity_ok']}", flush=True)
    if not all(pre["capacity_ok"].values()):
        raise SystemExit("preflight: a memory model cannot hold the pilot's tensors")
    with ThreadPoolExecutor(a.jobs) as ex:
        res = dict(zip(points, ex.map(lambda p: simulate(*dirs[p], yamls[points[p][1]], a.timeout), points)))
    after = frozen(list(yamls.values()))
    rows, checks = [], [("inputs unchanged during the runs (binary, every loaded library, yaml, config, scripts)", before == after)]
    for p, (w, m) in points.items():
        r = parse(dirs[p][0], dirs[p][1], w)
        r.update(point=p, workload=next(k for k, v in WORKLOADS.items() if v is w), memory=m, wall_s=res[p][0], peak_rss_kb=res[p][1],
                 latency_us=r["latency_cycles"] / ca.CLOCK_HZ * 1e6)
        rows.append(r)
        checks.append((f"{p}: the report names the selected yaml as the DRAM timing model", os.path.basename(yamls[m]) in r["dram_timing_model"]))
        if w["kind"] == "gemm":
            checks.append((f"{p}: MACs = M x K x N, DRAM payload = input + weights once, no output; mapping products = N, M, K",
                           r["scalar_macs"] == r["expected_macs"] and r["dram_payload_bytes"] == r["expected_payload_bytes"] and r["mapping_ok"]))
            checks.append((f"{p}: live tensors fit one MXU's declared VMEM share (activation, weight tile, accumulators)",
                           all(r["live_set_bytes"][k] <= r["glb_declared_bytes"][k] for k in ("Input", "Weight", "Output"))))
        else:
            checks.append((f"{p}: reported KV read / append bytes = geometry; K / V read backend-derived and exact",
                           r["kv_read_bytes_reported"] == r["expected_kv_read_bytes"] and r["kv_append_bytes_reported"] == r["expected_kv_append_bytes"]
                           and r["k_dram_read_cycles"] is not None and r["k_stream"] == r["v_stream"] == "exact"))
    by = {(r["workload"], r["memory"]): r for r in rows}
    for wn, w in WORKLOADS.items():
        b, h = by[wn, "hbm2e_refsb"], by[wn, "hbdram_example"]
        h["speedup_vs_hbm2e_refsb"], b["speedup_vs_hbm2e_refsb"] = b["latency_cycles"] / h["latency_cycles"], 1.0
        key = "busy_dram_cycles" if w["kind"] == "gemm" else "k_dram_read_cycles"
        checks.append((f"{wn}: the two memory models give different DRAM time (the backend processes the reads)", b[key] != h[key]))
        if w["kind"] == "gemm":
            checks.append((f"{wn}: same MACs, payload and compute busy on both memories (operation scope identical)",
                           all(b[k] == h[k] for k in ("scalar_macs", "dram_payload_bytes", "busy_pe_compute_cycles"))))
        else:
            checks.append((f"{wn}: QK / softmax / AV identical on both memories", all(b[k] == h[k] for k in ("qk_cycles", "softmax_cycles", "av_cycles"))))
    checks.append(("gemm_m32: both MOSAIC runs differ from the analytical-DRAM preflight in DRAM busy",
                   all(by["gemm_m32", m]["busy_dram_cycles"] != ref["busy_dram_cycles"] for m in MEMORIES)))
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(out, "pilot.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, keys)
        wr.writeheader()
        wr.writerows([{k: json.dumps(v) if isinstance(v, (dict, list)) else v for k, v in r.items()} for r in rows])
    json.dump(dict(frozen_after=after, checks=[dict(check=c, ok=bool(ok)) for c, ok in checks], rows=rows),
              open(os.path.join(out, "pilot.json"), "w"), indent=1)
    for m, f in mem.items():
        print(f"{m}: {f['impl']}, {f['channels']} channels x {f['dq_per_channel']} DQ x {f['rate_mbps_per_dq']} Mbps = "
              f"{f['nominal_interface_gbps']:.1f} GB/s nominal, {f['capacity_bytes'] / GIB:.0f} GiB, transaction {f['transaction_bytes']} B, tCK {f['tck_ns']} ns")
    for wn in WORKLOADS:
        b, h = by[wn, "hbm2e_refsb"], by[wn, "hbdram_example"]
        print(f"{wn:10s} {b['metric']}: {b['latency_us']:10.2f} us -> {h['latency_us']:10.2f} us, x{h['speedup_vs_hbm2e_refsb']:.3f}"
              + (f"   DRAM busy {b['busy_dram_pct']}% -> {h['busy_dram_pct']}%, compute busy {b['busy_pe_compute_pct']}% -> {h['busy_pe_compute_pct']}%, "
                 f"bottleneck {b['bottleneck']} -> {h['bottleneck']}" if "expected_macs" in b else
                 f"   K read {b['k_dram_read_cycles']} -> {h['k_dram_read_cycles']} cycles, KV exposed {b['kv_exposed_cycles']} -> {h['kv_exposed_cycles']}"))
    bad = [c for c, ok in checks if not ok]
    print(f"{len(checks) - len(bad)} of {len(checks)} checks passed" + "".join(f"\n   FAILED: {c}" for c in bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
