#!/usr/bin/env python3
"""GEMMs with M > 512 on NPUsim: the whole GEMM in one run against the row-tile pieces compare_v5e.py sums.

Per shape three runs on the SAME conventions (tpu_v5e.cfg as is -- isolated program: input from DRAM, output
materialized to DRAM -- with program_setup_cycle 0, every 128-column strip simulated, nothing extrapolated):
  whole  rows folded <= 512 per PE pass, the mt row tiles as the DRAM-level batch factor of the mapping
         (check_tiling.py's construction);
  first  one row tile (M / mt rows) as its own GEMM;
  inner  the same with layer_setup_cycle 0 (compare_v5e.py's `_inner` config).
  whole_glb (optional, --kinds)  as whole, but the row tiles as the GLB-level batch factor.
compare_v5e.py's estimate is first + (mt - 1) x inner (with strips extrapolated on top, not done here).
OUT/<run>/ keeps the config root, both reports, stdout / stderr and result.json (exit code, wall seconds,
counters, DRAM payload). Nothing is compared with TPU data here and nothing is fitted.

The dataset list: models/datasets/gemm/test.lst has 1024 (identical) lines; the front-end reads one per batch
row and runs off the list for batch > 1024 (SIGSEGV in nebula convolutional_t::load_data) -- that was the
"whole run fails at 4 row tiles". Each run gets a list of its own batch size; --list-lines 1024 reproduces it.

Usage: python3 row_tile_repro.py OUT_DIR M,K,N [...] [--kinds whole,first,inner] [-j 3] [--timeout S]
       python3 row_tile_repro.py OUT_DIR --table                                  (needs models/model)
"""
import argparse
import glob
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_v5e as c
palm = c.palm


def run(out, kind, M, K, N, timeout, list_lines=None):
    mt = math.ceil(M / c.MTILE)
    assert M % mt == 0 and N % (palm.MXUS * c.DIM) == 0, "whole row tiles and strips only"
    rows = M if kind.startswith("whole") else M // mt         # batch of this run
    name = f"{kind}_{M}x{K}x{N}"
    d = os.path.join(out, name)
    if os.path.exists(d):
        raise SystemExit(f"{d} exists; nothing was changed -- use a fresh OUT_DIR")
    cfgroot, rundir = os.path.join(d, "configs"), os.path.join(d, "run")
    for sub in ("accelerators", "networks", "mappings"):
        os.makedirs(os.path.join(cfgroot, sub))
    os.makedirs(os.path.join(rundir, "weights"))
    src = os.path.join(palm.ROOT, "models", "datasets")
    os.makedirs(os.path.join(rundir, "datasets", "gemm"))
    os.symlink(os.path.join(src, "imagenet"), os.path.join(rundir, "datasets", "imagenet"))
    shutil.copy(os.path.join(src, "gemm", "labels.lst"), os.path.join(rundir, "datasets", "gemm"))
    line = open(os.path.join(src, "gemm", "test.lst")).readline()
    lines = rows if list_lines is None else list_lines
    open(os.path.join(rundir, "datasets", "gemm", "test.lst"), "w").write(line * lines)
    cfg = open(os.path.join(palm.ROOT, "configs", "accelerators", f"{c.ACCEL}.cfg")).read()
    cfg = re.sub(r"(?m)^program_setup_cycle\s*=.*$", "program_setup_cycle = 0", cfg)
    if kind == "inner":
        cfg, n = re.subn(r"(?m)^layer_setup_cycle\s*=.*$", "layer_setup_cycle = 0", cfg)
        assert n == 1
    open(os.path.join(cfgroot, "accelerators", f"{c.ACCEL}.cfg"), "w").write(cfg)
    palm.gemm_files(cfgroot, rundir, name, c.ACCEL, rows, K, N)
    mf = os.path.join(cfgroot, "mappings", c.ACCEL, name, "ws.map")
    if kind.startswith("whole") and mt > 1:
        s = open(mf).read().replace(f"PE      = 1, {M},", f"PE      = 1, {M // mt},")
        if kind == "whole":                                   # row tiles at the DRAM level
            s, n = re.subn(r"DRAM    = (\d+), 1,", lambda m: f"DRAM    = {m.group(1)}, {mt},", s)
        else:                                                 # whole_glb: row tiles at the GLB level
            s, n = re.subn(r"GLB     = 1, 1,", f"GLB     = 1, {mt},", s)
        assert n == 1
        open(mf, "w").write(s)
    env = dict(os.environ, NPUSIM_CONFIG_ROOT=cfgroot)
    env.pop("NPUSIM_MOSAIC_DRAM", None)
    t0 = time.time()
    try:
        p = subprocess.run([os.path.join(palm.ROOT, "models", "model"), "run", c.ACCEL, name, "ws"], cwd=rundir,
                           env=env, capture_output=True, text=True, timeout=timeout)
        code, so, se = p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired as e:
        dec = lambda b: b.decode(errors="replace") if isinstance(b, bytes) else (b or "")
        code, so, se = "timeout", dec(e.stdout), dec(e.stderr)
    open(os.path.join(d, "stdout.txt"), "w").write(so)
    open(os.path.join(d, "stderr.txt"), "w").write(se)
    shutil.rmtree(os.path.join(rundir, "weights"))            # zero-filled placeholder, not evidence
    res = dict(kind=kind, M=M, K=K, N=N, row_tiles=mt, batch_rows=rows, list_lines=lines, exit=code,
               wall_s=round(time.time() - t0, 1),
               mapping=[l.strip() for l in open(mf) if re.match(r"\s*(DRAM|GLB|PE)\s*=", l)])
    rep = os.path.join(rundir, f"{c.ACCEL}_{name}_layer_0.txt")
    if code == 0 and os.path.exists(rep):
        t = open(rep).read()
        f = lambda rx: float(re.search(rx, t).group(1))
        blk = t[t.index("Multi-chip <-> DRAM transactions"):][:400]
        bits = [int(b) for b in re.search(r"(?m)^line_size\s*=\s*([\d:]+)", cfg[cfg.index("[dram]"):]).group(1).split(":")]
        tx = [int(re.search(rf"\* {k}\s*:\s*(\d+)", blk).group(1)) for k in ("Input data", "Weight", "Output data")]
        got = [x * b // 8 for x, b in zip(tx, bits)]
        want = [rows * K * 2, K * N * 2, rows * N * 2]        # bf16 input, weights (read once), output
        res.update(critical_path_cycles=f(r"Critical-path latency\s*:\s*([\d.]+)"),
                   dram_busy_cycles=f(r"\* DRAM\s*:\s*([\d.]+) cycles \("),
                   scalar_macs=f(r"Busy scalar-MAC cycles\s*:\s*([\d.]+)"), expected_macs=rows * K * N,
                   dram_payload_bytes=dict(zip(("input", "weight", "output"), got)),
                   expected_payload_bytes=dict(zip(("input", "weight", "output"), want)),
                   residency=re.search(r"Tensor residency\s*:\s*(.*)", t).group(1).strip())
        res["checks_ok"] = res["scalar_macs"] == res["expected_macs"] and got == want
    json.dump(res, open(os.path.join(d, "result.json"), "w"), indent=1)
    return res


def table(out):
    """Per shape: the GLB-level whole run (exact payload) against first + (mt - 1) x inner, and the DRAM-level
    whole run where it exists, with the DRAM payload each construction moves. Unrounded cycles in table.json."""
    rs = {(r["M"], r["K"], r["N"], r["kind"]): r for r in map(json.load, map(open, glob.glob(os.path.join(out, "*", "result.json"))))}
    us = lambda cyc: cyc / palm.CLOCK_HZ * 1e6
    mb = lambda r: "/".join(f"{r['dram_payload_bytes'][k] / r['expected_payload_bytes'][k]:g}" for k in ("input", "weight", "output"))
    rows = []
    print("payload x: DRAM payload moved / tensor size, input/weight/output (1/1/1 = each tensor once)")
    print("shape                 tiles | whole, GLB level us (payload x) | first + (mt-1) x inner us (weights x) | split - whole"
          " | per extra tile, cycles | whole, DRAM level us (payload x)")
    for M, K, N in sorted({k[:3] for k in rs}):
        g, f, i, d = (rs.get((M, K, N, k)) for k in ("whole_glb", "first", "inner", "whole"))
        ok = lambda r: bool(r) and r["exit"] == 0
        if not (ok(f) and ok(i) and (ok(g) or ok(d))):
            print(f"{M} x {K} x {N}: incomplete")
            continue
        mt = f["row_tiles"]
        split = f["critical_path_cycles"] + (mt - 1) * i["critical_path_cycles"]
        ref = g if ok(g) else d                                # N <= 512: both whole runs coincide
        row = dict(M=M, K=K, N=N, row_tiles=mt, split_cycles=split, split_weight_payload_x=mt,
                   whole_glb_cycles=ok(g) and g["critical_path_cycles"], whole_glb_payload_ok=ok(g) and g["checks_ok"],
                   whole_dram_cycles=ok(d) and d["critical_path_cycles"], whole_dram_payload_ok=ok(d) and d["checks_ok"],
                   pieces_payload_ok=f["checks_ok"] and i["checks_ok"],
                   split_minus_whole_pct=100 * (split / ref["critical_path_cycles"] - 1),
                   per_extra_tile_cycles=(split - ref["critical_path_cycles"]) / (mt - 1))
        rows.append(row)
        print(f"{M:5d} x {K:5d} x {N:5d}  {mt:3d}  | " + (f"{us(g['critical_path_cycles']):9.2f} ({mb(g)})" if ok(g) else "      not run   ")
              + f" | {us(split):9.2f} ({mt}) | {row['split_minus_whole_pct']:+7.2f}% | {row['per_extra_tile_cycles']:8.1f} | "
              + (f"{us(d['critical_path_cycles']):9.2f} ({mb(d)})" if ok(d) else "not run"))
    json.dump(rows, open(os.path.join(out, "table.json"), "w"), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("shapes", nargs="*")
    ap.add_argument("--kinds", default="whole,first,inner")
    ap.add_argument("--timeout", type=float, default=14400)
    ap.add_argument("--list-lines", type=int, help="lines of the dataset list (default: the batch); 1024 reproduces the crash")
    ap.add_argument("-j", "--jobs", type=int, default=3)
    ap.add_argument("--table", action="store_true")
    a = ap.parse_args()
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)
    todo = [(k, *(int(x) for x in s.split(","))) for s in a.shapes for k in a.kinds.split(",")]
    with ThreadPoolExecutor(a.jobs) as ex:
        for r in ex.map(lambda t: run(out, *t, a.timeout, a.list_lines), todo):
            print({k: v for k, v in r.items() if k != "mapping"}, flush=True)
    if a.table:
        table(out)


if __name__ == "__main__":
    main()
