#!/usr/bin/env python3
"""NPUsim on the loop cases of a `measure_tpu_v5e1.py --only balance` run: three ablations per case.

Each GEMM is simulated WHOLE (every 512-column strip a DRAM iteration; no strip extrapolation), in-program
(program_setup_cycle 0) with activations resident (activation_residency 2: the loop's output is reduced on chip;
its input read, M x K bf16 per iteration, is NOT modelled), on an unchanged DRAM model:
  A  the config as it is: off-chip stages follow the mapping's tiles
  B  offchip_window_count = the case's XLA window count (iteration_bounds of the run's own optimized HLO)
  C  B + window_switch_cycle = 590 (0.393 us at 1.5 GHz; fixed before the balance run, both dtypes)
B - A is the window schedule, C - B the per-window cost. Measured = the loop's fused op per iteration (device trace).

What the window schedule is and is not: COUNT-BASED, EQUAL-SHARE. Only the number of windows is taken from the HLO
(an observed compiler schedule, an input -- not predicted for unseen shapes); each window gets an equal share of
the layer's off-chip cost and a near-equal run of the simulator's tiles. XLA's window geometry (kernel_window,
partial edge windows) is not represented: `window_coverage` > 1 marks cases whose last windows are partial
(windows x window size / weight size, reading kernel_window as [8 K-rows, 128 columns] units -- an assumption).
Not modelled and not absorbed by w: the loop's input read and add, the output conversion and reduction.

The binary and library are hashed before and after the batch and the batch fails if they changed; the generated
configs, networks, mappings and NPUsim reports of every run are saved under OUT/generated_<dram>/. Every run has
its own network / weight / mapping name, and check_inputs() verifies network and mapping dimensions of EVERY
queued run after all inputs are generated and again after the runs.

Usage: python3 balance_npusim.py RUN_DIR OUT_DIR [--dram YAML] [--only dtype:M:K,...] [--variants ABC] [-j 6]
"""
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import statistics as st
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import subprocess

import balance_analyze
import compare_v5e as c
palm = c.palm

W_CYCLES = 590
ACCEL = {"bf16": "tpu_v5e", "int8": "tpu_v5e_int8w"}
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()


def window_config(path):
    """(iteration_bounds, kernel_window_bounds) of the program's one convolution, as written in the HLO."""
    m = re.findall(r'"iteration_bounds":\[([^\]]*)\],"kernel_window_bounds":\[([^\]]*)\]', open(path).read())
    assert len(m) == 1, (path, m)
    return tuple([int(x.strip('"')) for x in part.split(",")] for part in m[0])


def check_inputs(cfgroot, name, accel, M, K, N):
    """Fail unless the run's network AND its mapping describe an M x K x N GEMM. gemm_files() writes
    networks/<name>.cfg and weights/<name>.wgh shared by every accel: two runs with one name silently share
    the last shape written (the bug that invalidated the first batches)."""
    net = open(os.path.join(cfgroot, "networks", name + ".cfg")).read()
    g = lambda key: int(re.search(rf"(?m)^{key}=(\d+)", net).group(1))
    got = (g("batch"), g("height") * g("width") * g("channels"), g("output"))
    assert got == (M, K, N), f"{accel}/{name}: network is {got}, wanted {(M, K, N)}"
    assert f"weight=weights/{name}.wgh" in net, f"{accel}/{name}: network points at another weight file"
    rows = [[int(x) for x in line.split("=")[1].split(",") if x.strip()]
            for line in open(os.path.join(cfgroot, "mappings", accel, name, "ws.map")) if "=" in line]
    prod = lambda col: math.prod(r[col] for r in rows)
    got = (prod(1), prod(4), prod(0))                        # B, C, K columns: rows, reduction, output columns
    assert got == (M, K, N), f"{accel}/{name}: mapping covers {got}, wanted {(M, K, N)}"


def frozen():
    """Hashes of what a batch runs on."""
    root = palm.ROOT
    git = lambda *a: subprocess.run(["git", "-C", root, *a], capture_output=True, text=True).stdout.strip()
    return dict(binary_sha256=sha(os.path.join(root, "models", "model")),
                library_sha256=sha(os.path.join(root, "library", "libnpusim.so")),
                configs={n: sha(os.path.join(root, "configs", "accelerators", n + ".cfg")) for n in ACCEL.values()},
                git_head=git("rev-parse", "HEAD"), git_dirty=git("status", "--porcelain", "-uno").splitlines())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--dram", help="MOSAIC DRAM yaml (default: the config's analytical [dram])")
    ap.add_argument("--only", help="dtype:M:K,... (default: every loop case)")
    ap.add_argument("--variants", default="ABC")
    ap.add_argument("--w", type=int, default=W_CYCLES, help="window_switch_cycle of variant C (sensitivity only)")
    ap.add_argument("-j", "--jobs", type=int, default=6)
    ap.add_argument("--generate-only", action="store_true", help="write and check every run's inputs, run nothing")
    a = ap.parse_args()
    only = a.only and {tuple(x.split(":")) for x in a.only.split(",")}
    cases = []
    for r in csv.DictReader(open(os.path.join(a.run, "balance.csv"))):
        if r["schedule"] == "loop" and (not only or (r["dtype"], r["M"], r["K"]) in only):
            hlo = os.path.join(a.run, "hlo", f"balance_{r['dtype']}_{r['M']}x{r['K']}_loop.txt")
            it, kw = window_config(hlo)
            M, K, n = int(r["M"]), int(r["K"]), balance_analyze.windows(hlo)
            cases.append(dict(dtype=r["dtype"], M=M, K=K, windows=n, iteration_bounds="x".join(map(str, it)),
                              kernel_window_bounds="x".join(map(str, kw)),
                              window_coverage=round(n * kw[0] * 8 * kw[1] * 128 / (K * K), 3),
                              window_model="count-based equal-share",
                              measured_us=float(r["fused_ops_sum_per_iter_us"]), hlo_sha256=sha(hlo)))
    before = frozen()
    tmp, cfgroot, rundir = palm.workdir(tuple(ACCEL.values()))
    palm.MOSAIC_DRAM = a.dram and os.path.abspath(a.dram)
    todo = {}
    for x in cases:
        base = open(os.path.join(cfgroot, "accelerators", ACCEL[x["dtype"]] + ".cfg")).read()
        base = re.sub(r"(?m)^program_setup_cycle\s*=.*$", "program_setup_cycle = 0", base)
        base, k = re.subn(r"(?m)^activation_residency\s*=.*$", "activation_residency = 2", base)
        assert k == 1 and "offchip_window_count" not in base
        knobs = {"A": "", "B": f"offchip_window_count = {x['windows']}\n",
                 "C": f"offchip_window_count = {x['windows']}\nwindow_switch_cycle = {a.w}\n"}
        for v in a.variants:
            accel = f"bal_{x['dtype']}_{x['M']}x{x['K']}_{v}".lower()
            text = base.replace("[accelerator]\n", "[accelerator]\n" + knobs[v], 1)
            assert v == "A" or "offchip_window_count" in text
            open(os.path.join(cfgroot, "accelerators", accel + ".cfg"), "w").write(text)
            # the run's name is its accel: network, weights, mapping and report are this run's alone
            palm.gemm_files(cfgroot, rundir, accel, accel, x["M"], x["K"], x["K"])
            todo[x["dtype"], x["M"], x["K"], v] = (cfgroot, rundir, accel, accel)
    for (dt, M, K, v), (_, _, name, accel) in todo.items():   # after ALL generation: nothing was overwritten
        check_inputs(cfgroot, name, accel, M, K, K)
    assert len({t[2] for t in todo.values()}) == len(todo)
    os.makedirs(a.out, exist_ok=True)
    tag = (os.path.splitext(os.path.basename(a.dram))[0] if a.dram else "analytic") + (f"_w{a.w}" if a.w != W_CYCLES else "")
    gen = os.path.join(a.out, "generated_" + tag)            # what every run is given
    shutil.rmtree(gen, ignore_errors=True)
    for _, _, name, accel in todo.values():
        os.makedirs(os.path.join(gen, accel))
        shutil.copy(os.path.join(cfgroot, "accelerators", accel + ".cfg"), os.path.join(gen, accel))
        shutil.copy(os.path.join(cfgroot, "networks", name + ".cfg"), os.path.join(gen, accel, "network.cfg"))
        shutil.copytree(os.path.join(cfgroot, "mappings", accel), os.path.join(gen, accel, "mappings"))
    if a.generate_only:
        shutil.rmtree(tmp)
        print(f"generated and checked {len(todo)} runs -> {gen}")
        return
    print(f"NPUsim: {len(todo)} whole-GEMM runs on {a.jobs} threads ...", flush=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        sim = dict(zip(todo, ex.map(palm.run_gemm, todo.values())))
    after = frozen()
    for (dt, M, K, v), (_, _, name, accel) in todo.items():
        check_inputs(cfgroot, name, accel, M, K, K)           # still intact after the runs
        for report in (f"{accel}_{name}.txt", f"{accel}_{name}_layer_0.txt"):   # the run's own stage / timeline report
            shutil.copy(os.path.join(rundir, report), os.path.join(gen, accel))
    shutil.rmtree(tmp)
    if before != after:
        json.dump(dict(before=before, after=after), open(os.path.join(a.out, f"MIXED_BINARY_{tag}.json"), "w"), indent=1)
        raise SystemExit("binary, library or configs changed during the batch: results discarded")
    rows, err = [], {v: [] for v in a.variants}
    print(f"{'dtype':5} {'M':>4} {'K':>5} {'n':>3} | {'measured':>8} | " + " | ".join(f"{v + ' us':>8} {'err':>6}" for v in a.variants))
    for x in cases:
        row = dict(x, dram=tag, timing_kind="device_trace(fused op per loop iteration)")
        cells = []
        for v in a.variants:
            us = sim[x["dtype"], x["M"], x["K"], v] * 1e6
            e = us / x["measured_us"] - 1
            err[v].append(({x["dtype"], "full windows" if x["window_coverage"] <= 1 else "partial edge windows"}, e))
            row[f"npusim_{v}_us"], row[f"err_{v}_pct"] = round(us, 2), round(100 * e, 1)
            accel = todo[x["dtype"], x["M"], x["K"], v][3]       # what the run itself reports, not what was asked
            m = re.search(r"Tensor residency\s*:\s*(.*)", open(os.path.join(gen, accel, f"{accel}_{accel}_layer_0.txt")).read())
            row[f"residency_{v}"] = m.group(1).strip() if m else "not reported"
            cells.append(f"{us:8.2f} {100 * e:+5.1f}%")
        rows.append(row)
        print(f"{x['dtype']:5} {x['M']:4d} {x['K']:5d} {x['windows']:3d} | {x['measured_us']:8.2f} | " + " | ".join(cells))
    for g in ("all", "bf16", "int8", "full windows", "partial edge windows"):   # edge effects are not pooled away
        sel = {v: [e for d, e in err[v] if g == "all" or g in d] for v in a.variants}
        if sel[a.variants[0]]:
            print(f"{g:20s} n={len(sel[a.variants[0]]):2d} mean / MAPE: " + "  ".join(
                f"{v} {100 * st.mean(sel[v]):+5.1f} / {100 * st.mean(map(abs, sel[v])):4.1f}%" for v in a.variants))
    with open(os.path.join(a.out, f"balance_npusim_{tag}.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    env = json.load(open(os.path.join(a.run, "env.json")))
    json.dump(dict(frozen=before, dram=a.dram and {os.path.abspath(a.dram): sha(a.dram)},
                   balance_csv_sha256=sha(os.path.join(a.run, "balance.csv")),
                   measurement=dict(run=os.path.abspath(a.run), jax=env.get("jax"), device=env.get("device_kind"),
                                    date=env.get("date"), scope="streaming loop, fused op per iteration (device trace)"),
                   mapping="whole GEMM, 512-column strips as DRAM iterations, one <= 512-row tile",
                   window_model="count-based equal-share: offchip_window_count = HLO window count; geometry and "
                                "partial edge windows not represented",
                   variants=dict(A="config as is", B="offchip_window_count", C=f"B + window_switch_cycle {a.w}"),
                   window_switch_cycle=dict(value=a.w, clock_hz=palm.CLOCK_HZ, us=a.w / palm.CLOCK_HZ * 1e6,
                                            status="empirically fitted effective cost, not a measured hardware "
                                                   "constant; mechanism unidentified",
                                            source="2026-11-04 profile run, 4 int8 M = 128 cases (0.393 us), with "
                                                   "R 752 GB/s and F 2 us fitted on earlier stream data; fixed "
                                                   "before the 2026-11-12 balance run; one value for both dtypes",
                                            regime="XLA streaming loop on TPU v5e (JAX 0.11.2); not for staged / "
                                                   "isolated one-op programs or shapes without an HLO"),
                   assumptions=["program_setup_cycle = 0 (in-program)",
                                "activation_residency = 2: output reduced on chip; the loop's input read "
                                "(M x K bf16 per iteration) is not modelled",
                                "input add and output conversion / reduction of the fused op are not modelled"]),
              open(os.path.join(a.out, f"balance_npusim_{tag}.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
