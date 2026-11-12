#!/usr/bin/env python3
"""NPUsim (configs/accelerators/tpu_v5e.cfg) vs Cloud TPU v5e measurements from measure_tpu_v5e1.py.

  GEMM    every gemm.csv shape on NPUsim with validation/palm's mapping (output columns over the 4 MXUs
          in 128-column weight strips). Rows go in tiles of <= 512; one tile is simulated at 2 and 4
          strips (1 or 2 when that is all) and scaled, since the time is linear in strips from the
          second on. Per-run constants are charged once: program_setup_cycle once per program (not at all
          with --in-program), layer_setup_cycle by the first row tile only.
          ROW TILES ARE AN APPROXIMATION, NOT VALIDATED IN GENERAL: for M > 512 each <= 512-row tile is
          simulated as its own run and the tiles are summed, so every tile pays its own pipeline fill and
          drain. Against a direct run of the whole GEMM this was +1.5% at 1024 x 4096 x 512 (one early
          point) but +21.7% at 1024 x 128 x 1024 (results/2026-11-11_codex_audit/whole_gemm_small.log):
          the error grows as the per-tile compute shrinks. Unresolved (row-boundary fill / reuse); such
          rows are flagged in the log. check_tiling.py measures it for a given shape. Against the measured device time (wall time without a trace).
  decode  decode.csv per-scope times: the four projection GEMMs vs the rest (attention, norms, SwiGLU,
          psum) -- the time NPUsim's PaLM model leaves out -- and, at tp=1, NPUsim's projection GEMMs
          against the measured GEMM scopes.

Usage: python3 validation/tpu_v5e/compare_v5e.py RESULTS_DIR [-j 16] [--dram hbm2e_v5e.yaml]   (needs models/model)
Writes RESULTS_DIR/compare_<gemm csv>[_<accel>][_<yaml name>].csv.
"""
import argparse
import csv
import math
import os
import re
import shutil
import statistics
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "palm"))
import palm_decode as palm  # noqa: E402

palm.MXUS, palm.CLOCK_HZ = 4, 1.5e9        # TPU v5e; palm.gemm_files / run_gemm read these
ACCEL, DIM, MTILE = "tpu_v5e", 128, 512
IN_PROGRAM = False
CHAIN = False                              # per-kernel cost inside a chain: minus the inter-layer overlap


def run_sim(args):
    """(seconds, prologue s, tail s) of one NPUsim GEMM run (prologue/tail from its layer report)."""
    sec = palm.run_gemm(args)
    cfgroot, rundir, name, accel = args
    m = re.search(r"Overlap prologue/tail\s*:\s*([\d.]+)\s*/\s*([\d.]+)",
                  open(os.path.join(rundir, f"{accel}_{name}_layer_0.txt")).read())
    head, tail = (float(m.group(1)), float(m.group(2))) if m else (0.0, 0.0)
    return sec, head / palm.CLOCK_HZ, tail / palm.CLOCK_HZ
PEAK, HBM = 4 * DIM * DIM * 2 * 1.5e9, 819e9
D, HKV, DH = 4096, 8, 128                  # measure_tpu_v5e1.py's Llama-3-8B decode layer
PROJ = {"qkv": (4096, 6144), "o_proj": (4096, 4096), "ffn_up": (4096, 28672), "ffn_down": (14336, 4096)}


def num(r, k):
    return float(r[k]) if r.get(k) not in (None, "") else None


TIMING_KINDS = ("device_trace", "device_slope", "loop_slope", "wall")


def timing_kind(rows):
    """Best guess at what a GEMM csv's time column is, when --timing does not state it. A guess is marked
    "(inferred)"; anything not recognizable is "unknown" -- a device_us column alone proves nothing
    (gemm_int8.csv holds loop slopes in it).
    wall          no device_us on any row: host wall time of the call
    loop_slope    wall_us == device_us on every row: tpu_v5e1_stream_bench.py's exporter writes its per-iteration
                  wall slope into both columns (a loop iteration with the body's other fused ops, not a kernel)
    device_trace  every row has device_us and the op events of a parsed profiler trace (`ops`), as
                  measure_tpu_v5e1.py writes them"""
    has = lambda r, k: r.get(k) not in (None, "")
    if not any(has(r, "device_us") for r in rows):
        return "wall"
    if not all(has(r, "device_us") for r in rows):
        return "unknown"
    if all(has(r, "wall_us") and float(r["wall_us"]) == float(r["device_us"]) for r in rows):
        return "loop_slope(inferred)"
    if all(has(r, "ops") for r in rows):
        return "device_trace(inferred)"
    return "unknown"


def scale_time(sim, mt, first, inner, chain=False):
    """Kernel seconds of a GEMM from its simulated pieces. sim: {run key: (seconds, prologue s, tail s)};
    first / inner: (run keys, fit) of the first row tile (pays the layer setup) and of the other mt - 1 tiles
    (no layer setup; None when mt == 1); fit None for one run, else ((s0, s1), strips): the two runs have s0 and
    s1 strips and the time is linear in strips. chain: subtract the overlap with the previous kernel,
    min(prologue, tail) of the first run."""
    def tile(names, fit):
        if fit is None:
            return sim[names[0]][0]
        (s0, s1), strips = fit
        t0, t1 = sim[names[0]][0], sim[names[1]][0]
        return t0 + (strips - s0) * (t1 - t0) / (s1 - s0)
    hidden = min(sim[first[0][0]][1], sim[first[0][0]][2]) if chain else 0.0
    return tile(*first) + ((mt - 1) * tile(*inner) if inner else 0.0) - hidden


def simulate(shapes, jobs):
    """({(M, K, N): kernel seconds}, program setup seconds) on NPUsim. A GEMM is simulated in pieces (row tiles x
    column strips) and scaled, so the constants of a run are kept out of the pieces and charged once:
    program_setup_cycle is zeroed in every simulated run and returned for the caller to add once per program;
    layer_setup_cycle is charged by the first row tile only (the others run on an `_inner` config without it)."""
    INNER = ACCEL + "_inner"
    tmp, cfgroot, rundir = palm.workdir((ACCEL,))
    p = os.path.join(cfgroot, "accelerators", f"{ACCEL}.cfg")
    text = open(p).read()                     # read before reopening for write (which truncates)
    m = re.search(r"(?m)^program_setup_cycle\s*=\s*([\d.]+)", text)
    program_setup = float(m.group(1)) / palm.CLOCK_HZ if m else 0.0
    text = re.sub(r"(?m)^program_setup_cycle\s*=.*$", "program_setup_cycle = 0", text)
    if CHAIN:                                 # mid-chain kernel: its input and output stay on chip
        text = re.sub(r"(?m)^activation_residency\s*=.*$", "activation_residency = 2", text)
    open(p, "w").write(text)
    open(os.path.join(cfgroot, "accelerators", f"{INNER}.cfg"), "w").write(
        re.sub(r"(?m)^layer_setup_cycle\s*=.*$", "layer_setup_cycle = 0", text))
    plan, todo = {}, {}
    mxus = palm.MXUS

    def add(mtile, K, ncols, chips, ksplit=1, accel=None):
        """Queue a GEMM of ncols output columns over `chips` MXUs, each MXU group also splitting the reduction
        dimension K `ksplit` ways (partial sums); returns its run key."""
        accel = accel or ACCEL
        name = f"v5e_{mtile}x{K}x{ncols}_{chips}x{ksplit}mxu"
        if (accel, name) not in todo:
            palm.MXUS = chips                 # gemm_files reads it: CHIPS_Y = chips
            palm.gemm_files(cfgroot, rundir, name, accel, mtile, K, ncols)
            palm.MXUS = mxus
            if ksplit > 1:                    # CHIPS_Y also takes a C (= K) factor; the GLB keeps the rest
                p = os.path.join(cfgroot, "mappings", accel, name, "ws.map")
                s = open(p).read()
                s = s.replace(f"CHIPS_Y = {chips}, 1, 1, 1, 1,", f"CHIPS_Y = {chips}, 1, 1, 1, {ksplit},")
                s = re.sub(r"GLB     = 1, 1, 1, 1, (\d+),",
                           lambda m: f"GLB     = 1, 1, 1, 1, {math.ceil(int(m.group(1)) / ksplit)},", s)
                open(p, "w").write(s)
            todo[accel, name] = (cfgroot, rundir, name, accel)
        return accel, name

    for M, K, N in shapes:
        mt = math.ceil(M / MTILE)
        mtile = math.ceil(M / mt)
        nb = math.ceil(N / DIM)               # 128-column blocks
        strips = nb / mxus                    # blocks per MXU; fractional when the last round is partial
        # Fewer blocks than MXUs: nb MXUs take the columns and the rest split K (partial sums), so all 4 stream
        # weights. On one MXU its 128-element port caps an N = 128 GEMM at 384 GB/s while the TPU streams 621 GB/s
        # (8 x K x 128, K 4096 -> 16384); split 4 ways NPUsim streams 759 GB/s.
        # A partial last round is charged its column share: XLA tiles finer than whole 4-MXU rounds (K = N
        # 1024 -> 1100: measured 1.27x, column share 1.22x, whole rounds 1.6x). Whole rounds matched the absolute
        # K = N 1100-1280 times better only by covering the fixed per-kernel time NPUsim lacked.
        # 1 or 2 rounds are simulated as is, between 1 and 2 interpolated, more extrapolated from 2 and 4. Not from
        # 1 and 2: with Ramulator the first strip's DRAM fill is exposed and partly hidden from the second strip on,
        # so the 1 -> 2 step is below the steady per-strip time (8192^3 extrapolated under the compute roofline).
        def pieces(accel):
            if nb < mxus:
                ks = max(1, min(mxus // nb, math.ceil(K / DIM)))
                return [add(mtile, K, nb * DIM, nb, ks, accel)], None
            if strips in (1, 2):
                return [add(mtile, K, int(strips) * mxus * DIM, mxus, accel=accel)], None
            s = (1, 2) if strips < 2 else (2, 4)
            return [add(mtile, K, x * mxus * DIM, mxus, accel=accel) for x in s], (s, strips)
        # the first row tile pays the layer setup; the remaining mt - 1 tiles of the same kernel do not
        plan[M, K, N] = (mt, pieces(ACCEL), pieces(INNER) if mt > 1 else None)
    print(f"NPUsim: {len(todo)} runs on {jobs} threads ...", flush=True)
    with ThreadPoolExecutor(jobs) as ex:
        sim = dict(zip(todo, ex.map(run_sim, todo.values())))
    shutil.rmtree(tmp)

    return {k: scale_time(sim, *v, chain=CHAIN) for k, v in plan.items()}, program_setup


def roofline(M, K, N):
    wbytes = 1 if ACCEL.endswith("int8w") else 2           # int8 weights, bf16 activations
    return max(2 * M * K * N / PEAK, (2 * M * K + wbytes * K * N + 2 * M * N) / HBM)


def main():
    global ACCEL
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("-j", "--jobs", type=int, default=16)
    ap.add_argument("--dram", help="MOSAIC DRAM yaml (e.g. validation/tpu_v5e/hbm2e_v5e.yaml): NPUsim's DRAM tile "
                    "streams on MOSAIC's Ramulator2 instead of the analytical [dram] (needs a MOSAIC build)")
    ap.add_argument("--accel", default=ACCEL, help="configs/accelerators/<name>.cfg (tpu_v5e_int8w: int8 weights)")
    ap.add_argument("--gemm", default="gemm.csv", help="GEMM measurements in RESULTS_DIR")
    ap.add_argument("--in-program", action="store_true", help="the GEMMs ran inside a loop or chain of one "
                    "program (streamed, odd_gemm, per-kernel slopes): no program_setup_cycle per GEMM")
    ap.add_argument("--timing", choices=TIMING_KINDS, help="what the csv's time is (column timing_kind of the "
                    "output; default: guessed, see timing_kind()). loop_slope = a stream-bench loop iteration, "
                    "device_slope = slope of profiled device time over a chain")
    ap.add_argument("--chain", action="store_true", help="per-kernel slopes of a chain in one program: also "
                    "subtract the inter-layer overlap with the previous kernel (implies --in-program)")
    a = ap.parse_args()
    global IN_PROGRAM, CHAIN
    IN_PROGRAM, CHAIN = a.in_program or a.chain, a.chain
    ACCEL = a.accel
    palm.MOSAIC_DRAM = a.dram and os.path.abspath(a.dram)
    gemm = list(csv.DictReader(open(os.path.join(a.results, a.gemm))))
    kind = a.timing or timing_kind(gemm)
    print(f"measured time: {kind}" + ("" if a.timing else "  -- not stated; pass --timing"), flush=True)
    dpath = os.path.join(a.results, "decode.csv")
    decode = list(csv.DictReader(open(dpath))) if os.path.exists(dpath) else []
    shapes = {(int(r["M"]), int(r["K"]), int(r["N"])) for r in gemm}
    shapes |= {(int(r["batch"]), *kn) for r in decode if r["tp"] == "1" for kn in PROJ.values()}
    sim, program_setup = simulate(sorted(shapes), a.jobs)
    tiled = sorted(k for k in shapes if k[0] > MTILE)
    if tiled:
        print(f"WARNING: {len(tiled)} GEMM(s) with M > {MTILE} use the row-tiled estimate (one <= {MTILE}-row tile "
              f"simulated, summed over tiles): each tile pays its own pipeline fill. Observed against a direct "
              f"whole-GEMM run: +21.7% at 1024x128x1024, -4.4% at 1024x1024x1024 (check_tiling.py); there is no "
              f"general error bound. Not validated for these shapes: "
              + ", ".join("x".join(map(str, k)) for k in tiled), flush=True)
    if IN_PROGRAM:                            # kernels inside a loop/chain do not pay the program's start
        program_setup = 0.0

    # measured = the XLA module; copy = its copy-start/done ops (an operand staged HBM -> on-chip memory before
    # the fusion, not overlapped with it); kernel = measured - copy.
    # NPUsim = its kernel + the program setup once (a one-kernel program; 0 with --in-program), against both the
    # measured module (err) and the measured kernel (err_k): a row without copies has the same prediction for
    # both. err_k_noprog is a SENSITIVITY only (the program setup dropped from the prediction): how much of the
    # 0.69 us chain intercept is staging and how much other device-side start-up has not been measured.
    out, errs = [], {}
    print(f"\n{'M':>5} {'K':>6} {'N':>6} | {'measured':>9} {'copy':>7} {'kernel':>8} {'NPUsim':>9} {'roofline':>9} us"
          f" | {'err':>7} {'err_k':>7}  src")
    for r in gemm:
        M, K, N = int(r["M"]), int(r["K"]), int(r["N"])
        src = "device" if num(r, "device_us") is not None else "wall"
        meas = num(r, "device_us") if src == "device" else num(r, "wall_us")
        ops = dict(kv.rsplit("=", 1) for kv in (r.get("ops") or "").split(";") if "=" in kv)
        copy = sum(float(v) for k, v in ops.items() if "copy-" in k) if src == "device" else 0.0
        if num(r, "copy_us") is not None:     # smallk: the module's copy-start/done time recorded directly
            copy = num(r, "copy_us")
        kern = meas - copy
        sk, roof = sim[M, K, N] * 1e6, roofline(M, K, N) * 1e6
        s = sk + program_setup * 1e6
        err, errk, errk_noprog = 100 * (s / meas - 1), 100 * (s / kern - 1), 100 * (sk / kern - 1)
        grp = "decode (M<=128)" if M <= 128 else "other"
        errs.setdefault(grp, []).append(abs(err))
        errs.setdefault(grp + " kernel", []).append(abs(errk))
        out.append(dict(M=M, K=K, N=N, measured_us=meas, copy_us=round(copy, 2), kernel_us=round(kern, 2), src=src,
                        npusim_us=round(s, 2), roofline_us=round(roof, 2), err_pct=round(err, 1),
                        err_kernel_pct=round(errk, 1), npusim_noprog_us=round(sk, 2),
                        err_kernel_noprog_pct=round(errk_noprog, 1), timing_kind=kind))
        print(f"{M:5d} {K:6d} {N:6d} | {meas:9.2f} {copy:7.2f} {kern:8.2f} {s:9.2f} {roof:9.2f}    | {err:+6.1f}%"
              f" {errk:+6.1f}%  {src}")
    print("MAPE " + ", ".join(f"{k} {statistics.mean(v):.1f}% (n={len(v)})" for k, v in errs.items()))
    tag = ("_" + a.accel if a.accel != "tpu_v5e" else "") + \
        ("_" + os.path.splitext(os.path.basename(a.dram))[0] if a.dram else "") + ("_chain" if a.chain else "")
    with open(os.path.join(a.results, f"compare_{os.path.splitext(a.gemm)[0]}{tag}.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, list(out[0]))
        w.writeheader()
        w.writerows(out)

    if not decode:
        return
    if all(num(r, "device_us") is None for r in decode):
        print("\ndecode: no device trace (wall time only), per-scope breakdown unavailable")
        return
    print(f"\n{'tp':>2} {'B':>4} {'L':>5} | {'total':>8} {'GEMM':>8} {'unattr':>7} {'attn':>7} {'KVroof':>7} {'psum':>7}"
          f" {'other':>6} {'nonGEMM':>7} | {'NPUsim GEMM':>11} {'err':>7}   (us per layer)")
    for r in decode:
        tp, B, L, tot = int(r["tp"]), int(r["batch"]), int(r["context"]), num(r, "device_us")
        sc = {k[:-3]: float(v) for k, v in r.items() if k.endswith("_us") and k not in ("wall_us", "device_us") and v}
        g, un = sum(sc.get(s, 0) for s in PROJ), sc.get("-", 0)
        ps = sc.get("psum_attn", 0) + sc.get("psum_ffn", 0)
        att = sc.get("attn", 0)
        kvroof = B * L * HKV * DH * 2 * 2 / tp / HBM * 1e6
        line = (f"{tp:2d} {B:4d} {L:5d} | {tot:8.1f} {g:8.1f} {un:7.1f} {att:7.1f} {kvroof:7.1f} {ps:7.1f}"
                f" {tot - g - un - att - ps:6.1f} {100 * (1 - (g + un) / tot):6.1f}% |")
        if tp == 1:
            # one program per layer; the measured side is GEMM scopes + unattributed ops, so where the
            # program's fixed time lands in the trace decides whether this constant belongs here
            sg = (sum(sim[(B, *kn)] for kn in PROJ.values()) + program_setup) * 1e6
            line += f" {sg:11.1f} {100 * (sg / (g + un) - 1) if g else float('nan'):+6.1f}%"
        print(line)
    print("GEMM = qkv + o_proj + ffn_up + ffn_down scopes; unattr = ops outside any scope, mostly the weight copies\n"
          "(copy-done) of the GEMMs, so NPUsim GEMM is compared with GEMM + unattr; other = norms, SwiGLU, residual")


if __name__ == "__main__":
    main()
