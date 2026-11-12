#!/usr/bin/env python3
"""Where the int8-weight GEMM error is, from stored files only (no simulator, no TPU).

A  loop-slope scope (follow-up 2 / 3 stream bench: wall slope of one loop iteration, loop body included):
   t = F + weight_bytes / r over the four projections, per (run, M, dtype), M <= 32 -- TPU and NPUsim.
B  device-trace scope (2026-11-04 profile run: self time of the GEMM fusion): fusion / max(stream, compute)
   against compute / stream for every cold multi-window GEMM.
C  a HYPOTHESIS check, not a model: a compute-side cost w per XLA window, estimated on the int8 M = 128 cases
   only and applied to all other multi-window cases. Every constant here (752 GB/s, F = 2 us, w) is fitted to
   TPU data; none of them is in the simulator.

Usage: python3 int8_decompose.py RESULTS_DIR PROFILE_RUN_DIR
  RESULTS_DIR      validation/tpu_v5e/results (follow-up 2 / 3 and their compare_gemm*.csv)
  PROFILE_RUN_DIR  .../results/2026-11-04_colab_v5e1_profile (branch mosaic: xla_windows.csv, xprof/hlo_stats.csv)
"""
import collections
import csv
import json
import os
import statistics as st
import sys

PEAK = 4 * 128 * 128 * 2 * 1.5e9          # FLOP/s, 4 MXUs
R, F = 752e9, 2.0                          # B/s and us: the A fit's rate and intercept, used by B and C


def fit(pts):
    """Least squares t = F + x * s over (x, t); returns (F, s)."""
    n, sx, sy = len(pts), sum(x for x, _ in pts), sum(t for _, t in pts)
    s = (n * sum(x * t for x, t in pts) - sx * sy) / (n * sum(x * x for x, _ in pts) - sx * sx)
    return (sy - s * sx) / n, s


def pipeline(n, S, C, w):
    """n windows, each DMAed whole before its compute, the next one fetched meanwhile; w = compute-side cost
    per window."""
    d, c = S / n, C / n
    return F + d + (n - 1) * max(d, c + w) + c + w


def main(results, profile):
    print("A  loop-slope scope, M <= 32: t = F + weight_bytes / r   (F us, r GB/s)")
    print(f"{'':22s} {'M':>3} | {'bf16 F':>7} {'r':>5} | {'int8 F':>7} {'r':>5}")
    for run in ("2026-11-03_colab_v5e1_followup2", "2026-11-03_colab_v5e1_followup3"):
        d = json.load(open(os.path.join(results, run, "results.json")))
        t = {(r["M"], r["K"], r["N"], r["weights"]): r["us"] for r in d["stream_gemm"] if r["cold"]}
        for M in (1, 8, 32):
            cells = []
            for w, bpe in (("bf16", 2), ("int8", 1)):
                f0, s = fit([(K * N * bpe / 1e6, us) for (m, K, N, ww), us in t.items() if m == M and ww == w])
                cells.append(f"{f0:7.2f} {1e3 / s:5.0f}")
            print(f"TPU {run[-9:]:18s} {M:3d} | " + " | ".join(cells))
    f2 = os.path.join(results, "2026-11-03_colab_v5e1_followup2")
    for tag, name in (("_hbm2e_tpu_v5e_refsb", "refsb"), ("", "analytic")):
        for M in (1, 8, 32):
            cells = []
            for f, bpe in ((f"compare_gemm{tag}.csv", 2), (f"compare_gemm_int8_tpu_v5e_int8w{tag}.csv", 1)):
                f0, s = fit([(int(r["K"]) * int(r["N"]) * bpe / 1e6, float(r["npusim_us"]))
                             for r in csv.DictReader(open(os.path.join(f2, f)))
                             if int(r["M"]) == M and int(r["K"]) in (4096, 14336)])
                cells.append(f"{f0:7.2f} {1e3 / s:5.0f}")
            print(f"NPUsim {name:15s} {M:3d} | " + " | ".join(cells))

    win = {r["case"]: r for r in csv.DictReader(open(os.path.join(profile, "xla_windows.csv")))}
    ops = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(profile, "xprof", "hlo_stats.csv"))):
        if "fusion" in r["hlo_op_name"] or "convolution" in r["category"]:
            ops[r["case"]].append(float(r["avg_self_time"]))
    rows = []                              # (name, M, n windows, int8?, fusion us, stream us, compute us)
    for case, x in win.items():
        if case.endswith("_cold") and int(x["windows"]) >= 2 and ops[case]:
            M, K, N, i8 = int(x["M"]), int(x["K"]), int(x["N"]), "_int8_" in case
            rows.append((case[5:-5], M, int(x["windows"]), i8, max(ops[case]),
                         K * N * (1 if i8 else 2) / R * 1e6, 2 * M * K * N / PEAK * 1e6))
    print(f"\nB  device-trace scope: {len(rows)} cold multi-window GEMM fusions; S = weight bytes / {R / 1e9:.0f} GB/s,"
          " C = 2MKN / peak")
    print(f"{'case':26s} {'windows':>7} {'fusion':>8} {'S':>7} {'C':>7} {'C/S':>5} {'fusion/max':>10}")
    for name, M, n, i8, f, S, C in sorted(rows, key=lambda r: r[6] / r[5]):
        if 0.7 <= C / S <= 1.4 or (i8 and M >= 128) or (M == 512 and n >= 16):
            print(f"{name:26s} {n:7d} {f:8.2f} {S:7.2f} {C:7.2f} {C / S:5.2f} {f / max(S, C):10.3f}")

    train = [r for r in rows if r[3] and r[1] == 128]
    ws = [(f - F - S / n) / n - C / n for _, _, n, _, f, S, C in train]
    w = st.mean(ws)
    print(f"\nC  hypothesis: compute-side cost per window. w from the {len(train)} int8 M = 128 cases: "
          + ", ".join(f"{x:.3f}" for x in ws) + f" -> {w:.3f} us")
    big = [(name, (f - C - S / n) / n) for name, M, n, i8, f, S, C in rows if M == 512 and n >= 56]
    print("   the same quantity on compute-bound cases with >= 56 windows (F neglected): "
          + ", ".join(f"{name} {x:.3f}" for name, x in big))
    groups = collections.defaultdict(list)
    for name, M, n, i8, f, S, C in rows:
        if i8 and M == 128:
            continue
        g = ("int8" if i8 else "bf16") + (" C/S < 0.4" if C / S < 0.4 else " C/S 0.4-1.4" if C / S < 1.4 else " C/S > 1.4")
        groups[g].append((pipeline(n, S, C, 0.0) / f - 1, pipeline(n, S, C, w) / f - 1))
    groups["all held out"] = [e for v in groups.values() for e in v]
    print(f"   {'held-out group':18s} {'n':>3} | w = 0: mean / MAPE | w = {w:.3f}: mean / MAPE")
    for g, v in sorted(groups.items()):
        print(f"   {g:18s} {len(v):3d} | {100 * st.mean(e[0] for e in v):+6.1f} / {100 * st.mean(abs(e[0]) for e in v):4.1f}%"
              f"     | {100 * st.mean(e[1] for e in v):+6.1f} / {100 * st.mean(abs(e[1]) for e in v):4.1f}%")
    print("   w is not one number (0.17-0.41 us by context) and nothing shows what happens at a window boundary:"
          " not adopted.")


if __name__ == "__main__":
    assert abs(fit([(1, 5.0), (2, 7.0), (4, 11.0)])[0] - 3.0) < 1e-9 and abs(fit([(1, 5.0), (2, 7.0)])[1] - 2.0) < 1e-9
    assert abs(pipeline(4, 8.0, 4.0, 0.0) - (F + 2.0 + 3 * 2.0 + 1.0)) < 1e-9     # DMA-bound
    assert abs(pipeline(4, 4.0, 8.0, 0.5) - (F + 1.0 + 3 * 2.5 + 2.5)) < 1e-9     # compute-bound
    main(*sys.argv[1:3])
