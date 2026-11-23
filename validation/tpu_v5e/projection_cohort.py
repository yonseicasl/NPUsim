#!/usr/bin/env python3
"""Manifest and audit of the projection cohort of a stream-bench profile run (2026-11-04): the cold loop GEMMs of
the four Llama-3-8B projections, bf16 and int8 weights. No simulator. For every case of the full grid it records
provenance (HLO / trace file hashes), the compiler's window config, three separately scoped measurements and
three audits, and says whether the case can be compared with a whole-GEMM NPUsim run under fused_output traffic:

  compiler accounting  (static, from the op metadata in xplane.pb; NOT measured bus traffic): the loop's fused op
                       reads weight bytes + M*K*2 from memory space 1 and writes 4 bytes
  event audit          (xplane.pb device events, timestamps): every convolution-category event lies inside exactly
                       one "XLA Modules" interval and every module holds exactly `iters` of them
  reference time       conv_event_us = mean duration of those events (only when the event audit passes);
                       loop_slope_us = the bench's wall slope from results.json -- another scope, not mixed in

Calibration overlap is recorded per case: this run is where the per-window cost was estimated (int8, M = 128) and
where the window schedule was developed, so the cohort is a RETROSPECTIVE evaluation, not a holdout.

Usage: python3 projection_cohort.py PROFILE_RUN_DIR OUT_DIR
"""
import csv
import hashlib
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trace_bytes import breakdown, fields

PROJ = {"qkv": (4096, 6144), "o_proj": (4096, 4096), "ffn_up": (4096, 28672), "ffn_down": (14336, 4096)}
MS, ITERS = (1, 8, 32, 128, 512), 8
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()


def window_config(text):
    """[(iteration_bounds, kernel_window_bounds, input_window_bounds)] of the convolutions in an HLO text."""
    keys = ("iteration_bounds", "kernel_window_bounds", "input_window_bounds")
    found = [re.findall(rf'"{k}":\[([^\]]*)\]', text) for k in keys]
    ints = lambda t: tuple(int(x.strip('" ')) for x in t.split(",") if x.strip('" '))
    return [tuple(ints(f[i]) if i < len(f) else () for f in found) for i in range(len(found[0]))]


def device_plane(path):
    """(event metadata {id: (name, stats)}, {line name: [(start_ps, end_ps, metadata id)]}) of the TPU plane."""
    for f, _, plane in fields(open(path, "rb").read()):
        if f != 1:
            continue
        p = list(fields(plane))
        if "TPU" not in next((v.decode() for ff, _, v in p if ff == 2), ""):
            continue
        names = {}
        for ff, _, v in p:
            if ff == 5:
                e = dict((a, b) for a, _, b in fields(v))
                names[e[1]] = dict((a, b) for a, _, b in fields(e[2])).get(2, b"").decode()
        meta, lines = {}, {}
        for ff, _, v in p:
            if ff == 4:
                e = dict((a, b) for a, _, b in fields(v))
                md = list(fields(e[2]))
                stats = {}
                for a, _, b in md:
                    if a == 5:
                        d = dict((x, y) for x, _, y in fields(b))
                        stats[names.get(d.get(1), d.get(1))] = next((d[k] for k in (3, 4, 5, 6) if k in d), None)
                meta[e[1]] = (next((b.decode() for a, _, b in md if a == 2), ""), stats)
            elif ff == 3:
                ln = list(fields(v))
                t0 = next((b for a, _, b in ln if a == 3), 0) * 1000          # timestamp_ns -> ps
                evs = []
                for a, _, b in ln:
                    if a == 4:
                        d = dict((x, y) for x, _, y in fields(b))
                        evs.append((t0 + d.get(2, 0), t0 + d.get(2, 0) + d.get(3, 0), d.get(1)) if 3 in d else None)
                lines[next((b.decode() for a, _, b in ln if a == 2), "?")] = evs
        return meta, lines
    return {}, {}


def audit(trace_dir, M, K, N, wbytes):
    """Compiler accounting and event audit of one capture."""
    pbs = sorted(os.path.join(dp, f) for dp, _, fs in os.walk(trace_dir) for f in fs if f.endswith(".xplane.pb"))
    files = {os.path.relpath(os.path.join(dp, f), trace_dir): sha(os.path.join(dp, f))
             for dp, _, fs in os.walk(trace_dir) for f in fs}
    out = dict(trace_files=files, why=[], modules=0, conv_events=0, per_module_counts=[], conv_event_us=None,
               read_space1=None, write_space1=None, other_space=None, fused_ops=[])
    if len(pbs) != 1:
        out["why"].append(f"{len(pbs)} xplane.pb files under the capture directory")
        return out
    meta, lines = device_plane(pbs[0])
    is_conv = lambda mid: "convolution" in str(meta.get(mid, ("", {}))[1].get("hlo_category", b""))
    conv_ids = sorted({mid for mid in meta if is_conv(mid)})
    out["fused_ops"] = [re.match(r"%?([\w.\-]+)", meta[i][0]).group(1) for i in conv_ids]
    if len(conv_ids) != 1:
        out["why"].append(f"{len(conv_ids)} convolution-category ops in the program, expected 1")
    else:
        bd = meta[conv_ids[0]][1].get("memory_access_breakdown")
        bd = breakdown(bd) if isinstance(bd, (bytes, bytearray)) else []
        acc = lambda op, space: sum(b for o, s, b in bd if o == op and s == space)
        out.update(read_space1=acc(1, 1), write_space1=acc(2, 1),
                   other_space=sum(b for o, s, b in bd if s != 1))
        want = K * N * wbytes + M * K * 2
        if not (want <= out["read_space1"] <= want + 16 and out["write_space1"] <= 8 and not out["other_space"]):
            out["why"].append(f"compiler accounting is not weights + M*K*2 read / <= 8 B written: read "
                              f"{out['read_space1']} (want {want}), write {out['write_space1']}, other spaces {out['other_space']}")
    mods = [e for name, evs in lines.items() if "Module" in name for e in evs]
    convs = [e for name, evs in lines.items() if "Ops" in name and "Async" not in name for e in evs
             if e is None or is_conv(e[2])]
    if any(m is None for m in mods) or not mods:
        out["why"].append("module events missing or without duration")
        return out
    counts, stray = [0] * len(mods), 0
    for c in convs:
        inside = [i for i, m in enumerate(mods) if c is not None and c[0] >= m[0] and c[1] <= m[1]]
        if len(inside) == 1:
            counts[inside[0]] += 1
        else:
            stray += 1
    out.update(modules=len(mods), conv_events=len(convs), per_module_counts=sorted(set(counts)))
    if stray:
        out["why"].append(f"{stray} convolution events not inside exactly one module interval")
    if set(counts) != {ITERS}:
        out["why"].append(f"modules with a convolution-event count other than {ITERS}: {sorted(set(counts))}")
    if not out["why"]:
        out["conv_event_us"] = sum(c[1] - c[0] for c in convs) / len(convs) / 1e6
    return out


def main(run, outdir):
    res = json.load(open(os.path.join(run, "results.json")))
    slope = {(r["M"], r["K"], r["N"], r["weights"]): r["us"] for r in res["stream_gemm"] if r["cold"]}
    rows = []
    for proj, (K, N) in PROJ.items():
        for M in MS:
            for dt, wb in (("bf16", 2), ("int8", 1)):
                tag = f"gemm_{M}x{K}x{N}_{dt}_cold"
                note = ""
                if not os.path.isdir(os.path.join(run, "traces", tag)):       # the bench's one specially tagged case
                    alt = f"stream_gemm_{M}x{K}x{N}_{dt}"
                    if os.path.isdir(os.path.join(run, "traces", alt)):
                        tag, note = alt, "captured by the bench under its special tag stream_gemm_* (not gemm_*_cold)"
                hlo = os.path.join(run, "hlo", tag + ".txt")
                row = dict(proj=proj, dtype=dt, M=M, K=K, N=N, tag=tag, note=note, loop_slope_us=slope.get((M, K, N, dt)))
                why = []
                if os.path.exists(hlo):
                    wc = window_config(open(hlo).read())
                    row.update(hlo=os.path.relpath(hlo, run), hlo_sha256=sha(hlo), convolutions_in_hlo=len(wc))
                    if len(wc) == 1 and len(wc[0][0]) == 3 and len(wc[0][1]) == 2:
                        it, kw, iw = wc[0]
                        n = it[0] * it[1] * it[2]
                        row.update(windows=n, iteration_bounds="x".join(map(str, it)),
                                   kernel_window_bounds="x".join(map(str, kw)), input_window_bounds="x".join(map(str, iw)),
                                   window_coverage=round(n * kw[0] * 8 * kw[1] * 128 / (K * N), 4))
                    else:
                        why.append(f"window_config not one record of 3 iteration / 2 kernel bounds: {wc}")
                else:
                    why.append("no HLO text")
                tdir = os.path.join(run, "traces", tag)
                if os.path.isdir(tdir):
                    a = audit(tdir, M, K, N, wb)
                    why += a.pop("why")
                    row.update({k: (json.dumps(v) if isinstance(v, (dict, list)) else v) for k, v in a.items()})
                else:
                    why.append("no trace")
                row.update(w_calibration_case=(dt == "int8" and M == 128),
                           inspected_before=True,      # all of this run fed int8_decompose.py (tables B, C)
                           supported=not why, exclusion="; ".join(why))
                rows.append(row)
    os.makedirs(outdir, exist_ok=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(outdir, "manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    src = {f: sha(os.path.join(run, f)) for f in ("results.json", "gemm.csv", "xla_windows.csv", "env.json")
           if os.path.exists(os.path.join(run, f))}
    zipf = run.rstrip("/") + ".zip"
    json.dump(dict(run=os.path.abspath(run), env=res.get("env"), source_files=src,
                   raw_archive={os.path.basename(zipf): sha(zipf)} if os.path.exists(zipf) else None,
                   profiled_loop_iterations=ITERS,
                   scopes=dict(conv_event_us="mean duration of the loop's convolution-category device events (xplane.pb)",
                               loop_slope_us="bench wall slope (median wall(2n) - median wall(n)) / n, results.json",
                               read_space1="compiler static matrix-payload accounting of the fused op, not measured traffic"),
                   calibration="RETROSPECTIVE cohort: window_switch_cycle 590 (0.393 us) was estimated on this run's "
                               "four int8 M = 128 cases; the window schedule was developed on this run; R 752 GB/s and "
                               "F 2 us of the closed-form analyses came from follow-up 2 / 3 (same shapes, other "
                               "sessions) and are not used by the simulator",
                   not_modelled="the fused op's input add, output conversion and reduction, its 4-byte result"),
              open(os.path.join(outdir, "manifest.json"), "w"), indent=1)
    ok = [r for r in rows if r["supported"]]
    print(f"{len(rows)} grid cases, {len(ok)} supported, {sum(r['w_calibration_case'] for r in ok)} of them w-calibration cases")
    for r in rows:
        if not r["supported"]:
            print(f"  excluded {r['tag']}: {r['exclusion']}")
        elif r["note"]:
            print(f"  note {r['tag']}: {r['note']}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
