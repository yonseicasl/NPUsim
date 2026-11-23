#!/usr/bin/env python3
"""Next joint sweep, GENERATE ONLY: the screening grid of memory x TPU points with their derived quantities and
the accelerator configs of the non-baseline TPU points. Nothing is simulated.

Memory points are EXISTING MOSAIC yaml; their coupled quantities are derived, not chosen independently
(dse_pilot.memory_facts): total DQ = channels x DQ per channel, capacity = channels x organisation, nominal
interface rate = total DQ x rate. Two families that must not be mixed up:
  hbm2e_refsb            HBM stacks beside the die (the v5e baseline; 16 GiB)
  HBDRAM (on logic)      the repository's hybrid-bonded profiles, all 8 GiB:
      pin-count axis     cloud_{0.25, 0.5, 1, 2, 4}TBps -- capacity constant, channels x 1024 pins from 4 to 64:
                         RESOURCE-INCREASING (bond pins / interface width)
      organisation axis  hb_dram_16ch vs cloud_1TBps -- same DQ, capacity and rate, 128 vs 512 columns per row:
                         ISO-RESOURCE
  Stack depth is not a memory point: the performance yaml does not carry it (it lives in the thermal spec).
TPU points (accelerator cfg generated here; workloads must be REMAPPED per point at run time -- the mapping
generator takes the MXU count, a different peak-MAC number alone would not be an array sweep):
  mxu4       tpu_v5e.cfg as is (baseline)
  mxu8       num_chips and [multi_chip] height 8: RESOURCE-INCREASING (2 x MXUs)
  mxu8_glb2  mxu8 with the [separate] GLB input / weight / output sizes doubled: RESOURCE-INCREASING (SRAM)
  A 256 x 256 array is not generated: it needs new array geometry and mappings, not a text edit.
Clock: every point keeps [dram] frequency = compute clock = 1500 MHz.

Usage: python3 dse_manifest.py MOSAIC_ROOT OUT_DIR
"""
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dse_pilot as dp

MEM = [("hbm2e_refsb", "configs/architecture/dram/hbm2e_tpu_v5e_refsb.yaml", "HBM beside the die", "baseline"),
       ("hbdram_p4", "configs/architecture/dram/cloud/cloud_0.25TBps.yaml", "HBDRAM on logic", "pin-count axis"),
       ("hbdram_p8", "configs/architecture/dram/cloud/cloud_0.5TBps.yaml", "HBDRAM on logic", "pin-count axis"),
       ("hbdram_p16", "configs/architecture/dram/cloud/cloud_1TBps.yaml", "HBDRAM on logic", "pin-count axis (= the pilot's example profile)"),
       ("hbdram_p32", "configs/architecture/dram/cloud/cloud_2TBps.yaml", "HBDRAM on logic", "pin-count axis"),
       ("hbdram_p64", "configs/architecture/dram/cloud/cloud_4TBps.yaml", "HBDRAM on logic", "pin-count axis"),
       ("hbdram_p16_col128", "configs/architecture/dram/hb_dram_16ch.yaml", "HBDRAM on logic", "organisation axis (iso-resource with hbdram_p16)")]
GRID = [(m[0], "mxu4") for m in MEM] + [(m, "mxu8") for m in ("hbm2e_refsb", "hbdram_p16", "hbdram_p32", "hbdram_p64")] + \
       [("hbdram_p64", "mxu8_glb2")]


def tpu_cfg(base, name):
    """(cfg text, facts) of a TPU point derived from tpu_v5e.cfg."""
    mxus, glb = (8 if name.startswith("mxu8") else 4), (2 if name.endswith("glb2") else 1)
    cfg = base
    if mxus != 4:
        cfg, a = re.subn(r"(?m)^num_chips\s*=\s*4\s*$", f"num_chips={mxus}", cfg)
        sec = cfg.index("[multi_chip]")
        tail, b = re.subn(r"(?m)^height\s*=\s*4\s*$", f"height = {mxus}", cfg[sec:], count=1)
        assert a == 1 and b == 1
        cfg = cfg[:sec] + tail
    size = {}
    lo, hi = cfg.index("[separate]"), cfg.index("[multi_chip]")       # the GLB section; the PE buffers keep theirs
    sep = cfg[lo:hi]
    for k in ("input_size", "weight_size", "output_size"):
        v = int(re.search(rf"(?m)^{k}\s*=\s*(\d+)", sep).group(1)) * glb
        sep, n = re.subn(rf"(?m)^{k}\s*=\s*\d+", f"{k} = {v}", sep)
        assert n == 1
        size[k + "_kb"] = v
    cfg = cfg[:lo] + sep + cfg[hi:]
    return cfg, dict(tpu=name, mxus=mxus, array="128 x 128", peak_macs_per_cycle=mxus * 128 * 128,
                     strip_columns=mxus * 128, **size,
                     resource="baseline" if name == "mxu4" else "resource-increasing: " + ("MXUs" if glb == 1 else "MXUs and GLB SRAM"))


def main():
    root, out = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    if os.path.exists(os.path.join(out, "manifest.csv")):
        raise SystemExit(f"{out} already holds a manifest; nothing was changed")
    os.makedirs(os.path.join(out, "accelerators"), exist_ok=True)
    base = open(os.path.join(dp.ca.ROOT, "configs", "accelerators", f"{dp.ca.ACCEL}.cfg")).read()
    mem = {n: dict(dp.memory_facts(root, os.path.join(root, y)), family=fam, axis=axis) for n, y, fam, axis in MEM}
    tpu = {}
    for name in dict.fromkeys(t for _, t in GRID):
        cfg, tpu[name] = tpu_cfg(base, name)
        if name != "mxu4":
            open(os.path.join(out, "accelerators", f"tpu_v5e_{name}.cfg"), "w").write(cfg)
            tpu[name]["config"] = f"accelerators/tpu_v5e_{name}.cfg"
            tpu[name]["config_sha256"] = dp.sha(os.path.join(out, "accelerators", f"tpu_v5e_{name}.cfg"))
        else:
            tpu[name]["config"] = f"configs/accelerators/{dp.ca.ACCEL}.cfg (as is)"
    ref = mem["hbm2e_refsb"]
    rows = []
    for m, t in GRID:
        f, c = mem[m], tpu[t]
        w = dp.WORKLOADS["gemm_m512"]
        rows.append(dict(
            point=f"{m}__{t}", memory=m, family=f["family"], axis=f["axis"], yaml=f["yaml"], impl=f["impl"],
            channels=f["channels"], dq_per_channel=f["dq_per_channel"], total_dq=f["total_dq"], rate_mbps_per_dq=f["rate_mbps_per_dq"],
            capacity_gib=f["capacity_bytes"] / dp.GIB, nominal_interface_gbps=round(f["nominal_interface_gbps"], 1),
            transaction_bytes=f["transaction_bytes"], tck_ns=f["tck_ns"],
            dq_vs_baseline=f["total_dq"] / ref["total_dq"], capacity_vs_baseline=f["capacity_bytes"] / ref["capacity_bytes"],
            tpu=t, mxus=c["mxus"], peak_macs_per_cycle=c["peak_macs_per_cycle"], strip_columns=c["strip_columns"],
            glb_input_kb=c["input_size_kb"], glb_weight_kb=c["weight_size_kb"], glb_output_kb=c["output_size_kb"],
            tpu_resource=c["resource"],
            # the stream rate one strip of a K 4096 bf16 GEMM at M 512 would need to keep the MXUs busy
            gemm_m512_weight_gbps_to_match_compute=round(w["K"] * c["strip_columns"] * dp.ACT /
                                                         (w["M"] * w["K"] * c["strip_columns"] / c["peak_macs_per_cycle"] / dp.ca.CLOCK_HZ) / 1e9, 1),
            long_stream_rule="checked for hbm2e_refsb and hbdram_p16 only; run unittest/run_mosaic_kv_stream_test.sh "
                             "with this yaml before any extrapolated stream" if m not in ("hbm2e_refsb", "hbdram_p16") else "checked",
            status="generated; not simulated"))
    with open(os.path.join(out, "manifest.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    json.dump(dict(memories=mem, tpu_points=tpu, workloads=dp.WORKLOADS, grid=rows,
                   not_generated=["256 x 256 arrays (new geometry and mappings)", "HBM stack / channel expansion (needs a new yaml and "
                                  "a statement of what is scaled; architecture question open)", "stack depth as a performance axis "
                                  "(not in the performance yaml)", "DVFS / mixed clocks (cycle-domain conversion not verified)"]),
              open(os.path.join(out, "manifest.json"), "w"), indent=1)
    print(f"{len(rows)} points x {len(dp.WORKLOADS)} workloads = {len(rows) * len(dp.WORKLOADS)} simulations if run; nothing simulated")
    for r in rows:
        print(f"  {r['point']:32s} {r['channels']:3d} ch x {r['dq_per_channel']:5d} DQ = {r['total_dq']:6d} DQ, {r['capacity_gib']:4.0f} GiB, "
              f"{r['nominal_interface_gbps']:7.1f} GB/s nominal, tx {r['transaction_bytes']:4d} B | {r['mxus']} MXUs, GLB weight {r['glb_weight_kb']} KB | {r['tpu_resource']}")


if __name__ == "__main__":
    main()
