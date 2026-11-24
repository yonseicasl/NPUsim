#!/usr/bin/env python3
"""Consistency report of the fixed target profile (target_profile.json): per-die vs aggregate, geometry, derived
rates, and whether MOSAIC's DRAM model can represent the organisation as given. Nothing is simulated beyond three
short runtime probes of STAND-IN yaml, which are not the target (their capacity is wrong by construction).

Usage: python3 target_profile.py MOSAIC_ROOT OUT_DIR [--stream-probe PATH/hbm_stream_probe]
"""
import argparse
import json
import os
import re
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mosaic_root")
    ap.add_argument("out")
    ap.add_argument("--stream-probe", help="validation/tpu_v5e/hbm_stream_probe.cc built (see its header); optional")
    a = ap.parse_args()
    root, out = os.path.abspath(a.mosaic_root), os.path.abspath(a.out)
    if os.path.exists(os.path.join(out, "consistency_report.json")):
        raise SystemExit(f"{out} already holds a report; nothing was changed")
    os.makedirs(os.path.join(out, "standin_yaml"), exist_ok=True)
    p = json.load(open(os.path.join(HERE, "target_profile.json")))
    die, n = p["dram_die"], p["topology"]["dram_dies"]
    ch, banks = die["channels"], die["banks_per_channel"]
    dq = die["dq_bits_per_channel"]
    cap_b, bw = die["capacity"]["value"] * 1e9, die["nominal_peak_bandwidth"]["value"] * 1e12       # decimal, as stated
    fx, fy = die["footprint"]["x"], die["footprint"]["y"]
    lx, ly = p["logic_die"]["footprint"]["x"], p["logic_die"]["footprint"]["y"]
    agg = p["aggregate_as_stated_by_user"]
    d = dict(
        per_die=dict(channels=ch, banks=ch * banks, capacity_bytes=cap_b, nominal_peak_bytes_per_s=bw, footprint_mm2=fx * fy,
                     capacity_bits_per_channel=cap_b * 8 / ch, capacity_bits_per_bank=cap_b * 8 / ch / banks,
                     nominal_bytes_per_s_per_channel=bw / ch, nominal_bits_per_s_per_channel=bw * 8 / ch,
                     bit_density_Gb_per_mm2=cap_b * 8 / 1e9 / (fx * fy), seconds_to_read_whole_die_at_nominal=cap_b / bw),
        aggregate=dict(channels=n * ch, banks=n * ch * banks, capacity_bytes=n * cap_b, nominal_peak_bytes_per_s=n * bw,
                       dram_footprint_mm2=n * fx * fy),
        geometry=dict(logic_footprint_mm2=lx * ly, envelope_2x2_mm=[2 * fx, 2 * fy], envelope_fits_logic=2 * fx <= lx and 2 * fy <= ly,
                      slack_mm=[lx - 2 * fx, ly - 2 * fy], dram_over_logic_area_ratio=n * fx * fy / (lx * ly),
                      note="geometric envelope without gaps or keep-outs; says nothing about bonding, routing or thermal feasibility"))
    checks = [("aggregate channels as stated = dies x per die", agg["channels"] == n * ch),
              ("aggregate capacity as stated = dies x per die", agg["capacity_GB"] * 1e9 == n * cap_b),
              ("aggregate nominal bandwidth as stated = dies x per die", agg["nominal_peak_bandwidth_TB_per_s"] * 1e12 == n * bw),
              ("2 x 2 envelope fits the logic footprint (no gaps assumed)", d["geometry"]["envelope_fits_logic"])]
    # Fixed per-channel DQ; data rate is implied by the nominal peak, not a supplied clock.
    per_ch_bits = bw * 8 / ch
    d["interface_fixed_dq"] = dict(
        dq_bits_per_channel=dq,
        derived_required_rate_mbps_per_dq=per_ch_bits / dq / 1e6,
        data_dq_bits_per_die=ch * dq,
        data_dq_bits_aggregate=n * ch * dq,
        note="Data-interface width only, not total bond count. Rate derived from theoretical bandwidth; "
             "signaling mode/clock, protocol overhead and sustained efficiency remain unspecified.")
    # Representability in MOSAIC's HBDRAM model (python spec = source of the generated C++)
    spec = open(os.path.join(root, "simulator/src/dram/ramulator2/python/ramulator/dram/hbdram.py")).read()
    presets = {m.group(1): json.loads(m.group(2)) for m in re.finditer(r'"(HBDRAM_\w+)":\s*(\{"density"[^}]*\})', spec)}
    bank_counts = sorted({v["bankgroup"] * v["bank"] for v in presets.values()})
    densities = sorted({v["density"] for v in presets.values()})
    need_mb = cap_b * 8 / ch / 1e6
    rep = dict(
        channels=dict(target=ch, representable=None, how="runtime probe below"),
        banks_per_channel=dict(target=banks, representable=banks in bank_counts, model_has=bank_counts,
                               note="every HBDRAM preset has one bank per channel (stated in the spec's docstring); 4 banks need a new preset"),
        capacity_per_channel=dict(target_megabits_decimal=need_mb, preset_densities_megabits=densities,
                                  representable=False,
                                  note="187.5 Mb (decimal) per channel: with a power-of-two DQ the rows x columns per bank are not an integer; "
                                       "read as GiB it is 192 Mib = 3 x 2^26 bits, which needs a non-power-of-two row count and a density the "
                                       "model has no tRFC / tREFI entry for. NOT rounded to a preset here."),
        timings=dict(representable=False, note="target timings / refresh unknown; HBDRAM_500Mbps / 400Mbps are the repository's, not the target's"),
        four_domains=dict(representable=False, note="NPUsim has one [dram] endpoint and one staging link; four interface domains with local and "
                                                    "cross-domain routes are not modeled. A single 1024-channel endpoint is an IDEALIZED "
                                                    "aggregate-supply scenario only."))
    probes = []
    base = open(os.path.join(root, "configs/architecture/dram/thermal_sweep/hb3d_example.yaml")).read()
    for preset, chans, why in (("HBDRAM_2Gb_512pin", ch, "one die's channel count at 512 DQ x 500 Mbps (capacity 64 GB: NOT the target)"),
                               ("HBDRAM_2Gb_512pin", n * ch, "idealized single endpoint of all four dies (capacity 256 GB: NOT the target)"),
                               ("HBDRAM_1Gb_256pin", ch, "one die's channel count at 256 DQ (4.1 TB/s nominal: NOT the target rate)")):
        y = os.path.join(out, "standin_yaml", f"{preset}_x{chans}.yaml")
        t, k1 = re.subn(r"preset: HBDRAM_4Gb_1024pin_512col", f"preset: {preset}", base)
        t, k2 = re.subn(r"channel: 16", f"channel: {chans}", t)
        assert k1 == 1 and k2 == 1
        open(y, "w").write(t)
        r = subprocess.run([os.path.join(root, "simulator/build/bin/dram_bw_probe"), y, "1", "0"], capture_output=True, text=True).stdout
        m = re.search(r"channels=(\d+) tx=(\d+)B tCK=([\d.]+)ns", r)
        org = presets[preset]
        row = dict(standin=os.path.relpath(y, out), why=why, runs=bool(m), runtime_channels=m and int(m.group(1)),
                   transaction_bytes=m and int(m.group(2)), tck_ns=m and float(m.group(3)),
                   nominal_bytes_per_s=chans * org["dq"] * 500e6 / 8)
        if a.stream_probe and m:
            s = subprocess.run([a.stream_probe, y, str(64 << 20), str(64 * int(m.group(2))), "0"], capture_output=True, text=True).stdout
            g = re.search(r"([\d.]+) GB/s", s)
            row["sequential_read_GB_per_s_64MiB_64tx_requests"] = g and float(g.group(1))
        probes.append(row)
    rep["channels"]["representable"] = all(q["runs"] and q["runtime_channels"] in (ch, n * ch) for q in probes)
    d.update(checks=[dict(check=c, ok=bool(ok)) for c, ok in checks], representability=rep, standin_runtime_probes=probes,
             decode_bound_example=dict(
                 note="arithmetic only, no model selected: a batch-1 decode token of a model whose weights fill W bytes reads W once; "
                      "at the AGGREGATE THEORETICAL MAXIMUM that takes W / 32 TB/s. It is a lower bound on supply time, not a latency.",
                 examples={f"{w} GB of weights": f"{w * 1e9 / (n * bw) * 1e3:.3f} ms" for w in (8, 16, 24)}))
    json.dump(d, open(os.path.join(out, "consistency_report.json"), "w"), indent=1)
    pd, ag, ge = d["per_die"], d["aggregate"], d["geometry"]
    lines = [f"per die: {ch} channels x {banks} banks = {pd['banks']} banks; {cap_b / 1e9:g} GB; {bw / 1e12:g} TB/s nominal peak; {fx} x {fy} mm = {pd['footprint_mm2']} mm^2",
             f"   per channel: {pd['capacity_bits_per_channel'] / 1e6:g} Mb, {pd['nominal_bytes_per_s_per_channel'] / 1e9:g} GB/s = {pd['nominal_bits_per_s_per_channel'] / 1e9:g} Gb/s; per bank {pd['capacity_bits_per_bank'] / 1e6:g} Mb",
             f"   bit density {pd['bit_density_Gb_per_mm2']:.3f} Gb/mm^2; whole die read at the nominal peak in {pd['seconds_to_read_whole_die_at_nominal'] * 1e3:g} ms",
             f"aggregate (4 dies, 2 x 2, one tier): {ag['channels']} channels, {ag['banks']} banks, {ag['capacity_bytes'] / 1e9:g} GB, {ag['nominal_peak_bytes_per_s'] / 1e12:g} TB/s nominal peak, {ag['dram_footprint_mm2']} mm^2 of DRAM",
             f"geometry: logic {lx} x {ly} mm = {ge['logic_footprint_mm2']} mm^2; 2 x 2 envelope {ge['envelope_2x2_mm'][0]} x {ge['envelope_2x2_mm'][1]} mm (slack {ge['slack_mm'][0]} x {ge['slack_mm'][1]} mm, no gaps assumed); DRAM over {100 * ge['dram_over_logic_area_ratio']:.1f}% of the logic footprint",
             f"fixed interface: {dq} bits DQ/channel; derived required rate {d['interface_fixed_dq']['derived_required_rate_mbps_per_dq']:.4f} Mb/s per DQ (not a clock specification)",
             "representability in MOSAIC's HBDRAM model: " + ", ".join(f"{k} {'yes' if v['representable'] else 'NO'}" for k, v in rep.items())]
    lines += [f"   stand-in {q['standin']}: runs {q['runs']}, {q['runtime_channels']} channels, tx {q['transaction_bytes']} B, nominal {q['nominal_bytes_per_s'] / 1e12:.3f} TB/s"
              + (f", sequential read {q['sequential_read_GB_per_s_64MiB_64tx_requests']:.0f} GB/s" if q.get("sequential_read_GB_per_s_64MiB_64tx_requests") else "")
              + f"  [{q['why']}]" for q in probes]
    lines += [f"checks: {sum(c['ok'] for c in d['checks'])} of {len(d['checks'])} pass"]
    open(os.path.join(out, "consistency_report.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
