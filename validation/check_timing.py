#!/usr/bin/env python3
"""Acceptance gate for the committed Gemmini RTL and Eyeriss silicon baselines."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VALIDATION = ROOT / "validation"
RESULT = ROOT / "result"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as source:
        return list(csv.DictReader(source))


def metric(text: str, label: str, default: float | None = None) -> float:
    match = re.search(rf"{re.escape(label)}\s*:\s*([\d.]+)", text)
    if match is None:
        if default is not None:
            return default
        raise ValueError(f"missing {label!r} in simulator result")
    return float(match.group(1))


def compute_schedule_latency(text: str, source: str) -> float:
    """P0-4 identity gate: the official "Compute-schedule latency" metric must
    always equal "Computation cycle" + "Fold fill cycle" (the latter line is
    omitted from the printout when it is exactly zero). The tolerance accounts
    for both terms being independently rounded to one decimal at print time."""
    compute_schedule = metric(text, "Compute-schedule latency")
    computation = metric(text, "Computation cycle")
    fold_fill = metric(text, "Fold fill cycle", default=0.0)
    identity = computation + fold_fill
    tolerance = max(0.2, abs(identity) * 1e-6)
    if abs(compute_schedule - identity) > tolerance:
        raise AssertionError(
            f"{source}: Compute-schedule latency ({compute_schedule}) != "
            f"Computation cycle + Fold fill cycle ({identity})"
        )
    return compute_schedule


def busy_cycle_axes(text: str, source: str) -> list[tuple[float, float, float]]:
    """P0-2 network contract: parse the five (access / link / overlap) busy-cycle
    axis rows so callers can assert they were populated, not left at their
    zero-initialized default (the network-level regression this audit found)."""
    section = text.split("Busy-cycle axes", 1)[1].split("\n\n", 1)[0]
    axes = [
        tuple(float(value) for value in triple)
        for triple in re.findall(r":\s*([\d.]+)\s*/\s*([\d.]+)\s*/\s*([\d.]+)\s*cycles", section)
    ]
    if len(axes) != 5:
        raise AssertionError(f"{source}: expected 5 busy-cycle axis rows, found {len(axes)}")
    return axes


def stage_busy(text: str) -> list[float]:
    """The five per-stage busy values from the layer-timeline block."""
    timeline = text.split("Layer timeline", 1)[1].split("MAC result", 1)[0]
    busy = [float(v) for v in re.findall(r"\* [^:\n]+:\s*([\d.]+) cycles \(", timeline)]
    if len(busy) != 5:
        raise AssertionError(f"expected 5 stage busy rows, found {len(busy)}")
    return busy


def format_axis(text: str) -> float:
    timeline = text.split("Layer timeline", 1)[1].split("MAC result", 1)[0]
    return float(re.search(r"PE format-IP axis\s*:\s*([\d.]+)", timeline).group(1))


def check_network_timeline(path: Path) -> None:
    """P0-2/P0-4: the network-level rollup must carry the same self-consistent
    official metric and per-stage busy-cycle axes as the per-layer files it is
    aggregated from, not the zero-initialized defaults."""
    text = path.read_text(encoding="utf-8")
    compute_schedule = compute_schedule_latency(text, str(path))
    if compute_schedule <= 0.0:
        raise AssertionError(f"{path}: network Compute-schedule latency is non-positive ({compute_schedule})")
    axes = busy_cycle_axes(text, str(path))
    if all(value == 0.0 for row in axes for value in row):
        raise AssertionError(f"{path}: all network busy-cycle axes are zero")


def check_network_rollup(directory: Path) -> None:
    """L3: the network rollup contract. Network busy and network axes are TWO
    INDEPENDENT sums over the layers, and neither reconstructs the other:

        network busy[stage] == sum_layer(layer busy[stage])
        network axis[stage] == sum_layer(layer axis[stage])

    sum_layer(max(axes)) != max(sum_layer(axes)) in general, so the printed network axes
    must NOT be expected to reduce to network busy (the layer-scope T7 invariant in
    validation/traffic/check.py applies to layer files only). This check pins both sums so
    a future change cannot quietly redefine either side -- in particular, "fixing" network
    busy to max(network axes) would undercount serially executed layers and is caught here.
    The report states which contract each scope prints; assert that labelling too."""
    network_path = directory / "network.txt"
    network_text = network_path.read_text(encoding="utf-8")
    layer_paths = sorted(directory.glob("layer_*.txt"))
    if not layer_paths:
        raise AssertionError(f"{directory}: no layer_*.txt to roll up")

    if "[network: per-axis work summed over layers" not in network_text:
        raise AssertionError(f"{network_path}: network axis contract is not stated in the report")
    # L11/P4-14: the rollup must state its TIMING SCOPE next to the latency, and must say so
    # explicitly when layers were excluded. The count it declares has to match the number of
    # layer files that were actually folded in, or the scope statement is decorative.
    scope = re.search(r"Timing scope\s*:\s*(\d+) of (\d+) layers(.*)", network_text)
    if scope is None:
        raise AssertionError(f"{network_path}: network report does not state its timing scope")
    included, total, qualifier = int(scope.group(1)), int(scope.group(2)), scope.group(3)
    if included != len(layer_paths):
        raise AssertionError(
            f"{network_path}: declares {included} layers in scope but {len(layer_paths)} layer "
            f"result files were rolled up")
    if included < total and "PARTIAL" not in qualifier:
        raise AssertionError(
            f"{network_path}: {total - included} of {total} layers were excluded but the scope "
            f"line does not mark the rollup PARTIAL ({qualifier.strip()!r})")
    if included == total and "complete" not in qualifier:
        raise AssertionError(
            f"{network_path}: every layer is in scope but the report does not say so "
            f"({qualifier.strip()!r})")
    for layer_path in layer_paths:
        if "[layer: busy = max of these axes]" not in layer_path.read_text(encoding="utf-8"):
            raise AssertionError(f"{layer_path}: layer axis contract is not stated in the report")

    network_busy = stage_busy(network_text)
    network_axes = busy_cycle_axes(network_text, str(network_path))
    network_format = format_axis(network_text)
    summed_busy = [0.0]*5
    summed_axes = [[0.0]*3 for _ in range(5)]
    summed_format = 0.0
    for layer_path in layer_paths:
        layer_text = layer_path.read_text(encoding="utf-8")
        busy = stage_busy(layer_text)
        axes = busy_cycle_axes(layer_text, str(layer_path))
        summed_format += format_axis(layer_text)
        for stage in range(5):
            summed_busy[stage] += busy[stage]
            for axis in range(3):
                summed_axes[stage][axis] += axes[stage][axis]

    def close(measured: float, expected: float) -> bool:
        # Both sides are printed at one decimal, so the rollup of N layers can drift by
        # up to N*0.05 from the internally summed value.
        return abs(measured - expected) <= max(0.5, 0.05*len(layer_paths))

    for stage in range(5):
        if not close(network_busy[stage], summed_busy[stage]):
            raise AssertionError(
                f"{network_path}: network busy[{stage}] {network_busy[stage]} != "
                f"sum of layer busy {summed_busy[stage]}")
        for axis in range(3):
            if not close(network_axes[stage][axis], summed_axes[stage][axis]):
                raise AssertionError(
                    f"{network_path}: network axis[{stage}][{axis}] {network_axes[stage][axis]} "
                    f"!= sum of layer axis {summed_axes[stage][axis]}")
    if not close(network_format, summed_format):
        raise AssertionError(
            f"{network_path}: network format axis {network_format} != sum of layer "
            f"format axes {summed_format}")


def gemmini(check_baseline: bool) -> tuple[float, float]:
    baseline = {
        (int(row["m"]), int(row["k"]), int(row["n"])):
        (float(row["computation_cycles"]), float(row["fold_fill_cycles"]))
        for row in rows(VALIDATION / "phase2" / "npusim_baseline.csv")
    }
    errors: list[float] = []
    for row in rows(VALIDATION / "phase2" / "golden_rtl_cycles.csv"):
        point = int(row["m"]), int(row["k"]), int(row["n"])
        path = RESULT / "gemmini" / f"gemm_{point[0]}x{point[1]}x{point[2]}" / "ws" / "layer_0.txt"
        text = path.read_text(encoding="utf-8")
        measured = metric(text, "Computation cycle"), metric(text, "Fold fill cycle", default=0.0)
        if check_baseline and measured != baseline[point]:
            raise AssertionError(f"Gemmini NPUsim baseline changed at {point}: {measured} != {baseline[point]}")
        # P0-3: the golden comparison reads the official "Compute-schedule latency"
        # metric directly instead of re-deriving it from the two component fields.
        compute_schedule = compute_schedule_latency(text, str(path))
        check_network_timeline(path.with_name("network.txt"))
        check_network_rollup(path.parent)
        reference = float(row["rtl_cycles"])
        errors.append(abs((compute_schedule - reference) / reference * 100.0))
    return sum(errors) / len(errors), max(errors)


# L8/P4-13 external traffic scope.
#
# CORRECTION (found by re-running the phase-1/2/3 validation from scratch): this gate first
# scoped the traffic axes using a RAW point comparison of NPUsim's counters against Table V,
# reporting 30.5%/60.5% MAPE and a 133.94% max. That comparison is not apples-to-apples and the
# repository already held a better-founded one in validation/phase3/compare.py. Both axes now
# use that methodology, single-sourced in validation/phase3/traffic_reference.py:
#
#   DRAM -- Table V is measured AFTER the chip's run-length compression, so the reference is
#           converted to a dense equivalent before comparing. That is the difference between
#           "60.5% MAPE / 133.9% max" and the real 22.9% / 50.0%.
#   GLB  -- Table V depends on per-layer mapper parameters the paper does not publish, so no
#           point target exists. The verdict is whether the chip falls inside
#           [perfect retention, literal per-repetition streaming].
#
# Neither axis is inside the 15% milestone limit yet, so neither is accuracy-gated. What IS
# gated is that they do not get WORSE: each carries a regression ceiling pinned to today's
# measurement, so --check-traffic always means something while 15% stays the target.
TRAFFIC_AXES = {
    "DRAM": {
        "ceiling": 50.0,
        "milestone": 15.0,
        "reason": "input halo reuse is now capacity-gated and exact: traffic T3 is 1.00x the "
                  "full mapped input union in every validated CONV layer, while weight/output "
                  "remain fetch-once. The residual error has the opposite sign: the silicon "
                  "moves more DRAM data than this ideal dense-union schedule. Its unpublished "
                  "per-layer mapper/refetch behavior and the Table-V RLC dense-equivalent "
                  "conversion are not separable from public data, so this remains informational",
    },
}
# GLB has no point target; its verdict is the band. The gate requires that no layer's chip
# figure sits ABOVE NPUsim's literal-streaming upper bound by more than today's measurement --
# being above the upper bound means the chip streamed more through its GLB than THIS mapping does.
# CORRECTED 2026-08-20: an earlier note here called that a missing traffic SOURCE no mapping could
# reach. Measuring refuted it -- moving output factors up the hierarchy takes conv3 from 9.2 MB to
# 32.2 MB to 101.2 MB (configs/accelerators/eyeriss_psumprobe*.cfg), so the chip's 50.2 MB is
# reachable here. The real limit is that those mappings cost 4x-50x the measured computation
# cycles: psum spilling follows the tiling hierarchy rather than GLB capacity (E20-3).
# 2026-08-25 (E20-3b): 0. Every layer is now inside the band. The GLB psum LOAD used to be
# suppressed by a latch in global_buffer_t::data_transfer that never reset within a layer, so the
# array wrote its psum tile out once per revisit but read it back exactly once per output offset
# (Eyeriss conv3: 19656 write-backs against 312 loads -- a psum written and never read cannot be
# accumulated). Giving the load side the same psum_retention_valid condition the write-back side
# already carried took conv3's GLB upper bound 40.9 -> 73.1 MB with the compute-schedule latency
# unchanged. Keep this at 0: any layer drifting back above the upper bound is a regression.
GLB_BAND_CEILING_LAYERS = 0


def eyeriss(check_baseline: bool) -> tuple[float, float, dict]:
    baseline = {
        row["layer"]: (
            float(row["computation_cycles"]), float(row["fold_fill_cycles"]),
            float(row["glb_mb"]), float(row["dram_mb"]),
        )
        for row in rows(VALIDATION / "phase3" / "npusim_baseline.csv")
    }
    latency_errors: list[float] = []
    for row in rows(VALIDATION / "phase3" / "golden_eyeriss_silicon.csv"):
        path = RESULT / "eyeriss" / "alexnet" / "silicon" / f"layer_{row['layer_index']}.txt"
        text = path.read_text(encoding="utf-8")
        glb_section = text.split("Global buffer result", 1)[1]
        transaction_table = glb_section.split("serialized", 1)[1].split("# of request", 1)[0]
        # P1-B: the Eyeriss reference counts GLB *SRAM* traffic. A GLB-bypassed datatype
        # (Eyeriss streams filter weights from the chip ingress straight into the PE
        # spads) traverses the on-chip fabric -- and is charged for it -- but never
        # touches the SRAM, so it must be excluded here. The simulator now names the
        # bypassed datatypes explicitly instead of leaving them to be inferred from a
        # zero counter.
        bypassed = set(re.search(r"GLB-bypassed \(direct stream\)\s*:(.*)", glb_section).group(1).split())
        label_of = {"Input data": "input", "Weight": "weight", "Output data": "output"}
        glb_transactions = sum(
            int(serialized) for name, serialized in
            re.findall(r"(Input data|Weight|Output data)\s*:\s*\d+/\d+/(\d+)", transaction_table)
            if label_of[name] not in bypassed
        )
        dram_section = text.split("Multi-chip <-> DRAM transactions", 1)[1]
        dram_transactions = sum(int(value) for value in re.findall(r"/(\d+)\n", dram_section)[:3])
        measured = (
            metric(text, "Computation cycle"), metric(text, "Fold fill cycle", default=0.0),
            glb_transactions * 32 / 1e6, dram_transactions * 8 / 1e6,
        )
        name = row["layer"]
        if check_baseline and measured != baseline[name]:
            raise AssertionError(f"Eyeriss NPUsim baseline changed at {name}: {measured} != {baseline[name]}")
        # P0-3: the golden comparison reads the official "Compute-schedule latency"
        # metric directly instead of re-deriving it from the two component fields.
        compute_schedule = compute_schedule_latency(text, str(path))
        latency_ms = compute_schedule / 200_000.0
        latency_errors.append(abs((latency_ms - float(row["processing_ms"])) / float(row["processing_ms"]) * 100.0))
    # P0-2 network contract: verify the aggregated network.txt (not just per-layer
    # files) carries a self-consistent official metric and populated busy axes.
    check_network_timeline(RESULT / "eyeriss" / "alexnet" / "silicon" / "network.txt")
    check_network_rollup(RESULT / "eyeriss" / "alexnet" / "silicon")
    # P4-13: traffic verdicts come from the shared phase-3 methodology.
    sys.path.insert(0, str(VALIDATION / "phase3"))
    import traffic_reference
    traffic = traffic_reference.evaluate(RESULT / "eyeriss" / "alexnet" / "silicon")
    return sum(latency_errors) / len(latency_errors), max(latency_errors), traffic


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-baseline", action="store_true", help="require exact current NPUsim counters")
    parser.add_argument("--check-traffic", action="store_true", help="enable the milestone-3 traffic gate")
    arguments = parser.parse_args()

    gemmini_mape, gemmini_max = gemmini(arguments.check_baseline)
    eyeriss_mape, eyeriss_max, traffic_errors = eyeriss(arguments.check_baseline)
    print(f"Gemmini RTL: MAPE={gemmini_mape:.2f}% max={gemmini_max:.2f}% (limits 8%/8%)")
    print(f"Eyeriss silicon latency: MAPE={eyeriss_mape:.2f}% max={eyeriss_max:.2f}% (limits 5%/8%)")

    failed = gemmini_mape > 8.0 or gemmini_max > 8.0 or eyeriss_mape > 5.0 or eyeriss_max > 8.0

    # P4-13: external traffic scope, per axis, using the phase-3 methodology.
    print("Eyeriss external traffic (milestone-3 limit 15%; methodology: "
          "validation/phase3/traffic_reference.py):")

    # DRAM: measured against the RLC dense-equivalent chip figure.
    dram_errors = [(name, verdict["dram_error_pct"]) for name, verdict in traffic_errors.items()]
    worst_layer, worst = max(dram_errors, key=lambda entry: abs(entry[1]))
    dram_mape = sum(abs(value) for _, value in dram_errors)/len(dram_errors)
    scope = TRAFFIC_AXES["DRAM"]
    within = abs(worst) <= scope["milestone"]
    print(f" * DRAM MAPE={dram_mape:6.2f}%  max={worst:+7.2f}% ({worst_layer})  "
          f"ceiling {scope['ceiling']:.0f}%  "
          f"{'within milestone' if within else 'OUT OF GATE SCOPE'}")
    print("        vs RLC dense-equivalent chip DRAM (Table V is post-compression): " +
          ", ".join(f"{name} {value:+.1f}%" for name, value in dram_errors))
    if not within:
        print(f"        not gated because: {scope['reason']}")
    if abs(worst) > scope["ceiling"]:
        print(f"        REGRESSION: {abs(worst):.2f}% exceeds this axis's "
              f"{scope['ceiling']:.0f}% ceiling")
        if arguments.check_traffic:
            failed = True

    # GLB: a band, not a point. Report where the chip sits relative to it.
    above = [name for name, verdict in traffic_errors.items() if verdict["glb_above_upper"]]
    inside = [name for name, verdict in traffic_errors.items() if verdict["glb_in_band"]]
    print(f" * GLB  band [perfect retention, literal streaming]: {len(inside)} of "
          f"{len(traffic_errors)} layers inside")
    for name, verdict in traffic_errors.items():
        placement = ("inside" if verdict["glb_in_band"] else
                     "ABOVE upper" if verdict["glb_above_upper"] else "below lower")
        print(f"        {name}: lower {verdict['glb_lower']:.1f} <= chip "
              f"{verdict['glb_chip']:.1f} <= upper {verdict['glb_upper']:.1f} MB -> {placement}")
    if above:
        print(f"        {', '.join(above)}: the chip moved MORE through its GLB than the model "
              f"streams through its GLB in THIS mapping. Check the psum round trip first -- the "
              f"load and write-back sides must agree (E20-3b); a write-back without a matching "
              f"reload understates this bound without moving the compute-schedule latency")
    if len(above) > GLB_BAND_CEILING_LAYERS:
        print(f"        REGRESSION: {len(above)} layers above the upper bound exceeds the "
              f"{GLB_BAND_CEILING_LAYERS}-layer ceiling")
        if arguments.check_traffic:
            failed = True
    if arguments.check_traffic:
        print("   --check-traffic gates: regression ceilings only (DRAM "
              f"{TRAFFIC_AXES['DRAM']['ceiling']:.0f}%, GLB {GLB_BAND_CEILING_LAYERS} layers "
              "above band); neither axis is inside the 15% milestone yet")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
