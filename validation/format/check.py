#!/usr/bin/env python3
"""Format-IP axis regression with a NON-ZERO unit cost (L9 / P4-1).

The format-IP axis (payload/metadata packing and accumulator spill on the LB<->MAC path) is
wired into the PE busy combination, the network rollup and the report -- but every shipped
config leaves its unit cost at 0, so none of that wiring was ever exercised by a regression
and a slow format IP could have been silently free. These three configs differ from
configs/accelerators/gemmini.cfg in the format cost and the PE buffering only:

    gemmini_format         format_payload_cycle = 4,  PE double_buffer = 1
    gemmini_format_heavy   format_payload_cycle = 40, PE double_buffer = 1
    gemmini_format_single  format_payload_cycle = 4,  PE double_buffer = 0

Checks (asserted; non-zero exit on failure):
  FM1  The axis is actually non-zero in all three runs (the knob reaches the accounting).
  FM2  It is LINEAR in the unit cost: 10x the cost is exactly 10x the axis. That pins the
       axis to the format events rather than to anything else that moved.
  FM3  Double-buffered PE busy == max(compute, access, link, overlap, format). With the cheap
       cost the format hides behind the transfer makespan; with the heavy cost it DOMINATES
       and the busy value equals the format axis -- so a slow format IP is neither invisible
       nor silently fully hidden.
  FM4  Single-buffered PE busy == compute + max(access, link, overlap) + format, exactly.
       The format IP sits on the same non-overlappable tile-load path as the transfer.
  FM5  The network rollup sums the axis (it is a real busy axis, so network.txt's printed axes
       must account for it) -- checked here against a NON-ZERO value, which is what the
       zero-cost configs could never do.
  FM6  A slower format IP lengthens the critical path. Whatever else changes, the axis must
       not be free.

The METADATA stream and the format-IP ENERGY (FM7-FM11). Three of the four Format-IP unit costs --
`format_payload_energy`, `format_metadata_cycle`, `format_metadata_energy` -- were read by the code
and declared by NO config, so Format-IP energy was 0 in every run and the block-scale metadata
stream was never priced at all. The checks above validate the payload LATENCY and nothing else.
(Found by the undeclared-knob scan, validation/knobs KN9; the RE5 energy schema had to be completed
first, since it had been built from the keys the configs happened to declare and rejected two of
these as typos -- see KN10.)

gemmini_format_mxfp declares all four costs and switches the OPERANDS to mxfp8, leaving the output
int8. mxfp8's payload is 8 bits, exactly like int8, so the fixture differs from gemmini_format by
the metadata stream ALONE -- the payload traffic is untouched, which is what makes the two terms
separable.

  FM7  The metadata axis is exercised, and only where it should be: non-zero for the two mxfp8
       operands, exactly 0 for the int8 output (a format with no block scales has no metadata),
       and 0 everywhere in the all-int8 fixture.
  FM8  HAND IDENTITY by isolation: zeroing one unit cost isolates the other, and payload-only plus
       metadata-only must reconstruct the baseline exactly, on BOTH the cycle and the energy axis,
       for all three streams. The derived counts then follow -- payload and metadata are separate
       terms, not one blended "format" quantity.
  FM9  Each of the four costs is linear in its own term and moves no other.
  FM10 ENERGY/TIME SEPARATION, and RE8 CONTAINMENT: the cycle costs move no energy axis and the
       energy costs move no cycle axis, and the total Format-IP energy enters the MAC component
       subtotal by exactly its own value, once.
  FM11 GRANULARITY CONTRACT, stated rather than assumed. A block holds 32 elements, yet the
       metadata TRANSACTION count equals the payload count for a block-scaled operand: at MAC-tile
       granularity both ceil() to one transaction. That is an ACCESS count, not a volume ratio --
       the format IP must fetch a block's scale to decode any part of that block, so one scale
       access per tile is the intended reading. It is pinned here so the equality is a declared
       contract rather than an accident, and so the LIMITATION is on the record: this axis cannot
       express the 32:1 volume ratio between payload and metadata, exactly as RE1 found for the
       accumulator (which was moved to per-byte charging for that reason).
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKLOAD = ("gemm_64x64x64", "ws")
CHEAP, HEAVY, SINGLE = "gemmini_format", "gemmini_format_heavy", "gemmini_format_single"
VARIANTS = (CHEAP, HEAVY, SINGLE)
COST_RATIO = 10.0     # 40 / 4
# FM7-FM11: gemmini_format_mxfp declares all four Format-IP unit costs, with mxfp8 operands.
MXFP = "gemmini_format_mxfp"
PAYLOAD_CYCLE, METADATA_CYCLE = 4.0, 3.0
PAYLOAD_ENERGY, METADATA_ENERGY = 0.5, 0.25
STREAMS = ("Input", "Weight", "Output")


def run(target: str):
    directory = ROOT / "result" / target / WORKLOAD[0] / WORKLOAD[1]
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, WORKLOAD[0], WORKLOAD[1]],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    layer = directory / "layer_0.txt"
    network = directory / "network.txt"
    if not layer.exists() or not network.exists():
        raise SystemExit(f"missing simulator output under {directory}")
    return layer.read_text(encoding="utf-8"), network.read_text(encoding="utf-8")


def timeline(text: str) -> dict:
    block = text.split("Layer timeline", 1)[1].split("MAC result", 1)[0]
    axes = [[float(a), float(b), float(c)] for a, b, c in
            re.findall(r":\s*([\d.]+) /\s*([\d.]+) /\s*([\d.]+) cycles", block)]
    busy = [float(v) for v in re.findall(r"\* [^:\n]+:\s*([\d.]+) cycles \(", block)]
    return {
        "format": float(re.search(r"PE format-IP axis\s*:\s*([\d.]+)", block).group(1)),
        "schedule": float(re.search(r"Compute-schedule latency\s*:\s*([\d.]+)", block).group(1)),
        "critical": float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", block).group(1)),
        "double_buffered": "double-buffered" in re.search(
            r"PE local buffer\s*:\s*(\S+)", block).group(1),
        "pe_busy": busy[4],
        "pe_axes": axes[4],
    }


def main() -> int:
    measured = {}
    for name in VARIANTS:
        layer_text, network_text = run(name)
        measured[name] = timeline(layer_text)
        measured[name]["network_format"] = timeline(network_text)["format"]
    failures = []
    tol = 0.5

    # FM1
    for name in VARIANTS:
        if measured[name]["format"] <= 0.0:
            failures.append(f"FM1 {name}: format axis is {measured[name]['format']}, so the "
                            f"non-zero unit cost never reached the accounting")

    # FM2
    expected_heavy = measured[CHEAP]["format"]*COST_RATIO
    if abs(measured[HEAVY]["format"] - expected_heavy) > tol:
        failures.append(f"FM2: {COST_RATIO}x the unit cost gave {measured[HEAVY]['format']}, "
                        f"expected {expected_heavy} (the axis must be linear in the cost)")

    # FM3
    for name in (CHEAP, HEAVY):
        state = measured[name]
        if not state["double_buffered"]:
            failures.append(f"FM3 {name}: expected a double-buffered PE (config not picked up?)")
            continue
        expected = max(state["pe_axes"] + [state["schedule"], state["format"]])
        if abs(state["pe_busy"] - expected) > tol:
            failures.append(f"FM3 {name}: PE busy {state['pe_busy']} != "
                            f"max(compute, axes, format) {expected}")
    if abs(measured[HEAVY]["pe_busy"] - measured[HEAVY]["format"]) > tol:
        failures.append(f"FM3: with the heavy format cost the axis "
                        f"({measured[HEAVY]['format']}) should dominate PE busy "
                        f"({measured[HEAVY]['pe_busy']}), i.e. it must not stay hidden")
    if abs(measured[CHEAP]["pe_busy"] - measured[CHEAP]["format"]) <= tol:
        failures.append("FM3: with the cheap format cost the axis should still be hidden "
                        "behind the transfer makespan, so the two cases are distinguishable")

    # FM4
    single = measured[SINGLE]
    if single["double_buffered"]:
        failures.append("FM4: expected a single-buffered PE for the single-buffer variant")
    else:
        expected = single["schedule"] + max(single["pe_axes"]) + single["format"]
        if abs(single["pe_busy"] - expected) > tol:
            failures.append(f"FM4: single-buffered PE busy {single['pe_busy']} != compute + "
                            f"max(axes) + format {expected}")

    # FM5
    for name in VARIANTS:
        if abs(measured[name]["network_format"] - measured[name]["format"]) > tol:
            failures.append(f"FM5 {name}: network format axis "
                            f"{measured[name]['network_format']} != the sum of its layers' "
                            f"{measured[name]['format']}")
        if measured[name]["network_format"] <= 0.0:
            failures.append(f"FM5 {name}: network format axis is zero despite a non-zero cost")

    # FM6
    if measured[HEAVY]["critical"] <= measured[CHEAP]["critical"] + tol:
        failures.append(f"FM6: a {COST_RATIO}x slower format IP did not lengthen the critical "
                        f"path ({measured[CHEAP]['critical']} -> "
                        f"{measured[HEAVY]['critical']})")

    # ---- FM7-FM11: the metadata stream and the format-IP energy -----------------------------
    import shutil                                                          # noqa: E402
    sys.path.insert(0, str(ROOT / "validation"))
    import perturb_lib as pl                                               # noqa: E402

    def format_axes(text):
        mac = text.split("MAC result", 1)[1].split("PE result", 1)[0]
        cycles = {name: float(value) for name, value in re.findall(
            r"\*\s*(Input|Weight|Output)\s*:\s*([\d.]+) cycles",
            mac.split("Format-IP cycle", 1)[1])[:3]}
        energy = {name: float(value) for name, value in re.findall(
            r"\*\s*(Input|Weight|Output)\s*:\s*([\d.]+) \S+",
            mac.split("Format-IP energy", 1)[1])[:3]}
        mac_row = float(re.search(r"\* MAC \(compute\+format\)\s+([\d.]+)", text).group(1))
        return {"cycles": cycles, "energy": energy, "mac_row": mac_row}

    def mxfp_variant(label, edits):
        name = f"__fm_{label}"
        pl.variant(MXFP, name, edits)
        shutil.rmtree(pl.MAPPINGS / name, ignore_errors=True)
        shutil.copytree(pl.MAPPINGS / MXFP, pl.MAPPINGS / name)
        text = pl.run(name, *WORKLOAD)
        pl.discard(name)
        if text is None:
            raise SystemExit(f"the {label} format variant does not run")
        return format_axes(text)

    def zero(key):
        return lambda text: re.sub(rf"^{re.escape(key)} = .*$", f"{key} = 0:0:0", text,
                                   count=1, flags=re.M)

    def double(key, value):
        return lambda text: re.sub(rf"^{re.escape(key)} = .*$",
                                   f"{key} = {value*2}:{value*2}:{value*2}", text,
                                   count=1, flags=re.M)

    def compose(*edits):
        def apply(text):
            for edit in edits:
                text = edit(text)
            return text
        return apply

    mxfp = mxfp_variant("base", lambda text: text)
    payload_cyc = mxfp_variant("pcyc", zero("format_metadata_cycle"))
    metadata_cyc = mxfp_variant("mcyc", zero("format_payload_cycle"))
    payload_en = mxfp_variant("pen", zero("format_metadata_energy"))
    metadata_en = mxfp_variant("men", zero("format_payload_energy"))
    all_free = mxfp_variant("free", compose(zero("format_payload_cycle"),
                                            zero("format_metadata_cycle"),
                                            zero("format_payload_energy"),
                                            zero("format_metadata_energy")))
    int8_axes = format_axes((ROOT / "result" / CHEAP / WORKLOAD[0] / WORKLOAD[1] /
                             "layer_0.txt").read_text(encoding="utf-8"))

    # FM7
    for stream in ("Input", "Weight"):
        if metadata_cyc["cycles"][stream] <= 0.0:
            failures.append(f"FM7 {stream}: the mxfp8 operand carries no metadata cycle "
                            f"({metadata_cyc['cycles'][stream]}); the block-scale stream is not "
                            f"being priced at all")
        if metadata_en["energy"][stream] <= 0.0:
            failures.append(f"FM7 {stream}: the mxfp8 operand carries no metadata energy")
    if metadata_cyc["cycles"]["Output"] != 0.0 or metadata_en["energy"]["Output"] != 0.0:
        failures.append(f"FM7 Output: int8 has no block scales, so its metadata term must be 0, "
                        f"got cycle {metadata_cyc['cycles']['Output']} / energy "
                        f"{metadata_en['energy']['Output']}")
    for stream in STREAMS:
        if int8_axes["energy"][stream] != 0.0:
            failures.append(f"FM7 {stream}: the all-int8 fixture reports "
                            f"{int8_axes['energy'][stream]} of format energy but declares no "
                            f"format energy cost")

    # FM8 -- isolation reconstructs the baseline, and the derived counts follow
    counts = {}
    for stream in STREAMS:
        for axis, isolated_a, isolated_b, cost_a, cost_b in (
                ("cycles", payload_cyc, metadata_cyc, PAYLOAD_CYCLE, METADATA_CYCLE),
                ("energy", payload_en, metadata_en, PAYLOAD_ENERGY, METADATA_ENERGY)):
            total = isolated_a[axis][stream] + isolated_b[axis][stream]
            if abs(total - mxfp[axis][stream]) > 0.5:
                failures.append(f"FM8 {stream} {axis}: the isolated payload and metadata terms "
                                f"({isolated_a[axis][stream]} + {isolated_b[axis][stream]}) do not "
                                f"reconstruct the reported {mxfp[axis][stream]}")
            counts[(stream, axis, "payload")] = isolated_a[axis][stream]/cost_a
            counts[(stream, axis, "metadata")] = isolated_b[axis][stream]/cost_b
    for stream in STREAMS:
        if all_free["cycles"][stream] != 0.0 or all_free["energy"][stream] != 0.0:
            failures.append(f"FM8 {stream}: with all four costs at 0 the axis still reports "
                            f"cycle {all_free['cycles'][stream]} / energy "
                            f"{all_free['energy'][stream]}")

    # FM9 -- each cost linear in its own term only
    for key, value, axis, term, cost in (
            ("format_payload_cycle", PAYLOAD_CYCLE, "cycles", "payload", PAYLOAD_CYCLE),
            ("format_metadata_cycle", METADATA_CYCLE, "cycles", "metadata", METADATA_CYCLE),
            ("format_payload_energy", PAYLOAD_ENERGY, "energy", "payload", PAYLOAD_ENERGY),
            ("format_metadata_energy", METADATA_ENERGY, "energy", "metadata", METADATA_ENERGY)):
        doubled = mxfp_variant(f"x2_{key}", double(key, value))
        for stream in STREAMS:
            expected = mxfp[axis][stream] + counts[(stream, axis, term)]*cost
            if abs(doubled[axis][stream] - expected) > 0.5:
                failures.append(f"FM9 {key} {stream}: doubling it gave "
                                f"{doubled[axis][stream]}, expected {expected} "
                                f"({counts[(stream, axis, term)]:.0f} {term} events x {cost})")
            other = "energy" if axis == "cycles" else "cycles"
            if doubled[other][stream] != mxfp[other][stream]:
                failures.append(f"FM10 {key} {stream}: it also moved the {other} axis "
                                f"({mxfp[other][stream]} -> {doubled[other][stream]})")

    # FM10 -- RE8 containment of the format energy in the MAC subtotal
    energy_total = sum(mxfp["energy"][stream] for stream in STREAMS)
    row_delta = mxfp["mac_row"] - all_free["mac_row"]
    if abs(row_delta - energy_total) > 0.5:
        failures.append(f"FM10: the {energy_total} of Format-IP energy moved the MAC subtotal by "
                        f"{row_delta}; it must enter its own component exactly once")

    # FM11 -- the granularity contract, stated
    for stream in ("Input", "Weight"):
        payload_events = counts[(stream, "energy", "payload")]
        metadata_events = counts[(stream, "energy", "metadata")]
        if abs(payload_events - metadata_events) > 0.5:
            failures.append(f"FM11 {stream}: metadata events {metadata_events:.0f} != payload "
                            f"events {payload_events:.0f}. The declared contract is ONE block-scale "
                            f"access per tile -- a tile fits inside one 32-element block, so both "
                            f"round up to one transaction. If this ratio has become a volume ratio "
                            f"the axis changed meaning and the comment above needs rewriting")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"FM1 format axis non-zero in all three runs ok")
    print(f"FM2 axis linear in the unit cost: {measured[CHEAP]['format']} -> "
          f"{measured[HEAVY]['format']} ok")
    print(f"FM3 double-buffered PE busy == max(compute, axes, format); hidden at cost 4 "
          f"({measured[CHEAP]['pe_busy']}), dominant at cost 40 "
          f"({measured[HEAVY]['pe_busy']}) ok")
    print(f"FM4 single-buffered PE busy == compute + max(axes) + format = "
          f"{single['pe_busy']} ok")
    print(f"FM5 network format axis carries the layer sum ok")
    print(f"FM6 critical path {measured[CHEAP]['critical']} -> {measured[HEAVY]['critical']} ok")
    print(f"FM7 mxfp8 operands carry a metadata term; the int8 output carries none ok")
    print(f"FM8 payload-only + metadata-only == baseline on both axes; all costs 0 -> axis 0 ok")
    print(f"FM9 each of the four unit costs is linear in its own term ok")
    print(f"FM10 cycle costs move no energy and vice versa; the {energy_total:.0f} of format energy "
          f"enters the MAC subtotal once ok")
    print(f"FM11 granularity contract: "
          f"{counts[('Input', 'energy', 'metadata')]:.0f} metadata accesses == "
          f"{counts[('Input', 'energy', 'payload')]:.0f} payload transactions (one block scale per "
          f"tile; NOT the 32:1 volume ratio) ok")
    print("ALL FORMAT CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
