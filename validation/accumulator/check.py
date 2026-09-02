#!/usr/bin/env python3
"""Accumulator event and precision regression (E4 + RE1).

`accumulator_format` was parsed and printed and then consumed by nothing. An accumulator SPILL
was sized by the OUTPUT datatype, so `accumulator_format = fp32` and `= fp16` produced identical
traffic and identical energy -- the setting changed one string in the report and nothing else.

Two things were wrong, and both are fixed:
  * the spill moves ACCUMULATOR-precision partial sums, not output-precision values, and the
    final cast/pack is a SEPARATE event on the output precision;
  * charging either of them per LINK TRANSACTION cannot express precision at MAC-tile
    granularity -- ceil(elements * width / link_bits) rounds to the same 1 for fp32 and fp16
    alike. They are precision conversions on values, so they are charged per BYTE MOVED.

Two fixtures differ from configs/accelerators/gemmini.cfg only in `accumulator_format` (fp32 vs
fp16), with int8 operands and int8 output throughout, and distinct unit costs for the two events
(accumulator_spill_energy = 0.5, output_cast_energy = 0.25 on the output tensor).

RE1 -- the events themselves were wrong, not just their precision. `account_format_events(OUTPUT)`
fired on every output request and charged a spill AND a cast together, before it could even tell
whether the request was a fresh zero-initialized accumulator or a genuine reload. Consequences:
  * the "final cast" count followed the MAC count (262,144 on GEMM 64x64x64) instead of the final
    output element count (4,096) -- a 64x overcharge;
  * a zero-init paid reload energy it never spent;
  * with edge_accumulation the accumulator energy was billed to the PE even though the config says
    the accumulator sits at the array edge.
The four events now come from their own boundaries -- create, reload, spill, final cast -- and the
final cast is charged at the off-chip output store, the one boundary DR6/T1 establish fires exactly
once per output element. A ratio-only check could not see any of this, so the checks below pin
ABSOLUTE counts against the workload.

Checks (asserted; non-zero exit on failure):
  AC6  ABSOLUTE COUNT (RE1): the INT8 final cast equals the layer's final output element count --
       4,096 for GEMM 64x64x64 (64x64 outputs), not the 262,144 MAC issues. This is the check that
       distinguishes "cast once per output" from "cast once per multiply".
  AC7  The cast count is INVARIANT to the reduction depth: the accumulator format changes the spill
       and reload traffic but must not change how many outputs are produced.
  AC8  A fresh accumulator is a create, not a reload: the create count is non-zero and reload and
       spill are tracked separately, so a zero-init cannot be billed as a read-back.
  AC9  OWNERSHIP: with edge_accumulation the accumulator energy appears in the PE-array subtotal and
       NOT in the PE's own axis, matching the config's statement about where the accumulator lives.
  AC1  Both accumulator-side events follow the ACCUMULATOR width: reload and spill traffic each
       scale 2x from fp16 to fp32. (The old form compared spill against cast, which only worked
       because the two were produced by the same call -- exactly the bug RE1 describes.)
  AC2  Changing ONLY the accumulator format changes the spill: fp32 moves exactly 2x fp16.
  AC3  ... and leaves the output cast UNTOUCHED: the cast bytes are identical in both runs, since
       output_format did not change. Together with AC2 this shows the two events are separate
       rather than one merged "format" quantity.
  AC4  HAND IDENTITY, energy, ON THE OWNING COMPONENT:
         accumulator energy == (reload bytes + spill bytes) * accumulator_spill_energy
         cast energy        == cast bytes * output_cast_energy
       priced independently, at their own precisions, and reported by the component that performs
       each -- the array edge for the accumulator, the off-chip store for the cast.
  AC5  The lane->accumulator REDUCTION energy is reported on its own axis instead of being folded
       into the MAC total, so adder-tree work is attributable. (Its value is 0 here because
       mac_reduction_energy is not calibrated in this config; the axis must still exist.)

RE1's LATENCY HALF. Everything above prices the four events in ENERGY. Their TIME was modelled all
along -- `accumulator_spill_cycle` charges into Format-IP cycle, `output_cast_cycle` into the
multi-chip write-back -- but NO config declared either key, so a spill and a cast were free in
wall-clock in every fixture and the two knobs could have been wired to anything. That is invisible
to a value check and even to the dead-knob sweep, which perturbs only what a config declares; the
undeclared-knob scan (validation/knobs KN9) is what surfaced it. The accumulator fixtures now
declare both (spill 0.25, cast 0.5 cycles per byte), and the cast's cycle count is reported
(`Output cast cycle`) because it previously reached the report only as a MAX inside the fabric's
busy time -- a cast cheaper than the access/transfer axes was unobservable.

  AC10 SPILL/RELOAD TIME follows the ACCUMULATOR width, exactly as AC1's traffic does: the
       Format-IP output cycle is exactly 2x for fp32 vs fp16. This is the time mirror of AC1/AC2.
  AC11 CAST TIME, hand identity: Output cast cycle == cast bytes x output_cast_cycle
       (4,096 x 0.5 = 2,048), and INVARIANT to the accumulator format -- the time mirror of
       AC3/AC7.
  AC12 Each cycle knob scales its OWN axis linearly and only its own: doubling
       accumulator_spill_cycle doubles the Format-IP output cycle and leaves the cast cycle alone,
       and vice versa. With a knob undeclared its axis is exactly 0 -- the state that hid both.
RETAINED OUTPUT (AC14-AC16, E20-4). A reload is the physical READ-BACK of a partial sum over the
LB->MAC path. The accounting used to fire before the guard that skips that path when the output tile
is already resident in the MAC:

    account_accumulator_reload(tile);          // charged unconditionally
    if(!skip_transfer[OUTPUT]) { ... }         // the read-back that may not happen

Every fixture here missed it, because none of their mappings retains the output: injecting the
pre-fix order changes nothing in gemmini_accum_fp32 (1,048,560 bytes either way). gemmini_accum_retained
puts ONLY the reduction dimension C on the PE temporal level, so the MAC's output tile equals the
local buffer's and the read-back is skipped on 16,127 of the passes. There the pre-fix order charges
1,048,512 bytes of reload against 16,384 actually performed -- a 64x overcharge in both accumulator
traffic and accumulator energy, for read-backs that never happen.

  AC14 ABSOLUTE COUNT: on the retained mapping the reload is 16,384 bytes -- the read-backs that
       physically occur -- not the 1,048,512 the pre-fix accounting charged. Retention is visible as
       its own counter rather than as an absent event.
  AC15 Retention is a MAPPING property, not a precision one: the retained and create counts are
       IDENTICAL across fp32 and fp16 accumulators, while reload and spill bytes both scale exactly
       2x. A model that billed retained passes as reloads would make the reload count scale with
       retention too, which this separates.
  AC16 All four states are simultaneously distinguishable in one report -- create, retained, reload
       and spill are each non-zero and reload != spill, since a retained accumulator still spills
       when it finally leaves.
  AC13 ENERGY/TIME SEPARATION: doubling an energy unit cost moves no cycle axis, and doubling a
       cycle unit cost moves no energy axis. RE1 split these four events onto their own boundaries;
       this asserts the two dimensions of each event stay independent.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKLOAD = ("gemm_64x64x64", "ws")
# accumulator format -> (bits, target)
ACCUMULATORS = {"fp32": 32, "fp16": 16}
OUTPUT_BITS = 8          # output_format = int8
SPILL_ENERGY = 0.5       # accumulator_spill_energy on the output tensor
CAST_ENERGY = 0.25       # output_cast_energy on the output tensor
SPILL_CYCLE = 0.25       # accumulator_spill_cycle on the output tensor (RE1 latency half)
# E20-4: the retained-output mapping and its fp16 twin.
RETAINED = {"fp32": "gemmini_accum_retained", "fp16": "gemmini_accum_retained_fp16"}
RETAINED_RELOAD_BYTES = 16384      # the read-backs that physically occur, fp32
RETAINED_PREFIX_RELOAD = 1048512   # what the pre-E20-4 accounting charged instead
RETAINED_EVENTS = 16127            # passes where the accumulator stayed in the MAC
CAST_CYCLE = 0.5         # output_cast_cycle on the output tensor


def run(target: str) -> str:
    layer = ROOT / "result" / target / WORKLOAD[0] / WORKLOAD[1] / "layer_0.txt"
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, WORKLOAD[0], WORKLOAD[1]],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if not layer.exists():
        raise SystemExit(f"missing simulator output: {layer}")
    return layer.read_text(encoding="utf-8")


def measure(text: str) -> dict:
    mac = text.split("MAC result", 1)[1].split("PE result", 1)[0]
    summary = text.split("Energy summary", 1)[1].split("Power summary", 1)[0]
    # The MAC section labels these "Input"/"Weight"/"Output" (no "data" suffix).
    format_energy = [float(v) for v in re.findall(
        r"\*\s*(?:Input|Weight|Output)\s*:\s*([\d.]+) (?:pJ|normalized|uncalibrated)",
        mac.split("Format-IP energy", 1)[1])[:3]]
    reduction = re.search(r"Reduction energy\s*:\s*([\d.]+) (?:pJ|normalized|uncalibrated)", mac)
    def summary_row(component):
        row = re.search(rf"\* {re.escape(component)}\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", summary)
        return float(row.group(1)) if row else None
    return {
        "spill_bytes": int(re.search(r"Accumulator spill bytes\s*:\s*(\d+)", mac).group(1)),
        "reload_bytes": int(re.search(r"Accumulator reload bytes\s*:\s*(\d+)", mac).group(1)),
        "create_events": int(re.search(r"Accumulator create events\s*:\s*(\d+)", mac).group(1)),
        # RE1: the cast is reported by the component that commits it (the off-chip store), not
        # by the MAC section where the accumulator events live.
        "cast_bytes": int(re.search(r"Output cast bytes\s*:\s*(\d+)",
                                    text.split("Multi chip result", 1)[1]).group(1)),
        "macs": int(re.search(r"# of computations\s*:\s*(\d+)", mac).group(1)),
        # The PE-side accumulator energy is labelled "(PE)" to distinguish it from the edge one.
        "accumulator_energy": float(re.search(
            r"Accumulator energy \(PE\)\s*:\s*([\d.]+) (?:pJ|normalized|uncalibrated)", mac).group(1)),
        "mac_subtotal": summary_row("MAC (compute+format)"),
        "pe_array_subtotal": summary_row("PE array"),
        "edge_accumulator_energy": float(re.search(
            r"Accumulator energy \(edge\)\s*:\s*([\d.]+) (?:pJ|normalized|uncalibrated)",
            text.split("PE array result", 1)[1]).group(1)),
        "cast_energy": float(re.search(
            r"Output cast energy\s*:\s*([\d.]+) (?:pJ|normalized|uncalibrated)",
            text.split("Multi chip result", 1)[1]).group(1)),
        "format_energy_output": format_energy[2],
        # RE1 latency half: the two cycle axes the energy checks above have no view of.
        "format_cycle_output": float(re.search(
            r"Format-IP cycle.*?\* Output\s*:\s*([\d.]+) cycles", mac, re.S).group(1)),
        "cast_cycle": float(re.search(
            r"Output cast cycle\s*:\s*([\d.]+) cycles",
            text.split("Multi chip result", 1)[1]).group(1)),
        "reduction_axis_present": reduction is not None,
        "reduction_energy": float(reduction.group(1)) if reduction else None,
    }


def main() -> int:
    measured = {fmt: measure(run(f"gemmini_accum_{fmt}")) for fmt in ACCUMULATORS}
    failures = []
    tol = 0.5

    for fmt, bits in ACCUMULATORS.items():
        state = measured[fmt]
        # AC1: both accumulator-side events follow the accumulator width.
        for label, key in (("reload", "reload_bytes"), ("spill", "spill_bytes")):
            ratio = state[key]/measured["fp16"][key]
            expected_ratio = bits/ACCUMULATORS["fp16"]
            if abs(ratio - expected_ratio) > 1e-6:
                failures.append(f"AC1 {fmt} {label}: {state[key]} is {ratio}x the fp16 value, "
                                f"expected {expected_ratio}x ({bits}-bit vs 16-bit accumulator)")
        # AC4: energy identities on the components that own the events.
        expected_accumulator = (state["reload_bytes"] + state["spill_bytes"])*SPILL_ENERGY
        if abs(state["edge_accumulator_energy"] - expected_accumulator) > tol:
            failures.append(f"AC4 {fmt}: edge accumulator energy "
                            f"{state['edge_accumulator_energy']} != hand {expected_accumulator} "
                            f"(({state['reload_bytes']} reload + {state['spill_bytes']} spill) "
                            f"bytes x {SPILL_ENERGY})")
        expected_cast = state["cast_bytes"]*CAST_ENERGY
        if abs(state["cast_energy"] - expected_cast) > tol:
            failures.append(f"AC4 {fmt}: cast energy {state['cast_energy']} != hand "
                            f"{expected_cast} ({state['cast_bytes']} cast bytes x {CAST_ENERGY})")
        # AC5
        if not state["reduction_axis_present"]:
            failures.append(f"AC5 {fmt}: the reduction energy axis is not reported, so adder-tree "
                            f"work is folded back into the MAC total")

    # AC2
    expected_spill_ratio = ACCUMULATORS["fp32"]/ACCUMULATORS["fp16"]
    measured_spill_ratio = measured["fp32"]["spill_bytes"]/measured["fp16"]["spill_bytes"]
    if abs(measured_spill_ratio - expected_spill_ratio) > 1e-9:
        failures.append(f"AC2: fp32/fp16 spill byte ratio {measured_spill_ratio} != "
                        f"{expected_spill_ratio}; changing only the accumulator format must change "
                        f"the spill")

    # AC3
    if measured["fp32"]["cast_bytes"] != measured["fp16"]["cast_bytes"]:
        failures.append(f"AC3: the output cast changed with the accumulator format "
                        f"({measured['fp32']['cast_bytes']} vs "
                        f"{measured['fp16']['cast_bytes']}) even though output_format is fixed; "
                        f"the spill and the cast are separate events on separate precisions")

    # AC6/AC7: absolute final-cast count. GEMM 64x64x64 produces 64x64 = 4096 INT8 outputs.
    FINAL_OUTPUT_ELEMENTS = 64*64
    for fmt in ACCUMULATORS:
        state = measured[fmt]
        if state["cast_bytes"] != FINAL_OUTPUT_ELEMENTS:
            failures.append(
                f"AC6 {fmt}: final cast is {state['cast_bytes']} bytes, expected "
                f"{FINAL_OUTPUT_ELEMENTS} (one INT8 byte per final output element). "
                f"{state['macs']} would mean one cast per MAC issue instead of per output")
    if measured["fp32"]["cast_bytes"] != measured["fp16"]["cast_bytes"]:
        failures.append(f"AC7: the final cast count changed with the accumulator format "
                        f"({measured['fp32']['cast_bytes']} vs "
                        f"{measured['fp16']['cast_bytes']}); the accumulator precision changes the "
                        f"spill/reload traffic, not how many outputs are produced")

    # AC8: create vs reload are distinct events.
    for fmt in ACCUMULATORS:
        state = measured[fmt]
        if state["create_events"] <= 0:
            failures.append(f"AC8 {fmt}: no accumulator create events reported, so a zero-init "
                            f"cannot be distinguished from a reload")
        if state["reload_bytes"] <= 0:
            failures.append(f"AC8 {fmt}: no accumulator reload traffic reported")
        if state["reload_bytes"] == state["spill_bytes"] and state["create_events"] > 0:
            failures.append(f"AC8 {fmt}: reload ({state['reload_bytes']}) equals spill "
                            f"({state['spill_bytes']}) although {state['create_events']} "
                            f"accumulators were created from zero -- the create path is being "
                            f"charged as a reload")

    # AC9: with edge_accumulation the accumulator energy is the ARRAY's, not the PE's.
    for fmt in ACCUMULATORS:
        state = measured[fmt]
        if state["accumulator_energy"] != 0.0:
            failures.append(f"AC9 {fmt}: this config sets edge_accumulation, so the accumulator "
                            f"energy must not be charged to the PE, but the PE axis reports "
                            f"{state['accumulator_energy']}")
        if state["pe_array_subtotal"] is None or state["pe_array_subtotal"] <= 0.0:
            failures.append(f"AC9 {fmt}: the PE-array subtotal is "
                            f"{state['pe_array_subtotal']}, so the edge accumulator's energy did "
                            f"not reach the component that owns it")

    # ---- RE1 latency half (AC10-AC13) -------------------------------------------------------
    import shutil                                                          # noqa: E402
    sys.path.insert(0, str(ROOT / "validation"))
    import perturb_lib as pl                                               # noqa: E402

    def variant(label, base, edits):
        name = f"__ac_{label}"
        pl.variant(base, name, edits)
        shutil.rmtree(pl.MAPPINGS / name, ignore_errors=True)
        shutil.copytree(pl.MAPPINGS / base, pl.MAPPINGS / name)
        text = pl.run(name, *WORKLOAD)
        pl.discard(name)
        return measure(text) if text is not None else None

    def set_key(key, value):
        return lambda text: re.sub(rf"^{re.escape(key)} = .*$", f"{key} = {value}", text,
                                   count=1, flags=re.M)

    def drop_key(key):
        return lambda text: re.sub(rf"^{re.escape(key)} = .*$", f"#{key} removed", text,
                                   count=1, flags=re.M)

    fp32, fp16 = measured["fp32"], measured["fp16"]
    # AC10 -- spill/reload TIME follows the accumulator width, like AC1's traffic
    if fp16["format_cycle_output"] <= 0.0:
        failures.append(f"AC10: the Format-IP output cycle is {fp16['format_cycle_output']}, so "
                        f"accumulator_spill_cycle is not being charged at all -- the state that "
                        f"made a spill free in wall-clock everywhere")
    elif abs(fp32["format_cycle_output"] - 2*fp16["format_cycle_output"]) > 0.5:
        failures.append(f"AC10: spill/reload TIME {fp32['format_cycle_output']} (fp32) != 2x "
                        f"{fp16['format_cycle_output']} (fp16); the time must follow the "
                        f"accumulator width exactly as the traffic does (AC1)")
    # AC11 -- cast TIME hand identity, invariant to the accumulator format
    for fmt, state in measured.items():
        expected = state["cast_bytes"]*CAST_CYCLE
        if abs(state["cast_cycle"] - expected) > 1e-6:
            failures.append(f"AC11 {fmt}: cast cycle {state['cast_cycle']} != hand {expected} "
                            f"({state['cast_bytes']} bytes x {CAST_CYCLE})")
    if fp32["cast_cycle"] != fp16["cast_cycle"]:
        failures.append(f"AC11: the cast TIME changed with the accumulator format "
                        f"({fp32['cast_cycle']} vs {fp16['cast_cycle']}); the cast moves "
                        f"output-precision values, so it cannot depend on the accumulator width")
    # AC12 -- each cycle knob scales its own axis only; undeclared means exactly 0
    spill_2x = variant("spill2", "gemmini_accum_fp32", set_key("accumulator_spill_cycle", "0:0:0.5"))
    cast_2x = variant("cast2", "gemmini_accum_fp32", set_key("output_cast_cycle", "0:0:1.0"))
    spill_off = variant("spilloff", "gemmini_accum_fp32", drop_key("accumulator_spill_cycle"))
    cast_off = variant("castoff", "gemmini_accum_fp32", drop_key("output_cast_cycle"))
    energy_2x = variant("energy2", "gemmini_accum_fp32",
                        set_key("accumulator_spill_energy", "0:0:1.0"))
    for label, state in (("spill x2", spill_2x), ("cast x2", cast_2x), ("spill off", spill_off),
                         ("cast off", cast_off), ("energy x2", energy_2x)):
        if state is None:
            raise SystemExit(f"the {label} accumulator variant does not run")
    if abs(spill_2x["format_cycle_output"] - 2*fp32["format_cycle_output"]) > 0.5:
        failures.append(f"AC12: doubling accumulator_spill_cycle gave "
                        f"{spill_2x['format_cycle_output']}, expected "
                        f"{2*fp32['format_cycle_output']}")
    if spill_2x["cast_cycle"] != fp32["cast_cycle"]:
        failures.append(f"AC12: doubling accumulator_spill_cycle also moved the cast cycle "
                        f"({fp32['cast_cycle']} -> {spill_2x['cast_cycle']})")
    if abs(cast_2x["cast_cycle"] - 2*fp32["cast_cycle"]) > 1e-6:
        failures.append(f"AC12: doubling output_cast_cycle gave {cast_2x['cast_cycle']}, expected "
                        f"{2*fp32['cast_cycle']}")
    if cast_2x["format_cycle_output"] != fp32["format_cycle_output"]:
        failures.append(f"AC12: doubling output_cast_cycle also moved the spill cycle "
                        f"({fp32['format_cycle_output']} -> {cast_2x['format_cycle_output']})")
    if spill_off["format_cycle_output"] != 0.0:
        failures.append(f"AC12: with accumulator_spill_cycle undeclared the spill still costs "
                        f"{spill_off['format_cycle_output']} cycles; the default must be free")
    if cast_off["cast_cycle"] != 0.0:
        failures.append(f"AC12: with output_cast_cycle undeclared the cast still costs "
                        f"{cast_off['cast_cycle']} cycles; the default must be free")
    # AC13 -- energy and time are independent dimensions of the same event
    for label, state in (("accumulator_spill_cycle", spill_2x), ("output_cast_cycle", cast_2x)):
        for axis in ("edge_accumulator_energy", "cast_energy", "mac_subtotal",
                     "pe_array_subtotal"):
            if state[axis] != fp32[axis]:
                failures.append(f"AC13: doubling {label} moved the ENERGY axis {axis} "
                                f"({fp32[axis]} -> {state[axis]}); a cycle cost must not price "
                                f"anything")
    for axis in ("format_cycle_output", "cast_cycle"):
        if energy_2x[axis] != fp32[axis]:
            failures.append(f"AC13: doubling accumulator_spill_energy moved the TIME axis {axis} "
                            f"({fp32[axis]} -> {energy_2x[axis]}); an energy cost must not add "
                            f"latency")
    if abs(energy_2x["edge_accumulator_energy"] - 2*fp32["edge_accumulator_energy"]) > 0.5:
        failures.append(f"AC13: doubling accumulator_spill_energy gave "
                        f"{energy_2x['edge_accumulator_energy']}, expected "
                        f"{2*fp32['edge_accumulator_energy']}; the energy axis itself must respond")

    # ---- AC14-AC16: retained output -----------------------------------------------------------
    retained = {}
    for fmt, target in RETAINED.items():
        text = run(target)
        block = text.split("MAC result", 1)[1].split("PE result", 1)[0]
        def count(label):
            return int(re.search(rf"{re.escape(label)}\s*:\s*(\d+)", block).group(1))
        retained[fmt] = {"create": count("Accumulator create events"),
                         "retained": count("Accumulator retained events"),
                         "reload": count("Accumulator reload bytes"),
                         "spill": count("Accumulator spill bytes")}
    fp32_r, fp16_r = retained["fp32"], retained["fp16"]
    # AC14
    if fp32_r["retained"] <= 0:
        failures.append(f"AC14: the retained mapping records {fp32_r['retained']} retained passes, "
                        f"so it no longer exercises retention and the check below is vacuous")
    if fp32_r["reload"] != RETAINED_RELOAD_BYTES:
        failures.append(f"AC14: reload is {fp32_r['reload']} bytes on the retained mapping, expected "
                        f"{RETAINED_RELOAD_BYTES} -- the read-backs that physically occur. "
                        f"{RETAINED_PREFIX_RELOAD} would mean the reload is still charged before the "
                        f"skip_transfer guard, billing "
                        f"{RETAINED_PREFIX_RELOAD - RETAINED_RELOAD_BYTES} bytes of read-back that "
                        f"never happens")
    if fp32_r["retained"] != RETAINED_EVENTS:
        failures.append(f"AC14: {fp32_r['retained']} retained passes, expected {RETAINED_EVENTS}; "
                        f"the mapping's retention behaviour changed")
    # AC15
    for field in ("retained", "create"):
        if fp32_r[field] != fp16_r[field]:
            failures.append(f"AC15 {field}: {fp32_r[field]} (fp32) != {fp16_r[field]} (fp16). "
                            f"Retention is decided by the mapping, so it cannot depend on the "
                            f"accumulator width")
    for field in ("reload", "spill"):
        if fp32_r[field] != 2*fp16_r[field]:
            failures.append(f"AC15 {field}: {fp32_r[field]} (fp32) != 2x {fp16_r[field]} (fp16); "
                            f"the read-backs and evictions that DO happen move accumulator-precision "
                            f"values")
    # AC16
    for field in ("create", "retained", "reload", "spill"):
        if fp32_r[field] <= 0:
            failures.append(f"AC16 {field} is {fp32_r[field]} on the retained mapping; all four "
                            f"states must be simultaneously observable for them to be distinguishable")
    if fp32_r["reload"] == fp32_r["spill"]:
        failures.append(f"AC16: reload and spill are both {fp32_r['reload']}; a retained accumulator "
                        f"spills when it finally leaves but is not read back each pass, so the two "
                        f"counts must differ")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"{'accumulator':>12} {'reload B':>10} {'spill B':>10} {'cast B':>8} "
          f"{'edge acc E':>12} {'cast E':>9}")
    for fmt in ("fp32", "fp16"):
        state = measured[fmt]
        print(f"{fmt:>12} {state['reload_bytes']:>10} {state['spill_bytes']:>10} "
              f"{state['cast_bytes']:>8} {state['edge_accumulator_energy']:>12.2f} "
              f"{state['cast_energy']:>9.2f}")
    print(f"AC1 reload and spill both scale 2x with the accumulator width ok")
    print(f"AC2 fp32 spill == 2x fp16 spill ok")
    print(f"AC3 output cast unchanged by the accumulator format ok")
    print(f"AC4 edge accumulator energy == (reload+spill) x {SPILL_ENERGY}; cast energy == "
          f"cast bytes x {CAST_ENERGY}, each on its owning component ok")
    print(f"AC5 reduction energy reported on its own axis ok")
    print(f"AC6 final cast == {FINAL_OUTPUT_ELEMENTS} bytes == the final output element count "
          f"(not the {measured['fp32']['macs']} MAC issues) ok")
    print(f"AC7 the cast count is invariant to the accumulator precision ok")
    print(f"AC8 create ({measured['fp32']['create_events']}) / reload "
          f"({measured['fp32']['reload_bytes']}) / spill ({measured['fp32']['spill_bytes']}) are "
          f"distinct events ok")
    print(f"AC9 edge_accumulation puts the accumulator energy in the PE-array subtotal ok")
    print(f"AC10 spill/reload TIME follows the accumulator width: "
          f"{fp32['format_cycle_output']:.0f} (fp32) == 2x {fp16['format_cycle_output']:.0f} "
          f"(fp16) ok")
    print(f"AC11 cast TIME == {fp32['cast_bytes']} bytes x {CAST_CYCLE} = "
          f"{fp32['cast_cycle']:.0f} cycles, invariant to the accumulator format ok")
    print("AC12 each cycle knob scales only its own axis; undeclared means exactly 0 ok")
    print("AC13 energy and time are independent: neither knob moves the other's axis ok")
    print(f"AC14 retained mapping: {fp32_r['retained']} retained passes, reload "
          f"{fp32_r['reload']} B (pre-E20-4 charged {RETAINED_PREFIX_RELOAD} B) ok")
    print(f"AC15 retained/create identical across accumulator widths; reload and spill both 2x ok")
    print(f"AC16 create/retained/reload/spill are all non-zero and reload != spill ok")
    print("ALL ACCUMULATOR CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
