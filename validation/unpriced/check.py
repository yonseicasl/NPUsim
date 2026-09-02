#!/usr/bin/env python3
"""Active-event cost completeness (E20-1 / E20-2).

An energy axis reading 0.00 has four possible meanings and the report used to conflate two of them:
the event was modeled as free, or the event HAPPENED and nobody declared its price. RE5 separated
the config-side states (not modeled / partial / modeled zero / calibrated); this gate covers the
run-side one, which needs an event COUNT next to the declaration because both cases print 0.

`gemmini_cacti22` is what forced it. Its SRAM costs come from CACTI and its DRAM costs from
DRAMsim3, so it satisfied every earlier condition for absolute energy and reported 418.685 mW --
while the same report showed a layer setup, 1,048,560 bytes of accumulator reload, 1,048,576 of
spill and 4,096 of final cast, every one of them at 0.00 pJ because those keys were never declared.
The arithmetic was right; the numerator was missing charges of undeclared size.

The rule is uniform: an event is ACTIVE when its own counter is non-zero, and PRICED when the config
declares the key with ANY value, zero included. A declared zero is a statement ("this design
attributes no cost here"); an absent key is not a statement at all. Three activity counters had to
be added to make this checkable -- DRAM row activations, Format-IP payload/metadata transactions and
array reduction additions had energy axes but no counts.

Checks (asserted; non-zero exit on failure):
  UP1  The complete fixture reports ZERO unpriced active events and receives a wattage. Without this
       the rest would be satisfied by a gate that refuses everything.
  UP2  Removing ONE priced key at a time from the complete fixture is detected every time, and named
       -- so no active event is left without its own negative test.
  UP3  An INACTIVE feature's absent cost is NOT counted. Dense int8 does no Format-IP metadata work
       and a bus array has no reduction tree; flagging those would make the check useless noise.
  UP4  An explicit ZERO clears the flag while contributing no energy. This is the distinction the
       whole gate rests on: `key = 0` is a modeled zero, an absent key is an unpriced charge.
  UP5  gemmini_cacti22 -- externally calibrated on SRAM and DRAM but incomplete on the rest -- is
       refused a wattage and lists what is missing.
"""
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import perturb_lib as pl                                            # noqa: E402

COMPLETE = "gemmini_power_calibrated"
INCOMPLETE = "gemmini_cacti22"
WORKLOAD, MAPPING = "gemm_64x64x64", "ws"
# Every key the complete fixture declares for an ACTIVE event, with the section it lives in.
PRICED = (("systolic_array", "layer_setup_energy"),
          ("systolic_array", "weight_fold_fill_energy"),
          ("systolic_array", "accumulator_spill_energy"),
          ("systolic_array", "format_payload_energy"),
          ("systolic_array", "lb_static_energy"),
          ("multi_chip", "output_cast_energy"))
# Keys whose feature is INACTIVE in that fixture, so their absence must not be flagged.
INACTIVE = ("format_metadata_energy", "adder_energy", "row_miss_energy")


def measure(target, edits=None, label="base"):
    if edits is None:
        text = pl.run(target, WORKLOAD, MAPPING)
    else:
        name = f"__up_{label}"
        pl.variant(target, name, edits)
        shutil.rmtree(pl.MAPPINGS / name, ignore_errors=True)
        shutil.copytree(pl.MAPPINGS / target, pl.MAPPINGS / name)
        text = pl.run(name, WORKLOAD, MAPPING)
        pl.discard(name)
    if text is None:
        return None
    summary = text.split("Energy summary", 1)[1].split("MAC result", 1)[0]
    unpriced = re.findall(r"^  \* (.*) fired but \[(\w+)\] (\w+) is not declared$", summary, re.M)
    power = re.search(r"Average total power\s*:\s*([\d.]+) mW", summary)
    return {"unpriced": unpriced,
            "power_mW": float(power.group(1)) if power else None,
            "unit_note": re.search(r"Energy unit\s*:(.*)", summary).group(1).strip(),
            "note": re.search(r"Average power\s*:(.*)", summary).group(1).strip()
                    if re.search(r"Average power\s*:(.*)", summary) else "",
            "estimated": "ESTIMATED" in re.search(r"Layer total energy\s*:(.*)",
                                                  summary).group(1),
            "total": float(re.search(r"Layer total energy\s*:\s*([\d.]+)", summary).group(1))}


def drop(section, key):
    def apply(text):
        start = text.index(f"[{section}]")
        head, tail = text[:start], text[start:]
        return head + re.sub(rf"^{re.escape(key)}\s*=.*$", f"#{key} removed", tail,
                             count=1, flags=re.M)
    return apply


def set_zero(section, key, fields):
    def apply(text):
        start = text.index(f"[{section}]")
        head, tail = text[:start], text[start:]
        zero = ":".join(["0"]*fields)
        return head + re.sub(rf"^{re.escape(key)}\s*=.*$", f"{key} = {zero}", tail,
                             count=1, flags=re.M)
    return apply


def main() -> int:
    failures = []
    base = measure(COMPLETE)
    if base is None:
        raise SystemExit(f"{COMPLETE} does not run")

    # UP1
    if base["unpriced"]:
        failures.append(f"UP1 {COMPLETE}: {len(base['unpriced'])} active event(s) still unpriced "
                        f"({[u[0] for u in base['unpriced']]}); the fixture that proves the positive "
                        f"path must itself be complete")
    if base["power_mW"] is None:
        failures.append(f"UP1 {COMPLETE}: receives no wattage ({base['note']}), so a gate that "
                        f"refuses everything would pass the checks below")
    if base["estimated"]:
        failures.append(f"UP1 {COMPLETE}: its total is marked ESTIMATED although every active event "
                        f"is priced and every component cost declared")

    # UP2
    for index, (section, key) in enumerate(PRICED):
        state = measure(COMPLETE, drop(section, key), f"drop{index}")
        if state is None:
            failures.append(f"UP2 {section}.{key}: the variant does not run")
            continue
        found = [u for u in state["unpriced"] if u[2] == key]
        if not found:
            failures.append(f"UP2 {section}.{key}: removing it was NOT detected. Its event fires in "
                            f"this fixture, so an absent key leaves the total short by an "
                            f"undeclared amount and nothing says so "
                            f"(reported: {[u[2] for u in state['unpriced']]})")
        if state["power_mW"] is not None:
            failures.append(f"UP2 {section}.{key}: a wattage is still reported "
                            f"({state['power_mW']} mW) with that cost absent")
        if not state["estimated"]:
            failures.append(f"UP2 {section}.{key}: the total is not marked ESTIMATED with that cost "
                            f"absent")

    # UP3
    for key in INACTIVE:
        if any(u[2] == key for u in base["unpriced"]):
            failures.append(f"UP3 {key}: flagged as unpriced although its feature is inactive in "
                            f"{COMPLETE}; an inactive feature needs no cost")

    # UP4 -- an explicit zero is a statement, an absent key is not
    zeroed = measure(COMPLETE, set_zero("systolic_array", "layer_setup_energy", 1), "zero")
    removed = measure(COMPLETE, drop("systolic_array", "layer_setup_energy"), "gone")
    if zeroed is None or removed is None:
        failures.append("UP4: the explicit-zero or removed variant does not run")
    else:
        if any(u[2] == "layer_setup_energy" for u in zeroed["unpriced"]):
            failures.append("UP4: `layer_setup_energy = 0` is reported as unpriced; an explicit "
                            "zero is a modeled-zero STATEMENT and must clear the flag")
        if zeroed["power_mW"] is None:
            failures.append(f"UP4: declaring the cost as 0 still refuses a wattage "
                            f"({zeroed['note']})")
        if not any(u[2] == "layer_setup_energy" for u in removed["unpriced"]):
            failures.append("UP4: removing the key is not detected, so `= 0` and absent are still "
                            "the same state")
        if abs(zeroed["total"] - removed["total"]) > 0.5:
            failures.append(f"UP4: the explicit zero changed the total ({removed['total']} -> "
                            f"{zeroed['total']}); a modeled zero must cost nothing, it only states "
                            f"that it costs nothing")

    # UP5
    incomplete = measure(INCOMPLETE)
    if incomplete is None:
        failures.append(f"UP5 {INCOMPLETE} does not run")
    else:
        if not incomplete["unpriced"]:
            failures.append(f"UP5 {INCOMPLETE}: reports no unpriced event, but its MAC, "
                            f"accumulator, cast and control costs are not externally calibrated -- "
                            f"see validation/phase4/PROVENANCE.md")
        if incomplete["power_mW"] is not None:
            failures.append(f"UP5 {INCOMPLETE}: reports {incomplete['power_mW']} mW while its own "
                            f"provenance says several axes are uncalibrated")
        if "totals and power are meaningful" in incomplete["unit_note"]:
            failures.append(f"UP5 {INCOMPLETE}: the pre-run unit description still claims totals "
                            f"and power are meaningful even though the run finds unpriced events "
                            f"({incomplete['unit_note']!r})")
        if not incomplete["estimated"]:
            failures.append(f"UP5 {INCOMPLETE}: its total is not marked ESTIMATED")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"{COMPLETE:>28}: {len(base['unpriced'])} unpriced, {base['power_mW']} mW")
    print(f"{INCOMPLETE:>28}: {len(incomplete['unpriced'])} unpriced, wattage refused")
    for event, section, key in incomplete["unpriced"]:
        print(f"{'':>30}  {event} -> [{section}] {key}")
    print("UP1 the complete fixture has no unpriced event and does receive a wattage ok")
    print(f"UP2 removing each of the {len(PRICED)} priced keys is detected and refuses power ok")
    print(f"UP3 inactive features ({', '.join(INACTIVE)}) are not flagged ok")
    print("UP4 an explicit 0 clears the flag and costs nothing; an absent key does not ok")
    print(f"UP5 {INCOMPLETE} is refused a wattage and names what is missing ok")
    print("ALL UNPRICED-EVENT CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
