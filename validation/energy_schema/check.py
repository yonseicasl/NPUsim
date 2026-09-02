#!/usr/bin/env python3
"""Energy config SCHEMA and cost-state regression (RE3 / RE5 / RE8).

Two separate defects, both about what a number in the energy report MEANS.

RE3 -- the layer-setup EVENT count was derived from its unit energy. Gemmini declares
`layer_setup_cycle = 2270` and no setup energy, so the report said

    Layer-setup energy : 0.00 pJ over 0 setup event(s)

which reads as "no setup happened" for a schedule that pays 2270 cycles for one. The event now
comes from the setup executing; the energy is separately marked UNCALIBRATED when it is not priced.
The converse -- energy with no setup execution -- is rejected at config load, because it would put
the latency and energy event sources in different places.

RE5 -- `Uncosted components` was computed from the RESULT. A subtotal of 0 printed as "no energy
cost declared", conflating four situations a reader has to be able to tell apart:

    NOT MODELED    no energy key declared for the component at all
    PARTIAL        some costs declared, others missing -- the subtotal is an UNDERCOUNT, and this
                   is exactly what a key typo looks like
    modeled zero   every cost declared, all of them 0 -- a deliberate free
    NO ACTIVITY    every cost declared and non-zero, but the layer never exercised the component

The state is now derived from the DECLARATION (utils/energy_units.cc), so a typo is an error at
load time and a priced-but-idle component no longer reads as an unpriced one.

Fixtures (all gemmini.cfg with one thing changed):
  gemmini_cost_partial      mac_write_energy removed          -> MAC PARTIAL
  gemmini_cost_not_modeled  every [dram] energy key removed   -> DRAM NOT MODELED
  gemmini_cost_no_activity  bypass=1:1:1 on the GLB           -> GLB costed, no activity
  gemmini_setup_energy      layer_setup_energy = 7.5          -> a priced setup event
  gemmini                   the baseline                      -> modeled zeros, uncalibrated setup

Checks (asserted; non-zero exit on failure):
  ES1  A partially costed component is named as a PARTIAL UNDERCOUNT and names the missing key.
  ES2  A component with no declared cost is named NOT MODELED, distinctly from a modeled zero.
  ES3  A component whose costs are all declared 0 is named a modeled zero -- not "uncosted".
  ES4  A fully costed component with no activity says so, and is NOT counted as uncosted.
  ES5  `Uncosted components` counts DECLARATION gaps: the baseline has two all-zero rows and no
       warning, while the partial and not-modeled fixtures each raise exactly one.
  ES6  RE3: the baseline reports 1 setup event (not 0) with the energy marked UNCALIBRATED, and
       names the cycle cost that proves the setup runs.
  ES7  RE8 identity: 1 setup event x the declared 7.5 unit cost == the PE-array subtotal delta
       against the baseline, exactly once, and the same delta appears once in the layer total.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ROWS = ("MAC (compute+format)", "PE (local buffer)", "PE array", "Global buffer",
        "Multi-chip (NoP)", "DRAM")
FIXTURES = ("gemmini", "gemmini_cost_partial", "gemmini_cost_not_modeled",
            "gemmini_cost_no_activity", "gemmini_setup_energy")
SETUP_UNIT_ENERGY = 7.5      # gemmini_setup_energy.cfg: layer_setup_energy
SETUP_CYCLE = 2270           # gemmini.cfg: layer_setup_cycle


def run(target: str) -> str:
    layer = ROOT / "result" / target / "gemm_64x64x64" / "ws" / "layer_0.txt"
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, "gemm_64x64x64", "ws"],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if not layer.exists():
        raise SystemExit(f"missing simulator output: {layer}")
    return layer.read_text(encoding="utf-8")


def measure(text: str) -> dict:
    energy = text.split("Energy summary", 1)[1].split("Power summary", 1)[0]
    # The fold/setup/reduction event lines live in the MAC-result section, below the summary.
    events = text.split("MAC result", 1)[1].split("PE result", 1)[0]
    rows = {}
    for name in ROWS:
        found = re.search(rf" \* {re.escape(name)}\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)(.*)",
                          energy)
        if not found:
            raise SystemExit(f"missing energy row {name!r}")
        rows[name] = {"dynamic": float(found.group(1)), "static": float(found.group(2)),
                      "total": float(found.group(3)), "note": found.group(4).strip()}
    uncosted = re.search(r"Uncosted components\s*:\s*(\d+) of", energy)
    setup = re.search(r"Layer-setup energy\s*:\s*([\d.]+)\s+\S+\s+over\s+(\d+) setup event\(s\)(.*)",
                      events)
    return {"rows": rows,
            "uncosted": int(uncosted.group(1)) if uncosted else 0,
            "setup_energy": float(setup.group(1)),
            "setup_events": int(setup.group(2)),
            "setup_note": setup.group(3).strip(),
            "layer_total": float(re.search(r"Layer total energy\s*:\s*([\d.]+)", energy).group(1))}


def main() -> int:
    state = {name: measure(run(name)) for name in FIXTURES}
    failures = []
    base = state["gemmini"]

    # ES1
    mac = state["gemmini_cost_partial"]["rows"]["MAC (compute+format)"]["note"]
    if "PARTIAL" not in mac or "UNDERCOUNT" not in mac or "mac_write_energy" not in mac:
        failures.append(f"ES1: a half-priced MAC is not reported as a partial undercount naming "
                        f"the missing key ({mac!r})")
    # ES2
    dram = state["gemmini_cost_not_modeled"]["rows"]["DRAM"]["note"]
    if "NOT MODELED" not in dram:
        failures.append(f"ES2: a component with no declared energy key is not reported as NOT "
                        f"MODELED ({dram!r})")
    if dram == base["rows"]["Multi-chip (NoP)"]["note"]:
        failures.append("ES2: 'no cost declared' and 'declared zero' produce the same annotation, "
                        "so the report cannot distinguish them")
    # ES3
    for row in ("PE array", "Multi-chip (NoP)"):
        note = base["rows"][row]["note"]
        if "modeled zero" not in note:
            failures.append(f"ES3 {row}: an all-zero DECLARED cost is not reported as a modeled "
                            f"zero ({note!r})")
    # ES4
    glb = state["gemmini_cost_no_activity"]["rows"]["Global buffer"]
    if glb["total"] != 0.0:
        failures.append(f"ES4: the no-activity fixture still charges GLB energy ({glb['total']}); "
                        f"the fixture no longer isolates the state it is for")
    if "NO ACTIVITY" not in glb["note"]:
        failures.append(f"ES4: a fully costed but unexercised component is not reported as having "
                        f"no activity ({glb['note']!r})")
    if state["gemmini_cost_no_activity"]["uncosted"] != 0:
        failures.append(f"ES4: the no-activity fixture is counted as uncosted "
                        f"({state['gemmini_cost_no_activity']['uncosted']}); its costs ARE declared")
    # ES5
    zero_rows = sum(1 for row in ROWS if base["rows"][row]["total"] == 0.0)
    if zero_rows < 2:
        failures.append(f"ES5: the baseline no longer has >= 2 all-zero rows ({zero_rows}), so it "
                        f"cannot show that a zero subtotal is not an uncosted component")
    if base["uncosted"] != 0:
        failures.append(f"ES5: the baseline reports {base['uncosted']} uncosted component(s) even "
                        f"though every energy key it needs is declared -- the count is still being "
                        f"derived from the result rather than the declaration")
    for name in ("gemmini_cost_partial", "gemmini_cost_not_modeled"):
        if state[name]["uncosted"] != 1:
            failures.append(f"ES5 {name}: expected exactly 1 uncosted component, got "
                            f"{state[name]['uncosted']}")
    # ES6
    if base["setup_events"] != 1:
        failures.append(f"ES6: the baseline reports {base['setup_events']} setup event(s) for a "
                        f"schedule that pays {SETUP_CYCLE} setup cycles; the event count is still "
                        f"keyed off the unit energy")
    if base["setup_energy"] != 0.0:
        failures.append(f"ES6: the baseline declares no layer_setup_energy but reports "
                        f"{base['setup_energy']}")
    if "UNCALIBRATED" not in base["setup_note"] or str(SETUP_CYCLE) not in base["setup_note"]:
        failures.append(f"ES6: an unpriced setup does not say its unit cost is uncalibrated and "
                        f"name the {SETUP_CYCLE}-cycle execution that proves it runs "
                        f"({base['setup_note']!r})")
    # ES7
    priced = state["gemmini_setup_energy"]
    if priced["setup_events"] != 1:
        failures.append(f"ES7: expected 1 setup event, got {priced['setup_events']}")
    expected = priced["setup_events"]*SETUP_UNIT_ENERGY
    if abs(priced["setup_energy"] - expected) > 1e-9:
        failures.append(f"ES7: setup energy {priced['setup_energy']} != events x unit cost "
                        f"{expected}")
    array_delta = (priced["rows"]["PE array"]["total"] - base["rows"]["PE array"]["total"])
    if abs(array_delta - expected) > 1e-6:
        failures.append(f"ES7: the PE-array subtotal moved by {array_delta}, not by the "
                        f"{expected} the setup event costs -- the charge is missing from, or "
                        f"double counted in, its owning component")
    total_delta = priced["layer_total"] - base["layer_total"]
    if abs(total_delta - expected) > 1e-6:
        failures.append(f"ES7: the layer total moved by {total_delta}, not {expected}; the setup "
                        f"charge must appear in it exactly once")
    for row in ROWS:
        if row == "PE array":
            continue
        other = priced["rows"][row]["total"] - base["rows"][row]["total"]
        if abs(other) > 1e-6:
            failures.append(f"ES7: pricing the setup also moved the {row} subtotal by {other}")
    if priced["rows"]["PE array"]["note"]:
        failures.append(f"ES7: a PE array priced only through its optional setup cost is annotated "
                        f"{priced['rows']['PE array']['note']!r}; an optional-feature cost still "
                        f"makes the component costed")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"{'fixture':>26} {'uncosted':>9} {'setup ev':>9} {'setup E':>9}  annotated rows")
    for name in FIXTURES:
        notes = [f"{row.split(' ')[0]}={state[name]['rows'][row]['note'].split(':')[0][3:]}"
                 for row in ROWS if state[name]["rows"][row]["note"]]
        print(f"{name:>26} {state[name]['uncosted']:>9} {state[name]['setup_events']:>9} "
              f"{state[name]['setup_energy']:>9.2f}  {'; '.join(notes)}")
    print("ES1 a half-priced component is a PARTIAL UNDERCOUNT naming the missing key ok")
    print("ES2 a component with no declared cost is NOT MODELED, distinct from a modeled zero ok")
    print("ES3 an all-zero declared cost is a modeled zero, not an uncosted component ok")
    print("ES4 a costed but unexercised component says NO ACTIVITY and is not counted uncosted ok")
    print("ES5 Uncosted components counts declaration gaps (baseline 0 with 2 zero rows) ok")
    print(f"ES6 an unpriced setup reports 1 event over {SETUP_CYCLE} cycles, marked UNCALIBRATED ok")
    print(f"ES7 1 setup event x {SETUP_UNIT_ENERGY} == the PE-array delta == the layer-total delta,"
          f" and no other subtotal moved ok")
    print("ALL ENERGY SCHEMA CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
