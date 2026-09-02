#!/usr/bin/env python3
"""Fold/setup dynamic-energy regression (E5).

The model already charges LATENCY for two control activities: a weight residency (fold) reloads
the array, and a per-layer schedule setup runs once. Both were RTL-calibrated for Gemmini. But
neither had an energy counterpart, so the same activity was free on the energy side apart from
the leakage its extra cycles happened to accrue -- the more precise the latency model became, the
more asymmetric energy got.

A config now declares the energy per event:
    weight_fold_fill_energy   per weight-element residency
    layer_setup_energy        once per layer
and the SELECTION RULE mirrors the latency side exactly: a calibrated per-element fold cost is
charged per per-PE weight element, the analytical drain is charged per tile boundary, and the two
are never charged together -- the same source that supplies the latency supplies the event count.

Two fixtures differ from configs/accelerators/gemmini.cfg only in those two energy knobs
(3.0/500.0 and twice that), plus the shipped config as the uncalibrated reference.

Checks (asserted; non-zero exit on failure):
  FS1  HAND IDENTITY: fold energy == fold events * the declared per-fold cost, in both runs.
  FS2  HAND IDENTITY: setup energy == setup events * the declared per-layer cost, in both runs.
  FS3  Doubling ONLY the energy knobs doubles both energies and leaves the COMPUTE SCHEDULE
       bit-identical. Energy and latency are separate axes: an energy calibration must not move
       the validated latency metric, and this is what proves it.
  FS4  The shipped config -- which calibrates the latency bubbles but not their energy -- reports
       the EVENT COUNTS with zero energy. That is the asymmetry made visible: it used to be
       impossible to tell "modeled as free" from "not calibrated".
  FS5  CONTAINMENT: the fold/setup energy is inside the layer total. Doubling the two knobs must
       raise the layer total energy by exactly the energy delta and by nothing else.
  FS6  RE3: the SETUP event count comes from the setup executing, not from its energy being
       priced. FS4 above only asserted that the shipped config's setup ENERGY is zero -- and it
       reported "0.00 over 0 setup event(s)" for a schedule that pays 2270 setup cycles, because
       the event was keyed off the unit energy. An uncalibrated setup must report its 1 event, say
       its unit cost is uncalibrated, and give the same count as the calibrated fixtures, which
       run the identical schedule.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKLOAD = ("gemm_64x64x64", "ws")
RUNS = {
    "gemmini_fold_energy": (3.0, 500.0),
    "gemmini_fold_energy_2x": (6.0, 1000.0),
}
UNCALIBRATED = "gemmini"


def run(target: str) -> str:
    layer = ROOT / "result" / target / WORKLOAD[0] / WORKLOAD[1] / "layer_0.txt"
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, WORKLOAD[0], WORKLOAD[1]],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if not layer.exists():
        raise SystemExit(f"missing simulator output: {layer}")
    return layer.read_text(encoding="utf-8")


def measure(text: str) -> dict:
    fold = re.search(r"Weight-fold energy\s*:\s*([\d.]+) \S+ over (\d+) fold event", text)
    setup = re.search(r"Layer-setup energy\s*:\s*([\d.]+) \S+ over (\d+) setup event\(s\)(.*)",
                      text)
    if fold is None or setup is None:
        raise AssertionError("the report does not carry the fold/setup energy axes (E5)")
    summary = text.split("Energy summary", 1)[1].split("MAC result", 1)[0]
    return {
        "fold_energy": float(fold.group(1)),
        "fold_events": int(fold.group(2)),
        "setup_energy": float(setup.group(1)),
        "setup_events": int(setup.group(2)),
        "setup_note": setup.group(3).strip(),
        "schedule": float(re.search(r"Compute-schedule latency\s*:\s*([\d.]+)", text).group(1)),
        "layer_total": float(re.search(
            r"Layer total energy\s*:\s*([\d.]+) \S+", summary).group(1)),
    }


def main() -> int:
    measured = {name: measure(run(name)) for name in RUNS}
    uncalibrated = measure(run(UNCALIBRATED))
    failures = []
    tol = 0.5

    for name, (fold_cost, setup_cost) in RUNS.items():
        state = measured[name]
        # FS1
        expected_fold = state["fold_events"]*fold_cost
        if abs(state["fold_energy"] - expected_fold) > tol:
            failures.append(f"FS1 {name}: fold energy {state['fold_energy']} != hand "
                            f"{expected_fold} ({state['fold_events']} events x {fold_cost})")
        # FS2
        expected_setup = state["setup_events"]*setup_cost
        if abs(state["setup_energy"] - expected_setup) > tol:
            failures.append(f"FS2 {name}: setup energy {state['setup_energy']} != hand "
                            f"{expected_setup} ({state['setup_events']} events x {setup_cost})")

    base, doubled = measured["gemmini_fold_energy"], measured["gemmini_fold_energy_2x"]
    # FS3
    if abs(doubled["fold_energy"] - 2*base["fold_energy"]) > tol:
        failures.append(f"FS3: doubling the fold cost gave {doubled['fold_energy']}, expected "
                        f"{2*base['fold_energy']}")
    if abs(doubled["setup_energy"] - 2*base["setup_energy"]) > tol:
        failures.append(f"FS3: doubling the setup cost gave {doubled['setup_energy']}, expected "
                        f"{2*base['setup_energy']}")
    if abs(doubled["schedule"] - base["schedule"]) > tol:
        failures.append(f"FS3: an ENERGY calibration moved the compute schedule "
                        f"({base['schedule']} -> {doubled['schedule']}); energy and latency must "
                        f"stay separate axes")
    if abs(base["fold_events"] - doubled["fold_events"]) > 0 or \
       abs(base["setup_events"] - doubled["setup_events"]) > 0:
        failures.append("FS3: the event counts moved with the unit costs; the identity would then "
                        "not isolate the unit cost")

    # FS4
    if uncalibrated["fold_events"] <= 0:
        failures.append(f"FS4 {UNCALIBRATED}: calibrates a fold latency bubble but reports no fold "
                        f"events, so the asymmetry is invisible")
    if uncalibrated["fold_energy"] != 0.0 or uncalibrated["setup_energy"] != 0.0:
        failures.append(f"FS4 {UNCALIBRATED}: declares no fold/setup energy but reports "
                        f"{uncalibrated['fold_energy']}/{uncalibrated['setup_energy']}")

    # FS6 -- RE3: the event count is a property of the schedule, not of the price list.
    if uncalibrated["setup_events"] != 1:
        failures.append(f"FS6 {UNCALIBRATED}: reports {uncalibrated['setup_events']} setup "
                        f"event(s) for a schedule that pays a setup latency; the count is still "
                        f"derived from the unit energy rather than from the setup executing")
    if "UNCALIBRATED" not in uncalibrated["setup_note"]:
        failures.append(f"FS6 {UNCALIBRATED}: a setup that runs but is not priced does not say its "
                        f"unit cost is uncalibrated ({uncalibrated['setup_note']!r}), so 'modeled "
                        f"as free' and 'not calibrated' are still indistinguishable")
    for name in RUNS:
        if measured[name]["setup_events"] != uncalibrated["setup_events"]:
            failures.append(f"FS6 {name}: {measured[name]['setup_events']} setup event(s) vs "
                            f"{uncalibrated['setup_events']} for the same schedule; pricing the "
                            f"setup must not change how many times it happens")

    # FS5
    energy_delta = ((doubled["fold_energy"] + doubled["setup_energy"]) -
                    (base["fold_energy"] + base["setup_energy"]))
    total_delta = doubled["layer_total"] - base["layer_total"]
    if abs(total_delta - energy_delta) > tol:
        failures.append(f"FS5: the layer total moved {total_delta} while the fold/setup energy "
                        f"moved {energy_delta}; the control energy is not contained in the total "
                        f"(or something else moved with it)")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"{'config':>26} {'fold E':>9} {'events':>7} {'setup E':>9} {'sched':>8} {'total E':>14}")
    for name in (UNCALIBRATED, "gemmini_fold_energy", "gemmini_fold_energy_2x"):
        state = uncalibrated if name == UNCALIBRATED else measured[name]
        print(f"{name:>26} {state['fold_energy']:>9.2f} {state['fold_events']:>7} "
              f"{state['setup_energy']:>9.2f} {state['schedule']:>8.0f} "
              f"{state['layer_total']:>14.2f}")
    print("FS1 fold energy == events x declared cost ok")
    print("FS2 setup energy == events x declared cost ok")
    print(f"FS3 2x the energy knobs -> 2x the energy, compute schedule unchanged at "
          f"{base['schedule']:.0f} ok")
    print(f"FS4 the shipped config reports {uncalibrated['fold_events']} fold events at zero "
          f"energy (uncalibrated, not free) ok")
    print(f"FS5 layer total moved by exactly the control-energy delta ({total_delta:.2f}) ok")
    print(f"FS6 an unpriced setup reports {uncalibrated['setup_events']} event(s), marked "
          f"UNCALIBRATED, matching the priced fixtures' count ok")
    print("ALL FOLD/SETUP ENERGY CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
