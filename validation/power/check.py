#!/usr/bin/env python3
"""Average-power regression (E10 / Phase 5).

Power was unsupported: nothing converted cycles to seconds, so total energy and wall-clock were
never related. The contract added is deliberately small and fully stated rather than large and
partly founded:

    time  = critical-path cycles / authoritative frequency
    power = energy / time
    EDP   = energy x time      ED2P = energy x time^2

with three preconditions the report states rather than assumes:
  * ONE CLOCK -- the timeline is a single shared cycle axis, so cycles convert to seconds only if
    every modeled component runs on the same clock. A mixed-domain config gets "unsupported", not
    a number divided by an arbitrary domain.
  * NOTHING MISSING FROM THE NUMERATOR (E20-1) -- a component whose cost is only half declared, or
    an event that fired with no declared price, makes the total an undercount of unstated size.
    energy/time is then arithmetic, not a wattage. This is what `pJ + provenance + a precision MAC
    key` alone used to wave through: gemmini_cacti22 has externally derived SRAM and DRAM costs and
    reported 418.685 mW while its setup, accumulator and cast events were entirely unpriced.
  * ABSOLUTE ENERGY, QUALIFIED (RE2) -- watts require more than pJ-shaped arithmetic. The fixture
    must DECLARE `energy_unit = pJ`, declare WHERE that scale came from (`energy_reference`), and
    have a compute cost calibrated for the operand precision actually running. A config that
    declares pJ with no provenance, or that falls back to one MAC scalar shared by every
    precision, is reported as unsupported -- not as milliwatts.
  * AVERAGE ONLY -- peak power needs per-event concurrency, which this model does not resolve, so
    it is explicitly unsupported instead of approximated.

Fixtures: gemmini_power_calibrated (declared pJ + declared provenance + a MAC energy for the
int8 x int8 operand precision + non-zero leakage + single clock) and its slow-DRAM twin (the same
energy over a 10x longer window), plus three that must NOT get watts -- eyeriss_energy
(normalized), gemmini_power_no_provenance (pJ, no reference) and gemmini_power_precision_fallback
(pJ + reference, uncalibrated precision).

Checks (asserted; non-zero exit on failure):
  PW1  HAND IDENTITY, time: elapsed time == critical-path cycles / authoritative frequency.
  PW2  HAND IDENTITY, power: dynamic, leakage and total average power each equal the corresponding
       energy divided by that time, and total == dynamic + leakage.
  PW3  HAND IDENTITY, EDP/ED2P == total energy x time and x time^2.
  PW4  A LATENCY-ONLY change moves power the right way: the slow-DRAM run has the same dynamic
       energy (LK2) over a longer window, so its average dynamic power must fall by exactly the
       latency ratio -- while its LEAKAGE power stays constant, because leakage energy grows with
       the same window it is divided by. Getting both from one config change is what shows the
       time base is applied consistently to the two axes.
  PW5  A normalized-energy fixture reports power as UNSUPPORTED and says why. Watts from
       normalized costs would be meaningless.
  PW7  RE2 QUALIFICATION: each of the three non-qualifying fixtures reports NO wattage and names
       its specific disqualifying reason (no unit / no provenance / precision fallback), and its
       energy total is marked ESTIMATED. Only the fully calibrated fixture gets milliwatts. The
       precision cases are two: no per-precision MAC energy at all, and (RE4) one declared for the
       operands but not for the accumulator width -- both leave the compute term uncalibrated.
  PW8  E20-5: no SHIPPED architecture config receives a wattage. None of them declares an
       `energy_unit`, so their totals are relative event x unit-cost sums of unknown scale -- and
       the Phase-4 comparison against CACTI and DRAMsim3 showed their unit costs are not
       reproducible from the geometry and device they themselves declare. Refusing them is the
       correct behaviour and is checked here so it cannot regress quietly the way it would if it
       were only a property of the current configs.
  PW6  The report states what power excludes (DRAM background/refresh, clock network,
       controller/DMA, PHY) and that PEAK power is unsupported -- so a core-datapath figure is not
       mistaken for a chip-level one.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ABSOLUTE = {"gemmini_power_calibrated": ("gemm_64x64x64", "ws"),
            "gemmini_power_calibrated_slow_dram": ("gemm_64x64x64", "ws")}
# RE2: fixtures that must be refused watts, each with the substring its reason must contain.
UNQUALIFIED = {"eyeriss_energy": (("alexnet", "silicon"), "relative to the declared reference"),
               "gemmini_power_no_provenance": (("gemm_64x64x64", "ws"),
                                               "energy_reference is not"),
               "gemmini_power_precision_fallback": (("gemm_64x64x64", "ws"),
                                                    "not calibrated for the precision in use"),
               # RE4: operand widths priced, accumulator width not -- a PARTIAL calibration, which
               # is still not a calibrated compute cost.
               "gemmini_macacc_operand_only": (("gemm_64x64x64", "ws"),
                                               "not calibrated for the precision in use"),
               "gemmini": (("gemm_64x64x64", "ws"), "no energy_unit declared"),
               # E20-1: complete on EVENTS but a whole component's cost absent.
               "gemmini_power_missing_dram": (("gemm_64x64x64", "ws"),
                                              "1 component(s) with a missing cost"),
               # E20-1: externally calibrated on SRAM and DRAM, but its MAC, accumulator, cast and
               # control costs are not -- so several active events fire unpriced.
               "gemmini_cacti22": (("gemm_64x64x64", "ws"),
                                   "active event(s) with no declared price")}


def run(target: str, workload: str, mapping: str) -> str:
    layer = ROOT / "result" / target / workload / mapping / "layer_0.txt"
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, workload, mapping],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if not layer.exists():
        raise SystemExit(f"missing simulator output: {layer}")
    return layer.read_text(encoding="utf-8")


def measure(text: str) -> dict:
    energy = text.split("Energy summary", 1)[1].split("Power summary", 1)[0]
    power = text.split("Power summary", 1)[1].split("MAC result", 1)[0]

    def optional(label, source):
        found = re.search(rf"{re.escape(label)}\s*:\s*([\d.]+)", source)
        return float(found.group(1)) if found else None
    return {
        "dynamic_energy": optional("Layer dynamic energy", energy),
        "static_energy": optional("Layer static energy", energy),
        "total_energy": optional("Layer total energy", energy),
        "clock_note": re.search(r"Clock domain\s*:(.*)", power).group(1).strip(),
        "frequency": optional("Authoritative clock", power),
        "elapsed_ms": optional("Elapsed time", power),
        "dynamic_power": optional("Average dynamic power", power),
        "leakage_power": optional("Average leakage power", power),
        "total_power": optional("Average total power", power),
        "edp": optional("EDP", power),
        "ed2p": optional("ED2P", power),
        "power_note": (re.search(r"Average power\s*:(.*)", power).group(1).strip()
                       if re.search(r"Average power\s*:(.*)", power) else ""),
        "not_included": (re.search(r"Not included\s*:(.*)", power).group(1).strip()
                         if re.search(r"Not included\s*:(.*)", power) else ""),
        "peak": (re.search(r"Peak power\s*:(.*)", power).group(1).strip()
                 if re.search(r"Peak power\s*:(.*)", power) else ""),
        "critical": float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", text).group(1)),
        # RE2: the energy total must say so when its scale is not established.
        "total_estimated": "ESTIMATED" in (
            re.search(r"Layer total energy\s*:(.*)", energy).group(1)),
        "unit_note": re.search(r"Energy unit\s*:(.*)", energy).group(1).strip(),
    }


def main() -> int:
    measured = {name: measure(run(name, *workload)) for name, workload in ABSOLUTE.items()}
    unqualified = {name: (measure(run(name, *workload)), reason)
                   for name, (workload, reason) in UNQUALIFIED.items()}
    normalized = unqualified["eyeriss_energy"][0]
    failures = []

    # RE2: a qualified fixture must not be silently marked estimated, or the gate below would be
    # comparing two "estimated" reports and proving nothing.
    for name, state in measured.items():
        if state["total_estimated"]:
            failures.append(f"RE2 {name}: a fully calibrated fixture marks its total ESTIMATED "
                            f"({state['unit_note']!r})")

    for name, state in measured.items():
        if state["dynamic_power"] is None:
            failures.append(f"{name}: no average power reported ({state['power_note']})")
            continue
        seconds = state["elapsed_ms"]/1.0e3
        # PW1
        expected_seconds = state["critical"]/(state["frequency"]*1.0e6)
        if abs(seconds - expected_seconds) > max(1e-12, expected_seconds*1e-6):
            failures.append(f"PW1 {name}: elapsed {seconds}s != cycles/frequency "
                            f"{expected_seconds}s ({state['critical']} / "
                            f"{state['frequency']} MHz)")
        # PW2 -- energy in pJ over seconds, reported in mW
        for label, energy_key, power_key in (("dynamic", "dynamic_energy", "dynamic_power"),
                                             ("leakage", "static_energy", "leakage_power"),
                                             ("total", "total_energy", "total_power")):
            expected = state[energy_key]*1.0e-9/seconds
            if abs(state[power_key] - expected) > max(1e-3, abs(expected)*1e-5):
                failures.append(f"PW2 {name} {label}: {state[power_key]} mW != energy/time "
                                f"{expected} mW ({state[energy_key]} pJ / {seconds} s)")
        if abs(state["total_power"] - (state["dynamic_power"] + state["leakage_power"])) > 1e-2:
            failures.append(f"PW2 {name}: total power {state['total_power']} != dynamic + leakage "
                            f"{state['dynamic_power'] + state['leakage_power']}")
        # PW3
        for label, key, exponent in (("EDP", "edp", 1), ("ED2P", "ed2p", 2)):
            expected = state["total_energy"]*seconds**exponent
            if abs(state[key] - expected) > max(1e-6, abs(expected)*1e-5):
                failures.append(f"PW3 {name} {label}: {state[key]} != total energy x time^"
                                f"{exponent} {expected}")

    base = measured["gemmini_power_calibrated"]
    slow = measured["gemmini_power_calibrated_slow_dram"]
    if base["dynamic_power"] is not None and slow["dynamic_power"] is not None:
        latency_ratio = slow["critical"]/base["critical"]
        # PW4: same dynamic energy over a longer window -> power falls by the latency ratio.
        expected_dynamic = base["dynamic_power"]/latency_ratio
        if abs(slow["dynamic_power"] - expected_dynamic) > max(1e-3, expected_dynamic*1e-4):
            failures.append(f"PW4: dynamic power {slow['dynamic_power']} != "
                            f"{expected_dynamic} (the same dynamic energy over a "
                            f"{latency_ratio:.4f}x longer window)")
        # ... while leakage energy grows with the very window it is divided by, so leakage power
        # is invariant.
        if abs(slow["leakage_power"] - base["leakage_power"]) > max(1e-3,
                                                                   base["leakage_power"]*1e-4):
            failures.append(f"PW4: leakage power moved ({base['leakage_power']} -> "
                            f"{slow['leakage_power']}); leakage energy scales with the same window "
                            f"it is divided by, so its power must be invariant")

    # PW5
    if normalized["dynamic_power"] is not None:
        failures.append(f"PW5: the normalized fixture reports {normalized['dynamic_power']} mW; "
                        f"watts from normalized costs are meaningless")
    if "unsupported" not in normalized["power_note"] or \
       "normalized" not in normalized["power_note"]:
        failures.append(f"PW5: the normalized fixture does not explain why power is unavailable "
                        f"({normalized['power_note']!r})")

    # PW7 -- every non-qualifying fixture: no watts, the specific reason, ESTIMATED total.
    for name, (state, reason) in unqualified.items():
        if state["dynamic_power"] is not None:
            failures.append(f"PW7 {name}: reports {state['dynamic_power']} mW without qualifying "
                            f"for absolute power ({reason})")
        if "unsupported" not in state["power_note"] or reason not in state["power_note"]:
            failures.append(f"PW7 {name}: power note does not state the disqualifying reason "
                            f"{reason!r} ({state['power_note']!r})")
        if not state["total_estimated"]:
            failures.append(f"PW7 {name}: the energy total is not marked ESTIMATED even though "
                            f"its scale is not established ({reason})")

    # PW8 -- every shipped architecture config must be refused
    SHIPPED = {"gemmini": ("gemm_64x64x64", "ws"), "eyeriss": ("alexnet", "silicon"),
               "eyerissv2": ("gemm_64x64x64", "ws"), "tpu": ("gemm_64x64x64", "ws"),
               "tpuv3": ("gemm_64x64x64", "ws"), "simba": ("gemm_64x64x64", "ws"),
               "maeri": ("gemm_64x64x64", "ws"), "fsd": ("gemm_64x64x64", "ws")}
    for name, workload in SHIPPED.items():
        state = measure(run(name, *workload))
        if state["dynamic_power"] is not None or state["total_power"] is not None:
            failures.append(f"PW8 {name}: a shipped config reports a wattage "
                            f"({state['total_power']} mW). It declares no energy_unit, and Phase 4 "
                            f"shows its unit costs do not reproduce from its own declared geometry "
                            f"or DRAM device, so an absolute power from it is not meaningful")
        if not state["total_estimated"]:
            failures.append(f"PW8 {name}: its total is not marked ESTIMATED")

    # PW6
    for name, state in list(measured.items()) + [(n, s) for n, (s, _) in unqualified.items()]:
        for phrase in ("background", "clock network", "controller"):
            if phrase not in state["not_included"]:
                failures.append(f"PW6 {name}: the excluded-power scope does not mention "
                                f"{phrase!r} ({state['not_included']!r})")
        if "unsupported" not in state["peak"]:
            failures.append(f"PW6 {name}: peak power is not declared unsupported "
                            f"({state['peak']!r})")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"{'config':>28} {'cycles':>9} {'ms':>10} {'dyn mW':>9} {'leak mW':>9} {'total mW':>9}")
    for name, state in measured.items():
        print(f"{name:>28} {state['critical']:>9.0f} {state['elapsed_ms']:>10.4f} "
              f"{state['dynamic_power']:>9.3f} {state['leakage_power']:>9.3f} "
              f"{state['total_power']:>9.3f}")
    print("PW1 elapsed time == cycles / authoritative frequency ok")
    print("PW2 dynamic/leakage/total power == energy / time; total == dynamic + leakage ok")
    print("PW3 EDP/ED2P == total energy x time and x time^2 ok")
    print(f"PW4 a {slow['critical']/base['critical']:.2f}x longer window scales dynamic power down "
          f"by exactly that factor and leaves leakage power invariant ok")
    print("PW5 normalized-energy fixture reports power unsupported, with the reason ok")
    for name, (state, reason) in unqualified.items():
        print(f"PW7 {name:>32}: no watts; {reason} ok")
    print(f"PW8 all {len(SHIPPED)} shipped architecture configs are refused a wattage and marked "
          f"ESTIMATED ok")
    print("PW6 excluded scope stated; peak power declared unsupported ok")
    print("ALL POWER CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
