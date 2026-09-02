#!/usr/bin/env python3
"""Lane-reduction energy regression (E9).

The lane->accumulator adder tree has an energy model, but EVERY shipped config leaves
`mac_width` at 1 -- a single-lane MAC has no tree and no reduction work -- and none declares
`mac_reduction_energy`. So the knob had no regression whatsoever: it could have been off by any
factor, or wired to the wrong lane count, without a single check noticing. (E4 split the reduction
onto its own reported axis, which is what makes this checkable at all.)

Two fixtures give the MAC 4 lanes (number_of_macs = 1, mac_width = 4) with 4 reduction lanes live
at the MAC level, and declare mac_reduction_energy = 0.4 and 0.8. The mapping moves C from PE_Y to
the MAC level so the TOTAL work is identical to configs/mappings/gemmini.

Checks (asserted; non-zero exit on failure):
  RD1  The axis is actually non-zero -- i.e. the fixture really does exercise a multi-lane tree.
  RD2  HAND IDENTITY: an N-leaf binary adder tree performs N-1 additions, once per issue, so
         reduction energy == (scalar MACs / active lanes) * (active lanes - 1) * unit cost
       With 262144 scalar MACs over 4 lanes that is 65536 issues x 3 additions x the unit cost.
       This is the identity that distinguishes "N-1 additions" from "N" or from "one tree per
       scalar MAC", none of which anything checked before.
  RD3  Doubling ONLY the reduction unit cost doubles ONLY the reduction energy: the MAC count,
       the computation energy and the compute schedule must all stay bit-identical. Reduction is
       its own axis, not a rescaling of the MAC total.
  RD4  The shipped single-lane config has no reduction work at all (mac_width = 1), so its axis
       must be exactly 0 -- a one-leaf tree adds nothing. If it were non-zero the model would be
       charging a tree that does not exist.
  RD5  RE8 CONTAINMENT: reporting reduction on its own axis is only half the contract -- the
       energy also has to reach the totals exactly once. Doubling the unit cost must move the MAC
       component subtotal AND the layer total by precisely the same delta as the axis itself, and
       move no other component subtotal at all. Without this, an axis could be reported correctly
       and still be missing from (or double counted in) the total that gets published.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKLOAD = ("gemm_64x64x64", "ws")
# fixture -> declared mac_reduction_energy
RUNS = {"gemmini_reduction_energy": 0.4, "gemmini_reduction_energy_2x": 0.8}
ACTIVE_LANES = 4          # mac_width = 4 with C = 4 live at the MAC level
SINGLE_LANE = "gemmini"   # shipped config: mac_width = 1, no tree


def run(target: str) -> str:
    layer = ROOT / "result" / target / WORKLOAD[0] / WORKLOAD[1] / "layer_0.txt"
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, WORKLOAD[0], WORKLOAD[1]],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if not layer.exists():
        raise SystemExit(f"missing simulator output: {layer}")
    return layer.read_text(encoding="utf-8")


def measure(text: str) -> dict:
    mac = text.split("MAC result", 1)[1].split("PE result", 1)[0]
    return {
        "macs": int(re.search(r"# of computations\s*:\s*(\d+)", mac).group(1)),
        "reduction": float(re.search(
            r"Reduction energy\s*:\s*([\d.]+) (?:pJ|normalized|uncalibrated)", mac).group(1)),
        "computation": float(re.search(
            r"Computation energy\s*:\s*([\d.]+) (?:pJ|normalized|uncalibrated)", mac).group(1)),
        "schedule": float(re.search(r"Compute-schedule latency\s*:\s*([\d.]+)", text).group(1)),
        # RE8: the component subtotals and the layer total the axis has to land in.
        "rows": component_rows(text),
        "layer_total": float(re.search(r"Layer total energy\s*:\s*([\d.]+)", text).group(1)),
    }


ROWS = ("MAC (compute+format)", "PE (local buffer)", "PE array", "Global buffer",
        "Multi-chip (NoP)", "DRAM")


def component_rows(text: str) -> dict:
    energy = text.split("Energy summary", 1)[1].split("Power summary", 1)[0]
    rows = {}
    for name in ROWS:
        found = re.search(rf" \* {re.escape(name)}\s+[\d.]+\s+[\d.]+\s+([\d.]+)", energy)
        if not found:
            raise SystemExit(f"missing energy row {name!r}")
        rows[name] = float(found.group(1))
    return rows


def main() -> int:
    measured = {name: measure(run(name)) for name in RUNS}
    single = measure(run(SINGLE_LANE))
    failures = []
    tol = 0.5

    for name, unit in RUNS.items():
        state = measured[name]
        # RD1
        if state["reduction"] <= 0.0:
            failures.append(f"RD1 {name}: reduction energy is {state['reduction']}, so the fixture "
                            f"is not exercising a multi-lane adder tree")
            continue
        # RD2
        issues = state["macs"]/ACTIVE_LANES
        expected = issues*(ACTIVE_LANES - 1)*unit
        if abs(state["reduction"] - expected) > tol:
            failures.append(f"RD2 {name}: reduction energy {state['reduction']} != hand {expected} "
                            f"({state['macs']} MACs / {ACTIVE_LANES} lanes = {issues:.0f} issues x "
                            f"{ACTIVE_LANES - 1} additions x {unit})")

    base, doubled = measured["gemmini_reduction_energy"], measured["gemmini_reduction_energy_2x"]
    # RD3
    if abs(doubled["reduction"] - 2*base["reduction"]) > tol:
        failures.append(f"RD3: doubling the reduction unit cost gave {doubled['reduction']}, "
                        f"expected {2*base['reduction']}")
    for field in ("macs", "computation", "schedule"):
        if abs(doubled[field] - base[field]) > tol:
            failures.append(f"RD3: {field} moved with the reduction unit cost "
                            f"({base[field]} -> {doubled[field]}); reduction must be its own axis")

    # RD5 -- RE8 containment: the same delta, in the axis, in its component, in the layer total.
    axis_delta = doubled["reduction"] - base["reduction"]
    mac_delta = doubled["rows"]["MAC (compute+format)"] - base["rows"]["MAC (compute+format)"]
    if abs(mac_delta - axis_delta) > 1e-6:
        failures.append(f"RD5: the reduction axis moved by {axis_delta} but the MAC component "
                        f"subtotal moved by {mac_delta}; the axis is not contained in its own "
                        f"component exactly once")
    total_delta = doubled["layer_total"] - base["layer_total"]
    if abs(total_delta - axis_delta) > 1e-6:
        failures.append(f"RD5: the reduction axis moved by {axis_delta} but the layer total moved "
                        f"by {total_delta}; the axis must reach the published total exactly once")
    for row in ROWS:
        if row == "MAC (compute+format)":
            continue
        other = doubled["rows"][row] - base["rows"][row]
        if abs(other) > 1e-6:
            failures.append(f"RD5: doubling the reduction unit cost also moved the {row} subtotal "
                            f"by {other}")

    # RD4
    if single["reduction"] != 0.0:
        failures.append(f"RD4 {SINGLE_LANE}: a single-lane MAC has no adder tree but reports "
                        f"{single['reduction']} of reduction energy")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"{'config':>30} {'reduction E':>12} {'compute E':>12} {'MACs':>9} {'sched':>7}")
    for name in (SINGLE_LANE, "gemmini_reduction_energy", "gemmini_reduction_energy_2x"):
        state = single if name == SINGLE_LANE else measured[name]
        print(f"{name:>30} {state['reduction']:>12.2f} {state['computation']:>12.2f} "
              f"{state['macs']:>9} {state['schedule']:>7.0f}")
    print(f"RD1 multi-lane fixture exercises the tree ok")
    print(f"RD2 reduction energy == issues x (lanes-1) x unit cost ok")
    print(f"RD3 2x the unit cost -> 2x reduction only; MACs/compute/schedule unchanged ok")
    print(f"RD4 single-lane config has no reduction work ok")
    print(f"RD5 the {axis_delta:.2f} axis delta appears once in the MAC subtotal and once in the "
          f"layer total, and nowhere else ok")
    print("ALL REDUCTION ENERGY CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
