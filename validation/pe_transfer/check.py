#!/usr/bin/env python3
"""MAC <-> local-buffer transfer axis: the one axis mutation testing found unguarded.

WHY THIS EXISTS. The 2026-08-25 adjudication of the mutation campaign's survivors re-ran all nine
against the full suite (unit tests + 23 gates + check_timing). Eight survived. Seven of those are
explained -- an equivalent mutant, dead code, or code no shipped configuration reaches -- but one
was a real hole:

    components/pe.cc   transfer_cycle[type] += timing.link_transactions*u_transfer_cycle;
                                           -= ...

flipping that sign drives every PE's MAC<->local-buffer transfer cycle negative. stats_t aggregates
the field with std::max() from 0, so the report shows a clean 0.0 instead of the true positive
value, on eyeriss AND on gemmini. Nothing noticed:

  * the dead-knob sweep (KN4) passes -- the knob still MOVES other numbers, just in the wrong
    direction, and liveness has no sign;
  * latency monotonicity (KN7) only guards the critical path, which this axis does not set;
  * T7 combines the printed axes into stage busy with max(), and the PE stage is dominated by
    another axis, so a link axis collapsing to 0 changes no busy number;
  * no gate reads the "Cycle (MAC-Local buffer)" block at all.

So the axis was printed and never checked. These checks close that.

Checks (asserted; non-zero exit on failure):
  PT1  ACTIVITY CONSISTENCY. A datatype with serialized MAC<->LB transactions and a positive
       transfer_cycle_pe must report a POSITIVE transfer cycle. This is what the sign mutant
       violates: work happened, a cost was declared, and the reported cost is zero.
  PT2  HOMOGENEITY. Doubling transfer_cycle_pe exactly doubles every reported transfer cycle and
       moves no transaction count. The axis is linear and homogeneous in its one unit cost, so a
       hardcoded term or a second cost feeding it shows up here.
  PT3  AXIS IDENTITY. The PE stage's LINK busy axis equals the sum of the three transfer cycles,
       exactly. This is what ties the block to the timeline that T7 checks.
  PT4  ENERGY AXIS COVERAGE. transfer_energy_pe is declared 0 in ALL 68 shipped configs, so its
       multiplier is never exercised -- the same shape as the noc_energy defect (wrong by ~17x and
       invisible because every config priced it at zero). Pricing it here makes each transfer
       energy positive, exactly linear in the cost, and moves no cycle.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import perturb_lib as pl                                            # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
BASE, WORKLOAD = "gemmini_all_knobs", ("gemm_64x64x64", "ws")
STREAMS = ("Input data", "Weight", "Output data")
TRANSFER_CYCLE_PE = 1.0     # gemmini_all_knobs [spatial_arch]
EYERISS_RESULT = ROOT/"result"/"eyeriss"/"alexnet"/"silicon"/"layer_4.txt"


def set_key(key, value):
    return lambda text: re.sub(rf"^{re.escape(key)} = .*$", f"{key} = {value}", text,
                               count=1, flags=re.M)


def parse(text):
    """Transactions, transfer cycles/energies, and the PE stage's link busy axis."""
    block = text.split("MAC - Local buffer", 1)[1].split("PE array result", 1)[0]
    txn = {n: int(s) for n, s in
           re.findall(r"(Input data|Weight|Output data)\s*:\s*\d+/\d+/(\d+)", block)}
    cycles = {n: float(v) for n, v in re.findall(
        r"(Input data|Weight|Output data)\s*:\s*([\d.-]+) cycles",
        block.split("Transfer cycle", 1)[1].split("Total cycle", 1)[0])}
    energy = {n: float(v) for n, v in re.findall(
        r"(Input data|Weight|Output data)\s*:\s*([\d.-]+) \S+",
        block.split("Transfer energy", 1)[1].split("Static energy", 1)[0])}
    axis = re.search(r"\* PE \(compute\+LB\)\s*:\s*([\d.]+) /\s*([\d.]+) /\s*([\d.]+)",
                     text.split("Busy-cycle axes", 1)[1])
    return txn, cycles, energy, float(axis.group(2))


def measure(label, edits):
    name = f"__pt_{label}"
    pl.variant(BASE, name, edits)
    text = pl.run(name, *WORKLOAD)
    pl.discard(name)
    if text is None:
        raise SystemExit(f"the {label} variant does not run")
    return parse(text)


def main() -> int:
    failures = []
    base_txn, base_cycles, base_energy, base_axis = measure("base", lambda t: t)

    # PT1
    for s in STREAMS:
        if base_txn[s] > 0 and TRANSFER_CYCLE_PE > 0 and not base_cycles[s] > 0:
            failures.append(f"PT1 {s}: {base_txn[s]} MAC<->LB transactions at "
                            f"transfer_cycle_pe={TRANSFER_CYCLE_PE} report {base_cycles[s]} "
                            f"transfer cycles; a priced, active axis cannot cost nothing")

    # PT2
    d_txn, d_cycles, _, _ = measure("cyc2", set_key("transfer_cycle_pe", TRANSFER_CYCLE_PE*2))
    for s in STREAMS:
        if abs(d_cycles[s] - 2*base_cycles[s]) > 0.5:
            failures.append(f"PT2 {s}: doubling transfer_cycle_pe gave {d_cycles[s]}, "
                            f"expected exactly 2 x {base_cycles[s]}")
        if d_txn[s] != base_txn[s]:
            failures.append(f"PT2 {s}: a CYCLE cost changed the transaction count "
                            f"({base_txn[s]} -> {d_txn[s]})")

    # PT3 -- on this fixture and on the Eyeriss silicon run, which has a different shape
    cases = [("gemmini", base_cycles, base_axis)]
    if EYERISS_RESULT.exists():
        _, e_cycles, _, e_axis = parse(EYERISS_RESULT.read_text(encoding="utf-8"))
        cases.append(("eyeriss", e_cycles, e_axis))
    for name, cyc, axis in cases:
        total = sum(cyc.values())
        if abs(total - axis) > 0.5:
            failures.append(f"PT3 {name}: the PE link busy axis is {axis} but the three transfer "
                            f"cycles sum to {total}; the block and the timeline disagree")

    # PT4
    priced_txn, priced_cycles, priced_energy, _ = measure(
        "en1", set_key("transfer_energy_pe", "1"))
    twice_txn, twice_cycles, twice_energy, _ = measure(
        "en2", set_key("transfer_energy_pe", "2"))
    for s in STREAMS:
        if base_txn[s] > 0 and not priced_energy[s] > 0:
            failures.append(f"PT4 {s}: pricing transfer_energy_pe still reports "
                            f"{priced_energy[s]} transfer energy over {base_txn[s]} transactions")
        if abs(twice_energy[s] - 2*priced_energy[s]) > 0.5:
            failures.append(f"PT4 {s}: doubling transfer_energy_pe gave {twice_energy[s]}, "
                            f"expected exactly 2 x {priced_energy[s]}")
        if abs(priced_cycles[s] - base_cycles[s]) > 0.5:
            failures.append(f"PT4 {s}: an ENERGY cost moved the transfer CYCLE "
                            f"({base_cycles[s]} -> {priced_cycles[s]})")

    print(f"{'stream':>12} {'txns':>12} {'cycles':>12} {'energy@1':>10}")
    for s in STREAMS:
        print(f"{s:>12} {base_txn[s]:>12} {base_cycles[s]:>12.1f} {priced_energy[s]:>10.1f}")
    if failures:
        print("\nFAILED:")
        for f in failures:
            print(" - " + f)
        return 1
    print(f"PT1 every active, priced MAC<->LB stream reports a positive transfer cycle ok")
    print(f"PT2 transfer cycles are exactly linear in transfer_cycle_pe, counts unmoved ok")
    print(f"PT3 the PE link busy axis equals the sum of the transfer cycles on "
          f"{len(cases)} config(s) ok")
    print(f"PT4 transfer_energy_pe (0 in all 68 shipped configs) prices the axis exactly and "
          f"moves no cycle ok")
    print("ALL PE-TRANSFER CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
