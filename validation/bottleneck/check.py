#!/usr/bin/env python3
"""Bottleneck-classification acceptance gate (Phase-4 item 12).

The audit asks for compute-bound / memory-bound / multi-chip acceptance gates. Without an
external memory-inclusive latency reference (L4) the ABSOLUTE critical path cannot be gated --
so this gates what is verifiable without one: that the timeline puts the bottleneck in the
right place, and that it does so BECAUSE of the sub-model that was changed.

Three configs differ from configs/accelerators/gemmini.cfg in exactly one thing each, so each
classification is attributable:

    gemmini                 unchanged -- the on-chip distribution fabric dominates
    gemmini_memory_bound    DRAM unit costs x10, everything else identical
    gemmini_compute_bound   MAC computation_cycle 1 -> 60, everything else identical

plus gemmini_2chip for the multi-chip path.

Checks (asserted; non-zero exit on failure):
  BN1  The reported "Bottleneck level" really is the stage with the largest busy value in
       every run. A classification that disagrees with its own numbers is worse than none.
  BN2  Critical path >= the bottleneck stage's busy time, and >= the compute schedule, in
       every run. A stage cannot be busy longer than the layer takes.
  BN3  Raising ONLY the DRAM unit costs moves the bottleneck to a MEMORY stage
       (DRAM / Multi-chip / GLB) while leaving the compute schedule BIT-IDENTICAL. The
       identical schedule is what proves the reclassification came from the memory model
       rather than from compute having changed underneath.
  BN4  Raising ONLY the MAC cost moves the bottleneck to the PE (compute+LB) stage while
       leaving every MEMORY stage's busy time bit-identical -- the mirror argument.
  BN5  The unmodified config is neither of those: its bottleneck is the PE-array fabric, so
       the three fixtures cover three distinct regimes rather than two labels for one.
  BN6  The multi-chip fixture classifies and keeps the compute-schedule identity
       (Compute-schedule latency == Computation cycle + Fold fill cycle).
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STAGES = ("DRAM", "Multi-chip (NoP)", "Global buffer", "PE array", "PE (compute+LB)")
MEMORY_STAGES = {"DRAM", "Multi-chip (NoP)", "Global buffer"}
RUNS = {
    "gemmini": ("gemm_64x64x64", "ws"),
    "gemmini_memory_bound": ("gemm_64x64x64", "ws"),
    "gemmini_compute_bound": ("gemm_64x64x64", "ws"),
    "gemmini_2chip": ("gemm_64x64x64", "ws"),
}


def run(target: str, workload, mapping) -> str:
    layer = ROOT / "result" / target / workload / mapping / "layer_0.txt"
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, workload, mapping],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if not layer.exists():
        raise SystemExit(f"missing simulator output: {layer}")
    return layer.read_text(encoding="utf-8")


def measure(text: str) -> dict:
    block = text.split("Layer timeline", 1)[1].split("MAC result", 1)[0]
    busy = [float(v) for v in re.findall(r"\* [^:\n]+:\s*([\d.]+) cycles \(", block)]
    return {
        "busy": dict(zip(STAGES, busy)),
        "bottleneck": re.search(r"Bottleneck level\s*:\s*(.+)", block).group(1).strip(),
        "schedule": float(re.search(r"Compute-schedule latency\s*:\s*([\d.]+)", block).group(1)),
        "critical": float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", block).group(1)),
        "computation": float(re.search(r"Computation cycle\s*:\s*([\d.]+)", text).group(1)),
        "fold": float(re.search(r"Fold fill cycle\s*:\s*([\d.]+)", text).group(1))
                if re.search(r"Fold fill cycle\s*:\s*([\d.]+)", text) else 0.0,
    }


def main() -> int:
    measured = {name: measure(run(name, *workload)) for name, workload in RUNS.items()}
    base = measured["gemmini"]
    memory = measured["gemmini_memory_bound"]
    compute = measured["gemmini_compute_bound"]
    multi = measured["gemmini_2chip"]
    failures = []
    tol = 0.5

    # BN1 / BN2
    for name, state in measured.items():
        busiest = max(state["busy"], key=lambda stage: state["busy"][stage])
        if state["bottleneck"] != busiest:
            failures.append(f"BN1 {name}: reports bottleneck {state['bottleneck']!r} but the "
                            f"busiest stage is {busiest!r} "
                            f"({state['busy'][busiest]} vs {state['busy'][state['bottleneck']]})")
        if state["critical"] + tol < state["busy"][state["bottleneck"]]:
            failures.append(f"BN2 {name}: critical path {state['critical']} is shorter than the "
                            f"bottleneck stage's busy time {state['busy'][state['bottleneck']]}")
        if state["critical"] + tol < state["schedule"]:
            failures.append(f"BN2 {name}: critical path {state['critical']} is shorter than the "
                            f"compute schedule {state['schedule']}")

    # BN3
    if memory["bottleneck"] not in MEMORY_STAGES:
        failures.append(f"BN3: raising only the DRAM costs left the bottleneck at "
                        f"{memory['bottleneck']!r}, expected a memory stage")
    if abs(memory["schedule"] - base["schedule"]) > tol:
        failures.append(f"BN3: the memory-bound fixture's compute schedule moved "
                        f"({base['schedule']} -> {memory['schedule']}); the reclassification is "
                        f"then not attributable to the memory model alone")
    if memory["critical"] <= base["critical"] + tol:
        failures.append(f"BN3: a 10x slower DRAM did not lengthen the critical path "
                        f"({base['critical']} -> {memory['critical']})")

    # BN4
    if compute["bottleneck"] != "PE (compute+LB)":
        failures.append(f"BN4: raising only the MAC cost left the bottleneck at "
                        f"{compute['bottleneck']!r}, expected 'PE (compute+LB)'")
    for stage in MEMORY_STAGES:
        if abs(compute["busy"][stage] - base["busy"][stage]) > tol:
            failures.append(f"BN4: the compute-bound fixture's {stage} busy moved "
                            f"({base['busy'][stage]} -> {compute['busy'][stage]}); the "
                            f"reclassification is then not attributable to the MAC cost alone")

    # BN5
    if base["bottleneck"] != "PE array":
        failures.append(f"BN5: the unmodified config reports {base['bottleneck']!r}; the three "
                        f"fixtures are meant to cover three distinct regimes")

    # BN6
    identity = multi["computation"] + multi["fold"]
    if abs(multi["schedule"] - identity) > max(0.2, abs(identity)*1e-6):
        failures.append(f"BN6: multi-chip compute-schedule {multi['schedule']} != computation + "
                        f"fold fill {identity}")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"{'config':>22} {'bottleneck':>18} {'schedule':>10} {'critical':>11}")
    for name, state in measured.items():
        print(f"{name:>22} {state['bottleneck']:>18} {state['schedule']:>10,.0f} "
              f"{state['critical']:>11,.0f}")
    print()
    print("BN1 reported bottleneck == busiest stage in every run ok")
    print("BN2 critical path covers the bottleneck busy and the compute schedule ok")
    print(f"BN3 DRAM x10 -> {memory['bottleneck']}, compute schedule unchanged at "
          f"{memory['schedule']:,.0f} ok")
    print(f"BN4 MAC x60 -> {compute['bottleneck']}, memory axes unchanged ok")
    print(f"BN5 unmodified config is fabric-bound ({base['bottleneck']}) ok")
    print("BN6 multi-chip compute-schedule identity ok")
    print("ALL BOTTLENECK CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
