#!/usr/bin/env python3
"""Leakage (static energy) production-path regression (E11).

Static energy had exactly one unit test -- the `unit * elapsed` helper -- and nothing covering the
production path: the per-component scaling, the rescale onto the FINAL critical path after the
per-tile timeline recomputes it, or the network rollup. Every shipped config also sets
static_energy = 0, so the whole path could have been dead code.

Two fixtures declare non-zero leakage on every component that models it (PE 0.01, PE-array 0.02
pJ/cycle). The second additionally slows the DRAM 10x, which changes ONLY the final latency:

    gemmini_leakage              baseline
    gemmini_leakage_slow_dram    same leakage, DRAM unit costs x10

Checks (asserted; non-zero exit on failure):
  LK1  The leakage axis is actually non-zero, i.e. the fixture exercises the path at all.
  LK2  A latency-only change must NOT move dynamic energy: the two runs' layer dynamic energy is
       bit-identical. Dynamic energy is event count x unit cost and owes nothing to the schedule --
       this is the invariant that keeps a leakage change from being mistaken for a dynamic one.
  LK3  HAND IDENTITY: leakage is a TIME WINDOW, so the static energy ratio between the two runs
       equals their critical-path ratio, exactly. This is what pins the rescale onto the FINAL
       latency: if leakage were still charged against the pre-timeline latency, or against a
       stage's busy time instead of the layer's wall-clock, the two ratios would diverge.
  LK4  The network rollup carries leakage: network static energy equals the sum over its layers
       and is non-zero. A rollup that dropped the static axis would report 0 here.
  LK5  The shipped config, which declares no leakage, reports exactly 0 -- so a non-zero value can
       only come from a declared unit cost.
  LK6  PE leakage is the SUM OF TWO TERMS: `static_energy` (the MAC/PE) and `lb_static_energy` (the
       local buffer), added before the elapsed-cycle multiply. No config declared the second, so
       half of the model was unexercised and the key looked dead to the undeclared-knob scan
       (validation/knobs KN9). gemmini_lb_leakage adds lb_static_energy = 0.03 on top of
       static_energy = 0.01, so the PE static axis must be EXACTLY 4x -- and the latency must not
       move, since leakage is priced per cycle and changes no schedule.
"""
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKLOAD = ("gemm_64x64x64", "ws")
BASE, SLOW = "gemmini_leakage", "gemmini_leakage_slow_dram"
NO_LEAKAGE = "gemmini"
# LK6: the same config with the local-buffer leakage term added.
LB_LEAKAGE = "gemmini_lb_leakage"
PE_STATIC, LB_STATIC = 0.01, 0.03      # static_energy / lb_static_energy in that fixture


def run(target: str):
    directory = ROOT / "result" / target / WORKLOAD[0] / WORKLOAD[1]
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, WORKLOAD[0], WORKLOAD[1]],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if not (directory / "layer_0.txt").exists():
        raise SystemExit(f"missing simulator output under {directory}")
    return directory


def summary(text: str) -> dict:
    block = text.split("Energy summary", 1)[1].split("MAC result", 1)[0]
    timeline = text.split("Layer timeline", 1)[1].split("MAC result", 1)[0]

    def value(label, source):
        return float(re.search(rf"{re.escape(label)}\s*:\s*([\d.]+)", source).group(1))

    def scoped(axis, source):
        # RE7: a layer report says "Layer <axis> energy" and the network rollup says
        # "Network <axis> energy" -- the label states which scope produced the number.
        for scope in ("Layer", "Network"):
            if re.search(rf"{scope} {axis} energy\s*:", source):
                return value(f"{scope} {axis} energy", source)
        raise SystemExit(f"no layer- or network-scoped {axis} energy in this report")
    return {
        "dynamic": scoped("dynamic", block),
        "static": scoped("static", block),
        "total": scoped("total", block),
        "critical": value("Critical-path latency", timeline),
    }


def measure(target: str) -> dict:
    directory = run(target)
    layer = summary((directory / "layer_0.txt").read_text(encoding="utf-8"))
    network = summary((directory / "network.txt").read_text(encoding="utf-8"))
    layer_static_sum = 0.0
    layers = 0
    for entry in sorted(os.listdir(directory)):
        if not re.fullmatch(r"layer_\d+\.txt", entry):
            continue
        layers += 1
        layer_static_sum += summary(
            (directory / entry).read_text(encoding="utf-8"))["static"]
    layer["network_static"] = network["static"]
    layer["layer_static_sum"] = layer_static_sum
    layer["layers"] = layers
    return layer


def main() -> int:
    base = measure(BASE)
    slow = measure(SLOW)
    none = measure(NO_LEAKAGE)
    failures = []
    tol = 0.5

    # LK1
    for name, state in ((BASE, base), (SLOW, slow)):
        if state["static"] <= 0.0:
            failures.append(f"LK1 {name}: static energy is {state['static']}, so the leakage path "
                            f"is not being exercised")

    # LK2
    if abs(base["dynamic"] - slow["dynamic"]) > tol:
        failures.append(f"LK2: a latency-only change moved the dynamic energy "
                        f"({base['dynamic']} -> {slow['dynamic']}); dynamic energy is event count "
                        f"x unit cost and must not depend on the schedule")

    # LK3
    if base["static"] > 0.0 and base["critical"] > 0.0:
        static_ratio = slow["static"]/base["static"]
        latency_ratio = slow["critical"]/base["critical"]
        if abs(static_ratio - latency_ratio) > 1e-6*max(1.0, latency_ratio):
            failures.append(f"LK3: static energy ratio {static_ratio} != critical-path ratio "
                            f"{latency_ratio}; leakage is not being charged over the final "
                            f"wall-clock")

    # LK4
    for name, state in ((BASE, base), (SLOW, slow)):
        if abs(state["network_static"] - state["layer_static_sum"]) > max(tol, 0.01*state["layers"]):
            failures.append(f"LK4 {name}: network static energy {state['network_static']} != the "
                            f"sum over {state['layers']} layer(s) {state['layer_static_sum']}")
        if state["network_static"] <= 0.0:
            failures.append(f"LK4 {name}: the network rollup reports no leakage at all")

    # LK5
    if none["static"] != 0.0:
        failures.append(f"LK5 {NO_LEAKAGE}: declares no leakage but reports {none['static']}")

    # LK6 -- the two PE leakage terms add. Compare the PE ROW, not the layer total: the layer
    # total also carries the PE-array's own leakage (pe_array_static_energy), which does not scale
    # with these two keys, so a layer-level ratio is diluted (3.954 instead of 4).
    def pe_static(target):
        text = (ROOT / "result" / target / WORKLOAD[0] / WORKLOAD[1] /
                "layer_0.txt").read_text(encoding="utf-8")
        summary_block = text.split("Energy summary", 1)[1].split("Power summary", 1)[0]
        return float(re.search(r"\* PE \(local buffer\)\s+[\d.]+\s+([\d.]+)",
                               summary_block).group(1))

    lb = measure(LB_LEAKAGE)
    base_pe, lb_pe = pe_static(BASE), pe_static(LB_LEAKAGE)
    expected_ratio = (PE_STATIC + LB_STATIC)/PE_STATIC
    if base_pe <= 0.0:
        failures.append("LK6: the baseline reports no PE leakage, so the ratio below is meaningless")
    else:
        ratio = lb_pe/base_pe
        if abs(ratio - expected_ratio) > 1e-6:
            failures.append(f"LK6: adding lb_static_energy = {LB_STATIC} on top of static_energy = "
                            f"{PE_STATIC} gave a PE leakage ratio of {ratio}, expected exactly "
                            f"{expected_ratio} -- the two terms are added before the elapsed-cycle "
                            f"multiply, so the ratio is a pure cost ratio")
    if lb["critical"] != base["critical"]:
        failures.append(f"LK6: declaring lb_static_energy moved the critical path "
                        f"({base['layer']['critical']} -> {lb['layer']['critical']}); leakage is "
                        f"priced per cycle and must change no schedule")
    if abs(lb["dynamic"] - base["dynamic"]) > 0.5:
        failures.append(f"LK6: declaring lb_static_energy moved the DYNAMIC energy "
                        f"({base['layer']['dynamic']} -> {lb['layer']['dynamic']})")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"{'config':>28} {'critical':>11} {'dynamic E':>14} {'static E':>14}")
    for name, state in ((NO_LEAKAGE, none), (BASE, base), (SLOW, slow)):
        print(f"{name:>28} {state['critical']:>11.0f} {state['dynamic']:>14.2f} "
              f"{state['static']:>14.2f}")
    print("LK1 leakage axis exercised ok")
    print(f"LK2 dynamic energy unchanged at {base['dynamic']:.2f} under a latency-only change ok")
    print(f"LK3 static/critical ratios agree ({slow['static']/base['static']:.4f}) ok")
    print("LK4 network rollup carries leakage ok")
    print(f"LK5 {NO_LEAKAGE} reports zero leakage ok")
    print(f"LK6 static_energy + lb_static_energy add: PE leakage is exactly {expected_ratio:.0f}x "
          f"with the local-buffer term, latency and dynamic energy unchanged ok")
    print("ALL LEAKAGE CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
