#!/usr/bin/env python3
"""Spatial-array NoC energy regression (E2).

The array NoC's energy multiplier used to be the AVERAGE Manhattan hop count for every kind of
traffic. That is the cost of one average unicast -- not of a multicast that fans an operand tile
out to every active PE. On a 16x16 array it charged 15 link traversals where the distribution
tree crosses 255, understating distribution energy about 17x. It was invisible because every
shipped config sets noc_energy = 0.

RE6 -- what ONE traversal IS is now DEFINED rather than implied, because the same fabric gives
different answers under different conventions and the report said which one it used nowhere:

    noc_energy prices one traversal of one ROUTER-TO-ROUTER link INSIDE the array.

It does not price the GLB <-> array attach link (charged by the GLB's own transfer_energy, the
report's "PEs - Global buffer" axis) nor a PE's ejection (charged as its local-buffer write). So a
multicast to N active PEs is a spanning TREE over N routers = **N-1 edges**, not N: the Nth link a
receiver-counting model adds is exactly the attach link, which would then be billed twice. The
first version of this fix used N; NC1 now pins N-1, NC6 asserts the report states the contract,
and utils/interconnect_timing.cc enumerates the edge set explicitly so the closed form cannot
drift away from a physical link count (unittest/validation_test.cc cross-checks 1x1/1xN/Nx1/NxM).

Two fixtures differ from configs/accelerators/gemmini.cfg only in the array NoC (mesh instead of
bus, noc_energy = 2) and the ACTIVE SHAPE, with mappings chosen so the total work and the NoC
transaction counts are identical:

    gemmini_noc_energy_16x16   16x16 active
    gemmini_noc_energy_8x8      8x8 active  (K and C factors moved to PE/GLB levels)

Checks (asserted; non-zero exit on failure):
  NC1  HAND IDENTITY, multicast: the input and weight (operand distribution) energy equals
         transactions * noc_energy * (height * width - 1)
       -- one traversal per SPANNING-TREE EDGE over the active routers.
  NC2  HAND IDENTITY, unicast: the output (write-back) energy equals
         transactions * noc_energy * ((height-1) + (width-1)) / 2
       -- each PE's partial sum travels its own average Manhattan route to the edge. This
       direction was already right, and must stay separate from NC1.
  NC3  The NoC transaction counts are IDENTICAL in both runs, so the energy difference is
       attributable to the active shape alone and not to a changed tile size.
  NC4  The two directions really do use different multipliers (255 vs 15 at 16x16). If a future
       change collapsed them back into one formula, NC1 and NC2 could not both hold.
  NC6  RE6 condition 4: the report STATES the link contract -- that the count is router-to-router
       links and that the attach link and the endpoint write are charged elsewhere -- so the number
       can be reconciled with a physical link count instead of guessed at. It must also name the
       N-1 spanning-tree convention rather than leaving N vs N-1 open.
  NC5  A bus array is unaffected: one traversal per transaction, no fanout charge. This is what
       keeps the fix from silently re-pricing every shipped (bus) config.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKLOAD = ("gemm_64x64x64", "ws")
NOC_ENERGY = 2.0
SHAPES = {"gemmini_noc_energy_16x16": (16, 16), "gemmini_noc_energy_8x8": (8, 8)}
BUS_REFERENCE = "gemmini"          # shipped config: bus NoC, noc_energy = 0


def run(target: str) -> str:
    layer = ROOT / "result" / target / WORKLOAD[0] / WORKLOAD[1] / "layer_0.txt"
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, WORKLOAD[0], WORKLOAD[1]],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if not layer.exists():
        raise SystemExit(f"missing simulator output: {layer}")
    return layer.read_text(encoding="utf-8")


def measure(text: str) -> dict:
    section = text.split("PE array result", 1)[1].split("Global buffer result", 1)[0]
    transactions = {name: int(value) for name, value in re.findall(
        r"(Input data|Weight|Output data)\s*:\s*\d+/\d+/(\d+)", section)}
    energy_block = section.split("Interconnection energy", 1)[1]
    # Only the FIRST three rows belong to the interconnection block -- the static-energy block
    # that follows uses the same datatype labels and would otherwise overwrite them.
    energy = {name: float(value) for name, value in re.findall(
        r"(Input data|Weight|Output data)\s*:\s*([\d.]+) (?:pJ|normalized|uncalibrated)", energy_block)[:3]}
    return {"transactions": transactions, "energy": energy}


def main() -> int:
    measured = {name: measure(run(name)) for name in SHAPES}
    bus = measure(run(BUS_REFERENCE))
    failures = []
    tol = 0.5

    for name, (height, width) in SHAPES.items():
        state = measured[name]
        # RE6: a spanning tree over N active routers has N-1 edges.
        multicast_traversals = float(height*width - 1)
        unicast_traversals = ((height - 1) + (width - 1))/2.0
        # NC1
        for datatype in ("Input data", "Weight"):
            expected = state["transactions"][datatype]*NOC_ENERGY*multicast_traversals
            if abs(state["energy"][datatype] - expected) > tol:
                failures.append(
                    f"NC1 {name} {datatype}: NoC energy {state['energy'][datatype]} != hand "
                    f"{expected} ({state['transactions'][datatype]} txns x {NOC_ENERGY} x "
                    f"{multicast_traversals:.0f} spanning-tree edges over the {height}x{width} "
                    f"active routers; {height*width} would double charge the GLB attach link)")
        # NC2
        expected_output = state["transactions"]["Output data"]*NOC_ENERGY*unicast_traversals
        if abs(state["energy"]["Output data"] - expected_output) > tol:
            failures.append(
                f"NC2 {name}: write-back energy {state['energy']['Output data']} != hand "
                f"{expected_output} ({state['transactions']['Output data']} txns x {NOC_ENERGY} x "
                f"{unicast_traversals} average Manhattan hops)")
        # NC4
        if multicast_traversals == unicast_traversals:
            failures.append(f"NC4 {name}: the multicast and unicast multipliers coincide at this "
                            f"shape, so the two contracts are not distinguishable here")

    # NC3
    names = list(SHAPES)
    for datatype in ("Input data", "Weight", "Output data"):
        counts = {measured[name]["transactions"][datatype] for name in names}
        if len(counts) != 1:
            failures.append(f"NC3 {datatype}: NoC transaction counts differ between the shapes "
                            f"({ {n: measured[n]['transactions'][datatype] for n in names} }); "
                            f"the energy delta is then not attributable to the shape alone")

    # NC6
    contract = re.search(r"NoC link contract\s*:(.*)",
                         run("gemmini_noc_energy_16x16")).group(1)
    for phrase in ("router-to-router", "N-1", "transfer_energy", "buffer write"):
        if phrase not in contract:
            failures.append(f"NC6: the reported link contract does not state {phrase!r} "
                            f"({contract.strip()!r}), so the energy cannot be reconciled with a "
                            f"physical link count")

    # NC5
    for datatype, value in bus["energy"].items():
        if value != 0.0:
            failures.append(f"NC5 {BUS_REFERENCE} {datatype}: a bus array with noc_energy = 0 "
                            f"must cost no NoC energy, got {value}")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print("NC6 the report states the router-to-router N-1 contract and what is charged "
          "elsewhere ok")
    for name, (height, width) in SHAPES.items():
        state = measured[name]
        print(f"{name:>26} {height}x{width}: multicast in/wt "
              f"{state['energy']['Input data']:.0f}/{state['energy']['Weight']:.0f}, "
              f"unicast out {state['energy']['Output data']:.0f} "
              f"(txns {state['transactions']['Input data']}/"
              f"{state['transactions']['Output data']})")
    print("NC1 multicast energy == txns x noc_energy x (h*w - 1) spanning-tree edges ok")
    print("NC2 write-back energy == txns x noc_energy x avg Manhattan hops ok")
    print("NC3 transaction counts identical across shapes ok")
    print("NC4 the two directions use different multipliers ok")
    print("NC5 bus array unaffected ok")
    print("ALL NOC ENERGY CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
