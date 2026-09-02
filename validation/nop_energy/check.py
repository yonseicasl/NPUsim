#!/usr/bin/env python3
"""Multi-chip NoP energy regression (E12/E9).

The NoP separates two quantities on purpose: the BUSIEST link sets latency and serialized
transaction counts, while the TOTAL link traversals set energy. That contract had A/B coverage on
the traffic side (MC7) but none on the energy side, because every shipped config sets
nop_energy = 0.

Two fixtures place FOUR chips in a 1x4 MESH row with nop_energy = 5, differing only in
`nop_multicast`. A 2-chip mesh cannot discriminate the contracts -- its multicast tree and its
unicast routes cover the same two links -- so the row is long enough for them to diverge:

    unicast   routes cross 1 + 1 + 2 + 3 = 7 links; the ingress carries 4 copies
    multicast tree crosses the ingress + 3 row links = 4; the ingress carries 1 copy

The batch is split across the chips, so WEIGHT is the broadcast datatype and INPUT is partitioned.

Checks (asserted; non-zero exit on failure):
  NP1  HAND IDENTITY, from the report alone: the reported serialized transactions count the
       BOTTLENECK link's copies, and the energy counts ALL traversals, so
         energy == reported transactions * nop_energy * total_traversals / bottleneck_copies
       -- 4/1 for a multicast broadcast, 7/4 for anything unicast. A model that used one quantity
       for both latency and energy cannot satisfy this for both datatypes at once.
  NP2  The BROADCAST datatype's energy ratio between the two modes is exactly the traversal ratio
       7/4, while its transaction ratio is the bottleneck ratio 4/1. Two different ratios from one
       config change is the evidence that the two quantities are genuinely separate.
  NP3  The PARTITIONED datatype is unaffected: it needs a distinct chunk per chip either way, so
       both its transactions and its energy are identical in the two runs.
  NP4  Multicast is cheaper in energy AND in serialized traffic than unicast for the broadcast
       tensor -- a sanity direction check, so an inverted sign cannot pass NP1 by coincidence.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKLOAD = ("gemm_64x64x64", "ws")
MULTICAST, UNICAST = "gemmini_nop4_mc", "gemmini_nop4_uc"
NOP_ENERGY = 5.0
ACTIVE_CHIPS = 4
# 1x4 mesh from the ingress: hops 1,1,2,3
UNICAST_TRAVERSALS, UNICAST_BOTTLENECK = 7.0, 4.0
MULTICAST_TRAVERSALS, MULTICAST_BOTTLENECK = 4.0, 1.0
BROADCAST, PARTITIONED = "Weight", "Input data"


def run(target: str) -> str:
    layer = ROOT / "result" / target / WORKLOAD[0] / WORKLOAD[1] / "layer_0.txt"
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, WORKLOAD[0], WORKLOAD[1]],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if not layer.exists():
        raise SystemExit(f"missing simulator output: {layer}")
    return layer.read_text(encoding="utf-8")


def measure(text: str) -> dict:
    section = text.split("Multi chip result", 1)[1].split("DRAM result", 1)[0]
    transactions = {}
    for name, value in re.findall(
            r"(Input data|Weight|Output data)\s*:\s*\d+/\d+/(\d+)",
            section.split("NoP transactions", 1)[1]):
        transactions.setdefault(name, int(value))
    energy = {}
    for name, value in re.findall(
            r"(Input data|Weight|Output data)\s*:\s*([\d.]+) (?:pJ|normalized|uncalibrated)",
            section.split("Interconnection energy", 1)[1]):
        energy.setdefault(name, float(value))
    return {"transactions": transactions, "energy": energy}


def main() -> int:
    multicast = measure(run(MULTICAST))
    unicast = measure(run(UNICAST))
    failures = []
    tol = 0.5

    # NP1
    cases = [
        (MULTICAST, multicast, BROADCAST, MULTICAST_TRAVERSALS, MULTICAST_BOTTLENECK),
        (MULTICAST, multicast, PARTITIONED, UNICAST_TRAVERSALS, UNICAST_BOTTLENECK),
        (UNICAST, unicast, BROADCAST, UNICAST_TRAVERSALS, UNICAST_BOTTLENECK),
        (UNICAST, unicast, PARTITIONED, UNICAST_TRAVERSALS, UNICAST_BOTTLENECK),
    ]
    for name, state, datatype, traversals, bottleneck in cases:
        expected = state["transactions"][datatype]*NOP_ENERGY*traversals/bottleneck
        if abs(state["energy"][datatype] - expected) > tol:
            failures.append(
                f"NP1 {name} {datatype}: NoP energy {state['energy'][datatype]} != hand "
                f"{expected} ({state['transactions'][datatype]} reported txns x {NOP_ENERGY} x "
                f"{traversals:.0f} traversals / {bottleneck:.0f} bottleneck copies)")

    # NP2
    energy_ratio = unicast["energy"][BROADCAST]/multicast["energy"][BROADCAST]
    expected_energy_ratio = UNICAST_TRAVERSALS/MULTICAST_TRAVERSALS
    if abs(energy_ratio - expected_energy_ratio) > 1e-9:
        failures.append(f"NP2: broadcast energy ratio {energy_ratio} != the traversal ratio "
                        f"{expected_energy_ratio} (7/4)")
    txn_ratio = unicast["transactions"][BROADCAST]/multicast["transactions"][BROADCAST]
    expected_txn_ratio = UNICAST_BOTTLENECK/MULTICAST_BOTTLENECK
    if abs(txn_ratio - expected_txn_ratio) > 1e-9:
        failures.append(f"NP2: broadcast transaction ratio {txn_ratio} != the bottleneck ratio "
                        f"{expected_txn_ratio} (4/1)")
    if abs(energy_ratio - txn_ratio) < 1e-9:
        failures.append(f"NP2: the energy and transaction ratios coincide ({energy_ratio}); the "
                        f"fixture then cannot show that latency and energy use different "
                        f"quantities")

    # NP3
    for field in ("transactions", "energy"):
        if abs(multicast[field][PARTITIONED] - unicast[field][PARTITIONED]) > tol:
            failures.append(f"NP3 {field} {PARTITIONED}: changed with nop_multicast "
                            f"({multicast[field][PARTITIONED]} vs {unicast[field][PARTITIONED]}); "
                            f"a partitioned datatype needs its own chunk per chip either way")

    # NP4
    if not multicast["energy"][BROADCAST] < unicast["energy"][BROADCAST]:
        failures.append("NP4: multicast is not cheaper in energy than unicast for the broadcast "
                        "tensor")
    if not multicast["transactions"][BROADCAST] < unicast["transactions"][BROADCAST]:
        failures.append("NP4: multicast does not reduce the serialized traffic for the broadcast "
                        "tensor")

    # ---- NP5-NP7: the 2D multicast tree (E20-6) ----------------------------------------------
    import shutil
    sys.path.insert(0, str(ROOT / "validation"))
    import perturb_lib as pl

    # The chip grid of each fixture, and the spanning-tree edge count worked out by hand from the
    # model's own construction: the ingress link into (0,0), the row-0 horizontals up to the
    # rightmost used column, and the per-column verticals.
    #   1x4 : 1 + 3 + (0+0+0+0)      = 4   <- a single row, so the COLUMN term is 0
    #   2x2 : 1 + 1 + (1+1)          = 4
    #   3-wide, 4 active: 1 + 2 + (1+0+0) = 4   <- a non-rectangular active set
    # All three are 4, which is not a coincidence: a spanning tree over N routers plus the ingress
    # is N links whatever the grid's shape. What differs is HOW the 4 is composed, and only a grid
    # with more than one row exercises the column term at all -- which is why a max->min mutation of
    # `column_height[x]` was invisible to every fixture before these two.
    GRIDS = {"gemmini_nop4_mc": (4, 3, 4), "gemmini_nop2x2_mc": (4, 1, 4),
             "gemmini_nop2x3_mc": (4, 2, 4)}      # (active chips, row-only links, tree links)

    def weight_nop_energy(target):
        text = pl.run(target, *WORKLOAD)
        if text is None:
            return None, None
        block = text.split("Multi chip result", 1)[1]
        energy = [float(v) for _, v in re.findall(
            r"(Input data|Weight|Output data)\s*:\s*([\d.]+) \S+",
            block.split("Interconnection energy", 1)[1])[:3]][1]
        txns = re.findall(r"\*\s*(?:Input data|Weight|Output data)\s*:\d+/\d+/(\d+)",
                          text.split("GLB <-> multi-chip NoP transactions", 1)[1])[:3]
        return energy, int(txns[1])

    tree_counts = {}
    for target, (active, row_only, tree) in GRIDS.items():
        energy, txns = weight_nop_energy(target)
        if energy is None:
            failures.append(f"NP5 {target}: does not run")
            continue
        measured = energy/(txns*NOP_ENERGY)
        tree_counts[target] = measured
        # NP6 -- hand identity against the enumerated edge set
        if abs(measured - tree) > 1e-9:
            failures.append(f"NP6 {target}: {measured} tree links implied by the reported energy "
                            f"({energy} = {txns} txns x {NOP_ENERGY} x links), hand count {tree} "
                            f"(ingress + row-0 horizontals + per-column verticals)")
        # NP7 -- the COLUMN term must be exercised on a multi-row grid
        if row_only < tree and measured <= row_only + 1e-9:
            failures.append(f"NP7 {target}: {measured} links equals the row-only count {row_only}, "
                            f"so the per-column verticals contribute nothing. On a grid with more "
                            f"than one row they must -- this is exactly what a max->min mutation of "
                            f"column_height[x] produces, and what a single-row fixture cannot see")
    # NP5 -- the total is shape-invariant for a fixed active-chip count
    if len(tree_counts) == len(GRIDS) and len(set(round(v, 9) for v in tree_counts.values())) != 1:
        failures.append(f"NP5: the tree size differs across grid shapes with the same 4 active chips "
                        f"({ {k: v for k, v in tree_counts.items()} }); a spanning tree over N "
                        f"routers plus the ingress is N links whatever the shape")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"{'mode':>10} {'weight txns':>12} {'weight E':>10} {'input txns':>11} {'input E':>9}")
    for name, state in ((MULTICAST, multicast), (UNICAST, unicast)):
        label = "multicast" if name == MULTICAST else "unicast"
        print(f"{label:>10} {state['transactions'][BROADCAST]:>12} "
              f"{state['energy'][BROADCAST]:>10.0f} {state['transactions'][PARTITIONED]:>11} "
              f"{state['energy'][PARTITIONED]:>9.0f}")
    print("NP1 energy == reported txns x nop_energy x traversals / bottleneck copies ok")
    print(f"NP2 broadcast energy ratio {energy_ratio:.2f} (7/4) != transaction ratio "
          f"{txn_ratio:.2f} (4/1): latency and energy use different quantities ok")
    print("NP3 the partitioned datatype is unaffected by nop_multicast ok")
    print("NP4 multicast is cheaper in both energy and serialized traffic ok")
    print(f"NP5 the tree is {list(tree_counts.values())[0]:.0f} links in all three grid shapes with "
          f"4 active chips (shape-invariant) ok")
    print("NP6 each shape's tree matches its hand-enumerated edge set ok")
    print("NP7 the per-column verticals contribute on every multi-row grid ok")
    print("ALL NOP ENERGY CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
