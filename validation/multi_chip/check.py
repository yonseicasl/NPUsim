#!/usr/bin/env python3
"""Multi-chip live regression (P2-1/P2-2): two small 2-chip synthetic fixtures that need
no dataset download, so CE6 (first_active_pe initializes MIN statistics once across ALL
chips, not per chip), CE4 (entity/datatype aggregation order) and the NoP source-read
sharing contract have executable 2-chip regressions instead of none.

The two fixtures are the same hardware (configs/accelerators/gemmini_2chip.cfg and its
byte-for-byte copy gemmini_2chip_ksplit.cfg) running the same total work under MIRRORED
mappings:
  * gemmini_2chip        -- the 2-chip split is over B, so the WEIGHT tile is BROADCAST
                            (both chips need the same chunk) and input/output are
                            partitioned;
  * gemmini_2chip_ksplit -- the split is over K, so the weight tile is PARTITIONED and
                            the INPUT tile is broadcast.
Diffing them isolates the sharing contract from every other cost in the model (MC6/MC7).

NOTE on asymmetric fixtures. The audit asked for a fixture where the two chips carry
DIFFERENT per-datatype cycle distributions, to separate max_chip(combine_type(...)) from
combine_type(max_chip(...)). That is unreachable from configuration alone:
global_buffer_t::update_tile_size() assigns every chip the SAME
scheduler tile (m_scheduler->tile_size[GLOBAL_BUFFER]), so per-chip GLB access/fill
values are identical by construction and the two reduction orders always agree here.
The ordering itself is therefore covered where it CAN be made asymmetric -- as a
hand-computed unit test of the pure reduction (entity_combined_cycles() in
unittest/validation_test.cc, including the audit's own 100/100 two-chip example).
Making this fixture genuinely asymmetric requires per-chip mapping support, which is a
model feature, not a test fixture.

Checks (asserted; non-zero exit on failure):
  MC1  MIN <= AVG <= MAX computation cycle (CE6: a per-chip-reset MIN/AVG bug
       would corrupt this ordering or leave AVG un-averaged across chips)
  MC2  # of computations == analytic K*B*C*R*S (both chips' halves counted --
       a broken active-chip loop would silently drop one chip's work)
  MC3  Multi-chip busy-cycle axes (access/link) are non-zero (2 active chips
       genuinely exchange data over the NoP fabric)
  MC4  Compute-schedule latency == Computation cycle + Fold fill cycle
  MC5  network.txt busy-cycle axes are populated (P0 contract, re-asserted
       specifically against a multi-chip config)
  MC6  NoP source-read SHARING: the two mirrored mappings must report IDENTICAL
       multi-chip access cycles per datatype. A shared source read costs
       ceil(processor_tile/line) whichever tensor is broadcast, because a broadcast
       tensor's processor tile equals one chip's tile and is read once, while a
       partitioned tensor's is num_chips x one chip's tile and is read once per chip.
       Charging every chip a full read of a broadcast chunk breaks that symmetry: the
       B-split run's weight access would exceed the K-split run's.
  MC7  L7 MULTICAST link sharing, as an A/B on the nop_multicast flag alone
       (gemmini_2chip vs gemmini_2chip_unicast, same mapping, same hardware otherwise).
       A broadcast datatype -- one whose tile does not depend on the split dimension -- must
       cross the shared package ingress ONCE with multicast and once PER CHIP without it, so
       its NoP link transactions must scale by exactly the active chip count between the two
       runs. A partitioned datatype needs a distinct chunk per chip either way, so its count
       must be IDENTICAL in both. On a bus this is not an optimization but a correctness fix:
       one transmission is physically seen by every receiver.
  MC8  Both mappings do the same total work (identical computation counts), which is
       what makes the MC6 diff attributable to the split dimension alone.
  MC9  Mirror consistency: the B-split run's input link count equals the K-split run's weight
       count and vice versa, since the split dimension is what decides which tensor is
       broadcast. NOTE this is not a discriminating test on its own in this fixture -- under
       multicast a broadcast of a 2x-larger tile costs the same as a unicast of two 1x tiles,
       so the counts coincide numerically. MC7 is the check that separates the contracts.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULT = ROOT / "result" / "gemmini_2chip" / "gemm_64x64x64" / "ws"
KSPLIT_RESULT = ROOT / "result" / "gemmini_2chip_ksplit" / "gemm_64x64x64" / "ws"
UNICAST_RESULT = ROOT / "result" / "gemmini_2chip_unicast" / "gemm_64x64x64" / "ws"
TYPES = ("Input data", "Weight", "Output data")
ACTIVE_CHIPS = 2

# K=64, B=64 (32/32 split across 2 chips), P=Q=1, C=64, R=S=1 (single [connected] layer).
EXPECTED_COMPUTATIONS = 64 * 64 * 64 * 1 * 1


def multi_chip_section(text: str) -> str:
    return text.split("Multi chip result", 1)[1]


def per_type(section: str, header: str, pattern: str) -> dict:
    """First occurrence of each datatype row under `header`. Keeping the FIRST match
    matters: the tail of the section also carries the DRAM tables, whose rows use the
    same datatype labels, so a last-wins dict would silently read the wrong block."""
    body = section.split(header, 1)[1]
    found = {}
    for name, value in re.findall(rf"({'|'.join(TYPES)})\s*:\s*{pattern}", body):
        found.setdefault(name, value)
    missing = [t for t in TYPES if t not in found]
    if missing:
        raise AssertionError(f"missing {missing} under {header!r}")
    return {t: float(found[t]) for t in TYPES}


def metric(text: str, label: str, default: float = None) -> float:
    match = re.search(rf"{re.escape(label)}\s*:\s*([\d.]+)", text)
    if match is None:
        if default is not None:
            return default
        raise ValueError(f"missing {label!r} in simulator result")
    return float(match.group(1))


def busy_cycle_axes(text: str, source: str):
    section = text.split("Busy-cycle axes", 1)[1].split("\n\n", 1)[0]
    axes = [
        tuple(float(v) for v in triple)
        for triple in re.findall(r":\s*([\d.]+)\s*/\s*([\d.]+)\s*/\s*([\d.]+)\s*cycles", section)
    ]
    if len(axes) != 5:
        raise AssertionError(f"{source}: expected 5 busy-cycle axis rows, found {len(axes)}")
    return axes


def main() -> int:
    layer_path = RESULT / "layer_0.txt"
    network_path = RESULT / "network.txt"
    if not layer_path.exists() or not network_path.exists():
        print(f"missing fixture output; run: ./npusim.sh run gemmini_2chip gemm_64x64x64 ws", file=sys.stderr)
        return 1

    layer_text = layer_path.read_text(encoding="utf-8")
    network_text = network_path.read_text(encoding="utf-8")

    failures = []

    # MC1
    comp_min = metric(layer_text, "MIN")
    comp_max = metric(layer_text, "MAX")
    comp_avg = metric(layer_text, "AVG")
    if not (comp_min <= comp_avg + 1e-6 and comp_avg <= comp_max + 1e-6):
        failures.append(f"MC1: MIN({comp_min}) <= AVG({comp_avg}) <= MAX({comp_max}) violated")

    # MC2
    num_computations = int(re.search(r"# of computations\s*:\s*(\d+)", layer_text).group(1))
    if num_computations != EXPECTED_COMPUTATIONS:
        failures.append(f"MC2: # of computations {num_computations} != expected {EXPECTED_COMPUTATIONS} "
                        f"(a chip's work may have been dropped)")

    # MC3
    axes = busy_cycle_axes(network_text, str(network_path))
    multi_chip_access, multi_chip_link, _ = axes[1]
    if multi_chip_access <= 0.0 or multi_chip_link <= 0.0:
        failures.append(f"MC3: multi-chip busy axes are degenerate (access={multi_chip_access}, "
                        f"link={multi_chip_link}) despite 2 active chips")

    # MC4
    for label, text, path in (("layer", layer_text, layer_path), ("network", network_text, network_path)):
        compute_schedule = metric(text, "Compute-schedule latency")
        computation = metric(text, "Computation cycle")
        fold_fill = metric(text, "Fold fill cycle", default=0.0)
        identity = computation + fold_fill
        if abs(compute_schedule - identity) > max(0.2, abs(identity) * 1e-6):
            failures.append(f"MC4 ({label}): Compute-schedule latency ({compute_schedule}) != "
                            f"Computation cycle + Fold fill cycle ({identity}) [{path}]")

    # MC5
    if all(value == 0.0 for row in axes for value in row):
        failures.append(f"MC5: all network busy-cycle axes are zero [{network_path}]")

    # MC6/MC7/MC8: diff the mirrored mapping.
    ksplit_path = KSPLIT_RESULT / "layer_0.txt"
    if not ksplit_path.exists():
        print("missing mirrored fixture output; run: "
              "./npusim.sh run gemmini_2chip_ksplit gemm_64x64x64 ws", file=sys.stderr)
        return 1
    ksplit_text = ksplit_path.read_text(encoding="utf-8")
    b_split = multi_chip_section(layer_text)
    k_split = multi_chip_section(ksplit_text)

    access_header = "Access cycle (DRAM fill + NoP source reads)"
    b_access = per_type(b_split, access_header, r"([\d.]+) cycles")
    k_access = per_type(k_split, access_header, r"([\d.]+) cycles")
    b_link = per_type(b_split, "NoP transactions", r"\d+/\d+/(\d+)")
    k_link = per_type(k_split, "NoP transactions", r"\d+/\d+/(\d+)")

    # MC6
    for t in TYPES:
        if abs(b_access[t] - k_access[t]) > 0.5:
            failures.append(
                f"MC6: {t} multi-chip access cycles differ between the mirrored mappings "
                f"(B-split {b_access[t]} vs K-split {k_access[t]}); a broadcast tile's NoP "
                f"source read is not being shared across chips")

    # MC7: multicast A/B on the same mapping. In the B-split fixture WEIGHT is the broadcast
    # datatype (its tile does not depend on B) and INPUT is partitioned.
    unicast_path = UNICAST_RESULT / "layer_0.txt"
    if not unicast_path.exists():
        print("missing unicast fixture output; run: "
              "./npusim.sh run gemmini_2chip_unicast gemm_64x64x64 ws", file=sys.stderr)
        return 1
    unicast = multi_chip_section(unicast_path.read_text(encoding="utf-8"))
    u_link = per_type(unicast, "NoP transactions", r"\d+/\d+/(\d+)")
    expected_broadcast = b_link["Weight"]*ACTIVE_CHIPS
    if abs(u_link["Weight"] - expected_broadcast) > 0.5:
        failures.append(
            f"MC7: without multicast the broadcast weight should cross the ingress once per "
            f"chip ({b_link['Weight']} x {ACTIVE_CHIPS} = {expected_broadcast}), but the "
            f"unicast run reports {u_link['Weight']}")
    if abs(u_link["Input data"] - b_link["Input data"]) > 0.5:
        failures.append(
            f"MC7: the partitioned input needs a distinct chunk per chip either way, so its "
            f"NoP link count must not change with multicast "
            f"(multicast {b_link['Input data']} vs unicast {u_link['Input data']})")
    if abs(u_link["Weight"] - b_link["Weight"]) < 0.5:
        failures.append(
            "MC7: the multicast and unicast runs report the same broadcast weight link count, "
            "so the flag is not being applied and the A/B proves nothing")

    # MC8
    ksplit_computations = int(re.search(r"# of computations\s*:\s*(\d+)", ksplit_text).group(1))
    if ksplit_computations != num_computations:
        failures.append(f"MC8: mirrored mapping does the wrong amount of work "
                        f"({ksplit_computations} != {num_computations}); the MC6 diff is "
                        f"not attributable to the split dimension alone")

    # MC9
    if not (abs(b_link["Input data"] - k_link["Weight"]) < 0.5 and
            abs(b_link["Weight"] - k_link["Input data"]) < 0.5):
        failures.append(
            f"MC9: NoP link transactions do not mirror between the two mappings "
            f"(B-split in/wt {b_link['Input data']}/{b_link['Weight']}, "
            f"K-split in/wt {k_link['Input data']}/{k_link['Weight']}); the split dimension "
            f"decides which tensor is broadcast")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"MC1 MIN({comp_min})<=AVG({comp_avg})<=MAX({comp_max}) ok")
    print(f"MC2 computations={num_computations} ok")
    print(f"MC3 multi-chip axes access={multi_chip_access} link={multi_chip_link} ok")
    print("MC4 compute-schedule identity (layer+network) ok")
    print("MC5 network busy-cycle axes populated ok")
    print(f"MC6 NoP source-read sharing: access identical under mirrored splits "
          f"(in/wt/out {b_access['Input data']}/{b_access['Weight']}/{b_access['Output data']}) ok")
    print(f"MC7 multicast: broadcast weight {b_link['Weight']} -> {u_link['Weight']} without "
          f"multicast (x{ACTIVE_CHIPS}); partitioned input unchanged at {b_link['Input data']} ok")
    print(f"MC8 both mappings compute {ksplit_computations} MACs ok")
    print(f"MC9 link mirror: B-split in/wt {b_link['Input data']}/{b_link['Weight']} vs "
          f"K-split {k_link['Input data']}/{k_link['Weight']} ok")
    print("ALL MULTI-CHIP CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
