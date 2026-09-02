#!/usr/bin/env python3
"""Shipped-architecture smoke and cross-architecture invariant gate.

Six of the accelerator configs NPUsim ships -- tpu, tpuv3, simba, maeri, eyerissv2, fsd -- were
referenced by ZERO validation files. Nothing checked that they ran, let alone that their numbers
were self-consistent. That is not the same class of hole as an unexercised knob: entire component
implementations had no coverage, including maeri's [adder_tree] and simba's 36-chip NoP -- the very
NoP whose link contract RE6 changed.

The cause was a data dependency, not neglect: every one of them ships with mappings for real
networks (resnet50, bert, ...) whose weight blobs (weights/*.wgh) are not in the repo, so
`npusim.sh run tpu resnet50 cycle` fails at load. They now each carry a SYNTHETIC gemm_64x64x64
mapping instead, sized to that architecture's own array and chip grid, which needs no weight file.

THOSE MAPPINGS ARE FOR COVERAGE, NOT PERFORMANCE. They put spatial factors on the array and push
everything else to DRAM, which is the only shape that fits a 1-deep systolic PE buffer -- so e.g.
tpu's critical-path latency is ~260,000x its compute schedule, because every batch refetches
everything. Read them as "does this architecture compute the right thing self-consistently", never
as a mapping quality result.

Getting them to run surfaced two real defects, both fixed:
  * The multi-chip capacity check ignored the GLB bypass. global_buffer_t excuses a bypassed
    datatype from its own check, and multi_chip_t already consults chips[i]->bypass[type] when it
    accounts delivery traffic -- but its capacity check did not, so eyerissv2.cfg (weights streamed
    from the ingress straight to the spads, `weight_size = 0` at both levels) could not run at all.
  * multi_chip_t stored its buffer sizes as `unsigned` while global_buffer_t stored the same keys as
    `double`. eyerissv2 declares `input_size = 4.5` (KB), which failed the unsigned parse and
    SILENTLY left the buffer at 0 -- so every tile "exceeded" it. Nothing noticed because nothing
    ran the config.

Checks (asserted; non-zero exit on failure):
  SM1  Every shipped architecture runs and produces a parseable report.
  SM2  CROSS-ARCHITECTURE WORK IDENTITY: the MAC count is 262,144 -- K*B*C for GEMM 64x64x64 -- in
       EVERY architecture. The mapping decides how the work is spread, never how much of it there
       is, so a mapping that fails to cover the workload (or double-covers it) shows up here and
       nowhere else. This is the check that makes the rest meaningful.
  SM3  ENERGY SELF-CONSISTENCY: the layer total equals dynamic + static and equals the sum of the
       six component subtotals, in every architecture.
  SM4  LATENCY ORDERING: 0 < compute-schedule latency <= critical-path latency. The compute
       schedule is the memory-free lower bound, so it can never exceed the memory-inclusive path.
  SM5  TRAFFIC LOWER BOUND: DRAM traffic is at least the distinct tensor volume. No mapping can
       fetch less data than the workload contains.
  SM6  The mapping's chip factors match the reported chip utilization, so an architecture's own
       grid is genuinely exercised rather than collapsing to one chip. Derived from the mapping
       file, so it cannot drift.
  SM7  DETERMINISM across all of them: a second run is byte-identical.

MAERI's ADDER TREE (SM8-SM10). `adder_cycle` and `adder_energy` price the array-level reduction
network. They were read by the code and declared by no config, and the only architecture with an
[adder_tree] had no gate at all, so the network's own cost was unreachable twice over. maeri's
mapping now puts C SPATIALLY across 64 PEs -- MAERI's defining feature is that tree, and with C
temporal its fan-in is 1 and the whole path early-returns -- and maeri_adder declares the two costs.

  SM8  The tree is exercised: with the costs declared the output interconnection axis is non-zero,
       and with both at 0 it is exactly 0.
  SM9  HAND IDENTITY, energy: an N-leaf tree performs N-1 additions per output, so
         distinct outputs (K*B = 4,096) x (fan_in - 1 = 63) x adder_energy
       This is the array-level mirror of the per-MAC lane-tree identity in validation/reduction RD2,
       and it distinguishes N-1 from N just as that one does.
  SM10 HAND IDENTITY, cycle: the pipeline fill is ceil(log2(fan_in)) = 6 levels per write-back, and
       the model issues one write-back per element, so the delta over the zero-cost run is
         262,144 write-backs x 6 levels x adder_cycle
       linear in the knob, with adder_energy moving no cycle. NOTE the per-ELEMENT fill: a tree fill
       is charged for every element rather than amortized over the stream, which is an upper bound
       on the reduction latency. It is pinned here as the current contract, the same way the
       systolic fold bubble was after it was found to be charged per weight element.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import perturb_lib as pl                                            # noqa: E402

WORKLOAD, MAPPING = "gemm_64x64x64", "ws"
# gemmini is included as the reference architecture: it is the one with external validation, so if
# a cross-architecture invariant fails for it the invariant is wrong, not the architecture.
ARCHITECTURES = ("tpu", "tpuv3", "simba", "maeri", "eyerissv2", "fsd", "gemmini")
K, B, C = 64, 64, 64
EXPECTED_MACS = K*B*C
# Distinct tensor volume at 1 byte/element (every fixture uses int8 operands and output).
DISTINCT_BYTES = K*C + C*B + K*B
ROWS = ("MAC (compute+format)", "PE (local buffer)", "PE array", "Global buffer",
        "Multi-chip (NoP)", "DRAM")
# SM8-SM10: MAERI's array-level reduction tree.
ADDER_FIXTURE = "maeri_adder"
ADDER_CYCLE, ADDER_ENERGY = 2.0, 0.6      # [adder_tree] in maeri_adder.cfg
ADDER_FAN_IN = 64                         # C spatial across 64 PEs
ADDER_OUTPUTS = K*B                       # distinct outputs the tree reduces into
ADDER_DEPTH = 6                           # ceil(log2(64))
ADDER_WRITE_BACKS = K*B*C                 # one tree fill per element (see SM10)


def measure(text: str) -> dict:
    summary = text.split("Energy summary", 1)[1].split("Power summary", 1)[0]
    rows = []
    for row in ROWS:
        found = re.search(rf"\* {re.escape(row)}\s+[\d.]+\s+[\d.]+\s+([\d.]+)", summary)
        if not found:
            raise SystemExit(f"missing energy row {row!r}")
        rows.append(float(found.group(1)))
    dram = text.split("Multi-chip <-> DRAM transactions", 1)[1]
    return {
        "macs": int(re.search(r"# of computations\s*:\s*(\d+)", text).group(1)),
        "schedule": float(re.search(r"Compute-schedule latency\s*:\s*([\d.]+)", text).group(1)),
        "critical": float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", text).group(1)),
        "chip_utilization": float(re.search(r"Chip utilization\s*:\s*([\d.]+)", text).group(1)),
        "pe_utilization": float(re.search(r"PE array utilization\s*:\s*([\d.]+)", text).group(1)),
        "rows": rows,
        "dynamic": float(re.search(r"(?:Layer|Network) dynamic energy\s*:\s*([\d.]+)",
                                   summary).group(1)),
        "static": float(re.search(r"(?:Layer|Network) static energy\s*:\s*([\d.]+)",
                                  summary).group(1)),
        "total": float(re.search(r"(?:Layer|Network) total energy\s*:\s*([\d.]+)",
                                 summary).group(1)),
        "dram_transactions": sum(int(x) for x in re.findall(r"/(\d+)\n", dram)[:3]),
    }


def adder_axes(text: str):
    """The PE-array output interconnection cycle and energy -- where the tree's cost lands."""
    section = text.split("PE array result", 1)[1].split("Global buffer result", 1)[0]
    cycle = [float(v) for _, v in re.findall(
        r"(Input data|Weight|Output data)\s*:\s*([\d.]+) cycles",
        section.split("Interconnection cycle", 1)[1])[:3]]
    energy = [float(v) for _, v in re.findall(
        r"(Input data|Weight|Output data)\s*:\s*([\d.]+) \S+",
        section.split("Interconnection energy", 1)[1])[:3]]
    return {"cycle": cycle[2], "energy": energy[2]}


def mapping_chip_factors(target: str) -> int:
    """Active chips the mapping asks for: the product of the CHIPS_X and CHIPS_Y factors."""
    path = pl.MAPPINGS / target / WORKLOAD / f"{MAPPING}.map"
    active = 1
    for line in path.read_text(encoding="utf-8").splitlines():
        found = re.match(r"^(CHIPS_[XY])\s*=\s*(.+)$", line.split("#")[0].strip())
        if not found:
            continue
        # Only the seven loop dimensions are factors: K, B, P, Q, C, R, S. The H and W fields are
        # 0 placeholders for a GEMM and GROUP/STRIDE are not tiling factors, so multiplying the
        # whole line gives 0.
        fields = [f.strip() for f in found.group(2).split(",") if f.strip()]
        for field in fields[:7]:
            active *= int(field)
    return active


def dram_bitwidth(target: str) -> int:
    text = (pl.ACCELERATORS / f"{target}.cfg").read_text(encoding="utf-8")
    section = text.split("[dram]", 1)[1]
    return int(re.search(r"^bitwidth\s*=\s*(\d+)", section, re.M).group(1))


def main() -> int:
    failures = []
    state, repeated = {}, {}
    for target in ARCHITECTURES:
        text = pl.run(target, WORKLOAD, MAPPING)
        if text is None:
            failures.append(f"SM1 {target}: does not run on the synthetic {WORKLOAD} mapping")
            continue
        state[target] = measure(text)
        again = pl.run(target, WORKLOAD, MAPPING)
        repeated[target] = again

    for target, current in state.items():
        # SM2
        if current["macs"] != EXPECTED_MACS:
            failures.append(f"SM2 {target}: {current['macs']} MACs, expected {EXPECTED_MACS} "
                            f"(K*B*C = {K}*{B}*{C}); the mapping does not cover the workload "
                            f"exactly, so every other number for it is measuring the wrong problem")
        # SM3
        if abs(current["total"] - (current["dynamic"] + current["static"])) > 0.5:
            failures.append(f"SM3 {target}: total {current['total']} != dynamic + static "
                            f"{current['dynamic'] + current['static']}")
        if abs(current["total"] - sum(current["rows"])) > 0.5:
            failures.append(f"SM3 {target}: total {current['total']} != the sum of the six "
                            f"component subtotals {sum(current['rows'])}")
        # SM4
        if current["schedule"] <= 0.0:
            failures.append(f"SM4 {target}: compute-schedule latency is {current['schedule']}")
        if current["critical"] < current["schedule"] - 1e-9:
            failures.append(f"SM4 {target}: critical-path {current['critical']} < compute-schedule "
                            f"{current['schedule']}; the memory-free schedule is a LOWER bound and "
                            f"cannot exceed the memory-inclusive path")
        # SM5
        floor = -(-DISTINCT_BYTES*8 // dram_bitwidth(target))
        if current["dram_transactions"] < floor:
            failures.append(f"SM5 {target}: {current['dram_transactions']} DRAM transactions is "
                            f"below the {floor} needed to move the workload's {DISTINCT_BYTES} "
                            f"distinct bytes over a {dram_bitwidth(target)}-bit link")
        # SM6
        active = mapping_chip_factors(target)
        text = (pl.ACCELERATORS / f"{target}.cfg").read_text(encoding="utf-8")
        chips = int(re.search(r"^num_chips\s*=\s*(\d+)", text, re.M).group(1))
        expected = 100.0*active/chips
        if abs(current["chip_utilization"] - expected) > 0.15:
            failures.append(f"SM6 {target}: chip utilization {current['chip_utilization']}% != "
                            f"{expected:.1f}% ({active} of {chips} chips asked for by the mapping)")
        if current["pe_utilization"] <= 0.0:
            failures.append(f"SM6 {target}: PE-array utilization is 0, so the array is not being "
                            f"exercised at all")
        # SM7
        if repeated[target] != pl.run(target, WORKLOAD, MAPPING):
            failures.append(f"SM7 {target}: two runs of the same config differ")

    # ---- SM8-SM10: MAERI's adder tree -------------------------------------------------------
    import shutil                                                          # noqa: E402

    def adder_variant(label, edits):
        name = f"__sm_{label}"
        pl.variant(ADDER_FIXTURE, name, edits)
        shutil.rmtree(pl.MAPPINGS / name, ignore_errors=True)
        shutil.copytree(pl.MAPPINGS / ADDER_FIXTURE, pl.MAPPINGS / name)
        text = pl.run(name, WORKLOAD, MAPPING)
        pl.discard(name)
        if text is None:
            raise SystemExit(f"the {label} adder-tree variant does not run")
        return adder_axes(text)

    def set_key(key, value):
        return lambda text: re.sub(rf"^{re.escape(key)} = .*$", f"{key} = {value}", text,
                                   count=1, flags=re.M)

    def compose(*edits):
        def apply(text):
            for edit in edits:
                text = edit(text)
            return text
        return apply

    adder = adder_variant("base", lambda text: text)
    adder_free = adder_variant("free", compose(set_key("adder_cycle", "0"),
                                               set_key("adder_energy", "0")))
    adder_energy_2x = adder_variant("e2", set_key("adder_energy", str(ADDER_ENERGY*2)))
    adder_cycle_2x = adder_variant("c2", set_key("adder_cycle", str(ADDER_CYCLE*2)))

    # SM8
    if adder["energy"] <= 0.0:
        failures.append(f"SM8: the adder tree carries no energy ({adder['energy']}); with C "
                        f"temporal its fan-in is 1 and the whole path early-returns, so the mapping "
                        f"no longer exercises the reduction network")
    if adder_free["energy"] != 0.0:
        failures.append(f"SM8: with both adder costs at 0 the tree still charges "
                        f"{adder_free['energy']}")
    # SM9
    expected_energy = ADDER_OUTPUTS*(ADDER_FAN_IN - 1)*ADDER_ENERGY
    if abs(adder["energy"] - expected_energy) > 0.5:
        failures.append(f"SM9: adder energy {adder['energy']} != hand {expected_energy} "
                        f"({ADDER_OUTPUTS} outputs x {ADDER_FAN_IN - 1} additions x "
                        f"{ADDER_ENERGY}); {ADDER_OUTPUTS*ADDER_FAN_IN*ADDER_ENERGY} would mean N "
                        f"additions for an N-leaf tree instead of N-1")
    if abs(adder_energy_2x["energy"] - 2*adder["energy"]) > 0.5:
        failures.append(f"SM9: doubling adder_energy gave {adder_energy_2x['energy']}, expected "
                        f"{2*adder['energy']}")
    if adder_energy_2x["cycle"] != adder["cycle"]:
        failures.append(f"SM9: adder_energy moved the cycle axis ({adder['cycle']} -> "
                        f"{adder_energy_2x['cycle']})")
    # SM10
    expected_cycle = ADDER_WRITE_BACKS*ADDER_DEPTH*ADDER_CYCLE
    if abs((adder["cycle"] - adder_free["cycle"]) - expected_cycle) > 0.5:
        failures.append(f"SM10: the adder cycle delta {adder['cycle'] - adder_free['cycle']} != "
                        f"hand {expected_cycle} ({ADDER_WRITE_BACKS} write-backs x {ADDER_DEPTH} "
                        f"tree levels x {ADDER_CYCLE})")
    doubled_delta = adder_cycle_2x["cycle"] - adder_free["cycle"]
    if abs(doubled_delta - 2*expected_cycle) > 0.5:
        failures.append(f"SM10: doubling adder_cycle gave a delta of {doubled_delta}, expected "
                        f"{2*expected_cycle}")
    if adder_cycle_2x["energy"] != adder["energy"]:
        failures.append(f"SM10: adder_cycle moved the energy axis ({adder['energy']} -> "
                        f"{adder_cycle_2x['energy']})")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"{'architecture':>13} {'MACs':>9} {'sched':>8} {'critical':>11} {'PE%':>6} {'chip%':>6} "
          f"{'total E':>14} {'DRAM txns':>10}")
    for target in ARCHITECTURES:
        current = state[target]
        print(f"{target:>13} {current['macs']:>9} {current['schedule']:>8.0f} "
              f"{current['critical']:>11.0f} {current['pe_utilization']:>6.1f} "
              f"{current['chip_utilization']:>6.1f} {current['total']:>14.1f} "
              f"{current['dram_transactions']:>10}")
    print(f"SM1 all {len(ARCHITECTURES)} shipped architectures run on a synthetic GEMM mapping ok")
    print(f"SM2 every architecture performs exactly {EXPECTED_MACS} MACs ok")
    print("SM3 layer total == dynamic + static == the sum of the six component subtotals ok")
    print("SM4 0 < compute-schedule latency <= critical-path latency ok")
    print(f"SM5 DRAM traffic covers the workload's {DISTINCT_BYTES} distinct bytes ok")
    print("SM6 the mapping's chip factors match the reported chip utilization; the array is live ok")
    print("SM7 every architecture is deterministic across two runs ok")
    print(f"SM8 MAERI's adder tree is exercised ({adder['energy']:.1f}); free when unpriced ok")
    print(f"SM9 adder energy == {ADDER_OUTPUTS} outputs x {ADDER_FAN_IN - 1} additions x "
          f"{ADDER_ENERGY} = {expected_energy:.1f} (N-1, not N) ok")
    print(f"SM10 adder cycle delta == {ADDER_WRITE_BACKS} x {ADDER_DEPTH} levels x {ADDER_CYCLE} = "
          f"{expected_cycle:.0f}; the two costs are independent axes ok")
    print("ALL ARCHITECTURE CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
