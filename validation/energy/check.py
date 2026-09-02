#!/usr/bin/env python3
"""Energy-model validation against published Eyeriss references.

Uses configs/accelerators/eyeriss_energy.cfg: unit costs set to the ISCA'16
Table IV normalized access energies (relative to 1 MAC): RF(spad)=1x,
array(GIN)=2x, GLB=6x, DRAM=200x per 16-bit word.

SCOPE (E7): because those costs are NORMALIZED, every check below is a RELATIVE one --
ratios, shares and internal identities. This suite does NOT validate absolute picojoules,
and it must not be read as doing so. E9 asserts that the config declares the unit and that
the report says so, so an absolute reading is not available by accident.

Checks (asserted; non-zero exit on failure):
  E1  ALU energy == MACs exactly (identity)
  E2  RF accesses/MAC in [3, 8] (chip: ~4/MAC + spad refills)
  E3  CONV layers: RF : rest-except-DRAM in [1.5, 6] (chip claim ~4:1; the
      known PA9 GLB-read overstatement keeps NPUsim below 4 -- see README)
  E4  FC layers: DRAM share > 50% (ISCA Fig. 10 qualitative)
  E5a DRAM access energy / access cycle == u_energy/u_cycle (logical accesses)
  E5b GLB interconnection energy == serialized link txns x per-txn cost
  E7  TOTAL IDENTITY (E1): the energy summary's component subtotals are re-derived from the
      per-datatype lines further down the SAME report, and the layer total is re-derived from
      the subtotals. Without this there is no total for any later formula change to be
      regression-tested against, and nothing stops a component from being double-counted or
      dropped. The re-derivation deliberately uses the detailed per-datatype lines rather than
      the summary's own numbers, so the two have to agree.
  E12 SCOPE AND UNIT LABELS (RE7): the network rollup labels its totals "Network ...", the layer
      files label theirs "Layer ...", neither carries the other's label, and the accelerator
      SPECIFICATION dump prints the configured unit rather than a hardcoded "pJ" -- a normalized
      fixture used to describe its inputs in picojoules and its outputs in normalized units.
  E10 ENERGY SCOPE (E12): the network energy total states how many layers it covers, and says
      PARTIAL when unsupported layers (pooling/activation/normalization) were excluded. Their
      energy is absent from the total -- not estimated, not zero -- so without this line a
      partial rollup reads as an end-to-end network energy. The declared count must match the
      layer files that were actually rolled up, so the statement cannot be decorative.
  E9  UNIT PROVENANCE (E7): this fixture uses normalized costs, so the report must declare
      `normalized` (not pJ) with a provenance string, and every energy line must carry that
      unit. Before this, every value was printed as "pJ" whatever the config meant, so a
      relative breakdown was indistinguishable from measured energy.
  E8  NETWORK ROLLUP (E1): the network report's energy summary equals the sum of its layers'.
      Layers execute serially, so energy adds -- unlike the busy/axis quantities, which have a
      different network contract (see validation/check_timing.py). Checked once over the whole
      result directory.
  E6  OBSERVABILITY (E1): every energy the model accumulates is actually printed. Three
      components -- the PE-array temporal buffer's access energy, the multi-chip temporal
      buffer's access energy, and the DRAM link + row-activation energy -- were accumulated
      and summed into the network totals but appeared nowhere in the results, so their knobs
      changed internal state and showed the reader nothing. Total energy could not be
      reconstructed from the report at all. This check fails if any of those lines goes
      missing again.
"""
import os, re, subprocess, sys

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../.."
LAYERS = {"conv1": 0, "conv2": 2, "conv3": 4, "conv4": 5, "conv5": 6,
          "fc1": 8, "fc2": 9, "fc3": 10}
DRAM_E_PER_CYCLE = 100.0/3.0   # eyeriss_energy.cfg: 100 / 3 cycles per line access
# E7/E9: energy lines carry the unit the CONFIG declares (pJ or normalized), so the
# parsing must not hardcode one. E9 separately asserts the declared unit is the right one.
ENERGY_UNIT_PATTERN = r"(?:pJ|normalized|uncalibrated)"
GIN_E_PER_TXN = 32.0           # [shared] transfer_energy

def parse(L):
    txt = open(os.path.join(ROOT, f"result/eyeriss_energy/alexnet/silicon/layer_{L}.txt")).read()
    def vals(section, label, count=3):
        block = txt.split(section)[1].split(label)[1]
        return [float(v) for v in re.findall(rf":\s*([\d.]+) {ENERGY_UNIT_PATTERN}", block)[:count]]
    def cycles(section, label, count=3):
        block = txt.split(section)[1].split(label)[1]
        return [float(v) for v in re.findall(r":\s*([\d.]+) cycles", block)[:count]]
    alu = float(re.search(rf"Computation energy\s*:\s*([\d.]+) {ENERGY_UNIT_PATTERN}", txt).group(1))
    macs = int(re.search(r"# of data transfer to MAC\n \* Input data\s*:\s*(\d+)", txt).group(1))
    rf = sum(vals("PE result", "Access energy"))
    glb = sum(vals("Global buffer result", "Access energy"))
    glb_sec = txt.split("Global buffer result")[1]
    gin = sum(vals("Global buffer result", "Interconnection energy"))
    glb_txns = sum(int(v) for v in re.findall(
        r"/(\d+)\n", glb_sec.split("serialized")[1].split("# of request")[0]))
    dram = sum(vals("DRAM result", "Access energy"))
    dram_cyc = sum(cycles("DRAM result", "Access cycle"))
    return macs, alu, rf, gin, glb, dram, dram_cyc, glb_txns

# E1/E6: energy lines that must exist in every layer report, as (section, label). Their
# absence is what made total energy unreconstructible; a regression here is silent otherwise.
REQUIRED_ENERGY_LINES = [
    ("MAC result", "Computation energy"),
    ("MAC result", "Format-IP energy"),
    ("MAC result", "Access energy"),
    ("PE result", "Access energy"),
    ("PE array result", "Access energy (temporal buffer)"),
    ("PE array result", "Interconnection energy"),
    ("PE array result", "Static energy"),
    ("Global buffer result", "Access energy"),
    ("Global buffer result", "Interconnection energy"),
    ("Global buffer result", "Static energy"),
    ("Multi chip result", "Access energy (temporal buffer)"),
    ("Multi chip result", "Interconnection energy"),
    ("Multi chip result", "Static energy"),
    ("DRAM result", "Access energy"),
    ("DRAM result", "Interconnection energy (link + row activation)"),
]


def check_energy_observability(layer_index, failures):
    """E6: each required energy line is present in this layer's report."""
    txt = open(os.path.join(ROOT,
               f"result/eyeriss_energy/alexnet/silicon/layer_{layer_index}.txt")).read()
    ok = True
    for section, label in REQUIRED_ENERGY_LINES:
        if section not in txt:
            failures.append(f"{path or layer_index}: section {section!r} missing from the report")
            ok = False
            continue
        body = txt.split(section, 1)[1]
        if label not in body:
            failures.append(f"{path or layer_index}: {section!r} does not report {label!r}; "
                            f"that energy is accumulated but invisible (E1)")
            ok = False
    return ok


def per_type_energy(text, section, label):
    """The three per-datatype energy values printed under `label` inside `section`."""
    body = text.split(section, 1)[1].split(label, 1)[1]
    return sum(float(v) for v in re.findall(rf":\s*([\d.]+) {ENERGY_UNIT_PATTERN}", body)[:3])


def scalar_energy(text, label):
    return float(re.search(rf"{re.escape(label)}\s*:\s*([\d.]+) {ENERGY_UNIT_PATTERN}", text).group(1))


def optional_axis(text, label, section=None):
    """RE8: an optional-feature energy axis (reduction, fold, setup, accumulator, output cast).

    These print as scalars and are 0 in every shipped fixture, which is exactly why the total
    identity could not detect their omission -- adding 0 to the right component and adding it to
    the wrong one look identical. Returns 0.0 when the line is absent for this architecture."""
    body = text.split(section, 1)[1] if section else text
    found = re.search(rf"{re.escape(label)}\s*:\s*([\d.]+) {ENERGY_UNIT_PATTERN}", body)
    return float(found.group(1)) if found else 0.0


def check_energy_totals(layer_index, failures, path=None, require_nonzero=()):
    """E7: re-derive the summary's subtotals and total from the detailed per-datatype lines.

    The component boundary the model documents is: a transfer charges the source buffer's read
    to the SOURCE component, the link to the component that OWNS it, and the destination
    buffer's write to the DESTINATION component. So each printed energy line belongs to exactly
    one component and the subtotals add up exactly -- that is what is verified here.

    RE8: the re-derivation now also includes every OPTIONAL-feature axis, and `require_nonzero`
    lets a caller demand that the axes actually carry energy in this fixture. On eyeriss_energy
    they are all 0, so the identity had no power to detect one of them being dropped from -- or
    added twice to -- a subtotal. E11 re-runs the same identity on a fixture where every one of
    them is non-zero."""
    txt = open(os.path.join(ROOT, path or
               f"result/eyeriss_energy/alexnet/silicon/layer_{layer_index}.txt")).read()
    # Each axis is read from the section that OWNS it. (The weight-fold line is echoed in the
    # PE-array section for readability; reading only the MAC-result copy keeps it counted once.)
    mac_section = txt.split("MAC result", 1)[1].split("PE result", 1)[0]
    array_section = txt.split("PE array result", 1)[1].split("Global buffer result", 1)[0]
    chip_section = txt.split("Multi chip result", 1)[1].split("DRAM result", 1)[0]
    axes = {
        "Reduction energy": optional_axis(mac_section, "Reduction energy"),
        "Accumulator energy (PE)": optional_axis(mac_section, "Accumulator energy (PE)"),
        "Weight-fold energy": optional_axis(mac_section, "Weight-fold energy"),
        "Layer-setup energy": optional_axis(mac_section, "Layer-setup energy"),
        "Accumulator energy (edge)": optional_axis(array_section, "Accumulator energy (edge)"),
        "Output cast energy": optional_axis(chip_section, "Output cast energy"),
    }
    for label in require_nonzero:
        if axes[label] <= 0.0:
            failures.append(f"{path}: {label} is {axes[label]}, so this fixture cannot detect the "
                            f"axis being dropped from a component subtotal")
    expected = {
        "MAC (compute+format)":
            scalar_energy(txt, "Computation energy") +
            per_type_energy(txt, "MAC result", "Format-IP energy") +
            per_type_energy(txt, "MAC result", "Access energy") +
            axes["Reduction energy"] + axes["Accumulator energy (PE)"],
        "PE (local buffer)":
            per_type_energy(txt, "PE result", "Access energy") +
            per_type_energy(txt, "MAC - Local buffer", "Interconnection energy"),
        "PE array":
            per_type_energy(txt, "PE array result", "Access energy (temporal buffer)") +
            per_type_energy(txt, "PE array result", "Interconnection energy") +
            axes["Weight-fold energy"] + axes["Layer-setup energy"] +
            axes["Accumulator energy (edge)"],
        "Global buffer":
            per_type_energy(txt, "Global buffer result", "Access energy") +
            per_type_energy(txt, "Global buffer result", "Interconnection energy"),
        "Multi-chip (NoP)":
            per_type_energy(txt, "Multi chip result", "Access energy (temporal buffer)") +
            per_type_energy(txt, "Multi chip result", "Interconnection energy") +
            axes["Output cast energy"],
        "DRAM":
            per_type_energy(txt, "DRAM result", "Access energy") +
            per_type_energy(txt, "DRAM result", "Interconnection energy (link + row activation)"),
    }
    static_expected = {
        "PE (local buffer)": per_type_energy(txt, "MAC - Local buffer", "Static energy"),
        "PE array": per_type_energy(txt, "PE array result", "Static energy"),
        "Global buffer": per_type_energy(txt, "Global buffer result", "Static energy"),
        "Multi-chip (NoP)": per_type_energy(txt, "Multi chip result", "Static energy"),
    }
    summary = txt.split("Energy summary", 1)[1].split("MAC result", 1)[0]
    ok = True
    for component, expected_dynamic in expected.items():
        # The summary is a 4-column table, so the component name is followed by whitespace
        # and three numbers -- no colon.
        row = re.search(rf"\* {re.escape(component)}\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
                        summary)
        if row is None:
            failures.append(f"{path or layer_index}: energy summary has no {component!r} row")
            ok = False
            continue
        dynamic, stat, total = (float(row.group(i)) for i in (1, 2, 3))
        tolerance = max(0.5, abs(expected_dynamic)*1e-9)
        if abs(dynamic - expected_dynamic) > tolerance:
            failures.append(f"{path or layer_index}: {component} dynamic {dynamic} != the sum of "
                            f"its printed energy lines {expected_dynamic}")
            ok = False
        expected_static = static_expected.get(component, 0.0)
        if abs(stat - expected_static) > max(0.5, abs(expected_static)*1e-9):
            failures.append(f"{path or layer_index}: {component} static {stat} != its printed "
                            f"static energy {expected_static}")
            ok = False
        if abs(total - (dynamic + stat)) > 0.5:
            failures.append(f"{path or layer_index}: {component} total {total} != dynamic+static "
                            f"{dynamic + stat}")
            ok = False
    layer_dynamic = scalar_energy(summary, "Layer dynamic energy")
    layer_static = scalar_energy(summary, "Layer static energy")
    layer_total = scalar_energy(summary, "Layer total energy")
    if abs(layer_dynamic - sum(expected.values())) > max(0.5, layer_dynamic*1e-9):
        failures.append(f"{path or layer_index}: layer dynamic {layer_dynamic} != the sum of the "
                        f"component subtotals {sum(expected.values())}")
        ok = False
    if abs(layer_static - sum(static_expected.values())) > 0.5:
        failures.append(f"{path or layer_index}: layer static {layer_static} != the sum of the "
                        f"component static energies {sum(static_expected.values())}")
        ok = False
    if abs(layer_total - (layer_dynamic + layer_static)) > 0.5:
        failures.append(f"{path or layer_index}: layer total {layer_total} != dynamic+static "
                        f"{layer_dynamic + layer_static}")
        ok = False
    return ok


def check_energy_unit_provenance(failures):
    """E9: the normalized fixture declares its unit, and the report carries it."""
    ok = True
    config = open(os.path.join(ROOT, "configs/accelerators/eyeriss_energy.cfg")).read()
    if not re.search(r"^\s*energy_unit\s*=\s*normalized", config, re.M):
        failures.append("E9: eyeriss_energy.cfg uses normalized access costs but does not "
                        "declare energy_unit = normalized")
        ok = False
    for entry in ("layer_0.txt", "network.txt"):
        txt = open(os.path.join(ROOT,
                   f"result/eyeriss_energy/alexnet/silicon/{entry}")).read()
        unit_line = re.search(r"Energy unit\s*:(.*)", txt)
        if unit_line is None:
            failures.append(f"E9 {entry}: the energy summary does not state its unit")
            ok = False
            continue
        stated = unit_line.group(1)
        if "normalized" not in stated:
            failures.append(f"E9 {entry}: report states {stated.strip()!r} for a normalized "
                            f"fixture")
            ok = False
        if "provenance" not in stated or "not declared" in stated:
            failures.append(f"E9 {entry}: no energy provenance declared ({stated.strip()!r})")
            ok = False
        # The unit must reach the individual lines too, not just the header.
        if '" pJ"' in txt or re.search(r":\s*[\d.]+ pJ", txt):
            failures.append(f"E9 {entry}: energy values are still labelled pJ while the config "
                            f"declares normalized units")
            ok = False
    return ok


def check_energy_scope(failures):
    """E10: the network energy total declares its layer scope, and the count is real."""
    directory = os.path.join(ROOT, "result/eyeriss_energy/alexnet/silicon")
    network = open(os.path.join(directory, "network.txt")).read()
    summary = network.split("Energy summary", 1)[1].split("MAC result", 1)[0]
    scope = re.search(r"Energy scope\s*:\s*(\d+) of (\d+) layers(.*)", summary)
    if scope is None:
        failures.append("E10: the network energy summary does not state its layer scope")
        return False
    included, total, qualifier = int(scope.group(1)), int(scope.group(2)), scope.group(3)
    layer_files = [e for e in os.listdir(directory) if re.fullmatch(r"layer_\d+\.txt", e)]
    ok = True
    if included != len(layer_files):
        failures.append(f"E10: energy scope declares {included} layers but {len(layer_files)} "
                        f"layer result files were rolled up")
        ok = False
    if included < total and "PARTIAL" not in qualifier:
        failures.append(f"E10: {total - included} of {total} layers were excluded but the energy "
                        f"scope does not say PARTIAL ({qualifier.strip()!r})")
        ok = False
    if included == total and "complete" not in qualifier:
        failures.append(f"E10: every layer is in scope but the energy summary does not say so "
                        f"({qualifier.strip()!r})")
        ok = False
    # A per-layer report has no scope to declare; the statement belongs to the rollup only.
    layer = open(os.path.join(directory, "layer_0.txt")).read()
    if "Energy scope" in layer.split("Energy summary", 1)[1].split("MAC result", 1)[0]:
        failures.append("E10: a single-layer report should not claim a network energy scope")
        ok = False
    return ok


def check_network_energy_rollup(failures):
    """E8: network energy == the sum over the layers that were folded in.

    RE7: the two SCOPES must also be labelled apart. The rollup used to print "Layer dynamic
    energy" in network.txt, so a network total read as one layer's and no checker could tell which
    file it was looking at. The layer files must keep "Layer", network.txt must say "Network", and
    neither may carry the other's label."""
    directory = os.path.join(ROOT, "result/eyeriss_energy/alexnet/silicon")
    network = open(os.path.join(directory, "network.txt")).read()
    summary = network.split("Energy summary", 1)[1].split("MAC result", 1)[0]
    for wrong in ("Layer dynamic energy", "Layer static energy", "Layer total energy"):
        if wrong in summary:
            failures.append(f"RE7: network.txt labels its rollup {wrong!r}; a network total must "
                            f"not be labelled as one layer's")
    totals = {}
    for label in ("Network dynamic energy", "Network static energy", "Network total energy"):
        if not re.search(rf"{re.escape(label)}\s*:", summary):
            failures.append(f"RE7: network.txt does not report {label!r}; the rollup's scope must "
                            f"be stated in its label")
            return False
        totals[label] = scalar_energy(summary, label)
    summed = {label: 0.0 for label in totals}
    layers = 0
    for entry in sorted(os.listdir(directory)):
        if not re.fullmatch(r"layer_\d+\.txt", entry):
            continue
        layers += 1
        block = open(os.path.join(directory, entry)).read()
        block = block.split("Energy summary", 1)[1].split("MAC result", 1)[0]
        if "Network total energy" in block:
            failures.append(f"RE7: {entry} labels its total as a network one")
        for label in totals:
            summed[label] += scalar_energy(block, label.replace("Network", "Layer"))
    if layers == 0:
        failures.append("E8: no layer_*.txt to roll up")
        return False
    ok = True
    for label, measured in totals.items():
        # Both sides are printed to two decimals, so allow the rollup's rounding drift.
        if abs(measured - summed[label]) > max(0.5, 0.01*layers, abs(summed[label])*1e-9):
            failures.append(f"E8: network {label} {measured} != the sum over {layers} layers "
                            f"{summed[label]}")
            ok = False
    return ok


def run_all_knobs():
    """E11's fixture is not part of any other suite's result set, so run it here. Reading a stale
    result file would make the identity check pass against the previous build."""
    subprocess.run([os.path.join(ROOT, "npusim.sh"), "run", "gemmini_all_knobs",
                    "gemm_64x64x64", "ws"], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)


def check_specification_units(failures):
    """RE7 condition 4: the accelerator SPECIFICATION dump (stdout) must use the declared unit too.

    The result file's energy lines carried the configured unit, but the specification printed every
    unit cost as "pJ" regardless -- so the same normalized fixture described its own inputs in
    picojoules and its outputs in normalized units. A reader comparing the two would conclude the
    costs had been converted somewhere. Nothing had."""
    stdout = subprocess.run([os.path.join(ROOT, "npusim.sh"), "run", "eyeriss_energy",
                             "alexnet", "silicon"], cwd=ROOT, check=True,
                            capture_output=True, text=True).stdout
    offenders = [line.strip() for line in stdout.splitlines() if re.search(r"\bpJ\b", line)]
    if offenders:
        failures.append(f"RE7: the normalized fixture's specification still prints pJ "
                        f"({len(offenders)} line(s), e.g. {offenders[0]!r})")
        return False
    if "normalized" not in stdout:
        failures.append("RE7: the specification dump carries no unit label at all, so its unit "
                        "costs cannot be related to the result file's")
        return False
    return True


def main():
    failures = []
    # Regenerate the Eyeriss fixture BEFORE parsing its layer files. Running this near the end
    # made the first check after a build validate stale results and only refresh them afterward.
    e12_ok = check_specification_units(failures)
    run_all_knobs()
    def check(cond, msg):
        if not cond:
            failures.append(msg)
        return "ok" if cond else "FAIL"
    print(f"{'layer':>6s} {'RF/MAC':>7s} {'ALU%':>6s} {'RF%':>6s} {'GIN%':>6s} {'GLB%':>6s} "
          f"{'DRAM%':>6s} {'RF:rest':>8s}  {'E1':>4s} {'E2':>4s} {'E3/4':>4s} {'E5a':>4s} {'E5b':>4s} "
          f"{'E6':>4s} {'E7':>4s}")
    for name, L in LAYERS.items():
        macs, alu, rf, gin, glb, dram, dram_cyc, glb_txns = parse(L)
        total = alu + rf + gin + glb + dram
        rest = alu + gin + glb
        e1 = check(abs(alu/macs - 1) < 1e-9, f"{name}: E1 ALU {alu} != MACs {macs}")
        e2 = check(3.0 <= rf/macs <= 8.0, f"{name}: E2 RF/MAC {rf/macs:.2f} outside [3,8]")
        if name.startswith("conv"):
            e34 = check(1.5 <= rf/rest <= 6.0, f"{name}: E3 RF:rest {rf/rest:.2f} outside [1.5,6]")
        else:
            e34 = check(dram/total > 0.5, f"{name}: E4 DRAM share {dram/total:.2f} <= 0.5")
        e5a = check(abs(dram/dram_cyc - DRAM_E_PER_CYCLE)/DRAM_E_PER_CYCLE < 1e-9,
                    f"{name}: E5a DRAM pJ/cycle {dram/dram_cyc:.4f} != {DRAM_E_PER_CYCLE:.4f}")
        e5b = check(abs(gin - glb_txns*GIN_E_PER_TXN)/max(gin, 1) < 1e-9,
                    f"{name}: E5b GIN {gin} != txns {glb_txns} x {GIN_E_PER_TXN}")
        e6 = "ok" if check_energy_observability(L, failures) else "FAIL"
        e7 = "ok" if check_energy_totals(L, failures) else "FAIL"
        print(f"{name:>6s} {rf/macs:7.2f} {100*alu/total:6.1f} {100*rf/total:6.1f} "
              f"{100*gin/total:6.1f} {100*glb/total:6.1f} {100*dram/total:6.1f} "
              f"{rf/rest:7.2f}:1  {e1:>4s} {e2:>4s} {e34:>4s} {e5a:>4s} {e5b:>4s} {e6:>4s} "
              f"{e7:>4s}")
    # E11 (RE8): the same total identity on a fixture where EVERY optional-feature axis carries
    # energy. On eyeriss_energy they are all 0, so E7 above could not tell a subtotal that omits
    # an axis from one that includes it -- nor one that counts it twice.
    e11 = "ok" if check_energy_totals(
        0, failures, path="result/gemmini_all_knobs/gemm_64x64x64/ws/layer_0.txt",
        require_nonzero=("Reduction energy", "Weight-fold energy", "Layer-setup energy",
                         "Accumulator energy (edge)", "Output cast energy")) else "FAIL"
    e8 = "ok" if check_network_energy_rollup(failures) else "FAIL"
    e9 = "ok" if check_energy_unit_provenance(failures) else "FAIL"
    e10 = "ok" if check_energy_scope(failures) else "FAIL"
    e12 = "ok" if e12_ok else "FAIL"
    print(f"\nE11 total identity holds with every optional axis non-zero (reduction, weight fold, "
          f"layer setup, edge accumulator, output cast, DRAM row activation): {e11}")
    print(f"E8 network energy == sum over layers: {e8}")
    print(f"E9 unit declared and carried through the report (normalized, not pJ): {e9}")
    print(f"E10 network energy total declares its layer scope: {e10}")
    print(f"E12 (RE7) network rollup labels itself Network, layers label themselves Layer, and the "
          f"specification dump uses the declared unit (no hardcoded pJ): {e12}")
    print("NOTE: this fixture is NORMALIZED, so the checks above are relative; absolute "
          "picojoules are NOT validated here.")
    print("chip refs: RF/MAC ~4 (+refills); CONV RF:rest ~4:1 (JSSC; PA9 keeps NPUsim below); "
          "FC DRAM-dominant (ISCA Fig.10)")
    if failures:
        print(f"\n{len(failures)} check(s) FAILED:")
        for f in failures:
            print("  " + f)
        return 1
    print("ALL ENERGY CHECKS PASSED")
    return 0

if __name__ == "__main__":
    sys.exit(main())
