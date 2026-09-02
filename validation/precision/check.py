#!/usr/bin/env python3
"""Precision-dependent energy regression (E3).

MAC energy used to be ONE scalar multiplied by the active MAC count, with no dependence on the
operand formats. Switching input/weight precision therefore changed the memory traffic and left
the compute energy bit-identical -- an INT4 multiply was priced exactly like an FP16 one, and
nothing in the model or the report treated that as odd.

A config may now declare the cost for the precision it was measured at:
    mac_energy_<input format>_<weight format>_<accumulator format>   e.g. mac_energy_int4_int4_fp32
and when the running combination has no entry the report says UNCALIBRATED instead of presenting
the fallback scalar as a precision-aware number.

RE4 -- the key names all THREE datapath widths. The first version used only the two operands, so
an INT8 x INT8 MAC accumulating into FP32 and one accumulating into FP16 shared a single cost even
though the adder and the accumulator register are a large part of a MAC's energy. A config that
declares only the two-operand key is treated as a PARTIAL calibration: the number is used (it beats
the bare scalar) but the config does not qualify for absolute power (see validation/power).

Three fixtures differ from configs/accelerators/gemmini.cfg only in the operand formats and the
declared MAC energy for those formats (0.28 / 0.56 / 1.68 for int4 / int8 / fp16 -- synthetic
values monotone in operand width, not measured silicon numbers):

    gemmini_precision_int4   gemmini_precision_int8   gemmini_precision_fp16

Checks (asserted; non-zero exit on failure):
  PR1  Same work: the MAC count is IDENTICAL across the three precisions. This is what makes the
       two axes below separable -- neither difference can be attributed to doing more work.
  PR2  COMPUTE axis, hand identity: computation energy == MAC count * the cost DECLARED for that
       precision. Independent of any memory effect.
  PR3  MEMORY axis, per datatype: the INPUT and WEIGHT transaction counts scale with the operand
       bit width (int4 : int8 : fp16 = 1 : 2 : 4), while the OUTPUT count is UNCHANGED -- these
       fixtures move only input_format/weight_format and leave output_format at int8. That makes
       PR3 a stronger statement than a global ratio: precision is tracked PER DATATYPE, so a
       mixed-precision config is priced per tensor rather than by one global width.
  PR4  The two axes are genuinely independent: compute energy follows the declared cost even
       though it does NOT follow the bit-width ratio, so a single "scale everything by width"
       rule cannot satisfy PR2 and PR3 at once. (Here fp16's declared cost is 3x int8's while
       its traffic is 2x -- deliberately not proportional.)
  PR5  PROVENANCE: each run names the per-precision key it used and reports itself calibrated;
       the shipped gemmini config, which declares no per-precision cost, must report
       UNCALIBRATED and name the fallback scalar. Without this the fallback is indistinguishable
       from a measured value.
  PR6  RE4 ACCUMULATOR A/B: gemmini_macacc_fp32 and gemmini_macacc_fp16 differ ONLY in the
       accumulator width (INT8 x INT8 operands either way). The MAC count and every operand
       transaction count must be BIT-IDENTICAL, while the compute energy follows the cost declared
       for each accumulator -- the difference the two-operand key could not express.
  PR7  RE4 PARTIAL CALIBRATION: gemmini_macacc_operand_only declares the operand key but not the
       accumulator one. It must (a) still use that cost rather than the bare scalar, (b) report
       itself NOT calibrated and name the accumulator width that is missing, and (c) therefore
       produce the same compute energy as the fp32 fixture while being a different accumulator --
       which is exactly the blind spot, now visible in the report instead of hidden.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKLOAD = ("gemm_64x64x64", "ws")
# format -> (declared MAC energy, operand bits)
PRECISIONS = {"int4": (0.28, 4), "int8": (0.56, 8), "fp16": (1.68, 16)}
FALLBACK_TARGET = "gemmini"
# RE4: accumulator A/B -- INT8 x INT8 operands, accumulator width and its declared cost differ.
ACCUMULATORS = {"gemmini_macacc_fp32": ("fp32", 0.56), "gemmini_macacc_fp16": ("fp16", 0.34)}
OPERAND_ONLY = "gemmini_macacc_operand_only"   # operand key declared, accumulator key not


def run(target: str) -> str:
    layer = ROOT / "result" / target / WORKLOAD[0] / WORKLOAD[1] / "layer_0.txt"
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, WORKLOAD[0], WORKLOAD[1]],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if not layer.exists():
        raise SystemExit(f"missing simulator output: {layer}")
    return layer.read_text(encoding="utf-8")


def measure(text: str) -> dict:
    glb = text.split("Global buffer result", 1)[1]
    glb_table = glb.split("serialized", 1)[1].split("# of request", 1)[0]
    glb_per_type = {}
    for name, value in re.findall(
            r"(Input data|Weight|Output data)\s*:\s*\d+/\d+/(\d+)", glb_table):
        glb_per_type.setdefault(name, int(value))
    dram = text.split("Multi-chip <-> DRAM transactions", 1)[1]
    dram_per_type = {}
    for name, value in re.findall(
            r"(Input data|Weight|Output data)\s*:\s*\d+/\d+/(\d+)", dram):
        dram_per_type.setdefault(name, int(value))
    glb_txns = sum(glb_per_type.values())
    dram_txns = sum(dram_per_type.values())
    basis = re.search(r"MAC energy basis\s*:\s*(\S+)(.*)", text)
    return {
        "macs": int(re.search(r"# of computations\s*:\s*(\d+)", text).group(1)),
        "compute_energy": float(re.search(
            r"Computation energy\s*:\s*([\d.]+) (?:pJ|normalized|uncalibrated)", text).group(1)),
        "glb_txns": glb_txns,
        "dram_txns": dram_txns,
        "glb_per_type": glb_per_type,
        "dram_per_type": dram_per_type,
        "basis": basis.group(1),
        "basis_note": basis.group(2),
        "calibrated": "NOT calibrated" not in basis.group(2),
    }


def main() -> int:
    measured = {fmt: measure(run(f"gemmini_precision_{fmt}")) for fmt in PRECISIONS}
    fallback = measure(run(FALLBACK_TARGET))
    failures = []

    # PR1
    macs = {state["macs"] for state in measured.values()}
    if len(macs) != 1:
        failures.append(f"PR1: MAC counts differ across precisions "
                        f"({ {f: s['macs'] for f, s in measured.items()} }); the compute and "
                        f"memory axes are then not separable")
    mac_count = measured["int8"]["macs"]

    # PR2
    for fmt, (cost, _) in PRECISIONS.items():
        expected = mac_count*cost
        if abs(measured[fmt]["compute_energy"] - expected) > 0.5:
            failures.append(f"PR2 {fmt}: computation energy {measured[fmt]['compute_energy']} != "
                            f"hand {expected} ({mac_count} MACs x declared {cost})")

    # PR3: the operand tensors scale with the declared operand width; the output tensor, whose
    # format these fixtures do NOT change, must stay put.
    for boundary in ("glb_per_type", "dram_per_type"):
        for datatype in ("Input data", "Weight"):
            base = measured["int8"][boundary][datatype]
            for fmt, (_, bits) in PRECISIONS.items():
                expected = base*bits/8
                if abs(measured[fmt][boundary][datatype] - expected) > 0.5:
                    failures.append(
                        f"PR3 {fmt} {boundary} {datatype}: "
                        f"{measured[fmt][boundary][datatype]} != hand {expected:.0f} "
                        f"(int8 baseline {base} scaled by {bits}/8 bits)")
        output_counts = {measured[fmt][boundary]["Output data"] for fmt in PRECISIONS}
        if len(output_counts) != 1:
            failures.append(
                f"PR3 {boundary} Output data: changed with the operand precision "
                f"({ {f: measured[f][boundary]['Output data'] for f in PRECISIONS} }) even though "
                f"output_format is fixed at int8; precision must be tracked per datatype")

    # PR4
    compute_ratio = measured["fp16"]["compute_energy"]/measured["int8"]["compute_energy"]
    traffic_ratio = measured["fp16"]["dram_txns"]/measured["int8"]["dram_txns"]
    if abs(compute_ratio - traffic_ratio) < 1e-6:
        failures.append(f"PR4: the compute and traffic ratios coincide ({compute_ratio:.2f}); the "
                        f"fixture cannot then show that compute energy is an independent axis")

    # PR5
    for fmt in PRECISIONS:
        expected_key = f"mac_energy_{fmt}_{fmt}_fp32"
        if measured[fmt]["basis"] != expected_key:
            failures.append(f"PR5 {fmt}: MAC energy basis is {measured[fmt]['basis']!r}, expected "
                            f"{expected_key!r}")
        if not measured[fmt]["calibrated"]:
            failures.append(f"PR5 {fmt}: a declared per-precision cost is reported as "
                            f"uncalibrated ({measured[fmt]['basis_note'].strip()})")
    if fallback["basis"] != "computation_energy":
        failures.append(f"PR5 {FALLBACK_TARGET}: basis is {fallback['basis']!r}, expected the "
                        f"fallback scalar 'computation_energy'")
    if fallback["calibrated"]:
        failures.append(f"PR5 {FALLBACK_TARGET}: declares no per-precision MAC energy but does "
                        f"not report itself uncalibrated ({fallback['basis_note'].strip()})")

    # PR6 -- RE4 accumulator A/B
    accum = {name: measure(run(name)) for name in ACCUMULATORS}
    partial = measure(run(OPERAND_ONLY))
    for name, (fmt, cost) in ACCUMULATORS.items():
        state = accum[name]
        expected_key = f"mac_energy_int8_int8_{fmt}"
        if state["basis"] != expected_key:
            failures.append(f"PR6 {name}: MAC energy basis is {state['basis']!r}, expected "
                            f"{expected_key!r} -- the key must name the accumulator width")
        if not state["calibrated"]:
            failures.append(f"PR6 {name}: a declared three-part cost is reported as uncalibrated "
                            f"({state['basis_note'].strip()})")
        expected = state["macs"]*cost
        if abs(state["compute_energy"] - expected) > 0.5:
            failures.append(f"PR6 {name}: computation energy {state['compute_energy']} != hand "
                            f"{expected} ({state['macs']} MACs x declared {cost})")
    fp32, fp16 = accum["gemmini_macacc_fp32"], accum["gemmini_macacc_fp16"]
    if fp32["macs"] != fp16["macs"]:
        failures.append(f"PR6: the accumulator A/B changed the MAC count "
                        f"({fp32['macs']} vs {fp16['macs']}); the compute delta is then not "
                        f"attributable to the accumulator width")
    for boundary in ("glb_per_type", "dram_per_type"):
        for datatype in ("Input data", "Weight", "Output data"):
            if fp32[boundary][datatype] != fp16[boundary][datatype]:
                failures.append(
                    f"PR6 {boundary} {datatype}: operand traffic moved with the accumulator width "
                    f"({fp32[boundary][datatype]} vs {fp16[boundary][datatype]}); only the compute "
                    f"axis may respond to it")
    if fp32["compute_energy"] == fp16["compute_energy"]:
        failures.append(f"PR6: both accumulator widths give the same compute energy "
                        f"({fp32['compute_energy']}); the accumulator datapath is still invisible")

    # PR7 -- partial calibration
    if partial["basis"] != "mac_energy_int8_int8":
        failures.append(f"PR7 {OPERAND_ONLY}: basis is {partial['basis']!r}; the operand key must "
                        f"still be USED, not discarded in favour of the bare scalar")
    if partial["calibrated"]:
        failures.append(f"PR7 {OPERAND_ONLY}: an unpriced accumulator width is reported as "
                        f"calibrated ({partial['basis_note'].strip()})")
    if "fp16" not in partial["basis_note"]:
        failures.append(f"PR7 {OPERAND_ONLY}: the report does not name the accumulator width that "
                        f"is missing a cost ({partial['basis_note'].strip()})")
    if abs(partial["compute_energy"] - fp32["compute_energy"]) > 0.5:
        failures.append(f"PR7 {OPERAND_ONLY}: compute energy {partial['compute_energy']} != the "
                        f"fp32 fixture's {fp32['compute_energy']}; the operand-only key is the "
                        f"same 0.56, so this is the blind spot the report must declare")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"{'precision':>10} {'MACs':>9} {'compute E':>12} {'GLB txns':>9} {'DRAM txns':>10} "
          f"{'basis':>24}")
    for fmt in ("int4", "int8", "fp16"):
        state = measured[fmt]
        print(f"{fmt:>10} {state['macs']:>9} {state['compute_energy']:>12.2f} "
              f"{state['glb_txns']:>9} {state['dram_txns']:>10} {state['basis']:>24}")
    print(f"PR1 MAC count identical at {mac_count} ok")
    print("PR2 compute energy == MACs x declared per-precision cost ok")
    print("PR3 input/weight transactions scale 1:2:4 with operand width; output unchanged "
          "(output_format fixed) ok")
    print(f"PR4 compute ratio {compute_ratio:.2f} != traffic ratio {traffic_ratio:.2f}: the two "
          f"axes are independent ok")
    print("PR5 per-precision key named when declared; fallback reported UNCALIBRATED ok")
    print(f"PR6 fp32 vs fp16 accumulator: compute energy "
          f"{fp32['compute_energy']:.0f} vs {fp16['compute_energy']:.0f} with identical MAC count "
          f"and identical operand traffic ok")
    print(f"PR7 operand-only key is used ({partial['compute_energy']:.0f}) but reported NOT "
          f"calibrated, naming the unpriced fp16 accumulator ok")
    print("ALL PRECISION CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
