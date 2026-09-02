#!/usr/bin/env python3
"""Config-knob liveness and model-invariant regression (dead-knob sweep + metamorphic tests).

WHY THIS EXISTS. Every other gate in validation/ checks an axis that WORKS: it names a quantity,
derives it by hand, and asserts the model agrees. That shape cannot detect the defect this project
has hit once per assessment round -- a setting that is parsed, printed, and then consumed by
nothing. Such a knob produces a BIT-IDENTICAL report however it is set, so no amount of checking
the axes that do work will reveal it. The instances found so far were all found by hand, one at a
time:

    accumulator_format   parsed and printed; an accumulator spill was sized by the OUTPUT datatype
                         instead, so fp32 and fp16 produced identical traffic and energy   (E4)
    row_miss_energy      accumulated into internal state and printed nowhere               (E1a)
    noc_energy           zero in every shipped config, so its multiplier -- wrong by ~17x for
                         multicast -- was never exercised at all                        (E2/RE6)

This gate finds that class MECHANICALLY: perturb each declared setting alone, re-run, and require
at least one reported number to move. A knob that moves nothing is either dead, or not exercised by
the fixture -- and the second case must be declared, with the gate that does exercise it named.

The same perturbation machinery gives the metamorphic checks, which are statements about the whole
accounting graph rather than about one axis, and so cover code no hand identity names:

  KN1  DETERMINISM: the same config run twice produces an identical report.
  KN2  ENERGY SCALE INVARIANCE: multiplying EVERY energy unit cost by 10 multiplies every reported
       energy by exactly 10, and moves no cycle or latency by anything. Energy accounting is linear
       and homogeneous in its unit costs, so any energy that fails to scale is computed from a
       hardcoded constant, and any latency that moves means an energy cost is feeding the timing
       path. One check, whole-model coverage.
  KN3  ZERO COLLAPSE: with every energy unit cost at 0, every reported energy is EXACTLY 0 -- no
       tolerance. KN2 already catches most spurious additive terms, since it compares against
       10x the base rather than checking a ratio; KN3 adds the cases KN2's print-rounding slack can
       absorb (a term small enough that 9x it stays inside the slack) and states the property
       without any tolerance at all, which is worth having for one extra run.
  KN4  DEAD-KNOB SWEEP: see above. Exemptions are declared below, each with its reason and, where
       applicable, the gate that does cover the knob.
  KN5  ENERGY/TIMING SEPARATION, per knob: perturbing an ENERGY cost alone must not move any cycle
       or latency field. KN2 asserts this for a uniform scaling; this asserts it one knob at a time,
       so a single mis-wired cost cannot hide behind the others.
  KN6  ENERGY MONOTONICITY: raising one energy unit cost never LOWERS any reported energy.
  KN7  LATENCY MONOTONICITY: raising one cycle cost never LOWERS the critical-path latency.
  KN9  UNDECLARED-KNOB SCAN. KN4 can only perturb what a config DECLARES, so a knob no config
       mentions is outside its reach entirely -- and that is exactly where the biggest instance was
       hiding: `exist_temporal_buffer` gates the whole PE-array temporal buffer, defaults to OFF,
       and all 48 shipped configs shipped it COMMENTED OUT. Five cost knobs behind it had no
       coverage in any suite, and the sweep could not see them. So this check scans the
       complementary set: every key the code calls get_setting()/get_vector_setting() on, minus
       every key any config declares. Each survivor is classified below -- an alias that falls back
       to a declared key, an optional alternate with a default, or an UNCOVERED cost knob, which is
       a real backlog item rather than a passing check.
  KN11 A DECLARED-BUT-UNUSABLE VALUE IS REJECTED, not silently dropped. section_config_t's
       get_setting() used to return false whenever a present value would not parse into the
       requested C++ type, and every caller reads false as "not declared" -- so the value vanished
       and the default stood, with nothing printed. Two instances had shipped:
         * `[multi_chip] input_size = 4.5` -- the field was `unsigned`, the fractional value failed
           the eof() check, the buffer stayed 0, and eyerissv2.cfg could not run at all.
         * `[dram] transfer_cycle = 2.75:2.75:2.75` in maeri.cfg -- a per-datatype list for a key
           the model reads as a scalar, so the DRAM transfer cost was 0 for that architecture.
       Neither was noticed because neither config had a gate. A bool given "2" behaves the same way,
       which is also why a x2 perturbation of a boolean knob read as a dead knob in KN4. This check
       runs the simulator on a config carrying each such value and requires it to fail, naming the
       key.
  KN10 ENERGY SCHEMA COMPLETENESS. The RE5 schema (utils/energy_units.cc) rejects any key
       containing "energy" that it does not declare, on the grounds that such a key is a typo
       silently leaving a component at zero cost. But that table was first built from a survey of
       the keys the SHIPPED CONFIGS declared -- so it omitted every energy key the code reads that
       no config had ever set, and REJECTED two of them as typos when a new fixture finally tried
       to declare one. The schema has to be complete against what the CODE READS. This check
       compares the two directly.
  KN8  NO STALE EXEMPTIONS: every knob declared exempt below must actually behave as claimed. An
       exemption that has silently become live is a checker asserting nothing, so the tables cannot
       rot the way a hand-maintained allowlist usually does.

Fixture: gemmini_all_knobs -- the config with every optional-feature knob non-zero, which is what
makes a single sweep meaningful. A knob whose feature is switched off in the fixture is dead for an
uninteresting reason.
"""
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import perturb_lib as pl                                            # noqa: E402

BASE, WORKLOAD, MAPPING = "gemmini_all_knobs", "gemm_64x64x64", "ws"
SCALE = 10.0

# Boolean knobs: 0/1, so doubling is out of range and behaves like the original -- worse, a
# stringstream bool parse of "2" FAILS and silently leaves the default in place, which reads as a
# dead knob. The meaningful perturbation is a toggle. (double_buffer is 0 = single, 1 = double --
# NOT a depth. The two exist_temporal_buffer flags have DIFFERENT defaults: pe_array_t's is false,
# multi_chip_t's is true.)
TOGGLE = {("systolic_array", "double_buffer"), ("systolic_array", "edge_accumulation"),
          ("systolic_array", "exist_temporal_buffer"), ("multi_chip", "exist_temporal_buffer")}

# Knobs whose perturbation MUST leave the report untouched. A move here is a defect.
MUST_BE_INVARIANT = {
    ("accelerator", "num_threads"):
        "simulation parallelism, not a model parameter -- a result that depends on the thread "
        "count is non-deterministic, so this knob is REQUIRED to be dead",
}

# Knobs that are dead in THIS fixture for a stated reason. Each names where it is covered, or why
# nothing can cover it.
EXEMPT_DEAD = {
    ("systolic_array", "computation_energy"):
        "overridden by mac_energy_int8_int8_fp32 in this fixture, which is the RE4 precedence "
        "contract; the fallback scalar's own path is covered by validation/precision PR5",
    ("dram", "t_wtr_cycle"):
        "this schedule issues all of its loads and then its stores, so the bus flips read->write "
        "once and write->read never. Covered by validation/dram DR8, which runs the turnaround "
        "identity on gemm_512x512x512 where stores do precede later loads -- a gap this sweep "
        "found, since t_wtr_cycle previously had no coverage anywhere",
}
# `bandwidth` is declared in four sections and is dead in all of them. That is a real modelling
# fact worth stating rather than a bug: a reader may reasonably expect bandwidth to throttle.
_BANDWIDTH_REASON = ("bandwidth is only a SOURCE for deriving bitwidth when bitwidth is absent "
                     "(see derived_link_bitwidth); it is NOT a rate limit in this model, and every "
                     "section of this fixture declares bitwidth explicitly")
for _section in ("systolic_array", "separate", "multi_chip", "dram"):
    EXEMPT_DEAD[(_section, "bandwidth")] = _BANDWIDTH_REASON

# Knobs whose perturbed config cannot RUN, with the constraint that rejects it. A perturbation that
# produces an invalid machine is not evidence about the knob.
EXEMPT_INVALID = {
    ("accelerator", "num_chips"):
        "the mapping pins the chip count; doubling it leaves chip 1 with no work",
    ("systolic_array", "edge_accumulation"):
        "with edge accumulation off, partial sums must fit the PE local buffer (output 256 vs "
        "buffer 64 here); the A/B is covered by validation/accumulator AC9",
    ("multi_chip", "height"):           "the mapping pins the chip grid",
    ("multi_chip", "width"):            "the mapping pins the chip grid",
}

# Non-numeric settings (enums, formats, paths). A generic perturbation is not defined for them;
# each is covered by its own A/B fixture where it matters.
EXEMPT_NON_NUMERIC = {
    "energy_unit": "validation/power PW7", "energy_reference": "validation/power PW7",
    "compression_type": "validation/format", "input_format": "validation/precision",
    "weight_format": "validation/precision", "output_format": "validation/precision",
    "memory_type": "shipped-config A/B", "mac_stationary": "shipped-config A/B",
    "pe_stationary": "shipped-config A/B", "noc": "validation/noc NC5",
    "stationary_type": "shipped-config A/B", "dram_config": "external DRAM model selection",
    "output_dir": "output path, not a model parameter",
}


# KN9: keys the code reads that no config declares, each with what it actually is.
UNDECLARED_ALIAS = {
    "nop":            "falls back to `noc` (multi_chip.cc:154)",
    "nop_cycle":      "falls back to `noc_cycle` (multi_chip.cc:227)",
    "nop_energy":     "falls back to `noc_energy` (multi_chip.cc:230)",
    "num_processors": "alternative spelling of `num_chips`, which every config declares (npu.cc:70)",
}
UNDECLARED_OPTIONAL = {
    # 2026-08-28 (PyTorch frontend / SFU):
    "op_id":            "mapping-SECTION binding key for executable-IR runs -- it lives in mapping "
                        "files, not accelerator configs, so this sweep's declared-set can never "
                        "contain it by construction. Bound and validated (missing/duplicate/mixed/"
                        "wrong-kind all fail fast) by npu_frontend.cc bind_executable_mappings()",
    "strict_profiles":  "SFU validation-strictness toggle, default off: turns a profile-precision "
                        "mismatch into a config error instead of an UNCALIBRATED annotation. "
                        "Moves no modeled number by design, so the liveness sweep cannot see it",
    "profile_precision": "optional SFU profile provenance: names the format the declared "
                        "latency/energy profile describes; undeclared = every invocation is "
                        "marked UNCALIBRATED (the honest default)",
    "num_pes":                  "optional override; the mapping supplies the PE count",
    "data_format":              "shorthand that sets all three datatype formats at once; the "
                                "configs declare input/weight/output_format individually",
    "mxfp_metadata_layout":     "MXFP block-metadata placement; only meaningful for MXFP formats",
    "parameter_order":          "optional loop order override; the configs use stationary_type",
    "pe_array_parameter_order": "same, at the array level",
    "pe_parameter_order":       "same, at the PE level",
}
# Cost knobs that exist in the model and that NOTHING exercises. Not a passing state -- a recorded
# backlog. Each entry says what would close it.
UNDECLARED_UNCOVERED = {
    # OPEN 2026-08-28: softmax operand residency has a "glb" branch that changes the modeled
    # softmax operand traffic (operands staged in the GLB instead of streaming from DRAM), and
    # no config declares it -- every shipped SFU run exercises only the "dram" default. Closing
    # it needs a config that declares softmax_operand_residency = glb plus a gate pinning the
    # DRAM-vs-GLB traffic delta of a softmax layer.
    "softmax_operand_residency": "cost-affecting residency branch (dram|glb); only the dram "
                                 "default is exercised -- see the OPEN note above",
    # CLOSED 2026-08-20: accumulator_spill_cycle and output_cast_cycle used to sit here. The
    # accumulator fixtures now declare both and validation/accumulator AC10-AC13 pin their
    # latency identities, so they are declared knobs and the sweep above reaches them.
    # CLOSED 2026-08-20: format_payload_energy, format_metadata_cycle and format_metadata_energy
    # likewise -- gemmini_format_mxfp declares all four Format-IP costs and validation/format
    # FM7-FM11 separate the payload and metadata terms by measurement.
    # CLOSED 2026-08-20: row_miss_cycle is not a duplicate of t_ras+t_rp but the ALTERNATIVE flat
    # activation model dram_t selects when tRAS/tRP are absent -- gemmini_dram_flat exercises that
    # branch and validation/dram DR9 pins it. lb_static_energy is not an alternative to
    # static_energy either but a second ADDITIVE PE leakage term; gemmini_lb_leakage declares it and
    # validation/leakage LK6 pins the sum. (An earlier note here guessed row_miss_cycle was
    # redundant and should be deleted. Reading the code refuted that.)
    # CLOSED 2026-08-20: adder_cycle and adder_energy. The six shipped architecture configs now have
    # synthetic GEMM mappings and a gate (validation/architectures), and maeri_adder declares the
    # [adder_tree] costs -- SM8-SM10 pin the N-1 addition identity and the tree-fill latency.
}


def schema_energy_keys():
    """The energy keys utils/energy_units.cc declares (KN10)."""
    text = (pl.ROOT / "utils" / "energy_units.cc").read_text(encoding="utf-8")
    table = text.split("g_energy_keys[] = {", 1)[1].split("\n};", 1)[0]
    return set(re.findall(r'\{\s*"([a-z_0-9]+)"', table))


def code_read_keys():
    """Every key the code calls get_setting()/get_vector_setting() on."""
    read = set()
    for source in sorted((pl.ROOT / "components").glob("*.cc")) + \
            sorted((pl.ROOT / "scheduler").glob("*.cc")) + \
            sorted((pl.ROOT / "utils").glob("*.cc")):
        text = source.read_text(encoding="utf-8")
        read |= set(re.findall(r'get_(?:vector_)?setting\(\s*"([a-z_0-9]+)"', text))
    return read


def undeclared_keys():
    """Keys the code reads minus keys any config declares."""
    read = code_read_keys()
    declared = set()
    for config in sorted(pl.ACCELERATORS.glob("*.cfg")):
        for line in config.read_text(encoding="utf-8").splitlines():
            found = re.match(r"^([a-z_0-9]+)\s*=", line.split("#")[0].strip())
            if found:
                declared.add(found.group(1))
    return sorted(read - declared)


def energy_keys(text: str):
    return {(s, k) for s, k, _ in pl.settings(text)
            if "energy" in k and k not in ("energy_unit", "energy_reference")}


def is_cycle_knob(key: str) -> bool:
    return "cycle" in key or key.startswith("t_")


def perturbed(section: str, key: str, value: str, toggle: bool):
    """The perturbed value, or None when no numeric perturbation is defined."""
    if toggle:
        return "0" if value.strip() in ("1", "1.0") else "1"
    fields = value.split(":")
    try:
        numbers = [float(f) for f in fields]
    except ValueError:
        return None
    # Preserve integer-ness: an unsigned setting given "2.0" fails to parse and silently keeps its
    # default, which would read as a dead knob.
    integral = all(f.strip().lstrip("-").isdigit() for f in fields)

    def render(x):
        return str(int(x)) if integral else str(x)
    if all(n == 0.0 for n in numbers):
        return ":".join(render(1) for _ in numbers)
    return ":".join(render(n*2) for n in numbers)


def rewrite(section: str, key: str, new: str):
    def apply(text):
        out, current, done = [], None, False
        for raw in text.splitlines():
            body = raw.split("#")[0].strip()
            header = re.match(r"^\[([a-z_]+)\]$", body)
            if header:
                current = header.group(1)
            elif current == section and not done and re.match(re.escape(key) + r"\s*=", body):
                out.append(f"{key} = {new}")
                done = True
                continue
            out.append(raw)
        if not done:
            raise SystemExit(f"{section}.{key} not found for rewrite")
        return "\n".join(out) + "\n"
    return apply


def scale_energies(factor):
    def apply(text):
        out, current = [], None
        for raw in text.splitlines():
            body = raw.split("#")[0].strip()
            header = re.match(r"^\[([a-z_]+)\]$", body)
            if header:
                current = header.group(1)
                out.append(raw)
                continue
            found = re.match(r"^([a-z_0-9]+)\s*=\s*([0-9.:]+)\s*$", body)
            if found and "energy" in found.group(1):
                fields = found.group(2).split(":")
                out.append(f"{found.group(1)} = " +
                           ":".join(str(float(f)*factor) for f in fields))
                continue
            out.append(raw)
        return "\n".join(out) + "\n"
    return apply


def diff(baseline, other):
    """Fields that moved beyond print-rounding slack, and fields that appeared or vanished.

    A STRUCTURAL change is evidence of liveness just as much as a numeric one: raising one
    component's frequency breaks the single-clock-domain precondition, so the power block correctly
    disappears. Counting only numeric moves reported that knob as dead.
    """
    moved, structural = [], []
    for key, (value, unit, decimals) in baseline.items():
        if key not in other:
            structural.append((key, "field disappeared"))
            continue
        new = other[key][0]
        if abs(new - value) > pl.tolerance(decimals):
            moved.append((key, value, new, pl.classify(unit)))
    for key in other:
        if key not in baseline:
            structural.append((key, "field appeared"))
    return moved, structural


def main() -> int:
    failures = []
    base_text = (pl.ACCELERATORS / f"{BASE}.cfg").read_text(encoding="utf-8")
    first = pl.run(BASE, WORKLOAD, MAPPING)
    if first is None:
        raise SystemExit(f"the baseline fixture {BASE} does not run")
    baseline = pl.parse_report(first)

    # KN1 -- determinism
    again = pl.parse_report(pl.run(BASE, WORKLOAD, MAPPING))
    moved, structural = diff(baseline, again)
    if moved or structural:
        failures.append(f"KN1: the same config run twice differs in "
                        f"{len(moved) + len(structural)} field(s) "
                        f"(e.g. {(moved + structural)[0]}); the model is not deterministic, so "
                        f"every comparison below is unsound")

    # KN2 -- energy scale invariance
    name = pl.variant(BASE, "__knob_scale", scale_energies(SCALE))
    scaled = pl.run(name, WORKLOAD, MAPPING)
    if scaled is None:
        failures.append("KN2: the uniformly scaled config does not run")
    else:
        scaled_fields = pl.parse_report(scaled)
        energy_checked = time_checked = 0
        for key, (value, unit, decimals) in baseline.items():
            if key not in scaled_fields:
                continue
            new, kind = scaled_fields[key][0], pl.classify(unit)
            if kind in ("energy", "energy_derived"):
                energy_checked += 1
                if abs(new - value*SCALE) > pl.tolerance(decimals, SCALE):
                    failures.append(f"KN2 {key}: {value} -> {new}, expected {value*SCALE} "
                                    f"({SCALE}x every energy unit cost); an energy that does not "
                                    f"scale is computed from a hardcoded constant")
            elif kind == "time":
                time_checked += 1
                if abs(new - value) > pl.tolerance(decimals):
                    failures.append(f"KN2 {key}: {value} -> {new} when only ENERGY costs changed; "
                                    f"an energy cost is feeding the timing path")
        if energy_checked < 20 or time_checked < 20:
            failures.append(f"KN2: only {energy_checked} energy and {time_checked} time fields "
                            f"were compared; the parser is no longer seeing the report")
    pl.discard("__knob_scale")

    # KN3 -- zero collapse
    name = pl.variant(BASE, "__knob_zero", scale_energies(0.0))
    zeroed = pl.run(name, WORKLOAD, MAPPING)
    if zeroed is None:
        failures.append("KN3: the zero-energy config does not run")
    else:
        zero_fields = pl.parse_report(zeroed)
        for key, (value, unit, decimals) in baseline.items():
            if key not in zero_fields or pl.classify(unit) not in ("energy", "energy_derived"):
                continue
            if zero_fields[key][0] != 0.0:
                failures.append(f"KN3 {key}: {zero_fields[key][0]} with every energy unit cost at "
                                f"0; a non-zero energy here is a spurious ADDITIVE term, which a "
                                f"scaling test cannot see")
    pl.discard("__knob_zero")

    # KN4-KN7 -- per-knob sweep
    live, dead, invalid, skipped = [], [], [], []
    energy = energy_keys(base_text)
    critical = [k for k in baseline if k[1] == "Critical-path latency"]
    for section, key, value in pl.settings(base_text):
        knob = (section, key)
        new = perturbed(section, key, value, knob in TOGGLE)
        if new is None:
            skipped.append(knob)
            if key not in EXEMPT_NON_NUMERIC:
                failures.append(f"KN4 {section}.{key} = {value!r} is non-numeric and not declared "
                                f"in EXEMPT_NON_NUMERIC; a knob with no perturbation defined is "
                                f"unswept, which is indistinguishable from unchecked")
            continue
        name = pl.variant(BASE, "__knob_probe", rewrite(section, key, new))
        text = pl.run(name, WORKLOAD, MAPPING)
        if text is None:
            invalid.append(knob)
            if knob not in EXEMPT_INVALID:
                failures.append(f"KN4 {section}.{key}: {value} -> {new} produces a config that "
                                f"does not run, and no constraint is declared for it")
            pl.discard("__knob_probe")
            continue
        other = pl.parse_report(text)
        moved, structural = diff(baseline, other)
        if structural and knob in MUST_BE_INVARIANT:
            failures.append(f"KN4 {section}.{key} must not affect results but changed the report's "
                            f"STRUCTURE ({structural[:2]})")
        if not moved and not structural:
            dead.append(knob)
            if knob not in EXEMPT_DEAD and knob not in MUST_BE_INVARIANT:
                failures.append(f"KN4 {section}.{key}: {value} -> {new} changed NOTHING in the "
                                f"report. The knob is declared, parsed, and consumed by nothing -- "
                                f"or the fixture does not exercise it, in which case say so in "
                                f"EXEMPT_DEAD and name the gate that does")
        else:
            live.append(knob)
            if knob in MUST_BE_INVARIANT and moved:
                failures.append(f"KN4 {section}.{key} must not affect results "
                                f"({MUST_BE_INVARIANT[knob]}) but moved {len(moved)} field(s), "
                                f"e.g. {moved[0]}")
            # KN5 / KN6
            if knob in energy:
                for field, was, now, kind in moved:
                    if kind == "time":
                        failures.append(f"KN5 {section}.{key}: raising an ENERGY cost moved the "
                                        f"timing field {field} ({was} -> {now})")
                    if kind in ("energy", "energy_derived") and now < was:
                        failures.append(f"KN6 {section}.{key}: raising an energy unit cost LOWERED "
                                        f"{field} ({was} -> {now})")
            # KN7
            if is_cycle_knob(key):
                for field in critical:
                    if field in other and other[field][0] < baseline[field][0] - 1e-9:
                        failures.append(f"KN7 {section}.{key}: raising a cycle cost LOWERED the "
                                        f"critical-path latency ({baseline[field][0]} -> "
                                        f"{other[field][0]})")
        pl.discard("__knob_probe")

    # KN8 -- no stale exemptions
    for knob, reason in EXEMPT_DEAD.items():
        if knob in live:
            failures.append(f"KN8 {knob[0]}.{knob[1]} is declared dead ({reason[:60]}...) but is "
                            f"live now; remove the exemption")
        elif knob not in dead:
            failures.append(f"KN8 {knob[0]}.{knob[1]} is declared dead but was never swept "
                            f"(it is now non-numeric or invalid); update the tables")
    for knob in EXEMPT_INVALID:
        if knob not in invalid:
            failures.append(f"KN8 {knob[0]}.{knob[1]} is declared unrunnable when perturbed but "
                            f"now runs; move it out of EXEMPT_INVALID")
    for key in EXEMPT_NON_NUMERIC:
        if not any(k == key for _, k in skipped):
            failures.append(f"KN8 {key} is declared non-numeric but was not skipped; update the "
                            f"tables")

    # KN9 -- the complementary scan: keys the code reads that no config declares.
    undeclared = undeclared_keys()
    classified = set(UNDECLARED_ALIAS) | set(UNDECLARED_OPTIONAL) | set(UNDECLARED_UNCOVERED)
    for key in undeclared:
        if key not in classified:
            failures.append(f"KN9 '{key}' is read by the code but declared by NO config, and is "
                            f"not classified. A knob no config declares is invisible to the sweep "
                            f"above -- classify it as an alias, an optional alternate, or an "
                            f"UNCOVERED cost knob")
    for key in sorted(classified - set(undeclared)):
        failures.append(f"KN9 '{key}' is classified as undeclared but some config now declares it "
                        f"(or the code stopped reading it); remove it from the KN9 tables")

    # KN11 -- a declared value that cannot be used must be REJECTED, not silently dropped.
    for label, section, key, value, why in (
            ("fractional-unsigned", "dram", "num_banks", "8.5",
             "a fractional value for an integer setting"),
            ("vector-for-scalar", "dram", "transfer_cycle", "2.75:2.75:2.75",
             "a per-datatype list for a key the model reads as a scalar (the maeri.cfg bug)"),
            ("bool-out-of-range", "systolic_array", "exist_temporal_buffer", "2",
             "a boolean flag given something other than 0 or 1")):
        name = pl.variant(BASE, f"__kn11_{label}", rewrite(section, key, value))
        if pl.run(name, WORKLOAD, MAPPING) is not None:
            failures.append(f"KN11 {section}.{key} = {value!r} ({why}) was ACCEPTED. It cannot be "
                            f"parsed into the field's type, so the declared value is silently "
                            f"dropped and the default stands -- exactly the failure that made "
                            f"eyerissv2.cfg unrunnable and maeri.cfg's DRAM transfer free")
        pl.discard(name)

    # KN10 -- the schema must be complete against what the code reads, not against what the
    # configs happen to declare.
    schema = schema_energy_keys()
    prefixes = tuple(key for key in schema if key.endswith("_"))
    for key in sorted(code_read_keys()):
        if "energy" not in key or key in schema or key.startswith(prefixes):
            continue
        failures.append(f"KN10 '{key}' is read by the code but is not in the RE5 energy schema "
                        f"(utils/energy_units.cc). Any config declaring it would be REJECTED as a "
                        f"typo, so the knob can never be exercised -- add it to g_energy_keys[]")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"fixture {BASE}: {len(baseline)} reported numeric fields, "
          f"{len(pl.settings(base_text))} declared settings")
    print(f"  live knobs      : {len(live)}")
    print(f"  dead (declared) : {len(dead)}  {[f'{s}.{k}' for s, k in dead]}")
    print(f"  invalid when perturbed (declared): {len(invalid)}")
    print(f"  non-numeric (covered by A/B fixtures): {len(skipped)}")
    print("KN1 the same config run twice gives an identical report ok")
    print(f"KN2 {SCALE:.0f}x every energy unit cost -> every energy exactly {SCALE:.0f}x, "
          f"no cycle or latency moved ok")
    print("KN3 every energy unit cost at 0 -> every reported energy exactly 0 (no additive term) ok")
    print(f"KN4 every one of {len(live)} live knobs moves the report; every dead one is declared "
          f"with its reason ok")
    print("KN5 raising an energy cost never moves a timing field ok")
    print("KN6 raising an energy cost never lowers a reported energy ok")
    print("KN7 raising a cycle cost never lowers the critical-path latency ok")
    print("KN11 a fractional integer, a list for a scalar and an out-of-range bool are each "
          "rejected, not dropped ok")
    print(f"KN10 every energy key the code reads is in the RE5 schema ({len(schema)} declared) ok")
    print("KN8 no exemption has gone stale ok")
    print(f"KN9 {len(undeclared)} key(s) read but declared by no config, all classified: "
          f"{len(UNDECLARED_ALIAS)} alias, {len(UNDECLARED_OPTIONAL)} optional alternate, "
          f"{len(UNDECLARED_UNCOVERED)} UNCOVERED cost knob(s) ok")
    print("    uncovered backlog: " + ", ".join(sorted(UNDECLARED_UNCOVERED)))
    print("ALL KNOB CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
