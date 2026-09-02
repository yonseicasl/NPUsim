#!/usr/bin/env python3
"""Config perturbation and report parsing, shared by the knob-sweep gate.

Every check built on this follows the same shape: take a baseline config, change ONE thing (or
apply one uniform transform), re-run, and compare the two reports field by field. That makes a
class of defect checkable that per-axis hand identities cannot reach -- a knob that is parsed and
then consumed by nothing produces a bit-identical report, and no amount of checking the axes that
DO work will reveal it.
"""
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ACCELERATORS = ROOT / "configs" / "accelerators"
MAPPINGS = ROOT / "configs" / "mappings"

# Report units, and what a uniform energy scaling must do to each. "energy" scales with the energy
# unit costs; "time" must not move at all when only energy costs change.
ENERGY_UNITS = ("pJ", "normalized", "uncalibrated")
SCALES_WITH_ENERGY = ("mW", "pJ*s", "pJ*s^2")
TIME_UNITS = ("cycles", "ms", "MHz")


def parse_report(text: str) -> dict:
    """Every `label : number [unit]` line in a result file, keyed by (section, label).

    Keyed by section because the same label appears in several sections (`Access energy`,
    `static_energy`, ...) and a global dict would silently collapse them.
    """
    fields, section = {}, "(header)"
    for raw in text.splitlines():
        header = re.match(r"^=+\s*(.+?)\s*=+$", raw.strip())
        if header:
            section = header.group(1)
            continue
        line = raw.strip()
        found = re.match(r"^(?:\*\s*)?([A-Za-z][^:]*?)\s*:\s*(-?[\d.]+)\s*(\S*)", line)
        if not found:
            continue
        label, value, unit = found.group(1).strip(), found.group(2), found.group(3)
        if value in (".", "-"):
            continue
        key = (section, label)
        # A repeated label inside one section (MAX/MIN/AVG blocks reuse names) gets an index so
        # both copies are compared rather than one overwriting the other.
        index = 0
        while (key + (index,)) in fields:
            index += 1
        fields[key + (index,)] = (float(value), unit, len(value.split(".")[-1]) if "." in value else 0)
    return fields


def classify(unit: str) -> str:
    if unit in ENERGY_UNITS or unit.startswith("pJ*"):
        return "energy" if unit in ENERGY_UNITS else "energy_derived"
    if unit in SCALES_WITH_ENERGY:
        return "energy_derived"
    if unit in TIME_UNITS:
        return "time"
    return "other"


def tolerance(decimals: int, factor: float = 1.0) -> float:
    """Rounding slack for a value printed to `decimals` places and then multiplied by `factor`.

    Both sides of the comparison are rounded prints, so the base contributes factor*half-ulp and
    the scaled value another half-ulp. Deriving this rather than guessing a relative epsilon is
    what keeps a 6-decimal ED2P from reading as a scaling violation.
    """
    half_ulp = 0.5*10.0**(-decimals)
    return factor*half_ulp + half_ulp + 1e-9


def variant(base: str, name: str, edits) -> str:
    """Write configs/accelerators/<name>.cfg from <base>.cfg with `edits(text)` applied, and give
    it the base's mapping directory. Returns `name`."""
    text = (ACCELERATORS / f"{base}.cfg").read_text(encoding="utf-8")
    (ACCELERATORS / f"{name}.cfg").write_text(edits(text), encoding="utf-8")
    source, target = MAPPINGS / base, MAPPINGS / name
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(source, target)
    return name


def discard(name: str) -> None:
    (ACCELERATORS / f"{name}.cfg").unlink(missing_ok=True)
    shutil.rmtree(MAPPINGS / name, ignore_errors=True)
    shutil.rmtree(ROOT / "result" / name, ignore_errors=True)


def run(name: str, workload: str, mapping: str):
    """Run the simulator; return the layer report text, or None if the config does not run."""
    completed = subprocess.run([str(ROOT / "npusim.sh"), "run", name, workload, mapping],
                               cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                               text=True)
    layer = ROOT / "result" / name / workload / mapping / "layer_0.txt"
    if completed.returncode != 0 or not layer.exists():
        return None
    return layer.read_text(encoding="utf-8")


def settings(text: str):
    """Every (section, key, value) actually declared in a config, comments stripped."""
    out, section = [], None
    for raw in text.splitlines():
        line = raw.split("#")[0].strip()
        if not line:
            continue
        header = re.match(r"^\[([a-z_]+)\]$", line)
        if header:
            section = header.group(1)
            continue
        found = re.match(r"^([a-z_0-9]+)\s*=\s*(.+)$", line)
        if found and section:
            out.append((section, found.group(1), found.group(2).strip()))
    return out
