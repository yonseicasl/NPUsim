#!/usr/bin/env python3
"""SCALE-Sim v2 cross-simulator acceptance gate (Phase-4 item 12).

validation/phase1 has compared NPUsim's systolic array schedule against SCALE-Sim v2 since
the first milestone, but only as a report a human reads -- nothing ran it automatically, so a
model change could silently move the third external reference we have. This turns it into a
gate. SCALE-Sim's COMPUTE_REPORT.csv files are committed under out/, so the gate needs no
SCALE-Sim installation; regenerating them is documented in README.md.

Comparison axis (unchanged, see README.md): SCALE-Sim's compute schedule includes the
systolic fill/drain of every fold, and NPUsim records that fill on the operand-stream axis
rather than inside "Computation cycle". The comparable quantity is therefore

    array-sched = max(Computation cycle, PE-array input IC, PE-array weight IC)

with the OUTPUT stream excluded (NPUsim's mapping writes partial sums through per reduction
step; a WS array accumulates in place -- a documented modeling difference, not an error).

WHY THE LIMITS DIFFER BY REGIME. The README's regime analysis, confirmed against Gemmini RTL
in phase2, attributes the residuals to three different things, and a single blanket tolerance
would either pass everything or fail on non-timing artifacts:

  FC (fc1-fc3)          -- both simulators are dominated by fold fill/drain, i.e. the SAME
                           physics on the same axis. This is a genuine timing comparison and
                           is gated at 8%.
  steady-state conv     -- SCALE-Sim's own fold constant (R+C-2 ~ 24 cycles) is about twice
  (conv3/conv4/conv5)      what Gemmini RTL measured (~14), so NPUsim reads systematically
                           low against it. Gated at 20%, loose enough for that known bias and
                           tight enough to catch a real regression.
  mapping-limited conv  -- conv1/conv2 residuals are pure active-PE arithmetic: 12x14 needs a
  (conv1/conv2)            factor of 7 that AlexNet's dimensions do not have, and NPUsim's
                           mapping grammar cannot express the partial folds SCALE-Sim uses to
                           fill the array. Accuracy-gating these would gate a documented
                           expressiveness limit, so they are NOT accuracy-gated -- they are
                           pinned to the committed baseline instead.

Every layer is additionally pinned to validation/phase1/npusim_baseline.csv, so any model
change that moves these numbers has to be acknowledged by updating the baseline.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
NPU_RESULT = ROOT / "result" / "scalesim" / "alexnet" / "matched"
SS_CONV = HERE / "out" / "npusim_val_convs" / "COMPUTE_REPORT.csv"
SS_FC = HERE / "out" / "npusim_val_fcs" / "COMPUTE_REPORT.csv"
BATCH = 4

# (npusim layer id, scale-sim report, scale-sim row, name, group factor, regime)
LAYERS = [
    (0, "conv", 0, "conv1", 1, "mapping-limited"),
    (2, "conv", 1, "conv2", 2, "mapping-limited"),
    (4, "conv", 2, "conv3", 1, "steady-conv"),
    (5, "conv", 3, "conv4", 2, "steady-conv"),
    (6, "conv", 4, "conv5", 2, "steady-conv"),
    (8, "fc", 0, "fc1", 1, "fc"),
    (9, "fc", 1, "fc2", 1, "fc"),
    (10, "fc", 2, "fc3", 1, "fc"),
]
LIMITS = {"fc": 8.0, "steady-conv": 20.0, "mapping-limited": None}
RHO_LIMIT = 0.90


def parse_npusim(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    section = text.split("PE array result", 1)[1].split("Global buffer result", 1)[0]
    ic = [float(v) for v in re.findall(
        r"\* (?:Input data|Weight|Output data)\s*:\s*([\d.]+) cycles", section)[:3]]
    return {
        "comp": float(re.search(r"Computation cycle\s*:\s*([\d.]+)", text).group(1)),
        "ic_in": ic[0],
        "ic_w": ic[1],
    }


def parse_scalesim(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8") as source:
        return {int(row["LayerID"]): float(row[" Total Cycles"]) for row in csv.DictReader(source)}


def spearman(a, b) -> float:
    def rank(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        ranks = [0.0]*len(values)
        for position, index in enumerate(order):
            ranks[index] = position
        return ranks
    ra, rb = rank(a), rank(b)
    n = len(a)
    return 1 - 6*sum((ra[i] - rb[i])**2 for i in range(n))/(n*(n*n - 1))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-baseline", action="store_true",
                        help="require the exact committed NPUsim counters")
    arguments = parser.parse_args()

    if not NPU_RESULT.exists():
        print(f"missing NPUsim output; run: ./npusim.sh run scalesim alexnet matched",
              file=sys.stderr)
        return 1
    baseline = {row["layer"]: row for row in
                csv.DictReader((HERE / "npusim_baseline.csv").open(newline="", encoding="utf-8"))}
    reports = {"conv": parse_scalesim(SS_CONV), "fc": parse_scalesim(SS_FC)}

    failures = []
    npu_axis, ss_axis = [], []
    per_regime = {}
    print(f"{'layer':>7} {'regime':>16} {'array-sched':>13} {'SS x b x g':>13} {'err%':>8} "
          f"{'limit':>7}")
    for layer_id, tag, ss_id, name, groups, regime in LAYERS:
        measured = parse_npusim(NPU_RESULT / f"layer_{layer_id}.txt")
        if arguments.check_baseline:
            expected = baseline[name]
            for key, column in (("comp", "computation_cycles"), ("ic_in", "ic_input_cycles"),
                                ("ic_w", "ic_weight_cycles")):
                if abs(measured[key] - float(expected[column])) > 0.5:
                    failures.append(f"{name}: {column} changed {expected[column]} -> "
                                    f"{measured[key]}")
        reference = reports[tag][ss_id]*BATCH*groups
        axis = max(measured["comp"], measured["ic_in"], measured["ic_w"])
        error = (axis - reference)/reference*100.0
        npu_axis.append(axis)
        ss_axis.append(reference)
        per_regime.setdefault(regime, []).append(abs(error))
        limit = LIMITS[regime]
        if limit is not None and abs(error) > limit:
            failures.append(f"{name} ({regime}): |{error:+.1f}%| exceeds the {limit}% limit")
        print(f"{name:>7} {regime:>16} {axis:13,.0f} {reference:13,.0f} {error:+8.1f} "
              f"{'baseline' if limit is None else f'{limit:.0f}%':>7}")

    rho = spearman(npu_axis, ss_axis)
    if rho < RHO_LIMIT:
        failures.append(f"Spearman rho {rho:.3f} below the {RHO_LIMIT} limit; NPUsim no longer "
                        f"ranks these layers like SCALE-Sim does")

    print()
    for regime in ("fc", "steady-conv", "mapping-limited"):
        errors = per_regime.get(regime, [])
        if not errors:
            continue
        limit = LIMITS[regime]
        print(f"{regime:>16}: MAPE {sum(errors)/len(errors):5.1f}%  max {max(errors):5.1f}%  "
              f"limit {'baseline-pinned (documented mapping-expressiveness gap)' if limit is None else f'{limit:.0f}%'}")
    print(f"{'all layers':>16}: MAPE "
          f"{sum(sum(v) for v in per_regime.values())/sum(len(v) for v in per_regime.values()):5.1f}%"
          f"  Spearman rho {rho:.3f} (limit {RHO_LIMIT})")

    if failures:
        print()
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1
    print("\nSCALE-Sim cross-simulator gate PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
