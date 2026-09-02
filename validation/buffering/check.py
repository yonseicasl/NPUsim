#!/usr/bin/env python3
"""PE local-buffer double/single buffering A/B regression (P1-A).

Runs the SAME workload and mapping on two configs that differ in exactly one line --
configs/accelerators/gemmini.cfg (`double_buffer = 1`) and
configs/accelerators/gemmini_single_buffer.cfg (`double_buffer = 0`) -- and asserts the
hand-computable cycle contract, so a regression in either branch of the PE-busy
combination is caught by an executable test instead of by inspection.

Checks (asserted; non-zero exit on failure):
  LB-A  Compute-schedule latency is IDENTICAL in both runs. Buffering decides what can
        be hidden, never what is computed, so the externally validated metric must not
        move. (It is also the guard that the two runs really are the same workload.)
  LB-B  double-buffered PE busy == max(compute, access, link, overlap, format).
  LB-C  single-buffered PE busy == compute + max(access, link, overlap) + format.
        This is the P1-A fix: the overlap axis (cycle_mac_lb) is ALREADY the pipelined
        makespan of the access and link work, so the three axes are three views of one
        transfer, not three additive workloads.
  LB-D  single-buffered PE busy is STRICTLY LARGER than double-buffered (the
        serialization must be visible at all), and STRICTLY SMALLER than the pre-fix
        sum-of-all-axes value (the double counting must be gone). LB-D is skipped only
        if the two candidate values coincide, i.e. the config has no LB traffic to
        double count in the first place.
  LB-E  Critical-path latency does not shrink when buffering is removed.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKLOAD = ("gemm_64x64x64", "ws")
RUNS = {"double": "gemmini", "single": "gemmini_single_buffer"}


def run(target: str) -> str:
    layer = ROOT / "result" / target / WORKLOAD[0] / WORKLOAD[1] / "layer_0.txt"
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, WORKLOAD[0], WORKLOAD[1]],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if not layer.exists():
        raise SystemExit(f"missing simulator output: {layer}")
    return layer.read_text(encoding="utf-8")


def timeline(text: str) -> dict:
    tl = text.split("Layer timeline")[1].split("MAC result")[0]
    axes = [[float(a), float(b), float(c)] for a, b, c in
            re.findall(r":\s*([\d.]+) /\s*([\d.]+) /\s*([\d.]+) cycles", tl)]
    busy = [float(v) for v in re.findall(r"\* [^:\n]+:\s*([\d.]+) cycles \(", tl)]
    return {
        "schedule": float(re.search(r"Compute-schedule latency\s*:\s*([\d.]+)", tl).group(1)),
        "critical": float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", tl).group(1)),
        "format": float(re.search(r"PE format-IP axis\s*:\s*([\d.]+)", tl).group(1)),
        "mode": re.search(r"PE local buffer\s*:\s*(\S+)", tl).group(1),
        "pe_busy": busy[4],
        "pe_axes": axes[4],
    }


def main() -> int:
    result = {name: timeline(run(target)) for name, target in RUNS.items()}
    failures = []
    tol = 0.5

    for name, expected_mode in (("double", "double-buffered"), ("single", "single-buffered")):
        if result[name]["mode"] != expected_mode:
            failures.append(f"{name}: report says PE local buffer is {result[name]['mode']!r}, "
                            f"expected {expected_mode!r} (config not picked up?)")

    # LB-A
    if abs(result["double"]["schedule"] - result["single"]["schedule"]) > tol:
        failures.append(f"LB-A: compute schedule moved with buffering "
                        f"({result['double']['schedule']} -> {result['single']['schedule']})")

    # LB-B
    d = result["double"]
    expect_double = max(d["pe_axes"] + [d["schedule"], d["format"]])
    if abs(d["pe_busy"] - expect_double) > tol:
        failures.append(f"LB-B: double-buffered PE busy {d['pe_busy']} != "
                        f"max(compute, axes, format) {expect_double}")

    # LB-C
    s = result["single"]
    expect_single = s["schedule"] + max(s["pe_axes"]) + s["format"]
    if abs(s["pe_busy"] - expect_single) > tol:
        failures.append(f"LB-C: single-buffered PE busy {s['pe_busy']} != "
                        f"compute + max(axes) + format {expect_single}")

    # LB-D
    pre_fix_sum = s["schedule"] + sum(s["pe_axes"]) + s["format"]
    if s["pe_busy"] <= d["pe_busy"]:
        failures.append(f"LB-D: single-buffered PE busy {s['pe_busy']} is not larger than "
                        f"double-buffered {d['pe_busy']}; the serialization is invisible")
    if abs(pre_fix_sum - expect_single) <= tol:
        print("LB-D note: this config has no overlappable LB traffic, so the pre-fix "
              "sum-of-axes and the fixed makespan coincide; the double-counting half of "
              "LB-D is vacuous here")
    elif s["pe_busy"] >= pre_fix_sum - tol:
        failures.append(f"LB-D: single-buffered PE busy {s['pe_busy']} still equals the "
                        f"pre-fix sum-of-all-axes {pre_fix_sum}; the LB<->MAC transfer is "
                        f"charged more than once")

    # LB-E
    if result["single"]["critical"] + tol < result["double"]["critical"]:
        failures.append(f"LB-E: critical path shrank when buffering was removed "
                        f"({result['double']['critical']} -> {result['single']['critical']})")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"LB-A compute schedule identical: {d['schedule']} ok")
    print(f"LB-B double-buffered PE busy == max(compute, axes, format) = {d['pe_busy']} ok")
    print(f"LB-C single-buffered PE busy == compute + max(axes) + format = {s['pe_busy']} ok")
    print(f"LB-D single({s['pe_busy']}) > double({d['pe_busy']}) and < pre-fix sum({pre_fix_sum}) ok")
    print(f"LB-E critical path {d['critical']} -> {s['critical']} ok")
    print("ALL BUFFERING CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
