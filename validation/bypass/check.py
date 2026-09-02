#!/usr/bin/env python3
"""GLB bypass contract regression (P1-B accounting + L6 boundary depth).

Crosses the two settings that were previously conflated -- whether a datatype is
GLB-bypassed, and whether the GLB is double-buffered -- over four configs that differ from
configs/accelerators/gemmini.cfg in nothing else:

    gemmini_bypass_db    weight GLB-bypassed, GLB double_buffer = 1
    gemmini_bypass_sb    weight GLB-bypassed, GLB double_buffer = 0
    gemmini_resident_db  weight staged in GLB, GLB double_buffer = 1
    gemmini_resident_sb  weight staged in GLB, GLB double_buffer = 0

Crossing them is what makes the two contracts separable. The GLB double_buffer flag also
changes ACCOUNTING (a double-buffered destination hides its fill access, GB3), so comparing
latencies alone cannot tell "the bypass removed the decoupling" from "the flag hid the fill".
The boundary DEPTH the timeline ran at is therefore reported by the simulator and asserted
here directly.

Checks (asserted; non-zero exit on failure):
  BP1  P1-B accounting: a GLB-bypassed datatype costs ZERO GLB access cycles (no fill
       write, no read); the same datatype staged in the GLB costs more than zero.
  BP2  P1-B fabric: the GLB <-> PE-array weight link transactions are IDENTICAL in all four
       runs. Bypass skips the SRAM, not the wires -- the tile still has to be delivered to
       the array. (This is the half the pre-fix model got backwards: it kept the fill and
       deleted the delivery.)
  BP3  P1-B off-chip: the NoP weight transactions are IDENTICAL in all four runs. The tile
       still crosses the package either way.
  BP4  L6 depth contract: with a bypassed datatype the multi-chip -> GLB boundary runs at
       depth 1 (no tile-level decoupling) EVEN WHEN the GLB says double_buffer = 1, because
       a bypassed stream has no buffer there to decouple with. And the flag genuinely does
       decouple when every datatype is staged (resident + db runs at depth 2), so BP4 is not
       passing by ignoring the flag.
  BP5  Accounting ordering: resident+sb (fill charged) > resident+db (fill hidden) > 0, and
       bypass removes ALL of it. This pins BP1 to the amount actually removed rather than
       just its sign.
  BP6  Removing SRAM work never lengthens the layer: critical path (bypass) <= critical path
       (resident) at the same flag. Equality is expected whenever the GLB is not the
       bottleneck stage, which the report's bottleneck line shows.
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKLOAD = ("gemm_64x64x64", "ws")
VARIANTS = ("gemmini_bypass_db", "gemmini_bypass_sb",
            "gemmini_resident_db", "gemmini_resident_sb")
BYPASSED = {"gemmini_bypass_db", "gemmini_bypass_sb"}
DOUBLE_BUFFERED = {"gemmini_bypass_db", "gemmini_resident_db"}


def run(target: str) -> str:
    layer = ROOT / "result" / target / WORKLOAD[0] / WORKLOAD[1] / "layer_0.txt"
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, WORKLOAD[0], WORKLOAD[1]],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if not layer.exists():
        raise SystemExit(f"missing simulator output: {layer}")
    return layer.read_text(encoding="utf-8")


def measure(text: str) -> dict:
    glb = text.split("Global buffer result", 1)[1]
    access = glb.split("Access cycle", 1)[1].split("Energy", 1)[0]
    weight_access = float(re.search(r"Weight\s*:\s*([\d.]+) cycles", access).group(1))
    link = glb.split("serialized", 1)[1].split("# of request", 1)[0]
    weight_link = int(re.search(r"Weight\s*:\s*\d+/\d+/(\d+)", link).group(1))
    nop = text.split("Multi chip result", 1)[1].split("NoP transactions", 1)[1]
    nop_weight = int(re.search(r"Weight\s*:\s*\d+/\d+/(\d+)", nop).group(1))
    timeline = text.split("Layer timeline", 1)[1].split("MAC result", 1)[0]
    depths = [int(v) for v in re.findall(r":\s*(\d+) tiles", timeline)]
    if len(depths) != 4:
        raise AssertionError(f"expected 4 boundary depths, found {len(depths)}")
    return {
        "weight_access": weight_access,
        "weight_link": weight_link,
        "nop_weight": nop_weight,
        "depths": depths,
        "critical": float(re.search(r"Critical-path latency\s*:\s*([\d.]+)", timeline).group(1)),
        "bottleneck": re.search(r"Bottleneck level\s*:\s*(.+)", timeline).group(1).strip(),
        "bypassed_reported": "weight" in re.search(
            r"GLB-bypassed \(direct stream\)\s*:(.*)", glb).group(1).split(),
    }


def main() -> int:
    result = {name: measure(run(name)) for name in VARIANTS}
    failures = []

    for name in VARIANTS:
        if result[name]["bypassed_reported"] != (name in BYPASSED):
            failures.append(f"{name}: report says weight bypassed="
                            f"{result[name]['bypassed_reported']}, config says "
                            f"{name in BYPASSED} (config not picked up?)")

    # BP1
    for name in VARIANTS:
        access = result[name]["weight_access"]
        if name in BYPASSED and access != 0.0:
            failures.append(f"BP1 {name}: GLB-bypassed weight still costs {access} GLB "
                            f"access cycles")
        if name not in BYPASSED and access <= 0.0:
            failures.append(f"BP1 {name}: GLB-resident weight costs no GLB access cycles")

    # BP2 / BP3
    links = {result[name]["weight_link"] for name in VARIANTS}
    if len(links) != 1:
        failures.append(f"BP2: GLB<->PE-array weight link transactions differ across the "
                        f"four runs ({ {n: result[n]['weight_link'] for n in VARIANTS} }); "
                        f"bypass must skip the SRAM, not the fabric")
    nops = {result[name]["nop_weight"] for name in VARIANTS}
    if len(nops) != 1:
        failures.append(f"BP3: NoP weight transactions differ across the four runs "
                        f"({ {n: result[n]['nop_weight'] for n in VARIANTS} }); the tile "
                        f"crosses the package either way")

    # BP4
    for name in VARIANTS:
        depth = result[name]["depths"][1]
        expected = 2 if (name not in BYPASSED and name in DOUBLE_BUFFERED) else 1
        if depth != expected:
            failures.append(f"BP4 {name}: multi-chip -> GLB depth {depth}, expected "
                            f"{expected} (bypassed={name in BYPASSED}, "
                            f"double_buffer={name in DOUBLE_BUFFERED})")

    # BP5
    resident_sb = result["gemmini_resident_sb"]["weight_access"]
    resident_db = result["gemmini_resident_db"]["weight_access"]
    if not (resident_sb > resident_db > 0.0):
        failures.append(f"BP5: expected resident+sb ({resident_sb}) > resident+db "
                        f"({resident_db}) > 0; a double-buffered destination hides its fill "
                        f"access, so the single-buffered run must charge strictly more")

    # BP6
    for flag, suffix in (("db", "db"), ("sb", "sb")):
        bypass_critical = result[f"gemmini_bypass_{suffix}"]["critical"]
        resident_critical = result[f"gemmini_resident_{suffix}"]["critical"]
        if bypass_critical > resident_critical + 0.5:
            failures.append(f"BP6 ({flag}): bypassing the GLB lengthened the layer "
                            f"({bypass_critical} > {resident_critical})")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"BP1 bypassed weight GLB access = 0, resident = "
          f"{resident_db}/{resident_sb} (db/sb) ok")
    print(f"BP2 GLB<->PE-array weight link transactions identical = {links.pop()} ok")
    print(f"BP3 NoP weight transactions identical = {nops.pop()} ok")
    print("BP4 multi-chip -> GLB depth: bypass 1/1 (flag ignored), resident 2/1 (flag "
          "honoured) ok")
    print(f"BP5 fill accounting {resident_sb} > {resident_db} > 0 ok")
    print(f"BP6 critical path bypass <= resident at both flags "
          f"(bottleneck {result['gemmini_bypass_db']['bottleneck']}) ok")
    print("ALL BYPASS CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
