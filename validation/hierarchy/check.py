#!/usr/bin/env python3
"""Memory-hierarchy conservation gate.

Every other gate checks one boundary at a time. This one checks relations BETWEEN levels, which is
where a missing traffic SOURCE hides: the Eyeriss GLB-band gap (PA9) was found only because the
measured chip moved more through its GLB than the model accounts for at all, and it took an external
silicon reference to see it. A conservation law would find that class internally.

WHAT IS AND IS NOT DERIVABLE. A byte-level law across all five boundaries is NOT available from the
current report, and saying so is part of the result. Transactions are counted in each boundary's own
link width, and at fine granularity one transaction carries far less than a link word -- the
MAC<->LB boundary reports 262,144 transactions on a 128-bit link for 262,144 bytes of operand, so
multiplying by the width overstates the bytes 16x. This is the same granularity limit FM11 records
for the format metadata stream and RE1 hit for the accumulator.

The OUTPUT request chain likewise has no simple monotonicity, and three hypotheses for it were
refuted by measuring all seven architectures rather than one:
    hypothesis                                          refuted by
    array->GLB out == GLB->chip out - 1                 eyerissv2 (32,760 vs 2,048)
    monotone non-increasing outward                      tpu (63 < 64)
    0 exactly when edge_accumulation is on               simba and fsd (0 with it off)
So the output chain is recorded here as an OPEN question, not asserted. What the counter means
differs across architectures for reasons the report does not state; `Output tile residency` now
reports one of the contributing facts.

Checks (asserted; non-zero exit on failure), over every shipped architecture:
  HC1  LOAD REQUEST MONOTONICITY: for input and weight, the request count is non-increasing as you
       move OUTWARD (PE -> array -> GLB -> multi-chip -> DRAM). A load propagates outward only on a
       miss, so an outer level can never see more load requests than the inner level that generated
       them. An outer count that exceeds an inner one means traffic entering the hierarchy from
       nowhere -- exactly the missing-source shape.
  HC2  DRAM TRAFFIC LOWER BOUND: for every datatype, DRAM traffic is at least the distinct tensor
       volume. No mapping can move less data than the workload contains.
  HC3  REFETCH INTEGRALITY: DRAM traffic is an exact integer multiple of the distinct volume. The
       mapping factors divide the tensor dimensions (validation/architectures SM2 pins that), so a
       tensor is fetched a whole number of times; a fractional multiple would mean a partial tensor
       crossed the boundary.
  HC4  The `Output tile residency` statement is present, so a reader can tell an array-resident
       output tile from one written back during the layer.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import perturb_lib as pl                                            # noqa: E402

WORKLOAD, MAPPING = "gemm_64x64x64", "ws"
ARCHITECTURES = ("gemmini", "tpu", "tpuv3", "simba", "maeri", "eyerissv2", "fsd")
K = B = C = 64
# Distinct tensor volumes in bytes (int8 operands and output throughout).
VOLUME = {"Input data": C*B, "Weight": K*C, "Output data": K*B}
# Request counts, from innermost to outermost.
CHAIN = (("PE result", "# of data request to PE array"),
         ("PE array result", "# of request to Global buffer"),
         ("Global buffer result", "# of request to Chip-level processor"),
         ("Multi chip result", "# of request to off-chip memory"))
LOADS = ("Input data", "Weight")


def requests(text, section, label):
    body = text.split(section, 1)[1].split(label, 1)[1]
    found = re.findall(r"\*\s*(Input data|Weight|Output data)\s*:\s*(\d+)", body)[:3]
    return {name: int(value) for name, value in found}


def dram_bytes(target, text):
    section = (pl.ACCELERATORS / f"{target}.cfg").read_text(encoding="utf-8").split("[dram]", 1)[1]
    bits = int(re.search(r"^bitwidth\s*=\s*(\d+)", section, re.M).group(1))
    body = text.split("Multi-chip <-> DRAM transactions", 1)[1]
    found = re.findall(r"\*\s*(Input data|Weight|Output data)\s*:\s*\d+/\d+/(\d+)", body)[:3]
    return {name: int(value)*bits//8 for name, value in found}


def main() -> int:
    failures = []
    table = []
    for target in ARCHITECTURES:
        text = pl.run(target, WORKLOAD, MAPPING)
        if text is None:
            failures.append(f"HC1 {target}: does not run")
            continue
        chain = [requests(text, section, label) for section, label in CHAIN]
        traffic = dram_bytes(target, text)

        # HC1
        for datatype in LOADS:
            counts = [level[datatype] for level in chain]
            for inner, outer in zip(counts, counts[1:]):
                if outer > inner:
                    failures.append(
                        f"HC1 {target} {datatype}: an outer level sees {outer} load requests where "
                        f"the inner level generated {inner}. A load only propagates outward on a "
                        f"miss, so traffic is entering the hierarchy from nowhere "
                        f"(chain {counts})")
        # HC2 / HC3
        for datatype, volume in VOLUME.items():
            moved = traffic[datatype]
            if moved < volume:
                failures.append(f"HC2 {target} {datatype}: {moved} B of DRAM traffic is below the "
                                f"{volume} B the workload contains")
            elif moved % volume != 0:
                failures.append(f"HC3 {target} {datatype}: {moved} B is {moved/volume:.3f}x the "
                                f"{volume} B volume -- not a whole number of fetches, so a partial "
                                f"tensor crossed the boundary")
        # HC4
        if "Output tile residency" not in text:
            failures.append(f"HC4 {target}: the report does not state the output tile's residency, "
                            f"so an array-resident tile cannot be told from one written back")
        table.append((target, [level[d] for level in chain for d in ("Input data",)],
                      {d: traffic[d]//VOLUME[d] for d in VOLUME}))

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"{'architecture':>13} {'input load requests, inner -> outer':>38}   DRAM refetch (in/wt/out)")
    for target, counts, refetch in table:
        print(f"{target:>13} {str(counts):>38}   "
              f"{refetch['Input data']}x/{refetch['Weight']}x/{refetch['Output data']}x")
    print("HC1 load requests are non-increasing outward in every architecture ok")
    print("HC2 DRAM traffic covers the distinct tensor volume in every architecture ok")
    print("HC3 DRAM traffic is a whole number of tensor fetches ok")
    print("HC4 the output tile's residency is stated ok")
    print("OPEN: the OUTPUT request chain has no law this gate can assert -- see the module "
          "docstring for the three hypotheses measurement refuted")
    print("ALL HIERARCHY CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
