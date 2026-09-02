#!/usr/bin/env python3
"""PE-array temporal-buffer regression (found by the dead-knob sweep).

The PE array can hold an array-level temporal buffer between the GLB link and the PEs' local
buffers. It has five cost knobs -- `line_size`, `pe_array_read_cycle`, `pe_array_write_cycle`,
`pe_array_read_energy`, `pe_array_write_energy` -- and until this gate existed NONE of them had any
coverage in any suite. The cause was a single line: `exist_temporal_buffer` gates the whole path,
defaults to OFF, and all 48 shipped configs ship it COMMENTED OUT. So the buffer was pass-through
everywhere, its access energy printed 0.00 in every result file, and the five knobs could have been
wired to anything at all.

It was invisible to the dead-knob sweep too, which is worth stating: that sweep perturbs the
settings a config DECLARES, so a knob no config declares is outside its reach. What found this was
the complementary scan -- every key the code calls get_setting() on, minus every key any config
declares (validation/knobs KN9).

configs/accelerators/gemmini_all_knobs.cfg now turns the buffer on with distinct read (0.4) and
write (0.5) unit costs; gemmini.cfg remains the pass-through case.

Checks (asserted; non-zero exit on failure):
  AB1  The buffer is actually exercised: all three datatypes carry non-zero access energy. Without
       this the identities below would hold vacuously, which is exactly the state this gate ends.
  AB2  HAND IDENTITY, read and write separated by measurement rather than assumed. Setting one unit
       cost to 0 isolates the other, giving the access counts outright:
         read  : input 4096, weight 4096, output 0
         write : input 4096, weight 4096, output 65536
       and read + write must reconstruct the reported total exactly. The output READ count being 0
       is the contract, not an accident: with edge_accumulation the array never reads partial sums
       back out of the temporal buffer -- they stream to the edge accumulator -- but it does write
       them into it.
  AB3  The two directions are independent axes: doubling the read cost moves input and weight by
       `read accesses x cost` and leaves OUTPUT untouched (it has no read term); doubling the write
       cost moves all three.
  AB4  LINE-SIZE IDENTITY: doubling `line_size` halves every access count exactly -- the same bytes
       carried in half as many, twice as wide, lines. This is the only check that exercises
       `line_size` at the array level at all.
  AB5  PASS-THROUGH: with `exist_temporal_buffer` removed and the unit costs left non-zero, every
       temporal-buffer axis is exactly 0. This pins the default that made the knobs invisible, so a
       future change cannot start charging a buffer the config never asked for.
  AB6  RE8 CONTAINMENT: the buffer's energy reaches the PE-array component subtotal and the layer
       total by exactly the same delta, once each, and moves no other component.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import perturb_lib as pl                                            # noqa: E402

BASE, WORKLOAD, MAPPING = "gemmini_all_knobs", ("gemm_64x64x64", "ws"), None
READ_COST, WRITE_COST = 0.4, 0.5           # gemmini_all_knobs [systolic_array]
ARRAY_OUTPUT_TILE = 4096                   # array output tile in 8-bit line accesses
# CORRECTED 2026-08-25 (E20-3b): OUTPUT was 0 here, pinning a state that could not be physical --
# the array wrote output into its temporal buffer and never read it back out, yet the full 4096 B
# still arrived at DRAM. Data that leaves the array must be read out of the buffer it sat in. The
# array's finished output tile is 4096 elements at int8 over the 8-bit [spatial_arch] line, so the
# read-out is 4096 accesses, charged once because the psum is array-resident here (one output tile,
# no reduction loop with an output loop inside it).
READ_ACCESSES = {"Input data": 4096, "Weight": 4096, "Output data": 4096}
WRITE_ACCESSES = {"Input data": 4096, "Weight": 4096, "Output data": 65536}
STREAMS = ("Input data", "Weight", "Output data")
ROWS = ("MAC (compute+format)", "PE (local buffer)", "PE array", "Global buffer",
        "Multi-chip (NoP)", "DRAM")


def set_key(key, value):
    return lambda text: re.sub(rf"^{re.escape(key)} = .*$", f"{key} = {value}", text,
                               count=1, flags=re.M)


def drop_key(key):
    return lambda text: re.sub(rf"^{re.escape(key)} = .*$", f"#{key} removed", text,
                               count=1, flags=re.M)


def compose(*edits):
    def apply(text):
        for edit in edits:
            text = edit(text)
        return text
    return apply


def measure(label, edits):
    """Temporal-buffer access energy per stream, plus the component rows and the layer total."""
    name = f"__ab_{label}"
    pl.variant(BASE, name, edits)
    text = pl.run(name, *WORKLOAD)
    if text is None:
        pl.discard(name)
        return None
    array = text.split("PE array result", 1)[1].split("Global buffer result", 1)[0]
    block = array.split("Access energy (temporal buffer)", 1)[1]
    energy = {stream: float(value) for stream, value in
              re.findall(r"(Input data|Weight|Output data)\s*:\s*([\d.]+) \S+", block)[:3]}
    summary = text.split("Energy summary", 1)[1].split("Power summary", 1)[0]
    rows = {row: float(re.search(rf"\* {re.escape(row)}\s+[\d.]+\s+[\d.]+\s+([\d.]+)",
                                 summary).group(1)) for row in ROWS}
    total = float(re.search(r"Layer total energy\s*:\s*([\d.]+)", summary).group(1))
    pl.discard(name)
    return {"energy": energy, "rows": rows, "total": total}


def main() -> int:
    failures = []
    identity = (lambda text: text)
    base = measure("base", identity)
    if base is None:
        raise SystemExit(f"the baseline fixture {BASE} does not run")
    read_only = measure("ronly", set_key("pe_array_write_energy", "0"))
    write_only = measure("wonly", set_key("pe_array_read_energy", "0"))
    unpriced = measure("unpriced", compose(set_key("pe_array_read_energy", "0"),
                                           set_key("pe_array_write_energy", "0")))
    doubled_read = measure("r2", set_key("pe_array_read_energy", str(READ_COST*2)))
    doubled_write = measure("w2", set_key("pe_array_write_energy", str(WRITE_COST*2)))
    wide = measure("wide", set_key("line_size", "16:16:16"))
    passthrough = measure("passthrough", drop_key("exist_temporal_buffer"))
    for label, state in (("read-only", read_only), ("write-only", write_only),
                         ("unpriced", unpriced), ("2x read", doubled_read),
                         ("2x write", doubled_write), ("wide lines", wide),
                         ("pass-through", passthrough)):
        if state is None:
            raise SystemExit(f"the {label} variant does not run")

    # AB1
    for stream in STREAMS:
        if base["energy"][stream] <= 0.0:
            failures.append(f"AB1 {stream}: the temporal buffer carries no energy "
                            f"({base['energy'][stream]}), so every identity below is vacuous -- "
                            f"the state this gate exists to end")

    # AB2
    for stream in STREAMS:
        expected_read = READ_ACCESSES[stream]*READ_COST
        expected_write = WRITE_ACCESSES[stream]*WRITE_COST
        if abs(read_only["energy"][stream] - expected_read) > 0.5:
            failures.append(f"AB2 {stream} read term: {read_only['energy'][stream]} != hand "
                            f"{expected_read} ({READ_ACCESSES[stream]} accesses x {READ_COST})")
        if abs(write_only["energy"][stream] - expected_write) > 0.5:
            failures.append(f"AB2 {stream} write term: {write_only['energy'][stream]} != hand "
                            f"{expected_write} ({WRITE_ACCESSES[stream]} accesses x {WRITE_COST})")
        if abs(read_only["energy"][stream] + write_only["energy"][stream]
               - base["energy"][stream]) > 0.5:
            failures.append(f"AB2 {stream}: the isolated read and write terms "
                            f"({read_only['energy'][stream]} + {write_only['energy'][stream]}) do "
                            f"not reconstruct the reported {base['energy'][stream]}")

    # AB3
    for stream in STREAMS:
        read_delta = doubled_read["energy"][stream] - base["energy"][stream]
        if abs(read_delta - READ_ACCESSES[stream]*READ_COST) > 0.5:
            failures.append(f"AB3 {stream}: doubling the read cost moved it by {read_delta}, "
                            f"expected {READ_ACCESSES[stream]*READ_COST}")
        write_delta = doubled_write["energy"][stream] - base["energy"][stream]
        if abs(write_delta - WRITE_ACCESSES[stream]*WRITE_COST) > 0.5:
            failures.append(f"AB3 {stream}: doubling the write cost moved it by {write_delta}, "
                            f"expected {WRITE_ACCESSES[stream]*WRITE_COST}")
    # E20-3b: OUTPUT's read term is the tile read-out, and it must be exactly ONE tile -- the
    # partial sums themselves are accumulated at the array edge and never read back mid-reduction,
    # so a read term larger than one tile would mean psum traffic that should not exist here.
    output_read_delta = doubled_read["energy"]["Output data"] - base["energy"]["Output data"]
    if abs(output_read_delta - ARRAY_OUTPUT_TILE*READ_COST) > 0.5:
        failures.append(f"AB3 Output data: the read term is {output_read_delta/READ_COST:.0f} "
                        f"accesses, expected exactly one tile read-out ({ARRAY_OUTPUT_TILE}); with "
                        f"edge accumulation the array reads the FINISHED tile out once and never "
                        f"reads partial sums back mid-reduction")

    # AB4
    for stream in STREAMS:
        expected = base["energy"][stream]/2.0
        if abs(wide["energy"][stream] - expected) > 0.5:
            failures.append(f"AB4 {stream}: doubling line_size gave {wide['energy'][stream]}, "
                            f"expected exactly half of {base['energy'][stream]} = {expected} "
                            f"(the same bytes in half as many, twice as wide, lines)")

    # AB5
    for stream in STREAMS:
        if passthrough["energy"][stream] != 0.0:
            failures.append(f"AB5 {stream}: with exist_temporal_buffer removed the buffer still "
                            f"charges {passthrough['energy'][stream]}; the pass-through default "
                            f"must cost exactly nothing")

    # AB6
    delta = sum(base["energy"][stream] for stream in STREAMS)
    row_delta = base["rows"]["PE array"] - unpriced["rows"]["PE array"]
    if abs(row_delta - delta) > 0.5:
        failures.append(f"AB6: the buffer's {delta} of energy moved the PE-array subtotal by "
                        f"{row_delta}; it must land in its own component exactly once")
    total_delta = base["total"] - unpriced["total"]
    if abs(total_delta - delta) > 0.5:
        failures.append(f"AB6: the layer total moved by {total_delta}, not {delta}")
    for row in ROWS:
        if row == "PE array":
            continue
        other = base["rows"][row] - unpriced["rows"][row]
        if abs(other) > 0.5:
            failures.append(f"AB6: pricing the array buffer also moved the {row} subtotal by "
                            f"{other}")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"{'variant':>16} " + " ".join(f"{s:>12}" for s in STREAMS))
    for label, state in (("baseline", base), ("read only", read_only), ("write only", write_only),
                         ("2x line_size", wide), ("pass-through", passthrough)):
        print(f"{label:>16} " + " ".join(f"{state['energy'][s]:>12.1f}" for s in STREAMS))
    print("AB1 the array temporal buffer is exercised on all three datatypes ok")
    print(f"AB2 read {READ_ACCESSES} x {READ_COST} and write {WRITE_ACCESSES} x {WRITE_COST}, "
          f"isolated by measurement, reconstruct the total ok")
    print(f"AB3 read and write are independent axes; OUTPUT reads exactly one tile out "
          f"({ARRAY_OUTPUT_TILE} accesses) and never reads partial sums back ok")
    print("AB4 doubling line_size halves every access count exactly ok")
    print("AB5 without exist_temporal_buffer the buffer costs exactly 0 (pass-through default) ok")
    print(f"AB6 the {delta:.1f} lands once in the PE-array subtotal and once in the layer total, "
          f"and nowhere else ok")
    print("ALL ARRAY BUFFER CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
