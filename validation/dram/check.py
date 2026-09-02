#!/usr/bin/env python3
"""DRAM detail-knob regression (L8) and the analytical contract's reported scope.

The open-page row-buffer, JEDEC tRC, bank-count and bus-turnaround knobs existed but no
config enabled them, so every regression only exercised the flat model and the knobs could
have been arbitrarily wrong without anything noticing. This A/Bs three configs that differ
from each other in nothing but those knobs:

    gemmini              knobs off (flat per-access)
    gemmini_dram_detail  256 B rows, tRC = tRAS+tRP = 10, 8 banks, tWTR = tRTW = 2
    gemmini_dram_serial  the same with num_banks = 1

The 256 B row is deliberately small so one DRAM tile spans several rows on a fast workload;
the fixture exercises the knob arithmetic, it is not a calibrated DDR3 part.

Checks (asserted; non-zero exit on failure):
  DR1  Row activations and turnaround change TIMING ONLY: the Multi-chip <-> DRAM transaction
       counts are identical in all three runs.
  DR2  ... and they never touch the validated metric: the compute-schedule latency is
       identical in all three runs.
  DR3  Ordering: knobs off < 8 banks < 1 bank on the DRAM link axis. Activations cost time,
       and bank parallelism hides some of it.
  DR4  HAND IDENTITY for the bank term. Per datatype, a stream of `tile` elements spans
       ceil(tile_bytes / 256) rows; one bank serializes all of them while 8 banks serialize
       ceil(rows/8). The difference between the two runs must therefore be exactly
         sum_type( streams(type) * (rows(type) - ceil(rows(type)/8)) * tRC )
       with streams(type) recovered from the reported transaction count. The turnaround is
       identical in both runs and cancels out of this difference.
  DR5  HAND IDENTITY for the activation + turnaround term: the 8-bank run's excess over the
       knobs-off run must be exactly
         sum_type( streams(type) * ceil(rows(type)/8) * tRC ) + one bus turnaround
       (the single output write-back is the only direction flip in this workload).
  DR6  The report states which analytical contract is in force and what it does NOT capture,
       and those strings actually track the knobs -- the flat run must not claim a row model,
       and only the multi-bank run may claim bank interleaving.
  DR9  THE FLAT ROW-MISS MODEL. dram_t has two activation-latency models -- JEDEC tRC = tRAS+tRP
       when both are declared, and a flat `row_miss_cycle` per activation otherwise -- and the
       report names which one it used. Every config either left the row buffer off or declared
       tRAS/tRP, so the flat branch was read from the config and never executed. gemmini_dram_flat
       omits tRAS/tRP and declares row_miss_cycle, and the check requires the report to say "flat
       row-miss cost" and the activation latency to be exactly linear in that cost.
  DR8  WRITE-TO-READ TURNAROUND, on a workload that actually flips the bus that way. Every
       check above runs on gemm_64x64x64, where the schedule issues all of its loads and then
       its stores -- so the bus flips read->write exactly once and NEVER write->read, leaving
       `t_wtr_cycle` unexercised by this suite and, as the knob sweep found, by every other one.
       gemm_512x512x512 tiles the same machine enough times that stores precede later loads.
       The identity is exact and linear: 96 write->read flips, all charged to the stream that
       does the following read (weight here), so the weight transfer cycle must grow by
       96 x t_wtr_cycle and NOTHING else may move.
  DR7  E1/E9 HAND IDENTITY for row-activation ENERGY, and that it is OBSERVABLE. The energy
       accumulates into the DRAM interconnection (link) axis, which the report did not print at
       all until E1 -- so row_miss_energy = 20 pJ changed internal state and showed the reader
       nothing. Per datatype the printed value must equal
         streams(type) * rows(type) * row_miss_energy   (+ link transactions * transfer_energy,
                                                         which is 0 in this config)
       Unlike the latency term, energy does NOT benefit from bank parallelism: every activation
       costs energy whether or not it overlapped another bank's, so the same identity must hold
       for BOTH the 8-bank and the 1-bank run. With the knobs off it must be exactly 0.
"""
import math
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORKLOAD = ("gemm_64x64x64", "ws")
# DR8: a workload with enough DRAM tiles that stores precede later loads, so the write->read
# turnaround direction is exercised at all. Found by validation/knobs (the dead-knob sweep).
TURNAROUND_WORKLOAD = ("gemm_512x512x512", "ws")
TURNAROUND_FLIPS = 96              # hand-counted write->read bus flips for that schedule
TURNAROUND_STREAM = "Weight"       # the flip is charged to the stream doing the following read
VARIANTS = ("gemmini", "gemmini_dram_detail", "gemmini_dram_serial")
TYPES = ("Input data", "Weight", "Output data")

ROW_BYTES = 256          # row_buffer_size = 0.25 KB
TRC = 10.0               # t_ras_cycle 7 + t_rp_cycle 3
ROW_MISS_ENERGY = 20.0   # row_miss_energy, pJ per activation
BANKS = 8                # gemmini_dram_detail
TURNAROUND = 2.0         # t_wtr_cycle = t_rtw_cycle = 2
ELEMENT_BYTES = 1        # int8
DRAM_LINK_BYTES = 8      # [dram] bitwidth = 64


def run(target: str) -> str:
    layer = ROOT / "result" / target / WORKLOAD[0] / WORKLOAD[1] / "layer_0.txt"
    subprocess.run([str(ROOT / "npusim.sh"), "run", target, WORKLOAD[0], WORKLOAD[1]],
                   cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if not layer.exists():
        raise SystemExit(f"missing simulator output: {layer}")
    return layer.read_text(encoding="utf-8")


def measure(text: str) -> dict:
    dram = text.split("DRAM result", 1)[1]
    table = dram.split("Multi-chip <-> DRAM", 1)[1]
    transactions = {}
    for name, serialized in re.findall(
            r"(Input data|Weight|Output data)\s*:\s*\d+/\d+/(\d+)", table):
        transactions.setdefault(name, int(serialized))
    timeline = text.split("Layer timeline", 1)[1].split("MAC result", 1)[0]
    axes = re.findall(r":\s*([\d.]+) /\s*([\d.]+) /\s*([\d.]+) cycles", timeline)
    # E1: the DRAM link axis carries both the off-chip link energy and the row-activation
    # energy. This line only exists since E1 made it observable.
    link_energy_block = dram.split("Interconnection energy (link + row activation)", 1)
    if len(link_energy_block) < 2:
        raise AssertionError("DRAM interconnection (link + row activation) energy is not "
                             "reported; row_miss_energy would be invisible (E1)")
    link_energy = {name: float(value) for name, value in re.findall(
        # The unit label follows the config's declaration (E7/RE7), so it must not be hardcoded.
        r"(Input data|Weight|Output data)\s*:\s*([\d.]+) (?:pJ|normalized|uncalibrated)",
        link_energy_block[1])}
    processor = text.split("======= Processor ======", 1)[1].split("====", 1)[0]
    tile = {}
    for name, value in re.findall(r"(Input data|Weight|Output data)\s*:\s*(\d+)", processor):
        tile.setdefault(name, int(value))
    return {
        "transactions": transactions,
        "dram_link": float(axes[0][1]),
        "schedule": float(re.search(r"Compute-schedule latency\s*:\s*([\d.]+)", timeline).group(1)),
        "model": re.search(r"DRAM timing model\s*:(.*)", dram).group(1).strip(),
        "limits": re.search(r"not modeled here\s*:(.*)", dram).group(1).strip(),
        "tile": tile,
        "link_energy": link_energy,
    }


def turnaround_transfer_cycles(target: str, t_wtr: float):
    """DR8: the per-stream DRAM transfer cycles of `target` with t_wtr_cycle set to `t_wtr`."""
    import shutil
    sys.path.insert(0, str(ROOT / "validation"))
    import perturb_lib as pl
    name = f"__dr8_wtr{t_wtr:g}"
    pl.variant(target, name,
               lambda text: re.sub(r"^t_wtr_cycle\s*=.*$", f"t_wtr_cycle = {t_wtr:g}", text,
                                   count=1, flags=re.M))
    shutil.rmtree(pl.MAPPINGS / name, ignore_errors=True)
    shutil.copytree(pl.MAPPINGS / target, pl.MAPPINGS / name)
    text = pl.run(name, *TURNAROUND_WORKLOAD)
    if text is None:
        pl.discard(name)
        return None
    fields = pl.parse_report(text)
    out = {key[1]: value for key, (value, unit, decimals) in fields.items()
           if key[0] == "Multi chips - Off-chip memory" and key[2] == 1}
    pl.discard(name)
    return out


def main() -> int:
    result = {name: measure(run(name)) for name in VARIANTS}
    off, detail, serial = (result[name] for name in VARIANTS)
    failures = []

    # DR1
    for t in TYPES:
        counts = {result[name]["transactions"].get(t) for name in VARIANTS}
        if len(counts) != 1:
            failures.append(f"DR1 {t}: DRAM transaction counts differ across the knob variants "
                            f"({ {n: result[n]['transactions'].get(t) for n in VARIANTS} }); row "
                            f"activations and turnaround are timing, not traffic")

    # DR2
    schedules = {result[name]["schedule"] for name in VARIANTS}
    if len(schedules) != 1:
        failures.append(f"DR2: compute-schedule latency differs across the knob variants "
                        f"({ {n: result[n]['schedule'] for n in VARIANTS} }); DRAM timing must "
                        f"not move the validated metric")

    # DR3
    if not (off["dram_link"] < detail["dram_link"] < serial["dram_link"]):
        failures.append(f"DR3: expected knobs-off ({off['dram_link']}) < 8 banks "
                        f"({detail['dram_link']}) < 1 bank ({serial['dram_link']})")

    # Recover per-datatype stream counts and row spans from the report.
    streams, rows = {}, {}
    for t in TYPES:
        elements = off["tile"][t]
        per_stream = math.ceil(elements*ELEMENT_BYTES/DRAM_LINK_BYTES)
        total = off["transactions"][t]
        if per_stream == 0 or total % per_stream != 0:
            failures.append(f"DR4 {t}: {total} transactions is not a whole number of "
                            f"{per_stream}-transaction streams; cannot recover the stream count")
            streams[t], rows[t] = 0, 0
            continue
        streams[t] = total//per_stream
        rows[t] = math.ceil(elements*ELEMENT_BYTES/ROW_BYTES)

    # DR4
    expected_bank_saving = sum(
        streams[t]*(rows[t] - math.ceil(rows[t]/BANKS))*TRC for t in TYPES)
    measured_bank_saving = serial["dram_link"] - detail["dram_link"]
    if abs(measured_bank_saving - expected_bank_saving) > 0.5:
        failures.append(f"DR4: bank parallelism saved {measured_bank_saving} cycles, hand "
                        f"calculation says {expected_bank_saving} "
                        f"(streams { {t: streams[t] for t in TYPES} }, "
                        f"rows { {t: rows[t] for t in TYPES} }, tRC {TRC}, banks {BANKS})")

    # DR5
    expected_excess = sum(streams[t]*math.ceil(rows[t]/BANKS)*TRC for t in TYPES) + TURNAROUND
    measured_excess = detail["dram_link"] - off["dram_link"]
    if abs(measured_excess - expected_excess) > 0.5:
        failures.append(f"DR5: enabling the knobs added {measured_excess} cycles, hand "
                        f"calculation says {expected_excess} (activations + one "
                        f"{TURNAROUND}-cycle bus turnaround)")

    # DR6
    if "no row-buffer model" not in off["model"]:
        failures.append(f"DR6: the knobs-off run reports {off['model']!r}, which does not say "
                        f"the row model is disabled")
    for name in ("gemmini_dram_detail", "gemmini_dram_serial"):
        if "open-page" not in result[name]["model"]:
            failures.append(f"DR6 {name}: reports {result[name]['model']!r}, expected an "
                            f"open-page contract")
        if "turnaround" not in result[name]["model"]:
            failures.append(f"DR6 {name}: bus turnaround is enabled but not reported")
    if "bank interleaving" not in detail["model"]:
        failures.append(f"DR6: the 8-bank run does not report bank interleaving "
                        f"({detail['model']!r})")
    if "bank interleaving" in serial["model"]:
        failures.append(f"DR6: the 1-bank run claims bank interleaving ({serial['model']!r})")
    for name in VARIANTS:
        if not result[name]["limits"]:
            failures.append(f"DR6 {name}: the report does not state what the analytical model "
                            f"leaves out")

    # DR7: row-activation energy identity, on BOTH bank counts (energy does not benefit from
    # bank parallelism) and zero with the knobs off.
    row_energy = {}
    for t in TYPES:
        row_energy[t] = streams[t]*rows[t]*ROW_MISS_ENERGY
        for name in ("gemmini_dram_detail", "gemmini_dram_serial"):
            measured_energy = result[name]["link_energy"][t]
            if abs(measured_energy - row_energy[t]) > 0.5:
                failures.append(
                    f"DR7 {name} {t}: reported link+row energy {measured_energy} != hand "
                    f"{row_energy[t]} ({streams[t]} streams x {rows[t]} rows x "
                    f"{ROW_MISS_ENERGY} pJ)")
        if off["link_energy"][t] != 0.0:
            failures.append(f"DR7 gemmini {t}: knobs off must cost no link/row energy, got "
                            f"{off['link_energy'][t]}")

    # DR9 -- the flat row-miss model, and its linearity.
    def flat_cycles(row_miss):
        import shutil
        sys.path.insert(0, str(ROOT / "validation"))
        import perturb_lib as pl
        name = f"__dr9_{row_miss:g}"
        pl.variant("gemmini_dram_flat", name,
                   lambda text: re.sub(r"^row_miss_cycle\s*=.*$", f"row_miss_cycle = {row_miss:g}",
                                       text, count=1, flags=re.M))
        shutil.rmtree(pl.MAPPINGS / name, ignore_errors=True)
        shutil.copytree(pl.MAPPINGS / "gemmini_dram_flat", pl.MAPPINGS / name)
        text = pl.run(name, *WORKLOAD)
        pl.discard(name)
        if text is None:
            return None, None
        body = text.split("Multi chips - Off-chip memory", 1)[1]
        cycles = [float(v) for v in re.findall(
            r"\*\s*(?:Input data|Weight|Output data)\s*:\s*([\d.]+) cycles", body)[:3]]
        model = re.search(r"DRAM timing model\s*:(.*)", text).group(1)
        return cycles, model

    zero_cycles, flat_model = flat_cycles(0)
    if zero_cycles is None:
        failures.append("DR9: the flat row-miss fixture does not run")
    else:
        if "flat row-miss cost" not in flat_model:
            failures.append(f"DR9: the report names the model {flat_model.strip()!r}; with tRAS/tRP "
                            f"absent it must select and NAME the flat row-miss branch")
        base_cycles, _ = flat_cycles(12)
        double_cycles, _ = flat_cycles(24)
        for i, stream in enumerate(("input", "weight", "output")):
            per_unit = (base_cycles[i] - zero_cycles[i])/12.0
            if per_unit <= 0.0:
                failures.append(f"DR9 {stream}: row_miss_cycle changes nothing, so the flat branch "
                                f"is not being taken")
                continue
            expected = zero_cycles[i] + 24.0*per_unit
            if abs(double_cycles[i] - expected) > 0.5:
                failures.append(f"DR9 {stream}: doubling row_miss_cycle gave {double_cycles[i]}, "
                                f"expected {expected} ({per_unit:.0f} activations on the busiest "
                                f"bank x the flat cost)")

    # DR8 -- the write->read turnaround direction, on a workload that produces those flips.
    wtr_off = turnaround_transfer_cycles("gemmini_dram_detail", 0)
    for t_wtr in (2, 8):
        wtr_on = turnaround_transfer_cycles("gemmini_dram_detail", t_wtr)
        if wtr_off is None or wtr_on is None:
            failures.append(f"DR8: the {TURNAROUND_WORKLOAD[0]} turnaround run does not execute; "
                            f"t_wtr_cycle then has no coverage anywhere")
            break
        expected = TURNAROUND_FLIPS*t_wtr
        got = wtr_on[TURNAROUND_STREAM] - wtr_off[TURNAROUND_STREAM]
        if abs(got - expected) > 0.5:
            failures.append(f"DR8 t_wtr_cycle = {t_wtr}: the {TURNAROUND_STREAM} DRAM transfer "
                            f"cycles grew by {got}, expected {TURNAROUND_FLIPS} write->read flips "
                            f"x {t_wtr} = {expected}")
        for stream, value in wtr_off.items():
            if stream == TURNAROUND_STREAM:
                continue
            if abs(wtr_on[stream] - value) > 0.5:
                failures.append(f"DR8 t_wtr_cycle = {t_wtr}: the {stream} stream also moved "
                                f"({value} -> {wtr_on[stream]}); a write->read flip is charged to the "
                                f"stream that performs the following read, and only that one")

    if failures:
        for f in failures:
            print(f"FAIL {f}")
        print(f"{len(failures)} check(s) FAILED")
        return 1

    print(f"DR1 DRAM transaction counts identical across the three knob settings ok")
    print(f"DR2 compute-schedule latency identical at {off['schedule']} ok")
    print(f"DR3 DRAM link axis {off['dram_link']} < {detail['dram_link']} < "
          f"{serial['dram_link']} ok")
    print(f"DR4 bank parallelism saving {measured_bank_saving} == hand {expected_bank_saving} ok")
    print(f"DR5 activation + turnaround excess {measured_excess} == hand {expected_excess} ok")
    print(f"DR6 reported contract tracks the knobs; limits stated ok")
    print(f"DR9 the flat row-miss branch is selected, named, and linear in row_miss_cycle "
          f"({(base_cycles[0] - zero_cycles[0])/12.0:.0f} activations on the busiest bank) ok")
    print(f"DR8 write->read turnaround exercised: {TURNAROUND_FLIPS} flips x t_wtr_cycle on the "
          f"{TURNAROUND_STREAM} stream, nothing else moved ok")
    print(f"DR7 row-activation energy observable and == hand "
          f"{ {t: row_energy[t] for t in TYPES} } pJ on both bank counts, 0 with knobs off ok")
    print("ALL DRAM CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
