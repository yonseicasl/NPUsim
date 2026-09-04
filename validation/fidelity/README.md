# Fidelity measurement environment (evaluation.md Section B)

Measures how accurately NPUsim reproduces **execution cycles** (and DRAM traffic where the
reference has it) against RTL / silicon **ground truth**, for the configurations that have
committed golden reference data.

## Run

```bash
bash validation/fidelity/measure_fidelity.sh
```

Builds the model, runs the reference workloads, compares against golden, and prints a
consolidated table (also written to `fidelity_report.csv`).

## What is validated, and against what

| config | array | dataflow | reference (golden) | metric | result |
|---|---|---|---|---|---|
| **Gemmini** | 16×16 | weight-stationary | Chipyard `GemminiRocketConfig`, Verilator cycle-exact RTL, 6 GEMM points (`validation/phase2/`) | cycles | MAPE 4.40% / max 7.86% |
| **Eyeriss** | 12×14 | row-stationary | published silicon latency, AlexNet (`validation/phase3/`) | latency | MAPE 4.26% / max 6.39% |
| **NVDLA** | nv_small | — | official `dc_*` RTL traces, CRC-verified, cycle **holdout** + DRAM traffic (`validation/nvdla/`) | cycles | MAPE 0.55% / max 0.75% |

NVDLA additionally reproduces **DRAM traffic exactly** (input/weight/output bytes) on every
comparable trace — see `validation/nvdla/compare_full.py`.

The Gemmini number is in-sample (the fold-fill / setup constants were fit on these GEMMs);
the NVDLA number is a **holdout** (parameters fit only on `cal`-tagged traces, never on the
3 holdout traces) — the stronger out-of-sample evidence.

## What is NOT here, and why

- **Gemmini 128×128 (WS or OS)** — no golden. The committed Gemmini RTL is 16×16 WS only;
  there is no output-stationary Gemmini RTL trace at all. Validating a 128×128 or OS config
  requires new RTL runs (Chipyard/Verilator with `DIM=128` and the OS dataflow). The report
  leaves labelled slots so those rows can be filled once such traces exist. Until then,
  128×128/OS numbers would be NPUsim *predictions*, not fidelity.
- **Energy (Joules)** — NPUsim's energy is `UNCALIBRATED` (no pJ provenance declared) and no
  Joule-level golden exists for any of these designs, so an energy-fidelity figure would be
  meaningless. The evaluation.md Section B metrics are cycles / DRAM traffic / on-chip
  traffic — not energy — by the same reasoning.

## Extending (when 128×128 / OS RTL becomes available)

1. Produce golden cycles from the RTL run into a CSV (mirror `validation/phase2/golden_rtl_cycles.csv`).
2. Add an NPUsim accelerator config with `height=128 width=128` and `pe_stationary =
   weight_stationary` (or `output_stationary`).
3. Run it on the same GEMM shapes and add a comparison, following `validation/check_timing.py`.
