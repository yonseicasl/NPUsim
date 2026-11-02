# NPUsim vs Cloud TPU v5e (measured)

The PaLM comparison (`validation/palm/`) leaves NPUsim 41-60% under the published TPU v4 decode latency,
and public data cannot say where the rest goes. This directory measures it on a real TPU: per-GEMM
device times to check `configs/accelerators/tpu_v5e.cfg`, and a decode layer broken down by named scope
to size what NPUsim does not model (attention, vector work, collectives).

| file | runs on | does |
|---|---|---|
| `measure_tpu.py` | TPU (JAX only, one file) | GEMM sweep, vector ops, collectives, Llama-3-8B decode layer at tp=1 and tp=8 -> CSVs |
| `compare_v5e.py` | here (needs `models/model`) | NPUsim `tpu_v5e` on the measured GEMM shapes; decode GEMM vs non-GEMM share |
| `test_parse_trace.py` | anywhere with JAX | checks the profiler-trace parser on a synthetic TPU trace |
| `kernel-metadata.json` | Kaggle CLI | private script kernel on `TpuV5E8` |

`tpu_v5e.cfg` is `tpu_v4.cfg` with 4 MXUs (one TensorCore), 1500 MHz (197 TFLOPS / (4 x 128^2 x 2)) and
HBM2e 819 GB/s. Nothing in it is fitted to TPU data yet.

## Measuring

Kaggle (free, 20 TPU h/week, TPU v5e-8, no v6e; refs/tpu/README.md section 4). The account needs phone
verification for TPUs and internet.

1. Web: New Notebook -> Settings -> Accelerator `TPU v5e-8`, Internet on. Upload `measure_tpu.py` (Add
   Input -> upload, or paste it into a cell) and run `!python measure_tpu.py` (quick check:
   `--quick`). Results land in `/kaggle/working/tpu_v5e/`; download them from the Output tab.
2. CLI: `pip install kaggle`, put the API token in `~/.kaggle/kaggle.json`, replace `KAGGLE_USERNAME` in
   `kernel-metadata.json`, then
   ```
   kaggle kernels push -p validation/tpu_v5e
   kaggle kernels status KAGGLE_USERNAME/npusim-tpu-v5e-measure
   kaggle kernels output KAGGLE_USERNAME/npusim-tpu-v5e-measure -p validation/tpu_v5e/results/<date>
   ```
   `machine_shape` is what selects the v5e-8 per a 2026-10 report; check the run's `env.json`
   (`device_kind`) in case the CLI ignored it.

TRC or a paid Cloud TPU VM (`v5litepod-1` for gemm/vector/decode tp=1, `v5litepod-8` for all):
`python3 measure_tpu.py --out results/<date>` on the VM (installs `jax[tpu]` if missing).

A full run takes about 10-20 minutes. Device times come from the profiler (`jax.profiler.trace`,
perfetto export, XLA module/op events on `TPU:0`); if the trace has no TPU process the CSVs keep wall
time only and `compare_v5e.py` says so. `--keep-traces` keeps the raw traces for TensorBoard.

## Comparing

```
LD_LIBRARY_PATH=<npusim libs> python3 validation/tpu_v5e/compare_v5e.py validation/tpu_v5e/results/<date>
```

~15 minutes on 16 threads. Prints measured vs NPUsim vs roofline per GEMM (MAPE split into decode-shaped
M <= 128 and the rest) and writes `compare_gemm.csv`; then per decode case the measured layer time, its
projection GEMMs, attention against its KV-read roofline, psum, and the non-GEMM share -- at tp=1 also
NPUsim's GEMM sum against the measured GEMM scopes. Commit the result CSVs under `results/<date>/`.

Not covered: int8 weights (PaLM's low-latency setup), v6e (256x256 MXU; Kaggle has none), multi-chip
NPUsim. XLA may fuse ops across named scopes (e.g. a norm into the next matmul); `decode_ops.json`
keeps the per-op times for checking the attribution.
