# Item 5, step 2: K / V read of the attention consumer on the MOSAIC DRAM backend -- small gate (2026-11-23)

Source: `kv-integration`, this commit (the gate ran on a8e39fc + `core.diff`, which this commit contains
unchanged; `gate/kv_backend_probe.json` has the sha256 of `models/model`, libramulator, the yaml, `tpu_v5e.cfg`
and the scripts before and after the runs: equal. CORRECTION 2026-11-23 (found by Codex): libnpusim and libnebula
were NOT hashed in that gate -- the resolver only matched `name => path` lines of `ldd`, and those two are listed
by absolute path; libdramsim3 is not linked at all. An earlier version of this note claimed they were. The
resolver is fixed for later runs (every loaded shared object, unresolved = error); nothing is asserted
retroactively about this gate's library hashes). Reviewed scope:
Codex, 2026-11-23, "GO WITH CHANGES". No TPU data is used; nothing here is a fidelity result. The full
attention / decode cohort is NOT run.

## What changed in the core (`core.diff`, 5 files)

- `components/dram.cc` `ramulator_ns`: the process-static memo is keyed by the device configuration too -- the
  yaml's path and its bytes, read at every call. Before, two configurations in one process shared timings
  (found by Codex). Earlier batches ran one configuration per process and are not affected.
- `dram_t::mosaic_read_stream(bytes, &cycles, &extrapolated)`: ONE complete sequential read on a fresh device,
  issued through the same helper and DMA knobs as the tile streams (`mosaic_request_bytes`,
  `mosaic_max_inflight`, `mosaic_path_latency_ns`) but NOT amortized like them: the stream pays its own start.
  Result in [dram]-clock cycles (ns x `[dram] frequency` / 1e3, the unit of `u_read_cycle`). Exact up to 32 MiB;
  beyond, T(32 MiB) + (bytes - 32 MiB) x (T(32) - T(16)) / 16 MiB, flagged. Returns false without a selected yaml
  (analytical cost stands); with one selected a non-finite / non-positive time or slope ends the run -- no
  fallback. 0 bytes cost 0.
- `npu_t::apply_kv_cache_read`, attention branch: when the backend answers, its time REPLACES the analytical
  makespan including the analytical row-activation time. K and V are two fresh, independent streams (no DRAM
  state carried between the passes -- a modeling choice of the two-pass schedule, not a bound in either
  direction). Decoder share added once, as before. Energy, link beats, row-activation counts and the cache
  append stay the analytical `[dram]` bookkeeping; they are not backend counters.
- Report (only with a backend; without one the layer report is byte-identical to before): cost-basis line,
  `K / V DRAM read` (decoder excluded, bytes per pass, exact / extrapolated per stream), a line saying the stream
  total is spread uniformly over the tiles, and `(analytical)` on the append. `describe_timing_model` updated.

## Checks

`unittest/run_mosaic_kv_stream_test.sh` (`stream_diag_final/`, `stream_diag_final.log`; `stream_diag/` is the
same on the build before a no-op warning fix):
- **two configurations in one process** (refsb and `hbm2e_v5e.yaml`, 1 MiB read): 1412.500 / 1465.625 ns in fresh
  processes and A, B, A again in one process. With the old memo key the test fails (B returns 1412.500; checked
  by reverting the key and rebuilding).
- **long-stream rule, held out** -- acceptance fixed before the run: |predicted - exact| / exact <= 1% at 64 MiB.
  Predicted from 16 / 32 MiB 89188.750 ns, exact 89284.375 ns: **-0.107%**. (This helper's requests, not
  `dram_bw_probe`'s; it gates a computational approximation, not TPU accuracy.)
- **two fresh half streams vs one stream of both halves**, same helper: half 1 MiB -0.20%, 2 MiB +1.45%, 4 MiB
  -0.09%, 8 MiB -0.01%, 16 MiB +0.14%. No consistent sign: the choice of two fresh streams is worth at most
  ~1.5% of the read time at these sizes and is not a bound.
- a stream of 4 MiB + 13 B is accepted (5580.625 ns, as 4 MiB to the printed precision); 0 bytes cost 0.

Request pattern (decision (a): reported, not calibrated away). The helper issues 64-transaction requests (2048 B
on this yaml), unlimited in flight, no path latency -- `tpu_v5e.cfg` sets no `mosaic_*` DMA knob. Against
`dram_bw_probe` on the same yaml (`2026-11-11_codex_audit/kv_refsb.csv.json`): 2 MiB 2830.625 vs 2817.5 ns
(+0.47%), 8 MiB 11171.250 vs 11226.875 (-0.50%), 16 MiB 22345.000 vs 22375.000 (-0.13%), 32 MiB 44626.250 vs
44668.125 (-0.09%), 64 MiB 89284.375 vs 89251.875 (+0.04%).

`kv_backend_probe.py` (`gate/`: per run cfg + `[kvcache]`, network, mapping, both reports, stdout;
`kv_backend_probe.csv` / `.json`; `gate.log`): **39 of 39 checks**. Cycles at 1.5 GHz:

| B, L, schedule | KV read / append bytes (report = geometry) | K = V supply, no yaml -> refsb | stream | attention step, no yaml -> refsb | bytes / supply time |
|---|---|---|---|---|---|
| 1, 2048, streaming | 8 388 608 / 4096 | 7680.0 -> 8370.9 | exact (4 MiB) | 15759.2 -> 17141.1 | 819.2 -> 751.6 GB/s |
| 8, 2048, streaming | 67 108 864 / 32 768 | 61440.0 -> 66939.4 | exact (32 MiB) | 125625.4 -> 136624.2 | 819.2 -> 751.9 GB/s |
| 32, 8192, streaming | 1 073 741 824 / 131 072 | 983040.0 -> 1069595.6 | extrapolated (512 MiB) | 2008327.0 -> 2181438.3 | 819.2 -> 752.9 GB/s |
| 1, 2048, blocking | 8 388 608 / 4096 | 7680.0 -> 8370.9 | exact | 15951.2 -> 17333.1 | 819.2 -> 751.6 GB/s |

- Without a yaml: K / V supply 7680 / 61440 / 983040, the attention steps and the projection DRAM busy reproduce
  `2026-11-17_kv_backend_trace`, and the three layer reports are byte-identical to those of the unpatched binary.
- With the yaml: K and V DRAM read equal the stream helper for the same bytes (8370.9375, 66939.375, 1069595.625
  cycles); supply = DRAM read + half the decoder (0 here: uncompressed), KV supply = K + V -- nothing added twice;
  QK / softmax / AV, append cycles and bytes equal the no-yaml run; projection DRAM busy as before the patch
  (267.9 / 285.3 / 346.8).
- Regression: `run_validation.sh`, `run_stats_timeline_validation.sh` pass; knobs off, the six reference GEMMs
  (analytic and refsb) give the same critical-path cycles as before; bf16 32 x 4096 and int8 128 x 4096 A / B / C
  equal the saved 2026-11-12 batch on both DRAM models.
- Without MOSAIC support in the build: unchanged -- a set `NPUSIM_MOSAIC_DRAM` is refused at start-up; with no
  yaml the analytical path runs (`dram.cc` compiles without `NPUSIM_MOSAIC_DRAM`, syntax-checked).

Next to the external one-stream supply of `kv_streams.py` (Codex's, unchanged; NOT to be added to the internal
one): K + V here 11.161 / 89.253 / 1426.128 us against 11.227 / 89.254 / 1426.842 us (-0.58%, -0.00%, -0.05%) --
two fresh half streams on this helper's requests vs one stream on `dram_bw_probe`'s.

## Scope limits (to be stated with any result that uses this)

- The stream time is a fresh-stream TOTAL per pass. The streaming schedule in `stats.cc` divides it uniformly
  over `kv_tiles`, which spreads the first-access cost over the tiles: these are not per-tile DRAM arrival
  times, and no end-to-end overlap of supply and compute is validated by a total. Unchanged by this patch;
  blocking vs streaming differ by 192 cycles in the one case run both ways, with and without the yaml.
- QK / AV at `attention_macs_per_cycle = 65536` and softmax at 0.005 cycles per score are configured
  assumptions (`compare_attention.py`), not hardware-validated rates.
- The cache append, all energies and the link / row-activation counts are analytical.
- `compare_attention.py` still removes `NPUSIM_MOSAIC_DRAM` and keeps no runs; it is not used here and the
  attention / decode cohort waits for the review of this gate.
