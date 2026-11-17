# Item 5, step 1: which backend times the attention consumer's KV read? (2026-11-17)

Source: local integration branch `kv-integration` = `tpu-measure` f306df6 + `origin/mosaic` 46d117f (merge commit
ce9d9c9). Built in its own worktree; the item-4 binary and results in the `tpu-measure` worktree are untouched.
Integration checks before anything else: `run_validation.sh` and `run_stats_timeline_validation.sh` pass; with all
knobs off the six reference GEMM runs give the same critical-path cycles as the `tpu-measure` binary (analytic and
refsb); bf16 32x4096 and int8 128x4096 A / B / C (fused_output, analytic) equal the saved full batch.
Core delta the merge brings (vs f306df6): components/dram.cc (+108/-), dram.h (+24), kvcache.cc/.h (batch),
scheduler/stats.cc/.h (+11), unittest/validation_test.cc (+5), validation/tpu_v5e/hbm_stream_probe.cc; configs
gemmini_rtl_mem.cfg and two gemm networks.

## Answer: analytical, also when a MOSAIC DRAM yaml is selected

Code path (`kv-integration` ce9d9c9):

- `scheduler/npu.cc`, `npu_t::apply_kv_cache_read` (line 1947): with `kvcache->attention_enabled()` the K and V
  supplies are computed by the local lambda `read_stream` (line 1963) from the `[dram]` unit costs only:
  `max(accesses * dram->u_read_cycle[OUTPUT], beats * dram->u_transfer_cycle)` plus the open-page row-activation
  term when `row_buffer_bytes > 0`. It never calls `dram->mosaic_stream()` or `ramulator_ns()`, and no branch on
  `dram->mosaic_dram_config` exists in it. The cache append (`write_cycle`) is computed the same way.
- `components/dram.cc`, `dram_t::describe_timing_model` (line 383): with a yaml it reports "MOSAIC Ramulator2 (...)
  for tile streams; KV-cache/SFU operand streams analytical" -- the report says what the code does.
- `validation/tpu_v5e/compare_attention.py` line 66 removes `NPUSIM_MOSAIC_DRAM` from the environment and deletes
  its temporary runs: the existing attention comparison is analytical by construction and keeps no evidence.

Reproducer, evidence retained (`kv_backend_probe.py`; this directory: per run the cfg with its `[kvcache]`
section, network, mapping, both reports, stdout; `kv_backend_probe.csv`, `.log`). Each case runs twice on identical
inputs, without and with `NPUSIM_MOSAIC_DRAM=hbm2e_tpu_v5e_refsb.yaml`:

| B, L | KV bytes read | projection DRAM busy (cycles), unset -> yaml | K / V supply (cycles), unset -> yaml | attention step (cycles) |
|---|---|---|---|---|
| 1, 2048 | 8 388 608 | 249.5 -> 267.9 | 7680 / 7680 -> 7680 / 7680 | 15759.2 -> 15759.2 |
| 8, 2048 | 67 108 864 | 264.5 -> 285.3 | 61440 / 61440 -> 61440 / 61440 | 125625.4 -> 125625.4 |
| 32, 8192 | 1 073 741 824 | 320.8 -> 346.8 | 983040 / 983040 -> 983040 / 983040 | 2008327.0 -> 2008327.0 |

The yaml is in effect (the projection's tile streams change) and the KV supply does not move: KV bytes / KV supply
time = 819.2 GB/s in every run, i.e. `[dram] read_cycle = 1.875` cycles per 8192-bit line at 1500 MHz.

So no result of this code can be called a refsb KV result. Per the review, nothing further was run.

## Minimal change proposed (NOT implemented; for review)

Scope: only the attention consumer's K and V READ streams, only when a MOSAIC DRAM yaml is selected; default and
no-yaml behaviour bit-identical.

1. `dram_t`: one public function, e.g. `double mosaic_read_stream_cycles(size_t bytes)`, that times ONE contiguous
   sequential read of `bytes` on the yaml's device through the existing `ramulator_ns()` with the DMA settings the
   tile streams already use (`mosaic_request_bytes`, `mosaic_max_inflight`, `mosaic_path_latency_ns`); returns a
   negative value without a yaml.
2. `npu_t::apply_kv_cache_read`, attention branch: when that function returns a time, the stream's cycle is that
   time (it replaces the analytical device / link makespan and the row-activation term, which the device model
   covers); otherwise the present formula. K and V stay two sequential streams, the decoder half and all energy
   terms stay as they are; byte counts and link transactions are unchanged.
3. Long streams: a KV half is up to 2 GiB. Proposal: simulate exactly up to 64 MiB; beyond, time(64 MiB) + (bytes -
   64 MiB) x the steady slope between 32 and 64 MiB -- the rule of `kv_streams.py`, whose 64 MiB exact vs 16 / 32
   extrapolation differed by 0.003%. The report states per stream whether it was exact or extrapolated.
4. Reports: `describe_timing_model` and the attention section name the backend of the KV read streams; the cache
   append stays analytical and is said to be.

Expected traffic: unchanged (KV read bytes = B x L x 8 x 128 x 2 x 2; append as now).
Expected timing, derived independently of the change from MOSAIC's `dram_bw_probe` on the same yaml (one fresh
stream): 2 MiB 2.8175 us, 8 MiB 11.227 us, 16 MiB 22.375 us, 32 MiB 44.668 us -- against 2.56 / 10.24 / 20.48 /
40.96 us analytical. E.g. B 1, L 2048: K = V = 4 MiB each.
Tests proposed: (a) without a yaml the three cases above reproduce 7680 / 61440 / 983040 exactly; (b) with the yaml
the K supply of a case equals the probe's time for the same bytes (exact sizes) and the extrapolation rule beyond
64 MiB; (c) `kv_backend_probe.py`: K / V supply now differ with the yaml while bytes, MACs and the no-yaml run do
not; (d) run_validation, stats_timeline, knobs-off identity.

Not addressed by this change, and to be stated with any result: QK / AV at `attention_macs_per_cycle = 65536` and
softmax at 0.005 cycles per score are configured assumptions, not hardware-validated rates; the step's overlap
model (streaming schedule, buffer 2); the decode layer still adds four separately simulated GEMMs. The external
`kv_streams.py` supply must not be added on top of an internal KV supply.
