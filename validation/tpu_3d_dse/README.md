# TPU x directly bonded DRAM: decode DSE plan (updated 2026-11-24)

> **CURRENT PLAN:** [Fixed-memory decode DSE](#current-execution-plan-fixed-memory-decode-dse-updated-2026-11-24).
> Fixed four-die, 2 x 2, one-tier DRAM-on-logic target; decode latency first.
> Per die: 256 channels, **256-bit DQ per channel**, 4 banks/channel, 6 GB, 8 TB/s theoretical peak,
> 14 x 11 mm; aggregate 24 GB / 32 TB/s; logic 33 x 26 mm.
> The current plan adds DeepSeek-V4-Flash and Kimi K3 kernel/block workloads.
> Earlier pilot results and manifests below are historical; neither the 36-point memory grid nor the
> 40-point logic grid is the current execution plan. No expanded sweep has started.

Branch `tpu-3d-dse` (on the reviewed KV revision `kv-integration` 2ab9084). The TPU-fidelity phase is closed
here: no further measurement-comparison cohorts, probes or fitted-window experiments. Its unresolved limits
(below) travel with every number of this exploration.

**Question of this phase.** How much extra compute or buffering can use the bandwidth a DRAM stacked on the
logic die offers, and where does the limit move between memory, compute and on-chip buffers? It reports
performance / resource trade-offs. It does not name a physical optimum: no area, power or cooling budget has been
given, and NPUsim's energy numbers are unpriced.

## Current execution plan: fixed-memory decode DSE (updated 2026-11-24)

This section is the current plan. It supersedes the earlier workload selection and run proposals below,
including the 36-point memory sweep and the generated 40-point logic manifest. Those artifacts remain
historical records, not an execution queue. Status: **planning; no expanded workload sweep has started**.

### 1. Objective and fixed hardware contract

Minimize batch-1 decode service latency on one TPU-style logic die with four directly bonded DRAM dies.
Throughput and small batches are secondary. Report kernel, block and explicitly defined multi-block-window
latency; full-model token latency requires a complete, capacity-feasible model and system definition.

| Fixed quantity | Per DRAM die | Four dies |
|---|---|---|
| Arrangement | one vertical tier | lateral 2 x 2 above one logic die |
| Channels / banks | 256 channels, 4 banks/channel | 1024 channels, 4096 banks |
| DQ width | **256 bits per channel** | same per-channel width |
| Capacity | 6 GB | 24 GB |
| Theoretical peak bandwidth | 8 TB/s | 32 TB/s |
| Footprint | 14 x 11 mm | 28 x 22 mm envelope before gaps |

Logic footprint: 33 x 26 mm. GB and TB use decimal arithmetic. The latest user statement "DQ is 256bit"
is interpreted in the preceding per-channel interface context; it is not a 256-bit bus shared by the whole die.
The implied data rate is 8e12 * 8 / (256 * 256) = **976.5625 Mb/s per DQ**, or 0.9765625 Gb/s.
This is derived from the stated peak, not a measured signaling rate or a clock specification.
At 1 Gb/s per DQ the peak would be 8.192 TB/s/die, so that rounded rate must not silently replace 8 TB/s.

The 512-DQ stand-ins below are historical infrastructure probes and no longer candidate target interfaces.
DRAM timings, burst/address organization, refresh, read/write sharing, bonding overhead and actual operating
clock remain unresolved. Geometry fitting does not establish routing, area-budget or thermal feasibility.

For early screening, propose an explicitly idealized supply model capped at exactly 32e12 B/s aggregate;
later use four local domains capped at 8e12 B/s each, with finite cross-domain routes. Sustained efficiency
and timing uncertainty are sensitivity scenarios, not memory design axes. Exact 24 GB / 6 GB capacity
constraints are required even before a detailed DRAM organization is available.

### 2. Architecture-based workload suite

Keep the existing dense projection / conventional KV-attention microkernels as controls.
Add **DeepSeek-V4-Flash** and **Kimi K3** at their published dimensions; V4-Pro is a later shape sensitivity.
V4.1 and K3.1 are separate versions and are not automatic replacements. This is selection of architecture
workloads, not a claim that the full original models fit this chip.

| Workload family | Kernel scope | Design question |
|---|---|---|
| V4 local attention | shared K=V cache, local-window attention, partial RoPE, output projections | local-cache layout and matrix/vector feed |
| V4 CSA | compressor update, indexer, top-k, gather, local + compressed attention | sparse access, selection throughput, channel/bank utilization |
| V4 HCA | compressor update, attention over heavily compressed pool + local branch | compression overhead versus reduced traffic |
| V4 MoE + mHC | hash/learned routing, selected and shared experts, SwiGLU, combine, residual-stream mixing | expert placement, small matrix work, vector normalization, live buffers |
| K3 KDA | single-token recurrent update, projections, short convolution, gates, state read/write | persistent-state SRAM capacity, ports and spill cost |
| K3 gated MLA | latent-KV update/read, attention, projections and output gate | context-sensitive cache traffic and projection cost |
| K3 LatentMoE + AttnRes | latent down/up projection, routing, experts, SiTU-GLU, combination and depth-history access | expert parallelism, cross-domain traffic and activation storage |

Configuration anchors (official sources below; freeze exact revisions before runnable workloads):

- V4-Flash: hidden 4096; 43 decoder layers; 64 attention heads, head dimension 512; one shared KV head;
  local window 128; compression ratios 4 / 128 for CSA / HCA; index top-k 512; 256 routed experts,
  top-6 plus one shared expert, expert intermediate 2048; first three hash-routing layers;
  mHC multiplier 4. The compression schedule has an extra entry associated with model structure:
  confirm decoder versus next-token-prediction indexing in source before building the graph.
- K3: hidden 7168; 93 layers (69 KDA and 24 gated MLA); KDA 96 heads of dimension 128, short-convolution
  kernel 4; MLA KV rank 512; latent expert input dimension 3584 and intermediate 3072;
  896 routed experts, top-16 plus two shared experts; first layer dense.
  AttnRes block-size configuration is 12; infer its exact sublayer/boundary meaning from implementation,
  not from a guessed twelve-decoder-layer window.

KDA must use its recurrent decode path, not dense full-context QK/AV. V4 attention must not be substituted
with V3 MLA. LatentMoE must preserve latent projections instead of becoming a conventional hidden-size FFN.
Shared K=V storage is counted once; actual reads depend on scheduling and reuse.

### 3. Evaluation levels and scenarios

1. **Kernels:** actual tensor shapes, operation counts, payload versus physical traffic, padding, scratch,
   synchronization and allowed overlap. A selected-expert kernel may be evaluated conditionally, with its
   residency and route stated.
2. **Complete decode blocks:** five primary attention variants: V4 local / CSA / HCA and K3 KDA / gated MLA.
   Include the matching normalization, projections, MoE, mHC or AttnRes and state updates.
   Record layer ID and hash/dense/routed exceptions. Missing operations cannot have zero cost.
3. **Consecutive-layer windows:** preserve the actual attention/routing schedule and cross-layer state.
   Choose boundaries using the implementation, especially AttnRes. Supply and account for prior KV,
   recurrent and depth-history states; do not assume free initial loads or independent block execution.

Primary batch = 1; contexts = 2048 / 8192 / 32768. Use batch 4 / 8 and context 131072 only on shortlisted,
capacity-feasible cases. A KDA core's bounded recurrent state does not make the entire K3 context-independent:
its MLA layers still depend on context length.

For V4, sample ordinary tokens and positions before/on/after relevant compressor emission boundaries,
including cache rollover. Report periodic worst-phase and steady-state latency separately; these are not
empirical p99 statistics. For K3, compare resident versus spilled recurrent/depth state and measure the
corresponding transfer costs. Prefill is outside the timed scope, but initialized state and its footprint are not.

### 4. Precision and capacity contract

Use both explicit BF16 architectural controls and model-native mixed precision where supported.
V4 uses expert FP4 with other FP8 parameters; K3 publishes MXFP4 weights and MXFP8 activation QAT.
The K3 checkpoint configuration excludes several modules from weight quantization. Resolve precision per
operator from the pinned configuration/inference implementation; a top-level BF16 dtype or training
precision is not proof of every serving operator's execution format.

Count packed weights, scale metadata, padding, unpack/dequantization, accumulator and KV/state formats,
conversion throughput and scratch storage. INT8 is an optional separate sensitivity, not a proxy for FP4.
Unsupported formats receive an analytical-bound or unavailable label, not an assumed fourfold speedup.

Every complete block/window must fit:
resident weights + all KV/recurrent/depth state + activations + workspace + reserved memory <= 24 GB,
and <= 6 GB in each assigned DRAM domain. Record the reserved amount explicitly.
Apply corresponding peak live-set constraints to SRAM as well.

The original full models do not fit 24 GB. Even one K3 BF16 block containing its full routed-expert bank
can exceed 24 GB. Active parameters are not resident capacity. A selected-expert microbenchmark must not
stand in for arbitrary routing with free inactive experts. Native low-bit blocks are evaluated only after
their full capacity check; otherwise report infeasible. Reduced-size architecture proxies must have separate
names. Offload or multi-device execution requires a separate system model and is outside this plan.

### 5. Logic design space and implementation prerequisites

Explore at comparable total resources before increasing them:

- MXU array shape, count and partitioning across four memory domains; a GEMV/vector-oriented alternative
  needs an actual execution/dataflow model.
- Total SRAM and its allocation to weights, activations, recurrent state and depth history; banking,
  ports and feed bandwidth. Hold total SRAM fixed in one comparison: per-MXU VMEM otherwise grows with count.
- Staging/NoC bandwidth, DMA request size/concurrency, gather/scatter, top-k and vector/SFU capacity.
- Low-bit unpack/scale/conversion hardware and intermediate representations.
- Weight/expert/KV placement, local versus striped access, return/reduction traffic and synchronization.
  Copies count against capacity. Synthetic uniform, skewed and local expert routes are sensitivity cases,
  not measured routing traces. Batch size is not automatically tokens per expert.

Keep the clock fixed initially. Report compute, SRAM, vector, interface and interconnect resources together.
Process and component costs are missing, so do not turn the 858 mm^2 footprint into an invented area budget.

| Capability | Planning status / required action |
|---|---|
| Existing projection and KV-step paths | reuse after shape/mapping, batch-1 padding and live-set checks |
| V4 CSA/HCA, index/gather, mHC | operator and dependency coverage audit required; no verified full-block path yet |
| K3 recurrent KDA, gated MLA, LatentMoE, AttnRes | operator/state coverage audit required; no verified full-block path yet |
| Native low-bit mixed precision | audit storage and execution separately; unsupported costs must be explicit |
| Four local memory domains / cross-domain contention | model extension required; current aggregate endpoint is idealized |
| Exact four-bank DRAM scheduling / target timing | unsupported or unknown; preserve metadata without claiming fidelity |

### 6. Staged work and completion gates

| Stage | Work | Completion gate / deliverable |
|---|---|---|
| A. Workload specification | pin sources; extract shapes, layer schedules, precision and state; build operator DAGs | workload manifest, capacity feasibility, support-gap table |
| B. Minimal model coverage | add only operators/interfaces needed for the selected cases, ordered by bottleneck relevance | functional/traffic/state invariants and explicit model assumptions; no broad TPU calibration |
| C. Kernel screening | supported representative kernels across the existing five logic configurations; limited HBM anchors | identify compute, feed, vector, state or locality bottlenecks; mapping and padding checked |
| D. Block and window DSE | five block variants and selected consecutive windows; four-domain placement; shortlist 8-12 feasible logic candidates | critical-path latency and resources per scenario, with missing coverage visibly excluded |
| E. Robustness | 3-5 candidates: context, precision, batch, routing/locality, memory efficiency and periodic state updates | conditional ranking / rank reversals and 2-3 architecture families |

The previous estimate of approximately 24 initial runs and 24-36 block evaluations is superseded by the
expanded suite. After stage A, emit a bounded manifest listing every workload, candidate and scenario;
do not multiply all axes into a Cartesian sweep. Use concurrency <= 2 for the first execution batch.
The old 36- and 40-point manifests remain unstarted. This document update itself launches no simulations.

Latency follows the operator dependency graph with explicit compute/load/communication/serial-vector
terms and overlap rules. Do not blindly sum kernels or take their maximum. No full-model ITL, queueing
or p99 claims without the required model/system/arrival definition.

Stop exploring an axis when it does not materially change latency or rankings. Report rank reversals as
conditional results, without fitting simulator constants to force a winner. Physical power/thermal filtering
comes after process, energy, cooling and placement inputs; the stale full-footprint thermal scenario is unsuitable.

### 7. Reproducibility and deliverables

The workload manifest will record model ID, repository revision, configuration/source hashes, layer IDs,
operator shapes and dependencies, token position, batch/context, routing scenario, precision by tensor,
initial/final state, memory placement, peak live sets, traffic and support classification.
Use separate labels: **simulated**, **analytical estimate/bound**, **unsupported**, **capacity-infeasible**.
A schema requirement is not a claim that these new manifests or operators already exist.

First deliverables: pinned workload specifications and support-gap/feasibility matrix.
Final deliverables: per-kernel/block/window latency breakdown, utilization and traffic, resource tradeoffs,
capacity checks and a conditional shortlist. Broader simulator fidelity work stays closed.

Official sources checked 2026-11-24 (links are currently mutable; revision pinning is stage A):

- [DeepSeek-V4-Flash model card](https://huggingface.co/deepseek-ai/DeepSeek-V4-Flash)
- [DeepSeek-V4-Flash configuration](https://huggingface.co/deepseek-ai/DeepSeek-V4-Flash/raw/main/config.json)
- [Transformers V4 architecture implementation documentation](https://huggingface.co/docs/transformers/main/en/model_doc/deepseek_v4)
- [Kimi K3 official repository](https://github.com/MoonshotAI/Kimi-K3)
- [Kimi K3 configuration](https://huggingface.co/moonshotai/Kimi-K3/raw/main/config.json)
- [Kimi K3 technical blog](https://www.kimi.com/en/blog/kimi-k3)

---

## Historical pilot: what was compared

- **Chip (fixed in the pilot):** NPUsim `configs/accelerators/tpu_v5e.cfg` -- four 128 x 128 MXUs at 1.5 GHz with
  its declared buffers and the mappings of `validation/palm/palm_decode.py`. The `[dram]` clock equals the
  compute clock; MOSAIC stream times enter in that one cycle domain (a DVFS or mixed-clock axis needs a verified
  conversion first and is not enabled). Window knobs at their default 0 -- the 590-cycle window cost was fitted on
  v5e and is not carried to another architecture.
- **Memory, two architectures that are kept apart:**
  - `hbm2e_refsb` -- HBM2e stacks beside the die: the model compared with TPU v5e hardware; the paired baseline.
  - `hbdram_example` -- hybrid-bonded DRAM on logic: the repository's EXISTING example profile
    (`configs/architecture/dram/thermal_sweep/hb3d_example.yaml`). An existing model scenario, not a validated
    device, and nothing about its technology is added here.
  - **Target architecture (user, 2026-11-23, confirmed 2026-11-24): one DRAM die bonded directly on the xPU
    logic die -- 1-high.** HBM stack / channel expansion and 2 / 4 / 8-high stacks are out of scope; HBM2e stays
    only as the measured-against baseline. Channels, banks and DQ are organisation WITHIN the one die -- a channel
    count is not a number of dies. Capacity is what the selected organisation gives and is never multiplied by a
    layer count. No die area is known here, so no preset (including the example's 8 GiB) is shown to fit above
    the logic footprint: every point carries `footprint_checked = no`.
  - **1-high thermal scenario: prepared, not run** (`thermal_1high/dram_models_1high.yaml`, for MOSAIC's
    `tests/thermal_dram_sweep.py`). MOSAIC's examples start at 2 dies and are not results for this structure.
    Geometry = the script's defaults written out; `dram_power` 1.33 W/core is an ASSUMPTION (the repository's
    linear-with-depth rule), to be replaced by a justified power model before any thermal feasibility is claimed;
    the script's base chip is MOSAIC's 16-core flat-MAC chip, not the NPUsim model used here.
- **Stack depth is not a performance axis yet.** MOSAIC's `tests/thermal_dram_sweep.py` varies 2 / 4 / 8 stacked
  dies in the thermal model only; all three share one performance yaml. No stack-depth performance effect can be
  claimed until organisation / capacity / channels are tied to the depth. That sweep also uses MOSAIC's 16-core
  flat-MAC TPU-style chip, a different engine from this four-MXU NPUsim model: the two are not interchangeable.

Derived memory facts (`dse_pilot.memory_facts`: organisation from the yaml / preset tables, transaction size and
tCK from the runtime):

| memory | organisation | total DQ | capacity | nominal interface | transaction, tCK | sequential read, 1 MiB / 32 MiB |
|---|---|---|---|---|---|---|
| `hbm2e_refsb` | 16 channels x 128 DQ x 3200 Mbps | 2048 | 16 GiB | 819.2 GB/s | 32 B, 0.625 ns | 742 / 752 GB/s |
| `hbdram_example` | 16 channels x 1024 DQ x 500 Mbps | 16384 | 8 GiB | 1024.0 GB/s | 128 B, 2.0 ns | 1004 / 968 GB/s |

The two differ in capacity (16 vs 8 GiB) and in bond / pin count (8 x): this is NOT an iso-resource comparison.

## Pilot: six simulations (`dse_pilot.py` -> `results/2026-11-23_pilot/pilot/`)

Workloads -- synthetic and layer level. They are not end-to-end Llama throughput and do not show that a model
fits either memory (capacity was checked for the pilot's own tensors only):
- `gemm_m32`, `gemm_m512`: bf16 GEMM M x 4096 x 4096, `fused_output` scope (input and weights read from DRAM
  once, the output matrix not written back), `program_setup_cycle` 0 -- identical on both memories;
- `attn_b8`: bf16 KV attention step, B 8, L 2048, 8 KV heads x 128, two-pass, streaming schedule; the metric is
  the attention step appended after a small dummy projection.

| workload | HBM2e (baseline) | HBDRAM example | speed-up | regime, baseline -> HBDRAM |
|---|---|---|---|---|
| `gemm_m32` (32 x 4096 x 4096) | 45.79 us | 33.91 us | x1.351 | DRAM-bound (DRAM busy 99.5% -> 99.3%; compute busy 12.4% -> 16.8%) |
| `gemm_m512` (512 x 4096 x 4096) | 89.04 us | 88.99 us | x1.001 | compute-bound (compute busy 98.4%; DRAM busy 57.1% -> 42.2%) |
| `attn_b8` (B 8, L 2048) | 91.08 us | 71.19 us | x1.279 | KV-read bound (K = V read 66939 -> 52017 cycles; QK / softmax / AV 1024 / 2621 / 1024) |

Readings, with their scope:
- The memory-heavy kernels follow the memory; the compute-heavy one does not move (x1.001): at four MXUs it is
  already limited by the arrays, and its DRAM sits idle 43% (baseline) to 58% (HBDRAM) of the time. That idle
  share is the room "extra compute" would have to use -- the subject of the next sweep, not a result yet.
- **The size of the memory-bound gain depends on a modeling choice.** `gemm_m32` gains x1.351, which is the ratio
  of the two models' 1 MiB stream rates (1004 / 742 GB/s): tile streams are charged their share of a >= 1 MiB
  stream on a fresh device. `attn_b8`'s K / V read are complete 32 MiB streams and give x1.287 (752 -> 968 GB/s);
  the nominal interface ratio is x1.25. So the same hardware pair reads x1.25 - x1.35 depending on the stream
  length assumed. This is an approximation band to carry, not something to tune.
- The baseline numbers reproduce the fidelity phase on this build (45.79 us is variant A of the balance run;
  136624.2 cycles is the KV gate).

Checks (24 of 24, `pilot.json`): binary, every loaded shared object (56 files incl. libnpusim, libnebula,
libramulator), yaml, config and scripts hashed before and after -- equal; each report names the selected yaml as
its DRAM timing model, the two memories give different DRAM time for every workload, and both differ from an
analytical-DRAM preflight run (the backend processes GEMM and KV reads); MACs = M x K x N, DRAM payload = input +
weights once and no output, mapping products = N, M, K; reported KV read / append bytes = geometry, K / V read
exact (32 MiB, no extrapolation); compute busy, payload and QK / softmax / AV identical across the two memories.
Buffers: the live tensors fit the declared VMEM share. CORRECTED 2026-11-24: the pilot's check compared a whole
4-MXU weight strip (4 MiB) with the `[separate]` sizes, which the accelerator config declares PER MXU; per MXU
the live set is the activation (M x K x 2 B), one 128-column weight tile (K x 128 x 2 B = 1 MiB at K 4096, of
4 MiB) and one tile of fp32 accumulators (M x 128 x 4 B) -- independent of the MXU count (`dse_pilot.live_set`).
The earlier "weight strip equals the declared 4 MiB" flag is withdrawn; the saved pilot record keeps the
conservative numbers it was run with. (Arithmetic on declared sizes; the simulator does not enforce capacity.)
Numerical preflight for the HBDRAM profile (`stream_diag_hbdram_example/`): the long-stream rule predicts 64 MiB
to -0.18% (acceptance 1%); not used in the pilot (all streams exact).
Cost: 6 - 40 s per point, 9 min for M 512 (front-end load), 1.3 GB; two at a time.

## Historical memory sweep: superseded, generated and not run (`dse_manifest.py` -> `results/2026-11-24_next_sweep_manifest/`)

(`results/2026-11-23_next_sweep_manifest/` is the first version, kept: points only, no run directories.) The
2026-11-24 version records the 1-high architecture and the compute baseline (`tpu_v5e` for every point; the
unmeasured `tpu_v7.cfg` is a separate configuration, not in this grid) and holds all 36 run directories with the
mapping REGENERATED per point (`mxu8`: 8 x 128 = 1024-column strips, CHIPS_Y 8) and a `feasibility.json` each:
mapping products = N, M, K; whole strips; per-MXU VMEM share; the workload's tensors within the DRAM capacity; KV
stream within the exact range. All 36 pass. One 8-MXU point was simulated in a scratch copy to see that the
generated inputs run (memory-bound GEMM: latency as with 4 MXUs, compute busy 12.4% -> 6.4%); that is a smoke
test of the driver, not a sweep result. The grid itself is not started.

12 hardware points x the 3 workloads (36 simulations if run), every memory point an existing MOSAIC yaml with its
coupled quantities derived:
- pin-count axis, capacity constant at 8 GiB: `cloud_{0.25, 0.5, 1, 2, 4}TBps` = 4 .. 64 channels x 1024 pins
  (256 .. 4096 GB/s nominal) -- RESOURCE-INCREASING in bond pins;
- organisation axis: `hb_dram_16ch` vs `cloud_1TBps` (same DQ, capacity and rate; 128 vs 512 columns per row) --
  ISO-RESOURCE;
- TPU side: `mxu4` (as is), `mxu8` (2 x MXUs), `mxu8_glb2` (and 2 x GLB SRAM) -- RESOURCE-INCREASING; accelerator
  configs generated, workloads to be remapped per point (the mapping generator takes the MXU count).
Not generated, with the reason in `manifest.json`: 256 x 256 arrays, HBM expansion, stack depth, DVFS.
Before any of it runs: the long-stream check per new yaml (only `hbm2e_refsb` and `hbdram_p16` are checked).

## Limits carried from the closed fidelity phase

- NPUsim was compared with hardware on HBM2e TPU v5e-1 only. Streamed GEMM kernels: a few % (variant A:
  -2.8% mean / 3.2% MAPE on the 28 non-calibration projection cases); compute-bound kernels read ~6% short;
  the balanced regime (stream ~ compute) is unresolved; memory-bound bf16 is over-predicted by the opt-in window
  candidate, which stays off. A faster memory moves kernels toward exactly those less-certain regimes.
- The attention step has NOT been compared with hardware since its KV read moved to the DRAM backend. Its K / V
  reads are two fresh independent streams whose totals are spread uniformly over the supply tiles; QK / AV and
  softmax rates are configured assumptions; cache append, energy and link / row counts are analytical.
- Streams are fresh, contiguous and uncontended; no competing traffic, no real address layout.
- M > 512 GEMMs rely on an unvalidated row-tiling estimate; none is used here.
- No power, thermal or area result: NPUsim energy is unpriced and MOSAIC's example thermal powers are
  assumptions. Those are a later filter on a shortlist, with explicit cooling, geometry and limits.
- On a final shortlist only: repeat with the alternate documented supply assumptions (stream length, one vs two
  KV streams). A rank reversal is reported as inconclusive; the model is not tuned to produce a winner.

## TPU v7 (Ironwood) configuration (2026-11-23; no measurement behind it)

`configs/accelerators/tpu_v7.cfg` -- `tpu_v5e.cfg` re-parameterised from published figures; its header separates
what is published, what is derived and what is carried over. In short: clock 2200 MHz derived from the peak
(4 x 512^2 x 2 x f = 4614 TFLOPS fp8; bf16 at half that rate = 2307 TFLOPS), the chip modeled as 8 NPUsim arrays
of 256 x 256 (two per published "128K MAC/cycle" MXU), HBM3E 192 GiB at 7380 GB/s on the analytical DRAM model,
VMEM partitions x4, links scaled to the 256-wide array, and the two v5e-measured setup TIMES kept (an assumption;
set to 0 to drop it). Mappings need `palm.DIM = 256`, `palm.MXUS = 8`.

Parameterisation check (`check_tpu_v7.py` -> `results/2026-11-23_tpu_v7_config_check/`), bf16 GEMMs, isolated
scope, against the published rooflines -- a sanity check, not fidelity:

| GEMM | NPUsim | compute roofline | memory roofline | NPUsim / larger roofline | bottleneck |
|---|---|---|---|---|---|
| 512 x 4096 x 4096 | 8.46 us | 7.45 us | 5.68 us | x1.136 | compute |
| 128 x 4096 x 4096 | 5.24 us | 1.86 us | 4.83 us | x1.085 | DRAM |
| 32 x 8192 x 8192 | 18.80 us | 1.86 us | 18.33 us | x1.026 | DRAM |

MAC counts and DRAM payload are exact in all three. Open: no MOSAIC yaml of this HBM (the analytical model is
used), fp8 mode, the real MXU microarchitecture, and whether the v5e setup times apply.

## Fixed target and logic screening (2026-11-24)

**Fixed inputs** -- `target/target_profile.json`, units as the user stated them; not search axes:

| | per DRAM die | aggregate (4 dies, 2 x 2, one tier) |
|---|---|---|
| channels | 256 | 1024 |
| banks | 4 per channel = 1024 | 4096 |
| DQ width | 256 bits per channel | 256 bits per channel |
| capacity | 6 GB | 24 GB |
| bandwidth, theoretical maximum | 8 TB/s | 32 TB/s |
| footprint | 14 x 11 mm = 154 mm^2 | 616 mm^2, envelope 28 x 22 mm |

Logic die 33 x 26 mm = 858 mm^2. Objective: decode latency first (batch-1 token latency as the first screen), over
prefill and over throughput. V4-Flash and Kimi K3 architecture workloads are now selected in the current plan;
a deployed full model is not selected. DQ is now 256 bits per channel. Not given, and not invented: actual signaling clock/mode,
row / column organisation, timings and refresh, read / write sharing, die coordinates
and gaps, which logic attaches to which die, process node, area / power / cooling budgets.

**Derived** (`target/target_profile.py` -> `target/2026-11-24_consistency/`): per channel 187.5 Mb and 31.25 GB/s
(250 Gb/s) nominal; per bank 46.875 Mb; 0.312 Gb/mm^2; a whole die read at the nominal peak in 0.75 ms. The 2 x 2
envelope fits the logic footprint with 5 x 4 mm of slack before any gap or keep-out (71.8% of the logic area) --
a geometric statement only. The stated aggregates equal 4 x the per-die values. With the now-fixed 256-bit DQ, the stated 8 TB/s implies 976.5625 Mb/s per DQ.
The earlier 512-DQ / 500-Mb/s candidate is superseded. The original consistency reports are historical;
the updated generator reports the fixed width and derived rate.

**Can the performance model represent it? Partly.**

| aspect | representable now | note |
|---|---|---|
| 256 / 1024 channels | yes | runtime probe on stand-in yaml: runs; sequential read 7661 GB/s (256 ch x 512 DQ x 500 Mbps, 8.192 TB/s nominal) and 30897 GB/s (1024 ch) |
| 4 banks per channel | NO | every HBDRAM preset has one bank per channel |
| 6 GB per die | NO | 187.5 Mb per channel is not a power-of-two organisation, and no preset or tRFC / tREFI entry exists for it; NOT rounded |
| target timings | NO | unknown; the repository's 500 / 400 Mbps presets are not the target's |
| four interface domains | NO | NPUsim has one `[dram]` endpoint and one staging link |

These probes do not establish a faithful target-memory model. The historical 1024-channel stand-in has
512-bit DQ, a peak 2.4% above target, wrong capacity and no inter-domain locality/contending routes.
The current plan instead requires exact target capacity checks and proposes a separately labeled idealized
supply abstraction capped at 32e12 B/s before detailed four-domain timing is available.

**Decision variables (logic side):** MXU count and array shape; VMEM allocation, partition, banking, ports and
bandwidth; the on-chip feed from the DRAM interface to the MXUs; DMA request size and concurrency; mapping and
dataflow. Clock stays fixed until a power / thermal and clock-domain model exists. MXUs, SRAM, feed and interface
logic share one logic-area budget and no area model exists: every resource-changing point is exploratory.

**Historical screening manifest, generated and not run; superseded by the current plan** (`target/logic_manifest.py` -> `target/2026-11-24_logic_manifest/`):
5 logic points (MXUs 4 / 8 / 16, feed 1 / 4 / 32 KB per cycle) x 2 memory scenarios x 4 workload scenarios
(batch-1 matrix-vector 4096 x 4096 and 4096 x 14336, a batch-8 sensitivity, a batch-1 KV attention step) = 40
points, with accelerator configs and the first-order rates each configuration implies. Compute baseline
`tpu_v5e`; `tpu_v7.cfg` is separate and unmeasured.

What the configuration itself already says (arithmetic on declared rates, not a result): a batch-1 matrix-vector
product is all weight supply -- its compute bound is 0.17 us against 41 us of weight stream on HBM2e (a scratch
run of the existing model gives 41.8 us, DRAM-bound, compute busy 0.9%). With a 32 TB/s memory the limiter moves
on chip: `tpu_v5e.cfg`'s staging link carries 1 KB per cycle (1.54 TB/s) and each MXU takes one 128-element row
(256 B) per cycle, 1.54 TB/s for four. Widening the feed alone changes nothing (21.8 us); 8 / 16 MXUs give 10.9 /
5.5 us; reaching the memory's nominal peak would take about 21 x the weight-ingest rate of four MXUs. Whether a
real TPU's weight path is limited like the model's has never been tested -- the hardware comparison never
exceeded the HBM's 750 GB/s.

**Superseded:** the one-die description (2026-11-23); `results/2026-11-2{3,4}_next_sweep_manifest/` (memory
pin-count / channel grid: obsolete, not to be run); `thermal_1high/dram_models_1high.yaml` predates the four-die
footprint -- its memory die defaults to the full logic die and must be given the real placement before use.
