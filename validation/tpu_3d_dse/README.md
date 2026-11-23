# TPU x 3D-stacked DRAM: design-space exploration brief (first phase, 2026-11-23)

Branch `tpu-3d-dse` (on the reviewed KV revision `kv-integration` 2ab9084). The TPU-fidelity phase is closed
here: no further measurement-comparison cohorts, probes or fitted-window experiments. Its unresolved limits
(below) travel with every number of this exploration.

**Question of this phase.** How much extra compute or buffering can use the bandwidth a DRAM stacked on the
logic die offers, and where does the limit move between memory, compute and on-chip buffers? It reports
performance / resource trade-offs. It does not name a physical optimum: no area, power or cooling budget has been
given, and NPUsim's energy numbers are unpriced.

## What is being compared

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
  - **Target architecture (user, 2026-11-23): one DRAM die bonded directly on the xPU -- 1-high.** HBM stack /
    channel expansion is therefore out of scope and not prepared; HBM2e stays only as the measured-against
    baseline. Consequences still to be worked out, none assumed here: with a single DRAM layer the capacity is
    bounded by the area over the logic die (the example profile's 8 GiB in one layer is an organisation in a
    preset table, not a checked floorplan), and MOSAIC's thermal examples start at 2 stacked dies -- a 1-high
    entry does not exist yet.
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
Buffers: input, one weight strip and one output strip against the declared GLB -- fits; the weight strip (4 MiB)
equals the declared 4 MiB, so K > 4096 or more MXUs per strip needs a look at that buffer before it is run.
(Conservative arithmetic on declared sizes; the simulator does not enforce buffer capacity.)
Numerical preflight for the HBDRAM profile (`stream_diag_hbdram_example/`): the long-stream rule predicts 64 MiB
to -0.18% (acceptance 1%); not used in the pilot (all streams exact).
Cost: 6 - 40 s per point, 9 min for M 512 (front-end load), 1.3 GB; two at a time.

## Next sweep: generated, not run (`dse_manifest.py` -> `results/2026-11-23_next_sweep_manifest/`)

12 hardware points x the 3 workloads (36 simulations if run), every memory point an existing MOSAIC yaml with its
coupled quantities derived:
- pin-count axis, capacity constant at 8 GiB: `cloud_{0.25, 0.5, 1, 2, 4}TBps` = 4 .. 64 channels x 1024 pins
  (256 .. 4096 GB/s nominal) -- RESOURCE-INCREASING in bond pins;
- organisation axis: `hb_dram_16ch` vs `cloud_1TBps` (same DQ, capacity and rate; 128 vs 512 columns per row) --
  ISO-RESOURCE;
- TPU side: `mxu4` (as is), `mxu8` (2 x MXUs), `mxu8_glb2` (and 2 x GLB SRAM) -- RESOURCE-INCREASING; accelerator
  configs generated, workloads to be remapped per point (the mapping generator takes the MXU count).
Not generated, with the reason in `manifest.json`: 256 x 256 arrays, HBM expansion, stack depth, DVFS.
Before any of it runs: the long-stream check per new yaml, and the 8-MXU weight strip against the GLB.

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
