# PaLM decode on TPU v4 vs Pope et al. 2023

`python3 validation/palm/palm_decode.py` (needs `models/model` built) compares the decode latency of
PaLM 540B and 62B on TPU v4 with Pope et al., *Efficiently Scaling Transformer Inference*, MLSys 2023,
Tables 2-3 (`golden_pope2023_decode.csv`; refs/tpu/README.md). A run takes ~15 minutes.

## Model

- Chip: `configs/accelerators/tpu_v4.cfg` (`tpu_v4_int8w.cfg`: int8 weights). 8 MXUs of 128x128 WS as
  NPUsim chips, 1050 MHz (275 TFLOPS bf16), HBM2 1200 GB/s, VMEM shares as GLBs. No TPU golden calibrates
  it; the structural declarations come from the TPU's documented organization:
  - `[systolic_array] streams_pipelined = 1` (new knob, default off): the MXU's double-buffered weight
    registers let folds flow without draining, so neither the mesh route fill per transfer nor the
    analytical accumulation drain per weight tile (127 cycles on a 128x128 array) is charged.
  - Accumulator port 128 fp32 per cycle (`[separate] bitwidth = 4096`), DMA double-buffered
    (`[multi_chip] double_buffer = 1`).
- Partitioning as in the paper: 2D weight-stationary weights (E over X, F over YZ, X and YZ minimizing
  2B(E/X + F/YZ)), the parallel attention/FFN block fused into one input and one output matmul,
  multiquery attention sharded over batch. Per chip and step: the fused GEMMs and the unembedding on
  NPUsim (simulated at 1 and 2 weight strips per MXU and extrapolated: the time is linear in strips),
  the KV-cache read at HBM bandwidth (memory-bound), and ICI time from the paper's bandwidth formula at
  270 GB/s with bf16 activations. Serial = GEMM + KV + ICI; overlap = max(GEMM + KV, ICI).
- Not modeled: per-collective latency, vector work (layernorm, softmax, SwiGLU), int8 dequantization,
  XLA kernel efficiency, sampling and host.

## Result (2026-11-02)

ms per decoded token; GEMM, KV and ICI summed over the layers.

| case | X x YZ | GEMM | KV | ICI | serial | overlap | roofline | measured | serial err |
|---|---|---|---|---|---|---|---|---|---|
| 540B, 64 chips, B 64, int8 | 4 x 16 | 9.79 | 0.21 | 1.20 | 11.28 | 10.08 | 7.25 | 28.44 | -60.3% |
| 540B, 64 chips, B 512, bf16 | 4 x 16 | 38.06 | 1.68 | 9.62 | 49.66 | 40.04 | 33.09 | 93.75 | -47.0% |
| 62B, 16 chips, B 32, int8 | 2 x 8 | 4.11 | 0.23 | 0.31 | 4.78 | 4.47 | 3.48 | 11.41 | -58.1% |
| 62B, 8 chips, B 512, bf16 | 1 x 8 | 31.95 | 7.27 | 6.96 | 47.24 | 40.28 | 36.33 | 79.69 | -40.7% |

The per-chip GEMMs land at 1.2-1.4x the roofline (padding, pipeline fill). The measured times are
1.7-2.5x the serial model -- 104-145 us per layer unexplained at small batch, 374-507 us at batch 512
-- in line with the paper's MFU (8-37%). Public data does not attribute that time; fitting it with
per-layer constants on four points would be curve fitting, so it is not done. Before the structural
declarations above, NPUsim's per-fold drain and a 2048-bit accumulator port inflated the GEMMs ~2x and
made the batch-512 cases look within 2-10% -- by compensation, not by modeling.

Closing the gap needs measured component times (e.g. a TPU profile of one decode layer) or the
communication-latency and vector-unit terms with their own evidence.
