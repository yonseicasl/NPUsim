# NPUsim PyTorch frontend

`frontend/pytorch` captures a `torch.nn.Module`, lowers supported operations to
framework-neutral executable IR, and feeds that artifact to the existing NPUsim
timing/traffic/energy core. PyTorch is an export-time dependency only; compiling,
validating, and simulating an existing artifact do not import PyTorch.

```text
torch.nn.Module
  -> torch.export
  -> npusim.graph.v1 capture artifact
  -> strict lowering
  -> npusim.exec.v1 executable artifact
  -> C++ workload_graph_t
  -> transitional Nebula layer adapter
  -> scheduler + timing/traffic/energy models
```

## Quick start

Validate and compile the bundled capture fixture:

```bash
python3 -m frontend.pytorch.cli validate \
  frontend/fixtures/linear_relu_graph.json
python3 -m frontend.pytorch.cli compile \
  --graph frontend/fixtures/linear_relu_graph.json \
  --output /tmp/linear_relu.exec.json
python3 -m frontend.pytorch.cli validate-executable \
  /tmp/linear_relu.exec.json
```

Run it through the simulator from the repository root:

```bash
models/model run-ir \
  configs/accelerators/gemmini_sfu.cfg \
  /tmp/linear_relu.exec.json \
  frontend/fixtures/linear_relu.map \
  pytorch_linear_relu
```

The optional last argument is the safe result-name stem. Result files record the
capture graph hash, declared executable hash, frontend schema, model name, and the
fact that runtime datatypes come from the accelerator configuration.

Run frontend regression tests with:

```bash
unittest/run_pytorch_frontend_validation.sh   # lowering + loader/lifetime units
python3 validation/frontend/check.py          # acceptance gate (FE1-FE8)
```

The acceptance gate pins the whole-path claims: recompiling each checked-in
graph fixture reproduces its executable fixture bit-for-bit (FE1); the
linear_relu fixture -- which IS the RTL-validated gemm_64x64x64 -- produces
compute-schedule and DRAM traffic IDENTICAL to the nebula path, with the
critical-path delta equal to the SFU busy cycles exactly (FE2); the conv_relu
fixture's computation count, SFU scalar operations, and output-commit identity
all match its own geometry (FE3); missing-SFU / wrong-op_id runs fail with
their specific errors (FE4); the residual DAG pins and releases live tensors
while executing BatchNorm, add, pooling branches, and concat end to end (FE5);
the real-capture pool chain reproduces its max/avg scalar-op identities (FE6)
and pool SFU-contract refusals (FE7); and the lenet classifier -- a real
torch.export capture with a flatten between conv and FC -- runs end to end as
exactly 7 operations, with the flatten elided as an alias and the first Linear
reading the pool's storage GLB-resident through it (FE8).

## Exporting a model

Export requires a Python environment containing a compatible PyTorch version.
Create a factory such as `my_model.py`:

```python
import torch


def build():
    model = torch.nn.Linear(64, 64).eval()
    return model, (torch.randn(64, 64),), {}
```

Then run:

```bash
python3 -m frontend.pytorch.cli export \
  --factory my_model:build \
  --output /tmp/model.graph.json \
  --model-name model
python3 -m frontend.pytorch.cli compile \
  --graph /tmp/model.graph.json \
  --output /tmp/model.exec.json
```

## Contracts

`npusim.graph.v1` preserves tensor dependencies, concrete or symbolic shapes,
dtypes, logical bytes, tensor strides/storage offsets, parameter qualified names,
and positional/keyword literal arguments. The capture artifact does not contain
weight values.

`npusim.exec.v1` is static and executable. Every tensor dimension must be concrete,
every operation must be explicitly modeled, and geometry must agree with tensor
shapes. Unsupported operations fail during lowering; they are never silently
counted as zero-cycle layers. The C++ loader repeats schema, topology, dtype/byte,
and geometry validation before constructing timing layers.

Mapping sections may bind to executable operations with `op_id`:

```ini
[connected]
op_id = linear
MAC = ...
```

All sections must either declare `op_id` or use legacy ordinal binding. Mixed,
duplicate, missing, extra, or wrong-kind bindings fail before simulation.

## Current milestone scope

Supported executable lowering now includes:

- `aten.linear.default`, `aten.conv2d.default`, and non-transposed
  `aten.convolution.default`;
- single-consumer ReLU/Leaky ReLU fusion and last-axis Softmax;
- NCHW MaxPool2d and AvgPool2d (`ceil_mode=False`);
- equal-shape tensor add/multiply for residual paths (broadcasting is fail-fast);
- tensor concat on a static axis;
- inference BatchNorm with static channel vectors;
- element-preserving flatten/view/reshape of contiguous whole-storage tensors,
  elided during lowering as `alias_of` views of the producer's storage: they
  never become operations, cost zero cycles/bytes by construction, and the
  lifetime tracker keys residency on the storage tensor so a consumer reading
  through the view keeps the producer's buffer alive instead of forcing a
  fictitious DRAM round trip. Non-contiguous or offset views are rejected,
  because eliding those would hide a real copy.

Executable operations are consumed in validated topological order, so branches and
residual dependencies no longer need to form a linear tensor chain. A capacity-aware
lifetime tracker keeps live graph inputs and intermediate activations in the GLB,
reuses last-use storage, and materializes only graph outputs or tensors that do not
fit. The policy does not silently evict: an unretained tensor is explicitly read from
or written to DRAM by its consumer/producer model.

Pooling, elementwise, BatchNorm, concat, nonlinear activation, and Softmax timing use
the configured `[sfu]` datapath. Max/average pooling count valid boundary samples;
AvgPool additionally honors `count_include_pad`. BatchNorm is inference-only and is
costed as a per-element affine mul-add; training BatchNorm and LayerNorm remain
unsupported. Concat is a copy-only memory operation with zero arithmetic operations.
Every unsupported form fails during lowering or C++ validation.

Asymmetric Conv2d stride, dilation, transposed convolution, nonzero output padding,
pool `ceil_mode`, elementwise broadcasting, dynamic executable shapes, and tensor
values are not yet supported. Executable-IR runs remain timing-only; functional builds
fail fast. Nebula remains only as a temporary layer-allocation bridge, while operation
geometry, topology, residency, and non-MAC timing are owned by executable IR.
Removing that bridge, materialization semantics for non-contiguous views (only
whole-storage contiguous reshapes alias today), broader LLM
operators/normalization, and calibrated SFU profiles are later migration phases.
