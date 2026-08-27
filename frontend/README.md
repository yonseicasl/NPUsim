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
unittest/run_pytorch_frontend_validation.sh
```

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

Supported lowering is intentionally narrow:

- `aten.linear.default` -> `npusim.linear`;
- `aten.conv2d.default` and non-transposed `aten.convolution.default` ->
  `npusim.conv2d`;
- single-consumer ReLU/Leaky ReLU fused into Linear or Conv2d;
- last-axis Softmax -> terminal `npusim.softmax`.

The transitional C++ adapter currently accepts a linear tensor chain. It rejects
asymmetric Conv2d stride, dilation, transposed convolution, nonzero output padding,
and non-terminal Softmax. Nonlinear activation and Softmax require an accelerator
`[sfu]` section. Executable-IR runs are timing-only until a versioned tensor-value
artifact is implemented; functional builds fail fast instead of using random
Nebula tensors.

Nebula remains as a temporary layer-allocation adapter for differential parity.
Removing that final dependency, adding general DAG/view/materialization handling,
dynamic-shape bindings, tensor values, broader operator coverage, and calibrated
SFU profiles are later migration phases.
