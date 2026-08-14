# NPUsim frontend

`frontend/pytorch` is the first NPUsim PyTorch frontend milestone. It captures a
`torch.nn.Module` with `torch.export`, then writes a versioned, framework-neutral
JSON graph artifact. The artifact contains tensor shapes, dtypes, logical byte
sizes, data dependencies, operation names, and provenance. It deliberately does
not contain model-weight values.

This separates model capture from timing simulation: an exported graph is stable,
reviewable, and can be validated without importing PyTorch.

## Quick start

Validate the bundled artifact (PyTorch is not required):

```bash
python3 -m frontend.pytorch.cli validate frontend/fixtures/minimal_graph.json
```

For export, use a Python environment containing a compatible PyTorch version.
Create a factory such as `my_model.py`:

```python
import torch


def build():
    model = torch.nn.Linear(4, 2).eval()
    return model, (torch.randn(1, 4),), {}
```

Then export it:

```bash
python3 -m frontend.pytorch.cli export \
  --factory my_model:build \
  --output artifacts/linear.graph.json \
  --model-name linear-demo
python3 -m frontend.pytorch.cli validate artifacts/linear.graph.json
```

A factory must return `(model, example_args)` or
`(model, example_args, example_kwargs)`, where `example_args` is a tuple.

## Graph IR contract

The current schema version is `npusim.graph.v1`. Every artifact records:

- producer version and model structural hash;
- tensors with shape, dtype, layout, kind, and exact logical bytes;
- topologically ordered operations and their tensor dependencies;
- a content hash (`graph_sha256`) when written by the exporter.

The validator rejects invalid dtype/byte combinations, undeclared references,
duplicate IDs, redefined tensors, and non-topological dependencies before a
simulation can consume the artifact.

## Current scope

This is P0-A from `issue/frontend/pytorch_vllm.md`: model capture and a validated
interchange contract. NPUsim's C++ mapping/timing path does **not yet consume**
this graph IR, and it does not yet provide vLLM request-trace replay. P0-B will
lower supported graph operations into NPUsim layers and mappings; unsupported
operations must be reported explicitly rather than silently omitted.
