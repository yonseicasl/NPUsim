"""PyTorch-to-NPUsim graph IR export. PyTorch is imported only on export."""

from __future__ import annotations

import hashlib
import importlib
import json
from typing import Any

from .graph_ir import SCHEMA_VERSION, dump_graph_ir, validate_graph_ir

FRONTEND_VERSION = "0.1"


def _torch() -> Any:
    try:
        return importlib.import_module("torch")
    except ModuleNotFoundError as error:
        raise RuntimeError(
            "PyTorch is required only for export. Install a compatible torch package "
            "and run the frontend from that environment."
        ) from error


def _dtype_name(dtype: Any) -> str:
    return str(dtype).removeprefix("torch.")


def _shape_dimension(dimension: Any) -> int | str:
    try:
        return int(dimension)
    except (TypeError, ValueError):
        symbol = str(dimension)
        if not symbol:
            raise ValueError("empty symbolic dimension")
        return symbol


def _tensor_record(tensor_id: str, value: Any, kind: str) -> dict[str, Any]:
    torch = _torch()
    if isinstance(value, torch.Tensor):
        shape = [_shape_dimension(dimension) for dimension in value.shape]
        dtype = _dtype_name(value.dtype)
        logical_bytes = value.numel() * value.element_size()
    elif isinstance(value, bool):
        shape, dtype, logical_bytes = [], "bool", 1
    elif isinstance(value, int):
        shape, dtype, logical_bytes = [], "int64", 8
    elif isinstance(value, float):
        shape, dtype, logical_bytes = [], "float64", 8
    else:
        raise ValueError(f"node {tensor_id} has unsupported exported value type {type(value)!r}")
    return {
        "id": tensor_id,
        "shape": shape,
        "dtype": dtype,
        "layout": "contiguous",
        "kind": kind,
        "logical_bytes": logical_bytes,
    }


def _node_references(value: Any, fx_node_type: type) -> list[str]:
    if isinstance(value, fx_node_type):
        return [value.name]
    if isinstance(value, (list, tuple)):
        result: list[str] = []
        for item in value:
            result.extend(_node_references(item, fx_node_type))
        return result
    if isinstance(value, dict):
        result = []
        for item in value.values():
            result.extend(_node_references(item, fx_node_type))
        return result
    return []


def _structure_sha256(model: Any) -> str:
    entries = []
    for name, value in model.state_dict().items():
        entries.append({"name": name, "shape": list(value.shape), "dtype": _dtype_name(value.dtype)})
    encoded = json.dumps(entries, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _input_kind(input_spec: Any) -> str:
    name = getattr(getattr(input_spec, "kind", None), "name", "USER_INPUT").lower()
    if "parameter" in name:
        return "parameter"
    if "buffer" in name:
        return "buffer"
    if "constant" in name:
        return "constant"
    return "input"


def export_model(
    model: Any,
    example_args: tuple[Any, ...],
    example_kwargs: dict[str, Any] | None = None,
    *,
    model_name: str | None = None,
    dynamic_shapes: Any = None,
) -> dict[str, Any]:
    """Capture a PyTorch module with torch.export and return NPUsim graph IR.

    This timing frontend preserves graph shape/dtype/dependency metadata. It does
    not serialize weight values and does not imply that every exported op has a
    hardware timing model yet.
    """
    torch = _torch()
    if not isinstance(model, torch.nn.Module):
        raise TypeError("model must be a torch.nn.Module")
    if not isinstance(example_args, tuple):
        raise TypeError("example_args must be a tuple")
    if example_kwargs is None:
        example_kwargs = {}
    if not isinstance(example_kwargs, dict):
        raise TypeError("example_kwargs must be a dict")

    model.eval()
    exported = torch.export.export(
        model,
        args=example_args,
        kwargs=example_kwargs,
        dynamic_shapes=dynamic_shapes,
        strict=True,
    )
    graph_module = exported.graph_module
    fx_node_type = torch.fx.Node
    kinds = {
        input_spec.arg.name: _input_kind(input_spec)
        for input_spec in exported.graph_signature.input_specs
    }

    tensors: list[dict[str, Any]] = []
    nodes: list[dict[str, Any]] = []
    inputs: list[str] = []
    outputs: list[str] = []

    for node in graph_module.graph.nodes:
        if node.op == "output":
            outputs = _node_references(node.args, fx_node_type)
            continue
        value = node.meta.get("val")
        if value is None and node.op == "get_attr":
            value = getattr(graph_module, str(node.target))
        if isinstance(value, (list, tuple, dict)):
            raise ValueError(
                f"node {node.name} has a structured result; lower it to explicit tensor outputs first"
            )
        if value is None:
            raise ValueError(f"node {node.name} is missing torch.export value metadata")

        kind = kinds.get(node.name, "activation")
        tensors.append(_tensor_record(node.name, value, kind))
        if node.op == "placeholder":
            if kind == "input":
                inputs.append(node.name)
            continue
        nodes.append(
            {
                "id": node.name,
                "op": str(node.target),
                "inputs": _node_references((node.args, node.kwargs), fx_node_type),
                "outputs": [node.name],
                "attributes": {"fx_op": node.op},
                "source": {"target": str(node.target)},
            }
        )

    if not outputs:
        raise ValueError("torch.export graph has no tensor outputs")

    graph = {
        "schema_version": SCHEMA_VERSION,
        "producer": {"name": "npusim-pytorch-frontend", "version": FRONTEND_VERSION},
        "model": {
            "name": model_name or model.__class__.__qualname__,
            "structure_sha256": _structure_sha256(model),
        },
        "metadata": {
            "pytorch_version": torch.__version__,
            "export_mode": "torch.export",
            "dynamic_shapes_configured": dynamic_shapes is not None,
            "weight_values_included": False,
        },
        "inputs": inputs,
        "outputs": outputs,
        "tensors": tensors,
        "nodes": nodes,
    }
    validate_graph_ir(graph)
    return graph


def export_to_file(
    model: Any,
    example_args: tuple[Any, ...],
    output_path: str,
    example_kwargs: dict[str, Any] | None = None,
    **kwargs: Any,
) -> str:
    """Export a model and write a validated graph IR artifact."""
    return dump_graph_ir(export_model(model, example_args, example_kwargs, **kwargs), output_path)


def load_callable(specification: str) -> Any:
    """Resolve a user factory reference in the form ``module:attribute``."""
    module_name, separator, attribute = specification.partition(":")
    if not separator or not module_name or not attribute:
        raise ValueError("callable specification must use module:attribute")
    value = importlib.import_module(module_name)
    for part in attribute.split("."):
        value = getattr(value, part)
    if not callable(value):
        raise TypeError(f"{specification} is not callable")
    return value
