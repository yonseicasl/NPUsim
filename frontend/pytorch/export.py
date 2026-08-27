"""PyTorch-to-NPUsim capture graph export. PyTorch is imported only on export."""

from __future__ import annotations

import hashlib
import importlib
import json
from typing import Any

from .graph_ir import SCHEMA_VERSION, dump_graph_ir, validate_graph_ir

FRONTEND_VERSION = "0.2"


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


def _tensor_record(
    tensor_id: str, value: Any, kind: str, qualified_name: str | None = None
) -> dict[str, Any]:
    torch = _torch()
    if isinstance(value, torch.Tensor):
        shape = [_shape_dimension(dimension) for dimension in value.shape]
        dtype = _dtype_name(value.dtype)
        try:
            logical_bytes = int(value.numel()) * int(value.element_size())
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"tensor {tensor_id} has a symbolic byte size; provide a concrete export example"
            ) from error
        strides = [int(stride) for stride in value.stride()]
        storage_offset = int(value.storage_offset())
        layout = "contiguous" if value.is_contiguous() else "strided"
    elif isinstance(value, bool):
        shape, dtype, logical_bytes, strides, storage_offset, layout = [], "bool", 1, [], 0, "scalar"
    elif isinstance(value, int):
        shape, dtype, logical_bytes, strides, storage_offset, layout = [], "int64", 8, [], 0, "scalar"
    elif isinstance(value, float):
        shape, dtype, logical_bytes, strides, storage_offset, layout = [], "float64", 8, [], 0, "scalar"
    else:
        raise ValueError(f"node {tensor_id} has unsupported exported value type {type(value)!r}")
    record = {
        "id": tensor_id,
        "shape": shape,
        "dtype": dtype,
        "layout": layout,
        "strides": strides,
        "storage_offset": storage_offset,
        "kind": kind,
        "logical_bytes": logical_bytes,
    }
    if qualified_name:
        record["qualified_name"] = qualified_name
    return record


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


def _json_argument(value: Any, fx_node_type: type) -> Any:
    """Serialize literal arguments without losing their position around tensor operands."""
    if isinstance(value, fx_node_type):
        return {"tensor": value.name}
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, (list, tuple)):
        return [_json_argument(item, fx_node_type) for item in value]
    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise ValueError("exported keyword dictionaries must use string keys")
        return {key: _json_argument(item, fx_node_type) for key, item in value.items()}
    # torch dtype/device/layout and enum-like export literals have stable string forms.
    module_name = getattr(type(value), "__module__", "")
    if module_name.startswith("torch"):
        return str(value)
    raise ValueError(f"unsupported exported literal argument type {type(value)!r}")


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
    """Capture a PyTorch module with torch.export and return NPUsim capture IR."""
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
    qualified_names: dict[str, str] = {}
    qualified_names.update(getattr(exported.graph_signature, "inputs_to_parameters", {}))
    qualified_names.update(getattr(exported.graph_signature, "inputs_to_buffers", {}))

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
        tensors.append(_tensor_record(node.name, value, kind, qualified_names.get(node.name)))
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
                "arguments": [_json_argument(item, fx_node_type) for item in node.args],
                "keyword_arguments": {
                    key: _json_argument(item, fx_node_type) for key, item in node.kwargs.items()
                },
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
            "literal_arguments_preserved": True,
            "tensor_layout_preserved": True,
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
    return dump_graph_ir(export_model(model, example_args, example_kwargs, **kwargs), output_path)


def load_callable(specification: str) -> Any:
    module_name, separator, attribute = specification.partition(":")
    if not separator or not module_name or not attribute:
        raise ValueError("callable specification must use module:attribute")
    value = importlib.import_module(module_name)
    for part in attribute.split("."):
        value = getattr(value, part)
    if not callable(value):
        raise TypeError(f"{specification} is not callable")
    return value
