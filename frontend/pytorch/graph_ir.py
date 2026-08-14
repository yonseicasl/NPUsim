"""Versioned, framework-independent graph IR validation for NPUsim."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

SCHEMA_VERSION = "npusim.graph.v1"

_DTYPE_BYTES = {
    "bool": 1,
    "uint8": 1,
    "int8": 1,
    "int16": 2,
    "int32": 4,
    "int64": 8,
    "float16": 2,
    "bfloat16": 2,
    "float32": 4,
    "float64": 8,
}


class GraphIRError(ValueError):
    """Raised when a graph artifact violates the NPUsim frontend contract."""


def _require_mapping(value: Any, field: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise GraphIRError(f"{field} must be an object")
    return value


def _require_string(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value:
        raise GraphIRError(f"{field} must be a non-empty string")
    return value


def _require_string_list(value: Any, field: str) -> list[str]:
    if not isinstance(value, list) or any(not isinstance(item, str) or not item for item in value):
        raise GraphIRError(f"{field} must be a list of non-empty strings")
    return list(value)


def _shape_numel(shape: Sequence[Any]) -> int | None:
    numel = 1
    for dimension in shape:
        if isinstance(dimension, bool):
            raise GraphIRError("tensor shape dimensions cannot be boolean")
        if isinstance(dimension, int):
            if dimension <= 0:
                raise GraphIRError("static tensor shape dimensions must be positive")
            numel *= dimension
        elif isinstance(dimension, str) and dimension:
            return None
        else:
            raise GraphIRError("tensor shape dimensions must be positive integers or symbols")
    return numel


def _validate_tensor(tensor: Mapping[str, Any], tensor_ids: set[str]) -> None:
    tensor_id = _require_string(tensor.get("id"), "tensor.id")
    if tensor_id in tensor_ids:
        raise GraphIRError(f"duplicate tensor id: {tensor_id}")
    tensor_ids.add(tensor_id)

    shape = tensor.get("shape")
    if not isinstance(shape, list):
        raise GraphIRError(f"tensor {tensor_id} shape must be a list")
    numel = _shape_numel(shape)
    dtype = _require_string(tensor.get("dtype"), f"tensor {tensor_id} dtype")
    if dtype not in _DTYPE_BYTES:
        raise GraphIRError(f"tensor {tensor_id} has unsupported IR dtype: {dtype}")
    _require_string(tensor.get("kind"), f"tensor {tensor_id} kind")
    _require_string(tensor.get("layout", "contiguous"), f"tensor {tensor_id} layout")

    logical_bytes = tensor.get("logical_bytes")
    if not isinstance(logical_bytes, int) or logical_bytes < 0:
        raise GraphIRError(f"tensor {tensor_id} logical_bytes must be a non-negative integer")
    if numel is not None and logical_bytes != numel * _DTYPE_BYTES[dtype]:
        raise GraphIRError(f"tensor {tensor_id} logical_bytes disagrees with static shape and dtype")


def validate_graph_ir(graph: Mapping[str, Any]) -> None:
    """Validate schema and dependency invariants without importing PyTorch."""
    root = _require_mapping(graph, "graph")
    if root.get("schema_version") != SCHEMA_VERSION:
        raise GraphIRError(f"schema_version must be {SCHEMA_VERSION}")

    producer = _require_mapping(root.get("producer"), "producer")
    _require_string(producer.get("name"), "producer.name")
    _require_string(producer.get("version"), "producer.version")
    model = _require_mapping(root.get("model"), "model")
    _require_string(model.get("name"), "model.name")
    _require_string(model.get("structure_sha256"), "model.structure_sha256")

    tensors = root.get("tensors")
    if not isinstance(tensors, list) or not tensors:
        raise GraphIRError("tensors must be a non-empty list")
    tensor_ids: set[str] = set()
    for tensor in tensors:
        _validate_tensor(_require_mapping(tensor, "tensor"), tensor_ids)

    input_ids = _require_string_list(root.get("inputs"), "inputs")
    output_ids = _require_string_list(root.get("outputs"), "outputs")
    if any(tensor_id not in tensor_ids for tensor_id in input_ids + output_ids):
        raise GraphIRError("graph inputs and outputs must reference declared tensors")

    nodes = root.get("nodes")
    if not isinstance(nodes, list) or not nodes:
        raise GraphIRError("nodes must be a non-empty list")
    node_ids: set[str] = set()
    produced: set[str] = set(input_ids)
    for tensor in tensors:
        if tensor["kind"] in {"parameter", "buffer", "constant"}:
            produced.add(tensor["id"])

    for node in nodes:
        node_map = _require_mapping(node, "node")
        node_id = _require_string(node_map.get("id"), "node.id")
        if node_id in node_ids:
            raise GraphIRError(f"duplicate node id: {node_id}")
        node_ids.add(node_id)
        _require_string(node_map.get("op"), f"node {node_id} op")
        inputs = _require_string_list(node_map.get("inputs", []), f"node {node_id} inputs")
        outputs = _require_string_list(node_map.get("outputs"), f"node {node_id} outputs")
        if not outputs:
            raise GraphIRError(f"node {node_id} must produce at least one tensor")
        if any(tensor_id not in tensor_ids for tensor_id in inputs + outputs):
            raise GraphIRError(f"node {node_id} references an undeclared tensor")
        if any(tensor_id not in produced for tensor_id in inputs):
            raise GraphIRError(f"node {node_id} is not topologically ordered")
        if any(tensor_id in produced for tensor_id in outputs):
            raise GraphIRError(f"node {node_id} redefines tensor output")
        if not isinstance(node_map.get("attributes", {}), Mapping):
            raise GraphIRError(f"node {node_id} attributes must be an object")
        produced.update(outputs)

    if any(tensor_id not in produced for tensor_id in output_ids):
        raise GraphIRError("graph outputs must be produced by a node or declared input")


def canonical_json(graph: Mapping[str, Any]) -> str:
    validate_graph_ir(graph)
    return json.dumps(graph, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def graph_sha256(graph: Mapping[str, Any]) -> str:
    payload = dict(graph)
    payload.pop("graph_sha256", None)
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()


def load_graph_ir(path: str | Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as source:
        graph = json.load(source)
    validate_graph_ir(graph)
    return graph


def dump_graph_ir(graph: Mapping[str, Any], path: str | Path) -> str:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(graph)
    payload["graph_sha256"] = graph_sha256(payload)
    with destination.open("w", encoding="utf-8") as target:
        json.dump(payload, target, indent=2, sort_keys=True)
        target.write("\n")
    return payload["graph_sha256"]
