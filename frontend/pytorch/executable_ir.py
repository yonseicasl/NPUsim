"""Validated, framework-neutral executable IR for the NPUsim C++ core."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

EXECUTABLE_SCHEMA_VERSION = "npusim.exec.v1"
EXECUTABLE_PRODUCER_VERSION = "0.1"

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

_MODELED_KINDS = {"npusim.linear", "npusim.conv2d", "npusim.softmax"}
_ACTIVATIONS = {"linear", "relu", "leaky"}


class ExecutableIRError(ValueError):
    """Raised when an executable artifact cannot be simulated safely."""


def _mapping(value: Any, field: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ExecutableIRError(f"{field} must be an object")
    return value


def _string(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value:
        raise ExecutableIRError(f"{field} must be a non-empty string")
    return value


def _strings(value: Any, field: str) -> list[str]:
    if not isinstance(value, list) or any(not isinstance(item, str) or not item for item in value):
        raise ExecutableIRError(f"{field} must be a list of non-empty strings")
    return list(value)


def _positive_integer(value: Any, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ExecutableIRError(f"{field} must be a positive integer")
    return value


def _nonnegative_integer(value: Any, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ExecutableIRError(f"{field} must be a non-negative integer")
    return value


def _numel(shape: Sequence[Any], tensor_id: str) -> int:
    total = 1
    for dimension in shape:
        total *= _positive_integer(dimension, f"tensor {tensor_id} shape dimension")
    return total


def _validate_tensor(tensor: Mapping[str, Any], tensor_ids: set[str]) -> None:
    tensor_id = _string(tensor.get("id"), "tensor.id")
    if tensor_id in tensor_ids:
        raise ExecutableIRError(f"duplicate tensor id: {tensor_id}")
    tensor_ids.add(tensor_id)
    shape = tensor.get("shape")
    if not isinstance(shape, list):
        raise ExecutableIRError(f"tensor {tensor_id} shape must be a list")
    numel = _numel(shape, tensor_id)
    dtype = _string(tensor.get("dtype"), f"tensor {tensor_id} dtype")
    if dtype not in _DTYPE_BYTES:
        raise ExecutableIRError(f"tensor {tensor_id} has unsupported executable dtype: {dtype}")
    logical_bytes = _nonnegative_integer(
        tensor.get("logical_bytes"), f"tensor {tensor_id} logical_bytes"
    )
    if logical_bytes != numel * _DTYPE_BYTES[dtype]:
        raise ExecutableIRError(
            f"tensor {tensor_id} logical_bytes disagrees with its concrete shape and dtype"
        )
    _string(tensor.get("kind"), f"tensor {tensor_id} kind")
    _string(tensor.get("layout", "contiguous"), f"tensor {tensor_id} layout")


def _validate_geometry(operation: Mapping[str, Any]) -> None:
    operation_id = operation["id"]
    kind = operation["kind"]
    geometry = _mapping(operation.get("geometry"), f"operation {operation_id} geometry")
    if kind == "npusim.linear":
        for field in ("batch", "input_features", "output_features"):
            _positive_integer(geometry.get(field), f"operation {operation_id} geometry.{field}")
    elif kind == "npusim.conv2d":
        positive = (
            "batch",
            "input_channels",
            "input_height",
            "input_width",
            "output_channels",
            "output_height",
            "output_width",
            "filter_height",
            "filter_width",
            "stride_height",
            "stride_width",
            "dilation_height",
            "dilation_width",
            "groups",
        )
        for field in positive:
            _positive_integer(geometry.get(field), f"operation {operation_id} geometry.{field}")
        _nonnegative_integer(
            geometry.get("padding_height"), f"operation {operation_id} geometry.padding_height"
        )
        _nonnegative_integer(
            geometry.get("padding_width"), f"operation {operation_id} geometry.padding_width"
        )
        if geometry["input_channels"] % geometry["groups"]:
            raise ExecutableIRError(f"operation {operation_id} groups must divide input_channels")
        if geometry["output_channels"] % geometry["groups"]:
            raise ExecutableIRError(f"operation {operation_id} groups must divide output_channels")
    elif kind == "npusim.softmax":
        _positive_integer(geometry.get("rows"), f"operation {operation_id} geometry.rows")
        _positive_integer(
            geometry.get("row_length"), f"operation {operation_id} geometry.row_length"
        )


def validate_executable_ir(executable: Mapping[str, Any]) -> None:
    root = _mapping(executable, "executable")
    if root.get("schema_version") != EXECUTABLE_SCHEMA_VERSION:
        raise ExecutableIRError(f"schema_version must be {EXECUTABLE_SCHEMA_VERSION}")
    producer = _mapping(root.get("producer"), "producer")
    _string(producer.get("name"), "producer.name")
    _string(producer.get("version"), "producer.version")
    model = _mapping(root.get("model"), "model")
    _string(model.get("name"), "model.name")
    _string(model.get("structure_sha256"), "model.structure_sha256")
    source = _mapping(root.get("source_graph"), "source_graph")
    _string(source.get("schema_version"), "source_graph.schema_version")
    _string(source.get("sha256"), "source_graph.sha256")

    tensors = root.get("tensors")
    if not isinstance(tensors, list) or not tensors:
        raise ExecutableIRError("tensors must be a non-empty list")
    tensor_ids: set[str] = set()
    for tensor in tensors:
        _validate_tensor(_mapping(tensor, "tensor"), tensor_ids)

    inputs = _strings(root.get("inputs"), "inputs")
    outputs = _strings(root.get("outputs"), "outputs")
    if any(tensor_id not in tensor_ids for tensor_id in inputs + outputs):
        raise ExecutableIRError("inputs and outputs must reference declared tensors")

    operations = root.get("operations")
    if not isinstance(operations, list) or not operations:
        raise ExecutableIRError("operations must be a non-empty list")
    operation_ids: set[str] = set()
    produced = set(inputs)
    for tensor in tensors:
        if tensor["kind"] in {"parameter", "buffer", "constant"}:
            produced.add(tensor["id"])
    for operation_value in operations:
        operation = _mapping(operation_value, "operation")
        operation_id = _string(operation.get("id"), "operation.id")
        if operation_id in operation_ids:
            raise ExecutableIRError(f"duplicate operation id: {operation_id}")
        operation_ids.add(operation_id)
        kind = _string(operation.get("kind"), f"operation {operation_id} kind")
        status = _string(operation.get("status"), f"operation {operation_id} status")
        if status != "modeled" or kind not in _MODELED_KINDS:
            raise ExecutableIRError(
                f"operation {operation_id} is not executable: kind={kind}, status={status}"
            )
        operation_inputs = _strings(operation.get("inputs", []), f"operation {operation_id} inputs")
        operation_outputs = _strings(operation.get("outputs"), f"operation {operation_id} outputs")
        if not operation_outputs:
            raise ExecutableIRError(f"operation {operation_id} must produce at least one tensor")
        if any(tensor_id not in tensor_ids for tensor_id in operation_inputs + operation_outputs):
            raise ExecutableIRError(f"operation {operation_id} references an undeclared tensor")
        if any(tensor_id not in produced for tensor_id in operation_inputs):
            raise ExecutableIRError(f"operation {operation_id} is not topologically ordered")
        if any(tensor_id in produced for tensor_id in operation_outputs):
            raise ExecutableIRError(f"operation {operation_id} redefines a tensor")
        activation = _string(operation.get("activation", "linear"), f"operation {operation_id} activation")
        if activation not in _ACTIVATIONS:
            raise ExecutableIRError(
                f"operation {operation_id} has unsupported fused activation: {activation}"
            )
        source_nodes = _strings(
            operation.get("source_nodes", [operation_id]), f"operation {operation_id} source_nodes"
        )
        if not source_nodes:
            raise ExecutableIRError(f"operation {operation_id} source_nodes cannot be empty")
        _validate_geometry(operation)
        produced.update(operation_outputs)

    if any(tensor_id not in produced for tensor_id in outputs):
        raise ExecutableIRError("graph outputs are not produced by executable operations")


def canonical_executable_json(executable: Mapping[str, Any]) -> str:
    validate_executable_ir(executable)
    return json.dumps(executable, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def executable_sha256(executable: Mapping[str, Any]) -> str:
    payload = dict(executable)
    payload.pop("executable_sha256", None)
    return hashlib.sha256(canonical_executable_json(payload).encode("utf-8")).hexdigest()


def load_executable_ir(path: str | Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as source:
        executable = json.load(source)
    validate_executable_ir(executable)
    declared = executable.get("executable_sha256")
    if declared is not None and declared != executable_sha256(executable):
        raise ExecutableIRError("executable_sha256 does not match the executable artifact")
    return executable


def dump_executable_ir(executable: Mapping[str, Any], path: str | Path) -> str:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(executable)
    payload["executable_sha256"] = executable_sha256(payload)
    with destination.open("w", encoding="utf-8") as target:
        json.dump(payload, target, indent=2, sort_keys=True)
        target.write("\n")
    return payload["executable_sha256"]
