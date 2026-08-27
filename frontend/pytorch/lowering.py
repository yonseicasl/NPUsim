"""Strict lowering from PyTorch capture IR to NPUsim executable IR."""

from __future__ import annotations

import math
from typing import Any, Mapping

from .executable_ir import (
    EXECUTABLE_PRODUCER_VERSION,
    EXECUTABLE_SCHEMA_VERSION,
    validate_executable_ir,
)
from .graph_ir import graph_sha256, validate_graph_ir


class LoweringError(ValueError):
    """Raised when a captured operation has no safe executable lowering."""


_LINEAR_OPS = {"aten.linear.default"}
_CONV_OPS = {"aten.conv2d.default", "aten.convolution.default"}
_ACTIVATION_OPS = {
    "aten.relu.default": "relu",
    "aten.leaky_relu.default": "leaky",
}
_SOFTMAX_OPS = {"aten.softmax.int", "aten._softmax.default"}


def _product(values: list[int]) -> int:
    return math.prod(values) if values else 1


def _static_shape(tensor: Mapping[str, Any], context: str) -> list[int]:
    shape = tensor.get("shape")
    if not isinstance(shape, list) or any(
        isinstance(value, bool) or not isinstance(value, int) for value in shape
    ):
        raise LoweringError(f"{context} needs a concrete static shape binding")
    if any(value <= 0 for value in shape):
        raise LoweringError(f"{context} has a non-positive tensor dimension")
    return list(shape)


def _pair(value: Any, name: str, default: int | None = None) -> tuple[int, int]:
    if value is None:
        if default is None:
            raise LoweringError(f"missing required convolution attribute: {name}")
        return default, default
    if isinstance(value, bool):
        raise LoweringError(f"{name} cannot be boolean")
    if isinstance(value, int):
        return value, value
    if isinstance(value, (list, tuple)) and len(value) == 2:
        first, second = value
        if all(isinstance(item, int) and not isinstance(item, bool) for item in (first, second)):
            return first, second
    raise LoweringError(f"{name} must be an integer or a two-element integer list")


def _attributes(node: Mapping[str, Any]) -> Mapping[str, Any]:
    value = node.get("attributes", {})
    if not isinstance(value, Mapping):
        raise LoweringError(f"node {node.get('id')} attributes must be an object")
    return value


def _linear_operation(
    node: Mapping[str, Any],
    tensors: Mapping[str, Mapping[str, Any]],
    output_id: str,
    activation: str,
    source_nodes: list[str],
) -> dict[str, Any]:
    inputs = node.get("inputs", [])
    if not isinstance(inputs, list) or len(inputs) not in {2, 3}:
        raise LoweringError(f"linear node {node['id']} requires input, weight, and optional bias")
    input_shape = _static_shape(tensors[inputs[0]], f"linear node {node['id']} input")
    weight_shape = _static_shape(tensors[inputs[1]], f"linear node {node['id']} weight")
    output_shape = _static_shape(tensors[output_id], f"linear node {node['id']} output")
    if len(input_shape) < 1 or len(weight_shape) != 2 or len(output_shape) < 1:
        raise LoweringError(f"linear node {node['id']} has invalid tensor rank")
    input_features = input_shape[-1]
    output_features, weight_input_features = weight_shape
    batch = _product(input_shape[:-1])
    if input_features != weight_input_features:
        raise LoweringError(f"linear node {node['id']} input and weight dimensions disagree")
    if _product(output_shape[:-1]) != batch or output_shape[-1] != output_features:
        raise LoweringError(f"linear node {node['id']} output shape disagrees with its geometry")
    return {
        "id": str(node["id"]),
        "kind": "npusim.linear",
        "status": "modeled",
        "inputs": list(inputs),
        "outputs": [output_id],
        "activation": activation,
        "mapping_required": True,
        "geometry": {
            "batch": batch,
            "input_features": input_features,
            "output_features": output_features,
        },
        "source_nodes": source_nodes,
    }


def _conv_attributes(node: Mapping[str, Any]) -> dict[str, Any]:
    attributes = dict(_attributes(node))
    arguments = node.get("arguments")
    # torch.export may emit either the public conv2d schema or the lower-level
    # convolution schema. Their trailing positional arguments are different.
    # Capture-v1 fixtures may instead provide named attributes.
    if isinstance(arguments, list):
        if len(arguments) > 3:
            attributes.setdefault("stride", arguments[3])
        if len(arguments) > 4:
            attributes.setdefault("padding", arguments[4])
        if len(arguments) > 5:
            attributes.setdefault("dilation", arguments[5])
        if node.get("op") == "aten.convolution.default":
            if len(arguments) > 6:
                attributes.setdefault("transposed", arguments[6])
            if len(arguments) > 7:
                attributes.setdefault("output_padding", arguments[7])
            if len(arguments) > 8:
                attributes.setdefault("groups", arguments[8])
        elif len(arguments) > 6:
            attributes.setdefault("groups", arguments[6])
    return attributes


def _conv_operation(
    node: Mapping[str, Any],
    tensors: Mapping[str, Mapping[str, Any]],
    output_id: str,
    activation: str,
    source_nodes: list[str],
) -> dict[str, Any]:
    inputs = node.get("inputs", [])
    if not isinstance(inputs, list) or len(inputs) not in {2, 3}:
        raise LoweringError(f"conv2d node {node['id']} requires input, weight, and optional bias")
    input_shape = _static_shape(tensors[inputs[0]], f"conv2d node {node['id']} input")
    weight_shape = _static_shape(tensors[inputs[1]], f"conv2d node {node['id']} weight")
    output_shape = _static_shape(tensors[output_id], f"conv2d node {node['id']} output")
    if len(input_shape) != 4 or len(weight_shape) != 4 or len(output_shape) != 4:
        raise LoweringError(f"conv2d node {node['id']} currently requires NCHW/OIHW rank-4 tensors")
    batch, input_channels, input_height, input_width = input_shape
    output_channels, weight_channels, filter_height, filter_width = weight_shape
    attributes = _conv_attributes(node)
    if attributes.get("transposed", False) is not False:
        raise LoweringError(f"conv2d node {node['id']} does not support transposed convolution")
    output_padding = _pair(attributes.get("output_padding"), "output_padding", 0)
    if output_padding != (0, 0):
        raise LoweringError(f"conv2d node {node['id']} does not support output padding")
    stride_h, stride_w = _pair(attributes.get("stride"), "stride", 1)
    padding_h, padding_w = _pair(attributes.get("padding"), "padding", 0)
    dilation_h, dilation_w = _pair(attributes.get("dilation"), "dilation", 1)
    groups = attributes.get("groups", 1)
    if isinstance(groups, bool) or not isinstance(groups, int) or groups <= 0:
        raise LoweringError("groups must be a positive integer")
    if input_channels % groups or output_channels % groups:
        raise LoweringError(f"conv2d node {node['id']} groups do not divide channel dimensions")
    if weight_channels != input_channels // groups:
        raise LoweringError(f"conv2d node {node['id']} weight channels disagree with groups")
    expected_h = (
        input_height + 2 * padding_h - dilation_h * (filter_height - 1) - 1
    ) // stride_h + 1
    expected_w = (
        input_width + 2 * padding_w - dilation_w * (filter_width - 1) - 1
    ) // stride_w + 1
    if output_shape != [batch, output_channels, expected_h, expected_w]:
        raise LoweringError(f"conv2d node {node['id']} output shape disagrees with its geometry")
    return {
        "id": str(node["id"]),
        "kind": "npusim.conv2d",
        "status": "modeled",
        "inputs": list(inputs),
        "outputs": [output_id],
        "activation": activation,
        "mapping_required": True,
        "geometry": {
            "batch": batch,
            "input_channels": input_channels,
            "input_height": input_height,
            "input_width": input_width,
            "output_channels": output_channels,
            "output_height": expected_h,
            "output_width": expected_w,
            "filter_height": filter_height,
            "filter_width": filter_width,
            "stride_height": stride_h,
            "stride_width": stride_w,
            "padding_height": padding_h,
            "padding_width": padding_w,
            "dilation_height": dilation_h,
            "dilation_width": dilation_w,
            "groups": groups,
        },
        "source_nodes": source_nodes,
    }


def _softmax_operation(
    node: Mapping[str, Any], tensors: Mapping[str, Mapping[str, Any]]
) -> dict[str, Any]:
    inputs = node.get("inputs", [])
    outputs = node.get("outputs", [])
    if (
        not isinstance(inputs, list)
        or len(inputs) != 1
        or not isinstance(outputs, list)
        or len(outputs) != 1
    ):
        raise LoweringError(f"softmax node {node['id']} requires one input and one output")
    shape = _static_shape(tensors[inputs[0]], f"softmax node {node['id']} input")
    output_shape = _static_shape(tensors[outputs[0]], f"softmax node {node['id']} output")
    if shape != output_shape or not shape:
        raise LoweringError(f"softmax node {node['id']} input/output shape mismatch")
    attributes = _attributes(node)
    axis = attributes.get("axis", attributes.get("dim", -1))
    arguments = node.get("arguments")
    if isinstance(arguments, list) and len(arguments) > 1:
        axis = arguments[1]
    if not isinstance(axis, int) or isinstance(axis, bool):
        raise LoweringError(f"softmax node {node['id']} axis must be an integer")
    normalized_axis = axis if axis >= 0 else len(shape) + axis
    if normalized_axis != len(shape) - 1:
        raise LoweringError(f"softmax node {node['id']} currently supports only the last axis")
    return {
        "id": str(node["id"]),
        "kind": "npusim.softmax",
        "status": "modeled",
        "inputs": list(inputs),
        "outputs": list(outputs),
        "activation": "linear",
        "mapping_required": False,
        "geometry": {"rows": _product(shape[:-1]), "row_length": shape[-1]},
        "source_nodes": [str(node["id"])],
    }


def lower_graph(graph: Mapping[str, Any]) -> dict[str, Any]:
    """Lower a validated static capture graph, rejecting every unsupported node."""
    validate_graph_ir(graph)
    tensors = {str(tensor["id"]): tensor for tensor in graph["tensors"]}
    nodes = list(graph["nodes"])
    consumers: dict[str, list[Mapping[str, Any]]] = {}
    for node in nodes:
        for tensor_id in node.get("inputs", []):
            consumers.setdefault(str(tensor_id), []).append(node)

    fused_activation_by_node: dict[str, tuple[str, Mapping[str, Any]]] = {}
    fused_node_ids: set[str] = set()
    for node in nodes:
        if node.get("op") not in _LINEAR_OPS | _CONV_OPS:
            continue
        outputs = node.get("outputs", [])
        if not isinstance(outputs, list) or len(outputs) != 1:
            raise LoweringError(f"node {node.get('id')} must have exactly one output")
        output_consumers = consumers.get(str(outputs[0]), [])
        if len(output_consumers) != 1:
            continue
        candidate = output_consumers[0]
        activation = _ACTIVATION_OPS.get(str(candidate.get("op")))
        if activation is None or candidate.get("inputs") != outputs:
            continue
        fused_activation_by_node[str(node["id"])] = (activation, candidate)
        fused_node_ids.add(str(candidate["id"]))

    operations: list[dict[str, Any]] = []
    lowered_sources: list[str] = []
    for node in nodes:
        node_id = str(node["id"])
        if node_id in fused_node_ids:
            continue
        operation_name = str(node.get("op"))
        outputs = node.get("outputs", [])
        if not isinstance(outputs, list) or len(outputs) != 1:
            raise LoweringError(f"node {node_id} must have exactly one tensor output")
        activation = "linear"
        output_id = str(outputs[0])
        source_nodes = [node_id]
        if node_id in fused_activation_by_node:
            activation, activation_node = fused_activation_by_node[node_id]
            activation_outputs = activation_node.get("outputs", [])
            if not isinstance(activation_outputs, list) or len(activation_outputs) != 1:
                raise LoweringError(f"activation node {activation_node.get('id')} must have one output")
            output_id = str(activation_outputs[0])
            source_nodes.append(str(activation_node["id"]))
        if operation_name in _LINEAR_OPS:
            operation = _linear_operation(node, tensors, output_id, activation, source_nodes)
        elif operation_name in _CONV_OPS:
            operation = _conv_operation(node, tensors, output_id, activation, source_nodes)
        elif operation_name in _SOFTMAX_OPS:
            operation = _softmax_operation(node, tensors)
        elif operation_name in _ACTIVATION_OPS:
            raise LoweringError(
                f"standalone activation node {node_id} cannot yet be scheduled; "
                "fuse it with Linear/Conv2d"
            )
        else:
            raise LoweringError(f"unsupported PyTorch operation {operation_name} at node {node_id}")
        operations.append(operation)
        lowered_sources.extend(source_nodes)

    executable = {
        "schema_version": EXECUTABLE_SCHEMA_VERSION,
        "producer": {
            "name": "npusim-pytorch-lowering",
            "version": EXECUTABLE_PRODUCER_VERSION,
        },
        "model": dict(graph["model"]),
        "source_graph": {
            "schema_version": graph["schema_version"],
            "sha256": graph_sha256(graph),
        },
        "inputs": list(graph["inputs"]),
        "outputs": list(graph["outputs"]),
        "tensors": [dict(tensor) for tensor in graph["tensors"]],
        "operations": operations,
        "coverage": {
            "captured_nodes": len(nodes),
            "lowered_source_nodes": len(lowered_sources),
            "executable_operations": len(operations),
            "unsupported_nodes": 0,
        },
    }
    validate_executable_ir(executable)
    return executable
