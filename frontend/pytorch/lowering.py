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
_POOL_OPS = {
    "aten.max_pool2d.default": "max",
    "aten.avg_pool2d.default": "average",
}
_ELEMENTWISE_OPS = {"aten.add.Tensor": "add", "aten.mul.Tensor": "multiply"}
_CONCAT_OPS = {"aten.cat.default"}
_BATCH_NORM_OPS = {"aten.batch_norm.default"}
_LAYER_NORM_OPS = {"aten.layer_norm.default", "aten.native_layer_norm.default"}
_MATMUL_OPS = {"aten.matmul.default", "aten.bmm.default"}
# B-4: a last-two-dim transpose feeding a matmul's second operand (the q@k.transpose(-1,-2)
# attention idiom) is folded into the matmul's transpose_b flag and elided.
_TRANSPOSE_OPS = {"aten.transpose.int", "aten.t.default"}
# Element-preserving reshapes of contiguous whole-storage tensors are aliases, not work:
# they are elided during lowering and never become executable operations. The output
# tensor is annotated with alias_of so residency/lifetime shares the producer's storage.
_RESHAPE_OPS = {
    "aten.flatten.using_ints",
    "aten.view.default",
    "aten.reshape.default",
}


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


def _whole_contiguous_shape(tensor: Mapping[str, Any], context: str) -> list[int]:
    """Shape of a tensor that provably covers its whole storage in row-major order."""
    shape = _static_shape(tensor, context)
    if tensor.get("layout", "contiguous") != "contiguous":
        raise LoweringError(f"{context} must be contiguous to elide as an alias")
    if tensor.get("storage_offset", 0) != 0:
        raise LoweringError(f"{context} views storage at a nonzero offset; not a pure rename")
    strides = tensor.get("strides")
    if strides is not None:
        if not isinstance(strides, list) or len(strides) != len(shape):
            raise LoweringError(f"{context} strides disagree with its shape")
        expected = 1
        for index in range(len(shape) - 1, -1, -1):
            if shape[index] != 1 and strides[index] != expected:
                raise LoweringError(
                    f"{context} is not stored contiguously; eliding it would hide a real copy"
                )
            expected *= shape[index]
    return shape


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



def _pool_operation(
    node: Mapping[str, Any], tensors: Mapping[str, Mapping[str, Any]], mode: str
) -> dict[str, Any]:
    inputs, outputs = node.get("inputs", []), node.get("outputs", [])
    if not isinstance(inputs, list) or len(inputs) != 1 or not isinstance(outputs, list) or len(outputs) != 1:
        raise LoweringError(f"pool node {node['id']} requires one input and one output")
    input_shape = _static_shape(tensors[inputs[0]], f"pool node {node['id']} input")
    output_shape = _static_shape(tensors[outputs[0]], f"pool node {node['id']} output")
    if len(input_shape) != 4 or len(output_shape) != 4:
        raise LoweringError(f"pool node {node['id']} currently requires NCHW tensors")
    attributes = dict(_attributes(node))
    arguments = node.get("arguments")
    if isinstance(arguments, list):
        if len(arguments) > 1:
            attributes.setdefault("kernel_size", arguments[1])
        if len(arguments) > 2:
            attributes.setdefault("stride", arguments[2])
        if len(arguments) > 3:
            attributes.setdefault("padding", arguments[3])
        if mode == "max" and len(arguments) > 4:
            attributes.setdefault("dilation", arguments[4])
        if mode == "max" and len(arguments) > 5:
            attributes.setdefault("ceil_mode", arguments[5])
        if mode == "average" and len(arguments) > 4:
            attributes.setdefault("ceil_mode", arguments[4])
        if mode == "average" and len(arguments) > 5:
            attributes.setdefault("count_include_pad", arguments[5])
    if attributes.get("ceil_mode", False):
        raise LoweringError(f"pool node {node['id']} does not support ceil_mode")
    kernel_h, kernel_w = _pair(attributes.get("kernel_size"), "kernel_size")
    stride_value = attributes.get("stride")
    if stride_value is None or stride_value == []:
        stride_h, stride_w = kernel_h, kernel_w
    else:
        stride_h, stride_w = _pair(stride_value, "stride")
    padding_h, padding_w = _pair(attributes.get("padding"), "padding", 0)
    dilation_h, dilation_w = _pair(attributes.get("dilation"), "dilation", 1)
    if min(kernel_h, kernel_w, stride_h, stride_w, dilation_h, dilation_w) <= 0 or min(padding_h, padding_w) < 0:
        raise LoweringError(f"pool node {node['id']} has invalid kernel/stride/padding/dilation")
    batch, channels, input_height, input_width = input_shape
    effective_h = dilation_h * (kernel_h - 1) + 1
    effective_w = dilation_w * (kernel_w - 1) + 1
    expected_h = (input_height + 2 * padding_h - effective_h) // stride_h + 1
    expected_w = (input_width + 2 * padding_w - effective_w) // stride_w + 1
    if output_shape != [batch, channels, expected_h, expected_w]:
        raise LoweringError(f"pool node {node['id']} output shape disagrees with its geometry")
    return {
        "id": str(node["id"]), "kind": "npusim.pool2d", "status": "modeled",
        "inputs": list(inputs), "outputs": list(outputs), "activation": "linear",
        "mapping_required": False,
        "geometry": {
            "mode": mode, "batch": batch, "channels": channels,
            "input_height": input_height, "input_width": input_width,
            "output_height": expected_h, "output_width": expected_w,
            "kernel_height": kernel_h, "kernel_width": kernel_w,
            "stride_height": stride_h, "stride_width": stride_w,
            "padding_height": padding_h, "padding_width": padding_w,
            "dilation_height": dilation_h, "dilation_width": dilation_w,
            "count_include_pad": bool(attributes.get("count_include_pad", True)),
        },
        "source_nodes": [str(node["id"])],
    }


def _elementwise_operation(
    node: Mapping[str, Any], tensors: Mapping[str, Mapping[str, Any]], operator: str
) -> dict[str, Any]:
    inputs, outputs = node.get("inputs", []), node.get("outputs", [])
    if not isinstance(inputs, list) or len(inputs) != 2 or not isinstance(outputs, list) or len(outputs) != 1:
        raise LoweringError(f"elementwise node {node['id']} requires two tensors and one output")
    shapes = [_static_shape(tensors[value], f"elementwise node {node['id']} input") for value in inputs]
    output_shape = _static_shape(tensors[outputs[0]], f"elementwise node {node['id']} output")
    if shapes[0] != shapes[1] or shapes[0] != output_shape:
        raise LoweringError(f"elementwise node {node['id']} currently requires equal shapes (no broadcasting)")
    arguments = node.get("arguments")
    if operator == "add" and isinstance(arguments, list) and len(arguments) > 2 and arguments[2] != 1:
        raise LoweringError(f"elementwise add node {node['id']} requires alpha=1")
    return {
        "id": str(node["id"]), "kind": "npusim.elementwise", "status": "modeled",
        "inputs": list(inputs), "outputs": list(outputs), "activation": "linear",
        "mapping_required": False,
        "geometry": {"operator": operator, "elements": _product(output_shape)},
        "source_nodes": [str(node["id"])],
    }


def _concat_operation(node: Mapping[str, Any], tensors: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    inputs, outputs = node.get("inputs", []), node.get("outputs", [])
    if not isinstance(inputs, list) or len(inputs) < 2 or not isinstance(outputs, list) or len(outputs) != 1:
        raise LoweringError(f"concat node {node['id']} requires at least two inputs and one output")
    shapes = [_static_shape(tensors[value], f"concat node {node['id']} input") for value in inputs]
    output_shape = _static_shape(tensors[outputs[0]], f"concat node {node['id']} output")
    rank = len(output_shape)
    arguments = node.get("arguments")
    attributes = _attributes(node)
    axis = attributes.get("axis", attributes.get("dim", 0))
    if isinstance(arguments, list) and len(arguments) > 1:
        axis = arguments[1]
    if not isinstance(axis, int) or isinstance(axis, bool):
        raise LoweringError(f"concat node {node['id']} axis must be an integer")
    axis = axis if axis >= 0 else rank + axis
    if rank == 0 or axis < 0 or axis >= rank or any(len(shape) != rank for shape in shapes):
        raise LoweringError(f"concat node {node['id']} has an invalid axis or tensor rank")
    expected = list(shapes[0])
    expected[axis] = sum(shape[axis] for shape in shapes)
    for shape in shapes[1:]:
        if any(shape[index] != shapes[0][index] for index in range(rank) if index != axis):
            raise LoweringError(f"concat node {node['id']} non-concat dimensions disagree")
    if expected != output_shape:
        raise LoweringError(f"concat node {node['id']} output shape disagrees with its inputs")
    return {
        "id": str(node["id"]), "kind": "npusim.concat", "status": "modeled",
        "inputs": list(inputs), "outputs": list(outputs), "activation": "linear",
        "mapping_required": False, "geometry": {"axis": axis, "elements": _product(output_shape)},
        "source_nodes": [str(node["id"])],
    }


def _batch_norm_operation(node: Mapping[str, Any], tensors: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    inputs, outputs = node.get("inputs", []), node.get("outputs", [])
    if (
        not isinstance(inputs, list)
        or len(inputs) < 3
        or len(inputs) > 5
        or not isinstance(outputs, list)
        or len(outputs) != 1
    ):
        raise LoweringError(
            f"batch_norm node {node['id']} requires activation/stat tensors and one output"
        )
    input_shape = _static_shape(tensors[inputs[0]], f"batch_norm node {node['id']} input")
    output_shape = _static_shape(tensors[outputs[0]], f"batch_norm node {node['id']} output")
    if len(input_shape) < 2 or input_shape != output_shape:
        raise LoweringError(f"batch_norm node {node['id']} input/output shape mismatch")
    channels = input_shape[1]
    for tensor_id in inputs[1:]:
        if _static_shape(tensors[tensor_id], f"batch_norm node {node['id']} parameter") != [channels]:
            raise LoweringError(f"batch_norm node {node['id']} parameter/stat shape is invalid")
    arguments = node.get("arguments")
    training, epsilon = False, 1e-5
    if isinstance(arguments, list):
        if len(arguments) > 5:
            training = arguments[5]
        if len(arguments) > 7:
            epsilon = arguments[7]
    attributes = _attributes(node)
    training = attributes.get("training", training)
    epsilon = attributes.get("epsilon", attributes.get("eps", epsilon))
    if training is not False:
        raise LoweringError(f"batch_norm node {node['id']} supports inference mode only")
    if not isinstance(epsilon, (int, float)) or isinstance(epsilon, bool) or epsilon <= 0:
        raise LoweringError(f"batch_norm node {node['id']} epsilon must be positive")
    return {
        "id": str(node["id"]), "kind": "npusim.batch_norm", "status": "modeled",
        "inputs": list(inputs), "outputs": list(outputs), "activation": "linear",
        "mapping_required": False,
        "geometry": {
            "elements": _product(output_shape),
            "channels": channels,
            "epsilon": float(epsilon),
        },
        "source_nodes": [str(node["id"])],
    }


def _layer_norm_operation(node: Mapping[str, Any], tensors: Mapping[str, Any]) -> dict[str, Any]:
    # B-4: aten.layer_norm (input, normalized_shape, weight, bias, eps). Affine form only.
    inputs, outputs = node.get("inputs", []), node.get("outputs", [])
    if not isinstance(inputs, list) or len(inputs) < 3 or len(outputs) != 1:
        raise LoweringError(f"layer_norm node {node['id']} requires (input, weight, bias)")
    input_shape = _static_shape(tensors[inputs[0]], f"layer_norm node {node['id']} input")
    output_shape = _static_shape(tensors[outputs[0]], f"layer_norm node {node['id']} output")
    if not input_shape or input_shape != output_shape:
        raise LoweringError(f"layer_norm node {node['id']} input/output shape mismatch")
    # Data inputs are (input, weight, bias); weight/bias are the last two tensor operands.
    data_inputs = [t for t in inputs if isinstance(t, str) and t in tensors]
    weight, bias = data_inputs[-2], data_inputs[-1]
    normalized = _static_shape(tensors[weight], f"layer_norm node {node['id']} weight")
    if len(normalized) != 1 or normalized != _static_shape(tensors[bias], "layer_norm bias"):
        raise LoweringError(f"layer_norm node {node['id']} weight/bias must be 1-D of equal size")
    if input_shape[-1] != normalized[0]:
        raise LoweringError(f"layer_norm node {node['id']} normalizes the last axis only")
    attributes = _attributes(node)
    epsilon = attributes.get("eps", attributes.get("epsilon", 1e-5))
    arguments = node.get("arguments")
    if isinstance(arguments, list) and len(arguments) > 4 and isinstance(arguments[4], (int, float)):
        epsilon = arguments[4]
    if not isinstance(epsilon, (int, float)) or isinstance(epsilon, bool) or epsilon <= 0:
        raise LoweringError(f"layer_norm node {node['id']} eps must be positive")
    return {
        "id": str(node["id"]), "kind": "npusim.layer_norm", "status": "modeled",
        "inputs": [inputs[0], weight, bias], "outputs": list(outputs),
        "activation": "linear", "mapping_required": False,
        "geometry": {"elements": _product(output_shape), "normalized_size": normalized[0],
                     "epsilon": float(epsilon)},
        "source_nodes": [str(node["id"])],
    }


def _matmul_operation(node: Mapping[str, Any], tensors: Mapping[str, Any],
                      transpose_source: Mapping[str, str] | None = None) -> dict[str, Any]:
    # B-4: aten.matmul/bmm of two ACTIVATIONS. Supports 2-D [M,K]@[K,N] and batched
    # 3-D [Bt,M,K]@[Bt,K,N]; transpose_b when the second operand is [.,N,K] -- either an
    # explicitly [.,N,K]-shaped tensor or a folded transpose (transpose_source).
    inputs, outputs = list(node.get("inputs", [])), node.get("outputs", [])
    if len(inputs) != 2 or len(outputs) != 1:
        raise LoweringError(f"matmul node {node['id']} requires two operands and one output")
    folded_tb = False
    if transpose_source and str(inputs[1]) in transpose_source:
        inputs[1] = transpose_source[str(inputs[1])]
        folded_tb = True   # the folded transpose already presents B as [.,N,K]
    a = _static_shape(tensors[inputs[0]], f"matmul node {node['id']} A")
    b = _static_shape(tensors[inputs[1]], f"matmul node {node['id']} B")
    out = _static_shape(tensors[outputs[0]], f"matmul node {node['id']} out")
    if len(a) not in (2, 3) or len(b) != len(a) or len(out) != len(a):
        raise LoweringError(f"matmul node {node['id']} supports 2-D or batched 3-D matmul only")
    Bt = a[0] if len(a) == 3 else 1
    M, K = a[-2], a[-1]
    # With a folded transpose the source B is already [.,N,K]; otherwise pick the axis
    # that matches K. (A folded transpose whose source is [.,K,N] would be a no-op fold.)
    if folded_tb:
        if b[-1] != K:
            raise LoweringError(f"matmul node {node['id']} folded transpose inner dim disagrees")
        N, tb = b[-2], True
    elif b[-2] == K:
        N, tb = b[-1], False
    elif b[-1] == K:
        N, tb = b[-2], True
    else:
        raise LoweringError(f"matmul node {node['id']} inner dimensions disagree")
    if out[-2] != M or out[-1] != N or (len(out) == 3 and out[0] != Bt):
        raise LoweringError(f"matmul node {node['id']} output shape disagrees")
    return {
        "id": str(node["id"]), "kind": "npusim.matmul", "status": "modeled",
        "inputs": list(inputs), "outputs": list(outputs),
        "activation": "linear", "mapping_required": True,
        "geometry": {"matmul_batch": Bt, "matmul_m": M, "matmul_k": K, "matmul_n": N,
                     "matmul_transpose_b": tb},
        "source_nodes": [str(node["id"])],
    }


def _transpose_operation(node: Mapping[str, Any], tensors: Mapping[str, Any]) -> dict[str, Any]:
    # B-4: aten.transpose(dim0, dim1) / aten.permute swapping exactly two axes -> a physical
    # transpose op (multi-head split/merge). A transpose folded into a matmul's transpose_b
    # never reaches here (handled by the fold pre-pass).
    inputs, outputs = node.get("inputs", []), node.get("outputs", [])
    if len(inputs) != 1 or len(outputs) != 1:
        raise LoweringError(f"transpose node {node['id']} requires one input and one output")
    in_shape = _static_shape(tensors[inputs[0]], f"transpose node {node['id']} input")
    out_shape = _static_shape(tensors[outputs[0]], f"transpose node {node['id']} output")
    args = node.get("arguments")
    op_name = str(node.get("op"))
    if op_name == "aten.permute.default":
        perm = args[1] if isinstance(args, list) and len(args) > 1 else None
        if not isinstance(perm, list):
            raise LoweringError(f"permute node {node['id']} needs an explicit permutation")
        perm = [p % len(in_shape) for p in perm]
        swapped = [i for i, p in enumerate(perm) if p != i]
        if len(swapped) != 2 or perm[swapped[0]] != swapped[1] or perm[swapped[1]] != swapped[0]:
            raise LoweringError(f"permute node {node['id']} must swap exactly two axes")
        a0, a1 = sorted(swapped)
    else:  # aten.transpose.int(dim0, dim1)
        dims = [a for a in (args or []) if isinstance(a, int)]
        if len(dims) < 2:
            raise LoweringError(f"transpose node {node['id']} needs two dims")
        a0, a1 = sorted(d % len(in_shape) for d in dims[:2])
    expected = list(in_shape); expected[a0], expected[a1] = expected[a1], expected[a0]
    if expected != out_shape:
        raise LoweringError(f"transpose node {node['id']} output shape disagrees")
    return {
        "id": str(node["id"]), "kind": "npusim.transpose", "status": "modeled",
        "inputs": [str(inputs[0])], "outputs": [str(outputs[0])],
        "activation": "linear", "mapping_required": False,
        "geometry": {"axis0": a0, "axis1": a1},
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

    # B-4 transpose fold: a transpose node whose sole consumer is a matmul (as its second
    # operand) is absorbed into that matmul's transpose_b. Record {transpose_out: src} and
    # the set of matmul-operand ids that must be rewritten + transpose_b'd.
    transpose_source: dict[str, str] = {}
    folded_transpose_ids: set[str] = set()
    tensor_shapes = {str(t["id"]): t.get("shape") for t in graph["tensors"]}
    for node in nodes:
        op_name = str(node.get("op"))
        if op_name not in _TRANSPOSE_OPS:
            continue
        outs = node.get("outputs", [])
        ins = node.get("inputs", [])
        if not (isinstance(outs, list) and len(outs) == 1 and isinstance(ins, list) and ins):
            continue
        # Only fold a transpose of the LAST TWO axes (the matmul K/N axes). A head-split
        # transpose(0,1) is NOT a transpose_b -- it becomes a real WORKLOAD_TRANSPOSE.
        src_shape = tensor_shapes.get(str(ins[0]))
        nd = len(src_shape) if isinstance(src_shape, list) else 0
        args = node.get("arguments") or []
        dims = sorted(a % nd for a in args if isinstance(a, int)) if nd else []
        last_two = (op_name == "aten.t.default" and nd == 2) or (dims == [nd - 2, nd - 1])
        if not last_two:
            continue
        out_id = str(outs[0])
        cons = consumers.get(out_id, [])
        if (len(cons) == 1 and str(cons[0].get("op")) in _MATMUL_OPS
                and list(cons[0].get("inputs", []))[-1:] == [out_id]):
            transpose_source[out_id] = str(ins[0])
            folded_transpose_ids.add(str(node["id"]))

    # Outputs of real (non-folded) transpose/permute nodes: materialized contiguous in the
    # store, so a following reshape can alias them (head-merge idiom).
    materialized_transpose_ids = {
        str(n["outputs"][0]) for n in nodes
        if str(n.get("op")) in (_TRANSPOSE_OPS | {"aten.permute.default"})
        and str(n["id"]) not in folded_transpose_ids
        and isinstance(n.get("outputs"), list) and len(n["outputs"]) == 1
    }

    operations: list[dict[str, Any]] = []
    lowered_sources: list[str] = []
    tensor_aliases: dict[str, str] = {}
    for node in nodes:
        node_id = str(node["id"])
        if node_id in fused_node_ids or node_id in folded_transpose_ids:
            continue
        operation_name = str(node.get("op"))
        outputs = node.get("outputs", [])
        if not isinstance(outputs, list) or len(outputs) != 1:
            raise LoweringError(f"node {node_id} must have exactly one tensor output")
        if operation_name in _RESHAPE_OPS:
            node_inputs = node.get("inputs", [])
            if not isinstance(node_inputs, list) or len(node_inputs) != 1:
                raise LoweringError(f"reshape node {node_id} requires exactly one tensor input")
            source_id, view_id = str(node_inputs[0]), str(outputs[0])
            # A reshape whose source is a materialized transpose/permute (WORKLOAD_TRANSPOSE
            # physically reorders into contiguous storage) is a contiguous rename in our IR,
            # even though torch marks the transpose view "strided" (the head-merge idiom
            # ctx.transpose(0,1).reshape(T, D)).
            src_is_materialized = source_id in materialized_transpose_ids
            source_shape = ([int(d) for d in tensors[source_id]["shape"]]
                            if src_is_materialized
                            else _whole_contiguous_shape(tensors[source_id], f"reshape node {node_id} input"))
            view_shape = _whole_contiguous_shape(
                tensors[view_id], f"reshape node {node_id} output"
            )
            if tensors[source_id].get("dtype") != tensors[view_id].get("dtype"):
                raise LoweringError(f"reshape node {node_id} changes dtype; not a pure rename")
            if _product(source_shape) != _product(view_shape):
                raise LoweringError(f"reshape node {node_id} does not preserve the element count")
            # Chains (flatten of a view, ...) resolve to the single storage tensor.
            tensor_aliases[view_id] = tensor_aliases.get(source_id, source_id)
            lowered_sources.append(node_id)
            continue
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
        elif operation_name in _POOL_OPS:
            operation = _pool_operation(node, tensors, _POOL_OPS[operation_name])
        elif operation_name in _ELEMENTWISE_OPS:
            operation = _elementwise_operation(node, tensors, _ELEMENTWISE_OPS[operation_name])
        elif operation_name in _CONCAT_OPS:
            operation = _concat_operation(node, tensors)
        elif operation_name in _BATCH_NORM_OPS:
            operation = _batch_norm_operation(node, tensors)
        elif operation_name in _LAYER_NORM_OPS:
            operation = _layer_norm_operation(node, tensors)
        elif operation_name in _MATMUL_OPS:
            operation = _matmul_operation(node, tensors, transpose_source)
        elif operation_name in _TRANSPOSE_OPS or operation_name == "aten.permute.default":
            operation = _transpose_operation(node, tensors)
        elif operation_name in _ACTIVATION_OPS:
            raise LoweringError(
                f"standalone activation node {node_id} cannot yet be scheduled; "
                "fuse it with Linear/Conv2d"
            )
        else:
            raise LoweringError(f"unsupported PyTorch operation {operation_name} at node {node_id}")
        operations.append(operation)
        lowered_sources.extend(source_nodes)

    executable_tensors: list[dict[str, Any]] = []
    for tensor in graph["tensors"]:
        entry = dict(tensor)
        alias_target = tensor_aliases.get(str(tensor["id"]))
        if alias_target is not None:
            entry["alias_of"] = alias_target
        executable_tensors.append(entry)

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
        "tensors": executable_tensors,
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
