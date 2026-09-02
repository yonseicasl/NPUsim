#include "workload_graph.h"

#include <limits>
#include <stdexcept>

namespace {

void geometry_fail(const std::string &operation_id, const std::string &message) {
    throw std::runtime_error("executable IR: operation " + operation_id + " " + message);
}

size_t shape_product(const std::vector<size_t> &shape, size_t begin, size_t end,
                     const std::string &operation_id) {
    size_t result = 1;
    for(size_t index = begin; index < end; ++index) {
        if(result != 0 && shape[index] > std::numeric_limits<size_t>::max()/result) {
            geometry_fail(operation_id, "tensor shape overflows");
        }
        result *= shape[index];
    }
    return result;
}

} // namespace

void workload_graph_t::validate_operation_geometry(const workload_operation_t &operation) const {
    if(operation.outputs.size() != 1) {
        geometry_fail(operation.id, "must have exactly one output tensor");
    }
    const workload_tensor_t &output = tensor(operation.outputs[0]);

    if(operation.kind == WORKLOAD_LINEAR) {
        if(operation.inputs.size() != 2 && operation.inputs.size() != 3) {
            geometry_fail(operation.id, "linear needs input, weight, and optional bias");
        }
        const workload_tensor_t &input = tensor(operation.inputs[0]);
        const workload_tensor_t &weight = tensor(operation.inputs[1]);
        if(input.shape.empty() || output.shape.empty() || weight.shape.size() != 2) {
            geometry_fail(operation.id, "linear tensor ranks are invalid");
        }
        const size_t batch = shape_product(input.shape, 0, input.shape.size() - 1, operation.id);
        if(batch != operation.geometry.batch ||
           input.shape.back() != operation.geometry.input_features ||
           weight.shape[0] != operation.geometry.output_features ||
           weight.shape[1] != operation.geometry.input_features ||
           shape_product(output.shape, 0, output.shape.size() - 1, operation.id) != batch ||
           output.shape.back() != operation.geometry.output_features) {
            geometry_fail(operation.id, "linear geometry disagrees with tensor shapes");
        }
        if(operation.inputs.size() == 3) {
            const workload_tensor_t &bias = tensor(operation.inputs[2]);
            if(bias.shape.size() != 1 || bias.shape[0] != operation.geometry.output_features) {
                geometry_fail(operation.id, "linear bias shape is invalid");
            }
        }
        return;
    }

    if(operation.kind == WORKLOAD_CONV2D) {
        if(operation.inputs.size() != 2 && operation.inputs.size() != 3) {
            geometry_fail(operation.id, "conv2d needs input, weight, and optional bias");
        }
        const workload_tensor_t &input = tensor(operation.inputs[0]);
        const workload_tensor_t &weight = tensor(operation.inputs[1]);
        if(input.shape.size() != 4 || output.shape.size() != 4 || weight.shape.size() != 4) {
            geometry_fail(operation.id, "conv2d currently requires NCHW/OIHW rank-4 tensors");
        }
        const workload_geometry_t &g = operation.geometry;
        if(g.input_channels % g.groups || g.output_channels % g.groups ||
           input.shape[0] != g.batch || input.shape[1] != g.input_channels ||
           input.shape[2] != g.input_height || input.shape[3] != g.input_width ||
           weight.shape[0] != g.output_channels || weight.shape[1] != g.input_channels/g.groups ||
           weight.shape[2] != g.filter_height || weight.shape[3] != g.filter_width ||
           output.shape[0] != g.batch || output.shape[1] != g.output_channels ||
           output.shape[2] != g.output_height || output.shape[3] != g.output_width) {
            geometry_fail(operation.id, "conv2d geometry disagrees with tensor shapes");
        }
        const size_t effective_h = static_cast<size_t>(g.dilation_height)*(g.filter_height - 1) + 1;
        const size_t effective_w = static_cast<size_t>(g.dilation_width)*(g.filter_width - 1) + 1;
        const size_t padded_h = static_cast<size_t>(g.input_height) + 2ULL*g.padding_height;
        const size_t padded_w = static_cast<size_t>(g.input_width) + 2ULL*g.padding_width;
        if(padded_h < effective_h || padded_w < effective_w ||
           (padded_h - effective_h)/g.stride_height + 1 != g.output_height ||
           (padded_w - effective_w)/g.stride_width + 1 != g.output_width) {
            geometry_fail(operation.id, "conv2d output extent is invalid");
        }
        if(operation.inputs.size() == 3) {
            const workload_tensor_t &bias = tensor(operation.inputs[2]);
            if(bias.shape.size() != 1 || bias.shape[0] != g.output_channels) {
                geometry_fail(operation.id, "conv2d bias shape is invalid");
            }
        }
        return;
    }

    if(operation.kind == WORKLOAD_SOFTMAX) {
        if(operation.inputs.size() != 1) {
            geometry_fail(operation.id, "softmax needs exactly one input tensor");
        }
        const workload_tensor_t &input = tensor(operation.inputs[0]);
        if(input.shape.empty() || input.shape != output.shape ||
           input.shape.back() != operation.geometry.row_length ||
           shape_product(input.shape, 0, input.shape.size() - 1, operation.id) !=
               operation.geometry.rows) {
            geometry_fail(operation.id, "softmax geometry disagrees with tensor shapes");
        }
        return;
    }

    if(operation.kind == WORKLOAD_POOL2D) {
        if(operation.inputs.size() != 1) geometry_fail(operation.id, "pool2d needs one input");
        const workload_tensor_t &input = tensor(operation.inputs[0]);
        const workload_geometry_t &g = operation.geometry;
        if((g.mode != "max" && g.mode != "average") || input.shape.size() != 4 ||
           output.shape.size() != 4 || input.shape[0] != g.batch ||
           input.shape[1] != g.input_channels || input.shape[2] != g.input_height ||
           input.shape[3] != g.input_width || output.shape[0] != g.batch ||
           output.shape[1] != g.output_channels || output.shape[2] != g.output_height ||
           output.shape[3] != g.output_width) {
            geometry_fail(operation.id, "pool2d geometry disagrees with tensor shapes");
        }
        const size_t effective_h = static_cast<size_t>(g.dilation_height)*(g.kernel_height - 1) + 1;
        const size_t effective_w = static_cast<size_t>(g.dilation_width)*(g.kernel_width - 1) + 1;
        const size_t padded_h = static_cast<size_t>(g.input_height) + 2ULL*g.padding_height;
        const size_t padded_w = static_cast<size_t>(g.input_width) + 2ULL*g.padding_width;
        if(padded_h < effective_h || padded_w < effective_w ||
           (padded_h - effective_h)/g.stride_height + 1 != g.output_height ||
           (padded_w - effective_w)/g.stride_width + 1 != g.output_width) {
            geometry_fail(operation.id, "pool2d output extent is invalid");
        }
        return;
    }

    if(operation.kind == WORKLOAD_ELEMENTWISE) {
        if(operation.inputs.size() != 2 ||
           (operation.geometry.elementwise_operator != "add" &&
            operation.geometry.elementwise_operator != "multiply")) {
            geometry_fail(operation.id, "elementwise contract is invalid");
        }
        const workload_tensor_t &lhs = tensor(operation.inputs[0]);
        const workload_tensor_t &rhs = tensor(operation.inputs[1]);
        if(lhs.shape != rhs.shape || lhs.shape != output.shape ||
           output.elements() != operation.geometry.elements) {
            geometry_fail(operation.id, "elementwise shapes disagree (broadcasting is unsupported)");
        }
        return;
    }

    if(operation.kind == WORKLOAD_CONCAT) {
        if(operation.inputs.size() < 2 || operation.geometry.axis >= output.shape.size()) {
            geometry_fail(operation.id, "concat needs two inputs and a valid axis");
        }
        std::vector<size_t> expected = tensor(operation.inputs[0]).shape;
        if(expected.size() != output.shape.size()) geometry_fail(operation.id, "concat rank mismatch");
        expected[operation.geometry.axis] = 0;
        for(const std::string &input_id : operation.inputs) {
            const workload_tensor_t &input = tensor(input_id);
            if(input.shape.size() != output.shape.size()) geometry_fail(operation.id, "concat rank mismatch");
            for(size_t axis = 0; axis < output.shape.size(); ++axis) {
                if(axis != operation.geometry.axis && input.shape[axis] != output.shape[axis]) {
                    geometry_fail(operation.id, "concat non-axis dimensions disagree");
                }
            }
            expected[operation.geometry.axis] += input.shape[operation.geometry.axis];
        }
        if(expected != output.shape || output.elements() != operation.geometry.elements) {
            geometry_fail(operation.id, "concat output shape disagrees with inputs");
        }
        return;
    }

    if(operation.kind == WORKLOAD_BATCH_NORM) {
        if(operation.inputs.size() < 3 || operation.inputs.size() > 5) {
            geometry_fail(operation.id,
                          "batch_norm needs activation, running statistics, and optional affine tensors");
        }
        const workload_tensor_t &input = tensor(operation.inputs[0]);
        if(input.shape.size() < 2 || input.shape != output.shape ||
           output.elements() != operation.geometry.elements ||
           input.shape[1] != operation.geometry.output_channels) {
            geometry_fail(operation.id, "batch_norm geometry disagrees with tensor shapes");
        }
        for(size_t index = 1; index < operation.inputs.size(); ++index) {
            const workload_tensor_t &parameter = tensor(operation.inputs[index]);
            if(parameter.shape.size() != 1 ||
               parameter.shape[0] != operation.geometry.output_channels) {
                geometry_fail(operation.id, "batch_norm parameter/stat shape is invalid");
            }
        }
        return;
    }

    geometry_fail(operation.id, "has no geometry validator");
}
