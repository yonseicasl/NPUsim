#include "workload_graph.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {

using boost::property_tree::ptree;

void fail(const std::string &message) {
    throw std::runtime_error("executable IR: " + message);
}

std::string required_string(const ptree &tree, const std::string &path) {
    const boost::optional<std::string> value = tree.get_optional<std::string>(path);
    if(!value || value->empty()) fail(path + " must be a non-empty string");
    return *value;
}

size_t required_size(const ptree &tree, const std::string &path, bool allow_zero = false) {
    const boost::optional<std::string> text = tree.get_optional<std::string>(path);
    if(!text || text->empty()) fail(path + " must be an integer");
    if(std::find_if(text->begin(), text->end(), [](char value) {
           return !std::isdigit(static_cast<unsigned char>(value));
       }) != text->end()) {
        fail(path + " must be a concrete unsigned integer");
    }
    try {
        size_t parsed = 0;
        const unsigned long long value = std::stoull(*text, &parsed);
        if(parsed != text->size() || (!allow_zero && value == 0) ||
           value > std::numeric_limits<size_t>::max()) {
            fail(path + " is outside the supported range");
        }
        return static_cast<size_t>(value);
    } catch(const std::exception &) {
        fail(path + " is outside the supported range");
    }
    return 0;
}

unsigned required_unsigned(const ptree &tree, const std::string &path, bool allow_zero = false) {
    const size_t value = required_size(tree, path, allow_zero);
    if(value > std::numeric_limits<unsigned>::max()) fail(path + " exceeds unsigned range");
    return static_cast<unsigned>(value);
}

std::vector<std::string> string_array(const ptree &tree, const std::string &path,
                                      bool allow_empty = false) {
    const boost::optional<const ptree&> child = tree.get_child_optional(path);
    if(!child) fail(path + " must be an array");
    std::vector<std::string> result;
    for(ptree::const_iterator it = child->begin(); it != child->end(); ++it) {
        if(!it->first.empty()) fail(path + " must be a JSON array");
        const std::string value = it->second.get_value<std::string>();
        if(value.empty()) fail(path + " contains an empty string");
        result.push_back(value);
    }
    if(result.empty() && !allow_empty) fail(path + " cannot be empty");
    return result;
}

std::vector<size_t> size_array(const ptree &tree, const std::string &path) {
    const boost::optional<const ptree&> child = tree.get_child_optional(path);
    if(!child) fail(path + " must be an array");
    std::vector<size_t> result;
    for(ptree::const_iterator it = child->begin(); it != child->end(); ++it) {
        if(!it->first.empty()) fail(path + " must be a JSON array");
        ptree wrapper;
        wrapper.put("value", it->second.get_value<std::string>());
        result.push_back(required_size(wrapper, "value"));
    }
    return result;
}

size_t checked_multiply(size_t lhs, size_t rhs, const std::string &context) {
    if(lhs != 0 && rhs > std::numeric_limits<size_t>::max()/lhs) fail(context + " overflows");
    return lhs*rhs;
}

size_t dtype_bytes(const std::string &dtype) {
    if(dtype == "bool" || dtype == "uint8" || dtype == "int8") return 1;
    if(dtype == "int16" || dtype == "float16" || dtype == "bfloat16") return 2;
    if(dtype == "int32" || dtype == "float32") return 4;
    if(dtype == "int64" || dtype == "float64") return 8;
    fail("unsupported tensor dtype " + dtype);
    return 0;
}

workload_operation_kind_t parse_operation_kind(const std::string &kind) {
    if(kind == "npusim.linear") return WORKLOAD_LINEAR;
    if(kind == "npusim.conv2d") return WORKLOAD_CONV2D;
    if(kind == "npusim.softmax") return WORKLOAD_SOFTMAX;
    fail("unsupported operation kind " + kind);
    return WORKLOAD_UNDEFINED;
}

bool contains(const std::set<std::string> &values, const std::string &value) {
    return values.find(value) != values.end();
}

} // namespace

workload_geometry_t::workload_geometry_t() :
    batch(0), input_features(0), output_features(0), input_channels(0), input_height(0),
    input_width(0), output_channels(0), output_height(0), output_width(0), filter_height(0),
    filter_width(0), stride_height(0), stride_width(0), padding_height(0), padding_width(0),
    dilation_height(0), dilation_width(0), groups(0), rows(0), row_length(0) {}

workload_graph_t::workload_graph_t() {}

void workload_graph_t::load(const std::string &m_path) {
    ptree root;
    try {
        boost::property_tree::read_json(m_path, root);
    } catch(const boost::property_tree::json_parser::json_parser_error &error) {
        fail(std::string("cannot parse ") + m_path + ": " + error.message());
    }

    schema_version = required_string(root, "schema_version");
    if(schema_version != "npusim.exec.v1") fail("schema_version must be npusim.exec.v1");
    executable_sha256 = root.get<std::string>("executable_sha256", "not-declared");
    producer_name = required_string(root, "producer.name");
    producer_version = required_string(root, "producer.version");
    model_name = required_string(root, "model.name");
    model_structure_sha256 = required_string(root, "model.structure_sha256");
    source_schema_version = required_string(root, "source_graph.schema_version");
    source_sha256 = required_string(root, "source_graph.sha256");
    inputs = string_array(root, "inputs");
    outputs = string_array(root, "outputs");

    tensors.clear();
    operations.clear();
    std::set<std::string> tensor_ids;
    const boost::optional<ptree&> tensor_array = root.get_child_optional("tensors");
    if(!tensor_array || tensor_array->empty()) fail("tensors must be a non-empty array");
    for(ptree::const_iterator it = tensor_array->begin(); it != tensor_array->end(); ++it) {
        if(!it->first.empty()) fail("tensors must be a JSON array");
        const ptree &value = it->second;
        workload_tensor_t tensor_value;
        tensor_value.id = required_string(value, "id");
        if(!tensor_ids.insert(tensor_value.id).second) fail("duplicate tensor id " + tensor_value.id);
        tensor_value.shape = size_array(value, "shape");
        tensor_value.dtype = required_string(value, "dtype");
        tensor_value.layout = value.get<std::string>("layout", "contiguous");
        tensor_value.kind = required_string(value, "kind");
        tensor_value.logical_bytes = required_size(value, "logical_bytes", true);
        size_t elements = 1;
        for(size_t dimension : tensor_value.shape) {
            elements = checked_multiply(elements, dimension, "tensor " + tensor_value.id + " shape");
        }
        const size_t expected_bytes = checked_multiply(elements, dtype_bytes(tensor_value.dtype),
                                                       "tensor " + tensor_value.id + " bytes");
        if(tensor_value.logical_bytes != expected_bytes) {
            fail("tensor " + tensor_value.id + " logical_bytes disagrees with shape/dtype");
        }
        tensors.push_back(tensor_value);
    }

    for(const std::string &id : inputs) if(!contains(tensor_ids, id)) fail("unknown graph input " + id);
    for(const std::string &id : outputs) if(!contains(tensor_ids, id)) fail("unknown graph output " + id);

    std::set<std::string> produced(inputs.begin(), inputs.end());
    for(const workload_tensor_t &tensor_value : tensors) {
        if(tensor_value.kind == "parameter" || tensor_value.kind == "buffer" ||
           tensor_value.kind == "constant") {
            produced.insert(tensor_value.id);
        }
    }
    std::set<std::string> operation_ids;
    const boost::optional<ptree&> operation_array = root.get_child_optional("operations");
    if(!operation_array || operation_array->empty()) fail("operations must be a non-empty array");
    for(ptree::const_iterator it = operation_array->begin(); it != operation_array->end(); ++it) {
        if(!it->first.empty()) fail("operations must be a JSON array");
        const ptree &value = it->second;
        workload_operation_t operation;
        operation.id = required_string(value, "id");
        if(!operation_ids.insert(operation.id).second) fail("duplicate operation id " + operation.id);
        if(required_string(value, "status") != "modeled") fail("operation " + operation.id + " is not modeled");
        operation.kind = parse_operation_kind(required_string(value, "kind"));
        operation.inputs = string_array(value, "inputs", true);
        operation.outputs = string_array(value, "outputs");
        operation.activation = value.get<std::string>("activation", "linear");
        if(operation.activation != "linear" && operation.activation != "relu" &&
           operation.activation != "leaky") {
            fail("operation " + operation.id + " has unsupported activation " + operation.activation);
        }
        operation.mapping_required = value.get<bool>("mapping_required", false);
        operation.source_nodes = string_array(value, "source_nodes");
        for(const std::string &id : operation.inputs) {
            if(!contains(tensor_ids, id)) fail("operation " + operation.id + " uses unknown tensor " + id);
            if(!contains(produced, id)) fail("operation " + operation.id + " is not topologically ordered");
        }
        for(const std::string &id : operation.outputs) {
            if(!contains(tensor_ids, id)) fail("operation " + operation.id + " produces unknown tensor " + id);
            if(contains(produced, id)) fail("operation " + operation.id + " redefines tensor " + id);
        }

        workload_geometry_t &geometry = operation.geometry;
        if(operation.kind == WORKLOAD_LINEAR) {
            geometry.batch = required_unsigned(value, "geometry.batch");
            geometry.input_features = required_unsigned(value, "geometry.input_features");
            geometry.output_features = required_unsigned(value, "geometry.output_features");
            if(!operation.mapping_required) fail("linear operation " + operation.id + " requires a mapping");
        } else if(operation.kind == WORKLOAD_CONV2D) {
            geometry.batch = required_unsigned(value, "geometry.batch");
            geometry.input_channels = required_unsigned(value, "geometry.input_channels");
            geometry.input_height = required_unsigned(value, "geometry.input_height");
            geometry.input_width = required_unsigned(value, "geometry.input_width");
            geometry.output_channels = required_unsigned(value, "geometry.output_channels");
            geometry.output_height = required_unsigned(value, "geometry.output_height");
            geometry.output_width = required_unsigned(value, "geometry.output_width");
            geometry.filter_height = required_unsigned(value, "geometry.filter_height");
            geometry.filter_width = required_unsigned(value, "geometry.filter_width");
            geometry.stride_height = required_unsigned(value, "geometry.stride_height");
            geometry.stride_width = required_unsigned(value, "geometry.stride_width");
            geometry.padding_height = required_unsigned(value, "geometry.padding_height", true);
            geometry.padding_width = required_unsigned(value, "geometry.padding_width", true);
            geometry.dilation_height = required_unsigned(value, "geometry.dilation_height");
            geometry.dilation_width = required_unsigned(value, "geometry.dilation_width");
            geometry.groups = required_unsigned(value, "geometry.groups");
            if(!operation.mapping_required) fail("conv2d operation " + operation.id + " requires a mapping");
        } else {
            geometry.rows = required_unsigned(value, "geometry.rows");
            geometry.row_length = required_unsigned(value, "geometry.row_length");
            if(operation.mapping_required) fail("softmax operation " + operation.id + " cannot use a MAC mapping");
        }
        validate_operation_geometry(operation);
        produced.insert(operation.outputs.begin(), operation.outputs.end());
        operations.push_back(operation);
    }
    for(const std::string &id : outputs) if(!contains(produced, id)) fail("unproduced graph output " + id);
}

const workload_tensor_t &workload_graph_t::tensor(const std::string &m_id) const {
    for(const workload_tensor_t &value : tensors) if(value.id == m_id) return value;
    fail("unknown tensor " + m_id);
    return tensors.front();
}

std::string workload_graph_t::operation_kind_name(workload_operation_kind_t m_kind) const {
    if(m_kind == WORKLOAD_LINEAR) return "linear";
    if(m_kind == WORKLOAD_CONV2D) return "conv2d";
    if(m_kind == WORKLOAD_SOFTMAX) return "softmax";
    return "undefined";
}

std::vector<std::string> workload_graph_t::mapping_operation_ids() const {
    std::vector<std::string> result;
    for(const workload_operation_t &operation : operations) {
        if(operation.mapping_required) result.push_back(operation.id);
    }
    return result;
}

std::string workload_graph_t::legacy_network_config() const {
    if(operations.empty()) fail("cannot adapt an empty workload");
    const workload_operation_t &first = operations.front();
    std::ostringstream config;
    config << "# Generated from " << schema_version << " model " << model_name << "\n";
    config << "[net]\n";
    if(first.kind == WORKLOAD_LINEAR) {
        config << "height=1\nwidth=" << first.geometry.input_features << "\nchannels=1\n";
        config << "batch=" << first.geometry.batch << "\n";
    } else if(first.kind == WORKLOAD_CONV2D) {
        config << "height=" << first.geometry.input_height << "\n";
        config << "width=" << first.geometry.input_width << "\n";
        config << "channels=" << first.geometry.input_channels << "\n";
        config << "batch=" << first.geometry.batch << "\n";
    } else {
        fail("the transitional adapter requires Linear/Conv2d as the first operation");
    }
    config << "num_threads=1\n\n[data]\n";
    config << "test=/dev/null\nlabels=/dev/null\ntop=1\n\n";

    std::string previous_output;
    unsigned network_batch = first.geometry.batch;
    for(size_t index = 0; index < operations.size(); ++index) {
        const workload_operation_t &operation = operations[index];
        if(operation.inputs.empty()) fail("operation " + operation.id + " has no data input");
        if(index > 0 && operation.inputs.front() != previous_output) {
            fail("transitional Nebula adapter supports only a linear tensor chain; operation " +
                 operation.id + " does not consume the previous output");
        }
        if(operation.kind == WORKLOAD_LINEAR) {
            if(operation.geometry.batch != network_batch) fail("linear batch changes inside one workload");
            config << "[connected]\n";
            config << "output=" << operation.geometry.output_features << "\n";
            config << "activation=" << operation.activation << "\n\n";
        } else if(operation.kind == WORKLOAD_CONV2D) {
            if(operation.geometry.batch != network_batch) fail("conv2d batch changes inside one workload");
            if(operation.geometry.stride_height != operation.geometry.stride_width) {
                fail("legacy timing adapter cannot represent asymmetric conv2d stride");
            }
            if(operation.geometry.dilation_height != 1 || operation.geometry.dilation_width != 1) {
                fail("legacy timing adapter cannot represent dilated conv2d");
            }
            config << "[convolutional]\n";
            config << "filters=" << operation.geometry.output_channels << "\n";
            config << "filter_height=" << operation.geometry.filter_height << "\n";
            config << "filter_width=" << operation.geometry.filter_width << "\n";
            config << "stride=" << operation.geometry.stride_height << "\n";
            config << "padding_height=" << operation.geometry.padding_height << "\n";
            config << "padding_width=" << operation.geometry.padding_width << "\n";
            config << "group=" << operation.geometry.groups << "\n";
            config << "activation=" << operation.activation << "\n\n";
        } else {
            if(index + 1 != operations.size()) {
                fail("transitional Nebula adapter supports softmax only as the final operation");
            }
            if(operation.geometry.rows % network_batch != 0) {
                fail("softmax rows must be divisible by the fixed network batch");
            }
            config << "[softmax]\n";
            config << "groups=" << operation.geometry.rows/network_batch << "\n\n";
        }
        previous_output = operation.outputs.front();
    }
    return config.str();
}
