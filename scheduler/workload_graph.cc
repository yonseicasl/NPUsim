#include "workload_graph.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cmath>
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
    if(kind == "npusim.pool2d") return WORKLOAD_POOL2D;
    if(kind == "npusim.elementwise") return WORKLOAD_ELEMENTWISE;
    if(kind == "npusim.concat") return WORKLOAD_CONCAT;
    if(kind == "npusim.batch_norm") return WORKLOAD_BATCH_NORM;
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
    dilation_height(0), dilation_width(0), groups(0), rows(0), row_length(0),
    kernel_height(0), kernel_width(0), axis(0), elements(0), count_include_pad(true),
    epsilon(0.0), mode(), elementwise_operator() {}

size_t workload_tensor_t::elements() const {
    size_t result = 1;
    for(size_t dimension : shape) result *= dimension;
    return result;
}

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
        tensor_value.alias_of = value.get<std::string>("alias_of", "");
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

    for(const workload_tensor_t &tensor_value : tensors) {
        if(tensor_value.alias_of.empty()) continue;
        if(!contains(tensor_ids, tensor_value.alias_of)) {
            fail("tensor " + tensor_value.id + " aliases unknown tensor " + tensor_value.alias_of);
        }
        const workload_tensor_t &target = tensor(tensor_value.alias_of);
        if(!target.alias_of.empty()) {
            fail("tensor " + tensor_value.id + " aliases " + target.id + ", which is itself an alias");
        }
        if(target.dtype != tensor_value.dtype || target.logical_bytes != tensor_value.logical_bytes) {
            fail("tensor " + tensor_value.id + " alias does not preserve its storage dtype and bytes");
        }
        if(tensor_value.kind != "activation") {
            fail("tensor " + tensor_value.id + " alias must be an activation view");
        }
    }

    for(const std::string &id : inputs) if(!contains(tensor_ids, id)) fail("unknown graph input " + id);
    for(const std::string &id : outputs) if(!contains(tensor_ids, id)) fail("unknown graph output " + id);
    for(const std::string &id : inputs) {
        if(!tensor(id).alias_of.empty()) fail("graph input " + id + " must be a storage tensor, not an alias");
    }

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
            const workload_tensor_t &view = tensor(id);
            const std::string &storage = view.alias_of.empty() ? view.id : view.alias_of;
            if(!contains(produced, storage)) fail("operation " + operation.id + " is not topologically ordered");
        }
        for(const std::string &id : operation.outputs) {
            if(!contains(tensor_ids, id)) fail("operation " + operation.id + " produces unknown tensor " + id);
            if(!tensor(id).alias_of.empty()) fail("operation " + operation.id + " cannot produce alias tensor " + id);
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
        } else if(operation.kind == WORKLOAD_SOFTMAX) {
            geometry.rows = required_unsigned(value, "geometry.rows");
            geometry.row_length = required_unsigned(value, "geometry.row_length");
            if(operation.mapping_required) fail("softmax operation " + operation.id + " cannot use a MAC mapping");
        } else if(operation.kind == WORKLOAD_POOL2D) {
            geometry.mode = required_string(value, "geometry.mode");
            geometry.batch = required_unsigned(value, "geometry.batch");
            geometry.input_channels = geometry.output_channels = required_unsigned(value, "geometry.channels");
            geometry.input_height = required_unsigned(value, "geometry.input_height");
            geometry.input_width = required_unsigned(value, "geometry.input_width");
            geometry.output_height = required_unsigned(value, "geometry.output_height");
            geometry.output_width = required_unsigned(value, "geometry.output_width");
            geometry.kernel_height = required_unsigned(value, "geometry.kernel_height");
            geometry.kernel_width = required_unsigned(value, "geometry.kernel_width");
            geometry.stride_height = required_unsigned(value, "geometry.stride_height");
            geometry.stride_width = required_unsigned(value, "geometry.stride_width");
            geometry.padding_height = required_unsigned(value, "geometry.padding_height", true);
            geometry.padding_width = required_unsigned(value, "geometry.padding_width", true);
            geometry.dilation_height = required_unsigned(value, "geometry.dilation_height");
            geometry.dilation_width = required_unsigned(value, "geometry.dilation_width");
            geometry.count_include_pad = value.get<bool>("geometry.count_include_pad", true);
            if(operation.mapping_required) fail("pool operation " + operation.id + " cannot use a MAC mapping");
        } else if(operation.kind == WORKLOAD_ELEMENTWISE) {
            geometry.elementwise_operator = required_string(value, "geometry.operator");
            geometry.elements = required_size(value, "geometry.elements");
            if(operation.mapping_required) fail("elementwise operation " + operation.id + " cannot use a MAC mapping");
        } else if(operation.kind == WORKLOAD_CONCAT) {
            geometry.axis = required_unsigned(value, "geometry.axis", true);
            geometry.elements = required_size(value, "geometry.elements");
            if(operation.mapping_required) fail("concat operation " + operation.id + " cannot use a MAC mapping");
        } else if(operation.kind == WORKLOAD_BATCH_NORM) {
            geometry.elements = required_size(value, "geometry.elements");
            geometry.output_channels = required_unsigned(value, "geometry.channels");
            geometry.epsilon = value.get<double>("geometry.epsilon", 0.0);
            if(!std::isfinite(geometry.epsilon) || geometry.epsilon <= 0.0) {
                fail("batch_norm operation " + operation.id + " requires a positive finite epsilon");
            }
            if(operation.mapping_required) fail("batch_norm operation " + operation.id + " cannot use a MAC mapping");
        }
        validate_operation_geometry(operation);
        produced.insert(operation.outputs.begin(), operation.outputs.end());
        operations.push_back(operation);
    }
    for(const std::string &id : outputs) {
        const workload_tensor_t &view = tensor(id);
        const std::string &storage = view.alias_of.empty() ? view.id : view.alias_of;
        if(!contains(produced, storage)) fail("unproduced graph output " + id);
    }
}

const workload_tensor_t &workload_graph_t::tensor(const std::string &m_id) const {
    for(const workload_tensor_t &value : tensors) if(value.id == m_id) return value;
    fail("unknown tensor " + m_id);
    return tensors.front();
}

const workload_tensor_t &workload_graph_t::storage_tensor(const std::string &m_id) const {
    const workload_tensor_t &view = tensor(m_id);
    return view.alias_of.empty() ? view : tensor(view.alias_of);
}

std::string workload_graph_t::operation_kind_name(workload_operation_kind_t m_kind) const {
    if(m_kind == WORKLOAD_LINEAR) return "linear";
    if(m_kind == WORKLOAD_CONV2D) return "conv2d";
    if(m_kind == WORKLOAD_SOFTMAX) return "softmax";
    if(m_kind == WORKLOAD_POOL2D) return "pool2d";
    if(m_kind == WORKLOAD_ELEMENTWISE) return "elementwise";
    if(m_kind == WORKLOAD_CONCAT) return "concat";
    if(m_kind == WORKLOAD_BATCH_NORM) return "batch_norm";
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
    if(first.inputs.empty()) fail("first operation has no data input");
    const workload_tensor_t &first_input = tensor(first.inputs.front());
    const size_t batch = first_input.shape.empty() ? 1 : first_input.shape.front();
    if(batch == 0 || batch > std::numeric_limits<unsigned>::max()) fail("network batch is invalid");
    const size_t per_batch = first_input.elements()/batch;
    std::ostringstream config;
    config << "# Generated from " << schema_version << " model " << model_name << "\n";
    config << "[net]\n";
    if(first_input.shape.size() == 4) {
        config << "height=" << first_input.shape[2] << "\nwidth=" << first_input.shape[3]
               << "\nchannels=" << first_input.shape[1] << "\n";
    } else {
        config << "height=1\nwidth=" << per_batch << "\nchannels=1\n";
    }
    config << "batch=" << batch << "\nnum_threads=1\n\n[data]\n";
    config << "test=/dev/null\nlabels=/dev/null\ntop=1\n\n";

    for(const workload_operation_t &operation : operations) {
        if(operation.inputs.empty()) fail("operation " + operation.id + " has no data input");
        if(operation.kind == WORKLOAD_LINEAR) {
            config << "[connected]\noutput=" << operation.geometry.output_features
                   << "\nactivation=" << operation.activation << "\n\n";
        } else if(operation.kind == WORKLOAD_CONV2D) {
            config << "[convolutional]\nfilters=" << operation.geometry.output_channels << "\n";
            config << "filter_height=" << operation.geometry.filter_height << "\n";
            config << "filter_width=" << operation.geometry.filter_width << "\n";
            config << "stride=" << operation.geometry.stride_height << "\n";
            config << "padding_height=" << operation.geometry.padding_height << "\n";
            config << "padding_width=" << operation.geometry.padding_width << "\n";
            config << "group=" << operation.geometry.groups << "\n";
            config << "activation=" << operation.activation << "\n\n";
        } else if(operation.kind == WORKLOAD_SOFTMAX) {
            config << "[softmax]\ngroups=1\n\n";
        } else {
            // Timing-only shape placeholder. The framework-neutral descriptor below is
            // authoritative; after Nebula allocates this harmless identity layer, npu_t
            // overwrites every geometry field from the executable IR.
            config << "[maxpool]\nsize=1\nstride=1\npadding=0\n\n";
        }
    }
    return config.str();
}


workload_lifetime_t::workload_lifetime_t(
    const workload_graph_t &m_graph, size_t m_capacity,
    const std::map<std::string, size_t> &m_runtime_bytes) :
    graph(m_graph), capacity(m_capacity), occupied(0), runtime_bytes(m_runtime_bytes) {
    for(const workload_tensor_t &tensor_value : graph.tensors) {
        remaining_consumers[tensor_value.id] = 0;
        if(runtime_bytes.find(tensor_value.id) == runtime_bytes.end()) {
            runtime_bytes[tensor_value.id] = tensor_value.logical_bytes;
        }
    }
    for(const workload_operation_t &operation : graph.operations) {
        for(const std::string &input : operation.inputs) ++remaining_consumers[canonical(input)];
    }
    for(const std::string &output : graph.outputs) graph_output_storage.insert(canonical(output));
}

const std::string &workload_lifetime_t::canonical(const std::string &m_tensor_id) const {
    const workload_tensor_t &view = graph.tensor(m_tensor_id);
    return view.alias_of.empty() ? view.id : view.alias_of;
}

bool workload_lifetime_t::resident(const std::string &m_tensor_id) const {
    return glb_resident.find(canonical(m_tensor_id)) != glb_resident.end();
}

workload_residency_plan_t workload_lifetime_t::plan(size_t m_operation_index) const {
    if(m_operation_index >= graph.operations.size()) fail("lifetime operation index is out of range");
    const workload_operation_t &operation = graph.operations[m_operation_index];
    workload_residency_plan_t result;
    result.capacity = capacity;
    result.occupied_before = occupied;
    result.retain_inputs.assign(operation.inputs.size(), false);
    std::map<std::string, size_t> uses_in_operation;
    size_t reclaimable = 0;
    for(const std::string &input : operation.inputs) {
        result.inputs.push_back(resident(input) ? WORKLOAD_RESIDENCY_GLB : WORKLOAD_RESIDENCY_DRAM);
        ++uses_in_operation[canonical(input)];
    }
    for(const std::pair<const std::string, size_t> &entry : uses_in_operation) {
        const std::map<std::string, size_t>::const_iterator remaining =
            remaining_consumers.find(entry.first);
        if(resident(entry.first) && remaining != remaining_consumers.end() &&
           remaining->second <= entry.second) {
            reclaimable += runtime_bytes.at(entry.first);
        }
    }

    size_t planned_occupied = occupied - std::min(occupied, reclaimable);
    const std::string &output = operation.outputs.front();
    // Storage-resolved membership: a graph output declared through an alias view
    // still pins its producing storage tensor as an external output.
    const bool graph_output = graph_output_storage.find(output) != graph_output_storage.end();
    const size_t output_bytes = runtime_bytes.at(output);
    const size_t future_uses = remaining_consumers.at(output);
    result.retain_output = !graph_output && future_uses > 0 &&
                           output_bytes <= capacity - std::min(capacity, planned_occupied);
    if(result.retain_output) planned_occupied += output_bytes;

    // A DRAM-backed graph/activation input that will be used again is already staged
    // through the GLB for this operation. Pin it after the access when capacity remains;
    // parameters and constants stay DRAM-backed. Output retention has priority because
    // it feeds the next producer-consumer edge and avoids immediate materialization.
    std::set<std::string> newly_retained;
    for(size_t index = 0; index < operation.inputs.size(); ++index) {
        const std::string &input_id = canonical(operation.inputs[index]);
        if(result.inputs[index] == WORKLOAD_RESIDENCY_GLB ||
           newly_retained.find(input_id) != newly_retained.end()) {
            continue;
        }
        const workload_tensor_t &input_tensor = graph.tensor(input_id);
        if(input_tensor.kind == "parameter" || input_tensor.kind == "buffer" ||
           input_tensor.kind == "constant") {
            continue;
        }
        const size_t uses_now = uses_in_operation.at(input_id);
        const size_t future_input_uses = remaining_consumers.at(input_id) - uses_now;
        const size_t input_bytes = runtime_bytes.at(input_id);
        if(future_input_uses > 0 &&
           input_bytes <= capacity - std::min(capacity, planned_occupied)) {
            result.retain_inputs[index] = true;
            newly_retained.insert(input_id);
            planned_occupied += input_bytes;
        }
    }
    result.occupied_after = planned_occupied;
    return result;
}

void workload_lifetime_t::commit(size_t m_operation_index, workload_residency_plan_t *m_plan) {
    if(m_plan == NULL || m_operation_index >= graph.operations.size() ||
       m_plan->retain_inputs.size() != graph.operations[m_operation_index].inputs.size()) {
        fail("invalid lifetime commit");
    }
    const workload_operation_t &operation = graph.operations[m_operation_index];
    for(const std::string &raw_input : operation.inputs) {
        const std::string &input = canonical(raw_input);
        std::map<std::string, size_t>::iterator remaining = remaining_consumers.find(input);
        if(remaining == remaining_consumers.end() || remaining->second == 0) {
            fail("tensor consumer underflow " + input);
        }
        --remaining->second;
        if(remaining->second == 0 && glb_resident.erase(input) != 0) {
            occupied -= runtime_bytes.at(input);
            m_plan->released.push_back(input);
        }
    }
    const std::string &output = operation.outputs.front();
    if(m_plan->retain_output) {
        if(!glb_resident.insert(output).second) fail("tensor is already resident " + output);
        occupied += runtime_bytes.at(output);
    }
    for(size_t index = 0; index < operation.inputs.size(); ++index) {
        if(!m_plan->retain_inputs[index]) continue;
        const std::string &input = canonical(operation.inputs[index]);
        if(remaining_consumers.at(input) == 0) fail("cannot retain a dead input tensor " + input);
        if(glb_resident.insert(input).second) occupied += runtime_bytes.at(input);
    }
    if(occupied > capacity || occupied != m_plan->occupied_after) {
        fail("GLB lifetime accounting exceeded or disagreed with its plan");
    }
    m_plan->occupied_after = occupied;
}
