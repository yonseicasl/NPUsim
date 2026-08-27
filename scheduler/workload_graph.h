#ifndef __WORKLOAD_GRAPH_H__
#define __WORKLOAD_GRAPH_H__

#include <cstddef>
#include <string>
#include <vector>

enum workload_operation_kind_t {
    WORKLOAD_LINEAR,
    WORKLOAD_CONV2D,
    WORKLOAD_SOFTMAX,
    WORKLOAD_UNDEFINED
};

struct workload_tensor_t {
    std::string id;
    std::vector<size_t> shape;
    std::string dtype;
    std::string layout;
    std::string kind;
    size_t logical_bytes;

    workload_tensor_t() : logical_bytes(0) {}
};

struct workload_geometry_t {
    unsigned batch;
    unsigned input_features;
    unsigned output_features;
    unsigned input_channels;
    unsigned input_height;
    unsigned input_width;
    unsigned output_channels;
    unsigned output_height;
    unsigned output_width;
    unsigned filter_height;
    unsigned filter_width;
    unsigned stride_height;
    unsigned stride_width;
    unsigned padding_height;
    unsigned padding_width;
    unsigned dilation_height;
    unsigned dilation_width;
    unsigned groups;
    unsigned rows;
    unsigned row_length;

    workload_geometry_t();
};

struct workload_operation_t {
    std::string id;
    workload_operation_kind_t kind;
    std::vector<std::string> inputs;
    std::vector<std::string> outputs;
    std::string activation;
    bool mapping_required;
    workload_geometry_t geometry;
    std::vector<std::string> source_nodes;

    workload_operation_t() : kind(WORKLOAD_UNDEFINED), mapping_required(false) {}
};

// Framework-neutral input consumed by the C++ simulator. The first migration milestone
// still adapts it to Nebula's layer allocation after validation, but operation geometry,
// mapping identity and provenance are owned here rather than inferred from Nebula.
class workload_graph_t {
public:
    workload_graph_t();

    void load(const std::string &m_path);
    std::string legacy_network_config() const;
    std::vector<std::string> mapping_operation_ids() const;

    const workload_tensor_t &tensor(const std::string &m_id) const;
    std::string operation_kind_name(workload_operation_kind_t m_kind) const;
    void validate_operation_geometry(const workload_operation_t &m_operation) const;

    std::string schema_version;
    std::string executable_sha256;
    std::string source_schema_version;
    std::string source_sha256;
    std::string producer_name;
    std::string producer_version;
    std::string model_name;
    std::string model_structure_sha256;
    std::vector<std::string> inputs;
    std::vector<std::string> outputs;
    std::vector<workload_tensor_t> tensors;
    std::vector<workload_operation_t> operations;
};

#endif
