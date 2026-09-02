#ifndef __WORKLOAD_GRAPH_H__
#define __WORKLOAD_GRAPH_H__

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <vector>

enum workload_operation_kind_t {
    WORKLOAD_LINEAR,
    WORKLOAD_CONV2D,
    WORKLOAD_SOFTMAX,
    WORKLOAD_POOL2D,
    WORKLOAD_ELEMENTWISE,
    WORKLOAD_CONCAT,
    WORKLOAD_BATCH_NORM,
    WORKLOAD_UNDEFINED
};

struct workload_tensor_t {
    std::string id;
    std::vector<size_t> shape;
    std::string dtype;
    std::string layout;
    std::string kind;
    // Non-empty for elided reshape views (flatten/view/reshape): this tensor is a
    // rename of alias_of's storage, never a separate allocation or transfer.
    std::string alias_of;
    size_t logical_bytes;

    workload_tensor_t() : logical_bytes(0) {}
    size_t elements() const;
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
    unsigned kernel_height;
    unsigned kernel_width;
    unsigned axis;
    size_t elements;
    bool count_include_pad;
    double epsilon;
    std::string mode;
    std::string elementwise_operator;

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
    // Resolves alias views to the tensor that owns the storage; identity for
    // canonical tensors. Use for residency, byte, and datatype classification;
    // use tensor() where the view's shape is the semantic contract.
    const workload_tensor_t &storage_tensor(const std::string &m_id) const;
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


enum workload_tensor_residency_t {
    WORKLOAD_RESIDENCY_DRAM,
    WORKLOAD_RESIDENCY_GLB
};

struct workload_residency_plan_t {
    std::vector<workload_tensor_residency_t> inputs;
    std::vector<bool> retain_inputs;
    bool retain_output;
    size_t occupied_before;
    size_t occupied_after;
    size_t capacity;
    std::vector<std::string> released;

    workload_residency_plan_t() : retain_output(false), occupied_before(0),
                                  occupied_after(0), capacity(0) {}
};

// Capacity-aware tensor liveness for executable-IR timing. Parameters/constants remain
// DRAM-backed; graph/activation inputs and outputs stay in the GLB while they have future
// consumers and capacity permits. A last-use input buffer may be reused by the operation's
// output, avoiding artificial materialization without inventing eviction traffic.
class workload_lifetime_t {
public:
    workload_lifetime_t(const workload_graph_t &m_graph, size_t m_capacity,
                        const std::map<std::string, size_t> &m_runtime_bytes);

    workload_residency_plan_t plan(size_t m_operation_index) const;
    void commit(size_t m_operation_index, workload_residency_plan_t *m_plan);
    bool resident(const std::string &m_tensor_id) const;
    size_t occupied_bytes() const { return occupied; }

private:
    // All bookkeeping is keyed by storage tensors: alias views resolve through
    // canonical() so a consumer reading through a reshape keeps the producer's
    // buffer alive instead of forcing a fictitious DRAM round trip.
    const std::string &canonical(const std::string &m_tensor_id) const;

    const workload_graph_t &graph;
    size_t capacity;
    size_t occupied;
    std::map<std::string, size_t> runtime_bytes;
    std::map<std::string, size_t> remaining_consumers;
    std::set<std::string> glb_resident;
    std::set<std::string> graph_output_storage;
};

#endif
