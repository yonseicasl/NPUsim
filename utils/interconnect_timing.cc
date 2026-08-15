#include <algorithm>
#include <cmath>
#include "interconnect_timing.h"
#include "datatype.h"

bool is_supported_spatial_noc(noc_type_t topology) {
    return topology == noc_type_t::BUS ||
           topology == noc_type_t::STORE_AND_FORWARD ||
           topology == noc_type_t::MESH ||
           topology == noc_type_t::CROSSBAR;
}

bool is_supported_multi_chip_nop(noc_type_t topology) {
    return topology == noc_type_t::BUS;
}

adder_tree_reduction_cost_t adder_tree_reduction_cost(unsigned leaves) {
    if(leaves <= 1) return {0, 0};
    const unsigned num_additions = leaves - 1;
    const unsigned depth = static_cast<unsigned>(std::ceil(std::log2(static_cast<double>(leaves))));
    return {num_additions, depth};
}

bool has_valid_active_shape(unsigned physical_height, unsigned physical_width,
                            unsigned active_height, unsigned active_width) {
    return physical_height != 0 && physical_width != 0 &&
           active_height != 0 && active_width != 0 &&
           active_height <= physical_height && active_width <= physical_width;
}

bool is_valid_memory_line_bits(size_t width_bits) {
    return width_bits >= 8 && (width_bits & (width_bits - 1)) == 0;
}

spatial_noc_cost_t spatial_noc_cost(noc_type_t topology,
                                    unsigned active_height, unsigned active_width) {
    spatial_noc_cost_t cost = {1.0, 1.0};
    if(active_height == 0 || active_width == 0 || topology != noc_type_t::MESH) {
        return cost;
    }
    const unsigned endpoints = active_height*active_width;
    const unsigned max_hops = (active_height - 1) + (active_width - 1);
    const double sum_hops = static_cast<double>(active_width)*active_height*(active_height - 1)/2.0 +
                            static_cast<double>(active_height)*active_width*(active_width - 1)/2.0;
    cost.latency_multiplier = static_cast<double>(std::max(1U, max_hops));
    cost.energy_multiplier = std::max(1.0, sum_hops/static_cast<double>(endpoints));
    return cost;
}

double pipelined_transfer_cycles(unsigned transactions, double source_stage,
                                 double link_stage, double destination_stage) {
    if(transactions == 0) return 0.0;
    // A single transaction fills and drains with no overlap: source -> link -> destination.
    if(transactions == 1) return source_stage + link_stage + destination_stage;
    const double second_stage      = std::max(source_stage, link_stage);
    const double other_stage       = std::max(source_stage, std::max(link_stage, destination_stage));
    const double last_before_stage = std::max(link_stage, destination_stage);
    return source_stage + second_stage +
           static_cast<double>(transactions - 2)*other_stage +
           last_before_stage + destination_stage;
}


datatype_transfer_timing_t datatype_transfer_timing(data_type_t type, size_t elements,
                                                     size_t source_line_bits,
                                                     size_t destination_line_bits,
                                                     size_t link_bits) {
    datatype_transfer_timing_t timing = {
        runtime_datatypes().storage_transactions(type, elements, source_line_bits),
        runtime_datatypes().storage_transactions(type, elements, destination_line_bits),
        runtime_datatypes().storage_transactions(type, elements, link_bits),
        0,
        runtime_datatypes().payload_transactions(type, elements, link_bits),
        runtime_datatypes().metadata_transactions(type, elements, link_bits)
    };
    timing.pipeline_transactions = std::max(timing.source_accesses,
                                            std::max(timing.destination_accesses,
                                                     timing.link_transactions));
    return timing;
}
double static_energy_for_cycles(double unit_energy, double elapsed_cycles) {
    return unit_energy * elapsed_cycles;
}
