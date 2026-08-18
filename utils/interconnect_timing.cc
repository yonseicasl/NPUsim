#include <algorithm>
#include <cmath>
#include <iostream>
#include "interconnect_timing.h"
#include "datatype.h"

bool is_supported_spatial_noc(noc_type_t topology) {
    return topology == noc_type_t::BUS ||
           topology == noc_type_t::STORE_AND_FORWARD ||
           topology == noc_type_t::MESH ||
           topology == noc_type_t::CROSSBAR;
}

bool is_supported_multi_chip_nop(noc_type_t topology) {
    return topology == noc_type_t::BUS || topology == noc_type_t::MESH;
}

unsigned nop_route_hops(unsigned chip_index, unsigned grid_width) {
    if(grid_width == 0) return 1;
    const unsigned x = chip_index % grid_width;
    const unsigned y = chip_index / grid_width;
    return std::max(1U, x + y);
}

unsigned derived_link_bitwidth(const char *component, double bandwidth, double frequency) {
    if(frequency <= 0.0) return 0;
    const double bits = 8.0*bandwidth/frequency;
    const unsigned truncated = static_cast<unsigned>(bits);
    if(bits > 0.0 && bits != static_cast<double>(truncated)) {
        std::cerr << "Warning: " << component << " link bitwidth derived from bandwidth/frequency"
                  << " truncates " << bits << " -> " << truncated
                  << " bits; set 'bitwidth' explicitly" << std::endl;
    }
    return truncated;
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
    spatial_noc_cost_t cost = {0.0, 1.0};
    if(active_height == 0 || active_width == 0 || topology != noc_type_t::MESH) {
        return cost;
    }
    const unsigned endpoints = active_height*active_width;
    const unsigned max_hops = (active_height - 1) + (active_width - 1);
    const double sum_hops = static_cast<double>(active_width)*active_height*(active_height - 1)/2.0 +
                            static_cast<double>(active_height)*active_width*(active_width - 1)/2.0;
    cost.latency_fill_hops = static_cast<double>(std::max(1U, max_hops)) - 1.0;
    cost.energy_multiplier = std::max(1.0, sum_hops/static_cast<double>(endpoints));
    return cost;
}

double pipelined_transfer_cycles(unsigned transactions, double source_stage,
                                 double link_stage, double destination_stage) {
    if(transactions == 0) return 0.0;
    // CE5: ideal 3-stage pipeline makespan -- one fill through every stage, then
    // the slowest stage issues the remaining transactions back to back.
    const double slowest = std::max(source_stage, std::max(link_stage, destination_stage));
    return source_stage + link_stage + destination_stage +
           static_cast<double>(transactions - 1)*slowest;
}

// CE5: stage transaction counts may differ when the source line, link width, and
// destination line differ -- each stage works its OWN count; the makespan is one
// fill through every stage plus the bottleneck stage's remaining service time.
double pipelined_transfer_cycles(size_t source_transactions, double source_stage,
                                 size_t link_transactions, double link_stage,
                                 size_t destination_transactions, double destination_stage) {
    if(source_transactions == 0 && link_transactions == 0 && destination_transactions == 0) return 0.0;
    double fill = 0.0, bottleneck = 0.0;
    if(source_transactions > 0) {
        fill += source_stage;
        bottleneck = std::max(bottleneck, static_cast<double>(source_transactions - 1)*source_stage);
    }
    if(link_transactions > 0) {
        fill += link_stage;
        bottleneck = std::max(bottleneck, static_cast<double>(link_transactions - 1)*link_stage);
    }
    if(destination_transactions > 0) {
        fill += destination_stage;
        bottleneck = std::max(bottleneck, static_cast<double>(destination_transactions - 1)*destination_stage);
    }
    return fill + bottleneck;
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
