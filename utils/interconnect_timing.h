#ifndef __INTERCONNECT_TIMING_H__
#define __INTERCONNECT_TIMING_H__

#include <cstddef>
#include "def.h"

// Timing-only contracts shared by interconnect hierarchy components.
// STORE_AND_FORWARD is modeled as a serialized unicast fabric until a routed
// mesh model is introduced; unsupported topologies must fail at initialization.
bool is_supported_spatial_noc(noc_type_t topology);
bool has_valid_active_shape(unsigned physical_height, unsigned physical_width,
                            unsigned active_height, unsigned active_width);
bool is_valid_memory_line_bits(size_t width_bits);

// Fill/drain latency of a serialized source->link->destination transfer pipelined over
// `transactions` items. Takes the three raw per-item stage costs; the internal pipeline
// stages are derived from them. A single transaction has no overlap (source+link+dest);
// the empty case is 0.
double pipelined_transfer_cycles(unsigned transactions, double source_stage,
                                 double link_stage, double destination_stage);

// Unit static energy is expressed in pJ/cycle.
double static_energy_for_cycles(double unit_energy, double elapsed_cycles);

// All widths below are expressed in bits.  The payload/metadata counters are
// reported independently even for interleaved metadata layouts; `link` is the
// physical serialized transaction count used for timing and energy.
struct datatype_transfer_timing_t {
    size_t source_accesses;
    size_t destination_accesses;
    size_t link_transactions;
    size_t pipeline_transactions;
    size_t payload_link_transactions;
    size_t metadata_link_transactions;
};

datatype_transfer_timing_t datatype_transfer_timing(data_type_t type, size_t elements,
                                                     size_t source_line_bits,
                                                     size_t destination_line_bits,
                                                     size_t link_bits);
struct spatial_noc_cost_t {
    double latency_multiplier;
    double energy_multiplier;
};

spatial_noc_cost_t spatial_noc_cost(noc_type_t topology,
                                    unsigned active_height, unsigned active_width);

// Multi-chip currently models a single serialized bus fabric only.
bool is_supported_multi_chip_nop(noc_type_t topology);

// Balanced 2-input adder-tree reduction over `leaves` active operands.
// depth is the number of sequential levels (parallel additions per level);
// num_additions is the total adder-operation count across the whole tree.
struct adder_tree_reduction_cost_t {
    unsigned num_additions;
    unsigned depth;
};

adder_tree_reduction_cost_t adder_tree_reduction_cost(unsigned leaves);

#endif
