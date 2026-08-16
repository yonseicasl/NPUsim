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
// Pipelined-hop contract (SP1): a stream of T transactions to the farthest active
// destination completes in (T + latency_fill_hops)*noc_cycle -- hops pipeline, so
// the route depth is a one-time fill, not a per-transaction multiplier. Energy stays
// per-transaction x avg Manhattan hops (every hop transfer burns link energy).
struct spatial_noc_cost_t {
    double latency_fill_hops;   // max(1, max Manhattan hops) - 1; 0 for single-hop fabrics
    double energy_multiplier;   // average Manhattan hops per delivered transaction
};

spatial_noc_cost_t spatial_noc_cost(noc_type_t topology,
                                    unsigned active_height, unsigned active_width);

// Multi-chip NoP: a serialized bus fabric, or a routed-unicast mesh where each chip's
// stream traverses its Manhattan distance from the package ingress at grid (0,0).
bool is_supported_multi_chip_nop(noc_type_t topology);

// Manhattan hop count from the package ingress (0,0) to chip `chip_index`, laid out
// row-major on a grid with `grid_width` columns. At least one link is traversed.
unsigned nop_route_hops(unsigned chip_index, unsigned grid_width);

// B12: derive a link bit-width (bits) from bandwidth (GB/s) and frequency (MHz),
// warning on stderr when the result truncates fractionally -- silent truncation
// biases every ceil(bits/bitwidth) transaction count. Returns 0 when frequency <= 0.
unsigned derived_link_bitwidth(const char *component, double bandwidth, double frequency);

// Balanced 2-input adder-tree reduction over `leaves` active operands.
// depth is the number of sequential levels (parallel additions per level);
// num_additions is the total adder-operation count across the whole tree.
struct adder_tree_reduction_cost_t {
    unsigned num_additions;
    unsigned depth;
};

adder_tree_reduction_cost_t adder_tree_reduction_cost(unsigned leaves);

#endif
