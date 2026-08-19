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

nop_delivery_cost_t nop_delivery_cost(noc_type_t topology, unsigned grid_width,
                                      unsigned active_chips, bool multicast) {
    nop_delivery_cost_t cost = {0.0, 0.0, 0.0};
    if(active_chips == 0) return cost;
    if(topology != noc_type_t::MESH) {
        // BUS: one shared medium, no hops. A multicast transmission is seen by every
        // receiver at once; distinct chunks serialize on the medium.
        cost.bottleneck_link_tiles = multicast ? 1.0 : static_cast<double>(active_chips);
        cost.total_link_traversals = cost.bottleneck_link_tiles;
        cost.fill_hops = 0.0;
        return cost;
    }
    // MESH: dimension-order routes from the ingress at (0,0). Count how many links each
    // route traverses, and (for multicast) how many DISTINCT links the union of the routes
    // covers -- a multicast tree forwards one copy per link it uses, not one per receiver.
    const unsigned columns = (grid_width == 0) ? 1 : grid_width;
    unsigned unicast_traversals = 0, max_hops = 1;
    for(unsigned chip = 0; chip < active_chips; ++chip) {
        const unsigned hops = nop_route_hops(chip, columns);
        unicast_traversals += hops;
        max_hops = std::max(max_hops, hops);
    }
    if(!multicast) {
        cost.bottleneck_link_tiles = static_cast<double>(active_chips);
        cost.total_link_traversals = static_cast<double>(unicast_traversals);
        cost.fill_hops = static_cast<double>(max_hops - 1);
        return cost;
    }
    // Union of the dimension-order routes: the ingress link, every horizontal link along
    // row 0 up to the rightmost used column, and every vertical link in each used column.
    unsigned tree_links = 1;   // the ingress link into (0,0)
    unsigned widest_column = 0, tallest_row = 0;
    std::vector<unsigned> column_height(columns, 0);
    for(unsigned chip = 0; chip < active_chips; ++chip) {
        const unsigned x = chip % columns;
        const unsigned y = chip / columns;
        widest_column = std::max(widest_column, x);
        tallest_row = std::max(tallest_row, y);
        column_height[x] = std::max(column_height[x], y);
    }
    tree_links += widest_column;                       // row-0 horizontal links
    for(unsigned x = 0; x < columns; ++x) tree_links += column_height[x];
    (void)tallest_row;
    cost.bottleneck_link_tiles = 1.0;
    cost.total_link_traversals = static_cast<double>(tree_links);
    cost.fill_hops = static_cast<double>(max_hops - 1);
    return cost;
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

namespace {

size_t ceil_div_bits(size_t value, size_t divisor) {
    return (divisor == 0) ? 0 : (value + divisor - 1)/divisor;
}

// Derive the tail packet's count so the per-packet counts sum EXACTLY to the aggregate the
// caller charges cycles from. `per_packet` may need lowering too when the aggregate is
// smaller than the full-packet extrapolation (only reachable when the widths are not
// integer multiples of each other, or when a metadata stream is counted separately).
void reconcile_packet_counts(size_t packets, size_t aggregate,
                             size_t *per_packet, size_t *tail) {
    if(packets == 0) { *per_packet = 0; *tail = 0; return; }
    const size_t full = (packets - 1)*(*per_packet);
    if(aggregate >= full) {
        *tail = aggregate - full;
        return;
    }
    *per_packet = aggregate/packets;
    *tail = aggregate - (packets - 1)*(*per_packet);
}

} // namespace

size_t transfer_packet_groups_t::source_transactions() const {
    return (packets == 0) ? 0 : (packets - 1)*source_per_packet + source_tail;
}

size_t transfer_packet_groups_t::link_transactions() const {
    return (packets == 0) ? 0 : (packets - 1)*link_per_packet + link_tail;
}

size_t transfer_packet_groups_t::destination_transactions() const {
    return (packets == 0) ? 0 : (packets - 1)*destination_per_packet + destination_tail;
}

transfer_packet_groups_t transfer_packet_groups_t::without_source() const {
    transfer_packet_groups_t dropped = *this;
    dropped.source_per_packet = 0;
    dropped.source_tail = 0;
    return dropped;
}

transfer_packet_groups_t transfer_packet_groups_t::without_destination() const {
    transfer_packet_groups_t dropped = *this;
    dropped.destination_per_packet = 0;
    dropped.destination_tail = 0;
    return dropped;
}

transfer_packet_groups_t transfer_packet_groups_t::with_destination_total(size_t transactions) const {
    transfer_packet_groups_t replaced = *this;
    if(packets == 0) return replaced;
    replaced.destination_per_packet = ceil_div_bits(transactions, packets);
    reconcile_packet_counts(packets, transactions,
                            &replaced.destination_per_packet, &replaced.destination_tail);
    return replaced;
}

transfer_packet_groups_t transfer_packet_groups(size_t payload_bits,
                                                size_t source_line_bits, size_t link_bits,
                                                size_t destination_line_bits,
                                                size_t source_transactions,
                                                size_t link_transactions,
                                                size_t destination_transactions) {
    transfer_packet_groups_t groups = {0, 0, 0, 0, 0, 0, 0};
    const size_t packet_bits = std::max(source_line_bits,
                                        std::max(link_bits, destination_line_bits));
    if(payload_bits == 0 || packet_bits == 0) return groups;

    groups.packets = ceil_div_bits(payload_bits, packet_bits);
    const size_t tail_bits = payload_bits - (groups.packets - 1)*packet_bits;

    groups.source_per_packet = ceil_div_bits(packet_bits, source_line_bits);
    groups.link_per_packet = ceil_div_bits(packet_bits, link_bits);
    groups.destination_per_packet = ceil_div_bits(packet_bits, destination_line_bits);
    groups.source_tail = ceil_div_bits(tail_bits, source_line_bits);
    groups.link_tail = ceil_div_bits(tail_bits, link_bits);
    groups.destination_tail = ceil_div_bits(tail_bits, destination_line_bits);

    reconcile_packet_counts(groups.packets, source_transactions,
                            &groups.source_per_packet, &groups.source_tail);
    reconcile_packet_counts(groups.packets, link_transactions,
                            &groups.link_per_packet, &groups.link_tail);
    reconcile_packet_counts(groups.packets, destination_transactions,
                            &groups.destination_per_packet, &groups.destination_tail);
    return groups;
}

namespace {

// One packet through the chain of PRESENT stages, advancing each stage's resource
// availability. Returns the packet's completion time.
double advance_packet(const size_t *count, const double *cost, double *available) {
    bool have_upstream = false;
    size_t upstream_count = 0;
    double upstream_start = 0.0, upstream_end = 0.0, upstream_cost = 0.0;
    double completion = 0.0;
    for(unsigned stage = 0; stage < 3; ++stage) {
        if(count[stage] == 0) continue;
        double start = available[stage];
        if(have_upstream) {
            // A gather boundary (fewer transactions here than upstream) packs several
            // upstream transactions into one of ours, so the whole upstream packet must
            // land first. A fan-out boundary only needs the first upstream transaction.
            const double dependency = (count[stage] < upstream_count)
                ? upstream_end : upstream_start + upstream_cost;
            start = std::max(start, dependency);
        }
        double finish = start + static_cast<double>(count[stage])*cost[stage];
        // Our LAST transaction cannot precede the last upstream transaction feeding it.
        if(have_upstream) finish = std::max(finish, upstream_end + cost[stage]);
        available[stage] = finish;
        upstream_count = count[stage];
        upstream_start = start;
        upstream_end = finish;
        upstream_cost = cost[stage];
        have_upstream = true;
        completion = std::max(completion, finish);
    }
    return completion;
}

} // namespace

double pipelined_transfer_cycles(const transfer_packet_groups_t &groups,
                                 double source_stage, double link_stage,
                                 double destination_stage) {
    if(groups.packets == 0) return 0.0;
    const double cost[3] = {source_stage, link_stage, destination_stage};
    const size_t full_count[3] = {groups.source_per_packet, groups.link_per_packet,
                                  groups.destination_per_packet};
    const size_t tail_count[3] = {groups.source_tail, groups.link_tail,
                                  groups.destination_tail};
    double available[3] = {0.0, 0.0, 0.0};
    double completion = 0.0;

    // Steady-state period of the identical full packets: the busiest stage's per-packet
    // occupancy. The recurrence above is a (max,+) system, so over identical packets it
    // becomes periodic with exactly this period once each stage has been visited; simulate
    // the transient explicitly and extrapolate the rest instead of looping over a packet
    // count that scales with the payload.
    double period = 0.0;
    for(unsigned stage = 0; stage < 3; ++stage) {
        period = std::max(period, static_cast<double>(full_count[stage])*cost[stage]);
    }

    // The transient lasts at most one packet per stage; this bound only guards against a
    // floating-point comparison never reporting periodicity, so that the loop can never
    // become O(payload).
    const size_t transient_bound = 8;
    const double tolerance = 1e-9*std::max(1.0, period);
    const size_t full_packets = groups.packets - 1;
    size_t simulated = 0;
    while(simulated < full_packets) {
        double before[3] = {available[0], available[1], available[2]};
        completion = std::max(completion, advance_packet(full_count, cost, available));
        ++simulated;
        if(simulated < 2 || simulated >= full_packets) continue;
        // Only stages that actually carry work in a full packet advance; an absent stage
        // (a dropped source or destination) stays put and must not veto periodicity.
        bool periodic = true;
        for(unsigned stage = 0; stage < 3; ++stage) {
            if(full_count[stage] == 0) continue;
            if(std::fabs((available[stage] - before[stage]) - period) > tolerance) periodic = false;
        }
        if(periodic || simulated >= transient_bound) {
            const double remaining = static_cast<double>(full_packets - simulated);
            for(unsigned stage = 0; stage < 3; ++stage) {
                if(full_count[stage] == 0) continue;
                available[stage] += remaining*period;
            }
            completion += remaining*period;
            simulated = full_packets;
        }
    }
    completion = std::max(completion, advance_packet(tail_count, cost, available));
    return completion;
}

systolic_pipeline_cost_t systolic_pipeline_cost(unsigned active_height, unsigned active_width) {
    systolic_pipeline_cost_t cost = {0.0, 0.0};
    if(active_height == 0 || active_width == 0) return cost;
    cost.skew_fill_hops = static_cast<double>(active_height - 1) +
                          static_cast<double>(active_width - 1);
    cost.drain_hops = static_cast<double>(active_height - 1);
    return cost;
}

double pipelined_transfer_cycles(unsigned transactions, double source_stage,
                                 double link_stage, double destination_stage) {
    // CE5: equal-width endpoints are one transaction per stage per packet.
    transfer_packet_groups_t groups = {transactions, 1, 1, 1, 1, 1, 1};
    if(transactions == 0) groups.packets = 0;
    return pipelined_transfer_cycles(groups, source_stage, link_stage, destination_stage);
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
        runtime_datatypes().metadata_transactions(type, elements, link_bits),
        {0, 0, 0, 0, 0, 0, 0}
    };
    timing.pipeline_transactions = std::max(timing.source_accesses,
                                            std::max(timing.destination_accesses,
                                                     timing.link_transactions));
    // L2: the packet boundaries come from the payload and the three widths, and are
    // reconciled against the counts above so the makespan and the charged work agree.
    timing.groups = transfer_packet_groups(
        runtime_datatypes().storage_bits(type, elements),
        source_line_bits, link_bits, destination_line_bits,
        timing.source_accesses, timing.link_transactions, timing.destination_accesses);
    return timing;
}
double temporal_pipeline_run_cycles(unsigned tiles, const std::vector<double> &stage_totals) {
    if(stage_totals.empty()) return 0.0;
    if(tiles <= 1) {
        return *std::max_element(stage_totals.begin(), stage_totals.end());
    }
    double fill = 0.0, peak_per_tile = 0.0;
    for(double total : stage_totals) {
        const double per_tile = total/static_cast<double>(tiles);
        fill += per_tile;
        peak_per_tile = std::max(peak_per_tile, per_tile);
    }
    return fill + static_cast<double>(tiles - 1)*peak_per_tile;
}

double combine_datatype_cycles(const std::vector<double> &per_type, bool serialized_types) {
    double combined = 0.0;
    for(unsigned i = 0; i < per_type.size(); ++i) {
        combined = serialized_types ? combined + per_type[i]
                                    : std::max(combined, per_type[i]);
    }
    return combined;
}

double entity_combined_cycles(const std::vector<std::vector<double>> &per_entity_type,
                              bool serialized_types) {
    double combined = 0.0;
    for(unsigned entity = 0; entity < per_entity_type.size(); ++entity) {
        combined = std::max(combined, combine_datatype_cycles(per_entity_type[entity],
                                                             serialized_types));
    }
    return combined;
}

double entity_combined_cycles(const std::vector<std::vector<double>> &a,
                              const std::vector<std::vector<double>> &b,
                              bool serialized_types) {
    double combined = 0.0;
    const size_t entities = std::max(a.size(), b.size());
    for(size_t entity = 0; entity < entities; ++entity) {
        const std::vector<double> empty;
        const std::vector<double> &left = (entity < a.size()) ? a[entity] : empty;
        const std::vector<double> &right = (entity < b.size()) ? b[entity] : empty;
        std::vector<double> per_type(std::max(left.size(), right.size()), 0.0);
        for(size_t t = 0; t < per_type.size(); ++t) {
            // L1: the two sides of the SAME datatype partition serialize on that
            // partition's port, so they add here -- before the datatype combination.
            if(t < left.size()) per_type[t] += left[t];
            if(t < right.size()) per_type[t] += right[t];
        }
        combined = std::max(combined, combine_datatype_cycles(per_type, serialized_types));
    }
    return combined;
}

double pipeline_timeline_cycles(const std::vector<std::vector<double>> &stage_tile_costs,
                               const std::vector<unsigned> &boundary_depths,
                               std::vector<double> *stall) {
    const size_t stages = stage_tile_costs.size();
    if(stall) stall->assign(stages, 0.0);
    if(stages == 0) return 0.0;
    const size_t tiles = stage_tile_costs[0].size();
    for(size_t s = 1; s < stages; ++s) {
        if(stage_tile_costs[s].size() != tiles) {
            std::cerr << "Error: pipeline timeline stages must share one tile clock"
                      << std::endl;
            exit(1);
        }
    }
    if(boundary_depths.size() + 1 != stages) {
        std::cerr << "Error: pipeline timeline needs one buffer depth per stage boundary"
                  << std::endl;
        exit(1);
    }
    if(tiles == 0) return 0.0;

    // Only the last `depth` completions of each stage are ever read (the back-pressure
    // term), so keep a small ring per stage instead of a full history: the tile clock is
    // the GLB repetition count and can be large, and nothing here needs to look further
    // back than the deepest buffer.
    size_t max_depth = 1;
    for(size_t b = 0; b < boundary_depths.size(); ++b) {
        max_depth = std::max(max_depth, static_cast<size_t>(std::max(1U, boundary_depths[b])));
    }
    const size_t ring = max_depth + 2;
    std::vector<std::vector<double>> finish(stages, std::vector<double>(ring, 0.0));
    double makespan = 0.0;
    for(size_t k = 0; k < tiles; ++k) {
        for(size_t s = 0; s < stages; ++s) {
            // The stage is one serial resource: its previous tile must be done.
            double start = (k > 0) ? finish[s][(k - 1) % ring] : 0.0;
            // Its input tile must have arrived (written earlier in THIS k iteration).
            if(s > 0) start = std::max(start, finish[s - 1][k % ring]);
            // And a downstream buffer slot must be free (written `depth` iterations ago,
            // which the ring still holds because ring > max depth).
            double back_pressure = 0.0;
            if(s + 1 < stages) {
                const size_t depth = std::max(1U, boundary_depths[s]);
                if(k >= depth) back_pressure = finish[s + 1][(k - depth) % ring];
            }
            const double unblocked = start;
            start = std::max(start, back_pressure);
            if(stall && start > unblocked) stall->at(s) += start - unblocked;
            finish[s][k % ring] = start + stage_tile_costs[s][k];
            makespan = std::max(makespan, finish[s][k % ring]);
        }
    }
    return makespan;
}

dram_row_activation_cost_t dram_row_activations(size_t stream_bytes, size_t row_buffer_bytes,
                                                unsigned num_banks) {
    dram_row_activation_cost_t cost = {0, 0};
    if(row_buffer_bytes == 0 || stream_bytes == 0) return cost;
    const unsigned banks = (num_banks == 0) ? 1 : num_banks;
    cost.activations = (stream_bytes + row_buffer_bytes - 1)/row_buffer_bytes;
    cost.busiest_bank = (cost.activations + banks - 1)/banks;
    return cost;
}

double static_energy_for_cycles(double unit_energy, double elapsed_cycles) {
    return unit_energy * elapsed_cycles;
}
