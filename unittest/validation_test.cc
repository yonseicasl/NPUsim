#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "datatype.h"
#include "mapping_table.h"
#include "pe_lane.h"
#include "interconnect_timing.h"
#include "energy_units.h"

namespace {

void fail(const std::string &message) {
    std::cerr << "validation test failed: " << message << std::endl;
    std::exit(1);
}

// L5 test helpers: per-tile cost vectors on the shared tile clock.
std::vector<double> uniform_tiles(size_t tiles, double cost) {
    return std::vector<double>(tiles, cost);
}

// A stage that serves only `served` of `tiles` tiles (an off-chip stage refetching at a
// lower rate than the on-chip stages consume), carrying its whole cost on those tiles.
std::vector<double> sparse_tiles(size_t tiles, double cost, size_t served) {
    std::vector<double> costs(tiles, 0.0);
    for(size_t k = 0; k < tiles; ++k) {
        if(k == 0 || (k*served)/tiles != ((k - 1)*served)/tiles) costs[k] = cost;
    }
    return costs;
}

// A run whose LAST tile is cheaper than the rest.
std::vector<double> tail_tiles(size_t tiles, double cost, double tail_cost) {
    std::vector<double> costs(tiles, cost);
    if(tiles > 0) costs[tiles - 1] = tail_cost;
    return costs;
}

// Two-stage back-pressure attribution: the reported stall must be the time each stage sat
// blocked by a full downstream buffer, which is what tells a reader WHERE it stalled.
bool stall_is(size_t producer_tiles, double producer_cost, size_t consumer_tiles,
              double consumer_cost, unsigned depth, double expected_producer_stall,
              double expected_consumer_stall) {
    std::vector<double> stall;
    pipeline_timeline_cycles({uniform_tiles(producer_tiles, producer_cost),
                              uniform_tiles(consumer_tiles, consumer_cost)}, {depth}, &stall);
    return stall.size() == 2 &&
           std::fabs(stall[0] - expected_producer_stall) < 1e-9 &&
           std::fabs(stall[1] - expected_consumer_stall) < 1e-9;
}

// L2 test helpers: state a width-conversion case the way a real config produces it --
// payload bits plus the source-line / link / destination-line widths -- so the packet
// boundaries under test are the physical ones, not an arbitrary count triple.
double packet_makespan(size_t payload_bits, size_t source_line_bits, size_t link_bits,
                       size_t destination_line_bits, size_t source_transactions,
                       size_t link_transactions, size_t destination_transactions,
                       double source_stage, double link_stage, double destination_stage) {
    return pipelined_transfer_cycles(
        transfer_packet_groups(payload_bits, source_line_bits, link_bits, destination_line_bits,
                               source_transactions, link_transactions, destination_transactions),
        source_stage, link_stage, destination_stage);
}

// Brute-force reference for the packet recurrence: simulate EVERY packet, with no
// steady-state extrapolation. pipelined_transfer_cycles() short-circuits the identical
// full packets once the (max,+) recurrence becomes periodic, so that shortcut has to be
// pinned against the unabridged simulation -- including large packet counts, payloads that
// are not a whole number of packets, and dropped stages (which do not advance and must not
// veto periodicity).
double brute_force_packet_makespan(const transfer_packet_groups_t &groups, double source_stage,
                                   double link_stage, double destination_stage) {
    if(groups.packets == 0) return 0.0;
    const double cost[3] = {source_stage, link_stage, destination_stage};
    double available[3] = {0.0, 0.0, 0.0};
    double completion = 0.0;
    for(size_t packet = 0; packet < groups.packets; ++packet) {
        const bool tail = (packet + 1 == groups.packets);
        const size_t count[3] = {
            tail ? groups.source_tail : groups.source_per_packet,
            tail ? groups.link_tail : groups.link_per_packet,
            tail ? groups.destination_tail : groups.destination_per_packet};
        bool have_upstream = false;
        size_t upstream_count = 0;
        double upstream_start = 0.0, upstream_end = 0.0, upstream_cost = 0.0;
        for(unsigned stage = 0; stage < 3; ++stage) {
            if(count[stage] == 0) continue;
            double start = available[stage];
            if(have_upstream) {
                start = std::max(start, (count[stage] < upstream_count)
                                        ? upstream_end : upstream_start + upstream_cost);
            }
            double finish = start + static_cast<double>(count[stage])*cost[stage];
            if(have_upstream) finish = std::max(finish, upstream_end + cost[stage]);
            available[stage] = finish;
            upstream_count = count[stage];
            upstream_start = start;
            upstream_end = finish;
            upstream_cost = cost[stage];
            have_upstream = true;
            completion = std::max(completion, finish);
        }
    }
    return completion;
}

// True when the extrapolating implementation matches the unabridged simulation for the
// whole transfer and for each dropped-stage variant.
bool packet_extrapolation_exact(size_t payload_bits, size_t source_line_bits, size_t link_bits,
                                size_t destination_line_bits, double source_stage,
                                double link_stage, double destination_stage) {
    const size_t source_transactions = (payload_bits + source_line_bits - 1)/source_line_bits;
    const size_t link_transactions = (payload_bits + link_bits - 1)/link_bits;
    const size_t destination_transactions =
        (payload_bits + destination_line_bits - 1)/destination_line_bits;
    const transfer_packet_groups_t groups = transfer_packet_groups(
        payload_bits, source_line_bits, link_bits, destination_line_bits,
        source_transactions, link_transactions, destination_transactions);
    const transfer_packet_groups_t variants[3] = {
        groups, groups.without_destination(), groups.without_source()};
    for(unsigned i = 0; i < 3; ++i) {
        const double fast = pipelined_transfer_cycles(variants[i], source_stage, link_stage,
                                                      destination_stage);
        const double reference = brute_force_packet_makespan(variants[i], source_stage,
                                                             link_stage, destination_stage);
        if(std::fabs(fast - reference) > 1e-6) return false;
    }
    return true;
}

// The per-packet counts must sum back to the aggregate counts the caller charges cycles
// and energy from; otherwise the makespan describes different work than was billed.
bool packet_counts_consistent(size_t payload_bits, size_t source_line_bits, size_t link_bits,
                              size_t destination_line_bits, size_t source_transactions,
                              size_t link_transactions, size_t destination_transactions) {
    const transfer_packet_groups_t groups = transfer_packet_groups(
        payload_bits, source_line_bits, link_bits, destination_line_bits,
        source_transactions, link_transactions, destination_transactions);
    return groups.source_transactions() == source_transactions &&
           groups.link_transactions() == link_transactions &&
           groups.destination_transactions() == destination_transactions;
}

void validate_parser_contract() {
    section_config_t section("test");
    if(!section.add_setting("vector", "1:2:3") ||
       !section.add_setting("broadcast", "7") ||
       !section.add_setting("negative", "-1") ||
       section.add_setting("vector", "4:5:6")) {
        fail("duplicate-key contract");
    }

    std::vector<unsigned> values(3, 0);
    if(!section.get_vector_setting("vector", &values) ||
       values[0] != 1 || values[1] != 2 || values[2] != 3) {
        fail("exact vector parsing");
    }
    if(!section.get_vector_setting("broadcast", &values) ||
       values[0] != 7 || values[1] != 7 || values[2] != 7) {
        fail("scalar vector broadcast");
    }

    // A string setting is the WHOLE value, spaces included -- stream extraction used to stop at
    // the first space and then fail the eof() check, dropping a multi-word provenance string.
    section.add_setting("text", "a b c");
    std::string text;
    if(!section.get_setting("text", &text) || text != "a b c") {
        fail("string setting keeps its spaces");
    }

    // The REJECTION cases -- a negative value for an unsigned setting, a short or long
    // per-datatype vector, a fractional value for an integer setting, a bool given anything but
    // 0/1 -- are no longer observable in-process: they used to return false (silently leaving the
    // caller's default, which is the defect) and now abort with a message naming the section, key
    // and value. They are verified end to end by validation/knobs KN11, which runs the simulator
    // on a config carrying each bad value and requires it to fail.
}

void validate_legacy_glb_mapping_contract() {
    section_config_t section("mapping");
    section.add_setting("glb", "2,3,1,1,1,1,1,0,0,2,1");
    mapping_table_t table(section);
    if(table.calculate_active_component(component_type_t::GLOBAL_BUFFER) != 3) {
        fail("legacy GLB temporal repetition mapping");
    }
}

void validate_input_halo_contract() {
    section_config_t conv3("mapping");
    conv3.add_setting("pe", "16,1,1,1,4,1,3,0,0,1,1");
    conv3.add_setting("pe_x", "1,1,13,1,1,1,1,0,0,1,1");
    conv3.add_setting("pe_y", "4,1,1,1,1,3,1,0,0,1,1");
    conv3.add_setting("glb", "6,4,1,13,1,1,1,0,0,1,1");
    conv3.add_setting("dram", "1,1,1,1,64,1,1,0,0,1,1");
    const input_halo_reuse_t q_only = mapping_table_t(conv3).input_halo_reuse();
    if(!q_only.active || q_only.replicated_elements != 599040 ||
       q_only.unique_elements != 230400 || q_only.working_set_elements != 11520) {
        fail("Q-only input halo union/working-set contract");
    }

    section_config_t conv2("mapping");
    conv2.add_setting("pe", "16,1,1,1,1,1,5,0,0,1,1");
    conv2.add_setting("pe_x", "1,1,14,1,1,1,1,0,0,1,1");
    conv2.add_setting("pe_y", "1,1,1,1,2,5,1,0,0,1,1");
    conv2.add_setting("glb", "8,4,2,27,1,1,1,0,0,1,1");
    conv2.add_setting("dram", "2,1,1,1,48,1,1,0,0,2,1");
    const input_halo_reuse_t pq = mapping_table_t(conv2).input_halo_reuse();
    if(!pq.active || pq.replicated_elements != 1866240 ||
       pq.unique_elements != 380928 || pq.working_set_elements != 18624) {
        fail("P/Q input halo union/working-set contract");
    }
}

void validate_runtime_datatypes() {
    struct dense_format_case_t {
        const char *name;
        size_t bits;
    };
    const dense_format_case_t dense_formats[] = {
        {"fp32", 32}, {"fp16", 16}, {"bf16", 16}, {"int8", 8},
        {"int4", 4}, {"int2", 2}, {"uint8", 8}
    };
    for(const dense_format_case_t &format : dense_formats) {
        section_config_t format_section("accelerator");
        format_section.add_setting("data_format", format.name);
        runtime_datatypes().configure(format_section);
        if(runtime_datatypes().describe(data_type_t::INPUT) != format.name ||
           runtime_datatypes().describe(data_type_t::WEIGHT) != format.name ||
           runtime_datatypes().describe(data_type_t::OUTPUT) != format.name ||
           runtime_datatypes().payload_bits(data_type_t::INPUT, 3) != 3*format.bits ||
           runtime_datatypes().metadata_bits(data_type_t::INPUT, 3) != 0 ||
           runtime_datatypes().storage_transactions(data_type_t::INPUT, 3, 8) !=
               (3*format.bits + 7)/8) {
            fail(std::string("dense datatype format contract: ") + format.name);
        }
    }

    section_config_t section("accelerator");
    section.add_setting("data_format", "bf16");
    section.add_setting("weight_format", "mxfp4");
    section.add_setting("output_format", "int4");
    section.add_setting("accumulator_format", "fp32");
    runtime_datatypes().configure(section);

    if(runtime_datatypes().describe(data_type_t::INPUT) != "bf16" ||
       runtime_datatypes().storage_bytes(data_type_t::INPUT, 32) != 64 ||
       runtime_datatypes().describe(data_type_t::WEIGHT) != "mxfp4_b32_e8m0" ||
       runtime_datatypes().payload_bits(data_type_t::WEIGHT, 32) != 128 ||
       runtime_datatypes().metadata_bits(data_type_t::WEIGHT, 32) != 8 ||
       runtime_datatypes().storage_bytes(data_type_t::WEIGHT, 32) != 17 ||
       runtime_datatypes().storage_bytes(data_type_t::WEIGHT, 3) != 3 ||
       runtime_datatypes().payload_transactions(data_type_t::WEIGHT, 33, 16) != 9 ||
       runtime_datatypes().metadata_transactions(data_type_t::WEIGHT, 33, 16) != 1 ||
       runtime_datatypes().storage_transactions(data_type_t::WEIGHT, 33, 32) != 6 ||
       runtime_datatypes().storage_transactions(data_type_t::WEIGHT, 33, 16) != 10 ||
       runtime_datatypes().accumulator_format().name != "fp32" ||
       runtime_datatypes().describe(data_type_t::OUTPUT) != "int4" ||
       runtime_datatypes().storage_bytes(data_type_t::OUTPUT, 3) != 2 ||
       runtime_datatypes().storage_transactions(data_type_t::OUTPUT, 3, 8) != 2) {
        fail("runtime datatype storage accounting");
    }

    section_config_t interleaved("accelerator");
    interleaved.add_setting("weight_format", "mxfp4");
    interleaved.add_setting("mxfp_metadata_layout", "interleaved");
    runtime_datatypes().configure(interleaved);
    if(runtime_datatypes().metadata_layout() != mxfp_metadata_layout_t::INTERLEAVED ||
       runtime_datatypes().storage_bytes(data_type_t::WEIGHT, 3) != 3 ||
       runtime_datatypes().storage_transactions(data_type_t::WEIGHT, 33, 32) != 5) {
        fail("interleaved MXFP metadata layout accounting");
    }

    runtime_datatypes().configure(section);
    const datatype_transfer_timing_t mxfp_transfer = datatype_transfer_timing(
        data_type_t::WEIGHT, 33, 16, 32, 64);
    const datatype_transfer_timing_t int4_transfer = datatype_transfer_timing(
        data_type_t::OUTPUT, 3, 8, 16, 32);
    if(mxfp_transfer.source_accesses != 10 ||
       mxfp_transfer.destination_accesses != 6 ||
       mxfp_transfer.link_transactions != 4 ||
       mxfp_transfer.pipeline_transactions != 10 ||
       mxfp_transfer.payload_link_transactions != 3 ||
       mxfp_transfer.metadata_link_transactions != 1 ||
       int4_transfer.source_accesses != 2 ||
       int4_transfer.destination_accesses != 1 ||
       int4_transfer.link_transactions != 1 ||
       int4_transfer.pipeline_transactions != 2 ||
       int4_transfer.payload_link_transactions != 1 ||
       int4_transfer.metadata_link_transactions != 0) {
        fail("global buffer bit-width transaction contract");
    }

}
void validate_hierarchy_datatype_transactions() {
    section_config_t section("accelerator");
    section.add_setting("input_format", "int2");
    section.add_setting("weight_format", "mxfp8");
    section.add_setting("output_format", "fp32");
    section.add_setting("mxfp_metadata_layout", "separate");
    runtime_datatypes().configure(section);

    const datatype_transfer_timing_t int2 = datatype_transfer_timing(
        data_type_t::INPUT, 17, 8, 32, 16);
    const datatype_transfer_timing_t mxfp8 = datatype_transfer_timing(
        data_type_t::WEIGHT, 33, 64, 32, 16);
    const datatype_transfer_timing_t fp32 = datatype_transfer_timing(
        data_type_t::OUTPUT, 3, 64, 128, 32);

    if(int2.source_accesses != 5 || int2.destination_accesses != 2 ||
       int2.link_transactions != 3 || int2.pipeline_transactions != 5 ||
       int2.payload_link_transactions != 3 || int2.metadata_link_transactions != 0 ||
       mxfp8.source_accesses != 6 || mxfp8.destination_accesses != 10 ||
       mxfp8.link_transactions != 18 || mxfp8.pipeline_transactions != 18 ||
       mxfp8.payload_link_transactions != 17 || mxfp8.metadata_link_transactions != 1 ||
       fp32.source_accesses != 2 || fp32.destination_accesses != 1 ||
       fp32.link_transactions != 3 || fp32.pipeline_transactions != 3 ||
       !is_valid_memory_line_bits(8) || !is_valid_memory_line_bits(256) ||
       is_valid_memory_line_bits(0) || is_valid_memory_line_bits(12)) {
        fail("hierarchy-wide datatype transaction contract");
    }
}

void validate_pe_lane_contract() {
    const mac_lane_state_t full = calculate_mac_lane_state(8, 8, 64);
    if(!full.valid || full.scalar_capacity != 64 ||
       full.active_accumulator_units != 8 || full.final_accumulator_lanes != 8 ||
       full.utilization != 1.0) {
        fail("full PE lane contract");
    }

    const mac_lane_state_t partial = calculate_mac_lane_state(8, 8, 20);
    if(!partial.valid || partial.scalar_capacity != 64 ||
       partial.active_accumulator_units != 3 || partial.final_accumulator_lanes != 4 ||
       partial.utilization != 0.3125) {
        fail("partial PE lane contract");
    }

    if(calculate_mac_lane_state(8, 8, 0).valid ||
       calculate_mac_lane_state(8, 8, 65).valid ||
       calculate_mac_lane_state(0, 8, 1).valid ||
       accumulate_issue_cycles(3, 0.5) != 1.5 ||
       accumulate_issue_cycles(0, 0.5) != 0.0 ||
       // Lane->accumulator reduction fill: ceil(log2(mac_width)) stages, width 1 free.
       lane_reduction_fill_cycles(1, 2.0) != 0.0 ||
       lane_reduction_fill_cycles(2, 0.5) != 0.5 ||
       lane_reduction_fill_cycles(5, 1.0) != 3.0 ||
       lane_reduction_fill_cycles(8, 1.0) != 3.0 ||
       calculate_time_based_mac_utilization(84.0, 168.0) != 0.5 ||
       calculate_time_based_mac_utilization(0.0, 168.0) != 0.0 ||
       calculate_time_based_mac_utilization(1.0, 0.0) != 0.0) {
        fail("PE lane capacity or fractional issue cycle");
    }
}

// RE5: cross-check the declared energy schema for self-consistency.
bool energy_schema_is_well_formed() {
    unsigned count = 0;
    const energy_key_schema_t *keys = energy_key_schema(count);
    if(count == 0) return false;
    bool has_unit = false, has_reference = false;
    for(unsigned i = 0; i < count; ++i) {
        const std::string name(keys[i].name ? keys[i].name : "");
        if(name.empty()) return false;
        if(name.find("energy") == std::string::npos) return false;   // else it is never consulted
        for(unsigned j = i + 1; j < count; ++j) {
            if(name == std::string(keys[j].name)) return false;       // duplicate entry
        }
        if(keys[i].kind == ENERGY_KEY_PREFIX_SCALAR && name[name.size()-1] != '_') return false;
        if(name == "energy_unit") has_unit = true;
        if(name == "energy_reference") has_reference = true;
    }
    return has_unit && has_reference;
}

void validate_spatial_interconnect_contract() {
    if(!is_supported_spatial_noc(noc_type_t::BUS) ||
       !is_supported_spatial_noc(noc_type_t::STORE_AND_FORWARD) ||
       !is_supported_spatial_noc(noc_type_t::MESH) ||
       !is_supported_multi_chip_nop(noc_type_t::BUS) ||
       // SP1 pipelined-hop contract: 2x3 mesh max hops = 3 -> fill 2. E2/RE6: the energy quantity
       // depends on the direction -- a multicast is a spanning TREE over the 6 active routers, so
       // 5 edges (NOT 6: the 6th link a receiver count adds is the GLB attach link, priced by the
       // GLB's own transfer_energy), and per-source unicast is the 1.5 average Manhattan distance.
       spatial_noc_cost(noc_type_t::MESH, 2, 3, true).latency_fill_hops != 2.0 ||
       spatial_noc_cost(noc_type_t::MESH, 2, 3, true).link_traversals != 5.0 ||
       spatial_noc_cost(noc_type_t::MESH, 2, 3, false).link_traversals != 1.5 ||
       // RE6 condition 3: the closed form must equal an ENUMERATED edge set, for degenerate and
       // rectangular shapes alike. 1x1 has no internal link at all (its data arrives entirely
       // over the attach link); a 1xN or Nx1 line has N-1; an NxM mesh has N*M-1.
       spatial_multicast_edge_count(1, 1) != 0 ||
       spatial_noc_cost(noc_type_t::MESH, 1, 1, true).link_traversals != 0.0 ||
       spatial_multicast_edge_count(1, 4) != 3 ||
       spatial_noc_cost(noc_type_t::MESH, 1, 4, true).link_traversals != 3.0 ||
       spatial_multicast_edge_count(4, 1) != 3 ||
       spatial_noc_cost(noc_type_t::MESH, 4, 1, true).link_traversals != 3.0 ||
       spatial_multicast_edge_count(2, 3) != 5 ||
       spatial_multicast_edge_count(16, 16) != 255 ||
       spatial_noc_cost(noc_type_t::MESH, 16, 16, true).link_traversals != 255.0 ||
       // ... and the enumerated tree must be a TREE: exactly one edge per non-root router.
       spatial_multicast_edge_count(3, 5) != 3*5 - 1 ||
       spatial_multicast_edge_count(5, 3) != 5*3 - 1 ||
       // The contract statement must exist and differ by direction, so the report cannot print a
       // number without saying which convention produced it.
       std::string(spatial_noc_link_contract(noc_type_t::MESH, true)) ==
           std::string(spatial_noc_link_contract(noc_type_t::MESH, false)) ||
       std::string(spatial_noc_link_contract(noc_type_t::MESH, true)).find("N-1") ==
           std::string::npos ||
       std::string(nop_link_contract(noc_type_t::MESH, true)).find("ingress") ==
           std::string::npos ||
       // RE5: the energy schema itself must be well formed -- a duplicate or empty entry would
       // make lookups depend on declaration order, and the provenance keys must be present or
       // every config's `energy_unit` line would be rejected as an unknown key.
       !energy_schema_is_well_formed() ||
       spatial_noc_cost(noc_type_t::CROSSBAR, 2, 3, true).latency_fill_hops != 0.0 ||
       // R2: mesh NoP is now a supported routed-unicast model; store-and-forward is not.
       !is_supported_multi_chip_nop(noc_type_t::MESH) ||
       is_supported_multi_chip_nop(noc_type_t::STORE_AND_FORWARD) ||
       // Manhattan ingress hops: the same 4 chips route differently on 1x4 vs 2x2.
       nop_route_hops(0, 4) != 1 ||   // ingress-adjacent chip still crosses one link
       nop_route_hops(3, 4) != 3 ||   // 1x4 grid: chip 3 at (3,0) -> 3 hops
       nop_route_hops(3, 2) != 2 ||   // 2x2 grid: chip 3 at (1,1) -> 2 hops
       nop_route_hops(5, 3) != 3 ||   // 2x3 grid: chip 5 at (2,1) -> 3 hops
       // L7 NoP delivery contract. BUS: distinct chunks serialize one copy per chip on the
       // shared medium; a multicast transmission is physically seen by every receiver, so it
       // crosses the medium ONCE however many chips receive it.
       nop_delivery_cost(noc_type_t::BUS, 1, 4, false).bottleneck_link_tiles != 4.0 ||
       nop_delivery_cost(noc_type_t::BUS, 1, 4, true).bottleneck_link_tiles != 1.0 ||
       nop_delivery_cost(noc_type_t::BUS, 1, 4, true).total_link_traversals != 1.0 ||
       nop_delivery_cost(noc_type_t::BUS, 1, 4, false).fill_hops != 0.0 ||
       // A single chip cannot benefit from multicast -- one receiver, one copy either way.
       nop_delivery_cost(noc_type_t::BUS, 1, 1, false).bottleneck_link_tiles != 1.0 ||
       nop_delivery_cost(noc_type_t::MESH, 1, 1, false).bottleneck_link_tiles != 1.0 ||
       nop_delivery_cost(noc_type_t::MESH, 4, 0, false).bottleneck_link_tiles != 0.0 ||
       // MESH, 1x4 line of chips from the ingress at (0,0). Unicast: the ingress carries all
       // 4 copies, and the routes traverse 1+1+2+3 = 7 links in total (chip 0 still crosses
       // the ingress link), with a 3-hop deepest route leaving a 2-hop pipeline fill.
       nop_delivery_cost(noc_type_t::MESH, 4, 4, false).bottleneck_link_tiles != 4.0 ||
       nop_delivery_cost(noc_type_t::MESH, 4, 4, false).total_link_traversals != 7.0 ||
       nop_delivery_cost(noc_type_t::MESH, 4, 4, false).fill_hops != 2.0 ||
       // Multicast over the same line: one copy at the ingress, and the tree is the ingress
       // link plus the 3 row-0 links it fans out along = 4 traversals (vs 7 unicast).
       nop_delivery_cost(noc_type_t::MESH, 4, 4, true).bottleneck_link_tiles != 1.0 ||
       nop_delivery_cost(noc_type_t::MESH, 4, 4, true).total_link_traversals != 4.0 ||
       nop_delivery_cost(noc_type_t::MESH, 4, 4, true).fill_hops != 2.0 ||
       // MESH, 2x2 grid: chips at (0,0),(1,0),(0,1),(1,1); hops 1,1,1,2 -> 5 unicast
       // traversals. The multicast tree is the ingress link + 1 row-0 horizontal + one
       // vertical per used column (1 each) = 4.
       nop_delivery_cost(noc_type_t::MESH, 2, 4, false).total_link_traversals != 5.0 ||
       nop_delivery_cost(noc_type_t::MESH, 2, 4, true).total_link_traversals != 4.0 ||
       nop_delivery_cost(noc_type_t::MESH, 2, 4, true).bottleneck_link_tiles != 1.0 ||
       // B12: derived link width truncates (with a warning) and is 0 without a clock.
       derived_link_bitwidth("test", 2.0, 1.0) != 16 ||
       derived_link_bitwidth("test", 0.9, 1.0) != 7 ||
       derived_link_bitwidth("test", 1.0, 0.0) != 0 ||
       // E2/RE6 spatial NoC energy: MULTICAST distribution is a spanning TREE over the active
       // routers, so N-1 internal links; per-source UNICAST write-back is the average Manhattan
       // distance per transaction. A 16x16 array charges 255 for distribution and 15 for
       // write-back -- the old single formula used 15 for BOTH, understating distribution energy
       // ~17x. (The GLB attach link is the one link NOT counted here; it is priced by the GLB's
       // own transfer_energy, so a receiver count of 256 would bill that wire twice.)
       spatial_noc_cost(noc_type_t::MESH, 16, 16, true).link_traversals != 255.0 ||
       spatial_noc_cost(noc_type_t::MESH, 16, 16, false).link_traversals != 15.0 ||
       // Route depth is a latency quantity and is the same either way.
       spatial_noc_cost(noc_type_t::MESH, 16, 16, true).latency_fill_hops != 29.0 ||
       spatial_noc_cost(noc_type_t::MESH, 16, 16, false).latency_fill_hops != 29.0 ||
       // 1xN row: the tree along the row has N-1 links; a unicast averages (N-1)/2 hops -- for
       // N=8 that is 3.5, and the fill is 7-1 = 6.
       spatial_noc_cost(noc_type_t::MESH, 1, 8, true).link_traversals != 7.0 ||
       spatial_noc_cost(noc_type_t::MESH, 1, 8, false).link_traversals != 3.5 ||
       spatial_noc_cost(noc_type_t::MESH, 1, 8, true).latency_fill_hops != 6.0 ||
       // Nx1 column is the mirror image.
       spatial_noc_cost(noc_type_t::MESH, 8, 1, true).link_traversals != 7.0 ||
       spatial_noc_cost(noc_type_t::MESH, 8, 1, false).link_traversals != 3.5 ||
       // NxM asymmetric: 4x2 has 8 routers -> 7 tree links; unicast average = (sum of Manhattan
       // hops)/8 = (2*4*3/2 + 4*2*1/2)/8 = (12 + 4)/8 = 2.0.
       spatial_noc_cost(noc_type_t::MESH, 4, 2, true).link_traversals != 7.0 ||
       spatial_noc_cost(noc_type_t::MESH, 4, 2, false).link_traversals != 2.0 ||
       // A single PE has NO internal link to cross for a multicast -- its data arrives entirely
       // over the attach link, which this axis does not price. The unicast branch keeps its floor
       // of 1 because a write-back still crosses its own router.
       spatial_noc_cost(noc_type_t::MESH, 1, 1, true).link_traversals != 0.0 ||
       spatial_noc_cost(noc_type_t::MESH, 1, 1, false).link_traversals != 1.0 ||
       spatial_noc_cost(noc_type_t::MESH, 1, 1, true).latency_fill_hops != 0.0 ||
       // A bus transmission reaches every receiver at once and a crossbar is a direct path:
       // one traversal per transaction, no route depth, same in both directions.
       spatial_noc_cost(noc_type_t::BUS, 16, 16, true).link_traversals != 1.0 ||
       spatial_noc_cost(noc_type_t::BUS, 16, 16, false).link_traversals != 1.0 ||
       spatial_noc_cost(noc_type_t::CROSSBAR, 16, 16, true).link_traversals != 1.0 ||
       spatial_noc_cost(noc_type_t::BUS, 16, 16, true).latency_fill_hops != 0.0 ||
       // A degenerate active shape costs nothing rather than underflowing.
       spatial_noc_cost(noc_type_t::MESH, 0, 8, true).link_traversals != 1.0 ||
       !is_supported_spatial_noc(noc_type_t::CROSSBAR) ||
       !has_valid_active_shape(4, 8, 4, 8) ||
       has_valid_active_shape(4, 8, 5, 8) ||
       has_valid_active_shape(4, 8, 4, 9) ||
       // (transactions, source, link, dest). One transaction has no overlap: source+link+dest.
       // CE5: N>=2 uses the ideal makespan src+link+dst+(N-1)*max -- (4,5,1,2)=8+3*5=23
       // (the pre-fix double-max formula gave 24), symmetric for any bottleneck stage.
       pipelined_transfer_cycles(0, 5.0, 1.0, 2.0) != 0.0 ||
       pipelined_transfer_cycles(1, 5.0, 1.0, 2.0) != 8.0 ||
       pipelined_transfer_cycles(4, 5.0, 1.0, 2.0) != 23.0 ||
       pipelined_transfer_cycles(4, 1.0, 5.0, 2.0) != 23.0 ||
       pipelined_transfer_cycles(4, 2.0, 1.0, 5.0) != 23.0 ||
       // L2 PHYSICAL packet decomposition. Every case below is stated as
       // (payload bits, source line, link, destination line) -- a combination some real
       // config could produce -- instead of a hand-picked count triple, because the packet
       // BOUNDARIES are what the makespan depends on and only the widths determine them.
       // Stage costs are (5, 1, 2) throughout.
       //
       // 8b source line -> 256b link -> 8b destination line. One 256b link transaction
       // packs 32 source reads and unpacks into 32 destination writes.
       //
       // 32 B = exactly one packet: the link transaction is a hard barrier in both
       // directions and nothing else is in flight, so 32*5 + 1 + 32*2 = 225.
       packet_makespan(256, 8, 256, 8, 32, 1, 32, 5.0, 1.0, 2.0) != 225.0 ||
       // 64 B = two full packets. Packet 2's source assembly overlaps packet 1's link and
       // destination work, so the second packet costs only its own bottleneck period
       // (32*5 = 160): 225 + 160 = 385. A whole-stream barrier would give 449.
       packet_makespan(512, 8, 256, 8, 64, 2, 64, 5.0, 1.0, 2.0) != 385.0 ||
       // 200 B = SIX full packets plus an 8 B TAIL: 32,32,32,32,32,32,8 -- never the
       // 29,29,29,29,28,28,28 an even split of the aggregate counts would invent. The tail
       // packet costs less than a full period AND is gated by the destination stage, which
       // is still draining packet 6: 6 full packets reach (src 960, link 961, dst 1025),
       // then the tail is src 960->1000, link 1000->1001, dst max(1025,1001)=1025->1041.
       // The pre-fix even-split "first latency + sum of periods" formula gave 1059.
       packet_makespan(1600, 8, 256, 8, 200, 7, 200, 5.0, 1.0, 2.0) != 1041.0 ||
       // A payload smaller than one packet is a lone tail packet: 13*5 + 1 + 13*2 = 92.
       packet_makespan(104, 8, 256, 8, 13, 1, 13, 5.0, 1.0, 2.0) != 92.0 ||
       // 256b source line -> 32b link -> 32b destination line: one wide read FANS OUT to 8
       // link transactions, and the destination streams behind them (no barrier either
       // way), so 5 + 1 + 2 + 7*2 = 22.
       packet_makespan(256, 256, 32, 32, 1, 8, 8, 5.0, 1.0, 2.0) != 22.0 ||
       // 256b source -> 32b link -> 256b destination: fan-out then GATHER. The 8 link
       // transactions stream (5 + 8*1 = 13) but the single gathering write must wait for
       // all of them: 13 + 2 = 15.
       packet_makespan(256, 256, 32, 256, 1, 8, 1, 5.0, 1.0, 2.0) != 15.0 ||
       // Equal widths degenerate to one transaction per stage per packet, i.e. the CE5
       // formula above: 4 packets -> 8 + 3*5 = 23.
       packet_makespan(256, 64, 64, 64, 4, 4, 4, 5.0, 1.0, 2.0) != 23.0 ||
       // The decomposition is reconciled against the aggregate counts the caller charges
       // cycles from, so the per-packet counts always sum back to them.
       packet_counts_consistent(1600, 8, 256, 8, 200, 7, 200) != true ||
       packet_counts_consistent(104, 8, 256, 8, 13, 1, 13) != true ||
       packet_counts_consistent(512, 8, 256, 8, 64, 2, 64) != true ||
       packet_counts_consistent(256, 256, 32, 256, 1, 8, 1) != true ||
       // An absent stage is DROPPED, not charged at zero cost: with no destination the
       // 200 B stream ends when its last link transaction lands (960 -> 1000 -> 1001), and
       // with no source it ends when the last (cheap, 8-write) tail packet drains at 401.
       transfer_packet_groups(1600, 8, 256, 8, 200, 7, 200).without_destination().destination_transactions() != 0 ||
       transfer_packet_groups(1600, 8, 256, 8, 200, 7, 200).without_source().source_transactions() != 0 ||
       pipelined_transfer_cycles(transfer_packet_groups(1600, 8, 256, 8, 200, 7, 200).without_destination(),
                                 5.0, 1.0, 2.0) != 1001.0 ||
       pipelined_transfer_cycles(transfer_packet_groups(1600, 8, 256, 8, 200, 7, 200).without_source(),
                                 5.0, 1.0, 2.0) != 401.0 ||
       // A destination that moves only its SHARE of the payload keeps the shared packet
       // boundaries and spreads its own total over them (25 writes over 7 packets = 4x4+3x3).
       transfer_packet_groups(1600, 8, 256, 8, 200, 7, 200).with_destination_total(25).destination_transactions() != 25 ||
       // Nothing to move is 0, not a fill cost.
       transfer_packet_groups(0, 8, 256, 8, 0, 0, 0).packets != 0 ||
       pipelined_transfer_cycles(transfer_packet_groups(0, 8, 256, 8, 0, 0, 0), 5.0, 1.0, 2.0) != 0.0 ||
       // The steady-state extrapolation must equal the unabridged per-packet simulation --
       // for every bottleneck stage, for payloads that are not a whole number of packets,
       // for large packet counts, and with a stage dropped.
       !packet_extrapolation_exact(1600, 8, 256, 8, 5.0, 1.0, 2.0) ||
       !packet_extrapolation_exact(1600, 8, 256, 8, 1.0, 5.0, 2.0) ||
       !packet_extrapolation_exact(1600, 8, 256, 8, 1.0, 1.0, 9.0) ||
       !packet_extrapolation_exact(1048576, 8, 256, 8, 1.0, 1.0, 1.0) ||
       !packet_extrapolation_exact(1048576, 8, 256, 8, 0.5, 3.0, 0.25) ||
       !packet_extrapolation_exact(1000003, 8, 256, 8, 5.0, 1.0, 2.0) ||
       !packet_extrapolation_exact(4194304, 16, 64, 32, 1.0, 1.0, 1.0) ||
       !packet_extrapolation_exact(999, 8, 64, 256, 2.0, 3.0, 4.0) ||
       !packet_extrapolation_exact(64, 256, 32, 32, 5.0, 1.0, 2.0) ||
       !packet_extrapolation_exact(8, 8, 8, 8, 5.0, 1.0, 2.0) ||
       // P3: temporal_pipeline_run_cycles(). Empty run -> 0.
       temporal_pipeline_run_cycles(4, {}) != 0.0 ||
       // <=1 tile: nothing to pipeline against -- flat max(), NOT the general
       // formula's degenerate sum (40+8+8=56 would silently contradict "these
       // stages overlap").
       temporal_pipeline_run_cycles(1, {40.0, 8.0, 8.0}) != 40.0 ||
       temporal_pipeline_run_cycles(0, {40.0, 8.0, 8.0}) != 40.0 ||
       // >1 tile: per-tile costs (10, 2, 2), fill = 10+2+2 = 14, bottleneck stage
       // (10/tile) repeats 3 more times: 14 + 3*10 = 44 (vs the old flat max of 40 --
       // strictly larger, the one-time fill cost of the non-bottleneck stages).
       temporal_pipeline_run_cycles(4, {40.0, 8.0, 8.0}) != 44.0 ||
       // L5 per-tile timeline. Nothing to run is 0.
       pipeline_timeline_cycles({}, {}) != 0.0 ||
       pipeline_timeline_cycles({{}}, {}) != 0.0 ||
       // A DOUBLE-buffered boundary (depth 2) reproduces the fill+bottleneck closed form
       // for uniform per-tile costs -- the engine is a generalization, not a new model:
       // per-tile (10, 2, 2) over 4 tiles = 14 fill + 3*10 = 44, exactly the value
       // temporal_pipeline_run_cycles(4, {40, 8, 8}) gives above.
       pipeline_timeline_cycles({uniform_tiles(4, 10.0), uniform_tiles(4, 2.0),
                                 uniform_tiles(4, 2.0)}, {2, 2}) != 44.0 ||
       // Two stages, depth 2: 4*10 + 3 = 43 (the consumer's last tile trails the producer).
       pipeline_timeline_cycles({uniform_tiles(4, 10.0), uniform_tiles(4, 3.0)}, {2}) != 43.0 ||
       // A SINGLE-buffered boundary (depth 1) lets no tile-level decoupling happen: the
       // producer may not start tile k+1 until the consumer has taken tile k, so a
       // two-stage chain serializes to 4*(10+3) = 52 -- exactly the sum of the two stage
       // totals, which is what the previous "this boundary does not overlap" branch
       // produced. Both old cases are depths of the same engine.
       pipeline_timeline_cycles({uniform_tiles(4, 10.0), uniform_tiles(4, 3.0)}, {1}) != 52.0 ||
       // A deeper queue than the rates need changes nothing (depth 3 == depth 2 here).
       pipeline_timeline_cycles({uniform_tiles(4, 10.0), uniform_tiles(4, 3.0)}, {3}) != 43.0 ||
       // BACK-PRESSURE, the term an average-per-tile closed form cannot express: a FAST
       // producer in front of a SLOW consumer stalls once the buffer fills instead of
       // running to completion in parallel. Per-tile (3, 10) over 4 tiles at depth 2 costs
       // 3 + 4*10 = 43, and the producer spends 14 cycles blocked.
       pipeline_timeline_cycles({uniform_tiles(4, 3.0), uniform_tiles(4, 10.0)}, {2}) != 43.0 ||
       !stall_is(4, 3.0, 4, 10.0, 2, 14.0, 0.0) ||
       // The same chain with a deep enough buffer never stalls the producer.
       !stall_is(4, 3.0, 4, 10.0, 8, 0.0, 0.0) ||
       // Single buffering stalls the producer even when it is the slow side.
       !stall_is(4, 10.0, 4, 3.0, 1, 9.0, 0.0) ||
       // A RATE-MISMATCHED stage: an off-chip stage that refetches on only 2 of the 4 tiles
       // (cost 20 on those, 0 on the rest) while the on-chip stages serve every tile at 5.
       // The previous closed form had to fall back to a flat max() for any run containing
       // such a stage, because dividing its total by the on-chip tile count is meaningless.
       pipeline_timeline_cycles({sparse_tiles(4, 20.0, 2), uniform_tiles(4, 5.0),
                                 uniform_tiles(4, 5.0)}, {2, 2}) != 60.0 ||
       // A cheaper TAIL tile shortens the run: (10,10,10,2) against a 3/tile consumer is
       // 36, not the 43 a uniform 10/tile producer would give.
       pipeline_timeline_cycles({tail_tiles(4, 10.0, 2.0), uniform_tiles(4, 3.0)}, {2}) != 36.0 ||
       // One tile has nothing to pipeline against: the stages simply chain (10+2+2).
       pipeline_timeline_cycles({uniform_tiles(1, 10.0), uniform_tiles(1, 2.0),
                                 uniform_tiles(1, 2.0)}, {2, 2}) != 14.0 ||
       // CE4/P1-D datatype combination: a shared port serializes its three streams,
       // separate per-datatype partitions serve them concurrently.
       combine_datatype_cycles({10.0, 20.0, 30.0}, true) != 60.0 ||
       combine_datatype_cycles({10.0, 20.0, 30.0}, false) != 30.0 ||
       combine_datatype_cycles({}, true) != 0.0 ||
       // CE4/P1-D entity/datatype ORDER, the asymmetric case that separates the correct
       // reduction from the broken one. Two chips share a GLB; chip 0 spends 100 cycles
       // filling inputs and nothing on weights, chip 1 the mirror image. Each chip's own
       // serialized total is 100, and the chips run in parallel, so the stage axis is 100.
       // Collapsing chips first would give max_chip(input)=100 plus max_chip(weight)=100
       // = 200 cycles of elapsed time that never happened.
       entity_combined_cycles({{100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}}, true) != 100.0 ||
       // Same matrix on separate partitions: each chip's own max is 100 -> 100.
       entity_combined_cycles({{100.0, 0.0, 0.0}, {0.0, 100.0, 0.0}}, false) != 100.0 ||
       // A genuinely busier entity does dominate: chip 1 serializes 40+40+40.
       entity_combined_cycles({{100.0, 0.0, 0.0}, {40.0, 40.0, 40.0}}, true) != 120.0 ||
       entity_combined_cycles({}, true) != 0.0 ||
       // L1 two-sided reduction (a buffer's read side + its fill side). SEPARATE buffer:
       // the read side peaks on input (100) and the fill side on weight (100), but those
       // are two different partitions running in parallel -- each partition's own total is
       // 100, so the axis is 100. Adding the type-combined sides would report 200.
       entity_combined_cycles({{100.0, 0.0, 0.0}}, {{0.0, 100.0, 0.0}}, false) != 100.0 ||
       // Same partition on both sides: the two sides serialize on that partition's port.
       entity_combined_cycles({{100.0, 0.0, 0.0}}, {{40.0, 0.0, 0.0}}, false) != 140.0 ||
       // SHARED buffer: everything serializes on one port, so both forms agree at 200
       // (which is why the separate-buffer error stayed hidden).
       entity_combined_cycles({{100.0, 0.0, 0.0}}, {{0.0, 100.0, 0.0}}, true) != 200.0 ||
       // The Gemmini gemm_64x64x64 live case: separate GLB with base (256,256,1024) and
       // fill (256,256,0) -> per-type totals (512,512,1024) -> axis 1024. The pre-fix
       // form max_type(base)+max_type(fill) = 1024+256 reported 1280.
       entity_combined_cycles({{256.0, 256.0, 1024.0}}, {{256.0, 256.0, 0.0}}, false) != 1024.0 ||
       // Entity max still applies after the per-datatype add: chip 1 wins with 90+100.
       entity_combined_cycles({{10.0, 0.0, 0.0}, {90.0, 0.0, 0.0}},
                              {{100.0, 0.0, 0.0}, {100.0, 0.0, 0.0}}, false) != 190.0 ||
       // A missing row on one side contributes zeros rather than dropping the entity.
       entity_combined_cycles({{10.0, 0.0, 0.0}, {500.0, 0.0, 0.0}},
                              {{100.0, 0.0, 0.0}}, false) != 500.0 ||
       entity_combined_cycles({}, {}, false) != 0.0 ||
       // L8 DRAM open-page row activations. 1024 B over 256 B rows opens 4 rows; with one
       // bank all 4 serialize, with 8 banks the busiest holds 1.
       dram_row_activations(1024, 256, 1).activations != 4 ||
       dram_row_activations(1024, 256, 1).busiest_bank != 4 ||
       dram_row_activations(1024, 256, 8).activations != 4 ||
       dram_row_activations(1024, 256, 8).busiest_bank != 1 ||
       // 4096 B over 256 B rows is 16 activations; 8 banks leaves 2 per bank. Energy tracks
       // the 16 activations either way -- bank parallelism hides latency, not energy.
       dram_row_activations(4096, 256, 8).activations != 16 ||
       dram_row_activations(4096, 256, 8).busiest_bank != 2 ||
       // A partial last row still counts as an activation.
       dram_row_activations(300, 256, 1).activations != 2 ||
       // More banks than activations cannot go below one activation on the busiest bank.
       dram_row_activations(256, 256, 8).busiest_bank != 1 ||
       // Disabled (no row size) or empty streams cost nothing; 0 banks is treated as 1.
       dram_row_activations(1024, 0, 8).activations != 0 ||
       dram_row_activations(0, 256, 8).activations != 0 ||
       dram_row_activations(1024, 256, 0).busiest_bank != 4 ||
       // SY2/L9 systolic fill and drain. A 16x16 active array skews operands (16-1)+(16-1)
       // = 30 hops to the farthest PE, and its partial sums drain 15 hops down the columns
       // to the edge accumulator. Gemmini's RTL-measured ~DIM cycles of per-residency bubble
       // is that same column depth, which is why the drain follows the HEIGHT and not the
       // diameter.
       systolic_pipeline_cost(16, 16).skew_fill_hops != 30.0 ||
       systolic_pipeline_cost(16, 16).drain_hops != 15.0 ||
       // A 1xN row has no accumulation depth to drain, only skew along the row.
       systolic_pipeline_cost(1, 8).skew_fill_hops != 7.0 ||
       systolic_pipeline_cost(1, 8).drain_hops != 0.0 ||
       // An Nx1 column is the mirror image: all depth, no lateral skew beyond it.
       systolic_pipeline_cost(8, 1).skew_fill_hops != 7.0 ||
       systolic_pipeline_cost(8, 1).drain_hops != 7.0 ||
       // A single PE has neither.
       systolic_pipeline_cost(1, 1).skew_fill_hops != 0.0 ||
       systolic_pipeline_cost(1, 1).drain_hops != 0.0 ||
       // A degenerate shape costs nothing rather than underflowing.
       systolic_pipeline_cost(0, 8).drain_hops != 0.0 ||
       systolic_pipeline_cost(8, 0).skew_fill_hops != 0.0 ||
       // L9/PE1-PE2 active-lane reduction. 8 accumulators x 8 lanes, all 64 lanes live:
       // 8 full trees reduce concurrently, so the DEPTH is one tree's ceil(log2(8)) = 3, and
       // the WORK is 8 trees x 7 additions = 56.
       lane_reduction_fill_cycles(calculate_mac_lane_state(8, 8, 64), 1.0) != 3.0 ||
       lane_reduction_energy(calculate_mac_lane_state(8, 8, 64), 1.0) != 56.0 ||
       // Only 8 lanes live -> one full accumulator: depth 3, work 7. The structural form
       // reported the same depth but only ever one tree's work, which is why several live
       // accumulators used to be undercharged (56 vs 7 above).
       lane_reduction_fill_cycles(calculate_mac_lane_state(8, 8, 8), 1.0) != 3.0 ||
       lane_reduction_energy(calculate_mac_lane_state(8, 8, 8), 1.0) != 7.0 ||
       // A PARTIAL last accumulator: 3 lanes live in one unit reduce through a 2-level tree
       // doing 2 additions -- not the full-width 3 levels / 7 additions the structural form
       // charged.
       lane_reduction_fill_cycles(calculate_mac_lane_state(8, 8, 3), 1.0) != 2.0 ||
       lane_reduction_energy(calculate_mac_lane_state(8, 8, 3), 1.0) != 2.0 ||
       lane_reduction_fill_cycles(8, 1.0) != 3.0 ||
       lane_reduction_energy(8, 1.0) != 7.0 ||
       // Several units with a partial last one: 9 lanes = one full 8-lane unit plus a
       // single-lane unit. Depth is the full tree's 3; work is 7 + 0 (a 1-leaf tree adds
       // nothing).
       lane_reduction_fill_cycles(calculate_mac_lane_state(8, 8, 9), 1.0) != 3.0 ||
       lane_reduction_energy(calculate_mac_lane_state(8, 8, 9), 1.0) != 7.0 ||
       // Width 1 has no reduction at all, however many units are live.
       lane_reduction_fill_cycles(calculate_mac_lane_state(8, 1, 8), 1.0) != 0.0 ||
       lane_reduction_energy(calculate_mac_lane_state(8, 1, 8), 1.0) != 0.0 ||
       // An invalid lane state costs nothing rather than charging a fabricated tree.
       lane_reduction_fill_cycles(calculate_mac_lane_state(8, 8, 65), 1.0) != 0.0 ||
       lane_reduction_energy(calculate_mac_lane_state(8, 8, 0), 1.0) != 0.0 ||
       static_energy_for_cycles(2.5, 4.0) != 10.0) {
        fail("spatial interconnect timing contract");
    }
}

void validate_accelerator(config_t &config, const std::string &path) {
    // E8: every energy unit cost in the file must be a finite, non-negative number. Only 3 of
    // the 22 energy keys guarded themselves, so a negative or malformed cost used to reach the
    // model and produce negative/NaN energy that looked like a result.
    const std::string energy_error = validate_energy_settings(config);
    if(!energy_error.empty()) {
        fail(path + ": invalid energy unit cost: " + energy_error);
    }
    unsigned accelerator_count = 0;
    unsigned pe_array_count = 0;
    unsigned global_buffer_count = 0;
    unsigned multi_chip_count = 0;
    unsigned dram_count = 0;
    unsigned num_chips = 0;
    unsigned height = 0;
    unsigned width = 0;

    for(section_config_t &section : config.sections) {
        if(section.name == "accelerator") {
            accelerator_count++;
            if(!section.get_setting("num_chips", &num_chips) || num_chips == 0) {
                fail(path + ": invalid num_chips");
            }
        } else if(section.name == "spatial_arch" ||
                  section.name == "adder_tree" ||
                  section.name == "systolic_array") {
            pe_array_count++;
        } else if(section.name == "shared" || section.name == "separate") {
            global_buffer_count++;
        } else if(section.name == "multi_chip") {
            multi_chip_count++;
            if(!section.get_setting("height", &height) ||
               !section.get_setting("width", &width) ||
               height == 0 || width == 0) {
                fail(path + ": invalid multi-chip dimensions");
            }
        } else if(section.name == "dram") {
            dram_count++;
        } else {
            fail(path + ": unknown accelerator section " + section.name);
        }
    }

    if(accelerator_count != 1 || pe_array_count != 1 ||
       global_buffer_count != 1 || multi_chip_count != 1 || dram_count != 1) {
        fail(path + ": component cardinality");
    }
    if(height > 0 && width > 0 && height * width != num_chips) {
        fail(path + ": num_chips does not match height*width");
    }
}
void validate_mapping(config_t &config, const std::string &path) {
    if(config.sections.empty()) {
        fail(path + ": empty mapping");
    }
    for(section_config_t &section : config.sections) {
        mapping_table_t table(section);
        (void)table.calculate_parameter_size(component_type_t::DRAM);
    }
}

bool contains(const std::string &value, const std::string &needle) {
    return value.find(needle) != std::string::npos;
}

} // namespace

int main(int argc, char **argv) {
    if(argc < 2) {
        std::cerr << "Usage: " << argv[0] << " CONFIG..." << std::endl;
        return 2;
    }

    validate_parser_contract();
    validate_runtime_datatypes();
    validate_hierarchy_datatype_transactions();
    validate_legacy_glb_mapping_contract();
    validate_input_halo_contract();

    validate_pe_lane_contract();
    validate_spatial_interconnect_contract();
    for(int i = 1; i < argc; i++) {
        const std::string path = argv[i];
        config_t config;
        config.parse(path);
        if(contains(path, "/accelerators/")) {
            validate_accelerator(config, path);
        } else if(contains(path, ".map")) {
            validate_mapping(config, path);
        }
    }

    std::cout << "validated " << argc - 1 << " configuration files" << std::endl;
    return 0;
}
