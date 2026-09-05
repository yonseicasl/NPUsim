#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include "stats.h"
#include <limits>
#include "interconnect_timing.h"
#include "energy_units.h"
#include "datatype.h"
#include "pe_lane.h"

namespace {
// SFU list merging: append a token only when no existing SEGMENT equals it exactly.
// Substring matching silently dropped distinct operations whose name is contained in a
// recorded one ("elu" inside "relu"/"gelu", "sigmoid" inside "hsigmoid").
void append_unique_segment(std::string &m_list, const std::string &m_token,
                           const std::string &m_separator) {
    if(m_token.empty()) return;
    if(m_list.empty()) { m_list = m_token; return; }
    size_t start = 0;
    while(true) {
        const size_t end = m_list.find(m_separator, start);
        const std::string segment = (end == std::string::npos)
            ? m_list.substr(start) : m_list.substr(start, end - start);
        if(segment == m_token) return;
        if(end == std::string::npos) break;
        start = end + m_separator.size();
    }
    m_list += m_separator + m_token;
}

// Merge a separator-joined source list into a destination list, segment by segment.
void merge_unique_segments(std::string &m_list, const std::string &m_source,
                           const std::string &m_separator) {
    size_t start = 0;
    while(start <= m_source.size()) {
        const size_t end = m_source.find(m_separator, start);
        const std::string segment = (end == std::string::npos)
            ? m_source.substr(start) : m_source.substr(start, end - start);
        append_unique_segment(m_list, segment, m_separator);
        if(end == std::string::npos) break;
        start = end + m_separator.size();
    }
}

void print_transaction_breakdown(std::ofstream &output, const char *title,
                                 const std::vector<size_t> &payload,
                                 const std::vector<size_t> &metadata,
                                 const std::vector<size_t> &storage) {
    static const char *labels[] = {"Input data", "Weight", "Output data"};
    output << title << " transactions (payload/metadata/serialized)" << std::endl;
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        output << " * " << std::left << std::setw(19) << labels[i] << std::right << ":"
               << payload[i] << "/" << metadata[i] << "/" << storage[i] << std::endl;
    }
    output << std::endl;
}
}

stats_t::stats_t() :
    layer_latency(0.0),
    busy_cycle_pe(0.0),
    busy_cycle_pe_array(0.0),
    busy_cycle_global_buffer(0.0),
    busy_cycle_multi_chip(0.0),
    busy_cycle_dram(0.0),
    local_buffer_type(memory_type_t::UNDEFINED_MEMORY),
    num_computation(0),
    computation_cycle(0.0),
    max_computation_cycle(0.0),
    min_computation_cycle(0.0),
    avg_computation_cycle(0.0),
    computation_energy(0.0),
    mac_busy_cycle(0.0),
    mac_available_cycle(0.0),
    reduction_energy(0.0),
    weight_fold_energy(0.0),
    layer_setup_energy(0.0),
    weight_fold_events(0.0),
    layer_setup_events(0.0),
    stripe_transition_energy(0.0),
    stripe_transition_events(0.0),
    accumulator_reload_bytes(0),
    accumulator_spill_bytes(0),
    accumulator_create_events(0),
    accumulator_retained_events(0),
    output_cast_bytes(0),
    accumulator_energy(0.0),
    output_cast_energy(0.0),
    output_cast_cycle(0.0),
    row_activation_events(0),
    format_payload_events(0),
    format_metadata_events(0),
    reduction_additions(0.0),
    pe_array_accumulator_energy(0.0),
    utilization_mac(0.0),
    total_utilization_local_buffer(0.0),
    utilization_pe_array(0.0),
    fold_fill_cycle_pe_array(0.0),
    layer_setup_cycle_pe_array(0.0),
    stripe_transition_cycle_pe_array(0.0),
    global_buffer_type(memory_type_t::UNDEFINED_MEMORY),
    total_utilization_global_buffer(0.0),
    utilization_multi_chip(0.0) {

    init();

}

stats_t::~stats_t() {
    for(unsigned i = 0; i < component_type_t::NUM_COMPONENT_TYPES; i++) {
        tile_size[i].clear();
    }
}

// Initialize the stats.
void stats_t::init() {
    
    // Reserve the memory for calculating tile size.
    // (Was reserve() + operator[] on an empty outer vector -- undefined behavior that
    // crashed intermittently depending on allocator state; resize() constructs the
    // inner vectors first.)
    tile_size.resize(component_type_t::NUM_COMPONENT_TYPES);
    for(unsigned i = 0; i < component_type_t::NUM_COMPONENT_TYPES; i++) {
        tile_size[i].reserve(data_type_t::NUM_DATA_TYPES);
    }

    /* Initialize PE stats */
    // Initialize the number of request to the local buffer
    num_request_pe.reserve(data_type_t::NUM_DATA_TYPES);
    num_request_pe.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the number of data transfer to MAC unit
    num_data_transfer_pe.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer_pe.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize access cycle to computing units
    access_cycle_mac.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    max_access_cycle_mac.reserve(data_type_t::NUM_DATA_TYPES);
    max_access_cycle_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    min_access_cycle_mac.reserve(data_type_t::NUM_DATA_TYPES);
    min_access_cycle_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    avg_access_cycle_mac.reserve(data_type_t::NUM_DATA_TYPES);
    avg_access_cycle_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize access energy to computing units
    access_energy_mac.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize local buffer access cycle
    access_cycle_lb.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    max_access_cycle_lb.reserve(data_type_t::NUM_DATA_TYPES);
    max_access_cycle_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    min_access_cycle_lb.reserve(data_type_t::NUM_DATA_TYPES);
    min_access_cycle_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    avg_access_cycle_lb.reserve(data_type_t::NUM_DATA_TYPES);
    avg_access_cycle_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    
    // Initialize local buffer access energy
    access_energy_lb.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize local buffer utilization 
    utilization_local_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    utilization_local_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);
        
    // Initialize overlapped cycle between the computing unit and local buffer
    cycle_mac_lb.reserve(data_type_t::NUM_DATA_TYPES);
    cycle_mac_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize transfer cycle between computing unit and local buffers
    transfer_cycle_pe.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle_pe.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize transfer energy between computing unit and local buffers
    transfer_energy_pe.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy_pe.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    payload_link_transactions_pe.assign(data_type_t::NUM_DATA_TYPES, 0);
    metadata_link_transactions_pe.assign(data_type_t::NUM_DATA_TYPES, 0);
    storage_link_transactions_pe.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize static energy at PE
    static_energy_pe.reserve(data_type_t::NUM_DATA_TYPES);
    static_energy_pe.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    static_energy_pe_array.reserve(data_type_t::NUM_DATA_TYPES);
    static_energy_pe_array.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    utilization_pe_array_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    utilization_pe_array_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    utilization_multi_chip_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    utilization_multi_chip_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    format_cycle_pe.reserve(data_type_t::NUM_DATA_TYPES);
    format_cycle_pe.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    format_energy_pe.reserve(data_type_t::NUM_DATA_TYPES);
    format_energy_pe.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    
    /* Initialize PE array stats */

    // Initialize the number of request to PE array
    num_request_pe_array.reserve(data_type_t::NUM_DATA_TYPES);
    num_request_pe_array.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the number of data transfer of PE array (to PE)
    num_data_transfer_pe_array.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer_pe_array.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize access cycle to PE array (if temporal buffer exist)
    access_cycle_pe_array.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle_pe_array.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize access energy to PE array (if temporal buffer exist)
    access_energy_pe_array.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy_pe_array.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize transfer cycle between PE array and PE (Interconnection)
    transfer_cycle_pe_array.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle_pe_array.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    cycle_temporal_pe_array.reserve(data_type_t::NUM_DATA_TYPES);
    cycle_temporal_pe_array.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize transfer energy between PE array and PE (Interconnection)
    transfer_energy_pe_array.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy_pe_array.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    payload_link_transactions_pe_array.assign(data_type_t::NUM_DATA_TYPES, 0);
    metadata_link_transactions_pe_array.assign(data_type_t::NUM_DATA_TYPES, 0);
    storage_link_transactions_pe_array.assign(data_type_t::NUM_DATA_TYPES, 0);

    /* Initialize Global buffer stats */

    // Initialize the number of request to the global buffer
    num_request_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    num_request_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the number of data transfer to the PE array
    num_data_transfer_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0);
    psum_writeback_events_global_buffer = 0;
    payload_link_transactions_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0);
    metadata_link_transactions_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0);
    storage_link_transactions_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0);
    
    // Initialize access cycle of the global buffer
    access_cycle_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    fill_access_cycle_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    fill_access_energy_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    for(unsigned s = 0; s < 5; ++s) {
        stage_axis_access[s] = stage_axis_link[s] = stage_axis_overlap[s] = 0.0;
    }
    stage_axis_compute = 0.0;
    stage_axis_format = 0.0;
    entity_combined_format = 0.0;
    entity_combined_access_global_buffer = 0.0;
    entity_combined_link_global_buffer = 0.0;
    entity_combined_overlap_global_buffer = 0.0;
    entity_combined_access_pe_array = 0.0;
    entity_combined_link_pe_array = 0.0;
    entity_combined_overlap_pe_array = 0.0;
    entity_combined_access_lb = 0.0;
    entity_combined_link_pe = 0.0;
    entity_combined_overlap_mac_lb = 0.0;
    stage_fill_access_global_buffer = 0.0;
    global_buffer_bypass.assign(data_type_t::NUM_DATA_TYPES, false);
    network_rollup = false;
    network_rollup_mapped_seeded = false;
    output_tile_array_resident = false;
    compute_energy_basis = "computation_energy";
    authoritative_frequency_mhz = 0.0;
    single_clock_domain = false;
    clock_domain_note.clear();
    compute_energy_precision_calibrated = false;
    operand_precision.clear();
    network_timing_layers = 0;
    excluded_timing_layers = 0;
    chip_access_cycle_global_buffer.clear();
    chip_fill_access_cycle_global_buffer.clear();
    pe_double_buffer = true;
    timeline_boundary_overlap[0] = timeline_boundary_overlap[1] = false;
    timeline_boundary_overlap[2] = timeline_boundary_overlap[3] = true;
    timeline_physical_macs = 0.0;
    offchip_repetition_tiles = 1;
    input_halo_overlap = false;
    input_halo_capacity_sufficient = false;
    input_halo_reuse_applied = false;
    input_halo_unique_elements = 0;
    input_halo_replicated_elements = 0;
    input_halo_working_set_bytes = 0;
    input_halo_pre_dram_transactions = 0;
    dram_input_link_bits = 0;
    dram_input_source_line_bits = 0;
    dram_input_read_cycle = 0.0;
    dram_input_read_energy = 0.0;
    dram_input_access_hidden = false;
    for(unsigned s = 0; s < 5; ++s) timeline_stall[s] = 0.0;
    for(unsigned b = 0; b < 4; ++b) timeline_boundary_depth[b] = 1;
    timeline_offchip_outstanding = 0;
    psum_array_resident = false;
    psum_store_access_cycle_pe_array = 0.0;
    psum_store_access_energy_pe_array = 0.0;
    temporal_repetition_tiles = 1;

    /* SFU */
    sfu_present = false;
    sfu_active = false;
    sfu_operation.clear();
    sfu_physical_units = 0;
    sfu_lanes = 0;
    sfu_valid_elements = 0;
    sfu_operations = 0;
    sfu_invocations = 0;
    sfu_chunks = 0;
    sfu_commit_events = 0;
    sfu_commit_note.clear();
    sfu_profile_reference.clear();
    sfu_busy_cycle = 0.0;
    sfu_serial_cycle = 0.0;
    sfu_hidden_cycle = 0.0;
    sfu_stall_cycle = 0.0;
    sfu_queue_depth = 0;
    sfu_tail_lane_utilization = 0.0;
    sfu_pending = false;
    sfu_pending_invocation = sfu_invocation_t();
    sfu_pending_static_energy_per_cycle = 0.0;
    output_repetition_tiles = 1;
    sfu_stream_active = false;
    sfu_stream_residency.clear();
    sfu_stream_ingress_bytes = 0;
    sfu_stream_egress_bytes = 0;
    sfu_stream_ingress_cycle = 0.0;
    sfu_stream_egress_cycle = 0.0;
    sfu_ingress_elements = 0;
    sfu_egress_elements = 0;
    sfu_ingress_bytes = 0;
    sfu_egress_bytes = 0;
    sfu_ingress_transactions = 0;
    sfu_egress_transactions = 0;
    sfu_op_energy = 0.0;
    sfu_read_energy = 0.0;
    sfu_write_energy = 0.0;
    sfu_setup_energy = 0.0;
    sfu_static_energy = 0.0;
    sfu_timing_calibrated = true;
    sfu_only_layer = false;
    sfu_unpriced_events.clear();
    sfu_contract_note.clear();
    graph_residency_applied = false;
    graph_residency_note.clear();
    activation_unmodeled = false;
    activation_unmodeled_note.clear();

    /* Weight decompression */
    decomp_present = false;
    decomp_active = false;
    decomp_bypassed = false;
    decomp_dense_weight_bytes = 0;
    decomp_compressed_weight_bytes = 0;
    decomp_effective_ratio = 1.0;
    decomp_tiles = 0;
    decomp_decoder_cycles = 0.0;
    decomp_exposed_cycle = 0.0;
    decomp_hidden_cycle = 0.0;
    decomp_stall_cycle = 0.0;
    decomp_dram_weight_saved_cycle = 0.0;
    decomp_queue_depth = 0;
    decomp_overlap = false;
    decomp_decoder_energy = 0.0;
    decomp_static_energy = 0.0;
    decomp_timing_calibrated = true;
    decomp_unpriced_events.clear();
    decomp_profile_reference.clear();
    decomp_breakdown_dram = 0.0;
    decomp_breakdown_decode = 0.0;
    decomp_breakdown_compute = 0.0;
    decomp_breakdown_stall = 0.0;
    decomp_pending = false;
    decomp_pending_static_energy_per_cycle = 0.0;
    decomp_decoder_ratio = 0.0;
    decomp_startup_cycles = 0.0;
    decomp_tile_supply_fraction.clear();
    decomp_bypassed_tiles = 0;
    decomp_tile_ratio_cv = 0.0;
    decomp_output_buffer_tiles = 1;
    decomp_sink_cycles = 0.0;

    kv_present = false;
    kv_pending = false;
    kv_read_bytes = 0;
    kv_dense_read_bytes = 0;
    kv_bytes_per_token = 0;
    kv_context_length = 0;
    kv_compression_ratio = 1.0;
    kv_bypassed = false;
    kv_read_cycles = 0.0;
    kv_dram_read_cycles = 0.0;
    kv_decoder_cycles = 0.0;
    kv_decoder_calibrated = true;
    kv_read_energy = 0.0;
    kv_priced = true;
    kv_unpriced.clear();
    kv_profile_reference.clear();
    kv_schedule.clear();
    kv_sched_pending = false;
    kv_tiles = 0;
    kv_buffer_tiles = 1;
    kv_exposed_cycle = 0.0;
    kv_hidden_cycle = 0.0;
    kv_stall_cycle = 0.0;

    attn_present = false;
    attn_pending = false;
    attn_algorithm.clear();
    attn_qk_macs = 0;
    attn_softmax_elements = 0;
    attn_qk_cycles = 0.0;
    attn_av_cycles = 0.0;
    attn_softmax_cycles = 0.0;
    attn_supply_k_cycle = 0.0;
    attn_supply_v_cycle = 0.0;
    attn_write_cycle = 0.0;
    attn_step_latency = 0.0;
    attn_exposed_cycle = 0.0;
    attn_hidden_cycle = 0.0;
    attn_stall_cycle = 0.0;
    attn_compute_energy = 0.0;
    attn_compute_calibrated = true;
    attn_kv_occupancy_bytes = 0;
    attn_kv_capacity_bytes = 0;

    // Initialize access energy of the global buffer
    access_energy_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize global buffer utilization
    utilization_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    utilization_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize overlapped cycle between the global buffer and the PE array
    cycle_pe_array_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    cycle_pe_array_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize transfer cycle between the global buffer and the PE array
    transfer_cycle_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize transfer energy between the global buffer and the PE array
    transfer_energy_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize global buffer static energy
    static_energy_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    static_energy_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    /* Initialize Multi chip stats */

    // Initialize the number of request to chip-level processor
    num_request_multi_chip.reserve(data_type_t::NUM_DATA_TYPES);
    num_request_multi_chip.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the number of data transfer to the global buffer
    num_data_transfer_multi_chip.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer_multi_chip.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    payload_link_transactions_multi_chip.assign(data_type_t::NUM_DATA_TYPES, 0);
    metadata_link_transactions_multi_chip.assign(data_type_t::NUM_DATA_TYPES, 0);
    storage_link_transactions_multi_chip.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize access cycle of the chip-level processor (if temporal buffer exist)
    access_cycle_multi_chip.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle_multi_chip.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize access energy of the chip-level processor (if temporal buffer exist)
    access_energy_multi_chip.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy_multi_chip.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize transfer cycle between the chip-level processor and the global buffer (Network on Package)
    transfer_cycle_multi_chip.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle_multi_chip.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize transfer energy between the chip-level processor and the global buffer (Network on Package)
    transfer_energy_multi_chip.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy_multi_chip.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize static (leakage) energy of the Multi-chip temporal buffer
    static_energy_multi_chip.reserve(data_type_t::NUM_DATA_TYPES);
    static_energy_multi_chip.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    /* Initialize off-chip memory stats */
    // Initialize the number of request to the off-chip memory
    num_request_dram.reserve(data_type_t::NUM_DATA_TYPES);
    num_request_dram.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the number of data transfer to chip-level processor
    num_data_transfer_dram.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer_dram.assign(data_type_t::NUM_DATA_TYPES, 0);
    payload_link_transactions_dram.assign(data_type_t::NUM_DATA_TYPES, 0);
    metadata_link_transactions_dram.assign(data_type_t::NUM_DATA_TYPES, 0);
    storage_link_transactions_dram.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize access cycle to the off-chip memory
    access_cycle_dram.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle_dram.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize access energy to the off-chip memory
    access_energy_dram.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy_dram.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize overlapped cycle between the off-chip memory and on-chip processor
    cycle_chip_dram.reserve(data_type_t::NUM_DATA_TYPES);
    cycle_chip_dram.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize transfer cycle between the off-chip memory and on-chip processor
    transfer_cycle_dram.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle_dram.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize transfer energy between the off-chip memory and on-chip processor
    transfer_energy_dram.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy_dram.assign(data_type_t::NUM_DATA_TYPES, 0.0);

}

// Print out the stats of accelerator and neural network.
void stats_t::print_stats() {
    std::cout << "* Tile size" << std::endl;
    std::cout << "========= MAC =========="<< std::endl;
    std::cout << "Input data  : " << std::setw(10) 
                                  << tile_size[component_type_t::MAC][data_type_t::INPUT]  << std::endl;
    std::cout << "Weight      : " << std::setw(10) 
                                  << tile_size[component_type_t::MAC][data_type_t::WEIGHT] << std::endl;
    std::cout << "Output data : " << std::setw(10) 
                                  << tile_size[component_type_t::MAC][data_type_t::OUTPUT] << std::endl;
    std::cout << "========================"<< std::endl;
    std::cout << std::endl;

    std::cout << "===== Local buffer =====" << std::endl;
    std::cout << "Input data  : " << std::setw(10) 
                                  << tile_size[component_type_t::PE][data_type_t::INPUT]  << std::endl;
    std::cout << "Weight      : " << std::setw(10) 
                                  << tile_size[component_type_t::PE][data_type_t::WEIGHT] << std::endl;
    std::cout << "Output data : " << std::setw(10) 
                                  << tile_size[component_type_t::PE][data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;

    std::cout << "======= PE array =======" << std::endl;
    std::cout << "Input data  : " << std::setw(10) 
                                  << tile_size[component_type_t::PE_Y][data_type_t::INPUT]  << std::endl;
    std::cout << "Weight      : " << std::setw(10) 
                                  << tile_size[component_type_t::PE_Y][data_type_t::WEIGHT] << std::endl;
    std::cout << "Output data : " << std::setw(10) 
                                  << tile_size[component_type_t::PE_Y][data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;

    std::cout << "===== Global buffer ====" << std::endl;
    std::cout << "Input data  : " << std::setw(10) 
                                  << tile_size[component_type_t::GLOBAL_BUFFER][data_type_t::INPUT]  << std::endl;
    std::cout << "Weight      : " << std::setw(10) 
                                  << tile_size[component_type_t::GLOBAL_BUFFER][data_type_t::WEIGHT] << std::endl;
    std::cout << "Output data : " << std::setw(10) 
                                  << tile_size[component_type_t::GLOBAL_BUFFER][data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;

    std::cout << "======= Processor ======" << std::endl;
    std::cout << "Input data  : " << std::setw(10) 
                                  << tile_size[component_type_t::CHIPS_Y][data_type_t::INPUT]  << std::endl;
    std::cout << "Weight      : " << std::setw(10) 
                                  << tile_size[component_type_t::CHIPS_Y][data_type_t::WEIGHT] << std::endl;
    std::cout << "Output data : " << std::setw(10) 
                                  << tile_size[component_type_t::CHIPS_Y][data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;

    std::cout << "========= DRAM =========" << std::endl;
    std::cout << "Input data  : " << std::setw(10) 
                                  << tile_size[component_type_t::DRAM][data_type_t::INPUT]  << std::endl;
    std::cout << "Weight      : " << std::setw(10) 
                                  << tile_size[component_type_t::DRAM][data_type_t::WEIGHT] << std::endl;
    std::cout << "Output data : " << std::setw(10) 
                                  << tile_size[component_type_t::DRAM][data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;

}

// Print out the stats of accelerator and neural network.
void stats_t::print_stats(std::ofstream &m_output_file) {
    m_output_file << "* Tile size" << std::endl;
    m_output_file << "========= MAC =========="<< std::endl;
    m_output_file << "Input data  : " << std::setw(10) 
                                      << tile_size[component_type_t::MAC][data_type_t::INPUT]  << std::endl;
    m_output_file << "Weight      : " << std::setw(10) 
                                      << tile_size[component_type_t::MAC][data_type_t::WEIGHT] << std::endl;
    m_output_file << "Output data : " << std::setw(10) 
                                      << tile_size[component_type_t::MAC][data_type_t::OUTPUT] << std::endl;
    m_output_file << "========================"<< std::endl;
    m_output_file << std::endl;

    m_output_file << "===== Local buffer =====" << std::endl;
    m_output_file << "Input data  : " << std::setw(10) 
                                      << tile_size[component_type_t::PE][data_type_t::INPUT]  << std::endl;
    m_output_file << "Weight      : " << std::setw(10) 
                                      << tile_size[component_type_t::PE][data_type_t::WEIGHT] << std::endl;
    m_output_file << "Output data : " << std::setw(10) 
                                      << tile_size[component_type_t::PE][data_type_t::OUTPUT] << std::endl;
    m_output_file << "========================" << std::endl;
    m_output_file << std::endl;

    m_output_file << "======= PE array =======" << std::endl;
    m_output_file << "Input data  : " << std::setw(10) 
                                      << tile_size[component_type_t::PE_Y][data_type_t::INPUT]  << std::endl;
    m_output_file << "Weight      : " << std::setw(10) 
                                      << tile_size[component_type_t::PE_Y][data_type_t::WEIGHT] << std::endl;
    m_output_file << "Output data : " << std::setw(10) 
                                      << tile_size[component_type_t::PE_Y][data_type_t::OUTPUT] << std::endl;
    m_output_file << "========================" << std::endl;
    m_output_file << std::endl;

    m_output_file << "===== Global buffer ====" << std::endl;
    m_output_file << "Input data  : " << std::setw(10) 
                                      << tile_size[component_type_t::GLOBAL_BUFFER][data_type_t::INPUT]  << std::endl;
    m_output_file << "Weight      : " << std::setw(10) 
                                      << tile_size[component_type_t::GLOBAL_BUFFER][data_type_t::WEIGHT] << std::endl;
    m_output_file << "Output data : " << std::setw(10) 
                                      << tile_size[component_type_t::GLOBAL_BUFFER][data_type_t::OUTPUT] << std::endl;
    m_output_file << "========================" << std::endl;
    m_output_file << std::endl;

    m_output_file << "======= Processor ======" << std::endl;
    m_output_file << "Input data  : " << std::setw(10) 
                                      << tile_size[component_type_t::CHIPS_Y][data_type_t::INPUT]  << std::endl;
    m_output_file << "Weight      : " << std::setw(10) 
                                      << tile_size[component_type_t::CHIPS_Y][data_type_t::WEIGHT] << std::endl;
    m_output_file << "Output data : " << std::setw(10) 
                                      << tile_size[component_type_t::CHIPS_Y][data_type_t::OUTPUT] << std::endl;
    m_output_file << "========================" << std::endl;
    m_output_file << std::endl;

    m_output_file << "========= DRAM =========" << std::endl;
    m_output_file << "Input data  : " << std::setw(10) 
                                      << tile_size[component_type_t::DRAM][data_type_t::INPUT]  << std::endl;
    m_output_file << "Weight      : " << std::setw(10) 
                                      << tile_size[component_type_t::DRAM][data_type_t::WEIGHT] << std::endl;
    m_output_file << "Output data : " << std::setw(10) 
                                      << tile_size[component_type_t::DRAM][data_type_t::OUTPUT] << std::endl;
    m_output_file << "========================" << std::endl;
    m_output_file << std::endl;
}

void stats_t::update_tile_size(scheduler_t *m_scheduler) {
    tile_size = m_scheduler->tile_size;
}

void stats_t::update_stats(std::vector<pe_array_t*> m_pe_array, std::vector<global_buffer_t*> m_global_buffer, multi_chip_t *m_multi_chip, dram_t *m_dram) {

    // LB7: AND-reduced across every active PE below; stays true only if every PE's
    // local buffer is double-buffered.
    pe_double_buffer = true;

    // CE4/P1-D/L1: per-chip, per-datatype GLB accumulators, filled in the chip loop below
    // and combined only after repetition scaling (see merge_global_buffer_fill()).
    chip_access_cycle_global_buffer.clear();
    chip_fill_access_cycle_global_buffer.clear();

    // static_energy in the PE/GLB config is leakage energy in pJ/cycle. Leakage is an
    // always-on power: every physical component of every physical chip (active or not)
    // leaks for the whole modeled duration of the layer, not once per data-transfer
    // callback. That duration is the layer critical path, which must include memory
    // (GLB/NoP/DRAM) busy time -- charging leakage over the PE-array compute window
    // alone undercounts it heavily for memory-bound layers.
    double pe_elapsed_cycles = 0.0;
    size_t physical_scalar_macs = 0;
    for(unsigned i = 0; i < m_multi_chip->get_number_of_active_chips(); ++i) {
        for(unsigned j = 0; j < m_pe_array[i]->get_number_of_pes(); ++j) {
            pe_elapsed_cycles = std::max(pe_elapsed_cycles,
                                         m_pe_array[i]->pes[j]->modeled_elapsed_cycles());
            const size_t scalar_macs = m_pe_array[i]->pes[j]->scalar_mac_capacity();
            if(scalar_macs > std::numeric_limits<size_t>::max() - physical_scalar_macs) {
                std::cerr << "Error: scalar-MAC capacity overflow" << std::endl;
                exit(1);
            }
            physical_scalar_macs += scalar_macs;
        }
    }

    // A1 (global cycle / compute-memory overlap): place each hierarchy level's total
    // busy time on a shared analytical timeline. Adjacent stages overlap (pipeline over
    // tiles) when the boundary's downstream side buffers tiles -- the same
    // destination-side double-buffer convention the access-cycle accounting uses; a
    // single-buffered boundary serializes its two sides. The GLB->PE-array->PE path is
    // one continuous stream, so those stages always overlap. The layer critical-path
    // latency is the segment-combined result: max within an overlapped run of stages,
    // summed across serialized boundaries.
    busy_cycle_dram = m_dram->modeled_elapsed_cycles();
    // Phase-5: establish the clock contract. Every modeled component must share one frequency
    // for the shared cycle axis to be convertible to time; otherwise power is unsupported and
    // says so rather than being computed against an arbitrary domain.
    {
        const double frequencies[5] = {
            m_pe_array[0]->pes[0]->clock_mhz(),
            m_pe_array[0]->clock_mhz(),
            m_global_buffer[0]->clock_mhz(),
            m_multi_chip->clock_mhz(),
            m_dram->clock_mhz()};
        const char *names[5] = {"PE", "PE array", "global buffer", "multi-chip", "DRAM"};
        single_clock_domain = true;
        authoritative_frequency_mhz = frequencies[0];
        std::string mismatch;
        for(unsigned f = 1; f < 5; ++f) {
            if(frequencies[f] != frequencies[0]) {
                single_clock_domain = false;
                mismatch += std::string(mismatch.empty() ? "" : ", ") + names[f] + " " +
                            std::to_string(frequencies[f]) + " MHz";
            }
        }
        if(frequencies[0] <= 0.0) {
            single_clock_domain = false;
            clock_domain_note = "no positive component frequency is declared";
        } else if(single_clock_domain) {
            clock_domain_note = "all modeled components share one clock";
        } else {
            clock_domain_note = "mixed clock domains (PE " + std::to_string(frequencies[0]) +
                                " MHz vs " + mismatch +
                                "); the timeline is one shared cycle axis, so cycles from"
                                " different domains are not comparable";
            authoritative_frequency_mhz = 0.0;
        }
    }
    // L8: record the analytical DRAM contract so the report can state its scope.
    row_activation_events = m_dram->row_activation_events;
    dram_input_link_bits = m_dram->get_bitwidth();
    dram_input_source_line_bits = m_dram->line_size[data_type_t::INPUT];
    dram_input_read_cycle = m_dram->u_read_cycle[data_type_t::INPUT];
    dram_input_read_energy = m_dram->u_read_energy[data_type_t::INPUT];
    dram_input_access_hidden = m_multi_chip->double_buffer;
    dram_timing_model = m_dram->describe_timing_model();
    // RE6: record which links each fabric's noc_energy prices, so the reported energy can be
    // checked against an actual edge count instead of guessed at.
    if(!m_pe_array.empty()) {
        array_noc_link_contract = m_pe_array[0]->describe_noc_link_contract();
        output_tile_array_resident = m_pe_array[0]->output_tile_is_array_resident();
    }
    nop_link_contract = m_multi_chip->describe_nop_link_contract();
    dram_timing_limits = m_dram->describe_timing_limits();
    busy_cycle_multi_chip = m_multi_chip->modeled_elapsed_cycles();
    // RE1: the final output cast, charged at the off-chip output store.
    output_cast_bytes = m_multi_chip->output_cast_bytes;
    output_cast_energy = m_multi_chip->output_cast_energy;
    output_cast_cycle = m_multi_chip->output_cast_cycle;
    busy_cycle_global_buffer = 0.0;
    busy_cycle_pe_array = 0.0;
    bool global_buffer_double_buffer = false;
    for(unsigned i = 0; i < m_multi_chip->get_number_of_chips(); ++i) {
        busy_cycle_global_buffer = std::max(busy_cycle_global_buffer,
                                            m_global_buffer[i]->modeled_elapsed_cycles());
        busy_cycle_pe_array = std::max(busy_cycle_pe_array,
                                       m_pe_array[i]->modeled_elapsed_cycles());
        global_buffer_double_buffer = global_buffer_double_buffer || m_global_buffer[i]->double_buffer;
    }
    busy_cycle_pe = pe_elapsed_cycles;

    const double stage_busy[5] = {busy_cycle_dram, busy_cycle_multi_chip,
                                  busy_cycle_global_buffer, busy_cycle_pe_array, busy_cycle_pe};
    const bool boundary_overlaps[4] = {m_multi_chip->double_buffer,   // DRAM | multi-chip
                                       global_buffer_double_buffer,   // multi-chip | GLB
                                       true,                          // GLB | PE array (same stream)
                                       true};                         // PE array | PE (same stream)
    // CE1: keep the boundary contract and MAC capacity for the post-scaling
    // timeline recomputation (finalize_layer_timeline).
    for(unsigned b = 0; b < 4; ++b) timeline_boundary_overlap[b] = boundary_overlaps[b];
    // MC2: carry the outstanding-request depth alongside the boolean contract, because
    // finalize_layer_timeline() recomputes the timeline after scaling and has no components.
    timeline_offchip_outstanding = m_multi_chip->max_outstanding_requests;
    // E20-3b: whether the array's output tile survives the reduction. Decides below whether the
    // GLB<->array OUTPUT leg repeats with the collapsed reduction loop or only with the distinct
    // output tiles.
    psum_array_resident = m_pe_array.empty() ? false : m_pe_array[0]->psum_retention_valid;
    timeline_physical_macs = static_cast<double>(physical_scalar_macs);
    layer_latency = 0.0;
    double segment = stage_busy[0];
    for(unsigned b = 0; b < 4; ++b) {
        if(boundary_overlaps[b]) {
            segment = std::max(segment, stage_busy[b + 1]);
        } else {
            layer_latency += segment;
            segment = stage_busy[b + 1];
        }
    }
    layer_latency += segment;

    // Leakage accrues over the layer critical-path latency (always-on components).
    for(unsigned i = 0; i < m_multi_chip->get_number_of_chips(); ++i) {
        for(unsigned j = 0; j < m_pe_array[i]->get_number_of_pes(); ++j) {
            m_pe_array[i]->pes[j]->update_static_energy(layer_latency);
        }
        m_pe_array[i]->update_static_energy(layer_latency);
    }
    for(unsigned i = 0; i < m_multi_chip->get_number_of_chips(); ++i) {
        m_global_buffer[i]->update_static_energy(layer_latency);
    }
    m_multi_chip->update_static_energy(layer_latency);

    unsigned num_active_pe = 0;
    bool first_active_pe = true;   // CE6: MIN counters initialize once across ALL chips
    for(unsigned i = 0; i < m_multi_chip->get_number_of_active_chips(); i++) {
        for(unsigned j = 0; j < m_pe_array[i]->get_number_of_active_pes(); j++) {
            /* Update PE stats */
            // Update the number of computation
            num_computation +=  m_pe_array[i]->pes[j]->num_computation;

            // Update computation cycle 
            computation_cycle = std::max(computation_cycle, m_pe_array[i]->pes[j]->computation_cycle);

            max_computation_cycle = std::max(max_computation_cycle, m_pe_array[i]->pes[j]->computation_cycle);
            min_computation_cycle = first_active_pe ? m_pe_array[i]->pes[j]->computation_cycle : std::min(min_computation_cycle, m_pe_array[i]->pes[j]->computation_cycle);
            avg_computation_cycle += m_pe_array[i]->pes[j]->computation_cycle;

            // Update computation energy
            computation_energy += m_pe_array[i]->pes[j]->computation_energy;
            // E4: the reduction axis and the two format event counts add across PEs.
            reduction_energy += m_pe_array[i]->pes[j]->reduction_energy;
            accumulator_reload_bytes += m_pe_array[i]->pes[j]->accumulator_reload_bytes;
            accumulator_spill_bytes += m_pe_array[i]->pes[j]->accumulator_spill_bytes;
            accumulator_create_events += m_pe_array[i]->pes[j]->accumulator_create_events;
            accumulator_retained_events += m_pe_array[i]->pes[j]->accumulator_retained_events;
            accumulator_energy += m_pe_array[i]->pes[j]->accumulator_energy;

            // Account for scalar-MAC cycles that actually execute operations.
            mac_busy_cycle += static_cast<double>(m_pe_array[i]->pes[j]->num_computation) *
                              m_pe_array[i]->pes[j]->u_computation_cycle;

            // CE4: this PE's own datatype combination, reduced across entities (below)
            // in the correct max_entity(sum/max_type(...)) order.
            double pe_access_lb_sum_types = 0.0, pe_access_lb_max_type = 0.0;
            double pe_link_pe_types = 0.0;
            double pe_overlap_mac_lb_types = 0.0;
            double pe_format_types = 0.0;
            for(unsigned k = 0; k < data_type_t::NUM_DATA_TYPES; k++) {
                // Update the number of request to local buffer in PE
                num_request_pe[k] += m_pe_array[i]->pes[j]->num_request_to_lb[k];

                // Update the number of data transfer to the computing unit in PE
                num_data_transfer_pe[k] += m_pe_array[i]->pes[j]->num_data_transfer_to_mac[k];
                
                // Update access cycle of the computing unit
                access_cycle_mac[k] = std::max(access_cycle_mac[k], m_pe_array[i]->pes[j]->access_cycle_mac[k]);
                
                max_access_cycle_mac[k] = std::max(max_access_cycle_mac[k], m_pe_array[i]->pes[j]->access_cycle_mac[k]);
                min_access_cycle_mac[k] = first_active_pe ? m_pe_array[i]->pes[j]->access_cycle_mac[k] : std::min(min_access_cycle_mac[k], m_pe_array[i]->pes[j]->access_cycle_mac[k]);
                avg_access_cycle_mac[k] += m_pe_array[i]->pes[j]->access_cycle_mac[k];

                // Update access energy of the computing units
                access_energy_mac[k] += m_pe_array[i]->pes[j]->access_energy_mac[k];
                
                // Update access cycle of the local buffer
                access_cycle_lb[k] = std::max(access_cycle_lb[k], m_pe_array[i]->pes[j]->access_cycle_lb[k]);
                pe_access_lb_sum_types += m_pe_array[i]->pes[j]->access_cycle_lb[k];
                pe_access_lb_max_type = std::max(pe_access_lb_max_type, m_pe_array[i]->pes[j]->access_cycle_lb[k]);

                max_access_cycle_lb[k] = std::max(max_access_cycle_lb[k], m_pe_array[i]->pes[j]->access_cycle_lb[k]);
                min_access_cycle_lb[k] = first_active_pe ? m_pe_array[i]->pes[j]->access_cycle_lb[k] : std::min(min_access_cycle_lb[k], m_pe_array[i]->pes[j]->access_cycle_lb[k]);
                avg_access_cycle_lb[k] += m_pe_array[i]->pes[j]->access_cycle_lb[k];

                // Update access energy of the local buffer
                access_energy_lb[k] += m_pe_array[i]->pes[j]->access_energy_lb[k];

                // Update transfer cost between local buffer and computing unit
                transfer_cycle_pe[k] = std::max(transfer_cycle_pe[k], m_pe_array[i]->pes[j]->transfer_cycle[k]);
                pe_link_pe_types += m_pe_array[i]->pes[j]->transfer_cycle[k];
                transfer_energy_pe[k] += m_pe_array[i]->pes[j]->transfer_energy[k];
                payload_link_transactions_pe[k] += m_pe_array[i]->pes[j]->payload_link_transactions[k];
                metadata_link_transactions_pe[k] += m_pe_array[i]->pes[j]->metadata_link_transactions[k];
                storage_link_transactions_pe[k] += m_pe_array[i]->pes[j]->storage_link_transactions[k];

                format_cycle_pe[k] = std::max(format_cycle_pe[k], m_pe_array[i]->pes[j]->format_cycle[k]);
                pe_format_types += m_pe_array[i]->pes[j]->format_cycle[k];
                // E20-2: sum the Format-IP transaction counts over every PE, so an active but
                // unpriced format path is visible even when its energy axis reads 0.
                format_payload_events += m_pe_array[i]->pes[j]->format_payload_events[k];
                format_metadata_events += m_pe_array[i]->pes[j]->format_metadata_events[k];
                format_energy_pe[k] += m_pe_array[i]->pes[j]->format_energy[k];

                // Update overlapped cycle between local buffer and computing unit
                cycle_mac_lb[k] = std::max(cycle_mac_lb[k], m_pe_array[i]->pes[j]->cycle_mac_lb[k]);
                pe_overlap_mac_lb_types += m_pe_array[i]->pes[j]->cycle_mac_lb[k];

                // Update local buffer utilization
                utilization_local_buffer[k] = std::max(utilization_local_buffer[k], m_pe_array[i]->pes[j]->utilization_local_buffer[k]);
            }

            local_buffer_type = m_pe_array[i]->pes[j]->get_memory_type();
            // E3: every PE reads the same config, so record the compute-energy basis once, with the
            // operand precision it applies to.
            compute_energy_basis = m_pe_array[i]->pes[j]->compute_energy_basis;
            compute_energy_precision_calibrated =
                m_pe_array[i]->pes[j]->compute_energy_precision_calibrated;
            // RE4: the precision a MAC's energy depends on is all THREE widths, so the label
            // names the accumulator too. Reporting only the operands made an INT8 x INT8 -> FP32
            // MAC and an INT8 x INT8 -> FP16 MAC look like the same calibrated case.
            operand_precision = runtime_datatypes().format(data_type_t::INPUT).name + " x " +
                                runtime_datatypes().format(data_type_t::WEIGHT).name + " -> " +
                                runtime_datatypes().accumulator_format().name;
            pe_double_buffer = pe_double_buffer && m_pe_array[i]->pes[j]->double_buffer;

            // CE4: reduce THIS PE's own datatype-combined value into the running
            // max across PE entities -- max_entity(sum/max_type(...)), not the
            // sum/max_type(max_entity(...)) the flat vectors above would give.
            const double pe_combined_access_lb = (local_buffer_type == memory_type_t::SHARED)
                ? pe_access_lb_sum_types : pe_access_lb_max_type;
            entity_combined_access_lb = std::max(entity_combined_access_lb, pe_combined_access_lb);
            entity_combined_link_pe = std::max(entity_combined_link_pe, pe_link_pe_types);
            entity_combined_overlap_mac_lb = std::max(entity_combined_overlap_mac_lb, pe_overlap_mac_lb_types);
            entity_combined_format = std::max(entity_combined_format, pe_format_types);

            num_active_pe++;
            first_active_pe = false;
        }

        // CE4: this chip's own datatype combination, reduced across chip entities
        // (below the loop) in the correct max_entity(sum/max_type(...)) order.
        double chip_access_pe_array_types = 0.0;
        double chip_link_pe_array_types = 0.0;
        double chip_overlap_pe_array_types = 0.0;
        // CE4/P1-D/L1: this chip's own per-datatype base and fill, kept per chip AND per
        // datatype so both the entity dimension and the datatype correspondence survive
        // the (differing) repetition scaling that follows.
        std::vector<double> chip_base_access(data_type_t::NUM_DATA_TYPES, 0.0);
        std::vector<double> chip_fill_access(data_type_t::NUM_DATA_TYPES, 0.0);
        double chip_link_global_buffer_types = 0.0;
        double chip_overlap_global_buffer_types = 0.0;
        for(unsigned j = 0; j < data_type_t::NUM_DATA_TYPES; j++) {
            /* Update PE array stats */

            // Update the number of request to PE array
            num_request_pe_array[j] += m_pe_array[i]->num_request[j];

            // Update the number of data transfer from PE array
            num_data_transfer_pe_array[j] += m_pe_array[i]->num_data_transfer[j];
            payload_link_transactions_pe_array[j] += m_pe_array[i]->payload_link_transactions[j];
            metadata_link_transactions_pe_array[j] += m_pe_array[i]->metadata_link_transactions[j];
            storage_link_transactions_pe_array[j] += m_pe_array[i]->storage_link_transactions[j];

            // Update access cost of PE array
            access_cycle_pe_array[j] = std::max(access_cycle_pe_array[j], m_pe_array[i]->access_cycle[j]);
            chip_access_pe_array_types += m_pe_array[i]->access_cycle[j];
            access_energy_pe_array[j] += m_pe_array[i]->access_energy[j];
            if(j == data_type_t::OUTPUT) {
                psum_store_access_cycle_pe_array = std::max(psum_store_access_cycle_pe_array,
                                                            m_pe_array[i]->psum_store_access_cycle);
                psum_store_access_energy_pe_array += m_pe_array[i]->psum_store_access_energy;
            }

            // Update transfer cycle from PE array to local buffers
            transfer_cycle_pe_array[j] = std::max(transfer_cycle_pe_array[j], m_pe_array[i]->transfer_cycle[j]);
            // CE4 rule per the declared fabric: one shared distribution medium serializes the
            // datatype streams (sum); separate per-datatype buses run them concurrently (max).
            // Mirrors the GLB shared/separate port rule.
            if(m_pe_array[i]->array_fabric_separate) {
                chip_link_pe_array_types = std::max(chip_link_pe_array_types,
                                                    m_pe_array[i]->transfer_cycle[j]);
            } else {
                chip_link_pe_array_types += m_pe_array[i]->transfer_cycle[j];
            }
            cycle_temporal_pe_array[j] = std::max(cycle_temporal_pe_array[j], m_pe_array[i]->cycle_temporal_pe[j]);
            chip_overlap_pe_array_types += m_pe_array[i]->cycle_temporal_pe[j];
            transfer_energy_pe_array[j] += m_pe_array[i]->transfer_energy[j];

            /* Update global buffer stats */

            // Update global buffer type
            global_buffer_type = m_global_buffer[i]->get_memory_type();
            // P1-B: record the direct-stream (GLB-storage-bypassed) datatypes.
            if(m_global_buffer[i]->bypass[j]) global_buffer_bypass[j] = true;

            // Update the number of request to the global buffer
            num_request_global_buffer[j] += m_global_buffer[i]->num_request[j];
            payload_link_transactions_global_buffer[j] += m_global_buffer[i]->payload_link_transactions[j];
            metadata_link_transactions_global_buffer[j] += m_global_buffer[i]->metadata_link_transactions[j];
            storage_link_transactions_global_buffer[j] += m_global_buffer[i]->storage_link_transactions[j];
            
            // Update the number of data transfer from the global buffer
            num_data_transfer_global_buffer[j] += m_global_buffer[i]->num_data_transfer[j];
            if(j == data_type_t::OUTPUT) psum_writeback_events_global_buffer += m_global_buffer[i]->psum_writeback_events;

            // Update access cost to the global buffer
            access_cycle_global_buffer[j] = std::max(access_cycle_global_buffer[j], m_global_buffer[i]->access_cycle[j]);
            access_energy_global_buffer[j] += m_global_buffer[i]->access_energy[j];
            fill_access_cycle_global_buffer[j] = std::max(fill_access_cycle_global_buffer[j], m_global_buffer[i]->fill_access_cycle[j]);
            fill_access_energy_global_buffer[j] += m_global_buffer[i]->fill_access_energy[j];
            chip_fill_access[j] = m_global_buffer[i]->fill_access_cycle[j];
            // CE4/P1-D/L1: this chip's BASE access, per datatype (fill excluded -- it
            // scales per-datatype, not uniformly, so the two sides are combined in
            // merge_global_buffer_fill() after each has been scaled by its own repetition
            // factor; mixing them here before scaling would apply the wrong (uniform)
            // factor to the fill).
            chip_base_access[j] = m_global_buffer[i]->access_cycle[j];

            // Update transfer cost between the global buffer to the PE array
            transfer_cycle_global_buffer[j] = std::max(transfer_cycle_global_buffer[j], m_global_buffer[i]->transfer_cycle[j]);
            { static int d=0; if(d<3) std::cerr << "[GLBF] sep=" << m_global_buffer[i]->fabric_separate << " tc=" << m_global_buffer[i]->transfer_cycle[j] << std::endl; d++; }
            // CE4 rule per the declared GLB fabric: shared medium serializes the datatype
            // streams (sum); separate per-datatype buses run them concurrently (max).
            if(m_global_buffer[i]->fabric_separate) {
                chip_link_global_buffer_types = std::max(chip_link_global_buffer_types,
                                                         m_global_buffer[i]->transfer_cycle[j]);
            } else {
                chip_link_global_buffer_types += m_global_buffer[i]->transfer_cycle[j];
            }
            transfer_energy_global_buffer[j] += m_global_buffer[i]->transfer_energy[j];

            // Update overlapped cycle between the global buffer and PE array
            cycle_pe_array_global_buffer[j] = std::max(cycle_pe_array_global_buffer[j], m_global_buffer[i]->cycle_pe_array_global_buffer[j]);
            // Same fabric rule for the per-transfer pipeline (overlap) axis: separate
            // per-datatype buses fill and stream independently, so they combine by MAX.
            if(m_global_buffer[i]->fabric_separate) {
                chip_overlap_global_buffer_types = std::max(chip_overlap_global_buffer_types,
                    m_global_buffer[i]->cycle_pe_array_global_buffer[j]);
            } else {
                chip_overlap_global_buffer_types += m_global_buffer[i]->cycle_pe_array_global_buffer[j];
            }

            // Update global buffer utilization
            utilization_global_buffer[j] = std::max(utilization_global_buffer[j], m_global_buffer[i]->utilization[j]);
        }

        // CE4: reduce THIS chip's own datatype-combined value into the running max
        // across chip entities -- max_entity(sum/max_type(...)).
        entity_combined_access_pe_array = std::max(entity_combined_access_pe_array, chip_access_pe_array_types);
        entity_combined_link_pe_array = std::max(entity_combined_link_pe_array, chip_link_pe_array_types);
        entity_combined_overlap_pe_array = std::max(entity_combined_overlap_pe_array, chip_overlap_pe_array_types);
        // CE4/P1-D/L1: record this chip's own per-datatype base and fill. Datatype and
        // entity combination both happen in merge_global_buffer_fill(), after each side
        // has been scaled by its own repetition factor.
        chip_access_cycle_global_buffer.push_back(chip_base_access);
        chip_fill_access_cycle_global_buffer.push_back(chip_fill_access);
        entity_combined_link_global_buffer = std::max(entity_combined_link_global_buffer, chip_link_global_buffer_types);
        entity_combined_overlap_global_buffer = std::max(entity_combined_overlap_global_buffer, chip_overlap_global_buffer_types);

        // Update PE array utilization
        utilization_pe_array = std::max(utilization_pe_array, m_pe_array[i]->utilization);
        for(unsigned j = 0; j < data_type_t::NUM_DATA_TYPES; j++) {
            utilization_pe_array_buffer[j] = std::max(utilization_pe_array_buffer[j], m_pe_array[i]->buffer_utilization[j]);
        }
        // V2: fold fill serializes on the compute schedule (max across parallel arrays).
        fold_fill_cycle_pe_array = std::max(fold_fill_cycle_pe_array, m_pe_array[i]->fold_fill_cycle);
        layer_setup_cycle_pe_array = std::max(layer_setup_cycle_pe_array, m_pe_array[i]->u_layer_setup_cycle);
        // Per-stripe schedule bubble: one-time total = unit bubble x legacy C*R*S*P count.
        // Serializes on the compute schedule like fold fill / setup (max across arrays).
        stripe_transition_cycle_pe_array = std::max(stripe_transition_cycle_pe_array,
            m_pe_array[i]->u_stripe_transition_cycle
            * static_cast<double>(m_pe_array[i]->stripe_transitions));
        // E5: fold energy accumulates across arrays like the events that produced it; the
        // per-layer schedule setup fires once per array.
        weight_fold_energy += m_pe_array[i]->weight_fold_energy;
        pe_array_accumulator_energy += m_pe_array[i]->accumulator_energy;
        weight_fold_events += m_pe_array[i]->weight_fold_events;
        reduction_additions += m_pe_array[i]->reduction_additions;
        // RE3: the setup EVENT comes from the setup actually executing -- a declared setup cycle
        // cost -- not from its energy being priced. Keying the event off the energy made an
        // uncalibrated setup report "0.00 pJ over 0 setup event(s)", which reads as "no setup
        // happened" when the schedule in fact pays 2270 cycles for one. The converse (energy with
        // no setup execution) is rejected at config load; see validate_energy_settings().
        if(m_pe_array[i]->u_layer_setup_cycle > 0.0) {
            layer_setup_energy += m_pe_array[i]->u_layer_setup_energy;
            layer_setup_events += 1.0;
        }
        // RE3 convention as above: the transition EVENTS exist when the model models them (a
        // declared cycle cost), not when their energy is priced.
        if(m_pe_array[i]->u_stripe_transition_cycle > 0.0 && m_pe_array[i]->stripe_transitions > 0) {
            stripe_transition_events += static_cast<double>(m_pe_array[i]->stripe_transitions);
            stripe_transition_energy += m_pe_array[i]->u_stripe_transition_energy
                                      * static_cast<double>(m_pe_array[i]->stripe_transitions);
        }
    }

    // Guard the per-PE averages against a mapping that activates zero PEs (avoids NaN).
    if(num_active_pe > 0) {
        avg_computation_cycle /= num_active_pe;
        for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
            avg_access_cycle_mac[i] /= num_active_pe;
            avg_access_cycle_lb[i] /= num_active_pe;
        }
    }

    // Available MAC cycles must span the whole layer wall-clock (compute + memory stalls),
    // the same window used for leakage; otherwise memory-bound layers report near-full MAC
    // utilization while the MACs actually idle waiting on GLB/NoP/DRAM.
    mac_available_cycle = static_cast<double>(physical_scalar_macs) * layer_latency;
    utilization_mac = calculate_time_based_mac_utilization(mac_busy_cycle, mac_available_cycle);
    if(utilization_mac > 1.0 + 1e-9) {
        std::cerr << "Error: MAC busy cycles exceed physical MAC capacity" << std::endl;
        exit(1);
    }
    utilization_mac = std::min(1.0, utilization_mac);

    // Update PE static energy (leakage) over all physical chips and all physical PEs
    // (always-on): inactive chips/PEs still leak for the layer duration.
    for(unsigned i = 0; i < m_multi_chip->get_number_of_chips(); i++) {
        for(unsigned j = 0; j < m_pe_array[i]->get_number_of_pes(); j++) {
            for(unsigned k = 0; k < data_type_t::NUM_DATA_TYPES; k++) {
                static_energy_pe[k] += m_pe_array[i]->pes[j]->static_energy[k];
            }
        }
    }

    // Update PE-array temporal-buffer static energy over all physical chips (always-on).
    for(unsigned i = 0; i < m_multi_chip->get_number_of_chips(); i++) {
        for(unsigned j = 0; j < data_type_t::NUM_DATA_TYPES; j++) {
            static_energy_pe_array[j] += m_pe_array[i]->static_energy[j];
        }
    }

    // Update global buffer static energy over all physical chips (always-on).
    for(unsigned i = 0; i < m_multi_chip->get_number_of_chips(); i++) {
        for(unsigned j = 0; j < data_type_t::NUM_DATA_TYPES; j++) {
            static_energy_global_buffer[j] += m_global_buffer[i]->static_energy[j];
        }
    }

    // Update stats of Multi chip and DRAM.
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        /* Update chip-level processors' stats */

        // Update the number of request to the chip-level processor
        num_request_multi_chip[i] = m_multi_chip->num_request[i];

        // Update the number of data transfer from the chip-level processor
        num_data_transfer_multi_chip[i] = m_multi_chip->num_data_transfer[i];
        payload_link_transactions_multi_chip[i] = m_multi_chip->payload_link_transactions[i];
        metadata_link_transactions_multi_chip[i] = m_multi_chip->metadata_link_transactions[i];
        storage_link_transactions_multi_chip[i] = m_multi_chip->storage_link_transactions[i];

        // Update access cost of the chip-level processor
        access_cycle_multi_chip[i] = m_multi_chip->access_cycle[i];
        access_energy_multi_chip[i] = m_multi_chip->access_energy[i];

        // Update transfer cost between the global buffer and chip-level processor
        transfer_cycle_multi_chip[i] = m_multi_chip->transfer_cycle[i];
        transfer_energy_multi_chip[i] = m_multi_chip->transfer_energy[i];

        // Update static (leakage) energy of the Multi-chip temporal buffer
        static_energy_multi_chip[i] = m_multi_chip->static_energy[i];

        /* Update off-chip memory stats */

        // Update the number of request to the off-chip memory
        num_request_dram[i] = m_dram->num_request[i];

        // Update the number of data transfer from the off-chip memory
        num_data_transfer_dram[i] = m_dram->num_data_transfer[i];
        payload_link_transactions_dram[i] = m_dram->payload_link_transactions[i];
        metadata_link_transactions_dram[i] = m_dram->metadata_link_transactions[i];
        storage_link_transactions_dram[i] = m_dram->storage_link_transactions[i];

        // Update access cost of the off-chip memory
        access_cycle_dram[i] = m_dram->access_cycle[i];
        access_energy_dram[i] = m_dram->access_energy[i];

        // Update transfer cost between the off-chip memory and chip-level processor
        transfer_cycle_dram[i] = m_dram->transfer_cycle[i];
        transfer_energy_dram[i] = m_dram->transfer_energy[i];

        // Update overlapped cycle between the off-chip memory and chip-level processor
        cycle_chip_dram[i] = std::max(cycle_chip_dram[i], m_dram->cycle_chip_dram[i]);

    }
    
    // Update utilization of chip-level processor
    utilization_multi_chip = std::max(utilization_multi_chip, m_multi_chip->utilization);
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        utilization_multi_chip_buffer[i] = std::max(utilization_multi_chip_buffer[i],
                                                   m_multi_chip->buffer_utilization[i]);
    }

}

namespace {

unsigned scale_counter(unsigned value, unsigned repetitions, const char *label) {
    if(value != 0 && repetitions > std::numeric_limits<unsigned>::max()/value) {
        std::cerr << "Error: unsigned overflow while scaling " << label << std::endl;
        exit(1);
    }
    return value * repetitions;
}

size_t scale_counter(size_t value, unsigned repetitions, const char *label) {
    if(value != 0 && repetitions > std::numeric_limits<size_t>::max()/value) {
        std::cerr << "Error: size overflow while scaling " << label << std::endl;
        exit(1);
    }
    return value * repetitions;
}

void scale_counters(std::vector<size_t> *values, unsigned repetitions, const char *label) {
    for(unsigned i = 0; i < values->size(); ++i) {
        values->at(i) = scale_counter(values->at(i), repetitions, label);
    }
}

void scale_counters(std::vector<unsigned> *values, unsigned repetitions, const char *label) {
    for(unsigned i = 0; i < values->size(); ++i) {
        values->at(i) = scale_counter(values->at(i), repetitions, label);
    }
}

void scale_costs(std::vector<double> *values, unsigned repetitions) {
    for(unsigned i = 0; i < values->size(); ++i) {
        values->at(i) *= repetitions;
    }
}
size_t rescale_counter_ratio(size_t value, size_t numerator, size_t denominator,
                             const char *label) {
    if(value == 0 || numerator == denominator) return value;
    if(denominator == 0) {
        std::cerr << "Error: zero denominator while scaling " << label << std::endl;
        exit(1);
    }
    const long double scaled = static_cast<long double>(value)*numerator/denominator;
    if(scaled > static_cast<long double>(std::numeric_limits<size_t>::max())) {
        std::cerr << "Error: size overflow while scaling " << label << std::endl;
        exit(1);
    }
    return static_cast<size_t>(scaled + 0.5L);

}
} // namespace

void stats_t::scale_serial_repetitions(unsigned m_repetitions,
                                       const std::vector<unsigned> &m_datatype_repetitions,
                                       const input_halo_reuse_t &m_input_halo,
                                       bool m_halo_capacity_sufficient) {
    if(m_repetitions == 0) {
        std::cerr << "Error: temporal repetition count must be non-zero" << std::endl;
        exit(1);
    }
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        if(m_datatype_repetitions[i] == 0 || m_datatype_repetitions[i] > m_repetitions) {
            std::cerr << "Error: per-datatype repetition factor out of range" << std::endl;
            exit(1);
        }
    }
    input_halo_overlap = m_input_halo.active;
    input_halo_capacity_sufficient = m_input_halo.active && m_halo_capacity_sufficient;
    input_halo_unique_elements = m_input_halo.unique_elements;
    input_halo_replicated_elements = m_input_halo.replicated_elements;
    input_halo_working_set_bytes = runtime_datatypes().storage_bytes(
        data_type_t::INPUT, m_input_halo.working_set_elements);
    input_halo_reuse_applied = false;
    input_halo_pre_dram_transactions = 0;

    // P3: record the tile count for finalize_layer_timeline()'s fill+bottleneck
    // combination, covering both branches below.
    temporal_repetition_tiles = m_repetitions;
    // L5: and the off-chip refetch count on that same tile clock (see
    // stats_t::offchip_repetition_tiles).
    offchip_repetition_tiles = 1;
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        offchip_repetition_tiles = std::max(offchip_repetition_tiles, m_datatype_repetitions[i]);
    }
    // Phase-5: how many tiles of the shared clock COMMIT final outputs -- what the
    // streaming SFU stage consumes (see stats_t::output_repetition_tiles).
    output_repetition_tiles = std::min(m_repetitions,
                                       std::max(1u, m_datatype_repetitions[data_type_t::OUTPUT]));
    if(m_repetitions == 1) {
        // V2: the one-time per-layer setup applies regardless of repetition count.
        fold_fill_cycle_pe_array += layer_setup_cycle_pe_array + stripe_transition_cycle_pe_array;
        merge_global_buffer_fill();
        finalize_layer_timeline();
        return;
    }

    num_computation = scale_counter(num_computation, m_repetitions, "computation count");
    layer_latency *= m_repetitions;
    busy_cycle_pe *= m_repetitions;
    busy_cycle_pe_array *= m_repetitions;
    busy_cycle_global_buffer *= m_repetitions;
    busy_cycle_multi_chip *= m_repetitions;
    busy_cycle_dram *= m_repetitions;
    computation_cycle *= m_repetitions;
    max_computation_cycle *= m_repetitions;
    min_computation_cycle *= m_repetitions;
    avg_computation_cycle *= m_repetitions;
    mac_busy_cycle *= m_repetitions;
    mac_available_cycle *= m_repetitions;
    computation_energy *= m_repetitions;
    reduction_energy *= m_repetitions;
    accumulator_reload_bytes = scale_counter(accumulator_reload_bytes, m_repetitions,
                                             "accumulator reload bytes");
    accumulator_spill_bytes = scale_counter(accumulator_spill_bytes, m_repetitions,
                                                   "accumulator spill bytes");
    accumulator_create_events = scale_counter(accumulator_create_events, m_repetitions,
                                              "accumulator create events");
    accumulator_energy *= m_repetitions;

    scale_counters(&num_request_pe, m_repetitions, "PE request count");
    scale_counters(&num_data_transfer_pe, m_repetitions, "PE transfer count");
    scale_costs(&format_cycle_pe, m_repetitions);
    scale_costs(&format_energy_pe, m_repetitions);
    entity_combined_format *= m_repetitions;
    scale_costs(&access_cycle_mac, m_repetitions);
    scale_costs(&max_access_cycle_mac, m_repetitions);
    scale_costs(&min_access_cycle_mac, m_repetitions);
    scale_costs(&avg_access_cycle_mac, m_repetitions);
    scale_costs(&access_energy_mac, m_repetitions);
    scale_costs(&access_cycle_lb, m_repetitions);
    scale_costs(&max_access_cycle_lb, m_repetitions);
    scale_costs(&min_access_cycle_lb, m_repetitions);
    scale_costs(&avg_access_cycle_lb, m_repetitions);
    scale_costs(&access_energy_lb, m_repetitions);
    scale_costs(&transfer_cycle_pe, m_repetitions);
    scale_costs(&transfer_energy_pe, m_repetitions);
    scale_counters(&payload_link_transactions_pe, m_repetitions, "PE payload transactions");
    scale_counters(&metadata_link_transactions_pe, m_repetitions, "PE metadata transactions");
    scale_counters(&storage_link_transactions_pe, m_repetitions, "PE storage transactions");
    scale_costs(&cycle_mac_lb, m_repetitions);
    scale_costs(&static_energy_pe, m_repetitions);
    scale_costs(&static_energy_pe_array, m_repetitions);
    // CE4: the entity-combined scalars scale uniformly like their flat-vector
    // siblings above (this axis has no per-datatype repetition factor).
    entity_combined_access_lb *= m_repetitions;
    entity_combined_link_pe *= m_repetitions;
    entity_combined_overlap_mac_lb *= m_repetitions;

    scale_counters(&num_request_pe_array, m_repetitions, "PE-array request count");
    scale_counters(&num_data_transfer_pe_array, m_repetitions, "PE-array transfer count");
    scale_counters(&payload_link_transactions_pe_array, m_repetitions, "PE-array payload transactions");
    scale_counters(&metadata_link_transactions_pe_array, m_repetitions, "PE-array metadata transactions");
    scale_counters(&storage_link_transactions_pe_array, m_repetitions, "PE-array storage transactions");
    // E20-3b: capture the PE-array OUTPUT column BEFORE the uniform scaling -- the GLB<->array
    // output transfer may need a different factor (see below), and both of its endpoints must
    // carry the same one.
    const double base_output_array_store_cycle = psum_store_access_cycle_pe_array;
    const double base_output_array_store_energy = psum_store_access_energy_pe_array;
    const double base_output_array_other_cycle =
        access_cycle_pe_array[data_type_t::OUTPUT] - base_output_array_store_cycle;
    const double base_output_array_other_energy =
        access_energy_pe_array[data_type_t::OUTPUT] - base_output_array_store_energy;
    scale_costs(&access_cycle_pe_array, m_repetitions);
    scale_costs(&access_energy_pe_array, m_repetitions);
    scale_costs(&transfer_cycle_pe_array, m_repetitions);
    scale_costs(&cycle_temporal_pe_array, m_repetitions);
    scale_costs(&transfer_energy_pe_array, m_repetitions);
    entity_combined_access_pe_array *= m_repetitions;
    entity_combined_link_pe_array *= m_repetitions;
    entity_combined_overlap_pe_array *= m_repetitions;
    // V2: fold fill repeats with every global-buffer repetition; the per-layer setup
    // (config/flush/DMA prologue) is a one-time cost added after scaling.
    fold_fill_cycle_pe_array = fold_fill_cycle_pe_array*m_repetitions + layer_setup_cycle_pe_array
                             + stripe_transition_cycle_pe_array;
    // E5: the same split on the energy side -- per-fold energy repeats with every GLB
    // repetition, the one-time layer setup does not.
    weight_fold_energy *= m_repetitions;
    pe_array_accumulator_energy *= m_repetitions;
    weight_fold_events *= m_repetitions;

    // E20-3b: the GLB<->array OUTPUT leg does not always repeat with the uniform GLB factor.
    // When no psum has to leave the array between reduction steps, what crosses this boundary is
    // the FINISHED tile, once per distinct output tile -- so it scales with the OUTPUT datatype's
    // own repetition factor, exactly the argument RE1 already applies to the off-chip output cast
    // ("a GLB repetition over a reduction dimension re-accumulates the SAME output rather than
    // producing a new one"). Scaling it uniformly multiplied the write-out by the reduction depth
    // (8x on the 512x512x512 GEMM). When psums DO leave, every repetition genuinely moves one and
    // the uniform factor is correct, which is every validated CONV layer.
    const unsigned output_repetitions = psum_array_resident
        ? m_datatype_repetitions[data_type_t::OUTPUT] : m_repetitions;
    const size_t base_output_transfer = num_data_transfer_global_buffer[data_type_t::OUTPUT];
    const size_t base_output_stores = psum_writeback_events_global_buffer;
    const size_t base_output_payload = payload_link_transactions_global_buffer[data_type_t::OUTPUT];
    const size_t base_output_metadata = metadata_link_transactions_global_buffer[data_type_t::OUTPUT];
    const size_t base_output_storage = storage_link_transactions_global_buffer[data_type_t::OUTPUT];
    const double base_output_access_cycle = access_cycle_global_buffer[data_type_t::OUTPUT];
    const double base_output_access_energy = access_energy_global_buffer[data_type_t::OUTPUT];
    const double base_output_pipeline = cycle_pe_array_global_buffer[data_type_t::OUTPUT];
    const double base_output_transfer_cycle = transfer_cycle_global_buffer[data_type_t::OUTPUT];
    const double base_output_transfer_energy = transfer_energy_global_buffer[data_type_t::OUTPUT];
    std::vector<double> base_output_chip_access(chip_access_cycle_global_buffer.size(), 0.0);
    for(unsigned chip = 0; chip < chip_access_cycle_global_buffer.size(); ++chip) {
        if(data_type_t::OUTPUT < chip_access_cycle_global_buffer[chip].size()) {
            base_output_chip_access[chip] = chip_access_cycle_global_buffer[chip][data_type_t::OUTPUT];
        }
    }

    scale_counters(&num_request_global_buffer, m_repetitions, "global-buffer request count");
    scale_counters(&num_data_transfer_global_buffer, m_repetitions, "global-buffer transfer count");
    psum_writeback_events_global_buffer *= m_repetitions;
    scale_counters(&payload_link_transactions_global_buffer, m_repetitions, "global-buffer payload transactions");
    scale_counters(&metadata_link_transactions_global_buffer, m_repetitions, "global-buffer metadata transactions");
    scale_counters(&storage_link_transactions_global_buffer, m_repetitions, "global-buffer storage transactions");
    scale_costs(&access_cycle_global_buffer, m_repetitions);
    scale_costs(&access_energy_global_buffer, m_repetitions);
    scale_costs(&cycle_pe_array_global_buffer, m_repetitions);
    scale_costs(&transfer_cycle_global_buffer, m_repetitions);
    scale_costs(&transfer_energy_global_buffer, m_repetitions);
    scale_costs(&static_energy_global_buffer, m_repetitions);
    // CE4/P1-D/L1: the per-chip/per-datatype BASE access scales uniformly
    // (entity_combined_access_global_buffer is recomputed from these and the fill matrix
    // in merge_global_buffer_fill(), once both sides are final).
    for(unsigned chip = 0; chip < chip_access_cycle_global_buffer.size(); ++chip) {
        scale_costs(&chip_access_cycle_global_buffer[chip], m_repetitions);
    }
    entity_combined_link_global_buffer *= m_repetitions;
    entity_combined_overlap_global_buffer *= m_repetitions;
    if(output_repetitions != m_repetitions) {
        num_data_transfer_global_buffer[data_type_t::OUTPUT] =
            scale_counter(base_output_transfer, output_repetitions, "global-buffer output transfer count");
        psum_writeback_events_global_buffer = base_output_stores*output_repetitions;
        payload_link_transactions_global_buffer[data_type_t::OUTPUT] =
            scale_counter(base_output_payload, output_repetitions, "global-buffer output payload transactions");
        metadata_link_transactions_global_buffer[data_type_t::OUTPUT] =
            scale_counter(base_output_metadata, output_repetitions, "global-buffer output metadata transactions");
        storage_link_transactions_global_buffer[data_type_t::OUTPUT] =
            scale_counter(base_output_storage, output_repetitions, "global-buffer output storage transactions");
        access_cycle_global_buffer[data_type_t::OUTPUT] = base_output_access_cycle*output_repetitions;
        access_energy_global_buffer[data_type_t::OUTPUT] = base_output_access_energy*output_repetitions;
        cycle_pe_array_global_buffer[data_type_t::OUTPUT] = base_output_pipeline*output_repetitions;
        transfer_cycle_global_buffer[data_type_t::OUTPUT] = base_output_transfer_cycle*output_repetitions;
        transfer_energy_global_buffer[data_type_t::OUTPUT] = base_output_transfer_energy*output_repetitions;
        // The PE-array side of this same transfer -- reading the finished tile OUT of the array's
        // temporal buffer -- has to carry the identical factor. Scaling the two endpoints of one
        // transfer differently is not a smaller error than scaling both wrong; array_buffer AB2
        // caught exactly that (a 4x read term against a 1x write term on the GLB side).
        access_cycle_pe_array[data_type_t::OUTPUT] =
            base_output_array_other_cycle*m_repetitions + base_output_array_store_cycle*output_repetitions;
        access_energy_pe_array[data_type_t::OUTPUT] =
            base_output_array_other_energy*m_repetitions + base_output_array_store_energy*output_repetitions;
        psum_store_access_cycle_pe_array = base_output_array_store_cycle*output_repetitions;
        psum_store_access_energy_pe_array = base_output_array_store_energy*output_repetitions;
        // The per-chip access matrix feeds the GLB access AXIS through merge_global_buffer_fill(),
        // so it has to carry the same factor -- otherwise the axis keeps the uniform scale while
        // the per-datatype rows it is supposed to combine do not, which T10 catches.
        for(unsigned chip = 0; chip < chip_access_cycle_global_buffer.size(); ++chip) {
            if(data_type_t::OUTPUT < chip_access_cycle_global_buffer[chip].size()) {
                chip_access_cycle_global_buffer[chip][data_type_t::OUTPUT] =
                    base_output_chip_access[chip]*output_repetitions;
            }
        }
    }

    // Off-chip traffic scales per datatype: a GLB repetition over a dimension the
    // tensor does not depend on (e.g. the Q loop for weights) revisits the SAME
    // tile, which the multi-chip/GLB buffers retain -- only the repetitions that
    // actually index the tensor re-fetch it. Requests and leakage remain uniform.
    scale_counters(&num_request_multi_chip, m_repetitions, "multi-chip request count");
    scale_counters(&num_request_dram, m_repetitions, "DRAM request count");
    scale_costs(&static_energy_multi_chip, m_repetitions);
    // Weight decompression: the decoder decodes exactly what the weight DRAM stream
    // delivers, so every weight REFETCH (the WEIGHT datatype's off-chip repetition factor
    // below) is decoded again -- decoder work, decode energy, tile count, the scratchpad
    // sink and the reported DRAM saving all carry that factor. The dense/compressed BYTE
    // fields stay unscaled: they are the layer's storage footprint, not traffic. The
    // relative-throughput mode (decomp_decoder_ratio) is exempt because its decoder window
    // is derived in finalize_layer_timeline() from the already-scaled compute-side busy.
    if(decomp_present && data_type_t::WEIGHT < m_datatype_repetitions.size()) {
        const unsigned weight_refetch =
            std::max(1u, m_datatype_repetitions[data_type_t::WEIGHT]);
        if(weight_refetch > 1) {
            if(decomp_decoder_ratio <= 0.0) {
                decomp_decoder_cycles *= weight_refetch;
            }
            decomp_decoder_energy *= weight_refetch;
            decomp_tiles *= weight_refetch;
            decomp_sink_cycles *= weight_refetch;
            decomp_dram_weight_saved_cycle *= weight_refetch;
        }
    }
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        const unsigned repetitions = m_datatype_repetitions[i];
        num_data_transfer_multi_chip[i] = scale_counter(num_data_transfer_multi_chip[i], repetitions, "multi-chip transfer count");
        payload_link_transactions_multi_chip[i] = scale_counter(payload_link_transactions_multi_chip[i], repetitions, "multi-chip payload transactions");
        metadata_link_transactions_multi_chip[i] = scale_counter(metadata_link_transactions_multi_chip[i], repetitions, "multi-chip metadata transactions");
        storage_link_transactions_multi_chip[i] = scale_counter(storage_link_transactions_multi_chip[i], repetitions, "multi-chip storage transactions");
        access_cycle_multi_chip[i] *= repetitions;
        access_energy_multi_chip[i] *= repetitions;
        transfer_cycle_multi_chip[i] *= repetitions;
        transfer_energy_multi_chip[i] *= repetitions;

        num_data_transfer_dram[i] = scale_counter(num_data_transfer_dram[i], repetitions, "DRAM transfer count");
        payload_link_transactions_dram[i] = scale_counter(payload_link_transactions_dram[i], repetitions, "DRAM payload transactions");
        metadata_link_transactions_dram[i] = scale_counter(metadata_link_transactions_dram[i], repetitions, "DRAM metadata transactions");
        storage_link_transactions_dram[i] = scale_counter(storage_link_transactions_dram[i], repetitions, "DRAM storage transactions");
        access_cycle_dram[i] *= repetitions;
        access_energy_dram[i] *= repetitions;
        cycle_chip_dram[i] *= repetitions;
        transfer_cycle_dram[i] *= repetitions;
        transfer_energy_dram[i] *= repetitions;

        // RE1: the final output cast is committed at the off-chip store, so it scales with the
        // OUTPUT datatype's repetition factor -- not the uniform GLB repetition. A GLB repetition
        // over a reduction dimension (the C loop here) re-accumulates the SAME output rather than
        // producing a new one, so scaling the cast by it would multiply the final element count by
        // the reduction depth (4x on the GEMM fixture).
        if(i == data_type_t::OUTPUT) {
            output_cast_bytes = scale_counter(output_cast_bytes, repetitions, "output cast bytes");
            output_cast_energy *= repetitions;
            output_cast_cycle *= repetitions;
        }
        // GLB fill (write) side mirrors the off-chip supply, so it scales with the
        // same per-datatype factor before folding into the GLB access totals.
        fill_access_cycle_global_buffer[i] *= repetitions;
        fill_access_energy_global_buffer[i] *= repetitions;
        // CE4/P1-D: same factor on the per-chip copies, which keep the entity dimension.
        for(unsigned chip = 0; chip < chip_fill_access_cycle_global_buffer.size(); ++chip) {
            chip_fill_access_cycle_global_buffer[chip][i] *= repetitions;
        }
    }
    // E20-3: the live descriptor pass charged every legacy-GLB P/Q window independently.
    // When the GLB can hold the ring-buffer working set, keep the logical requests but coalesce
    // their overlapping dense payload to the exact union footprint. Cycle and energy on every
    // physical off-chip side use the same ratio, so traffic cannot change without its costs.
    //
    // EXTENDED 2026-08-26: the correction runs in BOTH directions, with an asymmetric guard.
    // Coalescing DOWN (windows overlap, per-repetition streaming over-fetches) is a retention
    // claim, so it stays gated on the working set fitting the GLB. Correcting UP happens when a
    // legacy-row FILTER loop leaves the below-legacy tile without its kernel extent, so the
    // repetition-scaled charge falls SHORT of the input union -- and the union is a hard lower
    // bound on off-chip input traffic under ANY loop order and ANY buffer size (every byte a
    // window touches crosses the interface at least once), so that direction is unconditional.
    if(input_halo_overlap && dram_input_link_bits > 0) {
        const unsigned input = data_type_t::INPUT;
        const size_t before = storage_link_transactions_dram[input];
        const size_t target_payload = runtime_datatypes().payload_transactions(
            data_type_t::INPUT, input_halo_unique_elements, dram_input_link_bits);
        const size_t target_metadata = runtime_datatypes().metadata_transactions(
            data_type_t::INPUT, input_halo_unique_elements, dram_input_link_bits);
        const size_t target_storage = runtime_datatypes().storage_transactions(
            data_type_t::INPUT, input_halo_unique_elements, dram_input_link_bits);
        input_halo_pre_dram_transactions = before;
        const bool coalesce_down = target_storage < before;   // retention claim: capacity-gated
        const bool correct_up = target_storage > before;      // union lower bound: unconditional
        if(before > 0 && (correct_up || (coalesce_down && input_halo_capacity_sufficient))) {
            const double ratio = static_cast<double>(target_storage)/static_cast<double>(before);

            payload_link_transactions_multi_chip[input] = rescale_counter_ratio(
                payload_link_transactions_multi_chip[input], target_storage, before,
                "halo-aware multi-chip payload transactions");
            metadata_link_transactions_multi_chip[input] = rescale_counter_ratio(
                metadata_link_transactions_multi_chip[input], target_storage, before,
                "halo-aware multi-chip metadata transactions");
            storage_link_transactions_multi_chip[input] = rescale_counter_ratio(
                storage_link_transactions_multi_chip[input], target_storage, before,
                "halo-aware multi-chip storage transactions");

            payload_link_transactions_dram[input] = target_payload;
            metadata_link_transactions_dram[input] = target_metadata;
            storage_link_transactions_dram[input] = target_storage;

            access_cycle_multi_chip[input] *= ratio;
            access_energy_multi_chip[input] *= ratio;
            transfer_cycle_multi_chip[input] *= ratio;
            transfer_energy_multi_chip[input] *= ratio;
            // Rebuild DRAM source accesses from the exact coalesced union. Unlike multiplying the
            // descriptor aggregate by `ratio`, this preserves an integer physical access count
            // and the access-energy/access-cycle unit-cost identity (energy E5a).
            const size_t target_source_accesses = runtime_datatypes().storage_transactions(
                data_type_t::INPUT, input_halo_unique_elements, dram_input_source_line_bits);
            access_cycle_dram[input] = dram_input_access_hidden
                ? 0.0 : static_cast<double>(target_source_accesses)*dram_input_read_cycle;
            access_energy_dram[input] =
                static_cast<double>(target_source_accesses)*dram_input_read_energy;
            cycle_chip_dram[input] *= ratio;
            transfer_cycle_dram[input] *= ratio;
            transfer_energy_dram[input] *= ratio;
            fill_access_cycle_global_buffer[input] *= ratio;
            fill_access_energy_global_buffer[input] *= ratio;
            for(unsigned chip = 0; chip < chip_fill_access_cycle_global_buffer.size(); ++chip) {
                chip_fill_access_cycle_global_buffer[chip][input] *= ratio;
            }
            input_halo_reuse_applied = true;
        }
    }
    merge_global_buffer_fill();
    finalize_layer_timeline();
}

// CE1/CE2/CE4: the live pass ran one GLB repetition; scaling multiplied the cycle
// vectors per datatype and added fold/setup. Stage busy and the critical path must
// come from those FINAL vectors, not from a uniform blow-up of the pre-scale state.
// CE4 combination rule: a single physical link/port serializes the datatype
// streams (sum across types); separate per-type partitions serve them concurrently
// (max across types).
void stats_t::finalize_layer_timeline() {
    auto sum_types = [](const std::vector<double> &values) {
        double total = 0.0;
        for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) total += values[i];
        return total;
    };
    // DRAM: one device/channel and one off-chip link.
    stage_axis_access[0] = sum_types(access_cycle_dram);
    stage_axis_link[0] = sum_types(transfer_cycle_dram);
    stage_axis_overlap[0] = sum_types(cycle_chip_dram);
    // KV-cache read (evaluation.md Sec 4): the decode step's KV-cache DRAM read is extra
    // DRAM ACCESS traffic on this layer's device -- it serializes on the same channel as the
    // weight/input/output reads, so it adds to the DRAM access axis. Injected here (post
    // repetition scaling) so it is charged exactly once per decode step, and left out of the
    // per-datatype vectors so weight decompression never touches it (KV is incompressible).
    if(kv_pending) {
        stage_axis_access[0] += kv_read_cycles;
    }
    // Multi-chip: one NoP fabric.
    stage_axis_access[1] = sum_types(access_cycle_multi_chip);
    stage_axis_link[1] = sum_types(transfer_cycle_multi_chip);
    stage_axis_overlap[1] = 0.0;
    // GLB: a shared SRAM serializes the three streams; separate partitions do not.
    // CE4: axes 2-4 come from the per-entity combined values (max_entity(sum/max_type(...))),
    // not sum/max_types() over the flat vectors (which would be
    // sum/max_type(max_entity(...)) -- wrong for asymmetric per-entity mappings; see
    // entity_combined_* in stats.h). Axes 3-4 are combined in update_stats(); axis 2 is
    // combined in merge_global_buffer_fill(), because its two sides (GLB read vs
    // multi-chip->GLB fill) scale by different repetition factors and can only be added
    // together once both are final.
    // CE4/P1-D: entity_combined_access_global_buffer already carries base+fill combined
    // per chip and maxed across chips (merge_global_buffer_fill()), so do NOT add the
    // fill-only diagnostic on top -- that would both double-count the fill and re-apply
    // the wrong (independent-max) entity order.
    stage_axis_access[2] = entity_combined_access_global_buffer;
    stage_axis_link[2] = entity_combined_link_global_buffer;
    stage_axis_overlap[2] = entity_combined_overlap_global_buffer;
    // PE array: one temporal buffer and one distribution fabric.
    stage_axis_access[3] = entity_combined_access_pe_array;
    stage_axis_link[3] = entity_combined_link_pe_array;
    stage_axis_overlap[3] = entity_combined_overlap_pe_array;
    // PE: the compute schedule serializes with fold fill and setup (CE2); the local
    // buffer combines per its memory type; MAC<->LB is one per-PE bus.
    stage_axis_compute = computation_cycle + fold_fill_cycle_pe_array;
    stage_axis_access[4] = entity_combined_access_lb;
    stage_axis_link[4] = entity_combined_link_pe;
    stage_axis_overlap[4] = entity_combined_overlap_mac_lb;
    // P4-4: the format-IP axis (payload/metadata/spill packing on LB<->MAC
    // transactions) previously only fed pe_t::modeled_elapsed_cycles()'s inert
    // pre-scale estimate; it now participates in the authoritative busy_cycle_pe
    // combination like the other PE axes, so a slow format-IP is no longer invisible
    // to (or silently fully hidden behind) the critical path.
    stage_axis_format = entity_combined_format;

    busy_cycle_dram = std::max(stage_axis_access[0], std::max(stage_axis_link[0], stage_axis_overlap[0]));
    busy_cycle_multi_chip = std::max(stage_axis_access[1], std::max(stage_axis_link[1], stage_axis_overlap[1]));
    busy_cycle_global_buffer = std::max(stage_axis_access[2], std::max(stage_axis_link[2], stage_axis_overlap[2]));
    busy_cycle_pe_array = std::max(stage_axis_access[3], std::max(stage_axis_link[3], stage_axis_overlap[3]));
    // LB7/P1-A: the three PE transfer axes are three VIEWS of one LB<->MAC transfer,
    // not three additive workloads. access[4] (LB port), link[4] (LB<->MAC bus) and
    // overlap[4] (cycle_mac_lb) are all accumulated from the same
    // account_descriptor_dense_mac_transfer() call, and overlap[4] is already the
    // pipelined MAKESPAN of source access + link + destination access -- so it
    // dominates each individual axis by construction. Summing all three (the pre-fix
    // single-buffer formula) charged every transaction two or three times over.
    //
    // Single buffering does not change what the transfer costs; it removes the ability
    // to hide that transfer BEHIND COMPUTE (there is no second half to fill while the
    // current tile is consumed). So the serialization is
    //   compute + transfer makespan + format-IP,
    // where the transfer makespan is the overlap axis (max() with the per-axis totals
    // only guards a future accounting path that fails to populate it), and the
    // format-IP axis sits on the same non-overlappable tile-load path.
    // A double-buffered LB keeps the original max() over all axes.
    const double lb_transfer_makespan = std::max(stage_axis_overlap[4],
                                                 std::max(stage_axis_access[4], stage_axis_link[4]));
    busy_cycle_pe = pe_double_buffer
        ? std::max(std::max(stage_axis_compute, stage_axis_access[4]),
                   std::max(stage_axis_link[4], std::max(stage_axis_overlap[4], stage_axis_format)))
        : stage_axis_compute + lb_transfer_makespan + stage_axis_format;

    const double stage_busy[5] = {busy_cycle_dram, busy_cycle_multi_chip,
                                  busy_cycle_global_buffer, busy_cycle_pe_array, busy_cycle_pe};
    // L5: the critical path is a PER-TILE event timeline over the five stages, not a
    // closed form over stage totals. Each stage is one serial resource with a per-tile
    // cost; each boundary is a staging buffer of a known depth in tiles; a stage starts a
    // tile when its own previous tile is done, its input tile has arrived, AND a downstream
    // slot is free. That last term is the back-pressure the previous average-per-tile
    // formula could not express -- a fast producer in front of a slow consumer now stalls
    // when the queue fills instead of running to completion in parallel.
    //
    // Buffer depth is the generalization of the old boolean "does this boundary overlap":
    //   depth 1 (single buffer)  -- the producer cannot run ahead, so the two sides
    //                              serialize per tile, which sums to exactly the segment
    //                              sum the old non-overlap branch produced;
    //   depth 2 (double buffer)  -- fill one half while the other drains, which for uniform
    //                              per-tile costs reproduces the old fill+bottleneck value.
    // So both previous cases are special cases; what is NEW is (a) rate-mismatched stages
    // (DRAM/multi-chip refetch fewer tiles than the on-chip stages consume, and used to
    // force the whole run back to a flat max()), (b) per-tile cost variation, and (c) the
    // stall attribution reported below.
    const unsigned tiles = std::max(1U, temporal_repetition_tiles);
    const unsigned offchip_tiles = std::min(tiles, std::max(1U, offchip_repetition_tiles));
    std::vector<std::vector<double>> stage_tile_costs(5, std::vector<double>(tiles, 0.0));
    for(unsigned stage = 0; stage < 5; ++stage) {
        // Stages 0 and 1 are off-chip: they serve only `offchip_tiles` of the shared tile
        // clock, spread evenly, and carry their whole cost on those tiles. Stages 2-4 scale
        // uniformly with the tile clock.
        const bool offchip = stage < 2;
        const unsigned served = offchip ? offchip_tiles : tiles;
        const double per_served_tile = stage_busy[stage]/static_cast<double>(served);
        for(unsigned k = 0; k < tiles; ++k) {
            const bool serves = !offchip || k == 0 ||
                (static_cast<unsigned long long>(k)*served)/tiles !=
                (static_cast<unsigned long long>(k - 1)*served)/tiles;
            if(serves) stage_tile_costs[stage][k] = per_served_tile;
        }
    }
    std::vector<unsigned> boundary_depths(4, 1);
    for(unsigned b = 0; b < 4; ++b) {
        boundary_depths[b] = timeline_boundary_overlap[b] ? 2 : 1;
    }
    // L6: the multi-chip -> GLB boundary is only decoupled for the datatypes actually
    // STAGED in the GLB. A GLB-bypassed datatype has no buffer there at all -- it is
    // streamed from the chip ingress straight through to the PE array -- so that stream
    // cannot hold a tile back while the next one arrives. One bypassed datatype therefore
    // removes the boundary's tile-level decoupling, no matter what the GLB's own
    // double_buffer flag says. (Splitting the boundary per datatype would need the GLB
    // stage itself decomposed per datatype; the depth is the boundary-level contract.)
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        if(global_buffer_bypass[i]) boundary_depths[1] = 1;
    }
    // MC2: an explicit outstanding-request limit overrides the DRAM <-> multi-chip depth only.
    // The boolean double_buffer contract cannot express more than 2 tiles in flight, and a real
    // memory interface (NVDLA's DBBIF, and every AXI-class port) is specified by exactly this
    // number. Unset (0) leaves the derived depth alone, so declaring nothing changes nothing --
    // and an explicit 1 must reproduce the implicit 1 bit for bit, which is what MC1 checks.
    if(timeline_offchip_outstanding > 0) boundary_depths[0] = timeline_offchip_outstanding;
    for(unsigned b = 0; b < 4; ++b) timeline_boundary_depth[b] = boundary_depths[b];
    std::vector<double> stall;
    double final_latency = pipeline_timeline_cycles(stage_tile_costs, boundary_depths, &stall);
    for(unsigned stage = 0; stage < 5; ++stage) {
        timeline_stall[stage] = (stage < stall.size()) ? stall[stage] : 0.0;
    }

    // Phase-5 (plan_sfu.md): integrate the pending SFU activation window into the SAME
    // timeline the producer ran on.
    //
    //   queue_depth = 1  (serial contract) -- the SFU drains the final output after the
    //       producer pipeline completes; its whole window extends the critical path.
    //       This reproduces the pre-streaming model bit for bit.
    //   queue_depth >= 2 (streaming) -- the SFU is a SIXTH per-tile pipeline stage behind
    //       a bounded output-tile queue. Final outputs are released only on the tiles
    //       that COMMIT outputs: a repetition group over a reduction dimension commits
    //       its output tile on the group's LAST tile (the activation cannot run before
    //       the reduction finishes), so the SFU's cost lands there rather than being
    //       spread optimistically over the whole group. A fast SFU then hides behind the
    //       producer except for its drain tail; a slow one back-pressures the producer
    //       through the queue, and that stall is attributed to the PE stage.
    if(sfu_pending && sfu_pending_invocation.busy_cycle > 0.0) {
        const double producer_latency = final_latency;
        if(sfu_queue_depth >= 2) {
            const unsigned sfu_tiles = std::min(tiles, std::max(1u, output_repetition_tiles));
            std::vector<double> sfu_tile_costs(tiles, 0.0);
            const double per_output_tile =
                sfu_pending_invocation.busy_cycle/static_cast<double>(sfu_tiles);
            // Serve on the LAST tile of each output-commit group (mirror of the off-chip
            // stages' first-tile rule -- a commit CLOSES its reduction group): group g of
            // the sfu_tiles groups closes at tile floor((g+1)*tiles/sfu_tiles) - 1.
            for(unsigned g = 0; g < sfu_tiles; ++g) {
                const unsigned close_tile = static_cast<unsigned>(
                    (static_cast<unsigned long long>(g + 1)*tiles)/sfu_tiles) - 1u;
                sfu_tile_costs[close_tile] += per_output_tile;
            }
            std::vector<std::vector<double>> streaming_costs = stage_tile_costs;
            streaming_costs.push_back(sfu_tile_costs);
            std::vector<unsigned> streaming_depths = boundary_depths;
            streaming_depths.push_back(sfu_queue_depth);
            std::vector<double> streaming_stall;
            const double streaming_latency = pipeline_timeline_cycles(
                streaming_costs, streaming_depths, &streaming_stall);
            for(unsigned stage = 0; stage < 5; ++stage) {
                timeline_stall[stage] = (stage < streaming_stall.size())
                                        ? streaming_stall[stage] : 0.0;
            }
            final_latency = streaming_latency;
            sfu_serial_cycle = std::max(0.0, streaming_latency - producer_latency);
            sfu_hidden_cycle = std::max(0.0,
                sfu_pending_invocation.busy_cycle - sfu_serial_cycle);
            // The PE stage's stall in the six-stage run is time it sat blocked on a full
            // SFU queue -- the producer back-pressure the plan asks to attribute.
            sfu_stall_cycle = (4 < streaming_stall.size()) ? streaming_stall[4] : 0.0;
        } else {
            final_latency += sfu_pending_invocation.busy_cycle;
            sfu_serial_cycle = sfu_pending_invocation.busy_cycle;
            sfu_hidden_cycle = 0.0;
            // Serial contract: the output drain waits for the whole SFU window.
            sfu_stall_cycle = sfu_pending_invocation.busy_cycle;
        }
        sfu_busy_cycle = sfu_pending_invocation.busy_cycle;
    }

    // Weight decompression (evaluation.md Sec 4): the decoder is a stage on the WEIGHT
    // supply path (DRAM compressed load -> decoder -> scratchpad -> compute). The weight
    // DRAM saving is already in the reduced DRAM busy above; here the decoder's own work
    // is placed on the timeline. overlap OFF serializes the whole decoder window onto the
    // critical path; overlap ON pipelines it against the producer/compute so a fast
    // decoder hides (fill/drain only) and a slow one back-pressures and is exposed. The
    // 3-stage pipeline (DRAM-supply -> decoder -> compute) is modeled with
    // pipeline_timeline_cycles over per-tile costs.
    //
    // Relative-throughput mode (decoder_ratio): weight demand = dense weight / the on-chip
    // COMPUTE-SIDE latency -- the time the datapath takes to process the layer if memory
    // were free = max on-chip busy axis (PE array / PE / GLB / multi-chip), DRAM excluded.
    // The pure MAC schedule (stage_axis_compute) is wrong here: a fill-bound array (e.g.
    // decode M=1 reloads the whole array to produce one output row) consumes weight over
    // its entire busy window, not just its MAC cycles. With throughput = ratio x demand the
    // decoder window collapses to startup + T_compute / ratio, so ratio 1 just keeps pace
    // with the datapath, < 1 exposes the decoder, > 1 leaves slack. The busy axes were
    // recomputed just above, which is why this derivation is deferred to here.
    if(decomp_pending && decomp_decoder_ratio > 0.0) {
        const double compute_side =
            std::max(std::max(busy_cycle_pe_array, busy_cycle_pe),
                     std::max(busy_cycle_global_buffer, busy_cycle_multi_chip));
        if(compute_side > 0.0) {
            decomp_decoder_cycles = decomp_startup_cycles + compute_side/decomp_decoder_ratio;
        }
    }
    if(decomp_pending && decomp_decoder_cycles > 0.0) {
        const double decoder = decomp_decoder_cycles;
        const double producer_latency = final_latency;
        if(decomp_overlap && decomp_queue_depth >= 1) {
            // The decoder pipelines over the WEIGHT tiles it processes (decomp_tiles), not
            // the layer's temporal-reuse tiles: the compressed weight arrives from DRAM in
            // decomp_tiles pieces regardless of temporal reuse, and the decoder overlaps
            // with the producer across them. (Using `tiles` collapsed decode M=1 -- one
            // temporal tile -- to a serial decoder that overlap could never hide.) Capped so
            // the per-tile vectors stay small; the pipeline makespan is granularity-stable.
            const unsigned dtiles = std::max<unsigned>(1,
                std::min<unsigned>(2048u, static_cast<unsigned>(decomp_tiles)));
            // Producer = the memory/compute pipeline delivering compressed weight tiles
            // (the layer's current per-tile makespan); decoder = a following stage of
            // per-tile cost; consumer = a sink. A bounded queue of depth queue_depth sits
            // between producer and decoder, and the decoder's dense-output buffer bounds
            // the decoder -> scratchpad boundary.
            std::vector<double> supply(dtiles, producer_latency/static_cast<double>(dtiles));
            if(!decomp_tile_supply_fraction.empty()) {
                // Per-tile compression variation (evaluation discussion 2026-09-03 Sec 5):
                // the compressed stream arrives BURSTY -- each supply tile takes time in
                // proportion to the compressed bytes it actually carries, so well-compressed
                // tiles arrive fast (and pile into the bounded queue) while poorly
                // compressed ones starve the decoder. Same total supply time, same DRAM
                // bytes as the uniform model -- only the arrival pattern differs, which is
                // exactly what a layer-mean analytical model cannot represent.
                std::vector<double> binned(dtiles, 0.0);
                const size_t chunks = decomp_tile_supply_fraction.size();
                for(size_t j = 0; j < chunks; ++j) {
                    binned[j*dtiles/chunks] += decomp_tile_supply_fraction[j];
                }
                for(unsigned k = 0; k < dtiles; ++k) {
                    supply[k] = producer_latency*binned[k];
                }
            }
            // Decode cost per tile stays uniform: the decoder emits DENSE bytes at a fixed
            // throughput, and dense tile sizes are uniform regardless of compression.
            std::vector<double> decode(dtiles, decoder/static_cast<double>(dtiles));
            // Sink = the scratchpad absorbing the dense output at the GLB weight-write
            // rate. With a real (non-zero) sink cost the decoder -> scratchpad boundary
            // depth (output_buffer_tiles) participates in the back-pressure.
            std::vector<double> sink(dtiles, decomp_sink_cycles/static_cast<double>(dtiles));
            std::vector<unsigned> depths = {decomp_queue_depth,
                                            std::max(1u, decomp_output_buffer_tiles)};
            std::vector<double> dstall;
            const double piped = pipeline_timeline_cycles({supply, decode, sink}, depths,
                                                          &dstall);
            final_latency = piped;
            decomp_exposed_cycle = std::max(0.0, piped - producer_latency);
            decomp_hidden_cycle = std::max(0.0, decoder - decomp_exposed_cycle);
            decomp_stall_cycle = (0 < dstall.size()) ? dstall[0] : 0.0;
        } else {
            final_latency += decoder;
            decomp_exposed_cycle = decoder;
            decomp_hidden_cycle = 0.0;
            decomp_stall_cycle = decoder;   // serial: the supply waits on the decoder
        }
    }

    // KV supply scheduling (evaluation discussion 2026-09-03 Sec 7): with a non-aggregate
    // [kvcache] kv_schedule, the KV read is a tile-level supply stage integrated against
    // the layer window instead of extra DRAM-axis busy. The layer window stands in for the
    // attention compute that consumes the KV tiles (attention itself is outside the modeled
    // scope -- stated in the report, not hidden). Same KV traffic and energy in every mode;
    // only the schedule differs:
    //   blocking        -- the whole KV read completes before compute: fully exposed.
    //   streaming       -- KV tiles pipeline against the window behind a bounded buffer
    //                      (kv_buffer_tiles): the fill and any supply excess are exposed.
    //   double_buffered -- streaming with the fill prefetched away: only the supply excess
    //                      over the window is exposed (ideal prefetch bound).
    if(kv_sched_pending && kv_read_cycles > 0.0) {
        const double producer_latency = final_latency;
        if(kv_schedule == "blocking") {
            final_latency += kv_read_cycles;
            kv_exposed_cycle = kv_read_cycles;
            kv_hidden_cycle = 0.0;
            kv_stall_cycle = kv_read_cycles;   // compute waits out the whole KV read
        } else if(kv_schedule == "double_buffered") {
            final_latency = std::max(producer_latency, kv_read_cycles);
            kv_exposed_cycle = std::max(0.0, kv_read_cycles - producer_latency);
            kv_hidden_cycle = kv_read_cycles - kv_exposed_cycle;
            kv_stall_cycle = 0.0;
        } else {   // streaming
            const unsigned nt = std::max<unsigned>(1,
                std::min<unsigned>(2048u, static_cast<unsigned>(kv_tiles)));
            std::vector<double> supply(nt, kv_read_cycles/static_cast<double>(nt));
            std::vector<double> consume(nt, producer_latency/static_cast<double>(nt));
            std::vector<unsigned> depths = {std::max(1u, kv_buffer_tiles)};
            std::vector<double> kstall;
            const double piped = pipeline_timeline_cycles({supply, consume}, depths, &kstall);
            final_latency = piped;
            kv_exposed_cycle = std::max(0.0, piped - producer_latency);
            kv_hidden_cycle = std::max(0.0, kv_read_cycles - kv_exposed_cycle);
            kv_stall_cycle = (0 < kstall.size()) ? kstall[0] : 0.0;
        }
        kv_sched_pending = false;
    }

    // Attention consumer step (evaluation discussion 2026-09-03 Sec 7): one decode step's
    // QK/softmax/AV appended AFTER the mapped layer (Q depends on the projections, so the
    // step cannot start earlier). Inside the step the KV supply carries a real tile
    // dependency -- tile k's QK/AV cannot start before tile k arrives through the bounded
    // KV buffer -- and the algorithm decides the softmax structure: "online" interleaves a
    // single K/V pass with a running softmax (no barrier); "two_pass" scores every cached
    // token first (K pass), runs the softmax BARRIER over all scores, then streams V for
    // AV. double_buffered is the ideal-prefetch bound of each pass (the fill hides behind
    // the mapped layer's window). The step ends with the cache append write.
    if(attn_pending) {
        const double supply_k = attn_supply_k_cycle;
        const double supply_v = attn_supply_v_cycle;
        const double compute_total = attn_qk_cycles + attn_softmax_cycles + attn_av_cycles;
        double step = 0.0;
        attn_stall_cycle = 0.0;
        if(kv_schedule == "blocking") {
            step = supply_k + supply_v + compute_total;
        } else if(kv_schedule == "double_buffered") {
            step = (attn_algorithm == "two_pass")
                ? std::max(supply_k, attn_qk_cycles) + attn_softmax_cycles +
                  std::max(supply_v, attn_av_cycles)
                : std::max(supply_k + supply_v, compute_total);
        } else {   // streaming
            const unsigned nt = std::max<unsigned>(1,
                std::min<unsigned>(2048u, static_cast<unsigned>(kv_tiles)));
            std::vector<unsigned> depths = {std::max(1u, kv_buffer_tiles)};
            auto piped = [&](double m_supply, double m_consume) {
                std::vector<double> supply(nt, m_supply/static_cast<double>(nt));
                std::vector<double> consume(nt, m_consume/static_cast<double>(nt));
                std::vector<double> stall_v;
                const double out = pipeline_timeline_cycles({supply, consume}, depths,
                                                            &stall_v);
                attn_stall_cycle += (0 < stall_v.size()) ? stall_v[0] : 0.0;
                return out;
            };
            step = (attn_algorithm == "two_pass")
                ? piped(supply_k, attn_qk_cycles) + attn_softmax_cycles +
                  piped(supply_v, attn_av_cycles)
                : piped(supply_k + supply_v, compute_total);
        }
        step += attn_write_cycle;
        attn_step_latency = step;
        attn_exposed_cycle = std::max(0.0, step - attn_write_cycle - compute_total);
        attn_hidden_cycle = std::max(0.0, supply_k + supply_v - attn_exposed_cycle);
        final_latency += step;
        attn_pending = false;
    }

    // Leakage accrued linearly over the uniform-scaled pre-recompute latency;
    // rescale it to the final window.
    if(layer_latency > 0.0 && final_latency != layer_latency) {
        const double leakage_scale = final_latency/layer_latency;
        for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
            static_energy_pe[i] *= leakage_scale;
            static_energy_pe_array[i] *= leakage_scale;
            static_energy_global_buffer[i] *= leakage_scale;
            static_energy_multi_chip[i] *= leakage_scale;
        }
    }
    layer_latency = final_latency;
    mac_available_cycle = timeline_physical_macs*layer_latency;
    utilization_mac = calculate_time_based_mac_utilization(mac_busy_cycle, mac_available_cycle);
    utilization_mac = std::min(1.0, utilization_mac);
    // Plan energy model: every physical SFU is always-on over the layer's final window
    // (bypass layers included -- no power gating is modeled).
    if(sfu_pending) {
        sfu_static_energy = static_cast<double>(sfu_physical_units)*layer_latency*
                            sfu_pending_static_energy_per_cycle;
        sfu_pending = false;
    }
    // Decompression: cycle breakdown (evaluation.md 4.C) and always-on decoder leakage.
    if(decomp_present) {
        decomp_static_energy = layer_latency*decomp_pending_static_energy_per_cycle;
        // The four buckets: memory (DRAM incl. compressed weight), decoder exposure,
        // compute schedule, and pipeline stall. They partition the critical path.
        decomp_breakdown_dram = busy_cycle_dram;
        decomp_breakdown_decode = decomp_exposed_cycle;
        decomp_breakdown_compute = stage_axis_compute;
        decomp_breakdown_stall = std::max(0.0, layer_latency - std::max(busy_cycle_dram,
            std::max(stage_axis_compute, decomp_exposed_cycle)));
        decomp_pending = false;
    }
}

void stats_t::merge_global_buffer_fill() {
    // CE4/P1-D: the fill (mc->GLB write) side scales per-datatype rather than uniformly
    // (see scale_serial_repetitions), so both sides reach their final value only here --
    // which is also the only point where the GLB access axis can be combined in the
    // correct order. Per chip, first combine that chip's OWN datatypes (a shared GLB
    // serializes its three streams on one port, separate partitions do not), then take
    // the max across chip entities:
    //     max_chip(combine_type(base[chip][t]) + combine_type(fill[chip][t]))
    // NOT sum_type(max_chip(...)), which lets an input-bound chip 0 and a weight-bound
    // chip 1 report 200 cycles of GLB time for work the two chips run in parallel in 100.
    const bool serialized_types = global_buffer_type == memory_type_t::SHARED;
    // Diagnostic fill-only component (reported) and the authoritative stage axis, both
    // reduced by the shared CE4 helper so the ordering is covered by hand-computed unit
    // tests (unittest/validation_test.cc) rather than only by inspection here. L1: the
    // axis adds the two sides PER DATATYPE before combining datatypes -- on a separate
    // buffer, max_type(base) and max_type(fill) may be peaks of different partitions,
    // which run in parallel and must not be summed.
    stage_fill_access_global_buffer = entity_combined_cycles(chip_fill_access_cycle_global_buffer,
                                                            serialized_types);
    entity_combined_access_global_buffer = entity_combined_cycles(chip_access_cycle_global_buffer,
                                                                 chip_fill_access_cycle_global_buffer,
                                                                 serialized_types);
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        access_cycle_global_buffer[i] += fill_access_cycle_global_buffer[i];
        access_energy_global_buffer[i] += fill_access_energy_global_buffer[i];
        fill_access_cycle_global_buffer[i] = 0.0;
        fill_access_energy_global_buffer[i] = 0.0;
    }
}

void stats_t::set_sfu_activation(const sfu_invocation_t &m_invocation,
                                 unsigned m_physical_units, unsigned m_lanes,
                                 double m_static_energy_per_cycle, unsigned m_queue_depth) {
    sfu_present = true;
    // "Active" means an SFU event actually FIRED. A linear-bypass layer records the
    // operation name and accrues always-on leakage (no power gating is modeled), but it
    // must not claim a dynamic event happened.
    sfu_active = sfu_active || m_invocation.operations > 0 || m_invocation.busy_cycle > 0.0;
    append_unique_segment(sfu_operation, m_invocation.operation, ", ");
    sfu_physical_units = m_physical_units;
    sfu_lanes = m_lanes;
    sfu_queue_depth = m_queue_depth;
    sfu_valid_elements += m_invocation.valid_elements;
    sfu_operations += m_invocation.operations;
    sfu_invocations += m_invocation.invocations;
    sfu_chunks += m_invocation.chunks;
    sfu_tail_lane_utilization = std::max(sfu_tail_lane_utilization,
                                         m_invocation.tail_lane_utilization);
    sfu_ingress_elements += m_invocation.ingress_elements;
    sfu_egress_elements += m_invocation.egress_elements;
    sfu_ingress_bytes += m_invocation.ingress_bytes;
    sfu_egress_bytes += m_invocation.egress_bytes;
    sfu_ingress_transactions += m_invocation.ingress_transactions;
    sfu_egress_transactions += m_invocation.egress_transactions;
    sfu_op_energy += m_invocation.op_energy;
    sfu_read_energy += m_invocation.read_energy;
    sfu_write_energy += m_invocation.write_energy;
    sfu_setup_energy += m_invocation.setup_energy;
    sfu_timing_calibrated = sfu_timing_calibrated && m_invocation.timing_calibrated;
    for(unsigned i = 0; i < m_invocation.unpriced.size(); ++i) {
        if(std::find(sfu_unpriced_events.begin(), sfu_unpriced_events.end(),
                     m_invocation.unpriced[i]) == sfu_unpriced_events.end()) {
            sfu_unpriced_events.push_back(m_invocation.unpriced[i]);
        }
    }
    // Timing integration is owned by finalize_layer_timeline(): serial append under the
    // queue_depth = 1 contract, or the sixth pipeline stage under streaming.
    sfu_pending = true;
    sfu_pending_invocation = m_invocation;
    sfu_pending_static_energy_per_cycle = m_static_energy_per_cycle;
}

void stats_t::apply_sfu_activation(const sfu_invocation_t &m_invocation,
                                   unsigned m_physical_units, unsigned m_lanes,
                                   double m_static_energy_per_cycle) {
    // SFU-only path (standalone softmax): there is no producer timeline to stream
    // against -- the multi-pass window IS the layer's execution, serially.
    set_sfu_activation(m_invocation, m_physical_units, m_lanes, m_static_energy_per_cycle, 1);
    sfu_pending = false;   // consumed here; finalize_layer_timeline() never runs

    const double previous_latency = layer_latency;
    layer_latency += m_invocation.busy_cycle;
    sfu_busy_cycle += m_invocation.busy_cycle;
    sfu_serial_cycle += m_invocation.busy_cycle;
    sfu_stall_cycle += m_invocation.busy_cycle;
    if(previous_latency > 0.0 && layer_latency != previous_latency) {
        // Leakage accrues over wall-clock: rescale every component's leakage window to
        // the SFU-extended latency, exactly as finalize_layer_timeline() rescales it.
        const double leakage_scale = layer_latency/previous_latency;
        for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
            static_energy_pe[i] *= leakage_scale;
            static_energy_pe_array[i] *= leakage_scale;
            static_energy_global_buffer[i] *= leakage_scale;
            static_energy_multi_chip[i] *= leakage_scale;
        }
    }
    // MAC availability follows the final latency: a slower SFU lowers MAC utilization.
    mac_available_cycle = timeline_physical_macs*layer_latency;
    if(mac_available_cycle > 0.0) {
        utilization_mac = std::min(1.0,
            calculate_time_based_mac_utilization(mac_busy_cycle, mac_available_cycle));
    }
    // Plan energy model: every physical SFU in the config is always-on for the layer's
    // final wall-clock window (no power gating modeled yet).
    sfu_static_energy = static_cast<double>(m_physical_units)*layer_latency*
                        m_static_energy_per_cycle;
}

void stats_t::record_sfu_only_layer(const sfu_invocation_t &m_invocation,
                                    unsigned m_physical_units, unsigned m_lanes,
                                    double m_static_energy_per_cycle,
                                    double m_frequency_mhz, bool m_single_clock,
                                    const std::string &m_clock_note,
                                    const sfu_operand_stream_t &m_stream) {
    // The SFU is the only modeled resource of this layer: the contract fields that
    // update_stats() normally captures from the memory/compute components are stated
    // here so the layer report and the power gate stay well defined.
    authoritative_frequency_mhz = m_frequency_mhz;
    single_clock_domain = m_single_clock;
    clock_domain_note = m_clock_note;
    compute_energy_basis = "n/a (SFU-only layer; no MAC activity)";
    // Vacuously calibrated: no MAC event fires in this layer, so the MAC precision
    // qualification cannot be the reason its energy is not absolute.
    compute_energy_precision_calibrated = true;
    operand_precision = runtime_datatypes().accumulator_format().name;
    dram_timing_model = m_stream.active
        ? "analytical unit-cost streaming of standalone graph operands (no bank-conflict"
          " scheduling for this layer)"
        : "n/a (no DRAM activity is modeled for this SFU-only layer)";
    dram_timing_limits = "n/a";
    array_noc_link_contract = "n/a (SFU-only layer)";
    nop_link_contract = "n/a (SFU-only layer)";
    // Tile sizes exist only for mapped layers; zero-fill so the stat printer stays safe.
    for(unsigned i = 0; i < component_type_t::NUM_COMPONENT_TYPES; ++i) {
        tile_size[i].assign(data_type_t::NUM_DATA_TYPES, 0);
    }
    sfu_only_layer = true;
    apply_sfu_activation(m_invocation, m_physical_units, m_lanes, m_static_energy_per_cycle);

    // Phase-7: operand streaming between the memory hierarchy and the SFU. The transfer
    // costs land on the components that own them (DRAM device+link, GLB ports), so the
    // energy summary's DRAM/GLB rows carry them; the serial ingress/egress makespans
    // extend the layer's critical path around the SFU window
    // (ingress -> multi-pass SFU -> egress -- the conservative serial contract).
    if(m_stream.active) {
        sfu_stream_active = true;
        sfu_stream_residency = m_stream.residency;
        sfu_stream_ingress_bytes += m_stream.ingress_bytes;
        sfu_stream_egress_bytes += m_stream.egress_bytes;
        sfu_stream_ingress_cycle += m_stream.ingress_cycle;
        sfu_stream_egress_cycle += m_stream.egress_cycle;
        access_cycle_dram[data_type_t::OUTPUT] += m_stream.dram_access_cycle;
        access_energy_dram[data_type_t::OUTPUT] += m_stream.dram_access_energy;
        // Row activations land on the DRAM transfer axis, exactly where dram_t puts its
        // own (account_row_activations); the event count feeds the unpriced-active check
        // for row_miss_energy like any other layer's activations.
        transfer_cycle_dram[data_type_t::OUTPUT] += m_stream.dram_link_cycle +
                                                    m_stream.dram_row_activation_cycle;
        transfer_energy_dram[data_type_t::OUTPUT] += m_stream.dram_link_energy +
                                                     m_stream.dram_row_activation_energy;
        row_activation_events += m_stream.dram_row_activations;
        storage_link_transactions_dram[data_type_t::OUTPUT] += m_stream.dram_link_transactions;
        payload_link_transactions_dram[data_type_t::OUTPUT] += m_stream.dram_link_transactions;
        access_cycle_global_buffer[data_type_t::OUTPUT] += m_stream.glb_access_cycle;
        access_energy_global_buffer[data_type_t::OUTPUT] += m_stream.glb_access_energy;
        // Busy axes for the layer-timeline table (busy = max of a stage's axes).
        stage_axis_access[0] += m_stream.dram_access_cycle;
        stage_axis_link[0] += m_stream.dram_link_cycle + m_stream.dram_row_activation_cycle;
        busy_cycle_dram = std::max(stage_axis_access[0], stage_axis_link[0]);
        stage_axis_access[2] += m_stream.glb_access_cycle;
        entity_combined_access_global_buffer += m_stream.glb_access_cycle;
        busy_cycle_global_buffer = std::max(busy_cycle_global_buffer,
                                            entity_combined_access_global_buffer);
        // Serial critical path: operand in, multi-pass, operand out. Leakage windows
        // (the SFU's own) follow the extended latency.
        layer_latency += m_stream.ingress_cycle + m_stream.egress_cycle;
        sfu_static_energy = static_cast<double>(m_physical_units)*layer_latency*
                            m_static_energy_per_cycle;
    }
}

void stats_t::apply_graph_residency(bool m_input_in_glb, bool m_output_in_glb,
                                    const std::string &m_note) {
    graph_residency_applied = m_input_in_glb || m_output_in_glb;
    if(!m_note.empty()) graph_residency_note = m_note;
    const unsigned input = data_type_t::INPUT;
    const unsigned output = data_type_t::OUTPUT;
    if(m_input_in_glb) {
        num_request_multi_chip[input] = num_data_transfer_multi_chip[input] = 0;
        payload_link_transactions_multi_chip[input] = 0;
        metadata_link_transactions_multi_chip[input] = 0;
        storage_link_transactions_multi_chip[input] = 0;
        access_cycle_multi_chip[input] = access_energy_multi_chip[input] = 0.0;
        transfer_cycle_multi_chip[input] = transfer_energy_multi_chip[input] = 0.0;
        num_request_dram[input] = num_data_transfer_dram[input] = 0;
        payload_link_transactions_dram[input] = 0;
        metadata_link_transactions_dram[input] = 0;
        storage_link_transactions_dram[input] = 0;
        access_cycle_dram[input] = access_energy_dram[input] = 0.0;
        transfer_cycle_dram[input] = transfer_energy_dram[input] = 0.0;
        cycle_chip_dram[input] = 0.0;
        fill_access_cycle_global_buffer[input] = 0.0;
        fill_access_energy_global_buffer[input] = 0.0;
        for(std::vector<double> &chip : chip_fill_access_cycle_global_buffer) {
            if(input < chip.size()) chip[input] = 0.0;
        }
    }
    if(m_output_in_glb) {
        num_request_multi_chip[output] = num_data_transfer_multi_chip[output] = 0;
        payload_link_transactions_multi_chip[output] = 0;
        metadata_link_transactions_multi_chip[output] = 0;
        storage_link_transactions_multi_chip[output] = 0;
        access_cycle_multi_chip[output] = access_energy_multi_chip[output] = 0.0;
        transfer_cycle_multi_chip[output] = transfer_energy_multi_chip[output] = 0.0;
        num_request_dram[output] = num_data_transfer_dram[output] = 0;
        payload_link_transactions_dram[output] = 0;
        metadata_link_transactions_dram[output] = 0;
        storage_link_transactions_dram[output] = 0;
        access_cycle_dram[output] = access_energy_dram[output] = 0.0;
        transfer_cycle_dram[output] = transfer_energy_dram[output] = 0.0;
        cycle_chip_dram[output] = 0.0;
    }
}

void stats_t::mark_unmodeled_activation(const std::string &m_note) {
    activation_unmodeled = true;
    if(activation_unmodeled_note.empty()) {
        activation_unmodeled_note = m_note;
    }
}

void stats_t::apply_decompression(decomp_t *m_engine, size_t m_dense_weight_bytes,
                                  double m_dram_bytes_per_cycle, size_t m_tile_bytes,
                                  double m_sink_bytes_per_cycle) {
    // Called BEFORE scale_serial_repetitions(): the compression ratio is
    // repetition-independent, so scaling the compressed weight traffic here (pre-scale)
    // and letting scale_serial_repetitions() multiply it afterwards is exact. The
    // compute schedule is only for the reported decoder ratio.
    const double compute_cycles = stage_axis_compute;
    const decomp_invocation_t inv = m_engine->decompress(m_dense_weight_bytes,
                                                         m_dram_bytes_per_cycle,
                                                         compute_cycles, m_tile_bytes);
    decomp_present = true;
    // A decoder runs whenever the layer is not bypassed and a throughput is declared --
    // either an absolute one (inv.decoder_cycles set now) or a relative one (decoder_ratio,
    // whose cycles finalize_layer_timeline() fills in from the final compute schedule).
    decomp_active = inv.active && !inv.bypassed &&
                    (inv.decoder_cycles > 0.0 || m_engine->get_decoder_ratio() > 0.0);
    decomp_bypassed = inv.bypassed;
    decomp_dense_weight_bytes = inv.dense_weight_bytes;
    decomp_compressed_weight_bytes = inv.compressed_weight_bytes;
    decomp_effective_ratio = inv.effective_ratio;
    decomp_tiles = inv.tiles;
    decomp_bypassed_tiles = inv.bypassed_tiles;
    decomp_tile_supply_fraction = inv.tile_supply_fraction;
    decomp_tile_ratio_cv = m_engine->get_tile_ratio_cv();
    decomp_output_buffer_tiles = std::max(1u, m_engine->get_output_buffer_tiles());
    // The scratchpad absorbs the decoder's DENSE output at the GLB weight-write rate; that
    // pacing is the sink stage of the supply->decode->scratchpad pipeline (a zero-cost
    // sink left output_buffer_tiles without any observable effect).
    decomp_sink_cycles = (m_sink_bytes_per_cycle > 0.0 && !inv.bypassed)
        ? static_cast<double>(inv.dense_weight_bytes)/m_sink_bytes_per_cycle : 0.0;
    decomp_decoder_cycles = inv.decoder_cycles;
    decomp_queue_depth = m_engine->get_queue_depth();
    decomp_overlap = m_engine->get_overlap();
    decomp_decoder_energy = inv.decoder_energy;
    decomp_timing_calibrated = inv.timing_calibrated;
    decomp_profile_reference = m_engine->get_profile_reference();
    for(unsigned i = 0; i < inv.unpriced.size(); ++i) {
        decomp_unpriced_events.push_back(inv.unpriced[i]);
    }

    // Effect 1: scale the WEIGHT DRAM traffic to the compressed footprint. Only the
    // compressed bytes are fetched; the saving is what relieves a memory-bound layer.
    // Scale cycles, energy AND the transaction/transfer counters together so the reported
    // DRAM identity (weight transactions x line size == compressed bytes) still holds --
    // otherwise the report shows halved DRAM cycles beside an un-reduced dense transaction
    // count. The counters use the exact integer ratio compressed/dense to stay on a whole
    // number of DRAM lines; the cycle/energy fields use the equivalent double ratio.
    if(inv.dense_weight_bytes > 0 && !inv.bypassed) {
        const size_t dense_b = inv.dense_weight_bytes;
        const size_t comp_b  = inv.compressed_weight_bytes;
        const double ratio = static_cast<double>(comp_b)/static_cast<double>(dense_b);
        decomp_dram_weight_saved_cycle =
            (access_cycle_dram[data_type_t::WEIGHT] +
             transfer_cycle_dram[data_type_t::WEIGHT])*(1.0 - ratio);
        access_cycle_dram[data_type_t::WEIGHT] *= ratio;
        access_energy_dram[data_type_t::WEIGHT] *= ratio;
        transfer_cycle_dram[data_type_t::WEIGHT] *= ratio;
        transfer_energy_dram[data_type_t::WEIGHT] *= ratio;
        cycle_chip_dram[data_type_t::WEIGHT] *= ratio;
        // Volume counters: compressed weight occupies fewer DRAM lines, so fewer fills and
        // fewer off-chip link transactions carry it. Integer-exact so the byte identity holds.
        const data_type_t W = data_type_t::WEIGHT;
        num_data_transfer_dram[W]        = static_cast<unsigned>(
            static_cast<size_t>(num_data_transfer_dram[W]) * comp_b / dense_b);
        payload_link_transactions_dram[W]  = payload_link_transactions_dram[W]  * comp_b / dense_b;
        metadata_link_transactions_dram[W] = metadata_link_transactions_dram[W] * comp_b / dense_b;
        storage_link_transactions_dram[W]  = storage_link_transactions_dram[W]  * comp_b / dense_b;
    }

    // Effect 2: stage the decoder for finalize_layer_timeline() (values already copied
    // into the scalar decomp_* fields above). In relative-throughput mode the decoder
    // cycles are computed in finalize from the final compute schedule, so carry the ratio
    // and startup across.
    decomp_pending = inv.active && !inv.bypassed;
    decomp_pending_static_energy_per_cycle = m_engine->get_static_energy_per_cycle();
    decomp_decoder_ratio = m_engine->get_decoder_ratio();
    decomp_startup_cycles = m_engine->get_startup_cycles();
}

void stats_t::apply_kv_cache_read(kvcache_t *m_engine, size_t m_dense_weight_bytes) {
    // Called BEFORE apply_decompression() (so the weight-DRAM reference is still dense) and
    // BEFORE scale_serial_repetitions(); the KV read is a fixed per-decode-step cost that
    // must NOT be repetition-scaled, so it is NOT written into the DRAM traffic vectors here
    // -- it is staged and injected once into the DRAM axis by finalize_layer_timeline().
    kv_present = true;
    kv_bytes_per_token = m_engine->get_bytes_per_token();
    kv_context_length = m_engine->get_context_length();
    kv_dense_read_bytes = m_engine->dense_read_bytes();
    kv_read_bytes = m_engine->compressed_read_bytes();      // what DRAM actually fetches
    kv_bypassed = m_engine->bypassed();
    kv_compression_ratio = kv_read_bytes > 0
        ? static_cast<double>(kv_dense_read_bytes)/static_cast<double>(kv_read_bytes) : 1.0;
    kv_profile_reference = m_engine->get_profile_reference();

    // Per-byte DRAM cost from THIS layer's own measured (dense) weight-DRAM access: it
    // already carries the accelerator's real DRAM model -- row activations, transfer,
    // read_cycle -- and tracks the bandwidth knob.
    double cost_per_byte = 0.0;
    if(m_dense_weight_bytes > 0 && access_cycle_dram[data_type_t::WEIGHT] > 0.0) {
        cost_per_byte = access_cycle_dram[data_type_t::WEIGHT]/
                        static_cast<double>(m_dense_weight_bytes);
    }
    // Only the COMPRESSED bytes are fetched from DRAM; the decoder then reconstitutes dense
    // KV. Both are on the KV supply path this decode step, so both fold into the DRAM axis.
    kv_dram_read_cycles = static_cast<double>(kv_read_bytes)*cost_per_byte;
    kv_decoder_cycles = m_engine->decoder_cycles(kv_dense_read_bytes);
    kv_decoder_calibrated = kv_bypassed || m_engine->decoder_calibrated();
    kv_read_cycles = kv_dram_read_cycles + kv_decoder_cycles;

    // Read energy priced on the compressed bytes actually moved.
    kv_read_energy = m_engine->read_energy(kv_read_bytes, &kv_priced);
    if(!kv_priced) {
        kv_unpriced.push_back("KV-cache read fired but [kvcache] kvcache_read_energy is not"
                              " declared");
    }
    if(!kv_decoder_calibrated) {
        kv_unpriced.push_back("KV decoder ran at an undeclared (ideal) throughput --"
                              " [kvcache] kv_decoder_bytes_per_cycle is not declared");
    }
    // Timeline routing by schedule (evaluation discussion 2026-09-03 Sec 7): aggregate
    // folds the read into the DRAM access axis; the tile-level schedules stage it for
    // finalize_layer_timeline()'s KV supply integration instead (never both).
    kv_schedule = m_engine->get_schedule();
    if(kv_schedule == "aggregate") {
        kv_pending = kv_read_cycles > 0.0;
    } else {
        kv_pending = false;
        kv_sched_pending = kv_read_cycles > 0.0;
        const size_t tile_bytes = std::max<size_t>(1, m_engine->get_tile_bytes());
        kv_tiles = std::max<size_t>(1, (kv_dense_read_bytes + tile_bytes - 1)/tile_bytes);
        kv_buffer_tiles = std::max(1u, m_engine->get_buffer_tiles());
    }
}

void stats_t::apply_attention_step(kvcache_t *m_engine, const kv_stream_cost_t &m_cost) {
    // Traffic bookkeeping (shares the KV report plumbing with the proxy modes).
    kv_present = true;
    kv_bytes_per_token = m_engine->get_bytes_per_token();
    kv_context_length = m_engine->get_context_length();
    kv_dense_read_bytes = m_engine->dense_read_bytes();
    kv_read_bytes = m_engine->compressed_read_bytes();
    kv_bypassed = m_engine->bypassed();
    kv_compression_ratio = kv_read_bytes > 0
        ? static_cast<double>(kv_dense_read_bytes)/static_cast<double>(kv_read_bytes) : 1.0;
    kv_profile_reference = m_engine->get_profile_reference();
    kv_schedule = m_engine->get_schedule();
    const size_t tile_bytes = std::max<size_t>(1, m_engine->get_tile_bytes());
    // Tiles per PASS (K or V stream = half the read); blocking needs no granularity.
    kv_tiles = std::max<size_t>(1, (kv_dense_read_bytes/2 + tile_bytes - 1)/tile_bytes);
    kv_buffer_tiles = std::max(1u, m_engine->get_buffer_tiles());
    kv_decoder_cycles = m_engine->decoder_cycles(kv_dense_read_bytes);
    kv_decoder_calibrated = m_engine->bypassed() || m_engine->decoder_calibrated();
    // The dedicated stream is priced from the DRAM component's own declared unit costs --
    // this path never uses the weight-derived estimate or the kvcache_read_energy scalar.
    kv_dram_read_cycles = m_cost.k_supply_cycle + m_cost.v_supply_cycle - kv_decoder_cycles;
    kv_read_cycles = m_cost.k_supply_cycle + m_cost.v_supply_cycle;
    kv_read_energy = m_cost.dram_energy;
    kv_priced = true;
    kv_pending = false;
    kv_sched_pending = false;

    // Attention consumer: compute windows from the declared geometry and rates. An
    // undeclared (0) rate is an IDEAL unit -- zero cycles -- and must be flagged, exactly
    // like the decoder's ideal mode.
    attn_present = true;
    attn_algorithm = m_engine->get_attention_algorithm();
    attn_qk_macs = m_engine->attention_qk_macs();
    attn_softmax_elements = m_engine->attention_softmax_elements();
    const double mac_rate = m_engine->get_attention_macs_per_cycle();
    const double softmax_rate = m_engine->get_softmax_cycles_per_element();
    attn_compute_calibrated = mac_rate > 0.0 && softmax_rate > 0.0;
    attn_qk_cycles = mac_rate > 0.0 ? static_cast<double>(attn_qk_macs)/mac_rate : 0.0;
    attn_av_cycles = attn_qk_cycles;   // AV moves the same MAC volume as QK^T at M=1
    attn_softmax_cycles = static_cast<double>(attn_softmax_elements)*softmax_rate;
    if(mac_rate <= 0.0) {
        kv_unpriced.push_back("attention QK/AV ran at an undeclared (ideal) rate --"
                              " [kvcache] attention_macs_per_cycle is not declared");
    }
    if(softmax_rate <= 0.0) {
        kv_unpriced.push_back("attention softmax ran at an undeclared (ideal) rate --"
                              " [kvcache] softmax_cycles_per_element is not declared");
    }
    attn_supply_k_cycle = m_cost.k_supply_cycle;
    attn_supply_v_cycle = m_cost.v_supply_cycle;
    attn_write_cycle = m_cost.write_cycle;
    bool mac_priced = true, softmax_priced = true;
    attn_compute_energy = m_engine->attention_compute_energy(&mac_priced, &softmax_priced);
    if(!mac_priced) {
        kv_unpriced.push_back("attention QK/AV MACs fired but [kvcache]"
                              " kvcache_attention_mac_energy is not declared");
    }
    if(!softmax_priced) {
        kv_unpriced.push_back("attention softmax fired but [kvcache] kvcache_softmax_energy"
                              " is not declared");
    }
    attn_kv_occupancy_bytes = kv_bytes_per_token*(kv_context_length + 1);
    attn_kv_capacity_bytes = m_engine->get_kv_capacity_bytes();
    attn_pending = true;   // integrated by finalize_layer_timeline()
}

void stats_t::update_network_stats(stats_t *m_source) {
    // L3: mark this object as the network rollup so print_results() states the network
    // axis contract instead of the layer one (see stats_t::network_rollup). The
    // min-reductions below are seeded from the first MAPPED layer folded in
    // (network_rollup_mapped_seeded) -- a min seeded from the default would always
    // report the default, and an SFU-only layer carries no boundary contract at all.
    network_rollup = true;
    // L11: one call == one layer folded into the rollup.
    ++network_timing_layers;

    /* Update PE stats */
    
    // Update the number of computation 
    num_computation += m_source->num_computation;

    // Update computation cost
    layer_latency += m_source->layer_latency;
    busy_cycle_pe += m_source->busy_cycle_pe;
    busy_cycle_pe_array += m_source->busy_cycle_pe_array;
    busy_cycle_global_buffer += m_source->busy_cycle_global_buffer;
    busy_cycle_multi_chip += m_source->busy_cycle_multi_chip;
    busy_cycle_dram += m_source->busy_cycle_dram;
    computation_cycle += m_source->computation_cycle;

    // CE2/CE7: the network-level "Compute-schedule latency" and per-stage busy
    // axes must be the sum of the per-layer values that already feed busy_cycle_*
    // above (layers run serially, so axis totals add like layer_latency does).
    stage_axis_compute += m_source->stage_axis_compute;
    // P4-1: the PE format-IP axis is a real busy axis of the PE stage (it participates
    // in busy_cycle_pe), so it must be summed into the network totals like the other
    // axes -- otherwise network.txt's printed axes cannot reconstruct network PE busy
    // once a config gives the format IP a non-zero unit cost.
    stage_axis_format += m_source->stage_axis_format;
    // The network PE stage is single-buffered as soon as ANY layer's PEs are, so the
    // printed combination rule stays checkable at network scope too.
    pe_double_buffer = pe_double_buffer && m_source->pe_double_buffer;
    // P1-B: bypass is a config property, so it is identical across layers; OR-reduce
    // anyway so the network report states the same contract as its layers.
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        if(m_source->global_buffer_bypass[i]) global_buffer_bypass[i] = true;
    }
    // Phase-5/E3/L8: the clock, energy-basis and DRAM contracts are config properties,
    // identical across MAPPED layers. An SFU-only layer (standalone softmax) carries
    // placeholder contract fields and no boundary contract at all, so it must not
    // overwrite what a fully mapped layer established -- unless it is the only kind of
    // layer the rollup has seen.
    const bool first_mapped_layer = !network_rollup_mapped_seeded;
    // An SFU-only layer's contract fields are placeholders: adopt them only while NO
    // mapped layer has been folded, so a softmax folded first cannot leave "n/a" contract
    // text standing in a rollup that later contains fully mapped layers.
    if(!m_source->sfu_only_layer || first_mapped_layer) {
        authoritative_frequency_mhz = m_source->authoritative_frequency_mhz;
        single_clock_domain = m_source->single_clock_domain;
        clock_domain_note = m_source->clock_domain_note;
        compute_energy_basis = m_source->compute_energy_basis;
        compute_energy_precision_calibrated = m_source->compute_energy_precision_calibrated;
        operand_precision = m_source->operand_precision;
        dram_timing_model = m_source->dram_timing_model;
        array_noc_link_contract = m_source->array_noc_link_contract;
        output_tile_array_resident = m_source->output_tile_array_resident;
        nop_link_contract = m_source->nop_link_contract;
        dram_timing_limits = m_source->dram_timing_limits;
    }
    // L5: layers run serially, so their back-pressure stalls add like layer_latency does.
    for(unsigned stage = 0; stage < 5; ++stage) timeline_stall[stage] += m_source->timeline_stall[stage];
    // L5/L6: a boundary is only as decoupled as its least-decoupled layer, so the network
    // rollup reports the min depth rather than a sum (depths are a contract, not work).
    // SFU-only layers carry no boundary contract, so they do not take part in the min --
    // the seed comes from the first MAPPED layer, wherever it appears in the network.
    if(!m_source->sfu_only_layer) {
        for(unsigned b = 0; b < 4; ++b) {
            timeline_boundary_depth[b] = first_mapped_layer
                ? m_source->timeline_boundary_depth[b]
                : std::min(timeline_boundary_depth[b], m_source->timeline_boundary_depth[b]);
        }
        network_rollup_mapped_seeded = true;
    }
    for(unsigned s = 0; s < 5; ++s) {
        stage_axis_access[s] += m_source->stage_axis_access[s];
        stage_axis_link[s] += m_source->stage_axis_link[s];
        stage_axis_overlap[s] += m_source->stage_axis_overlap[s];
    }

    /* SFU: counters/energy/traffic are work and sum; busy windows sum because layers run
       serially; units/lanes are config properties. layer_latency already contains each
       layer's SFU window (apply_sfu_activation extends it before the rollup). */
    sfu_present = sfu_present || m_source->sfu_present;
    sfu_active = sfu_active || m_source->sfu_active;
    if(m_source->sfu_present) {
        merge_unique_segments(sfu_operation, m_source->sfu_operation, ", ");
    }
    sfu_physical_units = std::max(sfu_physical_units, m_source->sfu_physical_units);
    sfu_lanes = std::max(sfu_lanes, m_source->sfu_lanes);
    sfu_valid_elements += m_source->sfu_valid_elements;
    sfu_operations += m_source->sfu_operations;
    sfu_invocations += m_source->sfu_invocations;
    sfu_chunks += m_source->sfu_chunks;
    sfu_commit_events += m_source->sfu_commit_events;
    if(!m_source->sfu_commit_note.empty()) {
        merge_unique_segments(sfu_commit_note, m_source->sfu_commit_note, "; ");
    }
    if(sfu_profile_reference.empty()) {
        sfu_profile_reference = m_source->sfu_profile_reference;
    }
    sfu_busy_cycle += m_source->sfu_busy_cycle;
    sfu_serial_cycle += m_source->sfu_serial_cycle;
    sfu_hidden_cycle += m_source->sfu_hidden_cycle;
    sfu_stall_cycle += m_source->sfu_stall_cycle;
    // Queue depth is a config property; the max survives the 0-default of layers that
    // recorded no SFU event.
    sfu_queue_depth = std::max(sfu_queue_depth, m_source->sfu_queue_depth);
    sfu_stream_active = sfu_stream_active || m_source->sfu_stream_active;
    if(m_source->sfu_stream_active) {
        merge_unique_segments(sfu_stream_residency, m_source->sfu_stream_residency, "; ");
        sfu_stream_ingress_bytes += m_source->sfu_stream_ingress_bytes;
        sfu_stream_egress_bytes += m_source->sfu_stream_egress_bytes;
        sfu_stream_ingress_cycle += m_source->sfu_stream_ingress_cycle;
        sfu_stream_egress_cycle += m_source->sfu_stream_egress_cycle;
    }
    sfu_tail_lane_utilization = std::max(sfu_tail_lane_utilization,
                                         m_source->sfu_tail_lane_utilization);
    sfu_ingress_elements += m_source->sfu_ingress_elements;
    sfu_egress_elements += m_source->sfu_egress_elements;
    sfu_ingress_bytes += m_source->sfu_ingress_bytes;
    sfu_egress_bytes += m_source->sfu_egress_bytes;
    sfu_ingress_transactions += m_source->sfu_ingress_transactions;
    sfu_egress_transactions += m_source->sfu_egress_transactions;
    sfu_op_energy += m_source->sfu_op_energy;
    sfu_read_energy += m_source->sfu_read_energy;
    sfu_write_energy += m_source->sfu_write_energy;
    sfu_setup_energy += m_source->sfu_setup_energy;
    sfu_static_energy += m_source->sfu_static_energy;
    sfu_timing_calibrated = sfu_timing_calibrated && m_source->sfu_timing_calibrated;
    for(unsigned i = 0; i < m_source->sfu_unpriced_events.size(); ++i) {
        if(std::find(sfu_unpriced_events.begin(), sfu_unpriced_events.end(),
                     m_source->sfu_unpriced_events[i]) == sfu_unpriced_events.end()) {
            sfu_unpriced_events.push_back(m_source->sfu_unpriced_events[i]);
        }
    }
    if(!m_source->sfu_contract_note.empty()) {
        merge_unique_segments(sfu_contract_note, m_source->sfu_contract_note, "; ");
    }
    // Decompression: counters/energy/breakdown are work and sum; contract fields are
    // config properties (max survives the 0-default of layers with no decomp event).
    decomp_present = decomp_present || m_source->decomp_present;
    decomp_active = decomp_active || m_source->decomp_active;
    decomp_bypassed = decomp_bypassed || m_source->decomp_bypassed;
    decomp_dense_weight_bytes += m_source->decomp_dense_weight_bytes;
    decomp_compressed_weight_bytes += m_source->decomp_compressed_weight_bytes;
    decomp_tiles += m_source->decomp_tiles;
    decomp_decoder_cycles += m_source->decomp_decoder_cycles;
    decomp_exposed_cycle += m_source->decomp_exposed_cycle;
    decomp_hidden_cycle += m_source->decomp_hidden_cycle;
    decomp_stall_cycle += m_source->decomp_stall_cycle;
    decomp_dram_weight_saved_cycle += m_source->decomp_dram_weight_saved_cycle;
    decomp_decoder_energy += m_source->decomp_decoder_energy;
    decomp_static_energy += m_source->decomp_static_energy;
    decomp_breakdown_dram += m_source->decomp_breakdown_dram;
    decomp_breakdown_decode += m_source->decomp_breakdown_decode;
    decomp_breakdown_compute += m_source->decomp_breakdown_compute;
    decomp_breakdown_stall += m_source->decomp_breakdown_stall;
    decomp_queue_depth = std::max(decomp_queue_depth, m_source->decomp_queue_depth);
    decomp_overlap = decomp_overlap || m_source->decomp_overlap;
    decomp_bypassed_tiles += m_source->decomp_bypassed_tiles;
    decomp_tile_ratio_cv = std::max(decomp_tile_ratio_cv, m_source->decomp_tile_ratio_cv);
    decomp_output_buffer_tiles = std::max(decomp_output_buffer_tiles,
                                          m_source->decomp_output_buffer_tiles);
    decomp_sink_cycles += m_source->decomp_sink_cycles;
    // Effective ratio at network scope: total dense / total compressed.
    if(decomp_compressed_weight_bytes > 0) {
        decomp_effective_ratio = static_cast<double>(decomp_dense_weight_bytes)/
                                 static_cast<double>(decomp_compressed_weight_bytes);
    }
    decomp_timing_calibrated = decomp_timing_calibrated && m_source->decomp_timing_calibrated;
    for(unsigned i = 0; i < m_source->decomp_unpriced_events.size(); ++i) {
        if(std::find(decomp_unpriced_events.begin(), decomp_unpriced_events.end(),
                     m_source->decomp_unpriced_events[i]) == decomp_unpriced_events.end()) {
            decomp_unpriced_events.push_back(m_source->decomp_unpriced_events[i]);
        }
    }
    if(decomp_profile_reference.empty()) {
        decomp_profile_reference = m_source->decomp_profile_reference;
    }
    // KV-cache read: bytes/cycles/energy are per-layer work and sum; the per-token footprint
    // and context length are config properties (max survives the 0-default of layers with no
    // KV read).
    kv_present = kv_present || m_source->kv_present;
    if(m_source->kv_pending) kv_pending = true;
    kv_read_bytes += m_source->kv_read_bytes;
    kv_dense_read_bytes += m_source->kv_dense_read_bytes;
    kv_read_cycles += m_source->kv_read_cycles;
    kv_dram_read_cycles += m_source->kv_dram_read_cycles;
    kv_decoder_cycles += m_source->kv_decoder_cycles;
    kv_read_energy += m_source->kv_read_energy;
    kv_bytes_per_token = std::max(kv_bytes_per_token, m_source->kv_bytes_per_token);
    kv_context_length = std::max(kv_context_length, m_source->kv_context_length);
    kv_bypassed = kv_bypassed || m_source->kv_bypassed;
    if(kv_read_bytes > 0) {
        kv_compression_ratio = static_cast<double>(kv_dense_read_bytes)/
                               static_cast<double>(kv_read_bytes);
    }
    kv_priced = kv_priced && m_source->kv_priced;
    kv_decoder_calibrated = kv_decoder_calibrated && m_source->kv_decoder_calibrated;
    // Schedule is a config property (identical across layers); exposure/stall are work.
    if(kv_schedule.empty()) kv_schedule = m_source->kv_schedule;
    kv_tiles += m_source->kv_tiles;
    kv_buffer_tiles = std::max(kv_buffer_tiles, m_source->kv_buffer_tiles);
    kv_exposed_cycle += m_source->kv_exposed_cycle;
    kv_hidden_cycle += m_source->kv_hidden_cycle;
    kv_stall_cycle += m_source->kv_stall_cycle;
    // Attention consumer step: counters/cycles/energy are work and sum; the algorithm and
    // capacity are config properties.
    attn_present = attn_present || m_source->attn_present;
    if(attn_algorithm.empty()) attn_algorithm = m_source->attn_algorithm;
    attn_qk_macs += m_source->attn_qk_macs;
    attn_softmax_elements += m_source->attn_softmax_elements;
    attn_qk_cycles += m_source->attn_qk_cycles;
    attn_av_cycles += m_source->attn_av_cycles;
    attn_softmax_cycles += m_source->attn_softmax_cycles;
    attn_supply_k_cycle += m_source->attn_supply_k_cycle;
    attn_supply_v_cycle += m_source->attn_supply_v_cycle;
    attn_write_cycle += m_source->attn_write_cycle;
    attn_step_latency += m_source->attn_step_latency;
    attn_exposed_cycle += m_source->attn_exposed_cycle;
    attn_hidden_cycle += m_source->attn_hidden_cycle;
    attn_stall_cycle += m_source->attn_stall_cycle;
    attn_compute_energy += m_source->attn_compute_energy;
    attn_compute_calibrated = attn_compute_calibrated && m_source->attn_compute_calibrated;
    attn_kv_occupancy_bytes += m_source->attn_kv_occupancy_bytes;
    attn_kv_capacity_bytes = std::max(attn_kv_capacity_bytes,
                                      m_source->attn_kv_capacity_bytes);
    for(unsigned i = 0; i < m_source->kv_unpriced.size(); ++i) {
        if(std::find(kv_unpriced.begin(), kv_unpriced.end(), m_source->kv_unpriced[i]) ==
           kv_unpriced.end()) {
            kv_unpriced.push_back(m_source->kv_unpriced[i]);
        }
    }
    if(kv_profile_reference.empty()) {
        kv_profile_reference = m_source->kv_profile_reference;
    }
    if(m_source->graph_residency_applied) graph_residency_applied = true;
    if(!m_source->graph_residency_note.empty()) {
        merge_unique_segments(graph_residency_note, m_source->graph_residency_note, "; ");
    }
    if(m_source->activation_unmodeled) {
        mark_unmodeled_activation(m_source->activation_unmodeled_note);
    }

    max_computation_cycle += m_source->max_computation_cycle;
    min_computation_cycle += m_source->min_computation_cycle;
    avg_computation_cycle += m_source->avg_computation_cycle;
    computation_energy += m_source->computation_energy;
    reduction_energy += m_source->reduction_energy;
    weight_fold_energy += m_source->weight_fold_energy;
    pe_array_accumulator_energy += m_source->pe_array_accumulator_energy;
    layer_setup_energy += m_source->layer_setup_energy;
    stripe_transition_energy += m_source->stripe_transition_energy;
    weight_fold_events += m_source->weight_fold_events;
    layer_setup_events += m_source->layer_setup_events;
    stripe_transition_events += m_source->stripe_transition_events;
    accumulator_reload_bytes += m_source->accumulator_reload_bytes;
    accumulator_spill_bytes += m_source->accumulator_spill_bytes;
    accumulator_create_events += m_source->accumulator_create_events;
    accumulator_retained_events += m_source->accumulator_retained_events;
    output_cast_bytes += m_source->output_cast_bytes;
    accumulator_energy += m_source->accumulator_energy;
    output_cast_energy += m_source->output_cast_energy;
    output_cast_cycle += m_source->output_cast_cycle;
    row_activation_events += m_source->row_activation_events;
    format_payload_events += m_source->format_payload_events;
    format_metadata_events += m_source->format_metadata_events;
    reduction_additions += m_source->reduction_additions;
    mac_busy_cycle += m_source->mac_busy_cycle;
    mac_available_cycle += m_source->mac_available_cycle;
    utilization_mac = calculate_time_based_mac_utilization(mac_busy_cycle, mac_available_cycle);
    if(utilization_mac > 1.0 + 1e-9) {
        std::cerr << "Error: network MAC busy cycles exceed physical MAC capacity" << std::endl;
        exit(1);
    }
    utilization_mac = std::min(1.0, utilization_mac);


    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        // Update the number of request to the local buffer
        num_request_pe[i] += m_source->num_request_pe[i];

        // Update the number of transfer from the local buffer
        num_data_transfer_pe[i] += m_source->num_data_transfer_pe[i];

        // Update access cost of the computing unit
        access_cycle_mac[i] += m_source->access_cycle_mac[i];
        access_energy_mac[i] += m_source->access_energy_mac[i];

        max_access_cycle_mac[i] += m_source->max_access_cycle_mac[i];
        min_access_cycle_mac[i] += m_source->min_access_cycle_mac[i];
        avg_access_cycle_mac[i] += m_source->avg_access_cycle_mac[i];

        // Update access cost of the local buffer
        access_cycle_lb[i] += m_source->access_cycle_lb[i];
        access_energy_lb[i] += m_source->access_energy_lb[i];
        static_energy_pe[i] += m_source->static_energy_pe[i];
        static_energy_pe_array[i] += m_source->static_energy_pe_array[i];

        max_access_cycle_lb[i] += m_source->max_access_cycle_lb[i];
        min_access_cycle_lb[i] += m_source->min_access_cycle_lb[i];
        avg_access_cycle_lb[i] += m_source->avg_access_cycle_lb[i];


        format_cycle_pe[i] += m_source->format_cycle_pe[i];
        format_energy_pe[i] += m_source->format_energy_pe[i];
        // Update transfer cost between the local buffer and computing units
        transfer_cycle_pe[i] += m_source->transfer_cycle_pe[i];
        transfer_energy_pe[i] += m_source->transfer_energy_pe[i];
        payload_link_transactions_pe[i] += m_source->payload_link_transactions_pe[i];
        metadata_link_transactions_pe[i] += m_source->metadata_link_transactions_pe[i];
        storage_link_transactions_pe[i] += m_source->storage_link_transactions_pe[i];

        // Update overlapped cycle between the local buffer and computing units
        cycle_mac_lb[i] += m_source->cycle_mac_lb[i];

        /* Update PE array stats */

        // Update the number of request to the PE array
        num_request_pe_array[i] += m_source->num_request_pe_array[i];

        // Update the number of data transfer from the PE array
        num_data_transfer_pe_array[i] += m_source->num_data_transfer_pe_array[i];
        payload_link_transactions_pe_array[i] += m_source->payload_link_transactions_pe_array[i];
        metadata_link_transactions_pe_array[i] += m_source->metadata_link_transactions_pe_array[i];
        storage_link_transactions_pe_array[i] += m_source->storage_link_transactions_pe_array[i];

        payload_link_transactions_global_buffer[i] += m_source->payload_link_transactions_global_buffer[i];
        metadata_link_transactions_global_buffer[i] += m_source->metadata_link_transactions_global_buffer[i];
        storage_link_transactions_global_buffer[i] += m_source->storage_link_transactions_global_buffer[i];
        // Update access cost of the PE array
        access_cycle_pe_array[i] += m_source->access_cycle_pe_array[i];
        access_energy_pe_array[i] += m_source->access_energy_pe_array[i];

        // Update transfer cost between the PE array and the local buffers
        transfer_cycle_pe_array[i] += m_source->transfer_cycle_pe_array[i];
        cycle_temporal_pe_array[i] += m_source->cycle_temporal_pe_array[i];
        transfer_energy_pe_array[i] += m_source->transfer_energy_pe_array[i];
        if(i == 0) fold_fill_cycle_pe_array += m_source->fold_fill_cycle_pe_array;

        /* Update global buffer stats */

        // Update the number of request to the global buffer
        num_request_global_buffer[i] += m_source->num_request_global_buffer[i];
        
        // Update the number of data transfer from the global buffer
        num_data_transfer_global_buffer[i] += m_source->num_data_transfer_global_buffer[i];
        if(i == data_type_t::OUTPUT) psum_writeback_events_global_buffer += m_source->psum_writeback_events_global_buffer;

        // Update access cost of the global buffer
        access_cycle_global_buffer[i] += m_source->access_cycle_global_buffer[i];
        access_energy_global_buffer[i] += m_source->access_energy_global_buffer[i];
        static_energy_global_buffer[i] += m_source->static_energy_global_buffer[i];

        // Update transfer cost between the global buffer and PE array
        transfer_cycle_global_buffer[i] += m_source->transfer_cycle_global_buffer[i];
        transfer_energy_global_buffer[i] += m_source->transfer_energy_global_buffer[i];

        // Update overlapped cycle between the global buffer and PE array
        cycle_pe_array_global_buffer[i] += m_source->cycle_pe_array_global_buffer[i];

        /* Update Chip-level processor stats */

        // Update the number of request to the chip-level processor
        num_request_multi_chip[i] += m_source->num_request_multi_chip[i];

        // Update the number of data transfer from the chip-level processor
        num_data_transfer_multi_chip[i] += m_source->num_data_transfer_multi_chip[i];
        payload_link_transactions_multi_chip[i] += m_source->payload_link_transactions_multi_chip[i];
        metadata_link_transactions_multi_chip[i] += m_source->metadata_link_transactions_multi_chip[i];
        storage_link_transactions_multi_chip[i] += m_source->storage_link_transactions_multi_chip[i];

        // Update access cost of the chip-level processor
        access_cycle_multi_chip[i] += m_source->access_cycle_multi_chip[i];
        access_energy_multi_chip[i] += m_source->access_energy_multi_chip[i];

        // Update transfer cost between the chip-level processor to the global buffer
        transfer_cycle_multi_chip[i] += m_source->transfer_cycle_multi_chip[i];
        transfer_energy_multi_chip[i] += m_source->transfer_energy_multi_chip[i];
        static_energy_multi_chip[i] += m_source->static_energy_multi_chip[i];

        /* Update off-chip memory stats */

        // Update the number of request to the off-chip memory
        num_request_dram[i] += m_source->num_request_dram[i];

        // Update the number of data transfer from the off-chip memory
        num_data_transfer_dram[i] += m_source->num_data_transfer_dram[i];
        payload_link_transactions_dram[i] += m_source->payload_link_transactions_dram[i];
        metadata_link_transactions_dram[i] += m_source->metadata_link_transactions_dram[i];
        storage_link_transactions_dram[i] += m_source->storage_link_transactions_dram[i];

        // Update access cost of the off-chip memory
        access_cycle_dram[i] += m_source->access_cycle_dram[i];
        access_energy_dram[i] += m_source->access_energy_dram[i];

        // Update transfer cost between the off-chip memory and chip-level processor
        transfer_cycle_dram[i] += m_source->transfer_cycle_dram[i];
        transfer_energy_dram[i] += m_source->transfer_energy_dram[i];

        // Update overlapped cycle between the off-chip memory and chip-level processor
        cycle_chip_dram[i] += m_source->cycle_chip_dram[i];
    }

    // P0/CE3 identity gate: "Compute-schedule latency" must always equal
    // "Computation cycle" + "Fold fill cycle" so the golden gate can trust the
    // official metric printed to network.txt (not just layer_*.txt).
    const double compute_schedule_identity = computation_cycle + fold_fill_cycle_pe_array;
    const double identity_tolerance = 1e-6*std::max(1.0, std::fabs(compute_schedule_identity));
    if(std::fabs(stage_axis_compute - compute_schedule_identity) > identity_tolerance) {
        std::cerr << "Error: network Compute-schedule latency diverged from "
                  << "Computation cycle + Fold fill cycle (" << stage_axis_compute
                  << " != " << compute_schedule_identity << ")" << std::endl;
        exit(1);
    }
}

// E1: energy summary -- component subtotals and the layer total.
//
// WHY THIS EXISTS. Three energy accumulators were never printed at all, and even the printed
// ones could not be added up into a layer total: there was no statement of which accumulators
// belong to which component, so a reader could not tell whether summing them double-counted a
// shared event. Without a total there is also nothing for a checker to re-derive, which is why
// a formula change could not be regression-tested.
//
// NO-DOUBLE-COUNTING BOUNDARY (the contract that makes the sum well defined). Every transfer
// in this model charges three DISTINCT physical resources, and each one is billed to the
// component that OWNS it:
//   * the source buffer's read port  -> the SOURCE component's access energy
//   * the link/fabric crossed        -> the component that owns that link
//   * the destination buffer's write -> the DESTINATION component's access energy
// So a GLB->PE-array transfer charges the GLB's read, the GLB<->array link (owned by the GLB),
// and the PE-array temporal buffer's write -- three different resources, three different
// owners, no event counted twice. Summing the subtotals is therefore exact, not an estimate.
// Leakage is modeled per physical component and is kept on its own axis, because it is a
// function of the layer's wall-clock rather than of any event count.
void stats_t::print_energy_summary(std::ofstream &m_output_file) {
    auto sum_types = [](const std::vector<double> &values) {
        double total = 0.0;
        for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) total += values[i];
        return total;
    };
    // MAC: the compute itself, the MAC register file, and the format IP on the LB<->MAC path.
    // E4: the reduction tree is its own axis now, but it is still the MAC component's energy --
    // include it in the subtotal so the layer total stays complete.
    // RE1: the PE-side accumulator energy (when the accumulator is IN the PE) and the final
    // output cast belong to the MAC/PE datapath component. With edge_accumulation the accumulator
    // energy was handed to the PE array instead and appears in that row.
    const double mac_dynamic = computation_energy + reduction_energy +
                               sum_types(access_energy_mac) + sum_types(format_energy_pe) +
                               accumulator_energy;
    // PE: the local buffer's own ports and the LB<->MAC bus. Leakage is modeled for the whole
    // PE (MAC included), so it is reported on this row.
    const double pe_dynamic = sum_types(access_energy_lb) + sum_types(transfer_energy_pe);
    const double pe_static = sum_types(static_energy_pe);
    // PE array: its temporal buffer's ports and the array's own distribution fabric.
    // E5: the fold/setup control energy is the PE array's own activity (weight reload, schedule
    // setup), so it belongs to that component's subtotal rather than sitting outside the total.
    const double pe_array_dynamic = sum_types(access_energy_pe_array) +
                                    sum_types(transfer_energy_pe_array) +
                                    weight_fold_energy + layer_setup_energy +
                                    stripe_transition_energy +
                                    pe_array_accumulator_energy;
    const double pe_array_static = sum_types(static_energy_pe_array);
    // GLB: its SRAM ports (the multi-chip fill write is already folded in by
    // merge_global_buffer_fill()) and the GLB<->PE-array link.
    const double global_buffer_dynamic = sum_types(access_energy_global_buffer) +
                                         sum_types(fill_access_energy_global_buffer) +
                                         sum_types(transfer_energy_global_buffer);
    const double global_buffer_static = sum_types(static_energy_global_buffer);
    // Multi-chip: its temporal buffer's ports and the NoP fabric.
    // RE1: the final output cast/pack is committed at the off-chip store, so it belongs to the
    // multi-chip component's subtotal.
    const double multi_chip_dynamic = sum_types(access_energy_multi_chip) +
                                      sum_types(transfer_energy_multi_chip) + output_cast_energy;
    const double multi_chip_static = sum_types(static_energy_multi_chip);
    // DRAM: its device access energy and the off-chip link plus row activations, plus the
    // KV-cache read energy (zero without a [kvcache] section, so the dense baseline is
    // unchanged). The KV read is DRAM traffic, so it belongs to the DRAM component subtotal.
    const double dram_dynamic = sum_types(access_energy_dram) + sum_types(transfer_energy_dram)
                                + kv_read_energy;
    // SFU: the activation operation itself plus its internal ingress/egress and setup, all
    // generated by the same invocation event (plan principle 6); leakage covers every
    // physical SFU over the final layer window. Zero -- and the row absent -- when no
    // [sfu] section is configured, so legacy totals are bit-identical.
    const double sfu_dynamic = sfu_op_energy + sfu_read_energy + sfu_write_energy +
                               sfu_setup_energy;
    const double sfu_static = sfu_static_energy;
    // Decompression: decoder dynamic energy + always-on engine leakage. Zero -- and out
    // of the total -- when no [decomp] section is configured (dense baseline unchanged).
    const double decomp_dynamic = decomp_decoder_energy;
    const double decomp_static_e = decomp_static_energy;

    // Attention consumer step: QK/AV + softmax energy (zero -- and flagged UNPRICED --
    // when the [kvcache] cost keys are undeclared). Its KV DRAM energy is already inside
    // dram_dynamic via kv_read_energy.
    const double dynamic_total = mac_dynamic + pe_dynamic + pe_array_dynamic +
                                 global_buffer_dynamic + multi_chip_dynamic + dram_dynamic +
                                 sfu_dynamic + decomp_dynamic + attn_compute_energy;
    const double static_total = pe_static + pe_array_static + global_buffer_static +
                                multi_chip_static + sfu_static + decomp_static_e;

    const std::vector<std::string> unpriced = unpriced_active_events();
    m_output_file << "============= Energy summary ============" << std::endl;
    // E3: name the key the MAC energy came from, and say plainly when it is not precision
    // calibrated -- otherwise the same scalar silently prices INT4 and FP16 identically.
    m_output_file << "MAC energy basis      : " << compute_energy_basis;
    if(compute_energy_precision_calibrated) {
        m_output_file << " (declared for " << operand_precision << ")";
    } else {
        m_output_file << " -- NOT calibrated for " << operand_precision
                      << " (declare mac_energy_<input>_<weight>_<accumulator> for it)";
    }
    m_output_file << std::endl;
    // E7: state the unit and its provenance right above the numbers. A normalized fixture's
    // absolute total is not meaningful, and nothing else in the report said so.
    m_output_file << "Energy unit           : " << energy_units().describe() << std::endl;
    const char *rows[7] = {"MAC (compute+format)", "PE (local buffer)", "PE array",
                           "Global buffer", "Multi-chip (NoP)", "DRAM", "SFU (activation)"};
    const double row_dynamic[7] = {mac_dynamic, pe_dynamic, pe_array_dynamic,
                                   global_buffer_dynamic, multi_chip_dynamic, dram_dynamic,
                                   sfu_dynamic};
    const double row_static[7] = {0.0, pe_static, pe_array_static, global_buffer_static,
                                  multi_chip_static, 0.0, sfu_static};
    // The SFU row is opt-in: shown only for runs that model an SFU (or folded one in at
    // network scope), so a legacy config's six-row breakdown is untouched.
    const unsigned printed_rows = (energy_cost_schema().sfu_declared || sfu_present) ? 7 : 6;
    m_output_file << std::left << std::setw(24) << "Component" << std::right
                  << std::setw(17) << "dynamic" << std::setw(17) << "static"
                  << std::setw(17) << "total" << std::endl;
    for(unsigned i = 0; i < printed_rows; ++i) {
        const bool has_energy = row_dynamic[i] != 0.0 || row_static[i] != 0.0;
        m_output_file << " * " << std::left << std::setw(21) << rows[i] << std::right
                      << std::setw(17) << std::setprecision(2) << row_dynamic[i]
                      << std::setw(17) << row_static[i]
                      << std::setw(17) << row_dynamic[i] + row_static[i]
                      // RE5: the annotation now comes from the component's DECLARATION state, not
                      // from the subtotal being zero. A bare 0 used to mean any of "not modeled",
                      // "half-priced", "deliberately free" and "never exercised".
                      << energy_cost_schema().annotate(static_cast<energy_component_t>(i),
                                                       has_energy)
                      << std::endl;
    }
    // E12: at network scope, state how many layers the total actually covers. Unsupported
    // layers (pooling, activation, normalization) are excluded from accelerator accounting
    // entirely, so their energy is ABSENT from this total -- not estimated, not zero. Without
    // this line a partial rollup reads as an end-to-end network energy, exactly the way the
    // timing rollup did before L11.
    if(network_rollup) {
        const unsigned total_layers = network_timing_layers + excluded_timing_layers;
        m_output_file << "Energy scope          : " << network_timing_layers << " of "
                      << total_layers << " layers";
        if(excluded_timing_layers > 0) {
            m_output_file << "  (PARTIAL: " << excluded_timing_layers
                          << " unsupported layer(s) excluded; this is NOT an end-to-end network"
                          << " energy)";
        } else {
            m_output_file << "  (complete: every layer is supported)";
        }
        m_output_file << std::endl;
    }
    // RE7: at network scope these are network totals, not layer ones. The rollup printed the
    // same "Layer ..." labels as a single-layer report, so a network file's total read as one
    // layer's -- and any checker matching on the label could not tell the two scopes apart.
    const char *scope = network_rollup ? "Network" : "Layer";
    m_output_file << std::left << std::setw(22) << (std::string(scope) + " dynamic energy")
                  << std::right << ":" << std::setw(18) << std::setprecision(2)
                  << dynamic_total << " " << energy_units().label() << std::endl;
    m_output_file << std::left << std::setw(22) << (std::string(scope) + " static energy")
                  << std::right << ":" << std::setw(18) << std::setprecision(2)
                  << static_total << " " << energy_units().label() << std::endl;
    m_output_file << std::left << std::setw(22) << (std::string(scope) + " total energy")
                  << std::right << ":" << std::setw(18) << std::setprecision(2)
                  << dynamic_total + static_total << " " << energy_units().label();
    // RE2/E20-1: say plainly when the total's SCALE is not established, so it is never read as
    // measured. Three independent reasons, any one of which is enough: the unit is not declared
    // absolute with provenance, the compute cost is not calibrated for the running precision, or an
    // event fired whose price nobody declared (E20-2). The third is the one that let a fixture with
    // externally derived SRAM and DRAM costs still be missing its setup, accumulator and cast
    // charges while claiming absolute pJ.
    if(!energy_units().is_absolute() || !compute_energy_precision_calibrated ||
       !unpriced.empty() || energy_cost_schema().undercounted() > 0) {
        m_output_file << "   (ESTIMATED: not calibrated for absolute comparison)";
    }
    m_output_file << std::endl;
    // RE5: an UNDERCOUNT warning keyed off the declaration. A component priced 0 on purpose does
    // not make the total suspect; one whose cost is missing does, whether or not it came out at 0.
    if(!unpriced.empty()) {
        m_output_file << "Unpriced active events: " << unpriced.size()
                      << " -- these FIRED and the config declares no cost for them, so the total"
                      << " above is missing a charge of undeclared size" << std::endl;
        for(unsigned i = 0; i < unpriced.size(); ++i) {
            m_output_file << "  * " << unpriced[i] << std::endl;
        }
        m_output_file << "  (declare the key as 0 to state a modeled zero; an ABSENT key is not the"
                      << " same statement)" << std::endl;
    }
    const unsigned undercounted = energy_cost_schema().undercounted();
    if(undercounted > 0) {
        m_output_file << "Uncosted components   : " << undercounted << " of "
                      << energy_cost_schema().components_in_scope()
                      << " have a missing or absent energy cost, so this"
                      << " total is an UNDERCOUNT (see the rows marked above)" << std::endl;
    }
    // Plan compatibility policy: a nonlinear activation that executed with no [sfu]
    // section keeps the legacy numbers, but its energy is ABSENT from this total -- scope,
    // stated next to the total it qualifies.
    if(activation_unmodeled) {
        m_output_file << "Activation scope      : NOT MODELED -- " << activation_unmodeled_note
                      << std::endl;
    }
    m_output_file << std::endl;

    print_power_summary(m_output_file, dynamic_total, static_total);
}

// E20-1/E20-2: which events fired without a declared price.
//
// The rule is uniform: an event is ACTIVE if its own counter is non-zero, and PRICED if the config
// declares the key -- with any value, zero included. An absent key and a declared zero are
// different statements and are treated as such; that distinction is why a count is needed at all,
// since both produce an energy of 0.
//
// PE local-buffer leakage is the one conditional entry: it counts as active only when the PE's own
// `static_energy` is declared, i.e. when the config has put leakage in scope at all. A config that
// models no leakage is not missing a leakage charge.
std::vector<std::string> stats_t::unpriced_active_events() const {
    std::vector<std::string> out;
    const energy_cost_schema_t &schema = energy_cost_schema();
    struct entry_t { bool active; const char *event; const char *section; const char *key; };
    const double accumulator_bytes = static_cast<double>(accumulator_reload_bytes) +
                                     static_cast<double>(accumulator_spill_bytes);
    const entry_t entries[] = {
        { layer_setup_events > 0.0,   "layer setup",                "pe_array",
          "layer_setup_energy" },
        { stripe_transition_events > 0.0, "stripe transition",      "pe_array",
          "stripe_transition_energy" },
        { weight_fold_events > 0.0,   "weight fold",                "pe_array",
          "weight_fold_fill_energy" },
        { accumulator_bytes > 0.0,    "accumulator reload/spill",   "pe_array",
          "accumulator_spill_energy" },
        { output_cast_bytes > 0,      "final output cast",          "multi_chip",
          "output_cast_energy" },
        { row_activation_events > 0,  "DRAM row activation",        "dram",
          "row_miss_energy" },
        { format_payload_events > 0,  "Format-IP payload",          "pe_array",
          "format_payload_energy" },
        { format_metadata_events > 0, "Format-IP metadata",         "pe_array",
          "format_metadata_energy" },
        { reduction_additions > 0.0,  "array reduction additions",  "pe_array",
          "adder_energy" },
        { schema.is_declared("pe_array", "static_energy"),
                                      "PE local-buffer leakage",    "pe_array",
          "lb_static_energy" },
    };
    for(unsigned i = 0; i < sizeof(entries)/sizeof(entries[0]); ++i) {
        if(!entries[i].active) continue;
        if(schema.is_declared(entries[i].section, entries[i].key)) continue;
        out.push_back(std::string(entries[i].event) + " fired but [" + entries[i].section +
                      "] " + entries[i].key + " is not declared");
    }
    // SFU events carry per-operation keys (sfu_op_energy_<op>), so the SFU records its own
    // unpriced-active list at invocation time; an entry here blocks absolute energy/power
    // exactly like the static entries above (plan: unpriced active SFU events disqualify
    // absolute results).
    for(unsigned i = 0; i < sfu_unpriced_events.size(); ++i) {
        out.push_back(sfu_unpriced_events[i]);
    }
    // Decompression decoder events fired without a declared cost key block absolute
    // energy/power exactly like every other unpriced event.
    for(unsigned i = 0; i < decomp_unpriced_events.size(); ++i) {
        out.push_back(decomp_unpriced_events[i]);
    }
    // KV-cache read events fired without a declared cost key block absolute energy/power too.
    for(unsigned i = 0; i < kv_unpriced.size(); ++i) {
        out.push_back(kv_unpriced[i]);
    }
    return out;
}

// Phase-5: average power, EDP/ED2P, and an explicit statement of what power does NOT include.
//
// The contract is deliberately small and complete rather than large and partly founded:
//
//   time      = critical-path cycles / authoritative frequency
//   power     = energy / time
//   EDP       = energy x time,   ED2P = energy x time^2
//
// Three preconditions, each reported rather than assumed:
//   * ONE CLOCK. The timeline is a single shared cycle axis, so cycles convert to seconds only if
//     every modeled component runs on the same clock (see stats_t::single_clock_domain). A
//     mixed-domain config gets "unsupported", not a number divided by an arbitrary domain.
//   * ABSOLUTE ENERGY. A normalized energy fixture has no watts. Power is therefore only reported
//     when energy_unit = pJ (see E7); otherwise it says why not.
//   * AVERAGE ONLY. Peak power needs the concurrency of individual events, which this model does
//     not resolve -- the per-tile timeline places stage-level work, not per-event overlap. Peak is
//     explicitly unsupported instead of being approximated by a sum of maxima.
//
// static_energy is pJ/CYCLE in the config and is already integrated over the layer's final
// wall-clock (see the leakage rescale in finalize_layer_timeline()), so dividing it by the same
// time window gives leakage power without double counting the duration.
void stats_t::print_power_summary(std::ofstream &m_output_file, double dynamic_total,
                                 double static_total) {
    m_output_file << "============= Power summary =============" << std::endl;
    m_output_file << "Clock domain          : " << clock_domain_note << std::endl;
    if(single_clock_domain) {
        m_output_file << "Authoritative clock   :" << std::setw(15) << std::setprecision(1)
                      << authoritative_frequency_mhz << " MHz" << std::endl;
    }
    // What power never includes here, stated next to the numbers so a chip-level comparison is
    // not attempted against a core-datapath figure.
    m_output_file << "Not included          : DRAM background/refresh and I/O termination, clock"
                  << " network, controller/DMA, PHY; core datapath only" << std::endl;
    m_output_file << "Peak power            : unsupported (needs per-event concurrency; this model"
                  << " resolves stage-level overlap only)" << std::endl;

    if(!single_clock_domain) {
        m_output_file << "Average power         : unsupported for this config (see the clock domain"
                      << " above)" << std::endl;
        m_output_file << std::endl;
        return;
    }
    // RE2: absolute power requires the energy to BE absolute -- a declared pJ unit with declared
    // provenance -- and it requires the compute cost to be calibrated for the precision actually
    // running. A fallback scalar shared by every precision is not a calibrated compute cost, so a
    // config using it does not qualify however its unit is labelled.
    if(!energy_units().is_absolute()) {
        m_output_file << "Average power         : unsupported (" << energy_units().label()
                      << " energy) -- " << energy_units().calibration_note() << std::endl;
        m_output_file << std::endl;
        return;
    }
    if(!compute_energy_precision_calibrated) {
        m_output_file << "Average power         : unsupported -- the compute cost is not calibrated"
                      << " for the precision in use (" << operand_precision << "; basis "
                      << compute_energy_basis << ")" << std::endl;
        m_output_file << std::endl;
        return;
    }
    // E20-1: absolute power additionally requires that nothing in the numerator is missing. A
    // component whose cost is half declared (RE5's PARTIAL/NOT_MODELED) makes the total an
    // undercount; an event that fired unpriced makes it an undercount of unstated size. Either way
    // the division is arithmetic, not a wattage.
    const unsigned undercounted = energy_cost_schema().undercounted();
    const std::vector<std::string> unpriced = unpriced_active_events();
    if(undercounted > 0 || !unpriced.empty()) {
        m_output_file << "Average power         : unsupported -- the energy total is an UNDERCOUNT ("
                      << undercounted << " component(s) with a missing cost, " << unpriced.size()
                      << " active event(s) with no declared price), so energy/time is not a power"
                      << std::endl;
        if(!unpriced.empty()) {
            m_output_file << "  first unpriced        : " << unpriced[0] << std::endl;
        }
        m_output_file << std::endl;
        return;
    }

    const double seconds = layer_latency/(authoritative_frequency_mhz*1.0e6);
    if(seconds <= 0.0) {
        m_output_file << "Average power         : unsupported (non-positive elapsed time)"
                      << std::endl;
        m_output_file << std::endl;
        return;
    }
    // Energy is pJ; pJ/s is pW, so scale to mW for a readable figure.
    const double to_milliwatt = 1.0e-9;
    const double total_energy = dynamic_total + static_total;
    m_output_file << "Elapsed time          :" << std::setw(15) << std::setprecision(6)
                  << seconds*1.0e3 << " ms  (" << std::setprecision(1) << layer_latency
                  << " cycles / " << authoritative_frequency_mhz << " MHz)" << std::endl;
    m_output_file << "Average dynamic power :" << std::setw(15) << std::setprecision(3)
                  << dynamic_total*to_milliwatt/seconds << " mW" << std::endl;
    m_output_file << "Average leakage power :" << std::setw(15) << std::setprecision(3)
                  << static_total*to_milliwatt/seconds << " mW" << std::endl;
    m_output_file << "Average total power   :" << std::setw(15) << std::setprecision(3)
                  << total_energy*to_milliwatt/seconds << " mW" << std::endl;
    // EDP/ED2P in pJ-based units: pJ*s and pJ*s^2.
    m_output_file << "EDP                   :" << std::setw(15) << std::setprecision(6)
                  << total_energy*seconds << " pJ*s" << std::endl;
    m_output_file << "ED2P                  :" << std::setw(15) << std::setprecision(6)
                  << total_energy*seconds*seconds << " pJ*s^2" << std::endl;
    m_output_file << std::endl;
}

// Print out the result of simulation.
void stats_t::print_results(std::ofstream &m_output_file) {
    m_output_file << std::fixed;

    // A1: shared analytical timeline -- critical-path latency and per-level busy
    // ratios. The bottleneck level identifies compute- vs memory-bound execution.
    m_output_file << "============ Layer timeline =============" << std::endl;
    // CE3 contract: "Compute-schedule latency" (computation + fold fill/setup) is the
    // validated official metric (Gemmini RTL / Eyeriss silicon golden gates).
    // "Critical-path latency" is the memory-inclusive analytical timeline
    // (informational; reflects config unit costs, not externally validated).
    // L11/P4-14: state the timing scope BEFORE the latencies, so a partial rollup cannot be
    // read as an end-to-end network latency. Unsupported layers (pooling, activation,
    // normalization) are excluded from accelerator timing entirely, so their time is simply
    // absent from the numbers below -- not estimated, not zero-filled.
    if(network_rollup) {
        const unsigned total = network_timing_layers + excluded_timing_layers;
        m_output_file << "Timing scope          : " << network_timing_layers << " of " << total
                      << " layers";
        if(excluded_timing_layers > 0) {
            m_output_file << "  (PARTIAL: " << excluded_timing_layers
                          << " unsupported layer(s) excluded; this is NOT an end-to-end network"
                          << " latency)";
        } else {
            m_output_file << "  (complete: every layer is supported)";
        }
        m_output_file << std::endl;
    }
    m_output_file << "Compute-schedule latency :" << std::setw(11) << std::setprecision(1)
                                                  << stage_axis_compute << " cycles (validated metric)" << std::endl;
    m_output_file << "Critical-path latency :" << std::setw(11) << std::setprecision(1)
                                               << layer_latency << " cycles" << std::endl;
    if(!graph_residency_note.empty()) {
        m_output_file << "Tensor residency     : " << graph_residency_note << std::endl;
    }
    const char *stage_names[5] = {"DRAM", "Multi-chip (NoP)", "Global buffer", "PE array", "PE (compute+LB)"};
    const double stage_values[5] = {busy_cycle_dram, busy_cycle_multi_chip,
                                    busy_cycle_global_buffer, busy_cycle_pe_array, busy_cycle_pe};
    unsigned bottleneck = 0;
    m_output_file << "Busy cycles (ratio of critical path)"
                  << (network_rollup ? "  [network: summed over layers]" : "") << std::endl;
    for(unsigned s = 0; s < 5; ++s) {
        const double ratio = (layer_latency > 0.0) ? stage_values[s]/layer_latency*100.0 : 0.0;
        m_output_file << " * " << std::left << std::setw(19) << stage_names[s] << std::right
                      << ":" << std::setw(11) << std::setprecision(1) << stage_values[s]
                      << " cycles (" << std::setprecision(1) << ratio << " %)" << std::endl;
        if(stage_values[s] > stage_values[bottleneck]) bottleneck = s;
    }
    m_output_file << "Bottleneck level      :" << std::setw(17)
                                               << stage_names[bottleneck] << std::endl;
    // CE7/L3: at LAYER scope these are the axes each stage busy was computed from
    // (busy = max of its axes; the PE stage additionally includes the compute schedule and
    // the format-IP axis below). At NETWORK scope they are per-axis work totals summed over
    // layers, which do NOT reconstruct network busy -- see stats_t::network_rollup. The
    // header states which contract applies so a reader (or a checker) cannot mix them up.
    m_output_file << "Busy-cycle axes (access / link / overlap)"
                  << (network_rollup ? "  [network: per-axis work summed over layers;"
                                       " NOT reducible to busy]"
                                     : "  [layer: busy = max of these axes]") << std::endl;
    for(unsigned s = 0; s < 5; ++s) {
        m_output_file << " * " << std::left << std::setw(19) << stage_names[s] << std::right
                      << ":" << std::setw(11) << std::setprecision(1) << stage_axis_access[s]
                      << " /" << std::setw(11) << stage_axis_link[s]
                      << " /" << std::setw(11) << stage_axis_overlap[s] << " cycles" << std::endl;
    }
    // P4-1: the PE stage has a fourth axis (format IP) that does not fit the shared
    // access/link/overlap table but does participate in busy_cycle_pe; print it so the PE
    // busy value stays reconstructible from a LAYER report (L3: at network scope it is a
    // work total like the other axes, not a term of the printed busy value).
    m_output_file << "PE format-IP axis     :" << std::setw(11) << std::setprecision(1)
                                               << stage_axis_format << " cycles" << std::endl;
    // LB7/P1-A: which rule combined the PE axes into PE busy. double-buffered = max()
    // over {compute, access, link, overlap, format}; single-buffered =
    // compute + max(access, link, overlap) + format (the overlap axis IS the pipelined
    // makespan of the access/link work, so the axes are not additive among themselves).
    m_output_file << "PE local buffer       :" << std::setw(18)
                                               << (pe_double_buffer ? "double-buffered" : "single-buffered")
                                               << std::endl;
    // L5: where the per-tile timeline stalled. A stage's stall is time it was ready and
    // its input had arrived, but the buffer in front of it was full -- i.e. downstream
    // back-pressure, the term an average-per-tile closed form cannot produce. Zero across
    // the board means every boundary was deep enough for the rates on either side.
    // L5/L6: the staging depth each boundary ran at -- 1 tile in flight means the two
    // sides alternate per tile, 2 means they pipeline. This is the overlap contract itself,
    // so it is reported rather than left implicit in a config flag.
    const char *boundary_names[4] = {"DRAM -> Multi-chip", "Multi-chip -> GLB",
                                     "GLB -> PE array", "PE array -> PE"};
    // An SFU-only scope (standalone softmax layer, or a rollup that folded only such
    // layers) never exercised the memory pipeline, so its init-value depths are not a
    // measured contract and must not print as one.
    const bool boundary_contract_measured =
        network_rollup ? network_rollup_mapped_seeded : !sfu_only_layer;
    if(boundary_contract_measured) {
        m_output_file << "Buffer depth (tiles in flight across each boundary)" << std::endl;
        for(unsigned b = 0; b < 4; ++b) {
            m_output_file << " * " << std::left << std::setw(19) << boundary_names[b] << std::right
                          << ":" << std::setw(11) << timeline_boundary_depth[b] << " tiles" << std::endl;
        }
    } else {
        m_output_file << "Buffer depth          : n/a (no mapped pipeline stages in this scope)"
                      << std::endl;
    }
    m_output_file << "Back-pressure stall (blocked by a full downstream buffer)" << std::endl;
    for(unsigned s = 0; s < 5; ++s) {
        m_output_file << " * " << std::left << std::setw(19) << stage_names[s] << std::right
                      << ":" << std::setw(11) << std::setprecision(1) << timeline_stall[s]
                      << " cycles" << std::endl;
    }
    // L1: how the GLB axis combines its three datatype streams. A shared SRAM serializes
    // them on one port (sum); separate per-datatype partitions serve them concurrently
    // (max). Printed because it is the rule that makes the GLB access axis checkable
    // against the per-datatype access cycles reported further down: the axis can never
    // exceed that combination, and on a single-chip config it equals it. The pre-fix code
    // added the type-combined read side to the type-combined fill side, which on a
    // separate buffer summed the peaks of two DIFFERENT partitions and broke this bound.
    m_output_file << "GLB datatype rule     :" << std::setw(18)
                                               << (global_buffer_type == memory_type_t::SHARED
                                                   ? "sum (shared port)" : "max (partitions)")
                                               << std::endl;
    m_output_file << std::endl;

    // SFU (plan/plan_sfu.md): explicit activation-cost block. Printed for runs that model
    // an SFU, and -- as a scope statement only -- for runs where a nonlinear activation
    // executed without one. Absent for legacy linear-only runs, keeping their reports
    // unchanged.
    if(sfu_present || sfu_active || activation_unmodeled) {
        m_output_file << "============ SFU (activation) ===========" << std::endl;
        if(sfu_present || sfu_active) {
            m_output_file << "Operation             : " << (sfu_operation.empty() ? "none" : sfu_operation)
                          << std::endl;
            const bool sfu_streaming_contract = sfu_queue_depth >= 2;
            m_output_file << "Physical units x lanes:" << std::setw(11) << sfu_physical_units
                          << " x " << sfu_lanes << "  (post-accumulator, fused)" << std::endl;
            m_output_file << "Queue depth           :" << std::setw(11)
                          << std::max(1u, sfu_queue_depth)
                          << (sfu_streaming_contract
                              ? "  (STREAMING: bounded output-tile queue; the SFU overlaps"
                                " the producer)"
                              : "  (serial contract: the SFU drains after the producer)")
                          << std::endl;
            m_output_file << "Valid output elements :" << std::setw(18) << sfu_valid_elements
                          << "  (network dims clamped by mapping coverage; padding excluded)"
                          << std::endl;
            m_output_file << "Scalar operations     :" << std::setw(18) << sfu_operations << std::endl;
            m_output_file << "Invocations / chunks  :" << std::setw(11) << sfu_invocations
                          << " / " << sfu_chunks << std::endl;
            // Phase-2: the final_output_tile event contract in force for this scope.
            if(!sfu_commit_note.empty()) {
                m_output_file << "Output tile commits   :" << std::setw(11) << sfu_commit_events
                              << "  (" << sfu_commit_note << ")" << std::endl;
            }
            // Phase-5 visibility split: busy = exposed (critical-path extension) + hidden
            // (overlapped with the producer pipeline). Under the serial contract the
            // whole window is exposed by definition.
            m_output_file << "SFU busy cycles       :" << std::setw(11) << std::setprecision(1)
                          << sfu_busy_cycle << " cycles" << std::endl;
            m_output_file << " * on critical path   :" << std::setw(11) << std::setprecision(1)
                          << sfu_serial_cycle << " cycles"
                          << (sfu_streaming_contract ? " (fill/drain + bottleneck exposure)"
                                                     : " (fully exposed by the serial contract)")
                          << std::endl;
            m_output_file << " * hidden by producer :" << std::setw(11) << std::setprecision(1)
                          << sfu_hidden_cycle << " cycles" << std::endl;
            m_output_file << "Producer queue stall  :" << std::setw(11) << std::setprecision(1)
                          << sfu_stall_cycle
                          << (sfu_streaming_contract
                              ? " cycles (PE stage blocked on a full SFU queue)"
                              : " cycles (serial contract: the output drain waits on the SFU)")
                          << std::endl;
            if(sfu_streaming_contract && !network_rollup) {
                m_output_file << "SFU bottleneck        : "
                              << (sfu_stall_cycle > 0.0
                                  ? "YES -- SFU throughput limits the producer"
                                  : "no -- the SFU hides behind the producer except"
                                    " fill/drain")
                              << std::endl;
            }
            m_output_file << "Tail lane utilization :" << std::setw(11) << std::setprecision(2)
                          << sfu_tail_lane_utilization << std::endl;
            m_output_file << "Ingress (elem/bytes/txn):" << std::setw(14) << sfu_ingress_elements
                          << "/" << sfu_ingress_bytes << "/" << sfu_ingress_transactions
                          << "  (accumulator precision; internal only, no extra DRAM traffic)"
                          << std::endl;
            m_output_file << "Egress  (elem/bytes/txn):" << std::setw(14) << sfu_egress_elements
                          << "/" << sfu_egress_bytes << "/" << sfu_egress_transactions
                          << std::endl;
            // Phase-7: standalone-softmax operand streaming between the memory hierarchy
            // and the SFU. The DRAM/GLB costs sit on those components' energy rows; the
            // serial makespans below are part of this layer's critical path.
            if(sfu_stream_active) {
                m_output_file << "Operand streaming     : " << sfu_stream_residency << std::endl;
                m_output_file << " * ingress            :" << std::setw(14)
                              << sfu_stream_ingress_bytes << " bytes / "
                              << std::setprecision(1) << sfu_stream_ingress_cycle
                              << " cycles (serial, before the SFU passes)" << std::endl;
                m_output_file << " * egress             :" << std::setw(14)
                              << sfu_stream_egress_bytes << " bytes / "
                              << std::setprecision(1) << sfu_stream_egress_cycle
                              << " cycles (serial, after the SFU passes)" << std::endl;
            }
            m_output_file << "Timing profile        : "
                          << (sfu_timing_calibrated
                              ? "declared per-operation latency/II"
                              : "DEFAULTED (1/1) for some active operation -- UNCALIBRATED")
                          << std::endl;
            // Phase-8: a declared profile is still only calibration-grade with declared
            // provenance -- otherwise the numbers cannot ground an absolute cycle claim.
            m_output_file << "Timing provenance     : "
                          << (sfu_profile_reference.empty()
                              ? "NOT DECLARED -- the SFU cycle numbers are not"
                                " calibration-grade (plan_sfu.md Phase 8)"
                              : sfu_profile_reference)
                          << std::endl;
            m_output_file << "SFU energy (op/read/write/setup/static)" << std::endl;
            m_output_file << " * " << std::setw(15) << std::setprecision(2) << sfu_op_energy
                          << "/" << sfu_read_energy << "/" << sfu_write_energy
                          << "/" << sfu_setup_energy << "/" << sfu_static_energy
                          << " " << energy_units().label() << std::endl;
            if(!sfu_contract_note.empty()) {
                m_output_file << "Scope                 : " << sfu_contract_note << std::endl;
            }
        }
        if(activation_unmodeled) {
            m_output_file << "Activation scope      : NOT MODELED -- "
                          << activation_unmodeled_note << std::endl;
        }
        m_output_file << std::endl;
    }

    // Weight decompression block (evaluation.md Sec 4). Printed for runs that model a
    // decompression engine; absent otherwise (dense baseline unchanged).
    if(decomp_present) {
        m_output_file << "======= Decompression (weight) ==========" << std::endl;
        m_output_file << "Compression (dense/comp):" << std::setw(11) << std::setprecision(2)
                      << decomp_effective_ratio << "x"
                      << (decomp_bypassed ? "  (BYPASSED: weight did not shrink -> dense)"
                                          : "")
                      << std::endl;
        m_output_file << "Weight bytes dense/comp :" << std::setw(14) << decomp_dense_weight_bytes
                      << " / " << decomp_compressed_weight_bytes << std::endl;
        m_output_file << "Weight DRAM saved       :" << std::setw(11) << std::setprecision(1)
                      << decomp_dram_weight_saved_cycle << " cycles" << std::endl;
        m_output_file << "Decoder tiles / queue   :" << std::setw(11) << decomp_tiles
                      << " / " << decomp_queue_depth
                      << (decomp_overlap ? " (overlap ON)" : " (overlap OFF)") << std::endl;
        if(decomp_tile_ratio_cv > 0.0) {
            m_output_file << "Per-tile CR variation   :" << std::setw(11) << std::setprecision(2)
                          << decomp_tile_ratio_cv << " +-spread, " << decomp_bypassed_tiles
                          << " tile(s) bypassed dense (layer total pinned)" << std::endl;
        }
        m_output_file << "Decoder cycles          :" << std::setw(11) << std::setprecision(1)
                      << decomp_decoder_cycles << " cycles" << std::endl;
        m_output_file << " * on critical path     :" << std::setw(11) << std::setprecision(1)
                      << decomp_exposed_cycle << " cycles" << std::endl;
        m_output_file << " * hidden by pipeline   :" << std::setw(11) << std::setprecision(1)
                      << decomp_hidden_cycle << " cycles" << std::endl;
        m_output_file << "Cycle breakdown (dram/decode/compute/stall)" << std::endl;
        m_output_file << " * " << std::setw(12) << std::setprecision(1) << decomp_breakdown_dram
                      << " /" << std::setw(11) << decomp_breakdown_decode
                      << " /" << std::setw(11) << decomp_breakdown_compute
                      << " /" << std::setw(11) << decomp_breakdown_stall << " cycles" << std::endl;
        m_output_file << "Decoder energy/static   :" << std::setw(11) << std::setprecision(2)
                      << decomp_decoder_energy << " / " << decomp_static_energy << " "
                      << energy_units().label() << std::endl;
        m_output_file << "Timing profile          : "
                      << (decomp_timing_calibrated ? "declared decoder throughput"
                                                   : "UNCALIBRATED (ideal/undeclared throughput)")
                      << std::endl;
        m_output_file << "Provenance              : "
                      << (decomp_profile_reference.empty()
                          ? "NOT DECLARED -- not calibration-grade"
                          : decomp_profile_reference) << std::endl;
        m_output_file << std::endl;
    }

    if(kv_present) {
        if(attn_present) {
            m_output_file << "==== KV cache + attention step ====" << std::endl;
            m_output_file << "Model scope            : attention consumer (QK/softmax/AV, "
                          << attn_algorithm << "), dedicated KV DRAM stream (sequential K/V"
                          << " layout), cache append; one decode step (no multi-step growth)"
                          << std::endl;
        } else {
            m_output_file << "==== KV-cache read (PROXY) ====" << std::endl;
            // The scope statement travels with every result: this is a traffic/supply-
            // scheduling proxy. No attention consumer (QK/softmax/AV), no KV address
            // stream, no KV write or per-token cache state, no cross-layer state; the
            // per-byte DRAM cost is estimated from this layer's own weight-DRAM access.
            m_output_file << "Model scope            : traffic/supply-scheduling proxy --"
                          << " attention execution is NOT modeled" << std::endl;
        }
        m_output_file << "Per-token / context    :" << std::setw(11) << kv_bytes_per_token
                      << " B / " << kv_context_length << " tokens" << std::endl;
        m_output_file << "KV compression         :" << std::setw(11) << std::setprecision(2)
                      << std::fixed << kv_compression_ratio << "x"
                      << (kv_bypassed ? "  (uncompressed)" : "") << std::endl;
        m_output_file << "KV bytes dense/comp    :" << std::setw(11) << kv_dense_read_bytes
                      << " / " << kv_read_bytes << " B" << std::endl;
        m_output_file << "KV schedule            :" << std::setw(11)
                      << (kv_schedule.empty() ? "aggregate" : kv_schedule);
        if(kv_schedule == "streaming" || kv_schedule == "double_buffered") {
            m_output_file << "  (" << kv_tiles << " tiles, buffer " << kv_buffer_tiles << ")";
        }
        m_output_file << std::endl;
        m_output_file << "KV supply DRAM cycles   :" << std::setw(11) << std::setprecision(1)
                      << kv_read_cycles << " cycles ("
                      << (kv_schedule.empty() || kv_schedule == "aggregate"
                          ? "added to DRAM access axis"
                          : "tile-level supply stage vs the layer window")
                      << ")" << std::endl;
        if(attn_present) {
            m_output_file << " * cost basis           : dedicated stream on the DRAM model"
                          << " (device read/write, link, row activations)" << std::endl;
        } else {
            m_output_file << " * cost/byte basis      : this layer's measured weight-DRAM"
                          << " access (estimate; no dedicated KV stream on the DRAM model)"
                          << std::endl;
        }
        if(!attn_present && !kv_schedule.empty() && kv_schedule != "aggregate") {
            m_output_file << " * exposed on crit path :" << std::setw(11) << std::setprecision(1)
                          << kv_exposed_cycle << " cycles" << std::endl;
            m_output_file << " * hidden behind layer  :" << std::setw(11) << std::setprecision(1)
                          << kv_hidden_cycle << " cycles" << std::endl;
            m_output_file << " * buffer-full stall    :" << std::setw(11) << std::setprecision(1)
                          << kv_stall_cycle << " cycles" << std::endl;
            m_output_file << " * consumer proxy       : layer window (attention compute is"
                          << " outside the modeled scope)" << std::endl;
        }
        if(attn_present) {
            m_output_file << "Attention step          :" << std::setw(11) << std::setprecision(1)
                          << attn_step_latency << " cycles appended after the layer ("
                          << kv_schedule << ", " << attn_algorithm << ")" << std::endl;
            m_output_file << " * QK / softmax / AV    :" << std::setw(11) << std::setprecision(1)
                          << attn_qk_cycles << " / " << attn_softmax_cycles << " / "
                          << attn_av_cycles << " cycles"
                          << (attn_compute_calibrated ? ""
                              : "  (UNCALIBRATED: ideal rate for an undeclared unit)")
                          << std::endl;
            m_output_file << " * K / V supply         :" << std::setw(11) << std::setprecision(1)
                          << attn_supply_k_cycle << " / " << attn_supply_v_cycle
                          << " cycles (" << kv_tiles << " tiles/pass, buffer "
                          << kv_buffer_tiles << ")" << std::endl;
            m_output_file << " * KV exposed / hidden  :" << std::setw(11) << std::setprecision(1)
                          << attn_exposed_cycle << " / " << attn_hidden_cycle << " cycles"
                          << std::endl;
            m_output_file << " * buffer-full stall    :" << std::setw(11) << std::setprecision(1)
                          << attn_stall_cycle << " cycles" << std::endl;
            m_output_file << " * cache append write   :" << std::setw(11) << std::setprecision(1)
                          << attn_write_cycle << " cycles" << std::endl;
            m_output_file << " * cache occupancy      :" << std::setw(11)
                          << attn_kv_occupancy_bytes << " B after append";
            if(attn_kv_capacity_bytes > 0) {
                m_output_file << " / " << attn_kv_capacity_bytes << " B capacity";
            }
            m_output_file << std::endl;
            m_output_file << " * attention energy     :" << std::setw(11) << std::setprecision(2)
                          << attn_compute_energy << " " << energy_units().label()
                          << " (QK/AV + softmax)" << std::endl;
        }
        m_output_file << " * compressed read      :" << std::setw(11) << std::setprecision(1)
                      << kv_dram_read_cycles << " cycles" << std::endl;
        m_output_file << " * decoder reconstitute :" << std::setw(11) << std::setprecision(1)
                      << kv_decoder_cycles << " cycles"
                      << (kv_decoder_calibrated ? "" : "  (UNCALIBRATED: ideal decoder)")
                      << std::endl;
        m_output_file << "KV read energy         :" << std::setw(11) << std::setprecision(2)
                      << kv_read_energy << " " << energy_units().label()
                      << (kv_priced ? "" : "  (UNPRICED)") << std::endl;
        m_output_file << "Provenance              : "
                      << (kv_profile_reference.empty()
                          ? "NOT DECLARED -- not calibration-grade"
                          : kv_profile_reference) << std::endl;
        m_output_file << std::endl;
    }

    print_energy_summary(m_output_file);

    /* PE result */
    m_output_file << "============== MAC result ===============" << std::endl;
    m_output_file << "# of computations     :" << std::setw(18) 
                                               << num_computation << std::endl;

    m_output_file << "# of data request to Local buffer" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18)
                                               << num_request_pe[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18)
                                               << num_request_pe[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18)
                                               << num_request_pe[data_type_t::OUTPUT] << std::endl;
    m_output_file << std::endl;
    
    m_output_file << "Cycle" <<std::endl;
    m_output_file << "Computation cycle     :" << std::setw(11) << std::setprecision(1) 
                                               << computation_cycle << " cycles" << std::endl;
    m_output_file << "MAX                   :" << std::setw(11) << std::setprecision(1) 
                                               << max_computation_cycle << " cycles" << std::endl;
    m_output_file << "MIN                   :" << std::setw(11) << std::setprecision(1) 
                                               << min_computation_cycle << " cycles" << std::endl; 
    m_output_file << "AVG                   :" << std::setw(11) << std::setprecision(1)
                                               << avg_computation_cycle << " cycles" << std::endl;

    m_output_file << "Access cycle" << std::endl;
    m_output_file << " * Input register     :" << std::setw(11) << std::setprecision(1) 
                                               << access_cycle_mac[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight register    :" << std::setw(11) << std::setprecision(1) 
                                               << access_cycle_mac[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output register    :" << std::setw(11) << std::setprecision(1) 
                                               << access_cycle_mac[data_type_t::OUTPUT] << " cycles" << std::endl;

    m_output_file << "MAX access cycle" << std::endl;
    m_output_file << " * Input register     :" << std::setw(11) << std::setprecision(1) 
                                               << max_access_cycle_mac[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight register    :" << std::setw(11) << std::setprecision(1) 
                                               << max_access_cycle_mac[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output register    :" << std::setw(11) << std::setprecision(1) 
                                               << max_access_cycle_mac[data_type_t::OUTPUT] << " cycles" << std::endl;

    m_output_file << "MIN access cycle" << std::endl;
    m_output_file << " * Input register     :" << std::setw(11) << std::setprecision(1) 
                                               << min_access_cycle_mac[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight register    :" << std::setw(11) << std::setprecision(1) 
                                               << min_access_cycle_mac[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output register    :" << std::setw(11) << std::setprecision(1) 
                                               << min_access_cycle_mac[data_type_t::OUTPUT] << " cycles" << std::endl;

    m_output_file << "AVG access cycle" << std::endl;
    m_output_file << " * Input register     :" << std::setw(11) << std::setprecision(1) 
                                               << avg_access_cycle_mac[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight register    :" << std::setw(11) << std::setprecision(1) 
                                               << avg_access_cycle_mac[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output register    :" << std::setw(11) << std::setprecision(1) 
                                               << avg_access_cycle_mac[data_type_t::OUTPUT] << " cycles" << std::endl;
    m_output_file << std::endl;

    m_output_file << "Energy" << std::endl;
    m_output_file << "Computation energy    :" << std::setw(15) << std::setprecision(2) 
                                               << computation_energy << " " << energy_units().label() << std::endl;
    m_output_file << "Format-IP cycle (payload/metadata/spill)" << std::endl;
    m_output_file << " * Input               :" << std::setw(11) << std::setprecision(1)
                  << format_cycle_pe[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight              :" << std::setw(11) << std::setprecision(1)
                  << format_cycle_pe[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output              :" << std::setw(11) << std::setprecision(1)
                  << format_cycle_pe[data_type_t::OUTPUT] << " cycles" << std::endl;
    // E4: the reduction tree's own energy, and the two format event counts on their two
    // different precisions. A merged "format energy" number could not show that an fp32
    // accumulator spill moves 4x what an int8 output cast moves.
    // E5: the fold/setup control energy and the events it was charged from, so the identity
    // (events x unit cost) closes from the report alone.
    m_output_file << "Weight-fold energy    :" << std::setw(15) << std::setprecision(2)
                                               << weight_fold_energy << " "
                                               << energy_units().label() << " over "
                                               << std::setprecision(0) << weight_fold_events
                                               << " fold event(s)" << std::endl;
    m_output_file << "Layer-setup energy    :" << std::setw(15) << std::setprecision(2)
                                               << layer_setup_energy << " "
                                               << energy_units().label() << " over "
                                               << std::setprecision(0) << layer_setup_events
                                               << " setup event(s)";
    // RE3: separate "no setup in this schedule" from "a setup runs but its unit cost is not
    // priced". Both used to print as zero events.
    // 2026-08-26: and separate "not priced" from a DECLARED zero -- the E20-2 distinction.
    // This annotation used to test the VALUE (== 0.0), so `layer_setup_energy = 0` -- an
    // explicit modeled-zero statement -- still printed as "not declared" while the
    // unpriced-event table (which consults the schema) correctly did not list it.
    if(layer_setup_events > 0.0 && layer_setup_energy == 0.0) {
        if(energy_cost_schema().is_declared("pe_array", "layer_setup_energy")) {
            m_output_file << " [modeled zero: layer_setup_energy is declared 0 over "
                          << std::setprecision(0) << layer_setup_cycle_pe_array
                          << " setup cycle(s)]";
        } else {
            m_output_file << " [unit cost UNCALIBRATED: the setup executes for "
                          << std::setprecision(0) << layer_setup_cycle_pe_array
                          << " cycle(s) but no layer_setup_energy is declared]";
        }
    } else if(layer_setup_events == 0.0) {
        m_output_file << " [no layer setup is modeled for this architecture]";
    }
    m_output_file << std::endl;
    // Per-stripe schedule bubble (see pe_array.h): same declared/modeled-zero discipline.
    m_output_file << "Stripe-transition energy:" << std::setw(13) << std::setprecision(2)
                  << stripe_transition_energy << " " << energy_units().label() << " over "
                  << std::setprecision(0) << stripe_transition_events << " transition(s)";
    if(stripe_transition_events > 0.0 && stripe_transition_energy == 0.0) {
        if(energy_cost_schema().is_declared("pe_array", "stripe_transition_energy")) {
            m_output_file << " [modeled zero: stripe_transition_energy is declared 0]";
        } else {
            m_output_file << " [unit cost UNCALIBRATED: transitions cost "
                          << std::setprecision(0) << stripe_transition_cycle_pe_array
                          << " cycle(s) total but no stripe_transition_energy is declared]";
        }
    } else if(stripe_transition_events == 0.0) {
        m_output_file << " [not modeled: no stripe_transition_cycle declared]";
    }
    m_output_file << std::endl;
    m_output_file << "Reduction energy      :" << std::setw(15) << std::setprecision(2)
                                               << reduction_energy << " " << energy_units().label()
                                               << std::endl;
    // RE1: the four accumulator/output events, each from its own boundary. A create is a
    // zero-initialized accumulator and is deliberately free -- shown so that "0 energy" is
    // visibly a modeled decision rather than a missing charge.
    m_output_file << "Accumulator retained events   :" << std::setw(10)
                  << accumulator_retained_events
                  << "   (prior partial sum already in the MAC: no LB->MAC read-back, so no reload)"
                  << std::endl;
    m_output_file << "Accumulator create events     :" << std::setw(10)
                                               << accumulator_create_events << std::endl;
    m_output_file << "Accumulator reload bytes      :" << std::setw(10)
                                               << accumulator_reload_bytes << std::endl;
    m_output_file << "Accumulator spill bytes       :" << std::setw(10)
                                               << accumulator_spill_bytes << std::endl;
    // RE1: the PE-side accumulator energy. It is 0 when edge_accumulation moves the accumulator to
    // the array edge -- the energy then appears on the PE-array row, where the config says the
    // accumulator lives.
    m_output_file << "Accumulator energy (PE):" << std::setw(14) << std::setprecision(2)
                                               << accumulator_energy << " "
                                               << energy_units().label() << std::endl;
    m_output_file << "Format-IP energy" << std::endl;
    m_output_file << " * Input               :" << std::setw(15) << std::setprecision(2)
                  << format_energy_pe[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight              :" << std::setw(15) << std::setprecision(2)
                  << format_energy_pe[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Output              :" << std::setw(15) << std::setprecision(2)
                  << format_energy_pe[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;
    m_output_file << "Access energy" << std::endl;
    m_output_file << " * Input register     :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_mac[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight register    :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_mac[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Output register    :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_mac[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;
    m_output_file << std::endl;

    // MAC utilization over the PE-modeled layer duration.
    m_output_file << "Utilization" << std::endl;
    m_output_file << "MAC utilization (time):" << std::setw(16) << std::setprecision(1)
                                               << utilization_mac*100 << " %" << std::endl;
    m_output_file << " * Busy scalar-MAC cycles     :" << std::setw(11) << std::setprecision(1)
                  << mac_busy_cycle << std::endl;
    m_output_file << " * Available scalar-MAC cycles:" << std::setw(11) << std::setprecision(1)
                  << mac_available_cycle << std::endl;
    m_output_file << std::endl;

    m_output_file << "============== PE result ================" << std::endl;
    m_output_file << "# of data transfer to MAC" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18)
                                               << num_data_transfer_pe[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18)
                                               << num_data_transfer_pe[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18)
                                               << num_data_transfer_pe[data_type_t::OUTPUT] << std::endl;
    m_output_file << std::endl;

    m_output_file << "# of data request to PE array" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18)
                                               << num_request_pe_array[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18)
                                               << num_request_pe_array[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18)        
                                               << num_request_pe_array[data_type_t::OUTPUT] << std::endl;
    m_output_file << "Cycle" << std::endl;
    m_output_file << "Access cycle" << std::endl; 
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << access_cycle_lb[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << access_cycle_lb[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << access_cycle_lb[data_type_t::OUTPUT] << " cycles" << std::endl;

    m_output_file << "MAX access cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << max_access_cycle_lb[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << max_access_cycle_lb[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << max_access_cycle_lb[data_type_t::OUTPUT] << " cycles" << std::endl;

    m_output_file << "MIN access cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << min_access_cycle_lb[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << min_access_cycle_lb[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << min_access_cycle_lb[data_type_t::OUTPUT] << " cycles" << std::endl;
    
    m_output_file << "AVG access cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << avg_access_cycle_lb[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << avg_access_cycle_lb[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << avg_access_cycle_lb[data_type_t::OUTPUT] << " cycles" << std::endl;
    m_output_file << std::endl;

    m_output_file << "Energy" << std::endl;
    m_output_file << "Access energy" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_lb[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_lb[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_lb[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;
    m_output_file << std::endl;

    // Local buffer utilization
    m_output_file << "Local buffer utilization" << std::endl;
    if(local_buffer_type == memory_type_t::SEPARATE) {
        // Separate partitions: Average is the mean per-type partition occupancy.
        total_utilization_local_buffer = (utilization_local_buffer[data_type_t::INPUT] +
                                          utilization_local_buffer[data_type_t::WEIGHT] +
                                          utilization_local_buffer[data_type_t::OUTPUT]) / 3.0;
        m_output_file << " * Input data         :" << std::setw(16) << std::setprecision(1)
                                                   << utilization_local_buffer[data_type_t::INPUT]*100 << " %" << std::endl;
        m_output_file << " * Weight             :" << std::setw(16) << std::setprecision(1)
                                                   << utilization_local_buffer[data_type_t::WEIGHT]*100 << " %" << std::endl;
        m_output_file << " * Output data        :" << std::setw(16) << std::setprecision(1)
                                                   << utilization_local_buffer[data_type_t::OUTPUT]*100 << " %" << std::endl;
        m_output_file << " * Average            :" << std::setw(16) << std::setprecision(1)
                                                   << total_utilization_local_buffer*100 << " %" << std::endl;
    }
    else if(local_buffer_type == memory_type_t::SHARED) {
        // Shared buffer holds all types: total occupancy is the sum of per-type ratios.
        total_utilization_local_buffer = utilization_local_buffer[data_type_t::INPUT] +
                                         utilization_local_buffer[data_type_t::WEIGHT] +
                                         utilization_local_buffer[data_type_t::OUTPUT];
        m_output_file << "buffer utilization    :" << std::setw(16) << std::setprecision(1)
                                                   << total_utilization_local_buffer*100 << " %" << std::endl;
    }
    m_output_file << std::endl;

    m_output_file << "========== MAC - Local buffer ===========" << std::endl;
    m_output_file << "# of data transfer" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18)
                                               << num_data_transfer_pe[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18) 
                                               << num_data_transfer_pe[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18) 
                                               << num_data_transfer_pe[data_type_t::OUTPUT] << std::endl;
    m_output_file << std::endl;
    print_transaction_breakdown(m_output_file, "MAC <-> local-buffer",
                                payload_link_transactions_pe,
                                metadata_link_transactions_pe,
                                storage_link_transactions_pe);
    m_output_file << "Cycle (MAC-Local buffer)" << std::endl;
    m_output_file << "Transfer cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_pe[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_pe[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_pe[data_type_t::OUTPUT] << " cycles" << std::endl;

    m_output_file << "Total cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1)
                                               << cycle_mac_lb[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1)
                                               << cycle_mac_lb[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1)
                                               << cycle_mac_lb[data_type_t::OUTPUT] << " cycles" <<  std::endl;
    m_output_file << std::endl;

    m_output_file << "Energy" << std::endl;
    m_output_file << "Transfer energy" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2) 
                                               << transfer_energy_pe[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2) 
                                               << transfer_energy_pe[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << transfer_energy_pe[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;
    m_output_file << std::endl;

    m_output_file << "Static energy (leakage over layer elapsed cycles)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_pe[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_pe[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_pe[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;
    m_output_file << std::endl;

    m_output_file << "============ PE array result ============" << std::endl;
    m_output_file << "# of data transfer to PEs" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18) 
                                               << num_data_transfer_pe_array[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18) 
                                               << num_data_transfer_pe_array[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18) 
                                               << num_data_transfer_pe_array[data_type_t::OUTPUT] << std::endl;
    m_output_file << std::endl;
    print_transaction_breakdown(m_output_file, "PE-array NoC",
                                payload_link_transactions_pe_array,
                                metadata_link_transactions_pe_array,
                                storage_link_transactions_pe_array);
    // RE6 condition 4: state the link contract next to the traffic it prices.
    m_output_file << "NoC link contract     : " << array_noc_link_contract << std::endl;
    m_output_file << "Output tile residency : "
                  << (output_tile_array_resident
                      ? "ARRAY-RESIDENT -- the array's output tile equals the GLB's, so"
                        " intermediate write-backs are skipped until the drain"
                      : "written back to the GLB during the layer")
                  << std::endl;
    m_output_file << std::endl;

    m_output_file << "# of request to Global buffer" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18) 
                                               << num_request_global_buffer[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18) 
                                               << num_request_global_buffer[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18) 
                                               << num_request_global_buffer[data_type_t::OUTPUT] << std::endl;
    m_output_file << std::endl;

    m_output_file << "Cycle" << std::endl;
    m_output_file << "Interconnection cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_pe_array[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_pe_array[data_type_t::WEIGHT] << " cycles" <<  std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_pe_array[data_type_t::OUTPUT] << " cycles" <<  std::endl;
    if(fold_fill_cycle_pe_array > 0.0) {
        // V2: RTL-calibrated weight-residency fold fill (+ per-layer setup) that serializes
        // with computation on the array compute schedule.
        m_output_file << "Fold fill cycle       :" << std::setw(11) << std::setprecision(1)
                                                   << fold_fill_cycle_pe_array << " cycles" << std::endl;
    }
    m_output_file << "Overlapped cycle (temporal buffer)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1)
                                               << cycle_temporal_pe_array[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1)
                                               << cycle_temporal_pe_array[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1)
                                               << cycle_temporal_pe_array[data_type_t::OUTPUT] << " cycles" << std::endl;
    m_output_file << std::endl;

    m_output_file << "Energy" << std::endl;
    // E1: the PE-array temporal buffer's own access energy was accumulated and summed into
    // the network totals but never printed, so a reader could not reconstruct this
    // component's energy -- and any config that gave the temporal buffer a non-zero
    // read/write energy saw no effect in the results at all.
    m_output_file << "Access energy (temporal buffer)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_pe_array[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_pe_array[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_pe_array[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;
    m_output_file << std::endl;

    m_output_file << "Interconnection energy" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2) 
                                               << transfer_energy_pe_array[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2) 
                                               << transfer_energy_pe_array[data_type_t::WEIGHT] << " " << energy_units().label() <<  std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2) 
                                               << transfer_energy_pe_array[data_type_t::OUTPUT] << " " << energy_units().label() <<  std::endl;
    m_output_file << std::endl;

    // RE1: accumulator energy handed over by the PEs when edge_accumulation puts the accumulator at
    // the array edge. Reported here because this is the component that owns it.
    m_output_file << "Accumulator energy (edge):" << std::setw(12) << std::setprecision(2)
                                               << pe_array_accumulator_energy << " "
                                               << energy_units().label() << std::endl;
    m_output_file << "Weight-fold energy    :" << std::setw(15) << std::setprecision(2)
                                               << weight_fold_energy << " "
                                               << energy_units().label() << std::endl;

    m_output_file << "Static energy (leakage over layer elapsed cycles)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_pe_array[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_pe_array[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_pe_array[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;
    m_output_file << std::endl;

    // PE utilization
    m_output_file << "Utilization" << std::endl;
    m_output_file << "PE array utilization  :" << std::setw(16) << std::setprecision(1)
                                               << utilization_pe_array*100 << " %" << std::endl;
    m_output_file << "Temporal-buffer occupancy" << std::endl;
    m_output_file << " * Input data         :" << std::setw(16) << std::setprecision(1)
                                               << utilization_pe_array_buffer[data_type_t::INPUT]*100 << " %" << std::endl;
    m_output_file << " * Weight             :" << std::setw(16) << std::setprecision(1)
                                               << utilization_pe_array_buffer[data_type_t::WEIGHT]*100 << " %" << std::endl;
    m_output_file << " * Output data        :" << std::setw(16) << std::setprecision(1)
                                               << utilization_pe_array_buffer[data_type_t::OUTPUT]*100 << " %" << std::endl;
    m_output_file << std::endl;

    m_output_file << "========= Global buffer result ==========" << std::endl;
    m_output_file << " # of data transfer to PE array" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18) 
                                               << num_data_transfer_global_buffer[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18) 
                                               << num_data_transfer_global_buffer[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18) 
                                               << num_data_transfer_global_buffer[data_type_t::OUTPUT] << std::endl;
    m_output_file << std::endl;
    // E20-3b: the psum round trip, stated as a pair. The OUTPUT line above is the LOAD leg
    // (GLB -> PE array); this is the STORE leg (PE array -> GLB). Every psum the array writes
    // out must come back to be accumulated, so with a reduction split above the array the two
    // legs track each other. They used to disagree by the reduction depth (Eyeriss conv3: 19656
    // stores against 312 loads) because the load side latched a skip flag that never reset
    // within a layer -- a half-counted round trip, which understated the GLB traffic bound
    // without moving the compute-schedule latency. Reported so the pairing is checkable.
    m_output_file << " Psum round trip      :" << std::setw(18)
                  << num_data_transfer_global_buffer[data_type_t::OUTPUT] << " loads / "
                  << psum_writeback_events_global_buffer << " stores"
                  << (psum_writeback_events_global_buffer == 0 ? "  (array-resident: no psum leaves the array)" : "")
                  << std::endl;
    m_output_file << std::endl;

    m_output_file << "PE-array <-> GLB transactions (payload/metadata/serialized)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(8)
                                               << payload_link_transactions_global_buffer[data_type_t::INPUT] << "/"
                                               << metadata_link_transactions_global_buffer[data_type_t::INPUT] << "/"
                                               << storage_link_transactions_global_buffer[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(8)
                                               << payload_link_transactions_global_buffer[data_type_t::WEIGHT] << "/"
                                               << metadata_link_transactions_global_buffer[data_type_t::WEIGHT] << "/"
                                               << storage_link_transactions_global_buffer[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(8)
                                               << payload_link_transactions_global_buffer[data_type_t::OUTPUT] << "/"
                                               << metadata_link_transactions_global_buffer[data_type_t::OUTPUT] << "/"
                                               << storage_link_transactions_global_buffer[data_type_t::OUTPUT] << std::endl;
    // P1-B: which of the streams above bypassed the GLB SRAM (no fill write, no read
    // access) and were delivered straight from the chip ingress to the PE array. Their
    // transactions are physically real fabric traffic but are NOT SRAM accesses, so a
    // GLB-SRAM-traffic comparison must exclude exactly these datatypes.
    m_output_file << "GLB-bypassed (direct stream) :" << std::setw(6)
                  << (global_buffer_bypass[data_type_t::INPUT] ? " input" : "")
                  << (global_buffer_bypass[data_type_t::WEIGHT] ? " weight" : "")
                  << (global_buffer_bypass[data_type_t::OUTPUT] ? " output" : "")
                  << (global_buffer_bypass[data_type_t::INPUT] || global_buffer_bypass[data_type_t::WEIGHT] ||
                      global_buffer_bypass[data_type_t::OUTPUT] ? "" : " none")
                  << std::endl;
    m_output_file << std::endl;

    m_output_file << " # of request to Chip-level processor" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18) 
                                               << num_request_multi_chip[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18) 
                                               << num_request_multi_chip[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18) 
                                               << num_request_multi_chip[data_type_t::OUTPUT] << std::endl;
    m_output_file << std::endl;

    m_output_file << " Cycle" << std::endl;
    m_output_file << " Access cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << access_cycle_global_buffer[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << access_cycle_global_buffer[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << access_cycle_global_buffer[data_type_t::OUTPUT] << " cycles" << std::endl;
    m_output_file << std::endl;
    
    m_output_file << "Energy" << std::endl;
    m_output_file << "Access energy" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2) 
                                               << access_energy_global_buffer[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2) 
                                               << access_energy_global_buffer[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_global_buffer[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;
    m_output_file << "Interconnection energy" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << transfer_energy_global_buffer[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << transfer_energy_global_buffer[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << transfer_energy_global_buffer[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;
    m_output_file << std::endl;

    m_output_file << "Static energy (leakage over layer elapsed cycles)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_global_buffer[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_global_buffer[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_global_buffer[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;
    m_output_file << std::endl;

    // Global buffer utilization
    m_output_file << "Global buffer utilization" << std::endl;
    if(global_buffer_type == memory_type_t::SEPARATE) {
        // Separate partitions: Average is the mean per-type partition occupancy
        // (previously never computed in this branch -> printed as 0).
        total_utilization_global_buffer = (utilization_global_buffer[data_type_t::INPUT] +
                                          utilization_global_buffer[data_type_t::WEIGHT] +
                                          utilization_global_buffer[data_type_t::OUTPUT]) / 3.0;
        m_output_file << " * input data         :" << std::setw(16) << std::setprecision(1)
                                                   << utilization_global_buffer[data_type_t::INPUT]*100 << " %" << std::endl;
        m_output_file << " * Weight             :" << std::setw(16) << std::setprecision(1) 
                                                   << utilization_global_buffer[data_type_t::WEIGHT]*100 << " %" <<  std::endl;
        m_output_file << " * Output data        :" << std::setw(16) << std::setprecision(1) 
                                                   << utilization_global_buffer[data_type_t::OUTPUT]*100 << " %" <<  std::endl;
        m_output_file << " * Average            :" << std::setw(16) << std::setprecision(1) 
                                                   << total_utilization_global_buffer*100 << std::endl;
    }
    else if(global_buffer_type == memory_type_t::SHARED) {
        total_utilization_global_buffer = utilization_global_buffer[data_type_t::INPUT] + utilization_global_buffer[data_type_t::WEIGHT] + utilization_global_buffer[data_type_t::OUTPUT];
        total_utilization_global_buffer *= 100;
        m_output_file << "Buffer utilization    :" << std::setw(16) << std::setprecision(1) 
                                                   << total_utilization_global_buffer << " %" << std::endl;
    }
    m_output_file << std::endl;

    m_output_file << "========== PEs - Global buffer ==========" << std::endl;
    m_output_file << "# of data transfer" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18) 
                                               << num_data_transfer_global_buffer[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18) 
                                               << num_data_transfer_global_buffer[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18) 
                                               << num_data_transfer_global_buffer[data_type_t::OUTPUT] << std::endl;
    m_output_file << std::endl;

    m_output_file << "Cycle (PE array-Global buffer)" << std::endl;
    m_output_file << "Transfer cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_global_buffer[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_global_buffer[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_global_buffer[data_type_t::OUTPUT] << " cycles" << std::endl;
    m_output_file << std::endl;

    m_output_file << "Total cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << cycle_pe_array_global_buffer[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << cycle_pe_array_global_buffer[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << cycle_pe_array_global_buffer[data_type_t::OUTPUT] << " cycles" << std::endl;
    m_output_file << std::endl;
   
    m_output_file << "=========== Multi chip result ===========" << std::endl;
    m_output_file << "# of data transfer to global buffers (loads)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18) 
                                               << num_data_transfer_multi_chip[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18) 
                                               << num_data_transfer_multi_chip[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18) 
                                               << num_data_transfer_multi_chip[data_type_t::OUTPUT] << std::endl;
    m_output_file << std::endl;
    print_transaction_breakdown(m_output_file, "GLB <-> multi-chip NoP",
                                payload_link_transactions_multi_chip,
                                metadata_link_transactions_multi_chip,
                                storage_link_transactions_multi_chip);
    m_output_file << "NoP link contract     : " << nop_link_contract << std::endl;
    m_output_file << std::endl;

    m_output_file << "# of request to off-chip memory" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18) 
                                               << num_request_dram[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18) 
                                               << num_request_dram[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18) 
                                               << num_request_dram[data_type_t::OUTPUT] << std::endl;
    m_output_file << std::endl;

    m_output_file << "Cycle" << std::endl;
    // P2-2: the multi-chip temporal buffer's own access side was aggregated into the
    // stage axis but never reported per datatype, so the NoP source-read SHARING
    // contract -- a broadcast datatype (one distinct chunk, e.g. a weight tile that does
    // not depend on the CHIPS_X/CHIPS_Y-split dimension) pays its source read ONCE for
    // the whole package, while a partitioned datatype pays one per active chip -- was
    // not visible in the report and could not be checked. This line is the buffer's
    // TOTAL access cycle: the DRAM-side fill writes plus the NoP-side source reads that
    // feed the per-chip global buffers.
    m_output_file << "Access cycle (DRAM fill + NoP source reads)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1)
                                               << access_cycle_multi_chip[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1)
                                               << access_cycle_multi_chip[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1)
                                               << access_cycle_multi_chip[data_type_t::OUTPUT] << " cycles" << std::endl;
    m_output_file << "Interconnection cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_multi_chip[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_multi_chip[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_multi_chip[data_type_t::OUTPUT] << " cycles" << std::endl;
    m_output_file << std::endl;

    m_output_file << "Energy" << std::endl;
    // E1: same gap on the multi-chip temporal buffer -- accumulated, summed, never printed.
    m_output_file << "Access energy (temporal buffer)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_multi_chip[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_multi_chip[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_multi_chip[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;
    m_output_file << std::endl;

    m_output_file << "Interconnection energy" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2) 
                                               << transfer_energy_multi_chip[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2) 
                                               << transfer_energy_multi_chip[data_type_t::WEIGHT] << " " << energy_units().label() <<  std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << transfer_energy_multi_chip[data_type_t::OUTPUT] << " " << energy_units().label() <<  std::endl;
    m_output_file << std::endl;

    // RE1: the final output cast/pack is committed here (the off-chip output store fires exactly
    // once per output element), so both the count and its energy are reported on this component.
    m_output_file << "Output cast bytes     :" << std::setw(18) << output_cast_bytes << std::endl;
    m_output_file << "Output cast energy    :" << std::setw(15) << std::setprecision(2)
                                               << output_cast_energy << " "
                                               << energy_units().label() << std::endl;
    m_output_file << "Output cast cycle     :" << std::setw(15) << std::setprecision(1)
                                               << output_cast_cycle << " cycles"
                                               << "  [enters the fabric's busy time as a MAX"
                                               << " against its access/transfer axes, not a serial"
                                               << " addition]" << std::endl;

    m_output_file << "Static energy (leakage over layer elapsed cycles)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_multi_chip[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_multi_chip[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_multi_chip[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;
    m_output_file << std::endl;

    // Multi-chip utilization
    m_output_file << "Utilization" << std::endl;
    m_output_file << "Chip utilization      :" << std::setw(16) << std::setprecision(1)
                                               << utilization_multi_chip*100 << " %" << std::endl;
    m_output_file << "Temporal-buffer occupancy" << std::endl;
    m_output_file << " * Input data         :" << std::setw(16) << std::setprecision(1)
                                               << utilization_multi_chip_buffer[data_type_t::INPUT]*100 << " %" << std::endl;
    m_output_file << " * Weight             :" << std::setw(16) << std::setprecision(1)
                                               << utilization_multi_chip_buffer[data_type_t::WEIGHT]*100 << " %" << std::endl;
    m_output_file << " * Output data        :" << std::setw(16) << std::setprecision(1)
                                               << utilization_multi_chip_buffer[data_type_t::OUTPUT]*100 << " %" << std::endl;
    m_output_file << std::endl;


#ifndef DRAMSIM3
    m_output_file << "============== DRAM result ==============" << std::endl;
    // L8: the analytical DRAM contract and its explicit limits. Stated in the report because
    // the numbers below are only meaningful within that scope: this is an approximation for
    // sequential, conflict-free streams, not a cycle-accurate DRAM. Build with DRAMSIM3 for
    // request-level channel/rank/bank/row timing.
    if(!network_rollup) {
        m_output_file << "Input halo reuse      : ";
        if(input_halo_reuse_applied) {
            // Both directions share one format so T11 can pin the same identity. The clause
            // after the working set names WHY the correction was allowed to run: coalescing
            // down is a retention claim (needs the ring working set resident), raising up to
            // the union is a lower-bound identity (holds for any buffer size and loop order).
            m_output_file << "applied; " << input_halo_replicated_elements << " -> "
                          << input_halo_unique_elements << " input elements, working set "
                          << input_halo_working_set_bytes
                          << (input_halo_unique_elements < input_halo_replicated_elements
                              ? " B fits GLB, DRAM serialized "
                              : " B union lower bound, DRAM serialized ")
                          << input_halo_pre_dram_transactions << " -> "
                          << storage_link_transactions_dram[data_type_t::INPUT];
        } else if(input_halo_overlap) {
            if(input_halo_capacity_sufficient)
                m_output_file << "not applied; coalesced target did not reduce DRAM traffic";
            else
                m_output_file << "not applied; working set " << input_halo_working_set_bytes
                              << " B does not fit the resident GLB allocation";
        } else {
            m_output_file << "not needed; per-repetition streaming already moves exactly the union";
        }
        m_output_file << std::endl;
    }
    m_output_file << "DRAM timing model     : " << dram_timing_model << std::endl;
    // E20-5: one contract, stated. `dram_config` names a DRAMsim3 device file and is read ONLY in a
    // -DDRAMSIM3 build (npusim.sh defaults to DRAMSIM3=0). It never derives an energy or a latency
    // for the analytical path -- those come from the declared [dram] keys and from nowhere else. So
    // a config naming an HBM2 part and a config naming a DDR3 part are charged identically unless
    // their keys differ, and the device name must not be read as provenance for the numbers.
    m_output_file << "DRAM cost provenance  : the [dram] unit costs, verbatim. `dram_config` selects"
                  << " a DRAMsim3 device model in a -DDRAMSIM3 build only and derives no cost here"
                  << std::endl;
    m_output_file << "  not modeled here    : " << dram_timing_limits << std::endl;
    m_output_file << "# of data transfer to multi chip (loads)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18)
                                               << num_data_transfer_dram[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18) 
                                               << num_data_transfer_dram[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18) 
                                               << num_data_transfer_dram[data_type_t::OUTPUT] << std::endl;
    m_output_file << std::endl;
    print_transaction_breakdown(m_output_file, "Multi-chip <-> DRAM",
                                payload_link_transactions_dram,
                                metadata_link_transactions_dram,
                                storage_link_transactions_dram);

    m_output_file << "Cycle" << std::endl;
    m_output_file << "Access cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << access_cycle_dram[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << access_cycle_dram[data_type_t::WEIGHT] << " cycles" <<  std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << access_cycle_dram[data_type_t::OUTPUT] << " cycles" <<  std::endl;
    m_output_file << std::endl;
    
    m_output_file << "Energy" << std::endl;
    m_output_file << "Access energy" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2) 
                                               << access_energy_dram[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2) 
                                               << access_energy_dram[data_type_t::WEIGHT] << " " << energy_units().label() <<  std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2) 
                                               << access_energy_dram[data_type_t::OUTPUT] << " " << energy_units().label() <<  std::endl;
    m_output_file << std::endl;

    // E1: the off-chip link energy AND the open-page row-activation energy both accumulate
    // into transfer_energy_dram, which was never printed. That made the row_miss_energy knob
    // invisible: gemmini_dram_detail.cfg sets it to 20 pJ/activation and the results showed
    // no trace of it.
    m_output_file << "Interconnection energy (link + row activation)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << transfer_energy_dram[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << transfer_energy_dram[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << transfer_energy_dram[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;
    m_output_file << std::endl;

    m_output_file << "======Multi chips - Off-chip memory======" << std::endl;
    m_output_file << "# of data transfer" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18) 
                                               << num_data_transfer_dram[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18) 
                                               << num_data_transfer_dram[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18) 
                                               << num_data_transfer_dram[data_type_t::OUTPUT] << std::endl;
    m_output_file << std::endl;

    m_output_file << "Cycle (Multi chips - Off-chip memory)" << std::endl;
    m_output_file << "Transfer cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_dram[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_dram[data_type_t::WEIGHT] << " cycles" <<  std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_dram[data_type_t::OUTPUT] << " cycles" <<  std::endl;
    m_output_file << std::endl;

    m_output_file << "Total cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << cycle_chip_dram[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << cycle_chip_dram[data_type_t::WEIGHT] << " cycles" <<  std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << cycle_chip_dram[data_type_t::OUTPUT] << " cycles" <<  std::endl;
    m_output_file << std::endl;
#endif

}
