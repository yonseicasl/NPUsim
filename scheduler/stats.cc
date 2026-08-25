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
    reduction_energy(0.0),
    weight_fold_energy(0.0),
    layer_setup_energy(0.0),
    weight_fold_events(0.0),
    layer_setup_events(0.0),
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
    mac_busy_cycle(0.0),
    mac_available_cycle(0.0),
    utilization_mac(0.0),
    total_utilization_local_buffer(0.0),
    utilization_pe_array(0.0),
    fold_fill_cycle_pe_array(0.0),
    layer_setup_cycle_pe_array(0.0),
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
            chip_link_pe_array_types += m_pe_array[i]->transfer_cycle[j];
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
            chip_link_global_buffer_types += m_global_buffer[i]->transfer_cycle[j];
            transfer_energy_global_buffer[j] += m_global_buffer[i]->transfer_energy[j];

            // Update overlapped cycle between the global buffer and PE array
            cycle_pe_array_global_buffer[j] = std::max(cycle_pe_array_global_buffer[j], m_global_buffer[i]->cycle_pe_array_global_buffer[j]);
            chip_overlap_global_buffer_types += m_global_buffer[i]->cycle_pe_array_global_buffer[j];

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
    if(m_repetitions == 1) {
        // V2: the one-time per-layer setup applies regardless of repetition count.
        fold_fill_cycle_pe_array += layer_setup_cycle_pe_array;
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
    fold_fill_cycle_pe_array = fold_fill_cycle_pe_array*m_repetitions + layer_setup_cycle_pe_array;
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
    if(input_halo_overlap && input_halo_capacity_sufficient && dram_input_link_bits > 0) {
        const unsigned input = data_type_t::INPUT;
        const size_t before = storage_link_transactions_dram[input];
        const size_t target_payload = runtime_datatypes().payload_transactions(
            data_type_t::INPUT, input_halo_unique_elements, dram_input_link_bits);
        const size_t target_metadata = runtime_datatypes().metadata_transactions(
            data_type_t::INPUT, input_halo_unique_elements, dram_input_link_bits);
        const size_t target_storage = runtime_datatypes().storage_transactions(
            data_type_t::INPUT, input_halo_unique_elements, dram_input_link_bits);
        input_halo_pre_dram_transactions = before;
        if(before > 0 && target_storage < before) {
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

void stats_t::update_network_stats(stats_t *m_source) {
    // L3: mark this object as the network rollup so print_results() states the network
    // axis contract instead of the layer one (see stats_t::network_rollup). The flag also
    // tells us whether this is the FIRST layer folded in, which the min-reductions below
    // need (a min seeded from the default would always report the default).
    const bool first_layer = !network_rollup;
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
    // Phase-5: the clock contract is a config property, identical across layers.
    authoritative_frequency_mhz = m_source->authoritative_frequency_mhz;
    single_clock_domain = m_source->single_clock_domain;
    clock_domain_note = m_source->clock_domain_note;
    // E3: the compute-energy basis is a config property, identical across layers.
    compute_energy_basis = m_source->compute_energy_basis;
    compute_energy_precision_calibrated = m_source->compute_energy_precision_calibrated;
    operand_precision = m_source->operand_precision;
    // L8: the DRAM contract is a config property, identical across layers; carry it so the
    // network report states the same scope as its layers.
    dram_timing_model = m_source->dram_timing_model;
    array_noc_link_contract = m_source->array_noc_link_contract;
    output_tile_array_resident = m_source->output_tile_array_resident;
    nop_link_contract = m_source->nop_link_contract;
    dram_timing_limits = m_source->dram_timing_limits;
    // L5: layers run serially, so their back-pressure stalls add like layer_latency does.
    for(unsigned stage = 0; stage < 5; ++stage) timeline_stall[stage] += m_source->timeline_stall[stage];
    // L5/L6: a boundary is only as decoupled as its least-decoupled layer, so the network
    // rollup reports the min depth rather than a sum (depths are a contract, not work).
    for(unsigned b = 0; b < 4; ++b) {
        timeline_boundary_depth[b] = first_layer
            ? m_source->timeline_boundary_depth[b]
            : std::min(timeline_boundary_depth[b], m_source->timeline_boundary_depth[b]);
    }
    for(unsigned s = 0; s < 5; ++s) {
        stage_axis_access[s] += m_source->stage_axis_access[s];
        stage_axis_link[s] += m_source->stage_axis_link[s];
        stage_axis_overlap[s] += m_source->stage_axis_overlap[s];
    }

    max_computation_cycle += m_source->max_computation_cycle;
    min_computation_cycle += m_source->min_computation_cycle;
    avg_computation_cycle += m_source->avg_computation_cycle;
    computation_energy += m_source->computation_energy;
    reduction_energy += m_source->reduction_energy;
    weight_fold_energy += m_source->weight_fold_energy;
    pe_array_accumulator_energy += m_source->pe_array_accumulator_energy;
    layer_setup_energy += m_source->layer_setup_energy;
    weight_fold_events += m_source->weight_fold_events;
    layer_setup_events += m_source->layer_setup_events;
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
    // DRAM: its device access energy and the off-chip link plus row activations.
    const double dram_dynamic = sum_types(access_energy_dram) + sum_types(transfer_energy_dram);

    const double dynamic_total = mac_dynamic + pe_dynamic + pe_array_dynamic +
                                 global_buffer_dynamic + multi_chip_dynamic + dram_dynamic;
    const double static_total = pe_static + pe_array_static + global_buffer_static +
                                multi_chip_static;

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
    const char *rows[6] = {"MAC (compute+format)", "PE (local buffer)", "PE array",
                           "Global buffer", "Multi-chip (NoP)", "DRAM"};
    const double row_dynamic[6] = {mac_dynamic, pe_dynamic, pe_array_dynamic,
                                   global_buffer_dynamic, multi_chip_dynamic, dram_dynamic};
    const double row_static[6] = {0.0, pe_static, pe_array_static, global_buffer_static,
                                  multi_chip_static, 0.0};
    m_output_file << std::left << std::setw(24) << "Component" << std::right
                  << std::setw(17) << "dynamic" << std::setw(17) << "static"
                  << std::setw(17) << "total" << std::endl;
    for(unsigned i = 0; i < 6; ++i) {
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
                      << NUM_ENERGY_COMPONENTS << " have a missing or absent energy cost, so this"
                      << " total is an UNDERCOUNT (see the rows marked above)" << std::endl;
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
    m_output_file << "Buffer depth (tiles in flight across each boundary)" << std::endl;
    for(unsigned b = 0; b < 4; ++b) {
        m_output_file << " * " << std::left << std::setw(19) << boundary_names[b] << std::right
                      << ":" << std::setw(11) << timeline_boundary_depth[b] << " tiles" << std::endl;
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
    if(layer_setup_events > 0.0 && layer_setup_energy == 0.0) {
        m_output_file << " [unit cost UNCALIBRATED: the setup executes for "
                      << std::setprecision(0) << layer_setup_cycle_pe_array
                      << " cycle(s) but no layer_setup_energy is declared]";
    } else if(layer_setup_events == 0.0) {
        m_output_file << " [no layer setup is modeled for this architecture]";
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
            m_output_file << "applied; " << input_halo_replicated_elements << " -> "
                          << input_halo_unique_elements << " input elements, working set "
                          << input_halo_working_set_bytes << " B fits GLB, DRAM serialized "
                          << input_halo_pre_dram_transactions << " -> "
                          << storage_link_transactions_dram[data_type_t::INPUT];
        } else if(input_halo_overlap) {
            if(input_halo_capacity_sufficient)
                m_output_file << "not applied; coalesced target did not reduce DRAM traffic";
            else
                m_output_file << "not applied; working set " << input_halo_working_set_bytes
                              << " B does not fit the resident GLB allocation";
        } else {
            m_output_file << "not needed; no overlapping legacy-GLB P/Q windows";
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
