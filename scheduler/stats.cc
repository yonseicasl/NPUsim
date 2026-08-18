#include <iomanip>
#include <fstream>
#include "stats.h"
#include <limits>
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
    timeline_boundary_overlap[0] = timeline_boundary_overlap[1] = false;
    timeline_boundary_overlap[2] = timeline_boundary_overlap[3] = true;
    timeline_physical_macs = 0.0;

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
    busy_cycle_multi_chip = m_multi_chip->modeled_elapsed_cycles();
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

            // Account for scalar-MAC cycles that actually execute operations.
            mac_busy_cycle += static_cast<double>(m_pe_array[i]->pes[j]->num_computation) *
                              m_pe_array[i]->pes[j]->u_computation_cycle;

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

                max_access_cycle_lb[k] = std::max(max_access_cycle_lb[k], m_pe_array[i]->pes[j]->access_cycle_lb[k]);
                min_access_cycle_lb[k] = first_active_pe ? m_pe_array[i]->pes[j]->access_cycle_lb[k] : std::min(min_access_cycle_lb[k], m_pe_array[i]->pes[j]->access_cycle_lb[k]);
                avg_access_cycle_lb[k] += m_pe_array[i]->pes[j]->access_cycle_lb[k];

                // Update access energy of the local buffer
                access_energy_lb[k] += m_pe_array[i]->pes[j]->access_energy_lb[k];

                // Update transfer cost between local buffer and computing unit
                transfer_cycle_pe[k] = std::max(transfer_cycle_pe[k], m_pe_array[i]->pes[j]->transfer_cycle[k]);
                transfer_energy_pe[k] += m_pe_array[i]->pes[j]->transfer_energy[k];
                payload_link_transactions_pe[k] += m_pe_array[i]->pes[j]->payload_link_transactions[k];
                metadata_link_transactions_pe[k] += m_pe_array[i]->pes[j]->metadata_link_transactions[k];
                storage_link_transactions_pe[k] += m_pe_array[i]->pes[j]->storage_link_transactions[k];

                format_cycle_pe[k] = std::max(format_cycle_pe[k], m_pe_array[i]->pes[j]->format_cycle[k]);
                format_energy_pe[k] += m_pe_array[i]->pes[j]->format_energy[k];

                // Update overlapped cycle between local buffer and computing unit
                cycle_mac_lb[k] = std::max(cycle_mac_lb[k], m_pe_array[i]->pes[j]->cycle_mac_lb[k]);
                
                // Update local buffer utilization
                utilization_local_buffer[k] = std::max(utilization_local_buffer[k], m_pe_array[i]->pes[j]->utilization_local_buffer[k]);
            }

            local_buffer_type = m_pe_array[i]->pes[j]->get_memory_type();

            num_active_pe++;
            first_active_pe = false;
        }

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
            access_energy_pe_array[j] += m_pe_array[i]->access_energy[j];

            // Update transfer cycle from PE array to local buffers
            transfer_cycle_pe_array[j] = std::max(transfer_cycle_pe_array[j], m_pe_array[i]->transfer_cycle[j]);
            cycle_temporal_pe_array[j] = std::max(cycle_temporal_pe_array[j], m_pe_array[i]->cycle_temporal_pe[j]);
            transfer_energy_pe_array[j] += m_pe_array[i]->transfer_energy[j];

            /* Update global buffer stats */

            // Update global buffer type
            global_buffer_type = m_global_buffer[i]->get_memory_type();

            // Update the number of request to the global buffer
            num_request_global_buffer[j] += m_global_buffer[i]->num_request[j];
            payload_link_transactions_global_buffer[j] += m_global_buffer[i]->payload_link_transactions[j];
            metadata_link_transactions_global_buffer[j] += m_global_buffer[i]->metadata_link_transactions[j];
            storage_link_transactions_global_buffer[j] += m_global_buffer[i]->storage_link_transactions[j];
            
            // Update the number of data transfer from the global buffer
            num_data_transfer_global_buffer[j] += m_global_buffer[i]->num_data_transfer[j];

            // Update access cost to the global buffer
            access_cycle_global_buffer[j] = std::max(access_cycle_global_buffer[j], m_global_buffer[i]->access_cycle[j]);
            access_energy_global_buffer[j] += m_global_buffer[i]->access_energy[j];
            fill_access_cycle_global_buffer[j] = std::max(fill_access_cycle_global_buffer[j], m_global_buffer[i]->fill_access_cycle[j]);
            fill_access_energy_global_buffer[j] += m_global_buffer[i]->fill_access_energy[j];

            // Update transfer cost between the global buffer to the PE array
            transfer_cycle_global_buffer[j] = std::max(transfer_cycle_global_buffer[j], m_global_buffer[i]->transfer_cycle[j]);
            transfer_energy_global_buffer[j] += m_global_buffer[i]->transfer_energy[j];
            
            // Update overlapped cycle between the global buffer and PE array
            cycle_pe_array_global_buffer[j] = std::max(cycle_pe_array_global_buffer[j], m_global_buffer[i]->cycle_pe_array_global_buffer[j]);

            // Update global buffer utilization
            utilization_global_buffer[j] = std::max(utilization_global_buffer[j], m_global_buffer[i]->utilization[j]);
        }
        
        // Update PE array utilization
        utilization_pe_array = std::max(utilization_pe_array, m_pe_array[i]->utilization);
        for(unsigned j = 0; j < data_type_t::NUM_DATA_TYPES; j++) {
            utilization_pe_array_buffer[j] = std::max(utilization_pe_array_buffer[j], m_pe_array[i]->buffer_utilization[j]);
        }
        // V2: fold fill serializes on the compute schedule (max across parallel arrays).
        fold_fill_cycle_pe_array = std::max(fold_fill_cycle_pe_array, m_pe_array[i]->fold_fill_cycle);
        layer_setup_cycle_pe_array = std::max(layer_setup_cycle_pe_array, m_pe_array[i]->u_layer_setup_cycle);
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

} // namespace

void stats_t::scale_serial_repetitions(unsigned m_repetitions,
                                       const std::vector<unsigned> &m_datatype_repetitions) {
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

    scale_counters(&num_request_pe, m_repetitions, "PE request count");
    scale_counters(&num_data_transfer_pe, m_repetitions, "PE transfer count");
    scale_costs(&format_cycle_pe, m_repetitions);
    scale_costs(&format_energy_pe, m_repetitions);
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

    scale_counters(&num_request_pe_array, m_repetitions, "PE-array request count");
    scale_counters(&num_data_transfer_pe_array, m_repetitions, "PE-array transfer count");
    scale_counters(&payload_link_transactions_pe_array, m_repetitions, "PE-array payload transactions");
    scale_counters(&metadata_link_transactions_pe_array, m_repetitions, "PE-array metadata transactions");
    scale_counters(&storage_link_transactions_pe_array, m_repetitions, "PE-array storage transactions");
    scale_costs(&access_cycle_pe_array, m_repetitions);
    scale_costs(&access_energy_pe_array, m_repetitions);
    scale_costs(&transfer_cycle_pe_array, m_repetitions);
    scale_costs(&cycle_temporal_pe_array, m_repetitions);
    scale_costs(&transfer_energy_pe_array, m_repetitions);
    // V2: fold fill repeats with every global-buffer repetition; the per-layer setup
    // (config/flush/DMA prologue) is a one-time cost added after scaling.
    fold_fill_cycle_pe_array = fold_fill_cycle_pe_array*m_repetitions + layer_setup_cycle_pe_array;

    scale_counters(&num_request_global_buffer, m_repetitions, "global-buffer request count");
    scale_counters(&num_data_transfer_global_buffer, m_repetitions, "global-buffer transfer count");
    scale_counters(&payload_link_transactions_global_buffer, m_repetitions, "global-buffer payload transactions");
    scale_counters(&metadata_link_transactions_global_buffer, m_repetitions, "global-buffer metadata transactions");
    scale_counters(&storage_link_transactions_global_buffer, m_repetitions, "global-buffer storage transactions");
    scale_costs(&access_cycle_global_buffer, m_repetitions);
    scale_costs(&access_energy_global_buffer, m_repetitions);
    scale_costs(&cycle_pe_array_global_buffer, m_repetitions);
    scale_costs(&transfer_cycle_global_buffer, m_repetitions);
    scale_costs(&transfer_energy_global_buffer, m_repetitions);
    scale_costs(&static_energy_global_buffer, m_repetitions);

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

        // GLB fill (write) side mirrors the off-chip supply, so it scales with the
        // same per-datatype factor before folding into the GLB access totals.
        fill_access_cycle_global_buffer[i] *= repetitions;
        fill_access_energy_global_buffer[i] *= repetitions;
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
    auto max_types = [](const std::vector<double> &values) {
        double peak = 0.0;
        for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) peak = std::max(peak, values[i]);
        return peak;
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
    stage_axis_access[2] = (global_buffer_type == memory_type_t::SHARED)
        ? sum_types(access_cycle_global_buffer) : max_types(access_cycle_global_buffer);
    stage_axis_link[2] = sum_types(transfer_cycle_global_buffer);
    stage_axis_overlap[2] = sum_types(cycle_pe_array_global_buffer);
    // PE array: one temporal buffer and one distribution fabric.
    stage_axis_access[3] = sum_types(access_cycle_pe_array);
    stage_axis_link[3] = sum_types(transfer_cycle_pe_array);
    stage_axis_overlap[3] = sum_types(cycle_temporal_pe_array);
    // PE: the compute schedule serializes with fold fill and setup (CE2); the local
    // buffer combines per its memory type; MAC<->LB is one per-PE bus.
    stage_axis_compute = computation_cycle + fold_fill_cycle_pe_array;
    stage_axis_access[4] = (local_buffer_type == memory_type_t::SHARED)
        ? sum_types(access_cycle_lb) : max_types(access_cycle_lb);
    stage_axis_link[4] = sum_types(transfer_cycle_pe);
    stage_axis_overlap[4] = sum_types(cycle_mac_lb);

    busy_cycle_dram = std::max(stage_axis_access[0], std::max(stage_axis_link[0], stage_axis_overlap[0]));
    busy_cycle_multi_chip = std::max(stage_axis_access[1], std::max(stage_axis_link[1], stage_axis_overlap[1]));
    busy_cycle_global_buffer = std::max(stage_axis_access[2], std::max(stage_axis_link[2], stage_axis_overlap[2]));
    busy_cycle_pe_array = std::max(stage_axis_access[3], std::max(stage_axis_link[3], stage_axis_overlap[3]));
    busy_cycle_pe = std::max(std::max(stage_axis_compute, stage_axis_access[4]),
                             std::max(stage_axis_link[4], stage_axis_overlap[4]));

    const double stage_busy[5] = {busy_cycle_dram, busy_cycle_multi_chip,
                                  busy_cycle_global_buffer, busy_cycle_pe_array, busy_cycle_pe};
    double final_latency = 0.0;
    double segment = stage_busy[0];
    for(unsigned b = 0; b < 4; ++b) {
        if(timeline_boundary_overlap[b]) {
            segment = std::max(segment, stage_busy[b + 1]);
        } else {
            final_latency += segment;
            segment = stage_busy[b + 1];
        }
    }
    final_latency += segment;

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
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        access_cycle_global_buffer[i] += fill_access_cycle_global_buffer[i];
        access_energy_global_buffer[i] += fill_access_energy_global_buffer[i];
        fill_access_cycle_global_buffer[i] = 0.0;
        fill_access_energy_global_buffer[i] = 0.0;
    }
}

void stats_t::update_network_stats(stats_t *m_source) {
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

    max_computation_cycle += m_source->max_computation_cycle;
    min_computation_cycle += m_source->min_computation_cycle;
    avg_computation_cycle += m_source->avg_computation_cycle;
    computation_energy += m_source->computation_energy;
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
    m_output_file << "Compute-schedule latency :" << std::setw(11) << std::setprecision(1)
                                                  << stage_axis_compute << " cycles (validated metric)" << std::endl;
    m_output_file << "Critical-path latency :" << std::setw(11) << std::setprecision(1)
                                               << layer_latency << " cycles" << std::endl;
    const char *stage_names[5] = {"DRAM", "Multi-chip (NoP)", "Global buffer", "PE array", "PE (compute+LB)"};
    const double stage_values[5] = {busy_cycle_dram, busy_cycle_multi_chip,
                                    busy_cycle_global_buffer, busy_cycle_pe_array, busy_cycle_pe};
    unsigned bottleneck = 0;
    m_output_file << "Busy cycles (ratio of critical path)" << std::endl;
    for(unsigned s = 0; s < 5; ++s) {
        const double ratio = (layer_latency > 0.0) ? stage_values[s]/layer_latency*100.0 : 0.0;
        m_output_file << " * " << std::left << std::setw(19) << stage_names[s] << std::right
                      << ":" << std::setw(11) << std::setprecision(1) << stage_values[s]
                      << " cycles (" << std::setprecision(1) << ratio << " %)" << std::endl;
        if(stage_values[s] > stage_values[bottleneck]) bottleneck = s;
    }
    m_output_file << "Bottleneck level      :" << std::setw(17)
                                               << stage_names[bottleneck] << std::endl;
    // CE7: the axes each stage busy was computed from (busy = max of its axes;
    // the PE stage additionally includes the compute schedule above).
    m_output_file << "Busy-cycle axes (access / link / overlap)" << std::endl;
    for(unsigned s = 0; s < 5; ++s) {
        m_output_file << " * " << std::left << std::setw(19) << stage_names[s] << std::right
                      << ":" << std::setw(11) << std::setprecision(1) << stage_axis_access[s]
                      << " /" << std::setw(11) << stage_axis_link[s]
                      << " /" << std::setw(11) << stage_axis_overlap[s] << " cycles" << std::endl;
    }
    m_output_file << std::endl;

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
                                               << computation_energy << " pJ" << std::endl;
    m_output_file << "Format-IP cycle (payload/metadata/spill)" << std::endl;
    m_output_file << " * Input               :" << std::setw(11) << std::setprecision(1)
                  << format_cycle_pe[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight              :" << std::setw(11) << std::setprecision(1)
                  << format_cycle_pe[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output              :" << std::setw(11) << std::setprecision(1)
                  << format_cycle_pe[data_type_t::OUTPUT] << " cycles" << std::endl;
    m_output_file << "Format-IP energy" << std::endl;
    m_output_file << " * Input               :" << std::setw(15) << std::setprecision(2)
                  << format_energy_pe[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight              :" << std::setw(15) << std::setprecision(2)
                  << format_energy_pe[data_type_t::WEIGHT] << " pJ" << std::endl;
    m_output_file << " * Output              :" << std::setw(15) << std::setprecision(2)
                  << format_energy_pe[data_type_t::OUTPUT] << " pJ" << std::endl;
    m_output_file << "Access energy" << std::endl;
    m_output_file << " * Input register     :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_mac[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight register    :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_mac[data_type_t::WEIGHT] << " pJ" << std::endl;
    m_output_file << " * Output register    :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_mac[data_type_t::OUTPUT] << " pJ" << std::endl;
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
                                               << access_energy_lb[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_lb[data_type_t::WEIGHT] << " pJ" << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_lb[data_type_t::OUTPUT] << " pJ" << std::endl;
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
                                               << transfer_energy_pe[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2) 
                                               << transfer_energy_pe[data_type_t::WEIGHT] << " pJ" << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << transfer_energy_pe[data_type_t::OUTPUT] << " pJ" << std::endl;
    m_output_file << std::endl;

    m_output_file << "Static energy (leakage over layer elapsed cycles)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_pe[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_pe[data_type_t::WEIGHT] << " pJ" << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_pe[data_type_t::OUTPUT] << " pJ" << std::endl;
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
    m_output_file << "Interconnection energy" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2) 
                                               << transfer_energy_pe_array[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2) 
                                               << transfer_energy_pe_array[data_type_t::WEIGHT] << " pJ" <<  std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2) 
                                               << transfer_energy_pe_array[data_type_t::OUTPUT] << " pJ" <<  std::endl;
    m_output_file << std::endl;

    m_output_file << "Static energy (leakage over layer elapsed cycles)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_pe_array[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_pe_array[data_type_t::WEIGHT] << " pJ" << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_pe_array[data_type_t::OUTPUT] << " pJ" << std::endl;
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
                                               << access_energy_global_buffer[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2) 
                                               << access_energy_global_buffer[data_type_t::WEIGHT] << " pJ" << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_global_buffer[data_type_t::OUTPUT] << " pJ" << std::endl;
    m_output_file << "Interconnection energy" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << transfer_energy_global_buffer[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << transfer_energy_global_buffer[data_type_t::WEIGHT] << " pJ" << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << transfer_energy_global_buffer[data_type_t::OUTPUT] << " pJ" << std::endl;
    m_output_file << std::endl;

    m_output_file << "Static energy (leakage over layer elapsed cycles)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_global_buffer[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_global_buffer[data_type_t::WEIGHT] << " pJ" << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_global_buffer[data_type_t::OUTPUT] << " pJ" << std::endl;
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

    m_output_file << "# of request to off-chip memory" << std::endl;
    m_output_file << " * Input data         :" << std::setw(18) 
                                               << num_request_dram[data_type_t::INPUT] << std::endl;
    m_output_file << " * Weight             :" << std::setw(18) 
                                               << num_request_dram[data_type_t::WEIGHT] << std::endl;
    m_output_file << " * Output data        :" << std::setw(18) 
                                               << num_request_dram[data_type_t::OUTPUT] << std::endl;
    m_output_file << std::endl;

    m_output_file << "Cycle" << std::endl;
    m_output_file << "Interconnection cycle" << std::endl;
    m_output_file << " * Input data         :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_multi_chip[data_type_t::INPUT] << " cycles" << std::endl;
    m_output_file << " * Weight             :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_multi_chip[data_type_t::WEIGHT] << " cycles" << std::endl;
    m_output_file << " * Output data        :" << std::setw(11) << std::setprecision(1) 
                                               << transfer_cycle_multi_chip[data_type_t::OUTPUT] << " cycles" << std::endl;
    m_output_file << std::endl;

    m_output_file << "Energy" << std::endl;
    m_output_file << "Interconnection energy" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2) 
                                               << transfer_energy_multi_chip[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2) 
                                               << transfer_energy_multi_chip[data_type_t::WEIGHT] << " pJ" <<  std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << transfer_energy_multi_chip[data_type_t::OUTPUT] << " pJ" <<  std::endl;
    m_output_file << std::endl;

    m_output_file << "Static energy (leakage over layer elapsed cycles)" << std::endl;
    m_output_file << " * Input data         :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_multi_chip[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_multi_chip[data_type_t::WEIGHT] << " pJ" << std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2)
                                               << static_energy_multi_chip[data_type_t::OUTPUT] << " pJ" << std::endl;
    m_output_file << std::endl;

    // Multi-chip utilization
    m_output_file << "Utilization" << std::endl;
    m_output_file << "Chip utilization      :" << std::setw(16) << std::setprecision(1)
                                               << utilization_multi_chip*100 << " %" << std::endl;
    m_output_file << std::endl;


#ifndef DRAMSIM3
    m_output_file << "============== DRAM result ==============" << std::endl;
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
                                               << access_energy_dram[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight             :" << std::setw(15) << std::setprecision(2) 
                                               << access_energy_dram[data_type_t::WEIGHT] << " pJ" <<  std::endl;
    m_output_file << " * Output data        :" << std::setw(15) << std::setprecision(2) 
                                               << access_energy_dram[data_type_t::OUTPUT] << " pJ" <<  std::endl;
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
