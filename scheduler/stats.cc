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
    tile_size.reserve(component_type_t::NUM_COMPONENT_TYPES);
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

    // static_energy in the PE config is leakage energy in pJ/cycle. Every
    // physical PE leaks for the modeled duration of the layer, rather than
    // once per data-transfer callback.
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
    for(unsigned i = 0; i < m_multi_chip->get_number_of_active_chips(); ++i) {
        for(unsigned j = 0; j < m_pe_array[i]->get_number_of_pes(); ++j) {
            m_pe_array[i]->pes[j]->update_static_energy(pe_elapsed_cycles);
        }
    }
    for(unsigned i = 0; i < m_multi_chip->get_number_of_chips(); ++i) {
        m_global_buffer[i]->update_static_energy(pe_elapsed_cycles);
    }

    unsigned num_active_pe = 0;
    for(unsigned i = 0; i < m_multi_chip->get_number_of_active_chips(); i++) {
        for(unsigned j = 0; j < m_pe_array[i]->get_number_of_active_pes(); j++) {
            /* Update PE stats */
            // Update the number of computation
            num_computation +=  m_pe_array[i]->pes[j]->num_computation;

            // Update computation cycle 
            computation_cycle = std::max(computation_cycle, m_pe_array[i]->pes[j]->computation_cycle);

            max_computation_cycle = std::max(max_computation_cycle, m_pe_array[i]->pes[j]->computation_cycle);
            min_computation_cycle = (j == 0) ? m_pe_array[i]->pes[j]->computation_cycle : std::min(min_computation_cycle, m_pe_array[i]->pes[j]->computation_cycle);
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
                min_access_cycle_mac[k] = (j == 0) ? m_pe_array[i]->pes[j]->access_cycle_mac[k] : std::min(min_access_cycle_mac[k], m_pe_array[i]->pes[j]->access_cycle_mac[k]);
                avg_access_cycle_mac[k] += m_pe_array[i]->pes[j]->access_cycle_mac[k];

                // Update access energy of the computing units
                access_energy_mac[k] += m_pe_array[i]->pes[j]->access_energy_mac[k];
                
                // Update access cycle of the local buffer
                access_cycle_lb[k] = std::max(access_cycle_lb[k], m_pe_array[i]->pes[j]->access_cycle_lb[k]);

                max_access_cycle_lb[k] = std::max(max_access_cycle_lb[k], m_pe_array[i]->pes[j]->access_cycle_lb[k]);
                min_access_cycle_lb[k] = (j == 0) ? m_pe_array[i]->pes[j]->access_cycle_lb[k] : std::min(min_access_cycle_lb[k], m_pe_array[i]->pes[j]->access_cycle_lb[k]);
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
        }
        
        for(unsigned j = 0; j < m_pe_array[i]->get_number_of_pes(); j++) {
            for(unsigned k = 0; k < data_type_t::NUM_DATA_TYPES; k++) {
                // Update static energy of PE
                static_energy_pe[k] += m_pe_array[i]->pes[j]->static_energy[k];
            }
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
    }

    avg_computation_cycle /= num_active_pe;
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        avg_access_cycle_mac[i] /= num_active_pe;
        avg_access_cycle_lb[i] /= num_active_pe;
    }

    mac_available_cycle = static_cast<double>(physical_scalar_macs) * pe_elapsed_cycles;
    utilization_mac = calculate_time_based_mac_utilization(mac_busy_cycle, mac_available_cycle);
    if(utilization_mac > 1.0 + 1e-9) {
        std::cerr << "Error: MAC busy cycles exceed physical MAC capacity" << std::endl;
        exit(1);
    }
    utilization_mac = std::min(1.0, utilization_mac);

    // Update global buffer static energy
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

void stats_t::scale_serial_repetitions(unsigned m_repetitions) {
    if(m_repetitions == 0) {
        std::cerr << "Error: temporal repetition count must be non-zero" << std::endl;
        exit(1);
    }
    if(m_repetitions == 1) return;

    num_computation = scale_counter(num_computation, m_repetitions, "computation count");
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

    scale_counters(&num_request_pe_array, m_repetitions, "PE-array request count");
    scale_counters(&num_data_transfer_pe_array, m_repetitions, "PE-array transfer count");
    scale_counters(&payload_link_transactions_pe_array, m_repetitions, "PE-array payload transactions");
    scale_counters(&metadata_link_transactions_pe_array, m_repetitions, "PE-array metadata transactions");
    scale_counters(&storage_link_transactions_pe_array, m_repetitions, "PE-array storage transactions");
    scale_costs(&access_cycle_pe_array, m_repetitions);
    scale_costs(&access_energy_pe_array, m_repetitions);
    scale_costs(&transfer_cycle_pe_array, m_repetitions);
    scale_costs(&transfer_energy_pe_array, m_repetitions);

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

    scale_counters(&num_request_multi_chip, m_repetitions, "multi-chip request count");
    scale_counters(&num_data_transfer_multi_chip, m_repetitions, "multi-chip transfer count");
    scale_counters(&payload_link_transactions_multi_chip, m_repetitions, "multi-chip payload transactions");
    scale_counters(&metadata_link_transactions_multi_chip, m_repetitions, "multi-chip metadata transactions");
    scale_counters(&storage_link_transactions_multi_chip, m_repetitions, "multi-chip storage transactions");
    scale_costs(&access_cycle_multi_chip, m_repetitions);
    scale_costs(&access_energy_multi_chip, m_repetitions);
    scale_costs(&transfer_cycle_multi_chip, m_repetitions);
    scale_costs(&transfer_energy_multi_chip, m_repetitions);

    scale_counters(&num_request_dram, m_repetitions, "DRAM request count");
    scale_counters(&num_data_transfer_dram, m_repetitions, "DRAM transfer count");
    scale_counters(&payload_link_transactions_dram, m_repetitions, "DRAM payload transactions");
    scale_counters(&metadata_link_transactions_dram, m_repetitions, "DRAM metadata transactions");
    scale_counters(&storage_link_transactions_dram, m_repetitions, "DRAM storage transactions");
    scale_costs(&access_cycle_dram, m_repetitions);
    scale_costs(&access_energy_dram, m_repetitions);
    scale_costs(&cycle_chip_dram, m_repetitions);
    scale_costs(&transfer_cycle_dram, m_repetitions);
    scale_costs(&transfer_energy_dram, m_repetitions);
}

void stats_t::update_network_stats(stats_t *m_source) {
    /* Update PE stats */
    
    // Update the number of computation 
    num_computation += m_source->num_computation;

    // Update computation cost
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
        transfer_energy_pe_array[i] += m_source->transfer_energy_pe_array[i];

        /* Update global buffer stats */

        // Update the number of request to the global buffer
        num_request_global_buffer[i] += m_source->num_request_global_buffer[i];
        
        // Update the number of data transfer from the global buffer
        num_data_transfer_global_buffer[i] += m_source->num_data_transfer_global_buffer[i];

        // Update access cost of the global buffer
        access_cycle_global_buffer[i] += m_source->access_cycle_global_buffer[i];
        access_energy_global_buffer[i] += m_source->access_energy_global_buffer[i];

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
    /* PE result */
    m_output_file << std::fixed;
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

    // PE utilization
    m_output_file << "Utilization" << std::endl;
    m_output_file << "PE array utilization  :" << std::setw(16) << std::setprecision(1) 
                                               << utilization_pe_array*100 << " %" << std::endl;
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
    m_output_file << std::endl;

    // Global buffer utilization
    m_output_file << "Global buffer utilization" << std::endl;
    if(global_buffer_type == memory_type_t::SEPARATE) {
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
    m_output_file << "# of data transfer to global buffers" << std::endl;
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

    // Multi-chip utilization
    m_output_file << "Utilization" << std::endl;
    m_output_file << "Chip utilization      :" << std::setw(16) << std::setprecision(1)
                                               << utilization_multi_chip*100 << " %" << std::endl;
    m_output_file << std::endl;


#ifndef DRAMSIM3
    m_output_file << "============== DRAM result ==============" << std::endl;
    m_output_file << "# of data transfer multi chip" << std::endl;
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
