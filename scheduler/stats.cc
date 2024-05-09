#include "stats.h"

stats_t::stats_t() :
    local_buffer_type(memory_type_t::UNDEFINED_MEMORY),
    num_computation(0),
    computation_cycle(0.0),
    max_computation_cycle(0.0),
    min_computation_cycle(0.0),
    avg_computation_cycle(0.0),
    computation_energy(0.0),
    utilization_mac(0.0),
    total_utilization_local_buffer(0.0),
    utilization_pe_array(0.0),
    global_buffer_type(memory_type_t::UNDEFINED_MEMORY),
    total_utilization_global_buffer(0.0),
    utilization_multi_chip(0.0) {

    init();

}

stats_t::~stats_t() {
}

// Initialize the stats.
void stats_t::init() {
    
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

    // Initialize static energy at PE
    static_energy_pe.reserve(data_type_t::NUM_DATA_TYPES);
    static_energy_pe.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    
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

    /* Initialize Global buffer stats */

    // Initialize the number of request to the global buffer
    num_request_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    num_request_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the number of data transfer to the PE array
    num_data_transfer_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0);
    
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
    /*
    std::cout << "* Tile size" << std::endl;
    std::cout << "========= MAC =========="<< std::endl;
    std::cout << "Input data  : " << std::setw(10) 
                                  << tile_size_mac[data_type_t::INPUT]  << std::endl;
    std::cout << "Weight      : " << std::setw(10) 
                                  << tile_size_mac[data_type_t::WEIGHT] << std::endl;
    std::cout << "Output data : " << std::setw(10) 
                                  << tile_size_mac[data_type_t::OUTPUT] << std::endl;
    std::cout << "========================"<< std::endl;
    std::cout << std::endl;

    std::cout << "===== Local buffer =====" << std::endl;
    std::cout << "Input data  : " << std::setw(10) 
                                  << tile_size_local_buffer[data_type_t::INPUT]  << std::endl;
    std::cout << "Weight      : " << std::setw(10) 
                                  << tile_size_local_buffer[data_type_t::WEIGHT] << std::endl;
    std::cout << "Output data : " << std::setw(10) 
                                  << tile_size_local_buffer[data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;

    std::cout << "======= PE array =======" << std::endl;
    std::cout << "Input data  : " << std::setw(10) 
                                  << tile_size_pe[data_type_t::INPUT]  << std::endl;
    std::cout << "Weight      : " << std::setw(10) 
                                  << tile_size_pe[data_type_t::WEIGHT] << std::endl;
    std::cout << "Output data : " << std::setw(10) 
                                  << tile_size_pe[data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;

    std::cout << "===== Global buffer ====" << std::endl;
    std::cout << "Input data  : " << std::setw(10) 
                                  << tile_size_global_buffer[data_type_t::INPUT]  << std::endl;
    std::cout << "Weight      : " << std::setw(10) 
                                  << tile_size_global_buffer[data_type_t::WEIGHT] << std::endl;
    std::cout << "Output data : " << std::setw(10) 
                                  << tile_size_global_buffer[data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;

    std::cout << "======= Processor ======" << std::endl;
    std::cout << "Input data  : " << std::setw(10) 
                                  << tile_size_processor[data_type_t::INPUT]  << std::endl;
    std::cout << "Weight      : " << std::setw(10) 
                                  << tile_size_processor[data_type_t::WEIGHT] << std::endl;
    std::cout << "Output data : " << std::setw(10) 
                                  << tile_size_processor[data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;

    std::cout << "========= DRAM =========" << std::endl;
    std::cout << "Input data  : " << std::setw(10) 
                                  << tile_size_dram[data_type_t::INPUT]  << std::endl;
    std::cout << "Weight      : " << std::setw(10) 
                                  << tile_size_dram[data_type_t::WEIGHT] << std::endl;
    std::cout << "Output data : " << std::setw(10) 
                                  << tile_size_dram[data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;
    */

}

// Print out the stats of accelerator and neural network.
void stats_t::print_stats(std::ofstream &m_output_file) {
    /*
    m_output_file << "* Tile size" << std::endl;
    m_output_file << "========= MAC =========="<< std::endl;
    m_output_file << "Input data  : " << std::setw(10) 
                                      << tile_size_mac[data_type_t::INPUT]  << std::endl;
    m_output_file << "Weight      : " << std::setw(10)
                                      << tile_size_mac[data_type_t::WEIGHT] << std::endl;
    m_output_file << "Output data : " << std::setw(10)
                                      << tile_size_mac[data_type_t::OUTPUT] << std::endl;
    m_output_file << "========================"<< std::endl;
    m_output_file << std::endl;

    m_output_file << "===== Local buffer =====" << std::endl;
    m_output_file << "Input data  : " << std::setw(10) 
                                      << tile_size_local_buffer[data_type_t::INPUT]  << std::endl;
    m_output_file << "Weight      : " << std::setw(10) 
                                      << tile_size_local_buffer[data_type_t::WEIGHT] << std::endl;
    m_output_file << "Output data : " << std::setw(10) 
                                      << tile_size_local_buffer[data_type_t::OUTPUT] << std::endl;
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
    */
}

void stats_t::update_tile_size(scheduler_t *m_scheduler) {
    //tile_size = m_scheduler->tile_size;

    //tile_size_mac = m_scheduler->tile_size[component_type_t::MAC];
}

void stats_t::update_stats(std::vector<pe_array_t*> m_pe_array, std::vector<global_buffer_t*> m_global_buffer, multi_chip_t *m_multi_chip, dram_t *m_dram) {

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

            // Update utilization of computing units
            utilization_mac = std::max(utilization_mac, m_pe_array[i]->pes[j]->utilization_mac);

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

    utilization_mac += m_source->utilization_mac;

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

        // Update transfer cost between the local buffer and computing units
        transfer_cycle_pe[i] += m_source->transfer_cycle_pe[i];
        transfer_energy_pe[i] += m_source->transfer_energy_pe[i];

        // Update overlapped cycle between the local buffer and computing units
        cycle_mac_lb[i] += m_source->cycle_mac_lb[i];

        /* Update PE array stats */

        // Update the number of request to the PE array
        num_request_pe_array[i] += m_source->num_request_pe_array[i];

        // Update the number of data transfer from the PE array
        num_data_transfer_pe_array[i] += m_source->num_data_transfer_pe_array[i];

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
    m_output_file << "Access energy" << std::endl;
    m_output_file << " * Input register     :" << std::setw(15) << std::setprecision(2)
                                               << access_energy_mac[data_type_t::INPUT] << " pJ" << std::endl;
    m_output_file << " * Weight register    :" << std::setw(15) << std::setprecision(2) 
                                               << access_energy_mac[data_type_t::WEIGHT] << " pJ" << std::endl;
    m_output_file << " * Output register    :" << std::setw(15) << std::setprecision(2) 
                                               << access_energy_mac[data_type_t::OUTPUT] << " pJ" << std::endl;
    m_output_file << std::endl;                                           
    
    // MAC utilization.
    m_output_file << "Utilization" << std::endl;
    m_output_file << "MAC utilization       :" << std::setw(16) << std::setprecision(1) 
                                               << utilization_mac*100 << " %" << std::endl;
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

