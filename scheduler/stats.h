#ifndef __STATS_H__
#define __STATS_H__


#include <vector>
#include "mapping_table.h"

#include "pe_array.h"

class stats_t {
	
public:
    stats_t();
    ~stats_t();
    // Initialize the stats.
    void init();

    // Print out the stats of accelerator and neural network.
    void print_stats();
    void print_stats(std::ofstream &m_output_file);

    // Update tile size
    void update_tile_size(scheduler_t *m_scheduler);

    // Update layer-wise stats.
    void update_stats(std::vector<pe_array_t*> m_pe_array, std::vector<global_buffer_t*> m_global_buffer, multi_chip_t *m_multi_chip, dram_t *m_dram);

    // Update network stats.
    void update_network_stats(stats_t *m_source);

    // Print out the result of simulation.
    void print_results(std::ofstream &m_output_file);
	
    /* Tile size */
    std::vector<std::vector<unsigned>> tile_size;

    /* PE */
    memory_type_t local_buffer_type;
    size_t num_computation;                                             // Number of computations.
    double computation_cycle;                                           // Total computation cycle.
    double max_computation_cycle;
    double min_computation_cycle;
    double avg_computation_cycle;
    double computation_energy;                                          // Total computation energy.

    std::vector<unsigned> num_request_pe;                               // Number of request to local buffer of PE.
    std::vector<unsigned> num_data_transfer_pe;                         // Number of data transfer to MAC unit of PE.

    std::vector<double> access_cycle_mac;                               // Total access cycle to MAC unit of PE.
    std::vector<double> max_access_cycle_mac;
    std::vector<double> min_access_cycle_mac;
    std::vector<double> avg_access_cycle_mac;
    std::vector<double> access_energy_mac;                              // Total access energy to MAC unit of PE.
    double utilization_mac;                                             // Utilization of computing units

    std::vector<double> access_cycle_lb;                                // Total access cycle to local buffer of PE.
    std::vector<double> max_access_cycle_lb;
    std::vector<double> min_access_cycle_lb;
    std::vector<double> avg_access_cycle_lb;
    std::vector<double> access_energy_lb;                               // Total access energy to local buffer in PE.
    std::vector<double> utilization_local_buffer;                       // Utilization of the local buffer in PE
    double total_utilization_local_buffer;

    std::vector<double> transfer_cycle_pe;                              // Total data transfer cycle between MAC unit and local buffer of PE.
    std::vector<double> transfer_energy_pe;                             // Total data transfer energy between MAC unit and local buffer of PE.
    std::vector<double> cycle_mac_lb;                                   // Overlapped cycle between computing and local buffer

    std::vector<double> static_energy_pe;                               // Static energy of PE

    /* PE array */
    std::vector<unsigned> num_request_pe_array;                         // Number of request to PE array (from PE).
    std::vector<unsigned> num_data_transfer_pe_array;                   // Number of data transfer from PE array (to PE).

    std::vector<double> access_cycle_pe_array;                          // Total access cycle to PE array.
    std::vector<double> access_energy_pe_array;                         // Total access energy to PE array.
    double utilization_pe_array;

    std::vector<double> transfer_cycle_pe_array;                        // Total data transfer cycle between PE and PE array.
    std::vector<double> transfer_energy_pe_array;                       // Total data transfer energy between PE and PE array.

    /* Global buffer */
    memory_type_t global_buffer_type;
    std::vector<unsigned> num_request_global_buffer;                    // Number of request to global buffer (from PE array).
    std::vector<unsigned> num_data_transfer_global_buffer;              // Number of data transfer from global buffer (to PE array).

    std::vector<double> access_cycle_global_buffer;                     // Total access cycle to global buffer.
    std::vector<double> access_energy_global_buffer;                    // Total access energy to global buffer.
    std::vector<double> cycle_pe_array_global_buffer;                   // Overlapped cycle between the global buffer and the PE array
    std::vector<double> utilization_global_buffer;                      // Global buffer utilization
    double total_utilization_global_buffer;

    std::vector<double> transfer_cycle_global_buffer;                   // Total data transfer cycle between global buffer and PE array.
    std::vector<double> transfer_energy_global_buffer;                  // Total data transfer energy between global buffer and PE array.

    std::vector<double> static_energy_global_buffer;                    // Static energy the global buffer

    /* Multi chip */
    std::vector<unsigned> num_request_multi_chip;                       // Number of request to Multi-chip (from global buffer).
    std::vector<unsigned> num_data_transfer_multi_chip;                 // Number of data transfer from Multi-chip (to global buffer).

    std::vector<double> access_cycle_multi_chip;                        // Total access cycle to Multi-chip.
    std::vector<double> access_energy_multi_chip;                       // Total access energy to Multi-chip.
    double utilization_multi_chip;

    std::vector<double> transfer_cycle_multi_chip;                      // Total data transfer cycle between Multi-chip and global buffer.
    std::vector<double> transfer_energy_multi_chip;                     // Total data transfer energy between Multi-chip and global buffer.

    /* DRAM */
    std::vector<unsigned> num_request_dram;                             // Number of request to DRAM (from Multi-chip).
    std::vector<unsigned> num_data_transfer_dram;                       // Number of data transfer from DRAM (to Multi-chip).

    std::vector<double> access_cycle_dram;                              // Total access cycle to DRAM.
    std::vector<double> access_energy_dram;                             // Total access energy to DRAM.
    std::vector<double> cycle_chip_dram;

    std::vector<double> transfer_cycle_dram;                            // Total data transfer cycle between DRAM and Multi-chip.
    std::vector<double> transfer_energy_dram;                           // Total data transfer energy between DRAM and Multi-chip.

};

#endif
