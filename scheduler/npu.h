#ifndef __NPU_H__
#define __NPU_H__

#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>

#ifdef Pytorch
    #include <Python.h>
#endif

#include "convolutional.h"
#include "fully_connected.h"
#include "recurrent.h"

#include "adder_tree.h"
#include "spatial_arch.h"
#include "systolic_array.h"
#include "pe_array.h"

#include "global_buffer.h"

#include "multi_chip.h"
#include "dram.h"

#include "mapping_table.h"
#include "scheduler.h"
#include "stats.h"

class network_t;
class scheduler_t;

class npu_t {

public:
    npu_t();
    ~npu_t();

    // Initialize the simulation environment.
    void init(const std::string m_accelerator_config, const std::string m_network_config, const std::string m_mapping_config);
    // Connect components
    void connect();
    // Execute simulation 
    void run(const std::string m_accelerator_config, const std::string m_network_config);

    // Check if the accelerator is idle or not.
    bool is_idle();

    /* Operation at accelerator components. */

    // MAC operation at PEs (Multi threads).
    void execute_thread(unsigned begin, unsigned end);
    // Computation at PEs.
    void execute();
    // Data transfer from PE array to local buffers in PEs.
    void transfer_data_to_pe();
    // Transfer tiled data from global buffer to PE array.
    void transfer_data_to_pe_array();
    // Data transfer from chip to global buffers.
    void transfer_data_to_global_buffer();
    // Transfer tiled data from DRAM to Multi Chip
    void transfer_data_to_multi_chip();

    // Send data request from Multi Chip to DRAM.
    void request_to_dram();
    // Request from the global buffer to Multi Chip
    void request_to_multi_chip();
    // Send data request from PE array to the global buffer.
    void request_to_global_buffer();
    // Send data request from PE to PE array.
    void request_to_pe_array();


    // Print out the Accelerator specification.
    void print_accelerator_specification();
    // Print out DNN configuration.
    void print_network_configuration(unsigned m_index);
    // Print out the stats.
    void print_stats(const std::string m_accelerator_config, const std::string m_network_config, unsigned m_index);

    /* Print result of the simulation */

    // Print out the simulation result.
    void print_layerwise_results(const std::string m_accelerator_config, const std::string m_network_config, unsigned m_index);
    // Print out the simulation result.
    void print_total_result(const std::string m_accelerator_config, const std::string m_network_config);

    // Reset performance counters and stats.
    void reset();
    // Update tile size for executing next layer
    void update_tile_size();

protected:
    unsigned num_processors;                        // The number of on-chip processors.
    data_format_t data_format;                      // Data format : Convolution or GEMM (General Matrix Multiplication).
    compression_type_t compression_type;            // Compression type : Dense, CSR, CSC, SparseMap

    /* Accelerator components */
    std::vector<pe_array_t*> pe_arrays;             // PE array
    std::vector<global_buffer_t*> global_buffers;   // Global buffer
    multi_chip_t *multi_chip;                       // On-chip processors
    dram_t *dram;                                   // DRAM

	nebula::network_t *network;                     // DNN model obtained from the software framework (PyTorch and Nebula)
	std::vector<mapping_table_t*> mapping_tables;	// Mapping tables.

    nebula::layer_t *layer;                         // Neural layers obtained from the software framework (PyTorch and Nebula)

#ifdef Pytorch

#endif
    std::vector<scheduler_t*> schedulers;           // A set of schedulers.
    scheduler_t *scheduler;
    std::vector<stats_t*> layer_stats;
    stats_t *network_stats;


};

#endif
