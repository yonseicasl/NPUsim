#ifndef __DRAM_H__
#define __DRAM_H__

#include "def.h"
#include "dnn_model.h"
#include "scheduler.h"
#include "multi_chip.h"
#include "memory.h"
#include "memory_controller.h"

#include "convolutional.h"

class layer_t;
class multi_chip_t;

class dram_t {

public:
    dram_t(section_config_t m_section_config);
    ~dram_t();

    // Initialize the DRAM.
    void init(section_config_t m_section_config);

    // Connect DRAM to chip-level processor
    void connect(multi_chip_t *m_multi_chip);
    // Connect DRAM to DNN model
    //void connect_layer(nebula::layer_t *m_layer);
    void connect_layer(layer_t *m_layer);
    // Connect DRAM to DNN model
    void disconnect_layer();

    // Update the tile size of DRAM.
    void update_tile_size(scheduler_t *m_scheduler);
    
    void update_offset();

    void check_tile_size();

    /* Get the off-chip memory specifications */

    unsigned get_bitwidth();

    /* Check DRAM status */

    // Check whether the DRAM  is idle or not.
    bool is_idle();

    /* DRAM actions */

    // Transfer data from DRAM to Global buffer.
    void data_transfer(scheduler_t *m_scheduler);

    // Print out the stats of DRAM.
    void print_specification();

    // Reset the stats. 
    void reset();

#ifdef DRAMSIM3
    void print_result();
    void send_request(data_t *m_data, mapping_table_t *m_mapping_table, unsigned m_offset, data_type_t m_data_type, action_type_t m_action_type);
#endif


    // Tile size of DRAM. 
    std::vector<unsigned> tile_size;
    std::vector<size_t> offsets;
    std::vector<bool> skip_transfer;

    /* Off-chip memory stats */
    std::vector<unsigned> num_request;              // The number of request from Global buffer to DRAM
    std::vector<unsigned> num_data_transfer;        // The number of data transfer from DRAM to Global buffer

    std::vector<double> access_cycle;               // Total access cycle of DRAM
    std::vector<double> access_energy;              // Total access energy of DRAM

    std::vector<double> transfer_cycle;             // Total data transfer cycle between Global buffer and DRAM
    std::vector<double> transfer_energy;            // Total data transfer energy between Global buffer and DRAM

    std::vector<double> cycle_chip_dram;            // Overlapped cycle between the off-chip memory and chip-level processor

    /* Off-chip memory unit stats */

    double u_transfer_cycle;                        // Unit transfer cycle between the chip-level processor and the off-chip memory
    double u_transfer_energy;                       // Unit transfer energy between the chip-level processor and the off-chip memory

    std::vector<double> u_read_cycle;               // Unit read cycle to the off-chip memory
    std::vector<double> u_read_energy;              // Unit read energy to the off-chip memory

    std::vector<double> u_write_cycle;              // Unit write cycle to the off-chip memory
    std::vector<double> u_write_energy;             // Unit write energy to the off-chip memory


    /* off chip memory specification */

    std::vector<unsigned> line_size;                // Line size of off-chip memory
    std::vector<unsigned> mask_bits;                // Mast bits
    bool transfer_output;

    //nebula::layer_t *layer;
    layer_t *layer;
    data_t *data;

private:
    size_t size;
    /* DRAM specification */
    float frequency;                                // Frequency of DRAM.
    float bandwidth;                                // Bandwidth between chip-level processor and the off-chip memory
    unsigned bitwidth;                              // The bitwidth between the off-chip memory and chip-level processor

    // Indices to count the sequence of tiled data,
    unsigned input_index;                           // Input index in DRAM.
    unsigned weight_index;                          // Weight index in DRAM.
    unsigned output_index;                          // Output index in DRAM.

    bool done;
    
    multi_chip_t *multi_chip;


#ifdef DRAMSIM3
    std::string dram_config;
    std::string output_dir;
    memory_controller_t *memory;
#endif


};

#endif
