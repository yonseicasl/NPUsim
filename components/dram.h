#ifndef __DRAM_H__
#define __DRAM_H__

#include "def.h"
#include "scheduler.h"
#include "multi_chip.h"
#include "memory_controller.h"

#include "convolutional.h"

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
    void connect_layer(nebula::layer_t *m_layer);
    // Connect DRAM to DNN model
    void disconnect_layer();

    // Update the tile size of DRAM.
    void update_tile_size(scheduler_t *m_scheduler);

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

    // Modeled busy duration of the DRAM for the current layer (max of its cost axes).
    double modeled_elapsed_cycles() const;

    // Reset the stats.
    void reset();

#ifdef DRAMSIM3
    void print_result();
    void send_request(data_t *m_data, mapping_table_t *m_mapping_table, unsigned m_offset, data_type_t m_data_type, action_type_t m_action_type);
#endif


    // Tile size of DRAM. 
    std::vector<unsigned> tile_size;
    std::vector<bool> skip_transfer;

    /* Off-chip memory stats */
    std::vector<unsigned> num_request;              // The number of request from Global buffer to DRAM
    std::vector<unsigned> num_data_transfer;        // The number of data transfer from DRAM to Global buffer

    std::vector<double> access_cycle;               // Total access cycle of DRAM
    std::vector<double> access_energy;              // Total access energy of DRAM

    std::vector<double> transfer_cycle;             // Total data transfer cycle between Global buffer and DRAM
    std::vector<double> transfer_energy;            // Total data transfer energy between Global buffer and DRAM

    std::vector<double> cycle_chip_dram;            // Overlapped cycle between the off-chip memory and chip-level processor

    // Runtime-datatype link transactions. All widths are bits.
    std::vector<size_t> payload_link_transactions;
    std::vector<size_t> metadata_link_transactions;
    std::vector<size_t> storage_link_transactions;

    /* Off-chip memory unit stats */

    // DR2: open-page row-buffer model (sequential-stream approximation). A dense
    // stream of B bytes activates ceil(B/row_buffer_bytes) rows; each activation
    // costs u_row_miss_cycle / u_row_miss_energy. Disabled when row_buffer_bytes == 0.
    size_t row_buffer_bytes;
    double u_row_miss_cycle;
    double u_row_miss_energy;
    // Charge the row activations of one dense stream (used by load and write-back).
    void account_row_activations(data_type_t type, size_t elements);

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

    nebula::layer_t *layer;

private:
    void account_descriptor_dense_load(data_type_t type, size_t elements);

    /* DRAM specification */
    float frequency;                                // Frequency of DRAM.
    float bandwidth;                                // Bandwidth between chip-level processor and the off-chip memory
    unsigned bitwidth;                              // The bitwidth between the off-chip memory and chip-level processor

    // Indices to count the sequence of tiled data,
    unsigned input_index;                           // Input index in DRAM.
    unsigned weight_index;                          // Weight index in DRAM.
    unsigned output_index;                          // Output index in DRAM.

    data_t* input_data;                             // Input neural data
    data_t* weight;                                 // Weight neural data
    data_t* output_data;                            // Output neural data

    bool done;
    
    multi_chip_t *multi_chip;


#ifdef DRAMSIM3
    std::string dram_config;
    std::string output_dir;
    memory_controller_t *memory;
#endif


};

#endif
