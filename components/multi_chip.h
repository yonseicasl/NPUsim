#ifndef __MULTI_CHIP_H__
#define __MULTI_CHIP_H__

#include <vector>
#include "global_buffer.h"
#include "dram.h"
#include "scheduler.h"

class global_buffer_t;
class dram_t;

class multi_chip_t {

public:
    
    multi_chip_t(section_config_t m_section_config);
    ~multi_chip_t();

    // Initialize the chip-level processors
    void init(section_config_t m_section_config);

    // Connect chip-level processors and the global buffer
    void connect(std::vector<global_buffer_t*> m_global_buffers);
    // Connect to chip-level processors and the off-chip memory.
    void connect(dram_t *m_dram);

    // Update tile size of chip-level processors.
    void update_tile_size(scheduler_t *m_scheduler);
    // Update offsets of each data of neural network.
    void update_offset();

    /* Get the chip-level processor's specifications */

    // Get stationary type.
    stationary_type_t get_stationary_type();
    // Get parameter order.
    std::string get_parameter_order();
    // Get the number of chips at chip-level processors
    unsigned get_number_of_chips();
    // Get the number of active chips at chip-level processors
    unsigned get_number_of_active_chips();

    /* Check chip-level processor's status */

    // Check whether the chip-level processor is idle or not
    bool is_idle();
    // Check whether the chip-level processor receives request from the global buffer
    bool is_exist_request_at_global_buffer();
    // Check whether the request exist in the PE array.
    bool is_exist_request_at_buffer();
    // Check whether the data exist is the Multi Chip
    bool is_exist_data();

    /* Chip-level processor's actions */

    // Wait for the data.
    bool wait_data();
    // Exist data in PE array.
    void fill_data();
    // Request data to Global buffer.
    void request_data();
    // Transfer the data each PE.
    void data_transfer(scheduler_t *m_scheduler);
    // Flush data at temporal buffer in PE array.
    void flush_data();

    // Print out specification of the chip-level processors.
    void print_specification();

    // Reset stats.
    void reset();

    dram_t *dram;

    data_t *data;                               // Temporal buffer of chip-level processor
    bool   double_buffer;                       //
   
    std::vector<bool> skip_transfer;              
    bool              equal_output_tile;         

    std::vector<unsigned> tile_size;            // Tile size of chip-level processor
    std::vector<unsigned> offsets;              // Offset index of chip-level processor
    unsigned duplicated_input;

    /* Chip-level processor's specifications */

    /* Chip-level processors signals */
    bool exist_temporal_buffer;
    std::vector<bool> exist_data;               // Data exist in chip-level processor (in temporal buffer)
    std::vector<bool> request_to_dram;          // Request data from the chip-level processors to off-chip memory
    std::vector<bool> is_waiting_data;          // Chip-level processor is waiting data or not

    /* Unit stats */
    double u_transfer_cycle;                    // The unit data transfer cycle between global buffer and chip-level processor (unit NoP cycle)
    double u_transfer_energy;                   // The unit data transfer energy between global buffer and chip-level processor (unit NoP energy)

    std::vector<double> u_read_cycle;           // The unit chip-level processor read cycle (if temporal buffer exist)
    std::vector<double> u_read_energy;          // The unit chip-level processor read energy (if temporal buffer exist)

    std::vector<double> u_write_cycle;          // The unit chip-level processor write cycle (if temporal buffer exist)
    std::vector<double> u_write_energy;         // The unit chip-level processor write energy (if temporal buffer exist)

    /* Chip-level processors stats */
    std::vector<unsigned> num_request;          // The number of data request from the global buffer to chip-level processor
    std::vector<unsigned> num_data_transfer;    // The number of data transfer from the chip-level processor to global buffer

    std::vector<double> access_cycle;           // Total access cycles to the chip-level processor (if temporal buffer exist)
    std::vector<double> access_energy;          // Total access energies to the chip-level processor (if temporal buffer exist)

    std::vector<double> transfer_cycle;         // Total data transfer cycle between the global buffer and chip-level processor (NoP)
    std::vector<double> transfer_energy;        // Total data transfer energy between the global buffer and chip-level processor (NoP)

    std::vector<double> cycle_temporal_chips;

    double utilization;                         // Utilization of chip-level processor

    double write_back_cycle;                    // 
    double overlapped_transfer_cycle;           //

    std::vector<unsigned> line_size;            // Line size of temporal buffer
    std::vector<unsigned> mask_bits;            // Mask bits of temporal buffer

protected:

    std::vector<global_buffer_t*> chips;        // Global buffers to connect

    stationary_type_t stationary_type;          // Stationary type at the chip-level processor (equal to global buffer stationary type)
    std::string parameter_order;                // Order of parameters used in data transfer
    noc_type_t nop_type;                        // NoP type (Mesh)
    memory_type_t memory_type;                  // Memory type of global buffer

    unsigned height;                            // Y-dimension of chip-level processors
    unsigned width;                             // X-dimension of chip-level processors

    unsigned num_chips;                         // The number of processors in chip-level processors
    unsigned num_active_chips_x;                // The number of active processors in X dimension
    unsigned num_active_chips_y;                // The number of active processors in Y dimension

    float frequency;                            // Frequency at chip-level processor
    float bandwidth;                            // The bandwidth between the global buffer and chip-level processors (temporal buffer) in GB/s
    unsigned bitwidth;              

    // Buffer size
    unsigned memory_size;                       // Temporal buffer size
    unsigned input_size;                        // Temporal input buffer size
    unsigned weight_size;                       // Temporal weight buffer size
    unsigned output_size;                       // Temporal output buffer size

    bool     initial;                           // 
    double nop_cycle;                           // NoP cycle
    double nop_energy;                          // NoP energy

    double scatter_cycle;                       // Scatter network cycle
    double scatter_energy;                      // Scatter network energy

    double local_cycle;                         // Local network cycle
    double local_energy;                        // Local network energy

    double gather_cycle;                        // Gather network cycle
    double gather_energy;                       // Gather network energy
};

#endif

