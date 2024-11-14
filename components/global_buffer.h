#ifndef __GLOBAL_BUFFER_H__
#define __GLOBAL_BUFFER_H__

#include "def.h"
#include "scheduler.h"
#include "pe_array.h"
#include "multi_chip.h"

class pe_array_t;
class multi_chip_t;

class global_buffer_t {

public:
    global_buffer_t(section_config_t m_section_config);
    virtual ~global_buffer_t();

    // Initialize the Global buffer.
    virtual void init(section_config_t m_section_config) = 0;

    // Connect Global buffer and PE array
    void connect(pe_array_t *m_pe_array);
    // Connect Global buffer and chip-level processors
    void connect(multi_chip_t *m_multi_chip);

    // Update tile size at the global buffer.
    void update_tile_size(scheduler_t *m_scheduler);
    // Update the data offset.
    virtual void update_offset() = 0;
    // Check the tile size.
    virtual void check_tile_size() = 0;
    
    /* Get the global buffer's specifications */

    // Get stationary type
    stationary_type_t get_stationary_type();
    // Get memory type
    memory_type_t get_memory_type();
    // Get global buffer size
    double get_buffer_size();
    // Get global buffer bitwidth
    unsigned get_bitwidth();

    // Get dynamic power of global buffer
    virtual double get_dynamic_power() = 0;
    // Get static power of global buffer
    virtual double get_static_power() = 0;

    /* Check global buffer's status */

    // Return true when the Global buffer has no data
    bool is_idle();
    // Check whether all data exist in the Global buffer.
    bool is_exist_data();
    // Check whether at least one request exist in the Global buffer.
    bool is_exist_request();

    /* Global buffer's actions */

    // Wait for the data.
    void wait_data();
    // Exist the data.
    void fill_data();
    // Request data to DRAM.
    void request_data();
    // Transfer the data to PE array.
    void data_transfer(scheduler_t *m_scheduler);
    // Flush the data 
    void flush_data(scheduler_t *m_scheduler);

    // Print out the configuration of the Global buffer.
    virtual void print_specification() = 0;

    void reset();

    /* Global buffer specifications */

    multi_chip_t *multi_chip;
    data_t *data;
    bool              double_buffer;                        // Check whether the accelerator uses single buffer or double buffer
    std::vector<bool> bypass;                               // Check if the accelerator uses bypass at the global buffer
    unsigned index;

    std::vector<bool> skip_transfer;                        // Determine whether transfer data to PE array or not.

    std::vector<unsigned> tile_size;                        // The tile size of Global buffer.
    std::vector<unsigned> offsets;                          // Offset of input data in Global buffer
    unsigned duplicated_input;


    /* Global buffer signals */

    std::vector<bool> exist_data;                           // Data exist in the global buffer
    std::vector<bool> request_to_multi_chip;                // Request data from the global buffer to chip-level processor

    /* Unit stats */

    double u_transfer_cycle;                                // The unit data transfer cycle between PE array and Global buffer
    double u_transfer_energy;                               // The unit data transfer energy between PE array and Global buffer

    std::vector<double> u_read_cycle;                       // The unit global buffer read cycle
    std::vector<double> u_read_energy;                      // The unit global buffer read energy

    std::vector<double> u_write_cycle;                      // The unit global buffer write cycle
    std::vector<double> u_write_energy;                     // The unit global buffer write energy

    std::vector<double> u_dynamic_power;                    // Dynamic power of the global buffer
    std::vector<double> u_static_power;                     // Static power of the global buffer

    /* Global buffer stats */

    std::vector<unsigned> num_request;                      // The number of data request from the PE array to global buffer
    std::vector<unsigned> num_data_transfer;                // The number of data transfer from the global buffer to PE array

    std::vector<double> access_cycle;                       // Total access cycles to the global buffer
    std::vector<double> access_energy;                      // Total access energies to the global buffer

    std::vector<double> cycle_pe_array_global_buffer;       // Total cycle between PE array and global buffer

    std::vector<double> transfer_cycle;                     // Total data transfer cycle between PE array and Global buffer
    std::vector<double> transfer_energy;                    // Total data transfer energy between PE array and Global buffer;

    std::vector<double> dynamic_power;                      // Dynamic power of global buffer. 

    std::vector<double> utilization;                        // Utilization of the global buffer

    double write_back_cycle;                                //
    double overlapped_transfer_cycle;                       // 
    

    std::vector<unsigned> line_size;                        // Line size of the global buffer
    std::vector<unsigned> mask_bits;                        // Mask bits of the global buffer
    bool transfer_output;                                   //

protected:

    pe_array_t *pe_array;

    memory_type_t memory_type;                              // Memory type of the global buffer (i.e., Separated or Shared).
    stationary_type_t stationary_type;                      // Stationary type at the global buffer (i.e., input stationary, weight stationary, or output stationary)
    std::string parameter_order;

    double size;                                            // Total size of the global buffer
    float frequency;                                        // The frequency.
    float bandwidth;                                        // The bandwidth between PE array and Global buffer.
    unsigned bitwidth;                                      // The bitwidth between PE array and Global buffer

    unsigned input_index;                                   // Input index in Global buffer.
    unsigned weight_index;                                  // Weight index in Global buffer.
    unsigned output_index;                                  // Output index in Global buffer.

    // Flush counter
    unsigned input_flush_counter;                           // Flush counter of input data.
    unsigned weight_flush_counter;                          // Flush counter of weight.
    unsigned output_flush_counter;                          // Flush counter of output data.

    bool idle;                                              // A signal whether the Global buffer is in idle state or not.
    bool initial;

};

class separate_buffer_t : public global_buffer_t {

public:
    separate_buffer_t(section_config_t m_section_config);
    ~separate_buffer_t();

    // Initialize the separate Global buffer.
    void init(section_config_t m_section_config);
    // Update the offset.
    void update_offset();
    // Check the tile size.
    void check_tile_size();
    // Get dynamic power
    double get_dynamic_power();
    // Get static power
    double get_static_power();
    // Print out the configuration of the separate Global buffer.
    void print_specification();

private:
    double input_size;                // Input Global buffer size.
    double weight_size;               // Weight Global buffer size.
    double output_size;               // Output Global buffer size.
};

class shared_buffer_t : public global_buffer_t {

public:
    shared_buffer_t(section_config_t m_section_config);
    ~shared_buffer_t();

    // Initialize the shared Global buffer.
    void init(section_config_t m_section_config);
    // Update the offset.
    void update_offset();
    // Check the tile size.
    void check_tile_size();
    // Get dynamic power
    double get_dynamic_power();
    // Get static power
    double get_static_power();
    // Print out the configuration of the shared Global buffer.
    void print_specification();

private:

};

#endif

