#ifndef __GLOBAL_BUFFER_H__
#define __GLOBAL_BUFFER_H__

#include "def.h"
#include <cstddef>
#include "scheduler.h"
#include "pe_array.h"
#include "multi_chip.h"

class pe_array_t;
class multi_chip_t;

class global_buffer_t {

public:

    // Phase-5: the component's clock in MHz. Power needs a single authoritative clock
    // across the modeled components (the timeline is one shared cycle axis), so stats_t
    // compares these and reports power as unsupported when they disagree.
    double clock_mhz() const { return static_cast<double>(frequency); }
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
    // E20-3: whether the physical GLB can hold the ring-buffer working set required to reuse
    // overlapping convolution input windows. Shared buffers also reserve the currently resident
    // weight/output tiles; separate buffers use the input partition only.
    bool can_retain_input_halo(size_t input_working_set_elements) const;
    double get_buffer_size();
    // Get global buffer bitwidth
    unsigned get_bitwidth();


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
    // Account a dense transfer using runtime datatype transactions.
    void account_descriptor_dense_transfer(data_type_t type);
    // Account a GLB->multi-chip OUTPUT write-back using runtime datatype transactions.
    void account_output_writeback_link();
    // Flush the data
    void flush_data(scheduler_t *m_scheduler);

    // Print out the configuration of the Global buffer.
    virtual void print_specification() = 0;
    // Accumulate leakage energy once for the modeled layer duration.
    void update_static_energy(double elapsed_cycles);

    // Modeled busy duration of the global buffer for the current layer (max of its cost axes).
    double modeled_elapsed_cycles() const;

    void reset();

    /* Global buffer specifications */

    multi_chip_t *multi_chip;
    data_t *data;
    bool              double_buffer;                        // Check whether the accelerator uses single buffer or double buffer
    std::vector<bool> bypass;                               // Check if the accelerator uses bypass at the global buffer
    // RE1: the FINAL output cast/pack. This is charged where the completed accumulation is read
    // out of the chip, not inside the PE: a PE writes the same output tile back once per
    // reduction pass, so charging the cast there made its count follow the MAC count (262,144 on
    // GEMM 64x64x64) instead of the final output element count (4,096). With edge_accumulation the
    // GLB output partition IS the accumulator, so this is also its physical owner; without it the
    // output still leaves through here, so the pack happens at this boundary either way.
    size_t output_cast_bytes;
    double output_cast_energy;
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

    std::vector<double> u_static_energy;                    // The unit static energy of the global buffer

    /* Global buffer stats */

    std::vector<unsigned> num_request;                      // The number of data request from the PE array to global buffer
    std::vector<unsigned> num_data_transfer;                // The number of data transfer from the global buffer to PE array
    // Logical payload/metadata and physical link transactions, all in bit-width units.
    std::vector<size_t> payload_link_transactions;
    std::vector<size_t> metadata_link_transactions;
    std::vector<size_t> storage_link_transactions;

    std::vector<double> access_cycle;                       // Total access cycles to the global buffer
    // Fill (multi-chip -> GLB write) side, kept separate so off-chip supply can be
    // repetition-scaled per datatype like the DRAM traffic it mirrors.
    std::vector<double> fill_access_cycle;
    std::vector<double> fill_access_energy;
    std::vector<double> access_energy;                      // Total access energies to the global buffer

    std::vector<double> cycle_pe_array_global_buffer;       // Total cycle between PE array and global buffer

    std::vector<double> static_energy;                      // Static energy of the global buffer

    std::vector<double> transfer_cycle;                     // Total data transfer cycle between PE array and Global buffer
    std::vector<double> transfer_energy;                    // Total data transfer energy between PE array and Global buffer;

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
    std::vector<double> capacity_per_type;                  // Per-type capacity in bytes (partition for separate, total for shared)
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
    // Print out the configuration of the shared Global buffer.
    void print_specification();

private:

};

#endif

