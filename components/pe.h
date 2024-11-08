#ifndef __PE_H__
#define __PE_H__

#include "def.h"
#include "utils.h"
#include "data.h"

#include "scheduler.h"
#include "pe_array.h"

class pe_array_t;
class pe_t {

public: 
    pe_t(section_config_t m_section_config);
    virtual ~pe_t();
    
    // Initialize the PE.
    void init(section_config_t m_section_config);

    // Connect PE to PE array.
    void connect(pe_array_t *m_pe_array);

    // Update tile size of PE.
    void update_tile_size(scheduler_t *m_scheduler);
    // Update the data offset
    void update_offset();
    // Check tile size of PE
    void check_tile_size();

    /* Get PE specifications */

    // Get stationary type of MAC register
    stationary_type_t get_mac_stationary_type();
    // Get stationary type of local buffer.
    stationary_type_t get_local_buffer_stationary_type();
    // Return parameter order.
    std::string get_parameter_order();
    // Return memory type.
    memory_type_t get_memory_type();

    /* Check PE status */

    // Return true when the PE has no data. 
    bool is_idle();
    // Check whether all data exist in the Local buffer.
    bool is_exist_data();
    // Check whether at least one request exist in the local buffer.
    bool is_exist_request();

    double get_static_power();

    /* PE actions */
    // Wait for the data.
    void wait_data();
    // Exist data in PE.
    void fill_data();
    // Request data to PE array.
    void request_data();
    // Execute the MAC operation.
    void data_transfer_to_mac(scheduler_t *m_scheduler);
    // Flush the PE data.
    void flush_data(scheduler_t *m_scheduler);
    // Execute MAC operation
    virtual void computation(scheduler_t *m_scheduler) = 0;

    // MAC operation
    void mac_operation();
    void activation();
    void max_pooling();
    void avg_pooling();

    // Print out the configuration of PE.
    void print_specification();

    // Reset component and stats.
    void reset();

    /* PE specifications */

    // Data at MAC unit.
    data_t *input_data_mac;                                 // Input data in MAC unit
    data_t *weight_mac;                                     // Weight in MAC unit
    data_t *output_data_mac;                                // Output data in MAC unit

    // Data at local buffer.
    data_t *input_data_lb;                                  // input data in local buffer
    data_t *weight_lb;                                      // weight in local buffer
    data_t *output_data_lb;                                 // output data in local buffer
    
    unsigned input_size;                                    // Input local buffer size
    unsigned weight_size;                                   // Weight local buffer size
    unsigned output_size;                                   // Output local buffer size
    
    std::vector<bool> bypass;                               // Check if bypass is applied at the local buffer
    unsigned index;                                         // Index of PE in PE array.

    std::vector<bool> skip_transfer;                        // Check whether skip data transfer from local buffer to MAC unit

    // Data tile size.
    std::vector<unsigned> tile_size_mac;                    // Tile size of MAC unit
    std::vector<unsigned> tile_size_lb;                     // Tile size of local buffer.

    std::vector<unsigned> offsets_mac;                      //
    std::vector<unsigned> offsets_lb;                       //

    unsigned duplicated_input_mac;                          //
    unsigned duplicated_input_lb;                           //

    /* PE signals */
    std::vector<bool> exist_data_mac;                       // Data exist in MAC
    std::vector<bool> request_to_lb;                        // Request data from MAC unit to Local buffer.

    std::vector<bool> exist_data_lb;                        // Data exist in local buffer
    std::vector<bool> request_to_pe_array;                  // Request data from Local buffer to PE array.


    /* Unit costs of PE */

    double u_computation_cycle;                             // The unit MAC cycle
    double u_computation_energy;                            // The unit MAC energy

    double u_transfer_cycle;                                // The unit data transfer cycle between MAC and local buffer
    double u_transfer_energy;                               // The unit data transfer energy between MAC and local buffer

    double u_dynamic_power_mac;                             // MAC dynamic power
    double u_static_power_mac;                              // MAC static power

    std::vector<double> u_read_cycle_mac;                   // The unit MAC read cycle
    std::vector<double> u_read_energy_mac;                  // The unit MAC read energy

    std::vector<double> u_write_cycle_mac;                  // The unit MAC write cycle
    std::vector<double> u_write_energy_mac;                 // The unit MAC write energy

    std::vector<double> u_read_cycle_lb;                    // The unit local buffer read cycle
    std::vector<double> u_read_energy_lb;                   // The unit local buffer read energy

    std::vector<double> u_write_cycle_lb;                   // The unit local buffer write cycle
    std::vector<double> u_write_energy_lb;                  // The unit local buffer write energy

    std::vector<double> u_dynamic_power_lb;                 // Dynamic power of local buffer
    std::vector<double> u_static_power_lb;                  // Static power of local buffer

    /* PE stats */

    unsigned num_computation;                               // The number of computations

    double    computation_cycle;                            // Total computation cycle 
    double    computation_energy;                           // Total computation energy.
    
    std::vector<unsigned> num_request_to_lb;                // The number of data request from MAC unit to local buffer
    std::vector<unsigned> num_data_transfer_to_mac;         // The number of data transfer from local buffer to MAC unit

    std::vector<double> access_cycle_mac;                   // Total access cycles to MAC unit
    std::vector<double> access_energy_mac;                  // Total access energies to MAC unit

    std::vector<double> access_cycle_lb;                    // Total access cycles to Local buffer.
    std::vector<double> access_energy_lb;                   // Total access energies to Local buffer.

    std::vector<double> transfer_cycle;                     // Total transfer cycle between MAC unit and local buffer
    std::vector<double> transfer_energy;                    // Total transfer energy between MAC unit and local buffer

    std::vector<double> cycle_mac_lb;                       // Total cycle between MAC unit and local buffer.

    double utilization_mac;                                 // Utilization of MAC units.
    std::vector<double> utilization_local_buffer;           // Local buffer utilization.

    double write_back_cycle_mac;                            // Access cycle to MAC units when writing back output data.
    double write_back_cycle_lb;                             // Access cycle to local buffers when writing back output data.
    double overlapped_transfer_cycle;                       // Total transfer cycles between MAC unit and local buffer
    
    std::vector<unsigned> line_size_mac;                    // Data block size of MAC unit.
    std::vector<unsigned> line_size_lb;                     // Data block size of local buffer.

    std::vector<unsigned> mask_bits_mac;                    // Mask bits
    std::vector<unsigned> mask_bits_lb;                     // Mask bits

protected:

    pe_array_t *pe_array;

    stationary_type_t stationary_type_mac;                  // Stationary type of MAC register
    stationary_type_t stationary_type_local_buffer;         // stationary type of local buffer
    std::string parameter_order;                            //
    memory_type_t memory_type;                              // Memory type.
    mac_type_t mac_type;                                    // Mac type.

    // Number of MACS.
    unsigned num_macs;                                      // Number of MAC units (accumulator)
    unsigned mac_width;                                     // Number of multiplier per MAC units
    unsigned num_active_macs;                               // Number of active MAC units.
    unsigned active_mac_width;

    float frequency;                                        // The frequency (GHz)
    float bandwidth;                                        // The bandwidth between MAC and local buffer (GB/sec)
    unsigned bitwidth;                                      // The bitwidth between MAC and local buffer

    // Counters for the number of data transfer.
    unsigned input_index;                                   // Input index  in local buffer.
    unsigned weight_index;                                  // Weight index in local buffer.
    unsigned output_index;                                  // Output index in local buffer.

    // Flush counter 
    unsigned input_flush_counter;                           // Input data flush counter.
    unsigned weight_flush_counter;                          // Weight flush counter.
    unsigned output_flush_counter;                          // Output data flush counter.

    // Whether the PE is in idle state or not.
    bool idle;

};

class undefined_stationary_t : public pe_t {
public:
    undefined_stationary_t(section_config_t m_section_config);
    ~undefined_stationary_t();

    void computation(scheduler_t *m_sheduler);
};

class input_stationary_t : public pe_t {

public:
    input_stationary_t(section_config_t m_section_config);
    ~input_stationary_t();

    // Execute MAC operation
    void computation(scheduler_t *m_scheduler);
};

class weight_stationary_t : public pe_t {

public:
    weight_stationary_t(section_config_t m_section_config);
    ~weight_stationary_t();

    // Execute the MAC operation
    void computation(scheduler_t *m_scheduler);
};

class output_stationary_t : public pe_t {

public:
    output_stationary_t(section_config_t m_section_config);
    ~output_stationary_t();
    
    //Execute the MAC operation
    void computation(scheduler_t *m_scheduler);
};

#endif
