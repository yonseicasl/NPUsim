#ifndef __PE_ARRAY_H__
#define __PE_ARRAY_H__

#include <vector>
#include <cmath>
#include "pe.h"
#include "global_buffer.h"
#include "scheduler.h"

class global_buffer_t;
class pe_t;

class pe_array_t {

public:
    
    pe_array_t(section_config_t m_section_config);
    virtual ~pe_array_t();

    // Initialize the PE array.
    virtual void init(section_config_t m_section_config) = 0;

    // Connect PE array and Global buffer
    void connect(global_buffer_t *m_global_buffer);

    // Update tile size of PE array
    virtual void update_tile_size(scheduler_t *m_scheduler) = 0;

    void update_offset();

    /* Get the PE array's specifications */

    // Get stationary type (same as PE stationary type)
    stationary_type_t get_stationary_type();
    // Get parameter order.
    std::string get_parameter_order();
    // Get the number of PEs in the PE array.
    unsigned get_number_of_pes();
    // Get the number of active PEs
    unsigned get_number_of_active_pes();

    /* Check PE array's status */

    // Return true when the PE array has no data.
    bool is_idle();
    // Check whether the request exist 
    bool is_exist_request_at_pe();
    // Check whether the request exist in the PE array.
    bool is_exist_request_at_buffer();

    bool is_exist_data();

    /* PE array's actions */

    // Wait for the data.
    bool wait_data();
    // Exist data in PE array.
    void fill_data();
    // Request data to Global buffer.
    void request_data();
    // Transfer the data each PE.
    virtual void data_transfer(scheduler_t *m_scheduler) = 0;
    // Flush data at temporal buffer in PE array.
    void flush_data();

    // Print out the configuration of PE array.
    virtual void print_specification() = 0;

    // Reset the component and stat
    void reset();

    std::vector<pe_t*> pes;                                 // A group of PEs.
    data_t *data;
    data_t *input_data;
    data_t *weight;
    data_t *output_data;
	data_t *workspace;							            // Temporal workspace for GEMM used in PE array.

    std::vector<bool> skip_transfer;                        // Check whether the data should transferred or not.
    bool              equal_output_tile;
   
    std::vector<unsigned> tile_size;                        // Tile size of PE array
    std::vector<unsigned> offsets;                          // Offset index of PE array
    unsigned duplicated_input;                              //

    /* PE array specifications */

    unsigned index;                                         // Index of PE array

    /* PE array signals */
    // A signal to indicate the existence of data.
    bool exist_temporal_buffer;
    std::vector<bool> exist_data;                           // Data exist in PE array (in temporal buffer).
    std::vector<bool> request_to_global_buffer;             // Request data from PE array to Global buffer
    std::vector<bool> is_waiting_data;


    /* PE array unit costs */
    
    std::vector<double> u_read_cycle;                       // The unit PE array read cycle
    std::vector<double> u_read_energy;                      // The unit PE array read energy

    std::vector<double> u_write_cycle;                      // The unit PE array write cycle
    std::vector<double> u_write_energy;                     // The unit PE array write energy

    
    /* PE array stats */

    std::vector<unsigned> num_request;                      // The number of data request from PE to PE array
    std::vector<unsigned> num_data_transfer;                // The number of data transfer from PE array to local buffer

    std::vector<double> access_cycle;                       // Total access cycles to PE array
    std::vector<double> access_energy;                      // Total access energies to PE array

    std::vector<double> transfer_cycle;                     // Total data transfer cycle between PE and PE array
    std::vector<double> transfer_energy;                    // Total data transfer energy between PE and PE array

    std::vector<double> cycle_temporal_pe;

    double utilization;                                     // Utilization of PE array

    double write_back_cycle;                                //
    double overlapped_transfer_cycle;                       //

    std::vector<unsigned> line_size;
    std::vector<unsigned> mask_bits;
    

protected:

    global_buffer_t *global_buffer;                         // Global buffer to connect

    stationary_type_t stationary_type;                      // Stationary type at the PE array (equal to PE stationary type)
    std::string array_parameter_order;                      // Order of DNN parameters
    noc_type_t noc_type;                                    // NoC type
	data_format_t data_format;                              // Data format (General or GEMM)
    memory_type_t memory_type;

    unsigned height;                                        // Y-dimension of PE array
    unsigned width;                                         // X-dimension of PE array
    
    unsigned num_pes;                                       // The number of PEs in PE array
    unsigned num_active_pe_x;                               // The number of active PEs in PE array
    unsigned num_active_pe_y;                               // The number of active PEs in PE array

    // Buffer size
    unsigned input_size;                                    // Input buffer size
    unsigned weight_size;                                   // Weight buffer size
    unsigned output_size;                                   // Output buffer size

    float frequency;                                        // The frequency.
    float bandwidth;                                        // The bandwidth between the PE array and Global buffer in GB/s.
    unsigned bitwidth;


    bool     initial;

    double noc_cycle;                                       // NoC cycle
    double noc_energy;                                      // NoC energy

};


#endif

