#ifndef __SCHEDULER_H__
#define __SCHEDULER_H__

#include <list>
#include "layer.h"
#include "convolutional_layer.h"
#include "connected_layer.h"

#include "mapping_table.h"
#include "def.h"

class scheduler_t {

public:
	scheduler_t(mapping_table_t *m_mappig_table, stationary_type_t pe_stationary, std::string pe_parameter_order,
                                                 stationary_type_t pe_array_stationary, std::string pe_array_parameter_order,
                                                 stationary_type_t multi_chip_stationary, std::string multi_chip_parameter_order);
	~scheduler_t();

    // Initialize the scheduler.
	void init(mapping_table_t *m_mapping_table, stationary_type_t pe_stationary, std::string pe_parameter_order, 
                                                stationary_type_t pe_array_stationary, std::string pe_array_parameter_order,
                                                stationary_type_t multi_chip_stationary, std::string multi_chip_parameter_order);
    
    // Print out the network stats.
    void print_stats();

    // Print out the network stats to file
    void print_stats(std::ofstream &m_output_file);

#ifdef FUNCTIONAL
    void transfer_data(data_t *m_dest, data_t *m_source, 
                       unsigned m_dest_offset, unsigned m_source_offset, 
                       component_type_t m_dest_type, component_type_t m_source_type, 
                       data_type_t m_data_type, stationary_type_t m_stationary_type, action_type_t m_action_type);
    /* TODO : update for the next version */
    void transfer_data_ver2(data_t *m_dest, data_t *m_source, 
                            component_type_t m_destination_type, component_type_t m_source_type, 
                            data_type_t m_data_type, stationary_type_t m_stationary_type, action_type_t m_action_type, bool last_component);

    /* TODO: update for the next version */
    void transfer_data_network_on_chip(data_t *m_dest, data_t *m_source, unsigned m_dest_offset, unsigned m_source_offset,
                                       component_type_t m_destination_type, component_type_t m_source_type, 
                                       data_type_t m_data_type, stationary_type_t m_stationary_type, action_type_t m_action_type);
#endif

    std::vector<unsigned> calculate_parameter_size(component_type_t m_component_type);

    compression_type_t compression_type;                                                // Data format 

    // Tile size
    std::vector<std::vector<unsigned>> tile_size;                                       // Data tile size at the components

#ifdef FUNCTIONAL
    // The number of zero data
    std::vector<unsigned> num_zeros;                                                    // The number of zero value elements in the case of functional simulation
#endif

    /* PE instantiates for ver2 */
    std::vector<unsigned> parameters_pe;                                                // DNN parameters 
    std::vector<std::list<unsigned>> offset_pe;
    std::vector<unsigned> iterations_pe;                                                // The number of iterations at PE

    /* PE scheduler */
    std::vector<std::list<unsigned>> offset_size_pe;                                    // The number of data in PE
    std::list<unsigned> input_offset_pe;                                                // Input offsets for version 1
    std::list<unsigned> weight_offset_pe;                                               // Weight offsets for version 1
    std::list<unsigned> output_offset_pe;                                               // Output offsets for version 1
    std::map<unsigned, bool> output_read_pe;                                            // To check output read or initialize

    /* PE array scheduler */
    std::vector<unsigned> offset_size_pe_array;                                         // The number of data in PE array
    std::vector<unsigned> input_offset_pe_array;                                        // Offsets of input data
    std::vector<unsigned> weight_offset_pe_array;                                       // Offsets of weight
    std::vector<unsigned> output_offset_pe_array;                                       // Offsets of output data
    
    std::vector<bool> read_tile_granular_pe_input;                                      // Check tile granular input data in PE array
    std::vector<bool> read_tile_granular_pe_weight;                                     // Check tile granular weight in PE array
    std::vector<bool> read_tile_granular_pe_output;                                     // Check tile granular output data in PE array
    std::vector<unsigned> num_tile_granular_data_pe_array;                              // The number of different tile-granular data


    /* Global buffer instantiates for ver2 */
    std::vector<unsigned> parameters_global_buffer;                                     // DNN parameters at the global buffer
    std::vector<unsigned> offset_global_buffer;                                         // Data offset at the global buffer
    std::vector<unsigned> iterations_global_buffer;                                     // The number of iterations at the global buffer

    /* Global buffer scheduler */
    std::vector<std::list<unsigned>> offset_size_global_buffer;                         // The number of data in the global buffer
    std::list<unsigned> input_offset_global_buffer;
    std::list<unsigned> weight_offset_global_buffer;
    std::list<unsigned> output_offset_global_buffer;
    std::map<unsigned, bool> output_read_global_buffer;                                 // To check read or initialize the output data

    /* Multi-chip scheduler */
    std::vector<unsigned> offset_size_multi_chip;                                       // The number of data in the on-chip processors
    std::vector<unsigned> input_offset_multi_chip;                                      // Offsets of input data
    std::vector<unsigned> weight_offset_multi_chip;                                     // Offsets of weight
    std::vector<unsigned> output_offset_multi_chip;                                     // Offsets of output data

    std::vector<bool> read_tile_granular_chip_input;                                    // Check tile granular input data in multi chip
    std::vector<bool> read_tile_granular_chip_weight;                                   // Check tile granular weight in multi chip
    std::vector<bool> read_tile_granular_chip_output;                                   // Check tile granular output data in multi-chip
    std::vector<unsigned> num_tile_granular_data_chip;                                  // The number of different tile-granular data

    /* Off-chip memory instantiates for ver2 */
    std::vector<unsigned> parameters_dram;                                              // DNN parameters at the off-chip memory
    std::vector<unsigned> offset_dram;                                                  // Data offset at the off-chip memory
    std::vector<unsigned> iterations_dram;                                              // The number of iterations at the off-chip memory

    /* Off-chip memory scheduler */
    std::vector<std::list<unsigned>> offset_size_dram;                                  // The number of data in the off-chip memory
    std::list<unsigned> input_offset_dram;
    std::list<unsigned> weight_offset_dram;
    std::list<unsigned> output_offset_dram;
    std::map<unsigned, bool> output_read_dram;                                          // To check transfer output data to multi-chip processor or initialize at the multi-chip processor

    mapping_table_t *mapping_table;

    unsigned num_active_mac;                                                            // The number of active MACs
    unsigned num_active_pe_x;                                                           // The number of active PEs in X dimension
    unsigned num_active_pe_y;                                                           // The number of active PEs in Y dimension
    unsigned num_active_chips_x;                                                        // The number of active chips in X dimension
    unsigned num_active_chips_y;                                                        // The number of active chips in Y dimension
    layer_name_t layer_name;

private:

    std::map<unsigned, bool> update_output_read(std::list<unsigned>m_output_offset);
   
    std::vector<std::vector<unsigned>> calculate_tile_size();
    
    // Calculate counter and the offsets of data
    // Case 1. Undefined stationary
    std::vector<std::list<unsigned>> calculate_offset_undefined_stationary(std::list<unsigned> *m_input_offsets, std::list<unsigned> *m_weight_offsets, std::list<unsigned> *m_output_offsets, 
                                                                           component_type_t m_dest_type, component_type_t m_source_type);

    // Case 2. Input stationary
    std::vector<std::list<unsigned>> calculate_offset_input_stationary(std::list<unsigned> *m_input_offsets, std::list<unsigned> *m_weight_offsets, std::list<unsigned> *m_output_offsets, 
                                                                       component_type_t m_dest_type, component_type_t m_source_type);
    // Case 3. Weight stationary
    std::vector<std::list<unsigned>> calculate_offset_weight_stationary(std::list<unsigned> *m_input_offsets, std::list<unsigned> *m_weight_offsets, std::list<unsigned> *m_output_offsets, 
                                                                        component_type_t m_dest_type, component_type_t m_source_type);

    // Case 4. Output stationary
    std::vector<std::list<unsigned>> calculate_offset_output_stationary(std::list<unsigned> *m_input_offsets, std::list<unsigned> *m_weight_offsets, std::list<unsigned> *m_output_offsets, 
                                                                        component_type_t m_dest_type, component_type_t m_source_type);
    // Offset calculation of spatial components
    void calculate_offset_network_on_chip(std::vector<unsigned> *m_input_offsets, std::vector<unsigned> *m_weight_offsets, std::vector<unsigned> *m_output_offsets, 
                                          component_type_t m_destination_type, component_type_t m_source_type);

    /* Functions for ver2 */
    // Counter and Offset calculation
    // Case 1. Undefined stationary
    std::vector<std::list<unsigned>> calculate_counter_undefined_stationary_ver2(component_type_t m_destination_type, component_type_t m_source_type, std::list<unsigned> *m_output_offset);

    std::vector<std::list<unsigned>> calculate_offset_undefined_stationary_ver2(data_type_t m_data_type, component_type_t m_destination_type, component_type_t m_source_type,
                                                                                             std::vector<unsigned> *m_offsets, std::vector<unsigned> *m_params, std::vector<unsigned> *m_iteration, 
                                                                                             std::string m_parameter_order, bool m_last_component);
    // Case 2. Input stationary
    std::vector<std::list<unsigned>> calculate_counter_input_stationary_ver2(component_type_t m_destination_type, component_type_t m_source_type, std::list<unsigned> *m_output_offset);

    std::vector<std::list<unsigned>> calculate_offset_input_stationary_ver2(data_type_t m_data_type, component_type_t m_destination_type, component_type_t m_source_type,
                                                                                         std::vector<unsigned> *m_offsets, std::vector<unsigned> *m_params, std::vector<unsigned> *m_iteration, 
                                                                                         std::string m_parameter_order, bool m_last_component);
    // Case 3. Weight stationary
    std::vector<std::list<unsigned>> calculate_counter_weight_stationary_ver2(component_type_t m_destination_type, component_type_t m_source_type, std::list<unsigned> *m_output_offset);

    std::vector<std::list<unsigned>> calculate_offset_weight_stationary_ver2(data_type_t m_data_type, component_type_t m_destination_type, component_type_t m_source_type, 
                                                                                          std::vector<unsigned> *m_offsets, std::vector<unsigned> *m_params, std::vector<unsigned> *m_iteration, std::vector<unsigned> m_counter,
                                                                                          std::string m_parameter_order, bool m_last_component);
    // Case 4. Output stationary
    std::vector<std::list<unsigned>> calculate_counter_output_stationary_ver2(component_type_t m_destination_type, component_type_t m_source_type, std::list<unsigned> *m_output_offset);

    std::vector<std::list<unsigned>> calculate_offset_output_stationary_ver2(data_type_t m_data_type, component_type_t m_destination_type, component_type_t m_source_type, 
                                                                                          std::vector<unsigned> *m_offsets, std::vector<unsigned> *m_params, std::vector<unsigned> *m_iteration, 
                                                                                          std::string m_parameter_order, bool m_last_component);
    
#ifdef FUNCTIONAL

    // Transfer input data from lower to higher component hierarchy.
    void input_data_load(data_t *m_dest, data_t *m_source, unsigned m_dest_offset, unsigned m_source_offset, 
                         component_type_t m_destination_type, component_type_t m_source_type);

    // Transfer weight from lower to higher component hierarchy.
    void weight_data_load(data_t *m_dest, data_t *m_source, unsigned m_dest_offset, unsigned m_source_offset, 
                          component_type_t m_destination_type, component_type_t m_source_type);

    // Transfer output data from lower to higher component hierarchy.
    void output_data_load(data_t *m_dest, data_t *m_source, unsigned m_dest_offset, unsigned m_source_offset, 
                          component_type_t m_destination_type, component_type_t m_source_type);

    // Write back
    void output_data_store(data_t *m_dest, data_t *m_source, unsigned m_dest_offset, unsigned m_source_offset, 
                           component_type_t m_destination_type, component_type_t m_source_type);
#endif

};

#endif

    /* TODO : update for the next version */
    //std::vector<std::list<unsigned>> calculate_counter_undefined_stationary(component_type_t m_destination_type, component_type_t m_source_type, std::list<unsigned> *m_output_offset);
    //std::vector<std::list<unsigned>> calculate_offset_undefined_stationary(data_type_t m_data_type, component_type_t m_destination_type, component_type_t m_source_type, std::vector<unsigned> *m_offsets, std::vector<unsigned> *m_params, std::vector<unsigned> *m_iteration, std::string m_parameter_order);

    /* TODO : update for the next version */
    //std::vector<std::list<unsigned>> calculate_counter_input_stationary(component_type_t m_destination_type, component_type_t m_source_type, std::list<unsigned> *m_output_offset);
    //void calculate_offset_input_stationary(data_type_t m_data_type, component_type_t m_destination_type, component_type_t m_source_type, std::vector<unsigned> *m_offsets, std::vector<unsigned> *m_params, std::vector<unsigned> *m_iteration, std::string m_parameter_order, bool m_last_component);

    /* TODO : update for the next version */
    //std::vector<std::list<unsigned>> calculate_counter_weight_stationary(component_type_t m_destination_type, component_type_t m_source_type, std::list<unsigned> *m_output_offset);
    //void calculate_offset_weight_stationary(data_type_t m_data_type, component_type_t m_destination_type, component_type_t m_source_type, std::vector<unsigned> *m_offsets, std::vector<unsigned> *m_params, std::vector<unsigned> *m_iteration, std::vector<unsigned> m_counter, std::string m_parameter_order, bool m_last_component);

    /* TODO : update for the next version */
    //std::vector<std::list<unsigned>> calculate_counter_output_stationary(component_type_t m_destination_type, component_type_t m_source_type, std::list<unsigned> *m_output_offset);
    //void calculate_offset_output_stationary(data_type_t m_data_type, component_type_t m_destination_type, component_type_t m_source_type, std::vector<unsigned> *m_offsets, std::vector<unsigned> *m_params, std::vector<unsigned> *m_iteration, std::string m_parameter_order, bool m_last_component);
