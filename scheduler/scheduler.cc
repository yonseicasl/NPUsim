#include "scheduler.h"


scheduler_t::scheduler_t(mapping_table_t *m_mapping_table, stationary_type_t pe_stationary, std::string pe_parameter_order,
                         stationary_type_t pe_array_stationary, std::string pe_array_parameter_order,
                         stationary_type_t multi_chip_stationary, std::string multi_chip_parameter_order) :
    compression_type(compression_type_t::DENSE),
    num_active_mac(1),
    num_active_pe_x(1),
    num_active_pe_y(1),
    num_active_chips_x(1),
    num_active_chips_y(1),
    layer_name(UNDEFINED_LAYER) {
    
    init(m_mapping_table, pe_stationary, pe_parameter_order, 
                          pe_array_stationary, pe_array_parameter_order, 
                          multi_chip_stationary, multi_chip_parameter_order);
}

scheduler_t::~scheduler_t() {
}

// Initialize the scheduler
void scheduler_t::init(mapping_table_t *m_mapping_table, stationary_type_t pe_stationary, std::string pe_parameter_order, 
                                                         stationary_type_t pe_array_stationary, std::string pe_array_parameter_order,
                                                         stationary_type_t multi_chip_stationary, std::string multi_chip_parameter_order) {

    mapping_table = m_mapping_table;

    // Calculate time-granular data size 
    tile_size.reserve(component_type_t::NUM_COMPONENT_TYPES);
    tile_size = calculate_tile_size();

#ifdef FUNCTIONAL
    // Initialize the number of zero-value data for a compression.
    num_zeros.reserve(data_type_t::NUM_DATA_TYPES);
    num_zeros.assign(data_type_t::NUM_DATA_TYPES, 0);
#endif

    /* Initialize PE scheduler */

    /* update for NPUsim ver2 */
    parameters_pe.reserve(parameter_type_t::NUM_PARAMETER_TYPES);
    parameters_pe.assign(parameter_type_t::NUM_PARAMETER_TYPES, 0);

    offset_pe.reserve(data_type_t::NUM_DATA_TYPES);
    offset_pe.assign(data_type_t::NUM_DATA_TYPES, 0);

    iterations_pe.reserve(data_type_t::NUM_DATA_TYPES);
    iterations_pe.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the offsets at PE.
    if(pe_stationary == stationary_type_t::INPUT_STATIONARY) {
        /* Counter and offset calculator : Version 1 */
        offset_size_pe = calculate_offset_input_stationary(&input_offset_pe, &weight_offset_pe, &output_offset_pe, 
                                                           component_type_t::MAC, component_type_t::PE);
        /* Counter calculator : Version 2 */
        //offset_size_pe = calculate_counter_input_stationary_ver2(component_type_t::MAC, component_type_t::PE, &output_offset_pe);
    } else if(pe_stationary == stationary_type_t::WEIGHT_STATIONARY) {
        /* Counter and offset calculator : Version 1 */
        offset_size_pe = calculate_offset_weight_stationary(&input_offset_pe, &weight_offset_pe, &output_offset_pe, 
                                                            component_type_t::MAC, component_type_t::PE);
        /* Counter calculator : Version 2 */
        //offset_size_pe = calculate_counter_weight_stationary_ver2(component_type_t::MAC, component_type_t::PE, &output_offset_pe);
    } else if(pe_stationary == stationary_type_t::OUTPUT_STATIONARY) {
        /* Counter and offset calculator : Version 1 */
        offset_size_pe = calculate_offset_output_stationary(&input_offset_pe, &weight_offset_pe, &output_offset_pe, 
                                                            component_type_t::MAC, component_type_t::PE);
        /* Counter calculator : Version 2 */
        //offset_size_pe = calculate_counter_output_stationary_ver2(component_type_t::MAC, component_type_t::PE, &output_offset_pe);
    } else if(pe_stationary == stationary_type_t::UNDEFINED_STATIONARY) {
        std::cerr << "No Local Reuse is not available on the current version" << std::endl;
        exit(1);
    }

    output_read_pe = update_output_read(output_offset_pe);
    
    //Active MAC calculation
    num_active_mac = mapping_table->calculate_active_component(component_type_t::MAC);

    /* Initialize PE array scheduler */

    offset_size_pe_array.reserve(data_type_t::NUM_DATA_TYPES);

    // Calculate the offsets of PE array to spread out the data to PE.
    calculate_offset_network_on_chip(&input_offset_pe_array, &weight_offset_pe_array, &output_offset_pe_array, component_type_t::PE, component_type_t::PE_Y);

    // Check tile-granular input data is read or not
    read_tile_granular_pe_input.reserve(input_offset_pe_array.size());
    read_tile_granular_pe_input.assign(input_offset_pe_array.size(), false);

    // Check tile-granular weight is read or not
    read_tile_granular_pe_weight.reserve(weight_offset_pe_array.size());
    read_tile_granular_pe_weight.assign(weight_offset_pe_array.size(), false);

    // Check tile-granular output data is read or not
    read_tile_granular_pe_output.reserve(output_offset_pe_array.size());
    read_tile_granular_pe_output.assign(output_offset_pe_array.size(), false);

    // Active PE calculation
    num_active_pe_y = mapping_table->calculate_active_component(component_type_t::PE_Y);
    num_active_pe_x = mapping_table->calculate_active_component(component_type_t::PE_X);

    // The number of tile-granular data
    num_tile_granular_data_pe_array.reserve(data_type_t::NUM_DATA_TYPES);
    num_tile_granular_data_pe_array.assign(data_type_t::NUM_DATA_TYPES, 1);

    // Calculate the number of tile-granular data
    mapping_table->calculate_num_tile_granular_data(component_type_t::PE_X, &num_tile_granular_data_pe_array);
    mapping_table->calculate_num_tile_granular_data(component_type_t::PE_Y, &num_tile_granular_data_pe_array);

    /* Initialize Global buffer scheduler */

    /* Update for NPUsim ver2 */
    parameters_global_buffer.reserve(parameter_type_t::NUM_PARAMETER_TYPES);
    parameters_global_buffer.assign(parameter_type_t::NUM_PARAMETER_TYPES, 0);

    offset_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    offset_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0);

    iterations_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    iterations_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0);

    offset_size_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);

    // Calculate the offsets at the global buffer
    if(pe_array_stationary == stationary_type_t::INPUT_STATIONARY) {
        /* Counter and offset calculator : Version 1 */
        offset_size_global_buffer = calculate_offset_input_stationary(&input_offset_global_buffer, &weight_offset_global_buffer, &output_offset_global_buffer, 
                                                                      component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER);
        /* Counter calculator : Version 2 */
        //offset_size_global_buffer = calculate_counter_input_stationary_ver2(component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, &output_offset_global_buffer);
    } else if(pe_array_stationary == stationary_type_t::WEIGHT_STATIONARY) {
        /* Counter and offset calculator : Version 1 */
        offset_size_global_buffer = calculate_offset_weight_stationary(&input_offset_global_buffer, &weight_offset_global_buffer, &output_offset_global_buffer, 
                                                                       component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER);
        /* Counter calculator : Version 2 */
        //offset_size_global_buffer = calculate_counter_weight_stationary_ver2(component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, &output_offset_global_buffer);
    } else if(pe_array_stationary == stationary_type_t::OUTPUT_STATIONARY) {
        /* Counter and offset calculator : Version 1 */
        offset_size_global_buffer = calculate_offset_output_stationary(&input_offset_global_buffer, &weight_offset_global_buffer, &output_offset_global_buffer, 
                                                                        component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER);
        /* Counter calculator : Version 2 */
        //offset_size_global_buffer = calculate_counter_output_stationary_ver2(component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, &output_offset_global_buffer);   
    } else if(pe_array_stationary == stationary_type_t::UNDEFINED_STATIONARY) {
        std::cerr << "No Local Reuse is not available on the current version" << std::endl;
        exit(1);
    }

    output_read_global_buffer = update_output_read(output_offset_global_buffer);

    /* Initialize chip-level processor scheduler */

    offset_size_multi_chip.reserve(data_type_t::NUM_DATA_TYPES);

    // calculate the offsets of chip-level processor
    calculate_offset_network_on_chip(&input_offset_multi_chip, &weight_offset_multi_chip, &output_offset_multi_chip, component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y);

    // Check tile-granular input data is read or not
    read_tile_granular_chip_input.reserve(input_offset_multi_chip.size());
    read_tile_granular_chip_input.assign(input_offset_multi_chip.size(), false);

    // Check tile-granular weight is read or not
    read_tile_granular_chip_weight.reserve(weight_offset_multi_chip.size());
    read_tile_granular_chip_weight.assign(weight_offset_multi_chip.size(), false);

    // Check tile-granular output data is read or not
    read_tile_granular_chip_output.reserve(output_offset_multi_chip.size());
    read_tile_granular_chip_output.assign(output_offset_multi_chip.size(), false);

    // Active chips calculation
    num_active_chips_y = mapping_table->calculate_active_component(component_type_t::CHIPS_Y);
    num_active_chips_x = mapping_table->calculate_active_component(component_type_t::CHIPS_X);

    // The number of tile-granular data
    num_tile_granular_data_chip.reserve(data_type_t::NUM_DATA_TYPES);
    num_tile_granular_data_chip.assign(data_type_t::NUM_DATA_TYPES, 1);

    // Calculate the number of tile-granular data
    mapping_table->calculate_num_tile_granular_data(component_type_t::CHIPS_X, &num_tile_granular_data_chip);
    mapping_table->calculate_num_tile_granular_data(component_type_t::CHIPS_Y, &num_tile_granular_data_chip);

    /**** DRAM ****/

    /* Update for NPUsim ver2 */
    parameters_dram.reserve(parameter_type_t::NUM_PARAMETER_TYPES);
    parameters_dram.assign(parameter_type_t::NUM_PARAMETER_TYPES, 0);

    offset_dram.reserve(data_type_t::NUM_DATA_TYPES);
    offset_dram.assign(data_type_t::NUM_DATA_TYPES, 0);

    iterations_dram.reserve(data_type_t::NUM_DATA_TYPES);
    iterations_dram.assign(data_type_t::NUM_DATA_TYPES, 0);

    offset_size_dram.reserve(data_type_t::NUM_DATA_TYPES);

    // Initialize the offsets at the off-chip memory
    if(multi_chip_stationary == stationary_type_t::INPUT_STATIONARY) {
        /* Counter and offset calculator : Version 1 */
        offset_size_dram = calculate_offset_input_stationary(&input_offset_dram, &weight_offset_dram, &output_offset_dram,
                                                             component_type_t::CHIPS_Y, component_type_t::DRAM);
        /* Counter calculator : Version 2 */
        //offset_size_dram = calculate_counter_input_stationary_ver2(component_type_t::CHIPS_Y, component_type_t::DRAM, &output_offset_dram);   
    } else if(multi_chip_stationary == stationary_type_t::WEIGHT_STATIONARY) {
        /* Counter and offset calculator : Version 1 */
        offset_size_dram = calculate_offset_weight_stationary(&input_offset_dram, &weight_offset_dram, &output_offset_dram,
                                                              component_type_t::CHIPS_Y, component_type_t::DRAM);
        /* Counter calculator : Version 2 */
        //offset_size_dram = calculate_counter_weight_stationary_ver2(component_type_t::CHIPS_Y, component_type_t::DRAM, &output_offset_dram);
    } else if(multi_chip_stationary == stationary_type_t::OUTPUT_STATIONARY) {
        /* Counter and offset calculator : Version 1 */
        offset_size_dram = calculate_offset_output_stationary(&input_offset_dram, &weight_offset_dram, &output_offset_dram,
                                                              component_type_t::CHIPS_Y, component_type_t::DRAM);
        /* Counter calculator : Version 2 */
        //offset_size_dram = calculate_counter_output_stationary_ver2(component_type_t::CHIPS_Y, component_type_t::DRAM, &output_offset_dram);
    } else if(multi_chip_stationary == stationary_type_t::UNDEFINED_STATIONARY) {
        std::cerr << "No Local Reuse is not available on the current version" << std::endl;
        exit(1);
    }
    output_read_dram = update_output_read(output_offset_dram);
}

void scheduler_t::print_stats() {
    std::cout << "* Tile size" << std::endl;
    std::cout << "========== MAC =========" << std::endl;
    std::cout << "# Input data  : " << tile_size[component_type_t::MAC][data_type_t::INPUT]  << std::endl;
    std::cout << "# Weight      : " << tile_size[component_type_t::MAC][data_type_t::WEIGHT] << std::endl;
    std::cout << "# Output data : " << tile_size[component_type_t::MAC][data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;

    std::cout << "===== Local buffer =====" << std::endl;
    std::cout << "# Input data  : " << tile_size[component_type_t::PE][data_type_t::INPUT]  << std::endl;
    std::cout << "# Weight      : " << tile_size[component_type_t::PE][data_type_t::WEIGHT] << std::endl;
    std::cout << "# Output data : " << tile_size[component_type_t::PE][data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;

    std::cout << "======= PE array =======" << std::endl;
    std::cout << "# Input data  : " << tile_size[component_type_t::PE_Y][data_type_t::INPUT]  << std::endl;
    std::cout << "# Weight      : " << tile_size[component_type_t::PE_Y][data_type_t::WEIGHT] << std::endl;
    std::cout << "# Output data : " << tile_size[component_type_t::PE_Y][data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;

    std::cout << "===== Global buffer ====" << std::endl;
    std::cout << "# Input data  : " << tile_size[component_type_t::GLOBAL_BUFFER][data_type_t::INPUT]  << std::endl;
    std::cout << "# Weight      : " << tile_size[component_type_t::GLOBAL_BUFFER][data_type_t::WEIGHT] << std::endl;
    std::cout << "# Output data : " << tile_size[component_type_t::GLOBAL_BUFFER][data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;

    std::cout << "===== Multi Chip =======" << std::endl;
    std::cout << "# Input data  : " << tile_size[component_type_t::CHIPS_Y][data_type_t::INPUT]  << std::endl;
    std::cout << "# Weight      : " << tile_size[component_type_t::CHIPS_Y][data_type_t::WEIGHT] << std::endl;
    std::cout << "# Output data : " << tile_size[component_type_t::CHIPS_Y][data_type_t::OUTPUT] << std::endl;
    std::cout << std::endl;

    std::cout << "========= DRAM =========" << std::endl;
    std::cout << "# Input data  : " << tile_size[component_type_t::DRAM][data_type_t::INPUT]  << std::endl;
    std::cout << "# Weight      : " << tile_size[component_type_t::DRAM][data_type_t::WEIGHT] << std::endl;
    std::cout << "# Output data : " << tile_size[component_type_t::DRAM][data_type_t::OUTPUT] << std::endl;
    std::cout << "========================" << std::endl;
    std::cout << std::endl;
}

void scheduler_t::print_stats(std::ofstream &m_output_file) {
     m_output_file << "* Tile size" << std::endl;
     m_output_file << "========== MAC =========" << std::endl;
     m_output_file << "# Input data  : " << tile_size[component_type_t::MAC][data_type_t::INPUT]  << std::endl;
     m_output_file << "# Weight      : " << tile_size[component_type_t::MAC][data_type_t::WEIGHT] << std::endl;
     m_output_file << "# Output data : " << tile_size[component_type_t::MAC][data_type_t::OUTPUT] << std::endl;
     m_output_file << "========================" << std::endl;
     m_output_file << std::endl;

     m_output_file << "===== Local buffer =====" << std::endl;
     m_output_file << "# Input data  : " << tile_size[component_type_t::PE][data_type_t::INPUT]  << std::endl;
     m_output_file << "# Weight      : " << tile_size[component_type_t::PE][data_type_t::WEIGHT] << std::endl;
     m_output_file << "# Output data : " << tile_size[component_type_t::PE][data_type_t::OUTPUT] << std::endl;
     m_output_file << "========================" << std::endl;
     m_output_file << std::endl;

     m_output_file << "======= PE array =======" << std::endl;
     m_output_file << "# Input data  : " << tile_size[component_type_t::PE_Y][data_type_t::INPUT]  << std::endl;
     m_output_file << "# Weight      : " << tile_size[component_type_t::PE_Y][data_type_t::WEIGHT] << std::endl;
     m_output_file << "# Output data : " << tile_size[component_type_t::PE_Y][data_type_t::OUTPUT] << std::endl;
     m_output_file << "========================" << std::endl;
     m_output_file << std::endl;

     m_output_file << "===== Global buffer ====" << std::endl;
     m_output_file << "# Input data  : " << tile_size[component_type_t::GLOBAL_BUFFER][data_type_t::INPUT]  << std::endl;
     m_output_file << "# Weight      : " << tile_size[component_type_t::GLOBAL_BUFFER][data_type_t::WEIGHT] << std::endl;
     m_output_file << "# Output data : " << tile_size[component_type_t::GLOBAL_BUFFER][data_type_t::OUTPUT] << std::endl;
     m_output_file << "========================" << std::endl;
     m_output_file << std::endl;

     m_output_file << "===== Multi Chip =======" << std::endl;
     m_output_file << "# Input data  : " << tile_size[component_type_t::CHIPS_Y][data_type_t::INPUT]  << std::endl;
     m_output_file << "# Weight      : " << tile_size[component_type_t::CHIPS_Y][data_type_t::WEIGHT] << std::endl;
     m_output_file << "# Output data : " << tile_size[component_type_t::CHIPS_Y][data_type_t::OUTPUT] << std::endl;
     m_output_file << std::endl;

     m_output_file << "========= DRAM =========" << std::endl;
     m_output_file << "# Input data  : " << tile_size[component_type_t::DRAM][data_type_t::INPUT]  << std::endl;
     m_output_file << "# Weight      : " << tile_size[component_type_t::DRAM][data_type_t::WEIGHT] << std::endl;
     m_output_file << "# Output data : " << tile_size[component_type_t::DRAM][data_type_t::OUTPUT] << std::endl;
     m_output_file << "========================" << std::endl;
     m_output_file << std::endl;
}


#ifdef FUNCTIONAL
void scheduler_t::transfer_data(data_t *m_dest, data_t *m_source, 
                                unsigned m_dest_offset, unsigned m_source_offset, 
                                component_type_t m_dest_type, component_type_t m_source_type, 
                                data_type_t m_data_type, stationary_type_t m_stationary_type, action_type_t m_action_type) {
    // Transfer input data.
    if(m_data_type == data_type_t::INPUT) {
        input_data_load(m_dest, m_source, 
                        m_dest_offset, m_source_offset, 
                        m_dest_type, m_source_type);
    }
    // Transfer weight
    else if(m_data_type == data_type_t::WEIGHT) {
        weight_data_load(m_dest, m_source, 
                         m_dest_offset, m_source_offset, 
                         m_dest_type, m_source_type);
    }
    // Transfer output data.
    else if(m_data_type == data_type_t::OUTPUT) {
        // Output data transfer
        if(m_action_type == action_type_t::LOAD) {
            output_data_load(m_dest, m_source, 
                             m_dest_offset, m_source_offset, 
                             m_dest_type, m_source_type);
        }
        // Write back output data
        else if(m_action_type == action_type_t::STORE) {
            output_data_store(m_dest, m_source, 
                              m_dest_offset, m_source_offset, 
                              m_dest_type, m_source_type);
        }
    }
}

void scheduler_t::transfer_data_ver2(data_t *m_dest, data_t *m_source, 
                                     component_type_t m_destination_type, component_type_t m_source_type, 
                                     data_type_t m_data_type, stationary_type_t m_stationary_type, action_type_t m_action_type, bool last_component) {
    // Calculate data offsets.
    std::string parameter_order = "BCPQKRS";

    std::vector<unsigned> *t_offsets;
    std::vector<unsigned> *t_params;
    std::vector<unsigned> *t_iterations;
    std::vector<unsigned> t_counters(data_type_t::NUM_DATA_TYPES);
    if(m_destination_type == component_type_t::MAC) {
        t_offsets = &offset_pe;
        t_params = &parameters_pe;
        t_iterations = &iterations_pe;
        t_counters[data_type_t::INPUT] = offset_size_pe[data_type_t::INPUT].front();
        t_counters[data_type_t::WEIGHT] = offset_size_pe[data_type_t::WEIGHT].front();
        t_counters[data_type_t::OUTPUT] = offset_size_pe[data_type_t::OUTPUT].front();
    } else if(m_destination_type == component_type_t::PE_Y) {
        t_offsets = &offset_global_buffer;
        t_params = &parameters_global_buffer;
        t_iterations = &iterations_global_buffer;
        t_counters[data_type_t::INPUT] = offset_size_global_buffer[data_type_t::INPUT].front();
        t_counters[data_type_t::WEIGHT] = offset_size_global_buffer[data_type_t::WEIGHT].front();
        t_counters[data_type_t::OUTPUT] = offset_size_global_buffer[data_type_t::OUTPUT].front();
    } else if(m_destination_type == component_type_t::CHIPS_Y) {
        t_offsets = &offset_dram;
        t_params = &parameters_dram;
        t_iterations = &iterations_dram;
        t_counters[data_type_t::INPUT] = offset_size_dram[data_type_t::INPUT].front();
        t_counters[data_type_t::WEIGHT] = offset_size_dram[data_type_t::WEIGHT].front();
        t_counters[data_type_t::OUTPUT] = offset_size_dram[data_type_t::OUTPUT].front();
    } else {
        std::cerr << "No matching accelerator components" << std::endl;
        exit(1);
    }

    if(m_stationary_type == stationary_type_t::INPUT_STATIONARY) {
        calculate_offset_input_stationary_ver2(m_data_type, m_destination_type, m_source_type, 
                                               t_offsets, t_params, 
                                               t_iterations, t_counters,
                                               parameter_order, last_component);
    } else if(m_stationary_type == stationary_type_t::WEIGHT_STATIONARY) {
        calculate_offset_weight_stationary_ver2(m_data_type, m_destination_type, m_source_type, 
                                                t_offsets, t_params, 
                                                t_iterations, t_counters, 
                                                parameter_order, last_component);
    } else if(m_stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
        calculate_offset_output_stationary_ver2(m_data_type, m_destination_type, m_source_type, 
                                                t_offsets, t_params, 
                                                t_iterations, t_counters, 
                                                parameter_order, last_component);
    } 

    // Transfer input data.
    if(m_data_type == data_type_t::INPUT) {
        input_data_load(m_dest, m_source, 0, t_offsets->at(data_type_t::INPUT), m_destination_type, m_source_type);
    }
    // Transfer weight
    else if(m_data_type == data_type_t::WEIGHT) {
        weight_data_load(m_dest, m_source, 0, t_offsets->at(data_type_t::WEIGHT), m_destination_type, m_source_type);
    }
    // Transfer output data.
    else if(m_data_type == data_type_t::OUTPUT) {
        if(m_action_type == action_type_t::LOAD) {
            output_data_load(m_dest, m_source, 0, t_offsets->at(data_type_t::OUTPUT), m_destination_type, m_source_type);
        }
        else if(m_action_type == action_type_t::STORE) {
            output_data_store(m_dest, m_source, t_offsets->at(data_type_t::OUTPUT), 0, m_destination_type, m_source_type);
        }
    }
}

#endif

std::vector<unsigned> scheduler_t::calculate_parameter_size(component_type_t m_component_type) {

    std::vector<unsigned> parameters(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    parameters = mapping_table->calculate_parameter_size(m_component_type);
    return parameters;
}

// Initialize output read to false (indicate output data can be initialized without data transfer)
std::map<unsigned, bool> scheduler_t::update_output_read(std::list<unsigned>m_output_offset) {

    std::map<unsigned, bool> output_read;
    for(auto it = m_output_offset.begin(); it != m_output_offset.end(); ++it) {
        if(output_read.find(*it) == output_read.end()) {
            output_read.insert(std::pair<unsigned, bool>(*it, false));
        }
    }
    return output_read;
}

// Calculate tile-granular data size
std::vector<std::vector<unsigned>> scheduler_t::calculate_tile_size() {

    std::vector<unsigned> parameters(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> tile_size(data_type_t::NUM_DATA_TYPES, 1);
    std::vector<std::vector<unsigned>> tile_sizes(component_type_t::NUM_COMPONENT_TYPES, std::vector<unsigned>(data_type_t::NUM_DATA_TYPES, 1));

    for(unsigned i = 0; i < component_type_t::NUM_COMPONENT_TYPES; i++) {
        parameters = mapping_table->calculate_parameter_size(static_cast<component_type_t>(i));
        tile_size[data_type_t::INPUT] = parameters[parameter_type_t::BATCH_SIZE]*parameters[parameter_type_t::INPUT_CHANNEL]
                                       *parameters[parameter_type_t::INPUT_HEIGHT]*parameters[parameter_type_t::INPUT_WIDTH];
        tile_size[data_type_t::WEIGHT] = parameters[parameter_type_t::OUTPUT_CHANNEL]*parameters[parameter_type_t::INPUT_CHANNEL]
                                        *parameters[parameter_type_t::FILTER_HEIGHT]*parameters[parameter_type_t::FILTER_WIDTH]
                                        /parameters[parameter_type_t::GROUP];
        tile_size[data_type_t::OUTPUT] = parameters[parameter_type_t::BATCH_SIZE]*parameters[parameter_type_t::OUTPUT_CHANNEL]
                                        *parameters[parameter_type_t::OUTPUT_HEIGHT]*parameters[parameter_type_t::OUTPUT_WIDTH];
        tile_sizes[i] = tile_size;
        parameters.assign(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    }

    return tile_sizes;
}


// Case 1 : Offset calculation in undefined stationary
//std::vector<std::list<unsigned>> scheduler_t::calculate_offset_undefined_stationary(std::list<unsigned> *m_input_offsets, std::list<unsigned> *m_weight_offsets, std::list<unsigned> *m_output_offsets,
//                                                                                    component_type_t m_destination_type, component_type_t m_source_type) {
//                                                                                                                        
//    std::vector<std::list<unsigned>> offset_size;                                                                       
//                                                                                                                        
//    return offset_size;                                                                                                 
//                                                                                                                        
//}         

// Case 2. Offset calculation in Input stationary 
std::vector<std::list<unsigned>> scheduler_t::calculate_offset_input_stationary(std::list<unsigned> *m_input_offsets, std::list<unsigned> *m_weight_offsets, std::list<unsigned> *m_output_offsets,
                                                                                component_type_t m_destination_type, component_type_t m_source_type) {


    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    std::vector<unsigned> dram_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    dram_param = mapping_table->calculate_parameter_size(component_type_t::DRAM);

    std::vector<std::list<unsigned>> offset_size;
    offset_size.reserve(data_type_t::NUM_DATA_TYPES);

    // Calculate the parameter size of destination and source.
    // The parameter sizes are used to calculate tile offset.
    dest_param = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);

    std::list<unsigned> input_offset_sizes, weight_offset_sizes, output_offset_sizes;
    unsigned input_offset_size = 0;

    // Calculate input hops
    unsigned height_hop = (dest_param[parameter_type_t::FILTER_HEIGHT] == source_param[parameter_type_t::FILTER_HEIGHT] &&
                           dest_param[parameter_type_t::INPUT_HEIGHT] < source_param[parameter_type_t::INPUT_HEIGHT]) ?
                           dest_param[parameter_type_t::STRIDE]*dest_param[parameter_type_t::OUTPUT_HEIGHT] : 1;
    unsigned width_hop  = (dest_param[parameter_type_t::FILTER_WIDTH] == source_param[parameter_type_t::FILTER_WIDTH] &&
                           dest_param[parameter_type_t::INPUT_WIDTH] < source_param[parameter_type_t::INPUT_WIDTH])  ?
                           dest_param[parameter_type_t::STRIDE]*dest_param[parameter_type_t::OUTPUT_WIDTH] : 1;
    unsigned stride     = dest_param[parameter_type_t::STRIDE];

    // The sequence of input data : B->C->H->W
    for(unsigned b = 0; b < source_param[parameter_type_t::BATCH_SIZE]; b += dest_param[parameter_type_t::BATCH_SIZE]) {
        for(unsigned g = 0; g < source_param[parameter_type_t::GROUP]; g += dest_param[parameter_type_t::GROUP]) {
            for(unsigned c = 0; c < source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; c += dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
                for(unsigned h = 0; h <= source_param[parameter_type_t::INPUT_HEIGHT] - dest_param[parameter_type_t::INPUT_HEIGHT];) {
                    for(unsigned w = 0; w <= source_param[parameter_type_t::INPUT_WIDTH] - dest_param[parameter_type_t::INPUT_WIDTH];) {

                        // Calculate the offsets of input data.
                        unsigned input_offset = b*source_param[parameter_type_t::GROUP]*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                              + g*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                              + c*source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                              + h*source_param[parameter_type_t::INPUT_WIDTH] + w;
                        m_input_offsets->emplace_back(input_offset);
                        input_offset_size++;

                        unsigned weight_offset_size = 0, output_offset_size = 0;
                        // The sequence of weight : C->K->R->S
                        for(unsigned k = 0; k < source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; k += dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
                            for(unsigned r = h%stride; r < source_param[parameter_type_t::FILTER_HEIGHT]; r += dest_param[parameter_type_t::FILTER_HEIGHT]*dest_param[parameter_type_t::STRIDE]) {
                                if(h >= r && source_param[parameter_type_t::INPUT_HEIGHT] - h >= source_param[parameter_type_t::FILTER_HEIGHT] - r
                                          && (h-r)/height_hop*dest_param[parameter_type_t::OUTPUT_HEIGHT] < source_param[parameter_type_t::OUTPUT_HEIGHT]) {
                                    for(unsigned s = w%stride; s < source_param[parameter_type_t::FILTER_WIDTH]; s += dest_param[parameter_type_t::FILTER_WIDTH]*dest_param[parameter_type_t::STRIDE]) {
                                        if(w >= s && source_param[parameter_type_t::INPUT_WIDTH] - w >= source_param[parameter_type_t::FILTER_WIDTH] - s
                                                  && (w-s)/width_hop*dest_param[parameter_type_t::OUTPUT_WIDTH] < source_param[parameter_type_t::OUTPUT_WIDTH]) {
                                            // Calculate the offsets of weight.
                                            unsigned weight_offset = g*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                                      *source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                                      *source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                                                   + k*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                                      *source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                                                   + c*source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                                                   + r*source_param[parameter_type_t::FILTER_WIDTH] + s;
                                            m_weight_offsets->emplace_back(weight_offset);
                                            weight_offset_size++;

                                            // Calculate the output height and output width to calculate the output data's offsets.

                                            unsigned p = (h - r)/dest_param[parameter_type_t::STRIDE], q = (w - s)/dest_param[parameter_type_t::STRIDE];
                                            // Calculate the offsets of output data.
                                            unsigned output_offset = b*source_param[parameter_type_t::GROUP]*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                                      *source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                                                   + g*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                                      *source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                                                   + k*source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                                                   + p*source_param[parameter_type_t::OUTPUT_WIDTH] + q;
                                            m_output_offsets->emplace_back(output_offset);
                                            output_offset_size++;
                                        }
                                    }
                                }
                            }
                        }
                        // Calculate the number of reuses of weight and output data.
                        if(weight_offset_size) {weight_offset_sizes.emplace_back(weight_offset_size);}
                        if(output_offset_size) {output_offset_sizes.emplace_back(output_offset_size);}
                        w += width_hop;
                    }
                    h += height_hop;
                }
            }
        }
    }
    input_offset_sizes.emplace_back(input_offset_size);

    offset_size.emplace_back(input_offset_sizes);
    offset_size.emplace_back(weight_offset_sizes);
    offset_size.emplace_back(output_offset_sizes);

    return offset_size;
}

// Case 3. Offset calculation in weight stationary 
std::vector<std::list<unsigned>> scheduler_t::calculate_offset_weight_stationary(std::list<unsigned> *m_input_offsets, std::list<unsigned> *m_weight_offsets, std::list<unsigned> *m_output_offsets,
                                                                                 component_type_t m_destination_type, component_type_t m_source_type) {

    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    // Calculate the parameter size of destination and source.
    // The parameter sizes are used to calculate tile offset.
    dest_param = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);

    std::vector<std::list<unsigned>> offset_size;
    offset_size.reserve(data_type_t::NUM_DATA_TYPES);

    unsigned stride = dest_param[parameter_type_t::STRIDE];

    unsigned weight_offset_size = 0;
    std::list<unsigned> input_offset_sizes, weight_offset_sizes, output_offset_sizes;
    for(unsigned g = 0; g < source_param[parameter_type_t::GROUP]; g += dest_param[parameter_type_t::GROUP]) {
        for(unsigned k = 0; k < source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; k += dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
            for(unsigned c = 0; c < source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; c += dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
                unsigned h = 0, w = 0;
                for(unsigned r = 0; r < source_param[parameter_type_t::FILTER_HEIGHT]; r += dest_param[parameter_type_t::FILTER_HEIGHT]) {
                    for(unsigned s = 0; s < source_param[parameter_type_t::FILTER_WIDTH]; s += dest_param[parameter_type_t::FILTER_WIDTH]) {
                        // Initialize the offset of weight.
                        // The sequence of offset : K->C->R->S
                        unsigned weight_offset = g*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                  *source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                               + k*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                  *source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                               + c*source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                               + r*source_param[parameter_type_t::FILTER_WIDTH] + s;
                        m_weight_offsets->emplace_back(weight_offset);
                        weight_offset_size++;

                        unsigned input_offset_size = 0, output_offset_size = 0;
                        for(unsigned b = 0; b < source_param[parameter_type_t::BATCH_SIZE]; b += dest_param[parameter_type_t::BATCH_SIZE]) {
                            for(unsigned p = 0; p < source_param[parameter_type_t::OUTPUT_HEIGHT]; p += dest_param[parameter_type_t::OUTPUT_HEIGHT]) {
                                for(unsigned q = 0; q < source_param[parameter_type_t::OUTPUT_WIDTH]; q += dest_param[parameter_type_t::OUTPUT_WIDTH]) {
                                    // Calculate the offsets of input data.
                                    // The sequence of input data : C->B->H->W
                                    unsigned input_offset = b*source_param[parameter_type_t::GROUP]*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                             *source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                                          + g*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                             *source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                                          + c*source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                                          + h*source_param[parameter_type_t::INPUT_WIDTH] + w
                                                          // The offset of input data should be different according to the position of output data.
                                                          + p*source_param[parameter_type_t::INPUT_WIDTH]*stride + q*stride;
                                    m_input_offsets->emplace_back(input_offset);
                                    input_offset_size++;

                                    // Calculate the offsets of output data.
                                    // The sequence of output data : K->B->P->Q
                                    unsigned output_offset = b*source_param[parameter_type_t::GROUP]*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                              *source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                                           + g*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                              *source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                                           + k*source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                                           + p*source_param[parameter_type_t::OUTPUT_WIDTH] + q;
                                    m_output_offsets->emplace_back(output_offset);
                                    output_offset_size++;
                                }
                            }
                        }
                        // Increase the index of input width.
                        w += dest_param[parameter_type_t::FILTER_WIDTH];
                        // Store the number of input and output data tiles.
                        // To indicate the number of reuse of weight.
                        input_offset_sizes.emplace_back(input_offset_size);
                        output_offset_sizes.emplace_back(output_offset_size);
                    }
                    // Increase the index of input height.
                    // Initialize the index of input width 0.
                    h += dest_param[parameter_type_t::FILTER_HEIGHT];
                    w = 0;
                }
            }
        }
    }
    weight_offset_sizes.emplace_back(weight_offset_size);

    offset_size.emplace_back(input_offset_sizes);
    offset_size.emplace_back(weight_offset_sizes);
    offset_size.emplace_back(output_offset_sizes);

    return offset_size;
}

// Case 4. Offset calculation in output stationary 
std::vector<std::list<unsigned>> scheduler_t::calculate_offset_output_stationary(std::list<unsigned> *m_input_offsets, std::list<unsigned> *m_weight_offsets, std::list<unsigned> *m_output_offsets,
                                                                                 component_type_t m_destination_type, component_type_t m_source_type) {

    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    // Calculate the parameter size of destination and source.
    // The parameter sizes are used to calculate tile offset.
    dest_param = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);

    std::vector<std::list<unsigned>> offset_size;
    offset_size.reserve(data_type_t::NUM_DATA_TYPES);

    unsigned stride = source_param[parameter_type_t::STRIDE];

    std::list<unsigned> input_offset_sizes, weight_offset_sizes, output_offset_sizes;
    unsigned output_offset_size = 0;

    for(unsigned b = 0; b < source_param[parameter_type_t::BATCH_SIZE]; b += dest_param[parameter_type_t::BATCH_SIZE]) {
        for(unsigned g = 0; g < source_param[parameter_type_t::GROUP]; g += dest_param[parameter_type_t::GROUP]) {
            for(unsigned k = 0; k < source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; k += dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
                unsigned h = 0, w = 0;
                for(unsigned p = 0; p < source_param[parameter_type_t::OUTPUT_HEIGHT]; p += dest_param[parameter_type_t::OUTPUT_HEIGHT]) {
                    for(unsigned q = 0; q < source_param[parameter_type_t::OUTPUT_WIDTH]; q += dest_param[parameter_type_t::OUTPUT_WIDTH]) {
                        // Calculate the offset of output data
                        // The sequence of output data : B->K->P->Q
                        unsigned output_offset = b*source_param[parameter_type_t::GROUP]*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                  *source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                               + g*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                  *source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                               + k*source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                               + p*source_param[parameter_type_t::OUTPUT_WIDTH] + q;
                        m_output_offsets->emplace_back(output_offset);
                        output_offset_size++;

                        unsigned input_offset_size = 0, weight_offset_size = 0;
                        for(unsigned c = 0; c < source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; c += dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
                            for(unsigned r = 0; r < source_param[parameter_type_t::FILTER_HEIGHT]; r += dest_param[parameter_type_t::FILTER_HEIGHT]) {
                                for(unsigned s = 0; s < source_param[parameter_type_t::FILTER_WIDTH]; s += dest_param[parameter_type_t::FILTER_WIDTH]) {
                                    // Calculate the offsets of input data.
                                    unsigned input_offset = b*source_param[parameter_type_t::GROUP]*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                             *source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                                          + g*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                             *source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                                          + c*source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                                          + h*source_param[parameter_type_t::INPUT_WIDTH] + w
                                                          + r*source_param[parameter_type_t::INPUT_WIDTH] + s;
                                    m_input_offsets->emplace_back(input_offset);
                                    input_offset_size++;

                                    // Calculate the offsets of weight.
                                    unsigned weight_offset = g*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                              *source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                              *source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                                           + k*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                              *source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                                           + c*source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                                           + r*source_param[parameter_type_t::FILTER_WIDTH] + s;
                                    m_weight_offsets->emplace_back(weight_offset);
                                    weight_offset_size++;
                                }
                            }
                        }
                        input_offset_sizes.emplace_back(input_offset_size);
                        weight_offset_sizes.emplace_back(weight_offset_size);

                        w += dest_param[parameter_type_t::OUTPUT_WIDTH]*stride;
                    }
                    h += dest_param[parameter_type_t::OUTPUT_HEIGHT]*stride;
                    w = 0;
                }
            }
        }
    }

    output_offset_sizes.emplace_back(output_offset_size);

    offset_size.emplace_back(input_offset_sizes);
    offset_size.emplace_back(weight_offset_sizes);
    offset_size.emplace_back(output_offset_sizes);

    return offset_size;
}


void scheduler_t::calculate_offset_network_on_chip(std::vector<unsigned> *m_input_offsets, std::vector<unsigned> *m_weight_offsets, std::vector<unsigned> *m_output_offsets, 
                                                   component_type_t m_destination_type, component_type_t m_source_type) {

    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    // Calculate the parameter size of destination and source.
    // The parameter sizes are used to calculate tile offset.
    dest_param = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);

    // The offset of Input data.
    for(unsigned b = 0; b < source_param[parameter_type_t::BATCH_SIZE]; b += dest_param[parameter_type_t::BATCH_SIZE]) {
        for(unsigned g = 0; g < source_param[parameter_type_t::GROUP]; g += dest_param[parameter_type_t::GROUP]) {
            for(unsigned c = 0; c < source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; c += dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
                for(unsigned h = 0; h < source_param[parameter_type_t::INPUT_HEIGHT]; h += dest_param[parameter_type_t::INPUT_HEIGHT]) {
                    for(unsigned w = 0; w < source_param[parameter_type_t::INPUT_WIDTH]; w += dest_param[parameter_type_t::INPUT_WIDTH]) {
                        unsigned input_offset = b*source_param[parameter_type_t::GROUP]*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                              + g*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                              + c*source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                              + h*source_param[parameter_type_t::INPUT_WIDTH] + w;
                        m_input_offsets->emplace_back(input_offset);
                    }
                }
            }
        }
    }

    // The offset of weight.
    for(unsigned g = 0; g < source_param[parameter_type_t::GROUP]; g += source_param[parameter_type_t::GROUP]) {
        for(unsigned k = 0; k < source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; k += dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
            for(unsigned c = 0; c < source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; c += dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
                for(unsigned r = 0; r < source_param[parameter_type_t::FILTER_HEIGHT]; r += dest_param[parameter_type_t::FILTER_HEIGHT]) {
                    for(unsigned s = 0; s < source_param[parameter_type_t::FILTER_WIDTH]; s += dest_param[parameter_type_t::FILTER_WIDTH]) {
                        unsigned weight_offset = g*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                  *source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                  *source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                               + k*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                  *source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                               + c*source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                               + r*source_param[parameter_type_t::INPUT_WIDTH] + s;
                        m_weight_offsets->emplace_back(weight_offset);
                    }
                }
            }
        }
    }

    // The offset of output data.
    for(unsigned b = 0; b < source_param[parameter_type_t::BATCH_SIZE]; b += dest_param[parameter_type_t::BATCH_SIZE]) {
        for(unsigned g = 0; g < source_param[parameter_type_t::GROUP]; g += dest_param[parameter_type_t::GROUP]) {
            for(unsigned k = 0; k < source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; k += dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
                for(unsigned p = 0; p < source_param[parameter_type_t::OUTPUT_HEIGHT]; p += dest_param[parameter_type_t::OUTPUT_HEIGHT]) {
                    for(unsigned q = 0; q < source_param[parameter_type_t::OUTPUT_WIDTH]; q += dest_param[parameter_type_t::OUTPUT_WIDTH]) {
                        unsigned output_offset = b*source_param[parameter_type_t::GROUP]*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                  *source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                               + g*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                  *source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                               + k*source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                               + p*source_param[parameter_type_t::OUTPUT_WIDTH] + q;
                        m_output_offsets->emplace_back(output_offset);
                    }
                }
            }
        }
    }
}

//std::vector<std::list<unsigned>> scheduler_t::calculate_counter_undefined_stationary_ver2(component_type_t m_destination_type, component_type_t m_source_type, std::list<unsigned> *m_output_offset) {
//}

//std::vector<std::list<unsigned>> scheduler_t::calculate_offset_undefined_stationary_ver2(data_type_t m_data_type, component_type_t m_destination_type, component_type_t m_source_type,
//                                                                                         std::vector<unsigned> *m_offsets, std::vector<unsigned> *m_params, std::vector<unsigned> *m_iteration, 
//                                                                                         std::string m_parameter_order, bool m_last_component) {
//}

std::vector<std::list<unsigned>> scheduler_t::calculate_counter_input_stationary_ver2(component_type_t m_destination_type, component_type_t m_source_type, std::list<unsigned> *m_output_offset) {

    // Initialize and calculate parameter size of accelerator components
    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    
    dest_param = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);

    std::vector<std::list<unsigned>> offset_size;
    offset_size.reserve(data_type_t::NUM_DATA_TYPES);

    std::list<unsigned> input_offset_sizes, weight_offset_sizes, output_offset_sizes;
    unsigned input_offset_size = 0;

    unsigned height_hop = (dest_param[parameter_type_t::FILTER_HEIGHT] == source_param[parameter_type_t::FILTER_HEIGHT] &&
                           dest_param[parameter_type_t::INPUT_HEIGHT] < source_param[parameter_type_t::INPUT_HEIGHT]) ? 
                           dest_param[parameter_type_t::STRIDE]*dest_param[parameter_type_t::OUTPUT_HEIGHT] : 1;
    unsigned width_hop  = (dest_param[parameter_type_t::FILTER_WIDTH] == source_param[parameter_type_t::FILTER_WIDTH] &&
                           dest_param[parameter_type_t::INPUT_WIDTH] < source_param[parameter_type_t::INPUT_WIDTH])  ?
                           dest_param[parameter_type_t::STRIDE]*dest_param[parameter_type_t::OUTPUT_WIDTH] : 1;
    unsigned stride     = dest_param[parameter_type_t::STRIDE];

    // The sequence of input data : B->C->H->W
    for(unsigned b = 0; b < source_param[parameter_type_t::BATCH_SIZE]; b += dest_param[parameter_type_t::BATCH_SIZE]) {
        for(unsigned g = 0; g < source_param[parameter_type_t::GROUP]; g += dest_param[parameter_type_t::GROUP]) {
            for(unsigned c = 0; c < source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; 
                         c += dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
                for(unsigned h = 0; h <= source_param[parameter_type_t::INPUT_HEIGHT] - dest_param[parameter_type_t::INPUT_HEIGHT];) {
                    for(unsigned w = 0; w <= source_param[parameter_type_t::INPUT_WIDTH] - dest_param[parameter_type_t::INPUT_WIDTH];) {
                        input_offset_size++;
                        unsigned weight_offset_size = 0, output_offset_size = 0;
                        // The sequence of weight : C->K->R->S
                        for(unsigned k = 0; k < source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; 
                                     k += dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
                            for(unsigned r = h%stride; r < source_param[parameter_type_t::FILTER_HEIGHT]; 
                                         r += dest_param[parameter_type_t::FILTER_HEIGHT]*dest_param[parameter_type_t::STRIDE]) {
                                if(h >= r && source_param[parameter_type_t::INPUT_HEIGHT] - h >= source_param[parameter_type_t::FILTER_HEIGHT] - r
                                          && (h-r)/height_hop*dest_param[parameter_type_t::OUTPUT_HEIGHT] < source_param[parameter_type_t::OUTPUT_HEIGHT]) {
                                    for(unsigned s = w%stride; s < source_param[parameter_type_t::FILTER_WIDTH]; 
                                                 s += dest_param[parameter_type_t::FILTER_WIDTH]*dest_param[parameter_type_t::STRIDE]) {
                                        if(w >= s && source_param[parameter_type_t::INPUT_WIDTH] - w >= source_param[parameter_type_t::FILTER_WIDTH] - s 
                                                  && (w-s)/width_hop*dest_param[parameter_type_t::OUTPUT_WIDTH] < source_param[parameter_type_t::OUTPUT_WIDTH]) {
                                            weight_offset_size++;
                                            output_offset_size++;
                                        }
                                    }
                                }
                            }
                        }
                        // Calculate the number of reuses of weight and output data.
                        if(weight_offset_size) {weight_offset_sizes.emplace_back(weight_offset_size);}
                        if(output_offset_size) {output_offset_sizes.emplace_back(output_offset_size);}
                        w += width_hop;
                    }
                    h += height_hop;
                }
            }
        }
    }
    input_offset_sizes.emplace_back(input_offset_size);

    offset_size.emplace_back(input_offset_sizes);
    offset_size.emplace_back(weight_offset_sizes);
    offset_size.emplace_back(output_offset_sizes);

    return offset_size;

}

void scheduler_t::calculate_offset_input_stationary_ver2(data_type_t m_data_type, component_type_t m_destination_type, component_type_t m_source_type, 
                                                         std::vector<unsigned> *m_offsets, std::vector<unsigned> *m_params, 
                                                         std::vector<unsigned> *m_iteration, std::vector<unsigned> m_counter,
                                                         std::string m_parameter_order, bool m_last_component) {
    // Initialize parameter size at source and destination components
    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    dest_param = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);

    std::vector<bool> update_params(data_type_t::NUM_DATA_TYPES, false);

    unsigned height_hop = (dest_param[parameter_type_t::FILTER_HEIGHT] == source_param[parameter_type_t::FILTER_HEIGHT] &&
                           dest_param[parameter_type_t::INPUT_HEIGHT] < source_param[parameter_type_t::INPUT_HEIGHT]) ? 
                           dest_param[parameter_type_t::STRIDE]*dest_param[parameter_type_t::OUTPUT_HEIGHT] : 1;
    unsigned width_hop  = (dest_param[parameter_type_t::FILTER_WIDTH] == source_param[parameter_type_t::FILTER_WIDTH] &&
                           dest_param[parameter_type_t::INPUT_WIDTH] < source_param[parameter_type_t::INPUT_WIDTH])  ?
                           dest_param[parameter_type_t::STRIDE]*dest_param[parameter_type_t::OUTPUT_WIDTH] : 1;
    unsigned stride     = dest_param[parameter_type_t::STRIDE];

    // Offset calculation of input data
    if(m_data_type == data_type_t::INPUT) {
        unsigned input_offset = m_params->at(parameter_type_t::BATCH_SIZE)
                               *source_param[parameter_type_t::GROUP]
                               *source_param[parameter_type_t::INPUT_CHANNEL]
                               /source_param[parameter_type_t::GROUP]
                               *source_param[parameter_type_t::INPUT_HEIGHT]
                               *source_param[parameter_type_t::INPUT_WIDTH] + 
                                m_params->at(parameter_type_t::GROUP)
                               *source_param[parameter_type_t::INPUT_CHANNEL]
                               /source_param[parameter_type_t::GROUP]
                               *source_param[parameter_type_t::INPUT_HEIGHT]
                               *source_param[parameter_type_t::INPUT_WIDTH] +
                                m_params->at(parameter_type_t::INPUT_CHANNEL)
                               *source_param[parameter_type_t::INPUT_HEIGHT]
                               *source_param[parameter_type_t::INPUT_WIDTH] + 
                                m_params->at(parameter_type_t::INPUT_HEIGHT)
                               *source_param[parameter_type_t::INPUT_WIDTH] +
                                m_params->at(parameter_type_t::INPUT_WIDTH);

        m_offsets->at(data_type_t::INPUT) = input_offset;
        if(m_last_component) {m_iteration->at(data_type_t::INPUT) += 1;}
    }
    // Offset calculation of weight
    if(m_data_type == data_type_t::WEIGHT) {
        unsigned weight_offset = m_params->at(parameter_type_t::GROUP)
                                *source_param[parameter_type_t::OUTPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::INPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::FILTER_HEIGHT]
                                *source_param[parameter_type_t::FILTER_WIDTH] + 
                                 m_params->at(parameter_type_t::OUTPUT_CHANNEL)
                                *source_param[parameter_type_t::INPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::FILTER_HEIGHT]
                                *source_param[parameter_type_t::FILTER_WIDTH] + 
                                 m_params->at(parameter_type_t::INPUT_CHANNEL)
                                *source_param[parameter_type_t::FILTER_HEIGHT]
                                *source_param[parameter_type_t::FILTER_WIDTH] + 
                                 m_params->at(parameter_type_t::FILTER_HEIGHT)
                                *source_param[parameter_type_t::FILTER_WIDTH] + 
                                 m_params->at(parameter_type_t::FILTER_WIDTH);
        m_offsets->at(data_type_t::WEIGHT) = weight_offset;
        if(m_last_component) {m_iteration->at(data_type_t::WEIGHT) += 1;}
    }
    // Offset calculation of output data
    if(m_data_type == data_type_t::OUTPUT) {

        m_params->at(parameter_type_t::OUTPUT_HEIGHT) = (m_params->at(parameter_type_t::INPUT_HEIGHT)-m_params->at(parameter_type_t::FILTER_HEIGHT))/stride;
        m_params->at(parameter_type_t::OUTPUT_WIDTH) = (m_params->at(parameter_type_t::INPUT_WIDTH)-m_params->at(parameter_type_t::FILTER_WIDTH))/stride;

        unsigned output_offset = m_params->at(parameter_type_t::BATCH_SIZE)
                                *source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::OUTPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::OUTPUT_HEIGHT]
                                *source_param[parameter_type_t::OUTPUT_WIDTH] + 
                                 m_params->at(parameter_type_t::GROUP)
                                *source_param[parameter_type_t::OUTPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::OUTPUT_HEIGHT]
                                *source_param[parameter_type_t::OUTPUT_WIDTH] + 
                                 m_params->at(parameter_type_t::OUTPUT_CHANNEL)
                                *source_param[parameter_type_t::OUTPUT_HEIGHT]
                                *source_param[parameter_type_t::OUTPUT_WIDTH] + 
                                 m_params->at(parameter_type_t::OUTPUT_HEIGHT)
                                *source_param[parameter_type_t::OUTPUT_WIDTH] + 
                                 m_params->at(parameter_type_t::OUTPUT_WIDTH);

        m_offsets->at(data_type_t::OUTPUT) = output_offset;
        if(m_last_component) {m_iteration->at(data_type_t::OUTPUT) += 1;}
    }

    // Reset iterations 
    if(m_iteration->at(data_type_t::WEIGHT) == m_counter[data_type_t::WEIGHT] &&
       m_iteration->at(data_type_t::OUTPUT) == m_counter[data_type_t::OUTPUT]) {

        // Reset input-related parameters
        if(m_iteration->at(data_type_t::INPUT) == m_counter[data_type_t::INPUT]) {
            m_iteration->at(data_type_t::INPUT) = 0;
            m_params->at(parameter_type_t::BATCH_SIZE) = 0;
            m_params->at(parameter_type_t::GROUP) = 0;
            m_params->at(parameter_type_t::INPUT_CHANNEL) = 0;
            m_params->at(parameter_type_t::INPUT_HEIGHT) = 0;
            m_params->at(parameter_type_t::INPUT_WIDTH) = 0;
        }
        // Reset input-irrelevant parameters
        m_iteration->at(data_type_t::WEIGHT) = 0, m_iteration->at(data_type_t::OUTPUT) = 0;
        m_params->at(parameter_type_t::OUTPUT_CHANNEL) = 0;
        m_params->at(parameter_type_t::FILTER_HEIGHT) = m_params->at(parameter_type_t::INPUT_HEIGHT)/stride;
        m_params->at(parameter_type_t::FILTER_WIDTH) = m_params->at(parameter_type_t::INPUT_WIDTH)/stride;
        m_params->at(parameter_type_t::OUTPUT_HEIGHT) = 0;
        m_params->at(parameter_type_t::OUTPUT_WIDTH) = 0;
    }

    // Update input-irrelevant parameters
    if(m_data_type == data_type_t::WEIGHT) {
        m_params->at(parameter_type_t::FILTER_WIDTH) += dest_param[parameter_type_t::FILTER_WIDTH]*stride;
        if(m_params->at(parameter_type_t::FILTER_WIDTH) >= source_param[parameter_type_t::FILTER_WIDTH]) {
            m_params->at(parameter_type_t::FILTER_WIDTH) = m_params->at(parameter_type_t::INPUT_WIDTH)%stride;
            m_params->at(parameter_type_t::FILTER_HEIGHT) += dest_param[parameter_type_t::INPUT_HEIGHT]*stride;
            if(m_params->at(parameter_type_t::FILTER_HEIGHT) >= source_param[parameter_type_t::FILTER_HEIGHT]) {
                m_params->at(parameter_type_t::FILTER_HEIGHT) = m_params->at(parameter_type_t::INPUT_HEIGHT)%stride;
                m_params->at(parameter_type_t::OUTPUT_CHANNEL) += dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP];
                if(m_params->at(parameter_type_t::OUTPUT_CHANNEL) > source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]) {
                    m_params->at(parameter_type_t::OUTPUT_CHANNEL) = 0;
                }
            }
        }
    }
    if(update_params[data_type_t::INPUT]) {
        update_params[data_type_t::INPUT] = false;
        m_params->at(parameter_type_t::INPUT_WIDTH) += width_hop;
        if(m_params->at(parameter_type_t::INPUT_WIDTH) > source_param[parameter_type_t::INPUT_WIDTH] - dest_param[parameter_type_t::INPUT_WIDTH]) {
            m_params->at(parameter_type_t::INPUT_WIDTH) = 0;
            m_params->at(parameter_type_t::INPUT_HEIGHT) += height_hop;
            if(m_params->at(parameter_type_t::INPUT_HEIGHT) > source_param[parameter_type_t::INPUT_HEIGHT] - dest_param[parameter_type_t::INPUT_HEIGHT]) {
                m_params->at(parameter_type_t::INPUT_HEIGHT) = 0;
                m_params->at(parameter_type_t::INPUT_CHANNEL) += dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP];
                if(m_params->at(parameter_type_t::INPUT_CHANNEL) >= source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]) {
                    m_params->at(parameter_type_t::INPUT_CHANNEL) = 0;
                    m_params->at(parameter_type_t::GROUP) += dest_param[parameter_type_t::GROUP];
                    if(m_params->at(parameter_type_t::GROUP) >= source_param[parameter_type_t::GROUP]) {
                        m_params->at(parameter_type_t::GROUP) = 0;
                        m_params->at(parameter_type_t::BATCH_SIZE) += dest_param[parameter_type_t::BATCH_SIZE];
                        if(m_params->at(parameter_type_t::BATCH_SIZE) >= source_param[parameter_type_t::BATCH_SIZE]) {
                            m_params->at(parameter_type_t::BATCH_SIZE) = 0;
                        }
                    }
                }
            }
        }
    }
}

std::vector<std::list<unsigned>> scheduler_t::calculate_counter_weight_stationary_ver2(component_type_t m_destination_type, component_type_t m_source_type, std::list<unsigned> *m_output_offset) {

    // Initialize and calculate tile size (parameter size)
    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    dest_param = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);

    std::vector<std::list<unsigned>> offset_size;
    offset_size.reserve(data_type_t::NUM_DATA_TYPES);

    unsigned weight_offset_size = 0;
    std::list<unsigned> input_offset_sizes, weight_offset_sizes, output_offset_sizes;
    for(unsigned g = 0; g < source_param[parameter_type_t::GROUP]; g += dest_param[parameter_type_t::GROUP]) {
        for(unsigned k = 0; k < source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; 
                     k += dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
            for(unsigned c = 0; c < source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; 
                         c += dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
                for(unsigned r = 0; r < source_param[parameter_type_t::FILTER_HEIGHT]; 
                             r += dest_param[parameter_type_t::FILTER_HEIGHT]) {
                    for(unsigned s = 0; s < source_param[parameter_type_t::FILTER_WIDTH]; 
                                 s += dest_param[parameter_type_t::FILTER_WIDTH]) {
                        // Increment the number of different set of weight
                        weight_offset_size++;

                        unsigned input_offset_size = 0, output_offset_size = 0;
                        for(unsigned b = 0; b < source_param[parameter_type_t::BATCH_SIZE]; 
                                     b += dest_param[parameter_type_t::BATCH_SIZE]) {
                            for(unsigned p = 0; p < source_param[parameter_type_t::OUTPUT_HEIGHT]; 
                                         p += dest_param[parameter_type_t::OUTPUT_HEIGHT]) {
                                for(unsigned q = 0; q < source_param[parameter_type_t::OUTPUT_WIDTH]; 
                                             q += dest_param[parameter_type_t::OUTPUT_WIDTH]) {
                                    // Increase the number of input data and output data sets.
                                    input_offset_size++;
                                    output_offset_size++;
                                }
                            }
                        }
                        // Store the number of input and output data tiles.
                        input_offset_sizes.emplace_back(input_offset_size);
                        output_offset_sizes.emplace_back(output_offset_size);
                    }
                }
            }
        }
    }
    weight_offset_sizes.emplace_back(weight_offset_size);

    offset_size.emplace_back(input_offset_sizes);
    offset_size.emplace_back(weight_offset_sizes);
    offset_size.emplace_back(output_offset_sizes);

    return offset_size;
}

void scheduler_t::calculate_offset_weight_stationary_ver2(data_type_t m_data_type, component_type_t m_destination_type, component_type_t m_source_type, 
                                                          std::vector<unsigned> *m_offsets, std::vector<unsigned> *m_params, 
                                                          std::vector<unsigned> *m_iteration, std::vector<unsigned> m_counter,
                                                          std::string m_parameter_order, bool m_last_component) {

    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    // Calculate the parameter size of destination and source.
    // The parameter sizes are used to calculate tile offset.
    dest_param = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);

    std::vector<bool> update_params(data_type_t::NUM_DATA_TYPES, false);

    unsigned stride = source_param[parameter_type_t::STRIDE];
    // Calculate the input data offset. 
    if(m_data_type == data_type_t::INPUT) {
        // Calculate input offset in B, C, H, W order.
        unsigned input_offset = m_params->at(parameter_type_t::BATCH_SIZE)
                               *source_param[parameter_type_t::GROUP]
                               *source_param[parameter_type_t::INPUT_CHANNEL]
                               /source_param[parameter_type_t::GROUP]
                               *source_param[parameter_type_t::INPUT_HEIGHT]
                               *source_param[parameter_type_t::INPUT_WIDTH] +
                                m_params->at(parameter_type_t::GROUP)
                               *source_param[parameter_type_t::INPUT_CHANNEL]
                               /source_param[parameter_type_t::GROUP]
                               *source_param[parameter_type_t::INPUT_HEIGHT]
                               *source_param[parameter_type_t::INPUT_WIDTH] + 
                                m_params->at(parameter_type_t::INPUT_CHANNEL)
                               *source_param[parameter_type_t::INPUT_HEIGHT]
                               *source_param[parameter_type_t::INPUT_WIDTH] + 
                                m_params->at(parameter_type_t::INPUT_HEIGHT)
                               *source_param[parameter_type_t::INPUT_WIDTH] + 
                                m_params->at(parameter_type_t::INPUT_WIDTH) + 
                                // The offset of input data should be different according to the position of output data.
                                m_params->at(parameter_type_t::OUTPUT_HEIGHT)
                               *source_param[parameter_type_t::INPUT_WIDTH]*stride + 
                                m_params->at(parameter_type_t::OUTPUT_WIDTH)*stride;
        m_offsets->at(data_type_t::INPUT) = input_offset;
        if(m_last_component) {m_iteration->at(data_type_t::INPUT) += 1;}
    }
    // Calculate the weight offset.
    if(m_data_type == data_type_t::WEIGHT) {
        unsigned weight_offset = m_params->at(parameter_type_t::GROUP)
                                *source_param[parameter_type_t::OUTPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::INPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::FILTER_HEIGHT]
                                *source_param[parameter_type_t::FILTER_WIDTH] + 
                                 m_params->at(parameter_type_t::OUTPUT_CHANNEL)
                                *source_param[parameter_type_t::INPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::FILTER_HEIGHT]
                                *source_param[parameter_type_t::FILTER_WIDTH] +
                                 m_params->at(parameter_type_t::INPUT_CHANNEL)
                                *source_param[parameter_type_t::FILTER_HEIGHT]
                                *source_param[parameter_type_t::FILTER_WIDTH] + 
                                 m_params->at(parameter_type_t::FILTER_HEIGHT)
                                *source_param[parameter_type_t::FILTER_WIDTH] + 
                                 m_params->at(parameter_type_t::FILTER_WIDTH);
        m_offsets->at(data_type_t::WEIGHT) = weight_offset;
        if(m_last_component) {m_iteration->at(data_type_t::WEIGHT) += 1;}
    }
    // Calculate the output offset.
    if(m_data_type == data_type_t::OUTPUT) {
        unsigned output_offset = m_params->at(parameter_type_t::BATCH_SIZE)
                                *source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::OUTPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::OUTPUT_HEIGHT]
                                *source_param[parameter_type_t::OUTPUT_WIDTH] + 
                                 m_params->at(parameter_type_t::GROUP)
                                *source_param[parameter_type_t::OUTPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::OUTPUT_HEIGHT]
                                *source_param[parameter_type_t::OUTPUT_WIDTH] + 
                                 m_params->at(parameter_type_t::OUTPUT_CHANNEL)
                                *source_param[parameter_type_t::OUTPUT_HEIGHT]
                                *source_param[parameter_type_t::OUTPUT_WIDTH] +
                                 m_params->at(parameter_type_t::OUTPUT_HEIGHT)
                                *source_param[parameter_type_t::OUTPUT_WIDTH] + 
                                 m_params->at(parameter_type_t::OUTPUT_WIDTH);
        m_offsets->at(data_type_t::OUTPUT) = output_offset;
        if(m_last_component) {m_iteration->at(data_type_t::OUTPUT) += 1;}
    }

    // Reset iterations of DNN parameters.
    //std::cout << m_iteration->at(data_type_t::INPUT) << " " << m_counter[data_type_t::INPUT] << " " << m_iteration->at(data_type_t::OUTPUT) << " " << m_counter[data_type_t::OUTPUT] << std::endl;
    if(m_iteration->at(data_type_t::INPUT) == m_counter[data_type_t::INPUT] &&
       m_iteration->at(data_type_t::OUTPUT) == m_counter[data_type_t::OUTPUT]) {
        // Reset iterations of weight-related parameters.
        if(m_iteration->at(data_type_t::WEIGHT) == m_counter[data_type_t::WEIGHT]) {
            m_iteration->at(data_type_t::WEIGHT) = 0;
            m_params->at(parameter_type_t::GROUP) = 0;
            m_params->at(parameter_type_t::OUTPUT_CHANNEL) = 0;
            m_params->at(parameter_type_t::INPUT_CHANNEL) = 0;
            m_params->at(parameter_type_t::FILTER_HEIGHT) = 0;
            m_params->at(parameter_type_t::FILTER_WIDTH) = 0;
        } 
        // Reset iterations of weight-irrelevant parameters.
        m_iteration->at(data_type_t::INPUT) = 0, m_iteration->at(data_type_t::OUTPUT) = 0;
        update_params[data_type_t::WEIGHT] = true;
        m_params->at(parameter_type_t::BATCH_SIZE) = 0;
        m_params->at(parameter_type_t::OUTPUT_HEIGHT) = 0;
        m_params->at(parameter_type_t::OUTPUT_WIDTH) = 0;
    }

    // Update weight-irrelevant parameters.
    if(m_data_type == data_type_t::INPUT) {
        m_params->at(parameter_type_t::OUTPUT_WIDTH) += dest_param[parameter_type_t::OUTPUT_WIDTH];
        if(m_params->at(parameter_type_t::OUTPUT_WIDTH) >= source_param[parameter_type_t::OUTPUT_WIDTH]) {
            m_params->at(parameter_type_t::OUTPUT_WIDTH) = 0;
            m_params->at(parameter_type_t::OUTPUT_HEIGHT) += dest_param[parameter_type_t::OUTPUT_HEIGHT];
            if(m_params->at(parameter_type_t::OUTPUT_HEIGHT) >= source_param[parameter_type_t::OUTPUT_HEIGHT]) {
                m_params->at(parameter_type_t::OUTPUT_HEIGHT) = 0;
                m_params->at(parameter_type_t::BATCH_SIZE) += dest_param[parameter_type_t::BATCH_SIZE];
                if(m_params->at(parameter_type_t::BATCH_SIZE) >= source_param[parameter_type_t::BATCH_SIZE]) {
                    m_params->at(parameter_type_t::BATCH_SIZE) = 0;
                }
            }
        }
    }

    // Update weight-related parameters.
    if(update_params[data_type_t::WEIGHT]) {
        update_params[data_type_t::WEIGHT] = false;
        // K->C->R->S
        m_params->at(parameter_type_t::FILTER_WIDTH) += dest_param[parameter_type_t::FILTER_WIDTH];
        m_params->at(parameter_type_t::INPUT_WIDTH) += dest_param[parameter_type_t::FILTER_WIDTH];
        if(m_params->at(parameter_type_t::FILTER_WIDTH) >= source_param[parameter_type_t::FILTER_WIDTH]) {
            m_params->at(parameter_type_t::FILTER_WIDTH) = 0;
            m_params->at(parameter_type_t::INPUT_HEIGHT) += dest_param[parameter_type_t::FILTER_HEIGHT];
            m_params->at(parameter_type_t::INPUT_WIDTH) = 0;
            m_params->at(parameter_type_t::FILTER_HEIGHT) += dest_param[parameter_type_t::FILTER_HEIGHT];
            if(m_params->at(parameter_type_t::FILTER_HEIGHT) >= source_param[parameter_type_t::FILTER_HEIGHT] ) {
                m_params->at(parameter_type_t::FILTER_HEIGHT) = 0;
                m_params->at(parameter_type_t::INPUT_HEIGHT) = 0, m_params->at(parameter_type_t::INPUT_WIDTH) = 0;
                m_params->at(parameter_type_t::INPUT_CHANNEL) += dest_param[parameter_type_t::INPUT_CHANNEL]
                                                                /dest_param[parameter_type_t::GROUP];
                if(m_params->at(parameter_type_t::INPUT_CHANNEL) >= source_param[parameter_type_t::INPUT_CHANNEL]
                                                                   /source_param[parameter_type_t::GROUP]) {
                    m_params->at(parameter_type_t::INPUT_CHANNEL) = 0;
                    m_params->at(parameter_type_t::OUTPUT_CHANNEL) += dest_param[parameter_type_t::OUTPUT_CHANNEL]
                                                                     /dest_param[parameter_type_t::GROUP];
                    if(m_params->at(parameter_type_t::OUTPUT_CHANNEL) >= source_param[parameter_type_t::OUTPUT_CHANNEL]
                                                                        /source_param[parameter_type_t::GROUP]) {
                        m_params->at(parameter_type_t::OUTPUT_CHANNEL) = 0;
                        m_params->at(parameter_type_t::GROUP) += dest_param[parameter_type_t::GROUP];
                        if(m_params->at(parameter_type_t::GROUP) >= source_param[parameter_type_t::GROUP]) {
                            m_params->at(parameter_type_t::GROUP) = 0;
                        }
                    }
                }
            }
        }
    }
}

std::vector<std::list<unsigned>> scheduler_t::calculate_counter_output_stationary_ver2(component_type_t m_destination_type, component_type_t m_source_type, std::list<unsigned> *m_output_offset) {
    // Initialize and calculate parameter sizes.
    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    dest_param = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);

    std::vector<std::list<unsigned>> offset_size;
    offset_size.reserve(data_type_t::NUM_DATA_TYPES);

    std::list<unsigned> input_offset_sizes, weight_offset_sizes, output_offset_sizes;
    unsigned output_offset_size = 0;
    for(unsigned b = 0; b < source_param[parameter_type_t::BATCH_SIZE]; b += dest_param[parameter_type_t::BATCH_SIZE]) {
        for(unsigned g = 0; g < source_param[parameter_type_t::GROUP]; g += dest_param[parameter_type_t::GROUP]) {
            for(unsigned k = 0; k < source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; 
                         k += dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
                for(unsigned p = 0; p < source_param[parameter_type_t::OUTPUT_HEIGHT]; p += dest_param[parameter_type_t::OUTPUT_HEIGHT]) {
                    for(unsigned q = 0; q < source_param[parameter_type_t::OUTPUT_WIDTH]; q += dest_param[parameter_type_t::OUTPUT_WIDTH]) {
                        // Calculate the offset of output data
                        // The sequence of output data : B->K->P->Q
                        output_offset_size++;

                        unsigned input_offset_size = 0, weight_offset_size = 0;
                        for(unsigned c = 0; c < source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; 
                                     c += dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]) {
                            for(unsigned r = 0; r < source_param[parameter_type_t::FILTER_HEIGHT]; r += dest_param[parameter_type_t::FILTER_HEIGHT]) {
                                for(unsigned s = 0; s < source_param[parameter_type_t::FILTER_WIDTH]; s += dest_param[parameter_type_t::FILTER_WIDTH]) {
                                    input_offset_size++;
                                    weight_offset_size++;
                                }
                            }
                        }
                        input_offset_sizes.emplace_back(input_offset_size);
                        weight_offset_sizes.emplace_back(weight_offset_size);
                    }
                }
            }
        }
    }

    output_offset_sizes.emplace_back(output_offset_size);

    offset_size.emplace_back(input_offset_sizes);
    offset_size.emplace_back(weight_offset_sizes);
    offset_size.emplace_back(output_offset_sizes);

    return offset_size;

}

void scheduler_t::calculate_offset_output_stationary_ver2(data_type_t m_data_type, component_type_t m_destination_type, component_type_t m_source_type, 
                                                          std::vector<unsigned> *m_offsets, std::vector<unsigned> *m_params, 
                                                          std::vector<unsigned> *m_iteration, std::vector<unsigned> m_counter,
                                                          std::string m_parameter_order, bool m_last_component) {

    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);


    // Calculate parameter size of destination and source components.
    dest_param = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);

    std::vector<bool> update_params(data_type_t::NUM_DATA_TYPES, false);

    unsigned stride = source_param[parameter_type_t::STRIDE];

    // Calculate input data offset
    if(m_data_type == data_type_t::INPUT) {
        unsigned input_offset = m_params->at(parameter_type_t::BATCH_SIZE)
                               *source_param[parameter_type_t::GROUP]
                               *source_param[parameter_type_t::INPUT_CHANNEL]
                               /source_param[parameter_type_t::GROUP]
                               *source_param[parameter_type_t::INPUT_HEIGHT]
                               *source_param[parameter_type_t::INPUT_WIDTH] + 
                                m_params->at(parameter_type_t::GROUP)
                               *source_param[parameter_type_t::INPUT_CHANNEL]
                               /source_param[parameter_type_t::GROUP]
                               *source_param[parameter_type_t::INPUT_HEIGHT]
                               *source_param[parameter_type_t::INPUT_WIDTH] + 
                                m_params->at(parameter_type_t::INPUT_CHANNEL)
                               *source_param[parameter_type_t::INPUT_HEIGHT]
                               *source_param[parameter_type_t::INPUT_WIDTH] +
                                m_params->at(parameter_type_t::INPUT_HEIGHT)
                               *source_param[parameter_type_t::INPUT_WIDTH] + 
                                m_params->at(parameter_type_t::INPUT_WIDTH) +
                                m_params->at(parameter_type_t::FILTER_HEIGHT)
                               *source_param[parameter_type_t::INPUT_WIDTH] + 
                                m_params->at(parameter_type_t::FILTER_WIDTH);

        m_offsets->at(data_type_t::INPUT) = input_offset;
        if(m_last_component) {m_iteration->at(data_type_t::INPUT) += 1;}
    }
    // Calculate weight offset
    if(m_data_type == data_type_t::WEIGHT) {
        unsigned weight_offset = m_params->at(parameter_type_t::GROUP)
                                *source_param[parameter_type_t::OUTPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::INPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::FILTER_HEIGHT]
                                *source_param[parameter_type_t::FILTER_WIDTH] + 
                                 m_params->at(parameter_type_t::OUTPUT_CHANNEL)
                                *source_param[parameter_type_t::INPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::FILTER_HEIGHT]
                                *source_param[parameter_type_t::FILTER_WIDTH] + 
                                 m_params->at(parameter_type_t::INPUT_CHANNEL)
                                *source_param[parameter_type_t::FILTER_HEIGHT]
                                *source_param[parameter_type_t::FILTER_WIDTH] +
                                 m_params->at(parameter_type_t::FILTER_HEIGHT) 
                                *source_param[parameter_type_t::FILTER_WIDTH] + 
                                 m_params->at(parameter_type_t::FILTER_WIDTH);

        m_offsets->at(data_type_t::WEIGHT) = weight_offset;
        if(m_last_component) {m_iteration->at(data_type_t::WEIGHT) += 1;}
    }

    // Calculate output offset
    if(m_data_type == data_type_t::OUTPUT) {
        unsigned output_offset = m_params->at(parameter_type_t::BATCH_SIZE)
                                *source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::OUTPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::OUTPUT_HEIGHT]
                                *source_param[parameter_type_t::OUTPUT_WIDTH] + 
                                 m_params->at(parameter_type_t::GROUP)
                                *source_param[parameter_type_t::OUTPUT_CHANNEL]
                                /source_param[parameter_type_t::GROUP]
                                *source_param[parameter_type_t::OUTPUT_HEIGHT]
                                *source_param[parameter_type_t::OUTPUT_WIDTH] + 
                                 m_params->at(parameter_type_t::OUTPUT_CHANNEL)
                                *source_param[parameter_type_t::OUTPUT_HEIGHT]
                                *source_param[parameter_type_t::OUTPUT_WIDTH] +
                                 m_params->at(parameter_type_t::OUTPUT_HEIGHT)
                                *source_param[parameter_type_t::OUTPUT_WIDTH] + 
                                 m_params->at(parameter_type_t::OUTPUT_WIDTH);

        m_offsets->at(data_type_t::OUTPUT) = output_offset;
        if(m_last_component) {m_iteration->at(data_type_t::OUTPUT) += 1;}
    }

    // Reset iteration value
    if(m_iteration->at(data_type_t::INPUT) == m_counter[data_type_t::INPUT] &&
       m_iteration->at(data_type_t::WEIGHT) == m_counter[data_type_t::WEIGHT]) {

        // Reset output-related parameters
        if(m_iteration->at(data_type_t::OUTPUT) == m_counter[data_type_t::OUTPUT]) {
            m_iteration->at(data_type_t::OUTPUT) = 0;
            m_params->at(parameter_type_t::BATCH_SIZE) = 0;
            m_params->at(parameter_type_t::GROUP) = 0;
            m_params->at(parameter_type_t::OUTPUT_CHANNEL) = 0;
            m_params->at(parameter_type_t::OUTPUT_HEIGHT) = 0;
            m_params->at(parameter_type_t::OUTPUT_WIDTH) = 0;
        }
        // Reset output-irrelevant parameters
        m_iteration->at(data_type_t::INPUT) = 0, m_iteration->at(data_type_t::WEIGHT) = 0;
        update_params[data_type_t::OUTPUT] = true;
        m_params->at(parameter_type_t::INPUT_CHANNEL) = 0;
        m_params->at(parameter_type_t::FILTER_HEIGHT) = 0;
        m_params->at(parameter_type_t::FILTER_WIDTH) = 0;
    }

    // Update output-irrelevant parameters
    /* Need to update : Not only input data */
    if(m_data_type == data_type_t::INPUT) {
        m_params->at(parameter_type_t::FILTER_WIDTH) += dest_param[parameter_type_t::FILTER_WIDTH];
        if(m_params->at(parameter_type_t::FILTER_WIDTH) >= source_param[parameter_type_t::FILTER_WIDTH]) {
            m_params->at(parameter_type_t::FILTER_WIDTH) = 0;
            m_params->at(parameter_type_t::FILTER_HEIGHT) += dest_param[parameter_type_t::FILTER_HEIGHT];
            if(m_params->at(parameter_type_t::FILTER_HEIGHT) >= source_param[parameter_type_t::FILTER_HEIGHT]) {
                m_params->at(parameter_type_t::FILTER_HEIGHT) = 0;
                m_params->at(parameter_type_t::INPUT_CHANNEL) += dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP];
                if(m_params->at(parameter_type_t::INPUT_CHANNEL) >= source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]) {
                    m_params->at(parameter_type_t::INPUT_CHANNEL) = 0;
                }
            }
        }
    }

    // Update output-related parameters
    if(update_params[data_type_t::OUTPUT]) {
        update_params[data_type_t::OUTPUT] = false;
        // B<-K<-P<-Q (Reverse order)
        m_params->at(parameter_type_t::OUTPUT_WIDTH) += dest_param[parameter_type_t::OUTPUT_WIDTH];
        m_params->at(parameter_type_t::INPUT_WIDTH) += dest_param[parameter_type_t::OUTPUT_WIDTH]*stride;
        if(m_params->at(parameter_type_t::OUTPUT_WIDTH) >= source_param[parameter_type_t::OUTPUT_WIDTH]) {
            m_params->at(parameter_type_t::OUTPUT_WIDTH) = 0;
            m_params->at(parameter_type_t::INPUT_WIDTH) = 0;
            m_params->at(parameter_type_t::OUTPUT_HEIGHT) += dest_param[parameter_type_t::OUTPUT_HEIGHT];
            m_params->at(parameter_type_t::INPUT_HEIGHT) += dest_param[parameter_type_t::OUTPUT_HEIGHT]*stride;
            if(m_params->at(parameter_type_t::OUTPUT_HEIGHT) >= source_param[parameter_type_t::OUTPUT_HEIGHT]) {
                m_params->at(parameter_type_t::OUTPUT_HEIGHT) = 0;
                m_params->at(parameter_type_t::INPUT_HEIGHT) = 0;
                m_params->at(parameter_type_t::OUTPUT_CHANNEL) += dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP];
                if(m_params->at(parameter_type_t::OUTPUT_CHANNEL) >= source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]) {
                    m_params->at(parameter_type_t::OUTPUT_CHANNEL) = 0;
                    m_params->at(parameter_type_t::GROUP) += dest_param[parameter_type_t::GROUP];
                    if(m_params->at(parameter_type_t::GROUP) >= source_param[parameter_type_t::GROUP]) {
                        m_params->at(parameter_type_t::GROUP) = 0;
                        m_params->at(parameter_type_t::BATCH_SIZE) += dest_param[parameter_type_t::BATCH_SIZE];
                        if(m_params->at(parameter_type_t::BATCH_SIZE) >= source_param[parameter_type_t::BATCH_SIZE]) {
                            m_params->at(parameter_type_t::BATCH_SIZE) = 0;
                        }
                    }
                }
            }
        }
    }
}

#ifdef FUNCTIONAL

// Transfer input data.
void scheduler_t::input_data_load(data_t *m_dest, data_t *m_source, unsigned m_dest_offset, unsigned m_source_offset, 
                                                    component_type_t m_destination_type, component_type_t m_source_type) {
    
    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    dest_param = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);

    num_zeros[data_type_t::INPUT] = 0;

    for(unsigned b = 0; b < dest_param[parameter_type_t::BATCH_SIZE]; b++) {
        for(unsigned g = 0; g < dest_param[parameter_type_t::GROUP]; g++) {
            for(unsigned c = 0; c < dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]; c++) {
                for(unsigned h = 0; h < dest_param[parameter_type_t::INPUT_HEIGHT]; h++) {
                    for(unsigned w = 0; w < dest_param[parameter_type_t::INPUT_WIDTH]; w++) {

                        unsigned source_index = m_source_offset + b*source_param[parameter_type_t::GROUP]*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                              + g*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                              + c*source_param[parameter_type_t::INPUT_HEIGHT]*source_param[parameter_type_t::INPUT_WIDTH]
                                              + h*source_param[parameter_type_t::INPUT_WIDTH] + w;                                         

                        unsigned dest_index   = m_dest_offset + b*dest_param[parameter_type_t::GROUP]*dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]
                                                 *dest_param[parameter_type_t::INPUT_HEIGHT]*dest_param[parameter_type_t::INPUT_WIDTH]
                                              + g*dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]
                                                 *dest_param[parameter_type_t::INPUT_HEIGHT]*dest_param[parameter_type_t::INPUT_WIDTH]
                                              + c*dest_param[parameter_type_t::INPUT_HEIGHT]*dest_param[parameter_type_t::INPUT_WIDTH]
                                              + h*dest_param[parameter_type_t::INPUT_WIDTH] + w;

                        if(m_source[source_index] == 0.0) {num_zeros[data_type_t::INPUT]++;}

                        m_dest[dest_index] = m_source[source_index];
                    }
                }
            }
        }
    }
}

// Transfer weight.
void scheduler_t::weight_data_load(data_t *m_dest, data_t *m_source, unsigned m_dest_offset, unsigned m_source_offset, 
                                   component_type_t m_destination_type, component_type_t m_source_type) {

    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    dest_param = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);

    num_zeros[data_type_t::WEIGHT] = 0;

    for(unsigned g = 0; g < dest_param[parameter_type_t::GROUP]; g++) {
        for(unsigned k = 0; k < dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]; k++) {
            for(unsigned c = 0; c < dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]; c++) {
                for(unsigned r = 0; r < dest_param[parameter_type_t::FILTER_HEIGHT]; r++) {
                    for(unsigned s = 0; s < dest_param[parameter_type_t::FILTER_WIDTH]; s++) {
                        unsigned source_index = m_source_offset + g*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                              + k*source_param[parameter_type_t::INPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                              + c*source_param[parameter_type_t::FILTER_HEIGHT]*source_param[parameter_type_t::FILTER_WIDTH]
                                              + r*source_param[parameter_type_t::FILTER_WIDTH] + s;

                        unsigned dest_index   = m_dest_offset + g*dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]
                                                 *dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]
                                                 *dest_param[parameter_type_t::FILTER_HEIGHT]*dest_param[parameter_type_t::FILTER_WIDTH]
                                              + k*dest_param[parameter_type_t::INPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]
                                                 *dest_param[parameter_type_t::FILTER_HEIGHT]*dest_param[parameter_type_t::FILTER_WIDTH]
                                              + c*dest_param[parameter_type_t::FILTER_HEIGHT]*dest_param[parameter_type_t::FILTER_WIDTH]
                                              + r*dest_param[parameter_type_t::FILTER_WIDTH] + s;

                        if(m_source[source_index] == 0.0) {num_zeros[data_type_t::WEIGHT]++;}
                                              
                        m_dest[dest_index] = m_source[source_index];
                    }
                }
            }
        }
    }
}

// Transfer output data.
void scheduler_t::output_data_load(data_t *m_dest, data_t *m_source, unsigned m_dest_offset, unsigned m_source_offset, 
                                   component_type_t m_destination_type, component_type_t m_source_type) {

    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    dest_param = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);

    for(unsigned b = 0; b < dest_param[parameter_type_t::BATCH_SIZE]; b++) {
        for(unsigned g = 0; g < dest_param[parameter_type_t::GROUP]; g++) {
            for(unsigned k = 0; k < dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]; k++) {
                for(unsigned p = 0; p < dest_param[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                    for(unsigned q = 0; q < dest_param[parameter_type_t::OUTPUT_WIDTH]; q++) {
                        unsigned source_index = m_source_offset + b*source_param[parameter_type_t::GROUP]*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                              + g*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                              + k*source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                              + p*source_param[parameter_type_t::OUTPUT_WIDTH] + q;

                        unsigned dest_index   = m_dest_offset + b*dest_param[parameter_type_t::GROUP]*dest_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *dest_param[parameter_type_t::OUTPUT_HEIGHT]*dest_param[parameter_type_t::OUTPUT_WIDTH]
                                              * g*dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]
                                                 *dest_param[parameter_type_t::OUTPUT_HEIGHT]*dest_param[parameter_type_t::OUTPUT_WIDTH]
                                              + k*dest_param[parameter_type_t::OUTPUT_HEIGHT]*dest_param[parameter_type_t::OUTPUT_WIDTH]
                                              + p*dest_param[parameter_type_t::OUTPUT_WIDTH] + q;

                        m_dest[dest_index] = m_source[source_index];
                    }
                }
            }
        }
    }
}

void scheduler_t::output_data_store(data_t *m_dest, data_t *m_source, unsigned m_dest_offset, unsigned m_source_offset,
                                    component_type_t m_destination_type, component_type_t m_source_type) {

    std::vector<unsigned> dest_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> source_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    dest_param   = mapping_table->calculate_parameter_size(m_destination_type);
    source_param = mapping_table->calculate_parameter_size(m_source_type);
    
    for(unsigned b = 0; b < source_param[parameter_type_t::BATCH_SIZE]; b++) {
        for(unsigned g = 0; g < source_param[parameter_type_t::GROUP]; g++) {
            for(unsigned k = 0; k < source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]; k++) {
                for(unsigned p = 0; p < source_param[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                    for(unsigned q = 0; q < source_param[parameter_type_t::OUTPUT_WIDTH]; q++) {
                        unsigned source_index = m_source_offset + b*source_param[parameter_type_t::GROUP]*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                              + g*source_param[parameter_type_t::OUTPUT_CHANNEL]/source_param[parameter_type_t::GROUP]
                                                 *source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                              + k*source_param[parameter_type_t::OUTPUT_HEIGHT]*source_param[parameter_type_t::OUTPUT_WIDTH]
                                              + p*source_param[parameter_type_t::OUTPUT_WIDTH] + q;

                        unsigned dest_index   = m_dest_offset + b*dest_param[parameter_type_t::GROUP]*dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]
                                                 *dest_param[parameter_type_t::OUTPUT_HEIGHT]*dest_param[parameter_type_t::OUTPUT_WIDTH]
                                              + g*dest_param[parameter_type_t::OUTPUT_CHANNEL]/dest_param[parameter_type_t::GROUP]
                                                 *dest_param[parameter_type_t::OUTPUT_HEIGHT]*dest_param[parameter_type_t::OUTPUT_WIDTH]
                                              + k*dest_param[parameter_type_t::OUTPUT_HEIGHT]*dest_param[parameter_type_t::OUTPUT_WIDTH]
                                              + p*dest_param[parameter_type_t::OUTPUT_WIDTH] + q;

                        m_dest[dest_index] = m_source[source_index];
                    }
                }
            }
        }
    }
}

#endif

