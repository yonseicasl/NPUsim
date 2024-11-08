#include <string>
#include <cassert>
#include <cstring>
#include <cmath>
#include "pe.h"

pe_t::pe_t(section_config_t m_section_config) :
    input_data_mac(NULL),
    weight_mac(NULL),
    output_data_mac(NULL),
    input_data_lb(NULL),
    weight_lb(NULL),
    output_data_lb(NULL),
    input_size(0),
    weight_size(0),
    output_size(0),
    index(0),
    duplicated_input_mac(0),
    duplicated_input_lb(0),
    /* Unit stats */
    u_computation_cycle(0.0),
    u_computation_energy(0.0),
    u_transfer_cycle(0.0),
    u_transfer_energy(0.0),
    u_dynamic_power_mac(0.0),
    u_static_power_mac(0.0),
    /* Accelerator stats */
    num_computation(0),
    computation_cycle(0.0),
    computation_energy(0.0),
    utilization_mac(0.0),
    write_back_cycle_mac(0.0),
    write_back_cycle_lb(0.0),
    overlapped_transfer_cycle(0.0),
    pe_array(NULL),
    stationary_type_mac(stationary_type_t::UNDEFINED_STATIONARY),
    stationary_type_local_buffer(stationary_type_t::UNDEFINED_STATIONARY),
    parameter_order("kbpqcrs"),
    memory_type(memory_type_t::SEPARATE),
    mac_type(mac_type_t::UNDEFINED_MAC),
    num_macs(1),
    mac_width(1),
    num_active_macs(1),
    frequency(0.0),
    bandwidth(0.0),
    bitwidth(0),
    input_index(0),
    weight_index(0),
    output_index(0),
    input_flush_counter(0),
    weight_flush_counter(0),
    output_flush_counter(0),
    idle(false) {

    init(m_section_config);

}

pe_t::~pe_t() {

    // Free memory space at MAC unit.
    delete [] input_data_mac;
    delete [] weight_mac;
    delete [] output_data_mac;

    // Free memory space at local buffer.
    delete [] input_data_lb;
    delete [] weight_lb;
    delete [] output_data_lb;

}

// Initialize the PE.
void pe_t::init(section_config_t m_section_config) {

    /* Initialize PE specifications */

    // Initialize the size of input, weight, and output buffer in Byte.
    m_section_config.get_setting("input_size", &input_size);
    m_section_config.get_setting("weight_size", &weight_size);
    m_section_config.get_setting("output_size", &output_size);

    // Initialize the number of elements in each buffer.
    unsigned num_input = input_size/sizeof(data_t);
    unsigned num_weight = weight_size/sizeof(data_t);
    unsigned num_output = output_size/sizeof(data_t);

    // Allocate the memory space for local buffer.
    input_data_lb  = new data_t[num_input]();
    weight_lb      = new data_t[num_weight]();
    output_data_lb = new data_t[num_output](); 
    memset(input_data_lb, 1.0, num_input * sizeof(data_t));
    memset(weight_lb, 1.0, num_weight * sizeof(data_t));
    memset(output_data_lb, 1.0, num_output * sizeof(data_t));

    // Initialize the size of MAC register
    m_section_config.get_setting("number_of_macs", &num_macs);
    m_section_config.get_setting("mac_width", &mac_width);

    input_data_mac  = new data_t[mac_width*num_macs]();
    weight_mac      = new data_t[mac_width*num_macs]();
    output_data_mac = new data_t[num_macs]();

    memset(input_data_mac, 1.0, mac_width*num_macs*sizeof(data_t));
    memset(weight_mac, 1.0, mac_width*num_macs*sizeof(data_t));
    memset(output_data_mac, 1.0, num_macs*sizeof(data_t));

    // Initialize the frequency (MHz) and bandwidth (GB/sec)
    m_section_config.get_setting("frequency", &frequency);
    m_section_config.get_setting("bandwidth", &bandwidth);
    bitwidth = 8*bandwidth/frequency;
    m_section_config.get_setting("bitwidth", &bitwidth);

    bypass.reserve(data_type_t::NUM_DATA_TYPES);
    bypass.assign(data_type_t::NUM_DATA_TYPES, 0);
    m_section_config.get_vector_setting("bypass", &bypass);

    // Initialize line size and mask bits of MAC register and local buffer
    line_size_mac.reserve(data_type_t::NUM_DATA_TYPES);
    line_size_mac.assign(data_type_t::NUM_DATA_TYPES, 8);

    mask_bits_mac.reserve(data_type_t::NUM_DATA_TYPES);
    mask_bits_mac.assign(data_type_t::NUM_DATA_TYPES, 0);

    line_size_lb.reserve(data_type_t::NUM_DATA_TYPES);
    line_size_lb.assign(data_type_t::NUM_DATA_TYPES, 8);

    mask_bits_lb.reserve(data_type_t::NUM_DATA_TYPES);
    mask_bits_lb.assign(data_type_t::NUM_DATA_TYPES, 0);

    m_section_config.get_vector_setting("mac_line_size", &line_size_mac);
    m_section_config.get_vector_setting("lb_line_size", &line_size_lb);

    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        while(line_size_mac[i] > 8) {
            line_size_mac[i] /= 2;
            mask_bits_mac[i]++;
        }

        while(line_size_lb[i] > 8) {
            line_size_lb[i] /= 2;
            mask_bits_lb[i]++;
        }
    }

    m_section_config.get_vector_setting("mac_line_size", &line_size_mac);
    m_section_config.get_vector_setting("lb_line_size", &line_size_lb); 

    // Initialize the stationary type.
    // Initialize the stationary type between MAC and local buffer.
    std::string stationary_str;
    if(m_section_config.get_setting("mac_stationary", &stationary_str)) {
        stationary_type_mac = (stationary_type_t)get_type(stationary_type_str, stationary_str);
    }
    
    // Initialize the order of parameters.
    if(stationary_type_mac == stationary_type_t::INPUT_STATIONARY) {
        parameter_order = "BCPQKRS";
    } else if(stationary_type_mac == stationary_type_t::WEIGHT_STATIONARY) {
        parameter_order = "KCRSBPQ";
    } else if(stationary_type_mac == stationary_type_t::OUTPUT_STATIONARY) {
        parameter_order = "BKPQKRS";
    }
    m_section_config.get_setting("pe_parameter_order", &parameter_order);

    // Initialize the stationary type between PE array and Global buffer.
    if(m_section_config.get_setting("pe_stationary", &stationary_str)) {
        stationary_type_local_buffer = (stationary_type_t)get_type(stationary_type_str, stationary_str);
    }

    std::string memory_str;
    if(m_section_config.get_setting("memory_type", &memory_str)) {
        memory_type = (memory_type_t)get_type(memory_type_str, memory_str);
    }

    skip_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);

    // Initialize the tile size of MAC unit.
    tile_size_mac.reserve(data_type_t::NUM_DATA_TYPES);
    tile_size_mac.assign(data_type_t::NUM_DATA_TYPES, 1);

    // Initialize the tile size of Local buffer.
    tile_size_lb.reserve(data_type_t::NUM_DATA_TYPES);
    tile_size_lb.assign(data_type_t::NUM_DATA_TYPES, 1);

    offsets_mac.reserve(data_type_t::NUM_DATA_TYPES);
    offsets_mac.assign(data_type_t::NUM_DATA_TYPES, 0);

    offsets_lb.reserve(data_type_t::NUM_DATA_TYPES);
    offsets_lb.assign(data_type_t::NUM_DATA_TYPES, 0);

    /* Initialize PE signals */

    // Exist data in MAC unit.
    exist_data_mac.reserve(data_type_t::NUM_DATA_TYPES);
    exist_data_mac.assign(data_type_t::NUM_DATA_TYPES, false);    

    // Exist data in Local buffer.
    exist_data_lb.reserve(data_type_t::NUM_DATA_TYPES);
    exist_data_lb.assign(data_type_t::NUM_DATA_TYPES, false);

    // Exist request to Local buffer (from MAC unit).
    request_to_lb.reserve(data_type_t::NUM_DATA_TYPES);
    request_to_lb.assign(data_type_t::NUM_DATA_TYPES, true);

    // Exist request to PE array (from Local buffer).
    request_to_pe_array.reserve(data_type_t::NUM_DATA_TYPES);
    request_to_pe_array.assign(data_type_t::NUM_DATA_TYPES, false);

    /* Initialize PE unit costs */

    // Initialize the MAC unit cycle and energy.
    m_section_config.get_setting("computation_cycle", &u_computation_cycle);
    m_section_config.get_setting("computation_energy", &u_computation_energy);

    // Initialize the transfer cycle/energy between PE and MAC unit and local buffer.
    m_section_config.get_setting("transfer_cycle_pe", &u_transfer_cycle);
    m_section_config.get_setting("transfer_energy_pe", &u_transfer_energy);

    m_section_config.get_setting("mac_static_power", &u_static_power_mac);
    m_section_config.get_setting("mac_dynamic_power", &u_dynamic_power_mac);

    // Initialize the unit cycle and energy of MAC unit
    u_read_cycle_mac.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_cycle_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("mac_read_cycle", &u_read_cycle_mac);

    u_read_energy_mac.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_energy_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("mac_read_energy", &u_read_energy_mac);

    u_write_cycle_mac.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_cycle_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("mac_write_cycle", &u_write_cycle_mac);

    u_write_energy_mac.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_energy_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("mac_write_energy", &u_write_energy_mac);

    // Initialize the unit cycle and energy of local buffer
    u_read_cycle_lb.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_cycle_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("lb_read_cycle", &u_read_cycle_lb);
    
    u_read_energy_lb.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_energy_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("lb_read_energy", &u_read_energy_lb);
    
    u_write_cycle_lb.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_cycle_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("lb_write_cycle", &u_write_cycle_lb);

    u_write_energy_lb.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_energy_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("lb_write_energy", &u_write_energy_lb);

    u_dynamic_power_lb.reserve(data_type_t::NUM_DATA_TYPES);
    u_dynamic_power_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("lb_dynamic_power", &u_dynamic_power_lb);

    u_static_power_lb.reserve(data_type_t::NUM_DATA_TYPES);
    u_static_power_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("lb_static_power", &u_static_power_lb);

    /* Initialize PE stats */

    // Total number of data request to local buffer
    num_request_to_lb.reserve(data_type_t::NUM_DATA_TYPES);
    num_request_to_lb.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Total number of data transfer to MAC unit
    num_data_transfer_to_mac.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer_to_mac.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Total access cycles to MAC unit
    access_cycle_mac.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Total access energies to MAC unit
    access_energy_mac.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Total access cycles to local buffer
    access_cycle_lb.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    
    // Total access energies to local buffer
    access_energy_lb.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Total transfer cycle between MAC unit and local buffer
    transfer_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Total transfer energies between MAC unit and local buffer
    transfer_energy.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Overlapped cycle between MAC units and local buffer
    cycle_mac_lb.reserve(data_type_t::NUM_DATA_TYPES);
    cycle_mac_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    utilization_local_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    utilization_local_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

}

// Connect PE to PE array.
void pe_t::connect(pe_array_t *m_pe_array) {
    pe_array = m_pe_array;
}

// Update tile size of local buffer and MAC unit.
void pe_t::update_tile_size(scheduler_t *m_scheduler) {

    num_active_macs = m_scheduler->num_active_mac;
    if(num_active_macs > num_macs*mac_width) {
        std::cerr << "The number of Active macs : " << num_active_macs 
                  << " is bigger than the number of macs : " << num_macs << std::endl;
        exit(1);
    }
    // Update tile sizes.
    tile_size_mac = m_scheduler->tile_size[component_type_t::MAC];
    tile_size_lb = m_scheduler->tile_size[component_type_t::PE];

}

void pe_t::update_offset() {
    // Update the offset of memory.
    if(memory_type == memory_type_t::SEPARATE) {
        offsets_lb[data_type_t::INPUT] = 0, offsets_lb[data_type_t::WEIGHT] = 0, offsets_lb[data_type_t::OUTPUT] = 0;
    }
    else if(memory_type == memory_type_t::SHARED) {
        offsets_lb[data_type_t::INPUT] = 0;
        offsets_lb[data_type_t::WEIGHT] = tile_size_lb[data_type_t::INPUT];
        offsets_lb[data_type_t::OUTPUT] = tile_size_lb[data_type_t::INPUT] + tile_size_lb[data_type_t::WEIGHT];
    }
}

void pe_t::check_tile_size() {
    // Check the validity.
    if(memory_type == memory_type_t::SEPARATE) {
        if((!bypass[data_type_t::INPUT] && tile_size_lb[data_type_t::INPUT]*sizeof(data_t) > input_size) || 
           (!bypass[data_type_t::WEIGHT] && tile_size_lb[data_type_t::WEIGHT]*sizeof(data_t) > weight_size) ||
           (!bypass[data_type_t::OUTPUT] && tile_size_lb[data_type_t::OUTPUT]*sizeof(data_t) > output_size)) {
            std::cerr << "The data size is bigger than local buffer size\n" 
                      << "Input data : " << tile_size_lb[data_type_t::INPUT]*sizeof(data_t) 
                      << " Input buffer : " << input_size << "\t"
                      << "Weight : " << tile_size_lb[data_type_t::WEIGHT]*sizeof(data_t)
                      << " Weight buffer : " << weight_size << "\t"
                      << "Output data : " << tile_size_lb[data_type_t::OUTPUT]*sizeof(data_t)
                      << " Output buffer : " << output_size << std::endl;
            exit(1);
        }
    }
    else if(memory_type == memory_type_t::SHARED) {
        unsigned data_size = tile_size_lb[data_type_t::INPUT] + tile_size_lb[data_type_t::WEIGHT] + tile_size_lb[data_type_t::OUTPUT];
        if(bypass[data_type_t::INPUT]) {
            data_size -= tile_size_lb[data_type_t::INPUT];
        }
        if(bypass[data_type_t::WEIGHT]) {
            data_size -= tile_size_lb[data_type_t::WEIGHT];
        }
        if(bypass[data_type_t::OUTPUT]) {
            data_size -= tile_size_lb[data_type_t::OUTPUT];
        }

        if(data_size*sizeof(data_t) > input_size + weight_size + output_size) {
            std::cerr << "The data size is bigger than local buffer size\n"
                      << "Data : " << data_size*sizeof(data_t)
                      << " Buffer : " << input_size + weight_size + output_size << std::endl;
         }
    }
}

// Get stationary_type of MAC register
stationary_type_t pe_t::get_mac_stationary_type() {
    return stationary_type_mac;
}

stationary_type_t pe_t::get_local_buffer_stationary_type() {
    return stationary_type_local_buffer;
}

// Return parameter order.
std::string pe_t::get_parameter_order() {
    return parameter_order;
}

memory_type_t pe_t::get_memory_type() {
    return memory_type;
}

// A signal to check whether the PE is idle state or not.
bool pe_t::is_idle() {
    return idle;
}

// A signal to check whether all data exist in the local buffer or not.
bool pe_t::is_exist_data() {
    if(exist_data_lb[data_type_t::INPUT] && exist_data_lb[data_type_t::WEIGHT] && exist_data_lb[data_type_t::OUTPUT]) {
        return true;
    }
    else {
        return false;
    }
}

// A signal to check whether at least one request exists from local buffer to PE array or not.
bool pe_t::is_exist_request() {
    if(request_to_pe_array[data_type_t::INPUT] || request_to_pe_array[data_type_t::WEIGHT] || request_to_pe_array[data_type_t::OUTPUT]) {
        return true;
    }
    else {
        return false;
    }
}

double pe_t::get_static_power() {
    double static_power = 0.0;
    if(get_memory_type() == memory_type_t::SEPARATE) { 
        static_power = u_static_power_mac 
                     + u_static_power_lb[data_type_t::INPUT] 
                     + u_static_power_lb[data_type_t::WEIGHT] 
                     + u_static_power_lb[data_type_t::OUTPUT];
    }
    else if(get_memory_type() == memory_type_t::SHARED) {
        static_power = u_static_power_mac + u_static_power_lb[data_type_t::INPUT];
    }
    return static_power;
}

// Wait for the data comes from Global buffer.
void pe_t::wait_data() {
    idle = true;
}

// A signal that the data exist in the PE.
void pe_t::fill_data() {
    idle = false;
}

// Request data to PE array.
void pe_t::request_data() {
    // If input data does not exist in Local buffer.
    if(!exist_data_lb[data_type_t::INPUT]) {
        // Request input data to PE array.
        request_to_pe_array[data_type_t::INPUT] = true;
        pe_array->exist_data[data_type_t::INPUT] = false;
        
        // Update Stats.
        pe_array->num_request[data_type_t::INPUT]++;
    }
    // If weight does not exist in Local buffer.
    if(!exist_data_lb[data_type_t::WEIGHT]) {
        // Request weight to PE array.
        request_to_pe_array[data_type_t::WEIGHT] = true;
        pe_array->exist_data[data_type_t::WEIGHT] = false;

        // Update Stats.
        pe_array->num_request[data_type_t::WEIGHT]++;
    }
    // If output data does not exist in Local buffer.
    if(!exist_data_lb[data_type_t::OUTPUT]) {
        // Request output data to PE array.
        request_to_pe_array[data_type_t::OUTPUT] = true;
        pe_array->exist_data[data_type_t::OUTPUT] = false;

        // Update stats.
        pe_array->num_request[data_type_t::OUTPUT]++;
    }
}

// Transfer data from local buffer to MAC unit.
// And execute MAC operation.
void pe_t::data_transfer_to_mac(scheduler_t *m_scheduler) {


    if(!bypass[data_type_t::INPUT]) {
        utilization_local_buffer[data_type_t::INPUT] = (float)(tile_size_lb[data_type_t::INPUT])/(float)(input_size);
    }
    if(bypass[data_type_t::WEIGHT]) {
        utilization_local_buffer[data_type_t::WEIGHT] = (float)(tile_size_lb[data_type_t::WEIGHT])/(float)(weight_size);
    }
    if(bypass[data_type_t::OUTPUT]) {
        utilization_local_buffer[data_type_t::OUTPUT] = (float)(tile_size_lb[data_type_t::OUTPUT])/(float)(output_size);
    }
    if(request_to_lb[data_type_t::INPUT]) {
#ifdef FUNCTIONAL       
        // Input data transfer 
        m_scheduler->transfer_data(input_data_mac, input_data_lb, 0, m_scheduler->input_offset_pe.front(),
                                   component_type_t::MAC, component_type_t::PE, 
                                   data_type_t::INPUT, get_mac_stationary_type(), action_type_t::LOAD);
        // Update for NPUsim ver2
        //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 && 
        //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
        //m_scheduler->transfer_data_ver2(input_data_mac, input_data_lb, 
        //                                component_type_t::MAC, component_type_t::PE, 
        //                                data_type_t::INPUT, get_mac_stationary_type(), action_type_t::LOAD, last_component);

        // Case 1. Dense format
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::INPUT]) {
                num_data_transfer_to_mac[data_type_t::INPUT]++;

                std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
                parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

                uint64_t address_mac = 0, address_lb = 0;
                unsigned num_access_mac = 0, num_access_lb = 0;
                for(unsigned b = 0; b < parameters_mac[parameter_type_t::BATCH_SIZE]; b++) {
                    for(unsigned c = 0; c < parameters_mac[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned h = 0; h < parameters_mac[parameter_type_t::INPUT_HEIGHT]; h++) {
                            for(unsigned w = 0; w < parameters_mac[parameter_type_t::INPUT_WIDTH]; w++) {
                                if(address_mac != ((uint64_t)&input_data_mac[b*parameters_mac[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_mac[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_mac[parameter_type_t::INPUT_WIDTH] +
                                                                             c*parameters_mac[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_mac[parameter_type_t::INPUT_WIDTH] + 
                                                                             h*parameters_mac[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                             mask_bits_mac[data_type_t::INPUT]) << mask_bits_mac[data_type_t::INPUT]) {
                                    // Update write cost to MAC unit
                                    access_energy_mac[data_type_t::INPUT] += u_write_energy_mac[data_type_t::INPUT];
                                    access_cycle_mac[data_type_t::INPUT] += u_write_cycle_mac[data_type_t::INPUT];
                                    num_access_mac++;

                                    // Update MAC address
                                    address_mac = ((uint64_t)&input_data_mac[b*parameters_mac[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_mac[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_mac[parameter_type_t::INPUT_WIDTH] +
                                                                             c*parameters_mac[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_mac[parameter_type_t::INPUT_WIDTH] + 
                                                                             h*parameters_mac[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                             mask_bits_mac[data_type_t::INPUT]) << mask_bits_mac[data_type_t::INPUT];
                                }
                                if(address_lb != ((uint64_t)&input_data_lb[m_scheduler->input_offset_pe.front() +
                                                                           b*parameters_lb[parameter_type_t::INPUT_CHANNEL]
                                                                            *parameters_lb[parameter_type_t::INPUT_HEIGHT]
                                                                            *parameters_lb[parameter_type_t::INPUT_WIDTH] +
                                                                           c*parameters_lb[parameter_type_t::INPUT_HEIGHT]
                                                                            *parameters_lb[parameter_type_t::INPUT_WIDTH] + 
                                                                           h*parameters_lb[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                           mask_bits_lb[data_type_t::INPUT]) << mask_bits_lb[data_type_t::INPUT]) {

                                    // Update read cost to local buffer.
                                    access_energy_lb[data_type_t::INPUT] += u_read_energy_lb[data_type_t::INPUT];
                                    access_cycle_lb[data_type_t::INPUT] += u_read_cycle_lb[data_type_t::INPUT];
                                    num_access_lb++;

                                    // Update local buffer address
                                    address_lb = ((uint64_t)&input_data_lb[m_scheduler->input_offset_pe.front() +
                                                                           b*parameters_lb[parameter_type_t::INPUT_CHANNEL]
                                                                            *parameters_lb[parameter_type_t::INPUT_HEIGHT]
                                                                            *parameters_lb[parameter_type_t::INPUT_WIDTH] +
                                                                           c*parameters_lb[parameter_type_t::INPUT_HEIGHT]
                                                                            *parameters_lb[parameter_type_t::INPUT_WIDTH] + 
                                                                           h*parameters_lb[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                           mask_bits_lb[data_type_t::INPUT]) << mask_bits_lb[data_type_t::INPUT];
                                }
                            }
                        }
                    }
                }
                
                // Update overlapped cycle at local buffer and MAC units.

                unsigned ratio = ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(line_size_mac[data_type_t::INPUT]));

                // At the 1, 2, before last, last stages
                double first_stage = u_read_cycle_lb[data_type_t::INPUT];
                double second_stage = std::max(u_read_cycle_lb[data_type_t::INPUT], 
                                                 u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth)));
                double last_before_stage = std::max(ratio*u_write_cycle_mac[data_type_t::INPUT], 
                                                      u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth)));
                double last_stage = ratio*u_write_cycle_mac[data_type_t::INPUT];

                // Remainder stages.
                double other_stage = std::max(u_read_cycle_lb[data_type_t::INPUT], 
                                       std::max(ratio*u_write_cycle_mac[data_type_t::INPUT], 
                                                u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth))));

                if(num_access_lb == 1) {
                    cycle_mac_lb[data_type_t::INPUT] += u_read_cycle_lb[data_type_t::INPUT] + 
                                                        u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth)) + 
                                                        ratio*u_write_cycle_mac[data_type_t::INPUT];
                } else {
                    cycle_mac_lb[data_type_t::INPUT] += first_stage + second_stage + 
                                                        (num_access_lb-2)*other_stage + 
                                                        last_before_stage + last_stage;
                } 
                // Update per data stats between MAC unit and local buffer
                transfer_cycle[data_type_t::INPUT] += num_access_lb*u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth));
                transfer_energy[data_type_t::INPUT] += num_access_lb*u_transfer_energy*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth));
            }
            if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                move_front(&m_scheduler->input_offset_pe);
            }
        }
        // Case 2. COO format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
            std::cerr << "Sparse format COO is not supported in this version" << std::endl;
            exit(1);
        }
        // Case 3. CSC format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::INPUT]) {
                // Row index calculation
                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::MAC);
                unsigned row = parameters[parameter_type_t::INPUT_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer_to_mac[data_type_t::INPUT]++;

                // Update local buffer access cost
                access_cycle_lb[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                      *u_read_cycle_lb[data_type_t::INPUT] + // Value
                                                       (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                      *u_read_cycle_lb[data_type_t::INPUT]
                                                      /(sizeof(data_t)*8/row_bit) + // Row index
                                                       parameters[parameter_type_t::BATCH_SIZE]
                                                      *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                      *(parameters[parameter_type_t::INPUT_WIDTH]+1)*u_read_cycle_lb[data_type_t::INPUT]
                                                      /(sizeof(data_t)*8/row_bit); // Column pointer

                access_energy_lb[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_read_energy_lb[data_type_t::INPUT] + // Non-zero data
                                                        (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_read_energy_lb[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/row_bit) + // Row index
                                                        parameters[parameter_type_t::BATCH_SIZE]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                       *u_read_energy_lb[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/row_bit); // Column pointer

                // Update MAC access cost
                access_cycle_mac[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_write_cycle_mac[data_type_t::INPUT] + // Non-zero data
                                                        (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_write_cycle_mac[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/row_bit) + // Row index
                                                        parameters[parameter_type_t::BATCH_SIZE]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                       *u_write_cycle_mac[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/row_bit); // Column pointer

                access_energy_mac[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                        *u_write_energy_mac[data_type_t::INPUT] + // DATA
                                                         (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                        *u_write_energy_mac[data_type_t::INPUT]
                                                        /(sizeof(data_t)*8/row_bit) + // Row index
                                                         parameters[parameter_type_t::BATCH_SIZE]
                                                        *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                        *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                        *u_write_energy_mac[data_type_t::INPUT]
                                                        /(sizeof(data_t)*8/row_bit); // Column pointer

                /* TODO */
                cycle_mac_lb[data_type_t::INPUT] += std::max((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_cycle_lb[data_type_t::INPUT] +  // Non-zero data
                                                             (tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_cycle_lb[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/row_bit)  + // Row index
                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                            *u_read_cycle_lb[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/row_bit), // Column pointer
                                                             (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_write_cycle_mac[data_type_t::INPUT] + // Non-zero data
                                                             (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_write_cycle_mac[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/row_bit) + // Row index
                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                            *u_write_cycle_mac[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/row_bit)); // Column pointer

                // Update transfer cost between MAC units and local buffers
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                      u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)/(double)bitwidth) + // Row index
                                                      u_transfer_cycle*ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *(parameters[parameter_type_t::INPUT_WIDTH+1])*row_bit)/(double)bitwidth);  // Column pointer
                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)/(double)bitwidth) + // Row index 
                                                       u_transfer_energy*ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                      *parameters[parameter_type_t::INPUT_CHANNEL]
                                                                                      /parameters[parameter_type_t::GROUP]
                                                                                      *(parameters[parameter_type_t::INPUT_WIDTH+1])*row_bit)/(double)bitwidth); // Column pointer
                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->input_offset_pe);
                }
            }
        }
        // Case 4. CSR format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::INPUT]) {
                // Column index calculation
                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::MAC);
                unsigned column = parameters[parameter_type_t::INPUT_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer_to_mac[data_type_t::INPUT]++;

                // Update local buffer access cost
                access_cycle_lb[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                      *u_read_cycle_lb[data_type_t::INPUT] + // Non-zero data
                                                       (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                      *u_read_cycle_lb[data_type_t::INPUT]
                                                      /(sizeof(data_t)*8/column_bit) + // Column index
                                                       parameters[parameter_type_t::BATCH_SIZE]
                                                      *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                      *(parameters[parameter_type_t::INPUT_HEIGHT]+1)*u_read_cycle_lb[data_type_t::INPUT]
                                                      /(sizeof(data_t)*8/column_bit); // Row pointer

                access_energy_lb[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_read_energy_lb[data_type_t::INPUT] + // Non-zero data
                                                        (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_read_energy_lb[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/column_bit) + // Column index
                                                        parameters[parameter_type_t::BATCH_SIZE]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                       *u_read_energy_lb[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/column_bit); // Row pointer

                // Update MAC access cost
                access_cycle_mac[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_write_cycle_mac[data_type_t::INPUT] + // Non-zero data
                                                        (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_write_cycle_mac[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/column_bit) + // Column index
                                                        parameters[parameter_type_t::BATCH_SIZE]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                       *u_write_cycle_mac[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/column_bit); // Row pointer

                access_energy_mac[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                        *u_write_energy_mac[data_type_t::INPUT] + // DATA
                                                         (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                        *u_write_energy_mac[data_type_t::INPUT]
                                                        /(sizeof(data_t)*8/column_bit) + // Column index
                                                         parameters[parameter_type_t::BATCH_SIZE]
                                                        *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                        *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                        *u_write_energy_mac[data_type_t::INPUT]
                                                        /(sizeof(data_t)*8/column_bit);

                /* TODO */
                cycle_mac_lb[data_type_t::INPUT] += std::max((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_cycle_lb[data_type_t::INPUT] +  // Non-zero data
                                                             (tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_cycle_lb[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/column_bit)  + // Column index
                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                            *u_read_cycle_lb[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/column_bit), // Row pointer
                                                             (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_write_cycle_mac[data_type_t::INPUT] + // Non-zero data
                                                             (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_write_cycle_mac[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/column_bit) + // Column index
                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                            *u_write_cycle_mac[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/column_bit)); // Row pointer

                // Update transfer cost between MAC units and local buffers
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                      u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)/(double)bitwidth) + // Column index
                                                      u_transfer_cycle*ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *(parameters[parameter_type_t::INPUT_HEIGHT+1])*column_bit)/(double)bitwidth);  // Row pointer
                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)/(double)bitwidth) + // Column index 
                                                       u_transfer_energy*ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                      *parameters[parameter_type_t::INPUT_CHANNEL]
                                                                                      /parameters[parameter_type_t::GROUP]
                                                                                      *(parameters[parameter_type_t::INPUT_HEIGHT+1])*column_bit)/(double)bitwidth); // Row pointer

                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->input_offset_pe);
                }
            }
        }
        // Case 5. SparseMap 
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::INPUT]) {
                num_data_transfer_to_mac[data_type_t::INPUT]++;

                // Update local buffer access cost
                access_cycle_lb[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                      *u_read_cycle_lb[data_type_t::INPUT] + // Non-zero value
                                                       tile_size_mac[data_type_t::INPUT]
                                                      *u_read_cycle_lb[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
                access_energy_lb[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_read_energy_lb[data_type_t::INPUT] + // Non-zero data
                                                        tile_size_mac[data_type_t::INPUT]
                                                       *u_read_energy_lb[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata

                // Update MAC access cost
                access_cycle_mac[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_write_cycle_mac[data_type_t::INPUT] + // Non-zero data
                                                        tile_size_mac[data_type_t::INPUT]
                                                       *u_write_cycle_mac[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
                access_energy_mac[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*u_write_energy_mac[data_type_t::INPUT]
                                                       + tile_size_mac[data_type_t::INPUT]*u_write_energy_mac[data_type_t::INPUT]/(sizeof(data_t)*8);
    
                // Update overlapped cycle between MAC and local buffers
                cycle_mac_lb[data_type_t::INPUT] += std::max((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *u_read_cycle_lb[data_type_t::INPUT] + // Non-zero data
                                                              tile_size_mac[data_type_t::INPUT]
                                                             *u_read_cycle_lb[data_type_t::INPUT]/(sizeof(data_t)*8), // Metadata
                                                              (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *u_write_cycle_mac[data_type_t::INPUT] +  // Non-zero data
                                                              tile_size_mac[data_type_t::INPUT]
                                                             *u_write_cycle_mac[data_type_t::INPUT]/(sizeof(data_t)*8)); // Metadata 

                // Transfer cycle between MAC unit and local buffer = transfer cycle for data + transfer cycle for index vector.
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(double)bitwidth) +  // Non-zero data
                                                      u_transfer_cycle*ceil((double)(tile_size_mac[data_type_t::INPUT])/(double)bitwidth); // Metadata

                // Transfer energy between MAC unit and local buffer = transfer energy for data + transfer energy for index vector
                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((double)(tile_size_mac[data_type_t::INPUT])/(double)bitwidth); // Metadata
                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->input_offset_pe);
                }
            }
        }
        else {
            std::cerr << "Undefined compression type" << std::endl;
            exit(1);
        }

#else   // Timing simulation
        // Update the number of data transfer from local buffer to MAC unit
        if(!skip_transfer[data_type_t::INPUT]) {
            num_data_transfer_to_mac[data_type_t::INPUT]++;

            std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
            parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

            uint64_t address_mac = 0, address_lb = 0;
            unsigned num_access_mac = 0, num_access_lb = 0;
            for(unsigned b = 0; b < parameters_mac[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned c = 0; c < parameters_mac[parameter_type_t::INPUT_CHANNEL]; c++) {
                    for(unsigned h = 0; h < parameters_mac[parameter_type_t::INPUT_HEIGHT]; h++) {
                        for(unsigned w = 0; w < parameters_mac[parameter_type_t::INPUT_WIDTH]; w++) {
                            if(address_mac != ((uint64_t)&input_data_mac[b*parameters_mac[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_mac[parameter_type_t::INPUT_HEIGHT]
                                                                          *parameters_mac[parameter_type_t::INPUT_WIDTH] +
                                                                         c*parameters_mac[parameter_type_t::INPUT_HEIGHT]
                                                                          *parameters_mac[parameter_type_t::INPUT_WIDTH] + 
                                                                         h*parameters_mac[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                         mask_bits_mac[data_type_t::INPUT]) << mask_bits_mac[data_type_t::INPUT]) {
                                // Update write cost to MAC unit
                                access_energy_mac[data_type_t::INPUT] += u_write_energy_mac[data_type_t::INPUT];
                                access_cycle_mac[data_type_t::INPUT] += u_write_cycle_mac[data_type_t::INPUT];
                                num_access_mac++;

                                // Update MAC address
                                address_mac = ((uint64_t)&input_data_mac[b*parameters_mac[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_mac[parameter_type_t::INPUT_HEIGHT]
                                                                          *parameters_mac[parameter_type_t::INPUT_WIDTH] +
                                                                         c*parameters_mac[parameter_type_t::INPUT_HEIGHT]
                                                                          *parameters_mac[parameter_type_t::INPUT_WIDTH] + 
                                                                         h*parameters_mac[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                         mask_bits_mac[data_type_t::INPUT]) << mask_bits_mac[data_type_t::INPUT];
                            }
                            if(address_lb != ((uint64_t)&input_data_lb[m_scheduler->input_offset_pe.front() +
                                                                       b*parameters_lb[parameter_type_t::INPUT_CHANNEL]
                                                                        *parameters_lb[parameter_type_t::INPUT_HEIGHT]
                                                                        *parameters_lb[parameter_type_t::INPUT_WIDTH] +
                                                                       c*parameters_lb[parameter_type_t::INPUT_HEIGHT]
                                                                        *parameters_lb[parameter_type_t::INPUT_WIDTH] + 
                                                                       h*parameters_lb[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                       mask_bits_lb[data_type_t::INPUT]) << mask_bits_lb[data_type_t::INPUT]) {

                                // Update read cost to local buffer.
                                access_energy_lb[data_type_t::INPUT] += u_read_energy_lb[data_type_t::INPUT];
                                access_cycle_lb[data_type_t::INPUT] += u_read_cycle_lb[data_type_t::INPUT];
                                num_access_lb++;

                                // Update local buffer address
                                address_lb = ((uint64_t)&input_data_lb[m_scheduler->input_offset_pe.front() + 
                                                                       b*parameters_lb[parameter_type_t::INPUT_CHANNEL]
                                                                        *parameters_lb[parameter_type_t::INPUT_HEIGHT]
                                                                        *parameters_lb[parameter_type_t::INPUT_WIDTH] +
                                                                       c*parameters_lb[parameter_type_t::INPUT_HEIGHT]
                                                                        *parameters_lb[parameter_type_t::INPUT_WIDTH] + 
                                                                       h*parameters_lb[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                       mask_bits_lb[data_type_t::INPUT]) << mask_bits_lb[data_type_t::INPUT];
                            }
                        }
                    }
                }
            }
            
            // Update overlapped cycle at local buffer and MAC units.

            unsigned ratio = ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(line_size_mac[data_type_t::INPUT]));

            // At the 1, 2, before last, last stages
            double first_stage = u_read_cycle_lb[data_type_t::INPUT];
            double second_stage = std::max(u_read_cycle_lb[data_type_t::INPUT], 
                                             u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth)));
            double last_before_stage = std::max(ratio*u_write_cycle_mac[data_type_t::INPUT], 
                                                  u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth)));
            double last_stage = ratio*u_write_cycle_mac[data_type_t::INPUT];

            // Remainder stages.
            double other_stage = std::max(u_read_cycle_lb[data_type_t::INPUT], 
                                   std::max(ratio*u_write_cycle_mac[data_type_t::INPUT], 
                                            u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth))));

            if(num_access_lb == 1) {
                cycle_mac_lb[data_type_t::INPUT] += u_read_cycle_lb[data_type_t::INPUT] + 
                                                    u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth)) + 
                                                    ratio*u_write_cycle_mac[data_type_t::INPUT];
            } else {
                cycle_mac_lb[data_type_t::INPUT] += first_stage + second_stage + 
                                                    (num_access_lb-2)*other_stage + 
                                                    last_before_stage + last_stage;
            } 
            // Update per data stats between MAC unit and local buffer
            transfer_cycle[data_type_t::INPUT] += num_access_lb*u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth));
            transfer_energy[data_type_t::INPUT] += num_access_lb*u_transfer_energy*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth));

            if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                move_front(&m_scheduler->input_offset_pe);
            }
        }
#endif

        // Increment the input index that finds the offset of input data in local buffer.
        input_index++;
        // Input data exists in MAC unit and request input data is not required.
        exist_data_mac[data_type_t::INPUT] = true, request_to_lb[data_type_t::INPUT] = false;

        if(tile_size_mac[data_type_t::INPUT] == tile_size_lb[data_type_t::INPUT]) { skip_transfer[data_type_t::INPUT] = true; }
    }
    // Transfer weight from local buffer to MAC unit.
    if(request_to_lb[data_type_t::WEIGHT]) {
#ifdef FUNCTIONAL

        // Weight data transfer 
        m_scheduler->transfer_data(weight_mac, weight_lb, 0, m_scheduler->weight_offset_pe.front(),
                                   component_type_t::MAC, component_type_t::PE, 
                                   data_type_t::WEIGHT, get_mac_stationary_type(), action_type_t::LOAD);
        // Update for NPUsim ver2
        //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 && 
        //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
        //m_scheduler->transfer_data_ver2(weight_mac, weight_lb, 
        //                                component_type_t::MAC, component_type_t::PE, 
        //                                data_type_t::WEIGHT, get_mac_stationary_type(), action_type_t::LOAD, last_component);

        // Case 1. Dense data format
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                num_data_transfer_to_mac[data_type_t::WEIGHT]++;

                std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
                parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

                uint64_t address_mac = 0, address_lb = 0;
                unsigned num_access_mac = 0, num_access_lb = 0;
                for(unsigned k = 0; k < parameters_mac[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned c = 0; c < parameters_mac[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned r = 0; r < parameters_mac[parameter_type_t::FILTER_HEIGHT]; r++) {
                            for(unsigned s = 0; s < parameters_mac[parameter_type_t::FILTER_WIDTH]; s++) {
                                if(address_mac != ((uint64_t)&weight_mac[k*parameters_mac[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_mac[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_mac[parameter_type_t::FILTER_WIDTH] +
                                                                         c*parameters_mac[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_mac[parameter_type_t::FILTER_WIDTH] + 
                                                                         r*parameters_mac[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                         mask_bits_mac[data_type_t::WEIGHT]) << mask_bits_mac[data_type_t::WEIGHT]) {

                                    // Update MAC cost.
                                    access_energy_mac[data_type_t::WEIGHT] += u_write_energy_mac[data_type_t::WEIGHT];
                                    access_cycle_mac[data_type_t::WEIGHT] += u_write_cycle_mac[data_type_t::WEIGHT];
                                    num_access_mac++;

                                    // Update MAC address
                                    address_mac = ((uint64_t)&weight_mac[k*parameters_mac[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_mac[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_mac[parameter_type_t::FILTER_WIDTH] +
                                                                         c*parameters_mac[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_mac[parameter_type_t::FILTER_WIDTH] + 
                                                                         r*parameters_mac[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                         mask_bits_mac[data_type_t::WEIGHT]) << mask_bits_mac[data_type_t::WEIGHT];
                                }

                                if(address_lb != ((uint64_t)&weight_lb[m_scheduler->weight_offset_pe.front() +
                                                                       k*parameters_lb[parameter_type_t::INPUT_CHANNEL]
                                                                        *parameters_lb[parameter_type_t::FILTER_HEIGHT]
                                                                        *parameters_lb[parameter_type_t::FILTER_WIDTH] +
                                                                       c*parameters_lb[parameter_type_t::FILTER_HEIGHT]
                                                                        *parameters_lb[parameter_type_t::FILTER_WIDTH] + 
                                                                       r*parameters_lb[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                       mask_bits_lb[data_type_t::WEIGHT]) << mask_bits_lb[data_type_t::WEIGHT]) {

                                    access_energy_lb[data_type_t::WEIGHT] += u_read_energy_lb[data_type_t::WEIGHT];
                                    access_cycle_lb[data_type_t::WEIGHT] += u_read_cycle_lb[data_type_t::WEIGHT];
                                    num_access_lb++;

                                    // Update Local buffer address
                                    address_lb = ((uint64_t)&weight_lb[m_scheduler->weight_offset_pe.front() +
                                                                       k*parameters_lb[parameter_type_t::INPUT_CHANNEL]
                                                                        *parameters_lb[parameter_type_t::FILTER_HEIGHT]
                                                                        *parameters_lb[parameter_type_t::FILTER_WIDTH] +
                                                                       c*parameters_lb[parameter_type_t::FILTER_HEIGHT]
                                                                        *parameters_lb[parameter_type_t::FILTER_WIDTH] + 
                                                                       r*parameters_lb[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                       mask_bits_lb[data_type_t::WEIGHT]) << mask_bits_lb[data_type_t::WEIGHT];
                                }
                            }
                        }
                    }
                }
                
                // Update overlapped cycle at local buffer and MAC units.
                unsigned ratio = ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(line_size_mac[data_type_t::WEIGHT]));

                // At the 1, 2, before last, last stages
                unsigned first_stage = u_read_cycle_lb[data_type_t::WEIGHT];
                unsigned second_stage = std::max(u_read_cycle_lb[data_type_t::WEIGHT], 
                                                 u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)));
                unsigned last_before_stage = std::max(ratio*u_write_cycle_mac[data_type_t::WEIGHT], 
                                                      u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)));
                unsigned last_stage = ratio*u_write_cycle_mac[data_type_t::WEIGHT];

                // Remainder stages
                unsigned other_stage = std::max(u_read_cycle_lb[data_type_t::WEIGHT],
                                       std::max(ratio*u_write_cycle_mac[data_type_t::WEIGHT], 
                                                u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth))));

                if(num_access_lb == 1) {
                    cycle_mac_lb[data_type_t::WEIGHT] += u_read_cycle_lb[data_type_t::WEIGHT] + 
                                                         u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)) +
                                                         ratio*u_write_cycle_mac[data_type_t::WEIGHT];
                } else {
                    cycle_mac_lb[data_type_t::WEIGHT] += first_stage + second_stage + 
                                                         (num_access_lb-2)*other_stage + 
                                                         last_before_stage + last_stage;
                }

                // Update per data stats between MAC unit and local buffer
                transfer_cycle[data_type_t::WEIGHT] += num_access_lb*u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)bitwidth);
                transfer_energy[data_type_t::WEIGHT] += num_access_lb*u_transfer_energy*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)bitwidth);
        
                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->weight_offset_pe);
                }
            }
        }
        // Case 2. COO data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
            std::cerr << "Data format COO is not supported in this version." << std::endl;
            exit(1);
        }
        // Case 3. CSC format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                // Row bit calculation
                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::MAC);
                unsigned row = parameters[parameter_type_t::FILTER_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer_to_mac[data_type_t::WEIGHT]++;

                // Update local buffer access cost
                access_cycle_lb[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle_lb[data_type_t::WEIGHT] + // Non-zero data
                                                        (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle_lb[data_type_t::WEIGHT]
                                                       /(sizeof(data_t)*8/row_bit) +  // Row index
                                                       parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                      *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                      *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                      *u_read_cycle_lb[data_type_t::WEIGHT]
                                                      /(sizeof(data_t)*8/row_bit); // Column pointer
                access_energy_lb[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_read_energy_lb[data_type_t::WEIGHT] + // Non-zero data
                                                         (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_read_energy_lb[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/row_bit) + // Row index
                                                         parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                        *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                        *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                        *u_read_energy_lb[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/row_bit); // Column pointer

                // Update MAC access cost
                access_cycle_mac[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_write_cycle_mac[data_type_t::WEIGHT] + // Non-zero data
                                                         (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_write_cycle_mac[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/row_bit) + // Row index
                                                         parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                        *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                        *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                        *u_write_cycle_mac[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/row_bit); // Column pointer
                access_energy_mac[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                         *u_write_energy_mac[data_type_t::WEIGHT] + // Non-zero data
                                                          (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                         *u_write_energy_mac[data_type_t::WEIGHT]
                                                         /(sizeof(data_t)*8/row_bit) + // Row index
                                                          parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                         /parameters[parameter_type_t::GROUP]
                                                         *parameters[parameter_type_t::INPUT_CHANNEL]
                                                         /parameters[parameter_type_t::GROUP]
                                                         *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                         *u_write_energy_mac[data_type_t::WEIGHT]
                                                         /(sizeof(data_t)*8/row_bit); // Column pointer

                /* TODO */
                // Update overlapped cycle between MAC units and local buffers
                cycle_mac_lb[data_type_t::WEIGHT] += std::max((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_cycle_lb[data_type_t::WEIGHT] + // Non-zero data
                                                              (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_cycle_lb[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/row_bit) + // Row index
                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                             *u_read_cycle_lb[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/row_bit), // Column pointer
                                                              (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_write_cycle_mac[data_type_t::WEIGHT] + // Non-zero data
                                                              (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_write_cycle_mac[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/row_bit) + // Row index
                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                             *u_write_cycle_mac[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/row_bit)); // Column pointer

                // Update transfer cost between MAC units and local buffers
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                       u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*row_bit)/(double)bitwidth) +  // Row index
                                                       u_transfer_cycle*ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                     *(parameters[parameter_type_t::FILTER_WIDTH+1])
                                                                                     *row_bit)/(double)bitwidth); // Column pointer
                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                        u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*row_bit)/(double)bitwidth) + // Row index
                                                        u_transfer_energy*ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                                                       /parameters[parameter_type_t::GROUP]
                                                                                       *parameters[parameter_type_t::INPUT_CHANNEL]
                                                                                       /parameters[parameter_type_t::GROUP]
                                                                                       *(parameters[parameter_type_t::FILTER_WIDTH+1])
                                                                                       *row_bit)/(double)bitwidth); // Column pointer
                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->weight_offset_pe);
                }
            }
        }
        // Case 4. CSR format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                // Column bit calculation
                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::MAC);
                unsigned column = parameters[parameter_type_t::FILTER_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer_to_mac[data_type_t::WEIGHT]++;

                // Update local buffer access cost
                access_cycle_lb[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle_lb[data_type_t::WEIGHT] + // Non-zero data
                                                        (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle_lb[data_type_t::WEIGHT]
                                                       /(sizeof(data_t)*8/column_bit) +  // Column index
                                                       parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                      /parameters[parameter_type_t::GROUP]
                                                      *parameters[parameter_type_t::INPUT_CHANNEL]
                                                      /parameters[parameter_type_t::GROUP]
                                                      *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                      *u_read_cycle_lb[data_type_t::WEIGHT]
                                                      /(sizeof(data_t)*8/column_bit); // Row pointer
                access_energy_lb[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_read_energy_lb[data_type_t::WEIGHT] + // Non-zero data
                                                         (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_read_energy_lb[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/column_bit) + // Column index
                                                         parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                        /parameters[parameter_type_t::GROUP]
                                                        *parameters[parameter_type_t::INPUT_CHANNEL]
                                                        /parameters[parameter_type_t::GROUP]
                                                        *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                        *u_read_energy_lb[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/column_bit); // Row pointer

                // Update MAC access cost
                access_cycle_mac[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_write_cycle_mac[data_type_t::WEIGHT] + // Non-zero data
                                                         (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_write_cycle_mac[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/column_bit) + // Column index
                                                         parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                        /parameters[parameter_type_t::GROUP]
                                                        *parameters[parameter_type_t::INPUT_CHANNEL]
                                                        /parameters[parameter_type_t::GROUP]
                                                        *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                        *u_write_cycle_mac[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/column_bit); // Row pointer
                access_energy_mac[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                         *u_write_energy_mac[data_type_t::WEIGHT] + // Non-zero data
                                                          (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                         *u_write_energy_mac[data_type_t::WEIGHT]
                                                         /(sizeof(data_t)*8/column_bit) + // Column index
                                                          parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                         /parameters[parameter_type_t::GROUP]
                                                         *parameters[parameter_type_t::INPUT_CHANNEL]
                                                         /parameters[parameter_type_t::GROUP]
                                                         *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                         *u_write_energy_mac[data_type_t::WEIGHT]
                                                         /(sizeof(data_t)*8/column_bit); // Row pointer

                /* TODO */
                // Update overlapped cycle between MAC units and local buffers
                cycle_mac_lb[data_type_t::WEIGHT] += std::max((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_cycle_lb[data_type_t::WEIGHT] + // Non-zero data
                                                              (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_cycle_lb[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/column_bit) + // Column index
                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                             /parameters[parameter_type_t::GROUP]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]
                                                             /parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                             *u_read_cycle_lb[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/column_bit), // Row pointer
                                                              (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_write_cycle_mac[data_type_t::WEIGHT] + // Non-zero data
                                                              (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_write_cycle_mac[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/column_bit) + // Column index
                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                             /parameters[parameter_type_t::GROUP]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]
                                                             /parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                             *u_write_cycle_mac[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/column_bit)); // Row pointer

                // Update transfer cost between MAC units and local buffers
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                       u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*column_bit)/(double)bitwidth) +  // Column index
                                                       u_transfer_cycle*ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                                                     /parameters[parameter_type_t::GROUP]
                                                                                     *parameters[parameter_type_t::INPUT_CHANNEL]
                                                                                     /parameters[parameter_type_t::GROUP]
                                                                                     *(parameters[parameter_type_t::FILTER_HEIGHT+1])
                                                                                     *column_bit)/(double)bitwidth); // Row pointer
                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                        u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*column_bit)/(double)bitwidth) + // Column index
                                                        u_transfer_energy*ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                                                       /parameters[parameter_type_t::GROUP]
                                                                                       *parameters[parameter_type_t::INPUT_CHANNEL]
                                                                                       /parameters[parameter_type_t::GROUP]
                                                                                       *(parameters[parameter_type_t::FILTER_HEIGHT+1])
                                                                                       *column_bit)/(double)bitwidth); // Row pointer
               
                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->weight_offset_pe);
                }
            }
        }
        // Case 4. SparseMap
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                num_data_transfer_to_mac[data_type_t::WEIGHT]++;

                // Update local buffer access cost
                access_cycle_lb[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle_lb[data_type_t::WEIGHT] + // Non-zero data
                                                        tile_size_mac[data_type_t::WEIGHT]
                                                       *u_read_cycle_lb[data_type_t::WEIGHT]/(sizeof(data_t)*8); // Metadata
                access_energy_lb[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_read_energy_lb[data_type_t::WEIGHT] + // Non-zero data
                                                         tile_size_mac[data_type_t::WEIGHT]
                                                        *u_read_energy_lb[data_type_t::WEIGHT]/(sizeof(data_t)*8); // Metadata
                // Update MAC access cost
                access_cycle_mac[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_write_cycle_mac[data_type_t::WEIGHT] + // Non-zero data
                                                         tile_size_mac[data_type_t::WEIGHT]
                                                        *u_write_cycle_mac[data_type_t::WEIGHT]/(sizeof(data_t)*8); // Metadata

                access_energy_mac[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*u_write_energy_mac[data_type_t::WEIGHT]
                                                        + tile_size_mac[data_type_t::WEIGHT]*u_write_energy_mac[data_type_t::WEIGHT]/(sizeof(data_t)*8);
                
                cycle_mac_lb[data_type_t::WEIGHT] += std::max((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*u_read_cycle_lb[data_type_t::WEIGHT]
                                                                         + tile_size_mac[data_type_t::WEIGHT]*u_read_cycle_lb[data_type_t::WEIGHT]/(sizeof(data_t)*8), 
                                                                         (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*u_write_cycle_mac[data_type_t::WEIGHT]
                                                                         + tile_size_mac[data_type_t::WEIGHT]*u_write_cycle_mac[data_type_t::WEIGHT]/(sizeof(data_t)*8));

                // Transfer cycle between MAC unit and Local buffer = transfer cycle for data + transfer cycle for index vector
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(double)bitwidth)
                                                     + u_transfer_cycle*ceil((double)(tile_size_mac[data_type_t::WEIGHT])/(double)bitwidth);
                // Transfer energy between MAC unit and local buffer = transfer energy for data + transfer energy for index vector
                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((double)(tile_size_mac[data_type_t::WEIGHT])/(double)bitwidth); // Metadata

                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->weight_offset_pe);
                }
            }
        }
        else {
            std::cerr << "Undefined compression type" << std::endl;
            exit(1);
        }
#else
        if(!skip_transfer[data_type_t::WEIGHT]) {
            num_data_transfer_to_mac[data_type_t::WEIGHT]++;

            std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
            parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

            uint64_t address_mac = 0, address_lb = 0;
            unsigned num_access_mac = 0, num_access_lb = 0;
            for(unsigned k = 0; k < parameters_mac[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                for(unsigned c = 0; c < parameters_mac[parameter_type_t::INPUT_CHANNEL]; c++) {
                    for(unsigned r = 0; r < parameters_mac[parameter_type_t::FILTER_HEIGHT]; r++) {
                        for(unsigned s = 0; s < parameters_mac[parameter_type_t::FILTER_WIDTH]; s++) {
                            if(address_mac != ((uint64_t)&weight_mac[k*parameters_mac[parameter_type_t::INPUT_CHANNEL]
                                                                      *parameters_mac[parameter_type_t::FILTER_HEIGHT]
                                                                      *parameters_mac[parameter_type_t::FILTER_WIDTH] +
                                                                     c*parameters_mac[parameter_type_t::FILTER_HEIGHT]
                                                                      *parameters_mac[parameter_type_t::FILTER_WIDTH] + 
                                                                     r*parameters_mac[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                     mask_bits_mac[data_type_t::WEIGHT]) << mask_bits_mac[data_type_t::WEIGHT]) {

                                // Update MAC cost.
                                access_energy_mac[data_type_t::WEIGHT] += u_write_energy_mac[data_type_t::WEIGHT];
                                access_cycle_mac[data_type_t::WEIGHT] += u_write_cycle_mac[data_type_t::WEIGHT];
                                num_access_mac++;

                                // Update MAC address
                                address_mac = ((uint64_t)&weight_mac[k*parameters_mac[parameter_type_t::INPUT_CHANNEL]
                                                                      *parameters_mac[parameter_type_t::FILTER_HEIGHT]
                                                                      *parameters_mac[parameter_type_t::FILTER_WIDTH] +
                                                                     c*parameters_mac[parameter_type_t::FILTER_HEIGHT]
                                                                      *parameters_mac[parameter_type_t::FILTER_WIDTH] + 
                                                                     r*parameters_mac[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                     mask_bits_mac[data_type_t::WEIGHT]) << mask_bits_mac[data_type_t::WEIGHT];
                            }

                            if(address_lb != ((uint64_t)&weight_lb[m_scheduler->weight_offset_pe.front() +
                                                                   k*parameters_lb[parameter_type_t::INPUT_CHANNEL]
                                                                    *parameters_lb[parameter_type_t::FILTER_HEIGHT]
                                                                    *parameters_lb[parameter_type_t::FILTER_WIDTH] +
                                                                   c*parameters_lb[parameter_type_t::FILTER_HEIGHT]
                                                                    *parameters_lb[parameter_type_t::FILTER_WIDTH] + 
                                                                   r*parameters_lb[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                   mask_bits_lb[data_type_t::WEIGHT]) << mask_bits_lb[data_type_t::WEIGHT]) {

                                access_energy_lb[data_type_t::WEIGHT] += u_read_energy_lb[data_type_t::WEIGHT];
                                access_cycle_lb[data_type_t::WEIGHT] += u_read_cycle_lb[data_type_t::WEIGHT];
                                num_access_lb++;

                                // Update Local buffer address
                                address_lb = ((uint64_t)&weight_lb[m_scheduler->weight_offset_pe.front() +
                                                                   k*parameters_lb[parameter_type_t::INPUT_CHANNEL]
                                                                    *parameters_lb[parameter_type_t::FILTER_HEIGHT]
                                                                    *parameters_lb[parameter_type_t::FILTER_WIDTH] +
                                                                   c*parameters_lb[parameter_type_t::FILTER_HEIGHT]
                                                                    *parameters_lb[parameter_type_t::FILTER_WIDTH] + 
                                                                   r*parameters_lb[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                   mask_bits_lb[data_type_t::WEIGHT]) << mask_bits_lb[data_type_t::WEIGHT];
                            }
                        }
                    }
                }
            }
            
            // Update overlapped cycle at local buffer and MAC units.

            unsigned ratio = ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(line_size_mac[data_type_t::WEIGHT]));

            // At the 1, 2, before last, last stages
            unsigned first_stage = u_read_cycle_lb[data_type_t::WEIGHT];
            unsigned second_stage = std::max(u_read_cycle_lb[data_type_t::WEIGHT], 
                                             u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)));
            unsigned last_before_stage = std::max(ratio*u_write_cycle_mac[data_type_t::WEIGHT], 
                                                  u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)));
            unsigned last_stage = ratio*u_write_cycle_mac[data_type_t::WEIGHT];

            // Remainder stages
            unsigned other_stage = std::max(u_read_cycle_lb[data_type_t::WEIGHT],
                                   std::max(ratio*u_write_cycle_mac[data_type_t::WEIGHT], 
                                            u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth))));

            if(num_access_lb == 1) {
                cycle_mac_lb[data_type_t::WEIGHT] += u_read_cycle_lb[data_type_t::WEIGHT] + 
                                                     u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)) +
                                                     ratio*u_write_cycle_mac[data_type_t::WEIGHT];
            } else {
                cycle_mac_lb[data_type_t::WEIGHT] += first_stage + second_stage + 
                                                     (num_access_lb-2)*other_stage + 
                                                     last_before_stage + last_stage;
            }

            // Update per data stats between MAC unit and local buffer
            transfer_cycle[data_type_t::WEIGHT] += num_access_lb*u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)bitwidth);
            transfer_energy[data_type_t::WEIGHT] += num_access_lb*u_transfer_energy*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)bitwidth);

            if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                move_front(&m_scheduler->weight_offset_pe);
            }
        }
#endif
        /* Stats */

        // Increase the weight index that finds the offset of weight in local buffer.
        weight_index++;
        // Weight exists in MAC unit and request weight is not required.
        exist_data_mac[data_type_t::WEIGHT] = true, request_to_lb[data_type_t::WEIGHT] = false;

        if(tile_size_mac[data_type_t::WEIGHT] == tile_size_lb[data_type_t::WEIGHT]) { skip_transfer[data_type_t::WEIGHT] = true;}
    }
    // Transfer output data from local buffer to MAC unit.
    if(request_to_lb[data_type_t::OUTPUT]) {
        // Load output data from local buffer to MAC unit.
        if(m_scheduler->output_read_pe[m_scheduler->output_offset_pe.front()]) {
            if(!skip_transfer[data_type_t::OUTPUT]) {
#ifdef FUNCTIONAL

                // Output data transfer 
                m_scheduler->transfer_data(output_data_mac, output_data_lb, 0, m_scheduler->output_offset_pe.front(),
                                           component_type_t::MAC, component_type_t::PE, 
                                           data_type_t::OUTPUT, get_mac_stationary_type(), action_type_t::LOAD);
                // Update for NPUsim ver2
                //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 && 
                //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
                //m_scheduler->transfer_data_ver2(output_data_mac, output_data_lb, 
                //                                component_type_t::MAC, component_type_t::PE, 
                //                                data_type_t::OUTPUT, get_mac_stationary_type(), action_type_t::LOAD, last_component);
#endif

                // Update the number of output data transfer from local buffer to MAC unit
                num_data_transfer_to_mac[data_type_t::OUTPUT]++;

                std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
                parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

                uint64_t address_mac = 0, address_lb = 0;
                unsigned num_access_mac = 0, num_access_lb = 0;
                for(unsigned b = 0; b < parameters_mac[parameter_type_t::BATCH_SIZE]; b++) {
                    for(unsigned k = 0; k < parameters_mac[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                        for(unsigned p = 0; p < parameters_mac[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                            for(unsigned q = 0; q < parameters_mac[parameter_type_t::OUTPUT_WIDTH]; q++) {
                                if(address_mac != ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                               *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                               *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                              k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                               *parameters_mac[parameter_type_t::OUTPUT_WIDTH] + 
                                                                              p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                              mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT]) {

                                    // Update write cost of MAC unit.
                                    access_energy_mac[data_type_t::OUTPUT] += u_write_energy_mac[data_type_t::OUTPUT];
                                    access_cycle_mac[data_type_t::OUTPUT] += u_write_cycle_mac[data_type_t::OUTPUT];
                                    num_access_mac++;

                                    // Update MAC address
                                    address_mac = ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                               *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                               *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                              k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                               *parameters_mac[parameter_type_t::OUTPUT_WIDTH] + 
                                                                              p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                              mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT];
                                }

                                if(address_lb != ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                            b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                             *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                             *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                            k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                             *parameters_lb[parameter_type_t::OUTPUT_WIDTH] + 
                                                                            p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                            mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                    // Update read cost of Local buffer
                                    access_energy_lb[data_type_t::OUTPUT] += u_read_energy_lb[data_type_t::OUTPUT];
                                    access_cycle_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT];
                                    num_access_lb++;

                                    // Update local buffer address
                                    address_lb = ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() + 
                                                                            b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                             *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                             *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                            k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                             *parameters_lb[parameter_type_t::OUTPUT_WIDTH] + 
                                                                            p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                            mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];

                                }
                            }
                        }
                    }
                }

                // Update overlapped cycle at local buffer and MAC units.
                unsigned ratio = ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(line_size_mac[data_type_t::OUTPUT]));

                // at the 1, 2, before last, and last stages
                unsigned first_stage = u_read_cycle_lb[data_type_t::OUTPUT];
                unsigned second_stage = std::max(u_read_cycle_lb[data_type_t::OUTPUT],
                                                 u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
                unsigned last_before_stage = std::max(ratio*u_write_cycle_mac[data_type_t::OUTPUT],
                                                      u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
                unsigned last_stage = ratio*u_write_cycle_mac[data_type_t::OUTPUT];

                // Remainder stages
                unsigned other_stage = std::max(u_read_cycle_lb[data_type_t::OUTPUT],
                                       std::max(ratio*u_write_cycle_mac[data_type_t::OUTPUT],
                                                u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth))));

                if(num_access_lb == 1) {
                    cycle_mac_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT] +
                                                         u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)) +
                                                         ratio*u_write_cycle_mac[data_type_t::OUTPUT];
                } else {
                    cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage +
                                                         (num_access_lb-2)*other_stage +
                                                         last_before_stage + last_stage;
                }

                // Update data transfer cycle and energy between MAC and local buffer.
                transfer_cycle[data_type_t::OUTPUT] += num_access_lb*u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth));
                transfer_energy[data_type_t::OUTPUT] += num_access_lb*u_transfer_energy*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth));
            }

            if(tile_size_mac[data_type_t::OUTPUT] == tile_size_lb[data_type_t::OUTPUT]) { skip_transfer[data_type_t::OUTPUT] = true;}
            /* Stats */
        }
        // Initialize the MAC unit
        else {
            m_scheduler->output_read_pe[m_scheduler->output_offset_pe.front()] = true;
        }

        // Increment the output data index that finds the offset of output data in local buffer.
        output_index++;
        // Output data exists in MAC unit and request output data is not required.
        exist_data_mac[data_type_t::OUTPUT]  = true, request_to_lb[data_type_t::OUTPUT] = false;
    }

    // Execute MAC operation.
    computation(m_scheduler);
   
    // Flush the data in local buffer.
    // When all the data in local buffer are used.
    if(input_index == m_scheduler->offset_size_pe[data_type_t::INPUT].front() &&
       weight_index == m_scheduler->offset_size_pe[data_type_t::WEIGHT].front() &&
       output_index == m_scheduler->offset_size_pe[data_type_t::OUTPUT].front()) {
        flush_data(m_scheduler);
    }
}

// Flush the data in local buffer.
void pe_t::flush_data(scheduler_t *m_scheduler) {
    // Case 1. Input stationary
    if(stationary_type_local_buffer == stationary_type_t::INPUT_STATIONARY) {
        // Case 1) Reuse input data && Request weight and output data
        if(weight_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front() - 1 &&
           output_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front() - 1) {
            // Weight and output data do not exist in local buffer.
            exist_data_lb[data_type_t::WEIGHT] = false, exist_data_lb[data_type_t::OUTPUT] = false;

#ifdef FUNCTIONAL
            // Output data transfer 
            m_scheduler->transfer_data(pe_array->output_data, output_data_lb, m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()], 0,
                                       component_type_t::PE_Y, component_type_t::PE, 
                                       data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
            // Update for NPUsim ver2
            //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 && 
            //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
            //m_scheduler->transfer_data_ver2(pe_array->output_data, output_data_lb, 
            //                                component_type_t::PE_Y, component_type_t::PE,
            //                                data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

            uint64_t address_pe = 0, address_pe_array = 0;
            unsigned num_access_pe = 0, num_access_pe_array = 0;
            for(unsigned b = 0; b < parameters_pe[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_pe[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_pe[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_pe[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            if(address_pe != ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                // Update read cost of Local buffer
                                access_energy_lb[data_type_t::OUTPUT] += u_read_energy_lb[data_type_t::OUTPUT];
                                access_cycle_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT];
                                num_access_pe++;

                                // Update local buffer address
                                address_pe = ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];
                            }

                            if(address_pe_array != ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT]) {

                                // Update write cost of accumulator in PE array.
                                pe_array->access_energy[data_type_t::OUTPUT] += pe_array->u_write_energy[data_type_t::OUTPUT];
                                pe_array->access_cycle[data_type_t::OUTPUT] += pe_array->u_write_cycle[data_type_t::OUTPUT];
                                num_access_pe_array++;

                                // Update accumulator address
                                address_pe_array = ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }
            
            // Increase flush counter of weight and output data.
            weight_flush_counter++;
            output_flush_counter++;

            // Waiting for weight and output data.
            // PE is in idle state.
            wait_data();
        }
        // Case 2) Request all data types 
        else {
            // Input data, weight, and output data do not exist in local buffer.
            exist_data_lb[data_type_t::INPUT] = false; exist_data_lb[data_type_t::WEIGHT] = false; exist_data_lb[data_type_t::OUTPUT] = false;

#ifdef FUNCTIONAL
            // Write back output data
            m_scheduler->transfer_data(pe_array->output_data, output_data_lb, m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()], 0,
                                       component_type_t::PE_Y, component_type_t::PE, 
                                       data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
            // Update for NPUsim ver2
            //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 && 
            //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
            //m_scheduler->transfer_data_ver2(pe_array->output_data, output_data_lb, 
            //                                component_type_t::PE_Y, component_type_t::PE,
            //                                data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

            uint64_t address_pe = 0, address_pe_array = 0;
            unsigned num_access_pe = 0, num_access_pe_array = 0;
            for(unsigned b = 0; b < parameters_pe[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_pe[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_pe[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_pe[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            if(address_pe != ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                // Update read cost of Local buffer
                                access_energy_lb[data_type_t::OUTPUT] += u_read_energy_lb[data_type_t::OUTPUT];
                                access_cycle_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT];
                                num_access_pe++;

                                // Update local buffer address
                                address_pe = ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];
                            }

                            if(address_pe_array != ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT]) {

                                // Update write cost of accumulator in PE array.
                                pe_array->access_energy[data_type_t::OUTPUT] += pe_array->u_write_energy[data_type_t::OUTPUT];
                                pe_array->access_cycle[data_type_t::OUTPUT] += pe_array->u_write_cycle[data_type_t::OUTPUT];
                                num_access_pe_array++;

                                // Update accumulator address
                                address_pe_array = ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }
            /*
            // Update local buffer read cycle and energy.
            access_cycle_lb[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*u_read_cycle_lb[data_type_t::OUTPUT]/(line_size_lb[data_type_t::OUTPUT]/8/sizeof(data_t));
            write_back_cycle_lb += tile_size_lb[data_type_t::OUTPUT]*u_read_cycle_lb[data_type_t::OUTPUT]/(line_size_lb[data_type_t::OUTPUT]/8/sizeof(data_t));
            access_energy_lb[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*u_read_energy_lb[data_type_t::OUTPUT]/(line_size_lb[data_type_t::OUTPUT]/8/sizeof(data_t));

            // Update PE array write cycle and energy.
            pe_array->access_cycle[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*pe_array->u_write_cycle[data_type_t::OUTPUT]/(pe_array->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            pe_array->write_back_cycle += tile_size_lb[data_type_t::OUTPUT]*pe_array->u_write_cycle[data_type_t::OUTPUT]/(pe_array->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            pe_array->access_energy[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*pe_array->u_write_energy[data_type_t::OUTPUT]/(pe_array->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            // Update data transfer cycle between PE and PE array.
            */

            // Set flush counter of weight and output data as zero.
            weight_flush_counter = 0;
            output_flush_counter = 0;

            // Waiting for input data, weight, and output data.
            // PE is in idle state.
            wait_data();
        }
    }
    // Case 2. Weight stationary
    else if(stationary_type_local_buffer == stationary_type_t::WEIGHT_STATIONARY) {
        // Case 1) Reuse weight && Request new input and output data.
        if(input_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front() - 1 &&
           output_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front() - 1) {
            // Input data and output data do not exist in local buffer.
            exist_data_lb[data_type_t::INPUT] = false, exist_data_lb[data_type_t::OUTPUT] = false;

#ifdef FUNCTIONAL
            // Write back output data
            m_scheduler->transfer_data(pe_array->output_data, output_data_lb, m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()], 0,
                                       component_type_t::PE_Y, component_type_t::PE, 
                                       data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
            // Update for NPUsim ver2
            //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 && 
            //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
            //m_scheduler->transfer_data_ver2(pe_array->output_data, output_data_lb, 
            //                                component_type_t::PE_Y, component_type_t::PE,
            //                                data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

            uint64_t address_pe = 0, address_pe_array = 0;
            unsigned num_access_pe = 0, num_access_pe_array = 0;
            for(unsigned b = 0; b < parameters_pe[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_pe[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_pe[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_pe[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            if(address_pe != ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                // Update read cost of Local buffer
                                access_energy_lb[data_type_t::OUTPUT] += u_read_energy_lb[data_type_t::OUTPUT];
                                access_cycle_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT];
                                num_access_pe++;

                                // Update local buffer address
                                address_pe = ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];
                            }

                            if(address_pe_array != ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT]) {

                                // Update write cost of accumulator in PE array.
                                pe_array->access_energy[data_type_t::OUTPUT] += pe_array->u_write_energy[data_type_t::OUTPUT];
                                pe_array->access_cycle[data_type_t::OUTPUT] += pe_array->u_write_cycle[data_type_t::OUTPUT];
                                num_access_pe_array++;

                                // Update accumulator address
                                address_pe_array = ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }

            /* TODO */
            // Update overlapped cycle between the local buffer and PE array.
            // Overlap read cycle of local buffer, write_back_cycle, and write cycle at the accumulator.
            access_cycle_lb[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*u_read_cycle_lb[data_type_t::OUTPUT]/(line_size_lb[data_type_t::OUTPUT]/8/sizeof(data_t));
            write_back_cycle_lb += tile_size_lb[data_type_t::OUTPUT]*u_read_cycle_lb[data_type_t::OUTPUT]/(line_size_lb[data_type_t::OUTPUT]/8/sizeof(data_t));
            access_energy_lb[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*u_read_energy_lb[data_type_t::OUTPUT]/(line_size_lb[data_type_t::OUTPUT]/8/sizeof(data_t));

            // Update PE array write cycle and energy.
            pe_array->access_cycle[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*pe_array->u_write_cycle[data_type_t::OUTPUT]/(pe_array->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            pe_array->write_back_cycle += tile_size_lb[data_type_t::OUTPUT]*pe_array->u_write_cycle[data_type_t::OUTPUT]/(pe_array->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            pe_array->access_energy[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*pe_array->u_write_energy[data_type_t::OUTPUT]/(pe_array->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            // Update data transfer cycle between PE and PE array.

            // Increase flush counter of input data and output data.
            input_flush_counter++;
            output_flush_counter++;

            // Waiting for input and output data.
            // PE is in idle state.
            wait_data();
        }
        // Case 2) Request all data types
        else {
            // Input data, weight, and output data do not exist in local buffer.
            exist_data_lb[data_type_t::INPUT] = false, exist_data_lb[data_type_t::WEIGHT] = false, exist_data_lb[data_type_t::OUTPUT] = false;
#ifdef FUNCTIONAL
            // Write back output data
            m_scheduler->transfer_data(pe_array->output_data, output_data_lb, m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()], 0,
                                       component_type_t::PE_Y, component_type_t::PE, 
                                       data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
            // Update for NPUsim ver2
            //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 && 
            //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
            //m_scheduler->transfer_data_ver2(pe_array->output_data, output_data_lb, 
            //                                component_type_t::PE_Y, component_type_t::PE,
            //                                data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

            uint64_t address_pe = 0, address_pe_array = 0;
            unsigned num_access_pe = 0, num_access_pe_array = 0;
            for(unsigned b = 0; b < parameters_pe[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_pe[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_pe[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_pe[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            if(address_pe != ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                // Update read cost of Local buffer
                                access_energy_lb[data_type_t::OUTPUT] += u_read_energy_lb[data_type_t::OUTPUT];
                                access_cycle_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT];
                                num_access_pe++;

                                // Update local buffer address
                                address_pe = ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];
                            }

                            if(address_pe_array != ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT]) {

                                // Update write cost of accumulator in PE array.
                                pe_array->access_energy[data_type_t::OUTPUT] += pe_array->u_write_energy[data_type_t::OUTPUT];
                                pe_array->access_cycle[data_type_t::OUTPUT] += pe_array->u_write_cycle[data_type_t::OUTPUT];
                                num_access_pe_array++;

                                // Update accumulator address
                                address_pe_array = ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }

            /*
            // Update stats.
            // Update local buffer read cycle and energy.
            access_cycle_lb[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*u_read_cycle_lb[data_type_t::OUTPUT]/(line_size_lb[data_type_t::OUTPUT]/8/sizeof(data_t));
            write_back_cycle_lb += tile_size_lb[data_type_t::OUTPUT]*u_read_cycle_lb[data_type_t::OUTPUT]/(line_size_lb[data_type_t::OUTPUT]/8/sizeof(data_t));
            access_energy_lb[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*u_read_energy_lb[data_type_t::OUTPUT]/(line_size_lb[data_type_t::OUTPUT]/8/sizeof(data_t));

            // Update PE array write cycle and energy.
            pe_array->access_cycle[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*pe_array->u_write_cycle[data_type_t::OUTPUT]/(pe_array->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            pe_array->write_back_cycle += tile_size_lb[data_type_t::OUTPUT]*pe_array->u_write_cycle[data_type_t::OUTPUT]/(pe_array->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            pe_array->access_energy[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*pe_array->u_write_energy[data_type_t::OUTPUT]/(pe_array->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            */

            // Set flush counter of input data and output data as zero.
            input_flush_counter = 0;
            output_flush_counter = 0;

            // Waiting for input data, weight, and output data.
            // PE is in idle state.
            wait_data();
        }
    }
    // Case 3. Output stationary
    else if(stationary_type_local_buffer == stationary_type_t::OUTPUT_STATIONARY) {
        // Case 1) Reuse output data && Request input data and weight
        if(input_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front() - 1 &&
           weight_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front() - 1) {
            // input data and weight do not exist in local buffer.
            exist_data_lb[data_type_t::INPUT] = false, exist_data_lb[data_type_t::WEIGHT] = false;

            // Increase flush counter of input data and weight.
            input_flush_counter++;
            weight_flush_counter++;

            // Waiting for input data and weight.
            // PE is in idle state
            wait_data();
        }
        // Case 2) Request all data types
        else {
            // Input data, weight, and output data are not in local buffer.
            exist_data_lb[data_type_t::INPUT] = false, exist_data_lb[data_type_t::WEIGHT] = false, exist_data_lb[data_type_t::OUTPUT] = false;
#ifdef FUNCTIONAL
            // Write back output data
            m_scheduler->transfer_data(pe_array->output_data, output_data_lb, m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()], 0,
                                       component_type_t::PE_Y, component_type_t::PE, 
                                       data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
            // Update for NPUsim ver2
            //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 && 
            //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
            //m_scheduler->transfer_data_ver2(pe_array->output_data, output_data_lb, 
            //                                component_type_t::PE_Y, component_type_t::PE,
            //                                data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

            uint64_t address_pe = 0, address_pe_array = 0;
            unsigned num_access_pe = 0, num_access_pe_array = 0;
            for(unsigned b = 0; b < parameters_pe[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_pe[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_pe[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_pe[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            if(address_pe != ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                // Update read cost of Local buffer
                                access_energy_lb[data_type_t::OUTPUT] += u_read_energy_lb[data_type_t::OUTPUT];
                                access_cycle_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT];
                                num_access_pe++;

                                // Update local buffer address
                                address_pe = ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];
                            }

                            if(address_pe_array != ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT]) {

                                // Update write cost of accumulator in PE array.
                                pe_array->access_energy[data_type_t::OUTPUT] += pe_array->u_write_energy[data_type_t::OUTPUT];
                                pe_array->access_cycle[data_type_t::OUTPUT] += pe_array->u_write_cycle[data_type_t::OUTPUT];
                                num_access_pe_array++;

                                // Update accumulator address
                                address_pe_array = ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }

            /*
            // Update local buffer read cycle and energy.
            access_cycle_lb[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*u_read_cycle_lb[data_type_t::OUTPUT]/(line_size_lb[data_type_t::OUTPUT]/8/sizeof(data_t));
            write_back_cycle_lb += tile_size_lb[data_type_t::OUTPUT]*u_read_cycle_lb[data_type_t::OUTPUT]/(line_size_lb[data_type_t::OUTPUT]/8/sizeof(data_t));
            access_energy_lb[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*u_read_energy_lb[data_type_t::OUTPUT]/(line_size_lb[data_type_t::OUTPUT]/8/sizeof(data_t));

            // Update PE array write cycle and energy.
            pe_array->access_cycle[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*pe_array->u_write_cycle[data_type_t::OUTPUT]/(pe_array->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            pe_array->write_back_cycle += tile_size_lb[data_type_t::OUTPUT]*pe_array->u_write_cycle[data_type_t::OUTPUT]/(pe_array->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            pe_array->access_energy[data_type_t::OUTPUT] += tile_size_lb[data_type_t::OUTPUT]*pe_array->u_write_energy[data_type_t::OUTPUT]/(pe_array->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            */

            // Set flush counter of input data and weight to zero.
            input_flush_counter = 0;
            weight_flush_counter = 0;

            // Waiting for input data, weight, and output data.
            // PE is in idle state.
            wait_data();
        }
    }
    // The counter of local buffer should be initialized to 0.
    input_index = 0, weight_index = 0, output_index = 0;
    //if(index == pe_array->get_number_of_active_pes()-1) {
    //    pe_array->update_pe_stats();
    //}
}

void pe_t::mac_operation() {
#if defined(USER_INTEGER) || defined(USER_FLOAT)
    for(unsigned i = 0; i < num_active_macs; i++) {
        output_data_mac[i].value += input_data_mac[i].value*weight_mac[i].value;
    }
#else
    for(unsigned i = 0; i < num_active_macs; i++) {
        output_data_mac[i] += input_data_mac[i]*weight_mac[i];
    }
#endif
}

void pe_t::activation() {
#if defined(USER_INTEGER) || defined(USER_FLOAT)
    for(unsigned i = 0; i < num_active_macs; i++) {
        output_data_mac[i].value = output_data_mac[i].value > 0.0 ? output_data_mac[i].value : 0.0;
        //output_data_mac[i] = output_data_mac[i] > 0.0 ? output_data_mac[i] : 0.0;
    }
#else
    for(unsigned i = 0; i < num_active_macs; i++) {
        output_data_mac[i] = output_data_mac[i] > 0.0 ? output_data_mac[i] : 0.0;
    }
#endif
}

/* TODO */
void pe_t::max_pooling() {

}

void pe_t::avg_pooling() {

}

            
// Print out the stats of the PE.
void pe_t::print_specification() {

    std::cout << std::fixed;
    std::cout << "============ MAC specification =============" << std::endl;
    std::cout << "# of MACs          :" << std::setw(24) 
                                        << num_macs*mac_width << std::endl;
    std::cout << "Computation energy :" << std::setw(21) << std::setprecision(2) 
                                        << u_computation_energy << " pJ" << std::endl;
    std::cout << "Computation cycle  :" << std::setw(17) << std::setprecision(1) 
                                        << u_computation_cycle << " cycles" << std::endl;

    std::cout << "Line size" << std::endl;
    std::cout << " * Input           :" << std::setw(19) << std::setprecision(0) 
                                        << line_size_mac[data_type_t::INPUT] << " bits" << std::endl;
    std::cout << " * Weight          :" << std::setw(19) << std::setprecision(0) 
                                        << line_size_mac[data_type_t::WEIGHT] << " bits" << std::endl;
    std::cout << " * Output          :" << std::setw(19) << std::setprecision(0) 
                                        << line_size_mac[data_type_t::OUTPUT] << " bits" << std::endl;

    std::cout << "Access energy (read/write)" << std::endl;
    std::cout << " * MAC units       :" << std::setw(16) << std::setprecision(2) 
                                        << u_read_energy_mac[data_type_t::INPUT] << "/" << u_write_energy_mac[data_type_t::INPUT] << " pJ" << std::endl;
    std::cout << "Access cycle (read/write)" << std::endl;
    std::cout << " * MAC units       :" << std::setw(13) << std::setprecision(1) 
                                        << u_read_cycle_mac[data_type_t::INPUT] << "/" << u_write_cycle_mac[data_type_t::INPUT] << " cycles" << std::endl;
    std::cout << std::endl;

    std::cout << "============= PE specification =============" << std::endl;
    // Print the size of local buffer.
    if(memory_type == memory_type_t::SEPARATE) {
        std::cout << "Local buffer type  :" << std::setw(24) 
                                            << "Separated" << std::endl;
        std::cout << " * Input           :" << std::setw(19) 
                                            << input_size << " Byte" << std::endl;
        std::cout << " * Weight          :" << std::setw(19) 
                                            << weight_size << " Byte" << std::endl;
        std::cout << " * Output          :" << std::setw(19) 
                                            << output_size << " Byte" << std::endl;
        
        std::cout << "Line size" << std::endl;
        std::cout << " * Input           :" << std::setw(19) << std::setprecision(0) 
                                            << line_size_lb[data_type_t::INPUT] << " bits" << std::endl;
        std::cout << " * Weight          :" << std::setw(19) << std::setprecision(0) 
                                            << line_size_lb[data_type_t::WEIGHT] << " bits" << std::endl;
        std::cout << " * Output          :" << std::setw(19) << std::setprecision(0) 
                                            << line_size_lb[data_type_t::OUTPUT] << " bits" << std::endl;
    }
    else if(memory_type == memory_type_t::SHARED) {
        std::cout << "Local buffer type  :" << std::setw(24) 
                                            << "Shared" << std::endl;
        std::cout << "Buffer size        :" << std::setw(19) 
                                            << input_size + weight_size + output_size << " Byte" << std::endl;
        std::cout << "Line size          :" << std::setw(19) << std::setprecision(0) 
                                            << line_size_lb[data_type_t::INPUT] << " bits" << std::endl;
    }
    // Print the stationary type.
    if(stationary_type_local_buffer == stationary_type_t::INPUT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24) 
                                            << "Input stationary" << std::endl;
    }
    else if(stationary_type_local_buffer == stationary_type_t::WEIGHT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24) 
                                            << "Weight stationary" << std::endl;
    } 
    else if(stationary_type_local_buffer == stationary_type_t::OUTPUT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24) 
                                            << "Output stationary" << std::endl;
    }
    // Print out the bandwidth
    std::cout << "Bandwidth          :" << std::setw(19) << std::setprecision(0) 
                                        << bandwidth << " GB/s" << std::endl;
    // Print out unit access cost of the local buffer. 
    if(memory_type == memory_type_t::SEPARATE) {
        std::cout << "Access energy (read/write)" << std::endl;
        std::cout << " * Input buffer    :" << std::setw(16) << std::setprecision(2) 
                                            << u_read_energy_lb[data_type_t::INPUT]  << "/" << u_write_energy_lb[data_type_t::INPUT]  << " pJ" << std::endl;
        std::cout << " * Weight buffer   :" << std::setw(16) << std::setprecision(2) 
                                            << u_read_energy_lb[data_type_t::WEIGHT] << "/" << u_write_energy_lb[data_type_t::WEIGHT] << " pJ" << std::endl;
        std::cout << " * Output buffer   :" << std::setw(16) << std::setprecision(2) 
                                            << u_read_energy_lb[data_type_t::OUTPUT] << "/" << u_write_energy_lb[data_type_t::OUTPUT] << " pJ" << std::endl;

        std::cout << "Access cycle (read/write)" << std::endl;
        std::cout << " * Input buffer    :" << std::setw(13) << std::setprecision(1) 
                                            << u_read_cycle_lb[data_type_t::INPUT] << "/" << u_write_cycle_lb[data_type_t::INPUT] << " cycles" << std::endl;
        std::cout << " * Weight buffer   :" << std::setw(13) << std::setprecision(1) 
                                            << u_read_cycle_lb[data_type_t::WEIGHT] << "/" << u_write_cycle_lb[data_type_t::WEIGHT] << " cycles" << std::endl;
        std::cout << " * Output buffer   :" << std::setw(13) << std::setprecision(1) 
                                            << u_read_cycle_lb[data_type_t::OUTPUT] << "/" << u_write_cycle_lb[data_type_t::OUTPUT] << " cycles" << std::endl;
    }
    else if(memory_type == memory_type_t::SHARED) {
        std::cout << "Access energy (read/write)" << std::endl;
        std::cout << " * Local buffer    :" << std::setw(16) << std::setprecision(2) 
                                            << u_read_energy_lb[data_type_t::INPUT]  << "/" << u_write_energy_lb[data_type_t::INPUT]  << " pJ" << std::endl;

        std::cout << "Access cycle (read/write)" << std::endl;
        std::cout << " * Local buffer    :" << std::setw(13) << std::setprecision(1) 
                                            << u_read_cycle_lb[data_type_t::INPUT] << "/" << u_write_cycle_lb[data_type_t::INPUT] << " cycles" << std::endl;
    }
	std::cout << std::endl;
}

// Reset the stats
void pe_t::reset() {
    
    memset(input_data_mac, 0.0, mac_width*num_macs*sizeof(data_t));
    memset(weight_mac, 0.0, mac_width*num_macs*sizeof(data_t));
    memset(output_data_mac, 0.0, num_macs*sizeof(data_t));

    memset(input_data_lb, 0.0, input_size);
    memset(weight_lb, 0.0, weight_size);
    memset(output_data_lb, 0.0, output_size);

    idle = false;

    num_computation = 0;

    computation_cycle = 0;
    computation_energy = 0;

    write_back_cycle_mac  = 0.0;
    write_back_cycle_lb = 0.0;
    overlapped_transfer_cycle = 0.0;

    utilization_mac = 0.0;

    num_request_to_lb.assign(data_type_t::NUM_DATA_TYPES, 0);
    num_data_transfer_to_mac.assign(data_type_t::NUM_DATA_TYPES, 0);
    
    // Reset access cycle and energy of MAC unit.
    access_cycle_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    access_energy_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Reset access cycle and energy of local buffer.
    access_cycle_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    access_energy_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    
    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Reset overlapped cycle between MAC unit and local buffer
    cycle_mac_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    
    utilization_local_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);

}

// All functions of undefined stationary
undefined_stationary_t::undefined_stationary_t(section_config_t m_section_config) :
    pe_t(m_section_config) {

}

undefined_stationary_t::~undefined_stationary_t() {

}

void undefined_stationary_t::computation(scheduler_t *m_scheduler) {

    if(exist_data_mac[data_type_t::INPUT] && exist_data_mac[data_type_t::WEIGHT] && exist_data_mac[data_type_t::OUTPUT]) {
#ifdef FUNCTIONAL

#endif

        num_computation++;
        computation_cycle += u_computation_cycle;
        computation_energy += u_computation_energy;

        exist_data_mac[data_type_t::INPUT] = false, request_to_lb[data_type_t::INPUT] = true;
        exist_data_mac[data_type_t::WEIGHT] = false, request_to_lb[data_type_t::WEIGHT] = true;
        exist_data_mac[data_type_t::OUTPUT] = false, request_to_lb[data_type_t::OUTPUT] = true;

        // Update stats (Memory operation)
        /* Stats */
        num_request_to_lb[data_type_t::INPUT]++;
        num_request_to_lb[data_type_t::WEIGHT]++;
        num_request_to_lb[data_type_t::OUTPUT]++;
        /* Stats */

#ifdef FUNCTIONAL
        // Write back output data
        m_scheduler->transfer_data(output_data_lb, output_data_mac, m_scheduler->output_offset_pe.front(), 0, 
                                   component_type_t::PE, component_type_t::MAC, 
                                   data_type_t::OUTPUT, get_mac_stationary_type(), action_type_t::STORE);
        // Update for NPUsim ver2
        //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 && 
        //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
        //m_scheduler->transfer_data_ver2(output_data_lb, output_data_mac,
        //                                component_type_t::PE, component_type_t::MAC,
        //                                data_type_t::OUTPUT, get_mac_stationary_type()), action_type_t::STORE);
#endif
        /* Stats */
        // Update MAC read cycle and energy.
        access_cycle_mac[data_type_t::OUTPUT] += tile_size_mac[data_type_t::OUTPUT]*u_read_cycle_mac[data_type_t::OUTPUT]/(line_size_mac[data_type_t::OUTPUT]/8/sizeof(data_t));
        write_back_cycle_mac += tile_size_mac[data_type_t::OUTPUT]*u_read_cycle_mac[data_type_t::OUTPUT]/(line_size_mac[data_type_t::OUTPUT]/8/sizeof(data_t));
        access_energy_mac[data_type_t::OUTPUT] += tile_size_mac[data_type_t::OUTPUT]*u_read_energy_mac[data_type_t::OUTPUT]/(line_size_mac[data_type_t::OUTPUT]/8/sizeof(data_t));

        // Update local buffer write cycle and energy.
        access_cycle_lb[data_type_t::OUTPUT] += tile_size_mac[data_type_t::OUTPUT]*u_write_cycle_lb[data_type_t::OUTPUT]/(line_size_mac[data_type_t::OUTPUT]/8/sizeof(data_t));
        write_back_cycle_lb += tile_size_mac[data_type_t::OUTPUT]*u_write_cycle_lb[data_type_t::OUTPUT]/(line_size_mac[data_type_t::OUTPUT]/8/sizeof(data_t));
        access_energy_lb[data_type_t::OUTPUT] += tile_size_mac[data_type_t::OUTPUT]*u_write_energy_lb[data_type_t::OUTPUT]/(line_size_mac[data_type_t::OUTPUT]/8/sizeof(data_t));

        // Update data transfer cycle and energy between MAC and local buffer.
        transfer_cycle[data_type_t::OUTPUT] += u_transfer_cycle*ceil((double)(tile_size_mac[data_type_t::OUTPUT]*8*sizeof(data_t))/(double)bitwidth);
        overlapped_transfer_cycle += u_transfer_cycle*ceil((double)(tile_size_mac[data_type_t::OUTPUT]*8*sizeof(data_t))/(double)bitwidth);
        transfer_energy[data_type_t::OUTPUT] += u_transfer_energy*ceil((double)(tile_size_mac[data_type_t::OUTPUT]*8*sizeof(data_t))/(double)bitwidth);
        /* Stats */


        if(input_index == m_scheduler->offset_size_pe[data_type_t::INPUT].front() &&
           weight_index == m_scheduler->offset_size_pe[data_type_t::WEIGHT].front()) {

            input_index = 0, weight_index = 0;
        }
    }
}


input_stationary_t::input_stationary_t(section_config_t m_section_config) :
    pe_t(m_section_config) {

}

input_stationary_t::~input_stationary_t() {

}

// Execute MAC operation 
void input_stationary_t::computation(scheduler_t *m_scheduler) {

    // When all the data exist, execute the MAC operation
    if(exist_data_mac[data_type_t::INPUT] && exist_data_mac[data_type_t::WEIGHT] && exist_data_mac[data_type_t::OUTPUT]) {
        if(m_scheduler->layer_name == layer_name_t::CONVOLUTIONAL_LAYER ||
           m_scheduler->layer_name == layer_name_t::CONNECTED_LAYER) {
#ifdef FUNCTIONAL
            mac_operation();
            activation();
            // Split active_macs and mac_width
            if(m_scheduler->compression_type == compression_type_t::DENSE) {
                for(unsigned i = 0; i < num_active_macs; i++) {
                    num_computation++;
                    computation_energy += u_computation_energy;
                }
                computation_cycle += u_computation_cycle;
            }
            else {
                for(unsigned i = 0; i < num_active_macs; i++) {
                    if(input_data_mac[i] != 0.0 && weight_mac[i] != 0.0) {
                        num_computation++;
                        computation_energy += u_computation_energy;
                    }
                }

                for(unsigned i = 0; i < num_active_macs; i++) {
                    if(input_data_mac[i] != 0.0 && weight_mac[i] != 0.0) {
                        computation_cycle += u_computation_cycle;
                        break;
                    }
                }
            }
#else
            for(unsigned i = 0; i < num_active_macs; i++) {
                num_computation++;
                computation_energy += u_computation_energy;
            }
            computation_cycle += u_computation_cycle;
#endif
        } else if(m_scheduler->layer_name == layer_name_t::MAXPOOL_LAYER) {
            max_pooling();
        } else if(m_scheduler->layer_name == layer_name_t::AVGPOOL_LAYER) {
            avg_pooling();
        }
#ifdef FUNCTIONAL
        // Write back output data
        m_scheduler->transfer_data(output_data_lb, output_data_mac, m_scheduler->output_offset_pe.front(), 0, 
                                   component_type_t::PE, component_type_t::MAC, 
                                   data_type_t::OUTPUT, get_mac_stationary_type(), action_type_t::STORE);
        // Update for NPUsim ver2
        //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 && 
        //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
        //m_scheduler->transfer_data_ver2(output_data_lb, output_data_mac,
        //                                component_type_t::PE, component_type_t::MAC,
        //                                data_type_t::OUTPUT, get_mac_stationary_type()), action_type_t::STORE);
#endif

        std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
        std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);

        parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
        parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

        uint64_t address_mac = 0, address_lb = 0;
        unsigned num_access_mac = 0, num_access_lb = 0;
        for(unsigned b = 0; b < parameters_mac[parameter_type_t::BATCH_SIZE]; b++) {
            for(unsigned k = 0; k < parameters_mac[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                for(unsigned p = 0; p < parameters_mac[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                    for(unsigned q = 0; q < parameters_mac[parameter_type_t::OUTPUT_WIDTH]; q++) {
                        if(address_mac != ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                      mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT]) {

                            // Update write cost of MAC unit.
                            access_energy_mac[data_type_t::OUTPUT] += u_read_energy_mac[data_type_t::OUTPUT];
                            access_cycle_mac[data_type_t::OUTPUT] += u_read_cycle_mac[data_type_t::OUTPUT];
                            num_access_mac++;

                            // Update MAC address
                            address_mac = ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                      mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT];
                        }

                        if(address_lb != ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                    b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                    mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                            // Update read cost of Local buffer
                            access_energy_lb[data_type_t::OUTPUT] += u_write_energy_lb[data_type_t::OUTPUT];
                            access_cycle_lb[data_type_t::OUTPUT] += u_write_cycle_lb[data_type_t::OUTPUT];
                            num_access_lb++;

                            // Update local buffer address
                            address_lb = ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                    b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                    mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];

                        }
                    }
                }
            }
        }

        // Update overlapped cycle at MAC units and local buffer
        unsigned ratio = ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(line_size_mac[data_type_t::OUTPUT]));

        unsigned first_stage = ratio*u_read_cycle_mac[data_type_t::OUTPUT];
        unsigned second_stage = std::max(ratio*u_read_cycle_mac[data_type_t::OUTPUT],
                                         u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
        unsigned last_before_stage = std::max(u_write_cycle_lb[data_type_t::OUTPUT],
                                              u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
        unsigned last_stage = u_write_cycle_lb[data_type_t::OUTPUT];

        unsigned other_stage = std::max(ratio*u_read_cycle_mac[data_type_t::OUTPUT],
                               std::max(u_write_cycle_lb[data_type_t::OUTPUT],
                                        u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth))));

        if(num_access_lb == 1) {
            cycle_mac_lb[data_type_t::OUTPUT] += ratio*u_read_cycle_mac[data_type_t::OUTPUT] +
                                                 u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)) +
                                                 u_write_cycle_lb[data_type_t::OUTPUT];
        } else if(num_access_lb == 2) {
            cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage + last_before_stage + last_stage;
        } else {
            cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage +
                                                 (num_access_lb - 2)*other_stage +
                                                 last_before_stage + last_stage;
        }

        transfer_energy[data_type_t::OUTPUT] += num_access_mac*u_transfer_energy*ceil((double)(line_size_mac[data_type_t::OUTPUT])/(double)(bitwidth));
        transfer_cycle[data_type_t::OUTPUT] += num_access_mac*u_transfer_cycle*ceil((double)(line_size_mac[data_type_t::OUTPUT])/(double)(bitwidth));

        // weight and output data should be requested from MAC to local buffer.
        exist_data_mac[data_type_t::WEIGHT] = false, request_to_lb[data_type_t::WEIGHT] = true;
        exist_data_mac[data_type_t::OUTPUT] = false, request_to_lb[data_type_t::OUTPUT] = true;

        // Update stats (Memory operation)
        num_request_to_lb[data_type_t::WEIGHT]++;
        num_request_to_lb[data_type_t::OUTPUT]++;

        // When all the weight and output data at local buffer are used.
        // A new input tile should be requested from MAC unit to local buffer.
        if(weight_index == m_scheduler->offset_size_pe[data_type_t::WEIGHT].front() &&
           output_index == m_scheduler->offset_size_pe[data_type_t::OUTPUT].front()) {

            // Request new input data.
            exist_data_mac[data_type_t::INPUT] = false, request_to_lb[data_type_t::INPUT] = true;

            // Update stats.
            num_request_to_lb[data_type_t::INPUT]++;

            // When not all input data in local buffer is transferred to MAC unit.
            if(input_index < m_scheduler->offset_size_pe[data_type_t::INPUT].front()) {
                weight_index = 0, output_index = 0;
            }
        }
    }
}

weight_stationary_t::weight_stationary_t(section_config_t m_section_config) :
    pe_t(m_section_config) {

}

weight_stationary_t::~weight_stationary_t() {

}

void weight_stationary_t::computation(scheduler_t *m_scheduler) {

    // Execute MAC operation when all DNN data exist in MAC register.
    if(exist_data_mac[data_type_t::INPUT] && exist_data_mac[data_type_t::WEIGHT] && exist_data_mac[data_type_t::OUTPUT]) {
        if(m_scheduler->layer_name == layer_name_t::CONVOLUTIONAL_LAYER || 
           m_scheduler->layer_name == layer_name_t::CONNECTED_LAYER) {
#ifdef FUNCTIONAL
            mac_operation();
            activation();
            // Split active_macs and mac_width
            if(m_scheduler->compression_type == compression_type_t::DENSE) {
                for(unsigned i = 0; i < num_active_macs; i++) {
                    num_computation++;
                    computation_energy += u_computation_energy;
                }
                computation_cycle += u_computation_cycle;
            }
            else {
                for(unsigned i = 0; i < num_active_macs; i++) {
                    if(input_data_mac[i] != 0.0 && weight_mac[i] != 0.0) {
                        num_computation++;
                        computation_energy += u_computation_energy;
                    }
                }

                for(unsigned i = 0; i < num_active_macs; i++) {
                    if(input_data_mac[i] != 0.0 && weight_mac[i] != 0.0) {
                        computation_cycle += u_computation_cycle;
                        break;
                    }
                }
            }
#else 
            for(unsigned i = 0; i < num_active_macs; i++) {
                num_computation++;
                computation_energy += u_computation_energy;
            }
            computation_cycle += u_computation_cycle;
#endif
        }
        else if(m_scheduler->layer_name == layer_name_t::MAXPOOL_LAYER) {
            max_pooling();
        }
        else if(m_scheduler->layer_name == layer_name_t::AVGPOOL_LAYER) {
            avg_pooling();
        }

#ifdef FUNCTIONAL
        // Write back output data
        m_scheduler->transfer_data(output_data_lb, output_data_mac, m_scheduler->output_offset_pe.front(), 0, 
                                   component_type_t::PE, component_type_t::MAC, 
                                   data_type_t::OUTPUT, get_mac_stationary_type(), action_type_t::STORE);
        // Update for NPUsim ver2
        //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 && 
        //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
        //m_scheduler->transfer_data_ver2(output_data_lb, output_data_mac,
        //                                component_type_t::PE, component_type_t::MAC,
        //                                data_type_t::OUTPUT, get_mac_stationary_type()), action_type_t::STORE);
#endif
        std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
        std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);

        parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
        parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

        uint64_t address_mac = 0, address_lb = 0;
        unsigned num_access_mac = 0, num_access_lb = 0;
        for(unsigned b = 0; b < parameters_mac[parameter_type_t::BATCH_SIZE]; b++) {
            for(unsigned k = 0; k < parameters_mac[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                for(unsigned p = 0; p < parameters_mac[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                    for(unsigned q = 0; q < parameters_mac[parameter_type_t::OUTPUT_WIDTH]; q++) {
                        if(address_mac != ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] + 
                                                                      p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                      mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT]) {

                            // Update write cost of MAC unit.
                            access_energy_mac[data_type_t::OUTPUT] += u_read_energy_mac[data_type_t::OUTPUT];
                            access_cycle_mac[data_type_t::OUTPUT] += u_read_cycle_mac[data_type_t::OUTPUT];
                            num_access_mac++;

                            // Update MAC address
                            address_mac = ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] + 
                                                                      p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                      mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT];
                        }

                        if(address_lb != ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                    b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] + 
                                                                    p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                    mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                            // Update read cost of Local buffer
                            access_energy_lb[data_type_t::OUTPUT] += u_write_energy_lb[data_type_t::OUTPUT];
                            access_cycle_lb[data_type_t::OUTPUT] += u_write_cycle_lb[data_type_t::OUTPUT];
                            num_access_lb++;

                            // Update local buffer address
                            address_lb = ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() + 
                                                                    b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] + 
                                                                    p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                    mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];

                        }
                    }
                }
            }
        }

        // Update overlapped cycle at MAC units and local buffer
        unsigned ratio = ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(line_size_mac[data_type_t::OUTPUT]));

        unsigned first_stage = ratio*u_read_cycle_mac[data_type_t::OUTPUT];
        unsigned second_stage = std::max(ratio*u_read_cycle_mac[data_type_t::OUTPUT],
                                         u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
        unsigned last_before_stage = std::max(u_write_cycle_lb[data_type_t::OUTPUT],
                                              u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
        unsigned last_stage = u_write_cycle_lb[data_type_t::OUTPUT];

        unsigned other_stage = std::max(ratio*u_read_cycle_mac[data_type_t::OUTPUT],
                               std::max(u_write_cycle_lb[data_type_t::OUTPUT],
                                        u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth))));

        if(num_access_lb == 1) {
            cycle_mac_lb[data_type_t::OUTPUT] += ratio*u_read_cycle_mac[data_type_t::OUTPUT] + 
                                                 u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)) +
                                                 u_write_cycle_lb[data_type_t::OUTPUT];
        } else if(num_access_lb == 2) {
            cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage + last_before_stage + last_stage;
        } else {
            cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage +
                                                 (num_access_lb - 2)*other_stage +
                                                 last_before_stage + last_stage;
        }
        
        transfer_energy[data_type_t::OUTPUT] += num_access_mac*u_transfer_energy*ceil((double)(line_size_mac[data_type_t::OUTPUT])/(double)(bitwidth));
        transfer_cycle[data_type_t::OUTPUT] += num_access_mac*u_transfer_cycle*ceil((double)(line_size_mac[data_type_t::OUTPUT])/(double)(bitwidth));

        // Input data and output data should be requested from MAC.
        exist_data_mac[data_type_t::INPUT] = false, request_to_lb[data_type_t::INPUT] = true;
        exist_data_mac[data_type_t::OUTPUT] = false, request_to_lb[data_type_t::OUTPUT] =true;

        // Update stats (Memory operation)
        num_request_to_lb[data_type_t::INPUT]++;
        num_request_to_lb[data_type_t::OUTPUT]++;

        // When all the input and output data at local buffer are used.
        // A new weight tile should be request from MAC to local buffer.
        if(input_index == m_scheduler->offset_size_pe[data_type_t::INPUT].front() &&
           output_index == m_scheduler->offset_size_pe[data_type_t::OUTPUT].front()) {

            // Request new weight.
            exist_data_mac[data_type_t::WEIGHT] = false, request_to_lb[data_type_t::WEIGHT] = true;

            // Update stats of memory operation.
            /* Stats */
            num_request_to_lb[data_type_t::WEIGHT]++;
            /* Stats */

            // The case when not all weight in local buffer is transferred to MAC.
            if(weight_index < m_scheduler->offset_size_pe[data_type_t::WEIGHT].front()) {
                input_index = 0, output_index = 0;
            }
        }
    }
}

output_stationary_t::output_stationary_t(section_config_t m_section_config) :
    pe_t(m_section_config) {

}

output_stationary_t::~output_stationary_t() {

}

void output_stationary_t::computation(scheduler_t *m_scheduler) {

    // Execute MAC operation
    // When all data exist.
    if(exist_data_mac[data_type_t::INPUT] && exist_data_mac[data_type_t::WEIGHT] && exist_data_mac[data_type_t::OUTPUT]) {
        if(m_scheduler->layer_name == layer_name_t::CONVOLUTIONAL_LAYER ||
           m_scheduler->layer_name == layer_name_t::CONNECTED_LAYER) {
#ifdef FUNCTIONAL
            mac_operation();
            activation();
            if(m_scheduler->compression_type == compression_type_t::DENSE) {
                for(unsigned i = 0; i < num_active_macs; i++) {
                    num_computation++;
                    computation_energy += u_computation_cycle;
                }
            }
            else {
                for(unsigned i = 0; i < num_active_macs; i++) {
                    if(input_data_mac[i] != 0.0 && weight_mac[i] != 0.0) {
                        num_computation++;
                        computation_energy += u_computation_cycle;
                    }
                }

                for(unsigned i = 0; i < num_active_macs; i++) {
                    if(input_data_mac[i] != 0.0 && weight_mac[i] != 0.0) {
                        computation_cycle += u_computation_cycle;
                        break;
                    }
                }
            }
#else
            for(unsigned i = 0; i < num_active_macs; i++) {
                num_computation++;
                computation_energy += u_computation_energy;
            }
            computation_cycle += u_computation_cycle;
            
#endif
        }

        else if(m_scheduler->layer_name == layer_name_t::MAXPOOL_LAYER) {
            max_pooling();
        }
        else if(m_scheduler->layer_name == layer_name_t::AVGPOOL_LAYER) {
            avg_pooling();
        }


        // After computation.
        // Input data and weight should be request from MAC to local buffer.
        exist_data_mac[data_type_t::INPUT] = false, request_to_lb[data_type_t::INPUT] = true;
        exist_data_mac[data_type_t::WEIGHT] = false, request_to_lb[data_type_t::WEIGHT] = true;

        // Update stats (memory operation).
        /* Stats */
        num_request_to_lb[data_type_t::INPUT]++;
        num_request_to_lb[data_type_t::WEIGHT]++;
        /* Stats */

        // When all the input and weight data at local buffer are used.
        // A new output data tile should be request from MAC to local buffer.
        if(input_index == m_scheduler->offset_size_pe[data_type_t::INPUT].front() &&
           weight_index == m_scheduler->offset_size_pe[data_type_t::WEIGHT].front()) {

            // Request new output data.
            exist_data_mac[data_type_t::OUTPUT] = false, request_to_lb[data_type_t::OUTPUT] = true;

            num_request_to_lb[data_type_t::OUTPUT]++;


#ifdef FUNCTIONAL
            // Write back output data
            m_scheduler->transfer_data(output_data_lb, output_data_mac, m_scheduler->output_offset_pe.front(), 0, 
                                       component_type_t::PE, component_type_t::MAC, 
                                       data_type_t::OUTPUT, get_mac_stationary_type(), action_type_t::STORE);
        // Update for NPUsim ver2
        //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 && 
        //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
        //m_scheduler->transfer_data_ver2(output_data_lb, output_data_mac,
        //                                component_type_t::PE, component_type_t::MAC,
        //                                data_type_t::OUTPUT, get_mac_stationary_type()), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
            parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

            uint64_t address_mac = 0, address_lb = 0;
            unsigned num_access_mac = 0, num_access_lb = 0;
            for(unsigned b = 0; b < parameters_mac[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_mac[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_mac[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_mac[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            if(address_mac != ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                           *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                           *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                          k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                           *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                          p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                          mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT]) {

                                // Update write cost of MAC unit.
                                access_energy_mac[data_type_t::OUTPUT] += u_read_energy_mac[data_type_t::OUTPUT];
                                access_cycle_mac[data_type_t::OUTPUT] += u_read_cycle_mac[data_type_t::OUTPUT];
                                num_access_mac++;

                                // Update MAC address
                                address_mac = ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                           *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                           *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                          k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                           *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                          p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                          mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT];
                            }

                            if(address_lb != ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                        b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                // Update read cost of Local buffer
                                access_energy_lb[data_type_t::OUTPUT] += u_write_energy_lb[data_type_t::OUTPUT];
                                access_cycle_lb[data_type_t::OUTPUT] += u_write_cycle_lb[data_type_t::OUTPUT];
                                num_access_lb++;

                                // Update local buffer address
                                address_lb = ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                        b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];

                            }
                        }
                    }
                }
            }

            // Update overlapped cycle at MAC units and local buffer
            unsigned ratio = ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(line_size_mac[data_type_t::OUTPUT]));

            unsigned first_stage = ratio*u_read_cycle_mac[data_type_t::OUTPUT];
            unsigned second_stage = std::max(ratio*u_read_cycle_mac[data_type_t::OUTPUT],
                                             u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
            unsigned last_before_stage = std::max(u_write_cycle_lb[data_type_t::OUTPUT],
                                                  u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
            unsigned last_stage = u_write_cycle_lb[data_type_t::OUTPUT];

            unsigned other_stage = std::max(ratio*u_read_cycle_mac[data_type_t::OUTPUT],
                                   std::max(u_write_cycle_lb[data_type_t::OUTPUT],
                                            u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth))));

            if(num_access_lb == 1) {
                cycle_mac_lb[data_type_t::OUTPUT] += ratio*u_read_cycle_mac[data_type_t::OUTPUT] + 
                                                     u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)) +
                                                     u_write_cycle_lb[data_type_t::OUTPUT];
            } else if(num_access_lb == 2) {
                cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage + last_before_stage + last_stage;
            } else {
                cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage +
                                                     (num_access_lb - 2)*other_stage +
                                                     last_before_stage + last_stage;
            }
            
            transfer_energy[data_type_t::OUTPUT] += num_access_mac*u_transfer_energy*ceil((double)(line_size_mac[data_type_t::OUTPUT])/(double)(bitwidth));
            transfer_cycle[data_type_t::OUTPUT] += num_access_mac*u_transfer_cycle*ceil((double)(line_size_mac[data_type_t::OUTPUT])/(double)(bitwidth));

            // When not all output in local buffer is transferred to MAC.
            if(output_index < m_scheduler->offset_size_pe[data_type_t::OUTPUT].front()) {
                input_index = 0, weight_index = 0;
            }
        }
    }
}

