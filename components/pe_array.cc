#include <iostream>
#include <cmath>
#include <cstring>
#include "pe_array.h"

pe_array_t::pe_array_t(section_config_t m_section_config) :
    input_data(NULL),
    weight(NULL),
    output_data(NULL),
	workspace(NULL),
    equal_output_tile(false),
    duplicated_input(0),
    index(0),
    exist_temporal_buffer(false),
    utilization(0.0),
    write_back_cycle(0.0),
    overlapped_transfer_cycle(0.0),
    global_buffer(NULL),
    stationary_type(stationary_type_t::UNDEFINED_STATIONARY),
    array_parameter_order("kbpqcrs"),
    noc_type(noc_type_t::UNDEFINED_NOC), 
    data_format(data_format_t::CONVOLUTION),
    memory_type(memory_type_t::SEPARATE),
    height(1),
    width(1),
    num_pes(1),
    num_active_pe_x(1),
    num_active_pe_y(1),
    input_size(0),
    weight_size(0),
    output_size(0),
    frequency(0.0),
    bandwidth(0.0),
    bitwidth(0),
    initial(true),
    noc_cycle(0),
    noc_energy(0.0) {
}

pe_array_t::~pe_array_t() {

}

void pe_array_t::connect(global_buffer_t *m_global_buffer) {
    global_buffer = m_global_buffer;
    index = global_buffer->index;

    for(unsigned i = 0; i < num_pes; i++) {
        pes[i]->connect(this);
    }
}
// Get the stationary type between the PE array and Global buffer.
stationary_type_t pe_array_t::get_stationary_type() { return stationary_type; }

// Return parameter order.
std::string pe_array_t::get_parameter_order() { return array_parameter_order; }

// Get the number of PEs.
unsigned pe_array_t::get_number_of_pes() { return num_pes; }

unsigned pe_array_t::get_number_of_active_pes() { return num_active_pe_x*num_active_pe_y; }

// A signal that checks whether the PE array is idle state or not.
bool pe_array_t::is_idle() {
    bool idle;
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        idle = pes[i]->is_idle();
        if(idle == false) {break;}
    }
    return idle;
}

// Check whether the request at PE exists.
bool pe_array_t::is_exist_request_at_pe() {
    bool request_input_ = false, request_weight_ = false, request_output_ = false;

    // At least one PE sends a request to PE array
    // Input request becomes true.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::INPUT]) {
            request_input_ = true;
            break;
        }
    }
    // At least one PE sends a request to PE array. 
    // Weight request becomes true.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::WEIGHT]) {
            request_weight_ = true;
            break;
        }
    }
    // At least one PE sends a request to PE array.
    // Output request becomes true.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::OUTPUT]) {
            request_output_ = true;
            break;
        }
    }

    if(request_input_ || request_weight_ || request_output_) { return true; }
    else { return false; }
}


// Check whether the request at temporal buffer exists.
bool pe_array_t::is_exist_request_at_buffer() {
    if(request_to_global_buffer[data_type_t::INPUT] || 
       request_to_global_buffer[data_type_t::WEIGHT] || 
       request_to_global_buffer[data_type_t::OUTPUT]) { 
        return true; 
    }
    else { return false; }
}

bool pe_array_t::is_exist_data() {
    if(exist_data[data_type_t::INPUT] && 
       exist_data[data_type_t::WEIGHT] && 
       exist_data[data_type_t::OUTPUT]) { 
        return true; 
    }
    else { return false; }
}


// Wait for the data comes from Global buffer.
bool pe_array_t::wait_data() {
    if(is_waiting_data[data_type_t::INPUT] || is_waiting_data[data_type_t::WEIGHT] || is_waiting_data[data_type_t::OUTPUT]) {
        return true;
    }
    else {
        return false;
    }
}

// A signal that the data exist in the PE array.
void pe_array_t::fill_data() {
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        pes[i]->fill_data();
    }
}

void pe_array_t::request_data() {
    // The case when at least on PE in PE array request data.
    // Request data from temporal buffer in PE array to Global buffer.
    
    // Request input data to Global buffer.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::INPUT]) {

            request_to_global_buffer[data_type_t::INPUT] = true;
            is_waiting_data[data_type_t::INPUT] = true;
            global_buffer->num_request[data_type_t::INPUT]++;
            break;
        }
    }
    // Request weight to Global buffer.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::WEIGHT]) {

            request_to_global_buffer[data_type_t::WEIGHT] = true;
            is_waiting_data[data_type_t::WEIGHT] = true;
            global_buffer->num_request[data_type_t::WEIGHT]++;
            break;
        }
    }
    // Request output data to Global buffer.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::OUTPUT]) {
            if(!initial && !equal_output_tile) {
#ifdef PRINT
                std::cout << "Write back output data from temporal buffer in PE array to Global buffer" << std::endl;
#endif
#ifdef FUNCTIONAL
#endif
                is_waiting_data[data_type_t::OUTPUT] = true;
                global_buffer->num_request[data_type_t::OUTPUT]++;

                access_cycle[data_type_t::OUTPUT] += tile_size[data_type_t::OUTPUT]*u_read_cycle[data_type_t::OUTPUT]/(line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
                write_back_cycle += tile_size[data_type_t::OUTPUT]*u_read_cycle[data_type_t::OUTPUT]/(line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
                access_energy[data_type_t::OUTPUT] += tile_size[data_type_t::OUTPUT]*u_read_energy[data_type_t::OUTPUT]/(line_size[data_type_t::OUTPUT]/8/sizeof(data_t));

                global_buffer->access_cycle[data_type_t::OUTPUT] += tile_size[data_type_t::OUTPUT]*global_buffer->u_write_cycle[data_type_t::OUTPUT]/(global_buffer->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
                global_buffer->write_back_cycle += tile_size[data_type_t::OUTPUT]*global_buffer->u_write_cycle[data_type_t::OUTPUT]/(global_buffer->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
                global_buffer->access_energy[data_type_t::OUTPUT] += tile_size[data_type_t::OUTPUT]*global_buffer->u_write_energy[data_type_t::OUTPUT]/(global_buffer->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));


                transfer_cycle[data_type_t::OUTPUT] += gather_cycle;
                overlapped_transfer_cycle += gather_cycle;
                transfer_energy[data_type_t::OUTPUT] += gather_energy;

                global_buffer->transfer_cycle[data_type_t::OUTPUT] += global_buffer->u_transfer_cycle*ceil((double)tile_size[data_type_t::OUTPUT]*8*sizeof(data_t)/(double)global_buffer->get_bitwidth());
                global_buffer->overlapped_transfer_cycle += global_buffer->u_transfer_cycle*ceil((double)tile_size[data_type_t::OUTPUT]*8*sizeof(data_t)/(double)global_buffer->get_bitwidth());
                global_buffer->transfer_energy[data_type_t::OUTPUT] += global_buffer->u_transfer_energy*ceil((double)tile_size[data_type_t::OUTPUT]*8*sizeof(data_t)/(double)global_buffer->get_bitwidth());
                if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
                    if(tile_size[data_type_t::OUTPUT] == global_buffer->tile_size[data_type_t::OUTPUT]) {equal_output_tile = true;}
                }
                else {
                    if(global_buffer->multi_chip->tile_size[data_type_t::OUTPUT] == global_buffer->multi_chip->dram->tile_size[data_type_t::OUTPUT]) {
                        if(tile_size[data_type_t::OUTPUT] == global_buffer->tile_size[data_type_t::OUTPUT]) {equal_output_tile = true;}
                    }
                    else {
                        if(tile_size[data_type_t::OUTPUT] == global_buffer->tile_size[data_type_t::OUTPUT] &&
                           tile_size[data_type_t::WEIGHT] != global_buffer->tile_size[data_type_t::WEIGHT]) {equal_output_tile = true;}
                    }

                }
            }
            else {
                initial = 0;
            }
            request_to_global_buffer[data_type_t::OUTPUT] = true;
            
            break;
        }
    }
}

// Flush data at temporal buffer in PE array.
void pe_array_t::flush_data() {
    // At least one PE sends an input request to PE array.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::INPUT]) {
            exist_data[data_type_t::INPUT] = false;
            break;
        }
    }
    // At least one PE sends a weight request to PE array.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::WEIGHT]) {
            exist_data[data_type_t::WEIGHT] = false;
            break;
        }
    }
    // At least one PE sends an output request to PE array.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::OUTPUT]) {
            exist_data[data_type_t::OUTPUT] = false;
            break;
        }
    }
}

void pe_array_t::reset() {
    for(unsigned i = 0; i < num_pes; i++) {
        pes[i]->reset();
    }

    memset(input_data, 0.0, input_size);
    memset(weight, 0.0, weight_size);
    memset(output_data, 0.0, output_size);

    initial = true;
    equal_output_tile = false;

    utilization = 0.0;
    write_back_cycle = 0.0, overlapped_transfer_cycle = 0.0;

    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);

    num_request.assign(data_type_t::NUM_DATA_TYPES, 0);
    num_data_transfer.assign(data_type_t::NUM_DATA_TYPES, 0);

    access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
}

