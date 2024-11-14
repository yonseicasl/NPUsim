#include <cmath>
#include <cstring>
#include "global_buffer.h"

global_buffer_t::global_buffer_t(section_config_t m_section_config) :
    multi_chip(NULL),
    data(NULL),
    double_buffer(false),
    index(0),
    duplicated_input(0),
    u_transfer_cycle(0.0),
    u_transfer_energy(0.0),
    write_back_cycle(0.0),
    overlapped_transfer_cycle(0.0),
    transfer_output(false),
    pe_array(NULL),
    memory_type(memory_type_t::UNDEFINED_MEMORY),
    stationary_type(stationary_type_t::UNDEFINED_STATIONARY),
    parameter_order("KBPQCRS"),
    size(0),
    frequency(0.0),
    bandwidth(0.0),
    bitwidth(0.0),
    input_index(0),
    weight_index(0),
    output_index(0),
    input_flush_counter(0),
    weight_flush_counter(0),
    output_flush_counter(0),
    idle(false),
    initial(true) {
}

global_buffer_t::~global_buffer_t() {

}

void global_buffer_t::connect(pe_array_t *m_pe_array) {
    pe_array = m_pe_array;
}

void global_buffer_t::connect(multi_chip_t *m_multi_chip) {
    multi_chip = m_multi_chip;
}

void global_buffer_t::update_tile_size(scheduler_t *m_scheduler) {
    tile_size = m_scheduler->tile_size[component_type_t::GLOBAL_BUFFER];
    
    update_offset();
    check_tile_size();
}

// Get global buffer's stationary type
stationary_type_t global_buffer_t::get_stationary_type() { return stationary_type; }

// Get global buffer's memory type
memory_type_t global_buffer_t::get_memory_type() { return memory_type;}

// Get global buffer's size
double global_buffer_t::get_buffer_size() { return size; }

// Get global buffer's bitwidth
unsigned global_buffer_t::get_bitwidth() { return bitwidth; }

// A signal to check whether the Global buffer is idle state or not.
bool global_buffer_t::is_idle() { return idle; }

// A signal to check whether all data exist in the Global buffer or not.
bool global_buffer_t::is_exist_data() {
    if(exist_data[data_type_t::INPUT] && exist_data[data_type_t::WEIGHT] && exist_data[data_type_t::OUTPUT]) {
        return true;
    }
    else { return false; }
}

// A signal to check whether at least one request exists from Global buffer to DRAM or not.
bool global_buffer_t::is_exist_request() {
    if(request_to_multi_chip[data_type_t::INPUT] || request_to_multi_chip[data_type_t::WEIGHT] || request_to_multi_chip[data_type_t::OUTPUT]) {
        return true;
    }
    else { return false; }
}

// Wait for the data comes from DRAM.
void global_buffer_t::wait_data() { idle = true; }

// A signal indicates that the data exists in the Global buffer.
void global_buffer_t::fill_data() { idle = false; }

// Send a request signal to DRAM.
void global_buffer_t::request_data() {
    // When the data does not exist in Global buffer
    // Send a request signal to DRAM.
    // Case 1. Input data
    if(!exist_data[data_type_t::INPUT]) {
        // Send a input data request signal to DRAM.
        request_to_multi_chip[data_type_t::INPUT] = true;
        multi_chip->exist_data[data_type_t::INPUT] = false;

        multi_chip->num_request[data_type_t::INPUT]++;
    }
    // Case 2. Weight data
    if(!exist_data[data_type_t::WEIGHT]) {
        // Send a weight request signal to DRAM 
        request_to_multi_chip[data_type_t::WEIGHT] = true;
        multi_chip->exist_data[data_type_t::WEIGHT] = false;

        multi_chip->num_request[data_type_t::WEIGHT]++;
    }
    // Case 3. Output data
    if(!exist_data[data_type_t::OUTPUT]) {
        if(initial) {
            initial = false;
        }
        request_to_multi_chip[data_type_t::OUTPUT] = true;
        multi_chip->exist_data[data_type_t::OUTPUT] = false;

        multi_chip->num_request[data_type_t::OUTPUT]++;
    }
}

// Transfer the data to temporal buffer of PE array.
void global_buffer_t::data_transfer(scheduler_t *m_scheduler) {

    dynamic_power = u_static_power;

    if(!bypass[data_type_t::INPUT]) {
        utilization[data_type_t::INPUT] = (float)(tile_size[data_type_t::INPUT]*sizeof(data_t))/(float)(size);
    }
    if(bypass[data_type_t::WEIGHT]) {
        utilization[data_type_t::WEIGHT] = (float)(tile_size[data_type_t::WEIGHT]*sizeof(data_t))/(float)(size);
    }
    if(bypass[data_type_t::OUTPUT]) {
        utilization[data_type_t::OUTPUT] = (float)(tile_size[data_type_t::OUTPUT]*sizeof(data_t))/(float)(size);
    }



    // Transfer input data from Global buffer to temporal buffer of PE array.
    if(pe_array->request_to_global_buffer[data_type_t::INPUT]) {
#ifdef FUNCTIONAL
        // Input data transfer
        m_scheduler->transfer_data(pe_array->input_data, data, 0, offsets[data_type_t::INPUT] + m_scheduler->input_offset_global_buffer.front(), 
                                   component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, 
                                   data_type_t::INPUT, pe_array->get_stationary_type(), action_type_t::LOAD);

        // Update for NPUsim ver2
        //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
        //m_scheduler->transfer_data_ver2(pe_array->input_data, data, 
        //                                component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, 
        //                                data_type_t::INPUT, pe_array->get_stationary_type(), action_type_t::LOAD, true);
                                   
        // Case 1. Dense data format
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::INPUT]) {
                dynamic_power[data_type_t::INPUT] = u_dynamic_power[data_type_t::INPUT];
                // Update PE array write cycle and energy.
                num_data_transfer[data_type_t::INPUT]++;

                std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);
                parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);

                uint64_t address_pe_array = 0, address_global_buffer = 0;
                unsigned num_access_array = 0, num_access_global_buffer = 0;
                for(unsigned b = 0; b < parameters_pe_array[parameter_type_t::BATCH_SIZE]; b++) {
                    for(unsigned c = 0; c < parameters_pe_array[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned h = 0; h < parameters_pe_array[parameter_type_t::INPUT_HEIGHT]; h++) {
                            for(unsigned w = 0; w < parameters_pe_array[parameter_type_t::INPUT_WIDTH]; w++) {
                                // Check temporal buffer address of PE array.
                                if(address_pe_array != ((uint64_t)&pe_array->input_data[b*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                         *parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                         *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                        c*parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                         *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                        h*parameters_pe_array[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                        pe_array->mask_bits[data_type_t::INPUT]) << pe_array->mask_bits[data_type_t::INPUT]) {

                                    // Update temporal buffer's costs of PE array
                                    pe_array->access_energy[data_type_t::INPUT] += pe_array->u_write_energy[data_type_t::INPUT];
                                    pe_array->access_cycle[data_type_t::INPUT] += pe_array->u_write_cycle[data_type_t::INPUT];
                                    num_access_array++;

                                    // Update temporal buffer address
                                    address_pe_array = ((uint64_t)&pe_array->input_data[b*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                         *parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                         *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                        c*parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                         *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                        h*parameters_pe_array[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                        pe_array->mask_bits[data_type_t::INPUT]) << pe_array->mask_bits[data_type_t::INPUT];
                                }
                                // Check global buffer address
                                if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::INPUT] + 
                                                                             m_scheduler->input_offset_global_buffer.front() +
                                                                             b*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] +
                                                                             c*parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                             h*parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                             mask_bits[data_type_t::INPUT] ) << mask_bits[data_type_t::INPUT]) {

                                    // Update global buffer cost
                                    access_energy[data_type_t::INPUT] += u_read_energy[data_type_t::INPUT];
                                    access_cycle[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT];
                                    num_access_global_buffer++;

                                    // Update global buffer address
                                    address_global_buffer = ((uint64_t)&data[offsets[data_type_t::INPUT] + 
                                                                             m_scheduler->input_offset_global_buffer.front() +
                                                                             b*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] +
                                                                             c*parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                             h*parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                             mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT];
                                }
                            }
                        }
                    }
                }

                // Update overlapped cost between PE array and global buffer
                unsigned ratio = ceil((double)(line_size[data_type_t::INPUT])/(double)(pe_array->line_size[data_type_t::INPUT]));

                // At the 1, 2, before last, and last stages
                unsigned first_stage = u_read_cycle[data_type_t::INPUT];
                unsigned second_stage = std::max(u_read_cycle[data_type_t::INPUT], 
                                                 u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)));

                unsigned last_before_stage = std::max(ratio*pe_array->u_write_cycle[data_type_t::INPUT],
                                                      u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)));
                unsigned last_stage = ratio*pe_array->u_write_cycle[data_type_t::INPUT];

                // Remainder stage
                unsigned other_stage = std::max(u_read_cycle[data_type_t::INPUT], 
                                       std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)),
                                                ratio*pe_array->u_write_cycle[data_type_t::INPUT]));

                if(num_access_global_buffer == 1) {
                    cycle_pe_array_global_buffer[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT] + 
                                                                        u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)) +
                                                                        pe_array->u_write_cycle[data_type_t::INPUT];
                } else {
                    cycle_pe_array_global_buffer[data_type_t::INPUT] += first_stage + second_stage +
                                                                        (num_access_global_buffer - 2)*other_stage +
                                                                        last_before_stage + last_stage;
                }

                // Update data transfer cycle and energy between Global buffer and PE array.
                transfer_cycle[data_type_t::INPUT] += num_access_global_buffer*u_transfer_cycle
                                                      *ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth));
                transfer_energy[data_type_t::INPUT] += num_access_global_buffer*u_transfer_energy
                                                      *ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth));

                pe_array->skip_transfer[data_type_t::INPUT] = false;

                //if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y -1) {
                //    move_front(&m_scheduler->input_offset_dram);
                //}
            }
        }
        // Case 2. COO data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
        }
        // Case 3. CSC data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::INPUT]) {
                // Column index calculation
                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::PE_Y);
                unsigned row = parameters[parameter_type_t::INPUT_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer[data_type_t::INPUT]++;

                // Update global buffer access cost
                access_cycle[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                    (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT]
                                                   /(sizeof(data_t)*8/row_bit) + // Row index
                                                    parameters[parameter_type_t::BATCH_SIZE]
                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                   *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                   *u_read_cycle[data_type_t::INPUT]
                                                   /(sizeof(data_t)*8/row_bit); // Column pointer
                access_energy[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT]  + // Non-zero data
                                                     (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT]/(sizeof(data_t)*8/row_bit) + // Row index
                                                     parameters[parameter_type_t::BATCH_SIZE]
                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                    *u_read_energy[data_type_t::INPUT]
                                                    /(sizeof(data_t)*8/row_bit); // Column pointer

                // Update PE array access cost (if temporal buffer exist)
                pe_array->access_cycle[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *pe_array->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                              (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *pe_array->u_write_cycle[data_type_t::INPUT]
                                                             /(sizeof(data_t)*8/row_bit) + // Row index
                                                              parameters[parameter_type_t::BATCH_SIZE]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                             *pe_array->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/row_bit); // Column pointer

                pe_array->access_energy[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                              *pe_array->u_write_energy[data_type_t::INPUT] + // Non-zero data
                                                               (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                              *pe_array->u_write_energy[data_type_t::INPUT]
                                                              /(sizeof(data_t)*8/row_bit) + // Row index
                                                               parameters[parameter_type_t::BATCH_SIZE]
                                                              *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                              *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                              *pe_array->u_write_energy[data_type_t::INPUT]
                                                              /(sizeof(data_t)*8/row_bit); // Column pointer

                // Update overlapped cycle between the Global buffer and temporal buffer in PE array
                cycle_pe_array_global_buffer[data_type_t::INPUT] += std::max((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *u_read_cycle[data_type_t::INPUT]
                                                                            /(sizeof(data_t)*8/row_bit) + // Row index
                                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                            *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                                            *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/row_bit), // Column pointer
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT]
                                                                            /(sizeof(data_t)*8/row_bit) + // Row index
                                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                            *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT]
                                                                            /(sizeof(data_t)*8/row_bit)); // Column pointer

                // Update transfer cost between the global buffer and temporal buffer in PE array
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                      u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)/(float)bitwidth) + // Row index
                                                      u_transfer_cycle*ceil((float)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                   *(parameters[parameter_type_t::INPUT_WIDTH+1])*row_bit)/(float)bitwidth); // Column pointer
                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)/(float)bitwidth) + // Row index
                                                       u_transfer_energy*ceil((float)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                     *(parameters[parameter_type_t::INPUT_WIDTH+1])*row_bit)/(float)bitwidth); // Column pointer
                pe_array->skip_transfer[data_type_t::INPUT] = false;

            }
        }
        // Case 4. CSR data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::INPUT]) {
                // Column index calculation
                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::PE_Y);
                unsigned column = parameters[parameter_type_t::INPUT_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer[data_type_t::INPUT]++;

                // Update global buffer access cost
                access_cycle[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                    (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT]
                                                   /(sizeof(data_t)*8/column_bit) + // Column index
                                                    parameters[parameter_type_t::BATCH_SIZE]
                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                   *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                   *u_read_cycle[data_type_t::INPUT]
                                                   /(sizeof(data_t)*8/column_bit); // Row pointer
                access_energy[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT]  + // Non-zero data
                                                     (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT]/(sizeof(data_t)*8/column_bit) + // Column index
                                                     parameters[parameter_type_t::BATCH_SIZE]
                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                    *u_read_energy[data_type_t::INPUT]
                                                    /(sizeof(data_t)*8/column_bit); // Row pointer

                // Update PE array access cost (if temporal buffer exist)
                pe_array->access_cycle[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *pe_array->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                              (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *pe_array->u_write_cycle[data_type_t::INPUT]
                                                             /(sizeof(data_t)*8/column_bit) + // Column index
                                                              parameters[parameter_type_t::BATCH_SIZE]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                             *pe_array->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/column_bit); // Row pointer

                pe_array->access_energy[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                              *pe_array->u_write_energy[data_type_t::INPUT] + // Non-zero data
                                                               (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                              *pe_array->u_write_energy[data_type_t::INPUT]
                                                              /(sizeof(data_t)*8/column_bit) + // Column index
                                                               parameters[parameter_type_t::BATCH_SIZE]
                                                              *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                              *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                              *pe_array->u_write_energy[data_type_t::INPUT]
                                                              /(sizeof(data_t)*8/column_bit); // Row pointer

                // Update overlapped cycle between the Global buffer and temporal buffer in PE array
                cycle_pe_array_global_buffer[data_type_t::INPUT] += std::max((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *u_read_cycle[data_type_t::INPUT]
                                                                            /(sizeof(data_t)*8/column_bit) + // Column index
                                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                            *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                                            *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/column_bit), // Row pointer
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT]
                                                                            /(sizeof(data_t)*8/column_bit) + // Column index
                                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                            *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT]
                                                                            /(sizeof(data_t)*8/column_bit)); // Row pointer

                // Update transfer cost between the global buffer and temporal buffer in PE array
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                      u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)/(float)bitwidth) + // Column index
                                                      u_transfer_cycle*ceil((float)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                   *(parameters[parameter_type_t::INPUT_HEIGHT+1])*column_bit)/(float)bitwidth); // Row pointer
                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)/(float)bitwidth) + // Column index
                                                       u_transfer_energy*ceil((float)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                     *(parameters[parameter_type_t::INPUT_HEIGHT+1])*column_bit)/(float)bitwidth); // Row pointer
                pe_array->skip_transfer[data_type_t::INPUT] = false;

                //if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y -1) {
                //    move_front(&m_scheduler->input_offset_dram);
                //}
            }
        }
        // Case 4. SparseMap
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::INPUT]) {

                num_data_transfer[data_type_t::INPUT]++;
                // Update Global buffer read cycle and energy.

                // Update global buffer cost
                access_cycle[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                    pe_array->tile_size[data_type_t::INPUT]
                                                   *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
                access_energy[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT] + // Non-zero data
                                                     pe_array->tile_size[data_type_t::INPUT]
                                                     *u_read_energy[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata

                // Update PE array cost
                pe_array->access_cycle[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *pe_array->u_write_cycle[data_type_t::INPUT] + // Non-zero value
                                                              pe_array->tile_size[data_type_t::INPUT]
                                                             *pe_array->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
                pe_array->access_energy[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                              *pe_array->u_write_energy[data_type_t::INPUT] + // Non-zero data
                                                               pe_array->tile_size[data_type_t::INPUT]
                                                              *pe_array->u_write_energy[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata

                // Update overlapped cycle between the global buffer and temporal buffer in PE array
                cycle_pe_array_global_buffer[data_type_t::INPUT] += std::max((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                                             pe_array->tile_size[data_type_t::INPUT]
                                                                            *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8), // Metadata
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT] + // Non-zero data 
                                                                             pe_array->tile_size[data_type_t::INPUT]
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8)); // Metadata
               
                // Update transfer cost between the global buffer and temporal buffer in PE array
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                      u_transfer_cycle*ceil((float)(pe_array->tile_size[data_type_t::INPUT])/(float)bitwidth); // Meta data

                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                     + u_transfer_energy*ceil((float)(pe_array->tile_size[data_type_t::INPUT])/(float)bitwidth); // Meta data

                pe_array->skip_transfer[data_type_t::INPUT] = false;

                //if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                //    move_front(&m_scheduler->input_offset_global_buffer);
                //}
            }
        }
        if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y -1) {
            move_front(&m_scheduler->input_offset_dram);
        }
#else
        /* Stats */
        if(!skip_transfer[data_type_t::INPUT]) {
            // Update PE array write cycle and energy.
            num_data_transfer[data_type_t::INPUT]++;

            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);
            parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);

            uint64_t address_pe_array = 0, address_global_buffer = 0;
            unsigned num_access_array = 0, num_access_global_buffer = 0;
            for(unsigned b = 0; b < parameters_pe_array[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned c = 0; c < parameters_pe_array[parameter_type_t::INPUT_CHANNEL]; c++) {
                    for(unsigned h = 0; h < parameters_pe_array[parameter_type_t::INPUT_HEIGHT]; h++) {
                        for(unsigned w = 0; w < parameters_pe_array[parameter_type_t::INPUT_WIDTH]; w++) {
                            // Check temporal buffer address of PE array.
                            if(address_pe_array != ((uint64_t)&pe_array->input_data[b*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                     *parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                    c*parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                    h*parameters_pe_array[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                    pe_array->mask_bits[data_type_t::INPUT]) << pe_array->mask_bits[data_type_t::INPUT]) {

                                // Update temporal buffer's costs of PE array
                                pe_array->access_energy[data_type_t::INPUT] += pe_array->u_write_energy[data_type_t::INPUT];
                                pe_array->access_cycle[data_type_t::INPUT] += pe_array->u_write_cycle[data_type_t::INPUT];
                                num_access_array++;

                                // Update temporal buffer address
                                address_pe_array = ((uint64_t)&pe_array->input_data[b*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                     *parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                    c*parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                    h*parameters_pe_array[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                    pe_array->mask_bits[data_type_t::INPUT]) << pe_array->mask_bits[data_type_t::INPUT];
                            }
                            // Check global buffer address
                            if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::INPUT] + 
                                                                         m_scheduler->input_offset_global_buffer.front() +
                                                                         b*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] +
                                                                         c*parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                         h*parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                         mask_bits[data_type_t::INPUT] ) << mask_bits[data_type_t::INPUT]) {

                                // Update global buffer cost
                                access_energy[data_type_t::INPUT] += u_read_energy[data_type_t::INPUT];
                                access_cycle[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT];
                                num_access_global_buffer++;

                                // Update global buffer address
                                address_global_buffer = ((uint64_t)&data[offsets[data_type_t::INPUT] + 
                                                                         m_scheduler->input_offset_global_buffer.front() +
                                                                         b*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] +
                                                                         c*parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                         h*parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                         mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT];
                            }
                        }
                    }
                }
            }

            // Update overlapped cost between PE array and global buffer
            unsigned ratio = ceil((double)(line_size[data_type_t::INPUT])/(double)(pe_array->line_size[data_type_t::INPUT]));

            // At the 1, 2, before last, and last stages
            unsigned first_stage = u_read_cycle[data_type_t::INPUT];
            unsigned second_stage = std::max(u_read_cycle[data_type_t::INPUT], 
                                             u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)));

            unsigned last_before_stage = std::max(ratio*pe_array->u_write_cycle[data_type_t::INPUT],
                                                  u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)));
            unsigned last_stage = ratio*pe_array->u_write_cycle[data_type_t::INPUT];

            // Remainder stage
            unsigned other_stage = std::max(u_read_cycle[data_type_t::INPUT], 
                                   std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)),
                                            ratio*pe_array->u_write_cycle[data_type_t::INPUT]));

            if(num_access_global_buffer == 1) {
                cycle_pe_array_global_buffer[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT] + 
                                                                    u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)) +
                                                                    pe_array->u_write_cycle[data_type_t::INPUT];
            } else {
                cycle_pe_array_global_buffer[data_type_t::INPUT] += first_stage + second_stage +
                                                                    (num_access_global_buffer - 2)*other_stage +
                                                                    last_before_stage + last_stage;
            }

            // Update data transfer cycle and energy between Global buffer and PE array.
            transfer_cycle[data_type_t::INPUT] += num_access_global_buffer*u_transfer_cycle
                                                  *ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth));
            transfer_energy[data_type_t::INPUT] += num_access_global_buffer*u_transfer_energy
                                                  *ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth));

            pe_array->skip_transfer[data_type_t::INPUT] = false;
            if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                move_front(&m_scheduler->input_offset_global_buffer);
            }
        }
#endif

        // Increment the input index that finds the offset of input data in Global buffer.
        input_index++;
        // Input data exists in PE array.
        // Request input data to Global buffer is not required.
        pe_array->exist_data[data_type_t::INPUT] = true, pe_array->request_to_global_buffer[data_type_t::INPUT] = false;
        if(pe_array->tile_size[data_type_t::INPUT] == tile_size[data_type_t::INPUT]) { skip_transfer[data_type_t::INPUT] = true; }
    }
    // Transfer weight from Global buffer to temporal buffer of PE array.
    if(pe_array->request_to_global_buffer[data_type_t::WEIGHT]) {
#ifdef FUNCTIONAL
        // Weight transfer
        m_scheduler->transfer_data(pe_array->weight, data, 0, offsets[data_type_t::WEIGHT] + m_scheduler->weight_offset_global_buffer.front(),
                                   component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, 
                                   data_type_t::WEIGHT, pe_array->get_stationary_type(), action_type_t::LOAD);
        // Update for NPUsim ver2
        //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
        //m_scheduler->transfer_data_ver2(pe_array->weight, data, 
        //                                component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, 
        //                                data_type_t::WEIGHT, pe_array->get_stationary_type(), action_type_t::LOAD, true);
        
        // Case 1. Dense
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                dynamic_power[data_type_t::WEIGHT] = u_dynamic_power[data_type_t::WEIGHT];
                num_data_transfer[data_type_t::WEIGHT]++;

                std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);
                parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);

                uint64_t address_pe_array = 0, address_global_buffer = 0;
                unsigned num_access_pe_array = 0, num_access_global_buffer = 0;

                for(unsigned k = 0; k < parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned c = 0; c < parameters_pe_array[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned r = 0; r < parameters_pe_array[parameter_type_t::FILTER_HEIGHT]; r++) {
                            for(unsigned s = 0; s < parameters_pe_array[parameter_type_t::FILTER_WIDTH]; s++) {
                                // Check temporal buffer address at PE array
                                if(address_pe_array != ((uint64_t)&pe_array->weight[k*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                     *parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                                    c*parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::FILTER_WIDTH] +
                                                                                    r*parameters_pe_array[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                    pe_array->mask_bits[data_type_t::WEIGHT]) << pe_array->mask_bits[data_type_t::WEIGHT]) {

                                    // Update PE array cost
                                    pe_array->access_energy[data_type_t::WEIGHT] += pe_array->u_write_energy[data_type_t::WEIGHT];
                                    pe_array->access_cycle[data_type_t::WEIGHT] += pe_array->u_write_cycle[data_type_t::WEIGHT];
                                    num_access_pe_array++;

                                    // Update PE array address (temporal buffer)
                                    address_pe_array = ((uint64_t)&pe_array->weight[k*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                     *parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                                    c*parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::FILTER_WIDTH] +
                                                                                    r*parameters_pe_array[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                    pe_array->mask_bits[data_type_t::WEIGHT]) << pe_array->mask_bits[data_type_t::WEIGHT];
                                }

                                // Check global buffer address
                                if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::WEIGHT] + 
                                                                             m_scheduler->weight_offset_global_buffer.front() +
                                                                             k*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] +
                                                                             c*parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                             r*parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                             mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT]) {

                                    // Update global buffer cost
                                    access_energy[data_type_t::WEIGHT] += u_read_energy[data_type_t::WEIGHT];
                                    access_cycle[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT];
                                    num_access_global_buffer++;

                                    // Update Global buffer address
                                    address_global_buffer = ((uint64_t)&data[offsets[data_type_t::WEIGHT] + 
                                                                             m_scheduler->weight_offset_global_buffer.front() +
                                                                             k*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] +
                                                                             c*parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                             r*parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                             mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT];
                                }
                            }
                        }
                    }
                }

                // Update overlapped cost between PE array and global buffer
                unsigned ratio = ceil((double)(line_size[data_type_t::WEIGHT])/(double)(pe_array->line_size[data_type_t::WEIGHT]));

                // At the 1, 2, before last, and last stages
                unsigned first_stage = u_read_cycle[data_type_t::WEIGHT];
                unsigned second_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                                 u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
                unsigned last_before_stage = std::max(ratio*pe_array->u_write_cycle[data_type_t::WEIGHT],
                                                      u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
                unsigned last_stage = ratio*pe_array->u_write_cycle[data_type_t::WEIGHT];

                // Remainder stages
                unsigned other_stage = std::max(u_read_cycle[data_type_t::WEIGHT], 
                                       std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)),
                                                ratio*pe_array->u_write_cycle[data_type_t::WEIGHT]));

                if(num_access_global_buffer == 1) {
                    cycle_pe_array_global_buffer[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT] +   
                                                                         u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)) +
                                                                         ratio*pe_array->u_write_cycle[data_type_t::WEIGHT];
                } else {
                    cycle_pe_array_global_buffer[data_type_t::WEIGHT] += first_stage + second_stage +
                                                                         (num_access_global_buffer - 2)*other_stage +
                                                                         last_before_stage + last_stage;
                }

                // Update data transfer cycle and energy between Global buffer and PE array.
                transfer_cycle[data_type_t::WEIGHT] += num_access_global_buffer*u_transfer_cycle
                                                       *ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth));
                transfer_energy[data_type_t::WEIGHT] += num_access_global_buffer*u_transfer_energy
                                                        *ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth));
                pe_array->skip_transfer[data_type_t::WEIGHT] = false;

                if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                    move_front(&m_scheduler->weight_offset_global_buffer);
                }
            }
        }
        // Case 2. COO data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
        }
        // Case 3. CSC format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                // Row bit calculation
                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::PE_Y);
                unsigned row = parameters[parameter_type_t::FILTER_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer[data_type_t::WEIGHT]++;

                // Update global buffer cost
                access_cycle[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                     (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT]/(sizeof(data_t)*8/row_bit) + // Row index
                                                     parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                    *u_read_cycle[data_type_t::WEIGHT]
                                                    /(sizeof(data_t)*8/row_bit); // Column pointer
                access_energy[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                      (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT]
                                                     /(sizeof(data_t)*8/row_bit) + // Row index
                                                      parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                     *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                     *u_read_energy[data_type_t::WEIGHT]
                                                     /(sizeof(data_t)*8/row_bit); // Column pointer
                                                     
                // Update PE array cost (if temporal buffer exist)
                pe_array->access_cycle[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                               (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                              /(sizeof(data_t)*8/row_bit) + // Row index
                                                               parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                              *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                              *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                              /(sizeof(data_t)*8/row_bit); // Column pointer
                pe_array->access_energy[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT] + // Non-zero data
                                                                (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT]
                                                               /(sizeof(data_t)*8/row_bit) + // Row index
                                                                parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT]
                                                               /(sizeof(data_t)*8/row_bit); // Column pointer

                // Update overlapped cycle between the global buffer and PE array
                cycle_pe_array_global_buffer[data_type_t::WEIGHT] += std::max((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                              (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *u_read_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/row_bit) + // Row index
                                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                                             *u_read_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/row_bit),  // Column pointer
                                                                              (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *pe_array->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                              (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/row_bit) + // Row index
                                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                                             *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/row_bit)); // Column pointer

                // Update transfer cost between the global buffer and PE array
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                       u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*row_bit)/(float)bitwidth) + // Row index
                                                       u_transfer_cycle*ceil((float)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *(parameters[parameter_type_t::FILTER_WIDTH+1])*row_bit)/(float)bitwidth); // Column pointer
                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                        u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*row_bit)/(float)bitwidth) + // Row index
                                                        u_transfer_energy*ceil((float)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                      *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                      *(parameters[parameter_type_t::FILTER_WIDTH+1])*row)/(float)bitwidth); // Column pointer
                pe_array->skip_transfer[data_type_t::WEIGHT] = false;
                if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                    move_front(&m_scheduler->weight_offset_global_buffer);
                }
            }
        }
        // Case 4. CSR format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                // Column bit calculation
                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::PE_Y);
                unsigned column = parameters[parameter_type_t::FILTER_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer[data_type_t::WEIGHT]++;

                // Update global buffer cost
                access_cycle[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                     (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT]/(sizeof(data_t)*8/column_bit) + // Column index
                                                     parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                    *u_read_cycle[data_type_t::WEIGHT]
                                                    /(sizeof(data_t)*8/column_bit); // Row pointer
                access_energy[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                      (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT]
                                                     /(sizeof(data_t)*8/column_bit) + // Column index
                                                      parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                     *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                     *u_read_energy[data_type_t::WEIGHT]
                                                     /(sizeof(data_t)*8/column_bit); // Row pointer
                                                     
                // Update PE array cost (if temporal buffer exist)
                pe_array->access_cycle[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                               (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                              /(sizeof(data_t)*8/column_bit) + // Column index
                                                               parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                              *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                              *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                              /(sizeof(data_t)*8/column_bit); // Row pointer
                pe_array->access_energy[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT] + // Non-zero data
                                                                (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT]
                                                               /(sizeof(data_t)*8/column_bit) + // Column index
                                                                parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT]
                                                               /(sizeof(data_t)*8/column_bit); // Row pointer

                // Update overlapped cycle between the global buffer and PE array
                cycle_pe_array_global_buffer[data_type_t::WEIGHT] += std::max((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                              (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *u_read_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/column_bit) + // Column index
                                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                                             *u_read_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/column_bit),  // Row pointer
                                                                              (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *pe_array->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                              (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/column_bit) + // Column index
                                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                                             *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/column_bit)); // row pointer

                // Update transfer cost between the global buffer and PE array
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                       u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*column_bit)/(float)bitwidth) + // Column index
                                                       u_transfer_cycle*ceil((float)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *(parameters[parameter_type_t::FILTER_HEIGHT+1])*column_bit)/(float)bitwidth); // Row pointer
                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                        u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*column_bit)/(float)bitwidth) + // Column index
                                                        u_transfer_energy*ceil((float)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                      *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                      *(parameters[parameter_type_t::FILTER_HEIGHT+1])*column_bit)/(float)bitwidth); // Row pointer
                pe_array->skip_transfer[data_type_t::WEIGHT] = false;
                if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                    move_front(&m_scheduler->weight_offset_global_buffer);
                }
            }
        }
        // Case 5. SparseMap
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                num_data_transfer[data_type_t::WEIGHT]++;

                // Update Global buffer access cost
                access_cycle[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                     pe_array->tile_size[data_type_t::WEIGHT]
                                                    *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
                access_energy[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                      pe_array->tile_size[data_type_t::WEIGHT]
                                                     *u_read_energy[data_type_t::WEIGHT]/(sizeof(data_t)*8); // Metadata

                // Update PE array access cost (if temporal buffer exist)
                pe_array->access_cycle[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                               pe_array->tile_size[data_type_t::WEIGHT]
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT]/(sizeof(data_t)*8); // Metadata
                pe_array->access_energy[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT] + // Non-zero data
                                                                pe_array->tile_size[data_type_t::WEIGHT]
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT]/(sizeof(data_t)*8); // Metadata

                // Update overlapped cycle between the global buffer and PE array
                cycle_pe_array_global_buffer[data_type_t::WEIGHT] += std::max((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                              *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                               pe_array->tile_size[data_type_t::WEIGHT]
                                                                              *u_read_cycle[data_type_t::WEIGHT]/(sizeof(data_t)*8), // Metadata
                                                                               (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                               pe_array->tile_size[data_type_t::WEIGHT]
                                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT]/(sizeof(data_t)*8)); // Metadata


                // Update transfer cost between the global buffer and PE array
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                       u_transfer_cycle*ceil((float)(pe_array->tile_size[data_type_t::WEIGHT])/(float)bitwidth); // Metadata

                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                        u_transfer_energy*ceil((float)(pe_array->tile_size[data_type_t::WEIGHT])/(float)bitwidth); // Meta data

                pe_array->skip_transfer[data_type_t::WEIGHT] = false;
                if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                    move_front(&m_scheduler->weight_offset_global_buffer);
                }
            }
        }
        else {
            std::cerr << "Undefined compressed data type" << std::endl;
            exit(1);
        }

        /* Stats */

#else
        if(!skip_transfer[data_type_t::WEIGHT]) {
            num_data_transfer[data_type_t::WEIGHT]++;

            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);
            parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);

            uint64_t address_pe_array = 0, address_global_buffer = 0;
            unsigned num_access_pe_array = 0, num_access_global_buffer = 0;

            for(unsigned k = 0; k < parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                for(unsigned c = 0; c < parameters_pe_array[parameter_type_t::INPUT_CHANNEL]; c++) {
                    for(unsigned r = 0; r < parameters_pe_array[parameter_type_t::FILTER_HEIGHT]; r++) {
                        for(unsigned s = 0; s < parameters_pe_array[parameter_type_t::FILTER_WIDTH]; s++) {
                            // Check temporal buffer address at PE array
                            if(address_pe_array != ((uint64_t)&pe_array->weight[k*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                 *parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                 *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                                c*parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                 *parameters_pe_array[parameter_type_t::FILTER_WIDTH] +
                                                                                r*parameters_pe_array[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                pe_array->mask_bits[data_type_t::WEIGHT]) << pe_array->mask_bits[data_type_t::WEIGHT]) {

                                // Update PE array cost
                                pe_array->access_energy[data_type_t::WEIGHT] += pe_array->u_write_energy[data_type_t::WEIGHT];
                                pe_array->access_cycle[data_type_t::WEIGHT] += pe_array->u_write_cycle[data_type_t::WEIGHT];
                                num_access_pe_array++;

                                // Update PE array address (temporal buffer)
                                address_pe_array = ((uint64_t)&pe_array->weight[k*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                 *parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                 *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                                c*parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                 *parameters_pe_array[parameter_type_t::FILTER_WIDTH] +
                                                                                r*parameters_pe_array[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                pe_array->mask_bits[data_type_t::WEIGHT]) << pe_array->mask_bits[data_type_t::WEIGHT];
                            }

                            // Check global buffer address
                            if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::WEIGHT] + 
                                                                         m_scheduler->weight_offset_global_buffer.front() +
                                                                         k*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] +
                                                                         c*parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                         r*parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                         mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT]) {

                                // Update global buffer cost
                                access_energy[data_type_t::WEIGHT] += u_read_energy[data_type_t::WEIGHT];
                                access_cycle[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT];
                                num_access_global_buffer++;

                                // Update Global buffer address
                                address_global_buffer = ((uint64_t)&data[offsets[data_type_t::WEIGHT] + 
                                                                         m_scheduler->weight_offset_global_buffer.front() + 
                                                                         k*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] +
                                                                         c*parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                         r*parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                         mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT];
                            }
                        }
                    }
                }
            }

            // Update overlapped cost between PE array and global buffer
            unsigned ratio = ceil((double)(line_size[data_type_t::WEIGHT])/(double)(pe_array->line_size[data_type_t::WEIGHT]));

            // At the 1, 2, before last, and last stages
            unsigned first_stage = u_read_cycle[data_type_t::WEIGHT];
            unsigned second_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                             u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
            unsigned last_before_stage = std::max(ratio*pe_array->u_write_cycle[data_type_t::WEIGHT],
                                                  u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
            unsigned last_stage = ratio*pe_array->u_write_cycle[data_type_t::WEIGHT];

            // Remainder stages
            unsigned other_stage = std::max(u_read_cycle[data_type_t::WEIGHT], 
                                   std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)),
                                            ratio*pe_array->u_write_cycle[data_type_t::WEIGHT]));

            if(num_access_global_buffer == 1) {
                cycle_pe_array_global_buffer[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT] +   
                                                                     u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)) +
                                                                     ratio*pe_array->u_write_cycle[data_type_t::WEIGHT];
            } else {
                cycle_pe_array_global_buffer[data_type_t::WEIGHT] += first_stage + second_stage +
                                                                     (num_access_global_buffer - 2)*other_stage +
                                                                     last_before_stage + last_stage;
            }

            // Update data transfer cycle and energy between Global buffer and PE array.
            transfer_cycle[data_type_t::WEIGHT] += num_access_global_buffer*u_transfer_cycle
                                                   *ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth));
            transfer_energy[data_type_t::WEIGHT] += num_access_global_buffer*u_transfer_energy
                                                    *ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth));

            pe_array->skip_transfer[data_type_t::WEIGHT] = false;

            if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                move_front(&m_scheduler->weight_offset_global_buffer);
            }
        }
#endif
        /* Stats */

        // Increase the weight index by one.
        weight_index++;

        // Update PE array signals
        pe_array->exist_data[data_type_t::WEIGHT] = true, pe_array->request_to_global_buffer[data_type_t::WEIGHT] = false;
        if(pe_array->tile_size[data_type_t::WEIGHT] == tile_size[data_type_t::WEIGHT]) { skip_transfer[data_type_t::WEIGHT] = true; }
    }
    // Transfer output data from Global buffer to temporal buffer of PE array.
    if(pe_array->request_to_global_buffer[data_type_t::OUTPUT]) {
        if(m_scheduler->output_read_global_buffer[m_scheduler->output_offset_global_buffer.front()]) {
#ifdef FUNCTIONAL
            // Output data transfer
            m_scheduler->transfer_data(pe_array->output_data, data, 0, offsets[data_type_t::OUTPUT] + m_scheduler->output_offset_global_buffer.front(), 
                                       component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, 
                                       data_type_t::OUTPUT, pe_array->get_stationary_type(), action_type_t::LOAD);

            // Update for NPUsim ver2
            //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
            //m_scheduler->transfer_data_ver2(pe_array->output_data, data, 
            //                                component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, 
            //                                data_type_t::OUTPUT, pe_array->get_stationary_type(), action_type_t::LOAD, true);
#endif
            if(!skip_transfer[data_type_t::OUTPUT]) {
                dynamic_power[data_type_t::OUTPUT] = u_dynamic_power[data_type_t::OUTPUT];
                num_data_transfer[data_type_t::OUTPUT]++;

                std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);
                parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);

                uint64_t address_pe_array = 0, address_global_buffer = 0;
                unsigned num_access_pe_array = 0, num_access_global_buffer = 0;

                for(unsigned b = 0; b < parameters_pe_array[parameter_type_t::BATCH_SIZE]; b++) {
                    for(unsigned k = 0; k < parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                        for(unsigned p = 0; p < parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                            for(unsigned q = 0; q < parameters_pe_array[parameter_type_t::OUTPUT_WIDTH]; q++) {
                                // Check address at PE array (temporal buffer)
                                if(address_pe_array != ((uint64_t)&pe_array->output_data[b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                          *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                          *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                         k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                          *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                         p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                         pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT]) {

                                    // Update PE array cost
                                    pe_array->access_energy[data_type_t::OUTPUT] += pe_array->u_write_energy[data_type_t::OUTPUT];
                                    pe_array->access_cycle[data_type_t::OUTPUT] += pe_array->u_write_cycle[data_type_t::OUTPUT];
                                    num_access_pe_array++;

                                    // Update PE array address (temporal buffer)
                                    address_pe_array = ((uint64_t)&pe_array->output_data[b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                          *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                          *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                         k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                          *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                         p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                         pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT];
                                }

                                // Check global buffer address 
                                if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                             m_scheduler->output_offset_global_buffer.front() + 
                                                                             b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                              *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                             k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                             p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                             mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT]) {
                                    // Update global buffer cost
                                    access_energy[data_type_t::OUTPUT] += u_read_energy[data_type_t::OUTPUT];
                                    access_cycle[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT];
                                    num_access_global_buffer++;

                                    // Update global buffer address 
                                    address_global_buffer = ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                             m_scheduler->output_offset_global_buffer.front() + 
                                                                             b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                              *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                             k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                             p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                             mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT];
                                }
                            }
                        }
                    }
                }
                // Update overlapped cost between PE array and global buffer
                unsigned ratio = ceil((double)(line_size[data_type_t::OUTPUT])/(double)(pe_array->line_size[data_type_t::OUTPUT]));

                // At the 1, 2, before last, and last stages
                unsigned first_stage = u_read_cycle[data_type_t::OUTPUT];
                unsigned second_stage = std::max(u_read_cycle[data_type_t::OUTPUT],
                                                 u_transfer_cycle*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth)));
                unsigned last_before_stage = std::max(ratio*pe_array->u_write_cycle[data_type_t::OUTPUT],
                                                      u_transfer_cycle*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth)));
                unsigned last_stage = ratio*pe_array->u_write_cycle[data_type_t::OUTPUT];

                // Remainder stages
                unsigned other_stage = std::max(u_read_cycle[data_type_t::OUTPUT],
                                       std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth)),
                                                ratio*pe_array->u_write_cycle[data_type_t::OUTPUT]));

                if(num_access_global_buffer == 1) {
                    cycle_pe_array_global_buffer[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT] +
                                                                         u_transfer_cycle*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth)) +
                                                                         ratio*pe_array->u_write_cycle[data_type_t::OUTPUT];
                } else {
                    cycle_pe_array_global_buffer[data_type_t::OUTPUT] += first_stage + second_stage +
                                                                         (num_access_global_buffer - 2)*other_stage +
                                                                         last_before_stage + last_stage;
                }

                // Update data transfer cycle and energy between Global buffer and PE array.
                transfer_cycle[data_type_t::OUTPUT] += num_access_global_buffer*u_transfer_cycle
                                                       *ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth));
                transfer_energy[data_type_t::OUTPUT] += num_access_global_buffer*u_transfer_energy
                                                        *ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth));


                pe_array->skip_transfer[data_type_t::OUTPUT] = false;
                pe_array->equal_output_tile = false;
                transfer_output = true;
            }
            if(pe_array->tile_size[data_type_t::OUTPUT] == tile_size[data_type_t::OUTPUT]) {skip_transfer[data_type_t::OUTPUT] = true;}
        }
        else {
            transfer_output = false;
            m_scheduler->output_read_global_buffer[m_scheduler->output_offset_global_buffer.front()] = true;
            for(auto it = m_scheduler->output_read_pe.begin(); it != m_scheduler->output_read_pe.end(); ++it) {
                it->second = false;
            }
        }

        if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y -1) {
            move_front(&m_scheduler->output_offset_global_buffer);
        }
        // Increment the output index.
        output_index++;

        // Update PE array signals
        pe_array->exist_data[data_type_t::OUTPUT] = true, pe_array->request_to_global_buffer[data_type_t::OUTPUT] = false;
    }
    pe_array->fill_data();

    // Check iterations of data
    // Case 1. Input stationary
    if(pe_array->get_stationary_type() == stationary_type_t::INPUT_STATIONARY) {
        // Reuse input data
        if(weight_index == m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front() &&
           output_index == m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front() &&
           input_index < m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front()) {

            weight_index = 0, output_index = 0;

            move_front(&m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT]);
            move_front(&m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT]);
        }
    }
    // Case 2. Weight stationary
    else if(pe_array->get_stationary_type() == stationary_type_t::WEIGHT_STATIONARY) {
        // Reuse weight
        if(input_index == m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front() &&
           output_index == m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front() &&
           weight_index < m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front()) {

            input_index = 0, output_index = 0;

            move_front(&m_scheduler->offset_size_global_buffer[data_type_t::INPUT]);
            move_front(&m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT]);
        }
    }
    // Case 3. Output stationary
    else if(pe_array->get_stationary_type() == stationary_type_t::OUTPUT_STATIONARY) {
        // Reuse output data
        if(input_index == m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front() &&
           weight_index == m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front() &&
           output_index < m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front()) {

            input_index = 0, weight_index = 0;

            move_front(&m_scheduler->offset_size_global_buffer[data_type_t::INPUT]);
            move_front(&m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT]);
        }
    }
    // Refresh all data
    if(input_index == m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front() &&
       weight_index == m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front() &&
       output_index == m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front()) {
        flush_data(m_scheduler);
    }
}


void global_buffer_t::flush_data(scheduler_t *m_scheduler) {
    // Case 1. Input stationary
    if(stationary_type == stationary_type_t::INPUT_STATIONARY) {
        // Case 1) Can reuse input data && request weight and output data
        if(weight_flush_counter < m_scheduler->offset_size_dram[data_type_t::WEIGHT].front() - 1 &&
           output_flush_counter < m_scheduler->offset_size_dram[data_type_t::OUTPUT].front() - 1) {
            // Weight and output data do not exist in Global buffer.
            exist_data[data_type_t::WEIGHT] = false, exist_data[data_type_t::OUTPUT] = false;

#ifdef FUNCTIONAL
            // Write back Output data 
            m_scheduler->transfer_data(multi_chip->data, data, m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()], offsets[data_type_t::OUTPUT],
                                       component_type_t::CHIPS_Y, component_type_t::GLOBAL_BUFFER, 
                                       data_type_t::OUTPUT, get_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
            parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);

            uint64_t address_global_buffer = 0, address_multi_chip = 0;
            unsigned num_access_global_buffer = 0, num_access_multi_chip = 0;

            for(unsigned b = 0; b < parameters_global_buffer[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            // Check global buffer address
                            if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                                   b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                   mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT]) {

                                // Update global buffer cost
                                access_energy[data_type_t::OUTPUT] += u_read_energy[data_type_t::OUTPUT];
                                access_cycle[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT];
                                num_access_global_buffer++;

                                // Update global buffer address
                                address_global_buffer = ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                         b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                         k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                         p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                         mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT];

                            }
                            // Check multi chip address (temporal buffer)
                            if(address_multi_chip != ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::OUTPUT] + 
                                                                                  m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()] + 
                                                                                  b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                  multi_chip->mask_bits[data_type_t::OUTPUT]) << multi_chip->mask_bits[data_type_t::OUTPUT]) {

                                // Update costs of chip-level processor (temporal buffer)
                                multi_chip->access_energy[data_type_t::OUTPUT] += multi_chip->u_write_energy[data_type_t::OUTPUT];
                                multi_chip->access_cycle[data_type_t::OUTPUT] += multi_chip->u_write_cycle[data_type_t::OUTPUT];
                                num_access_multi_chip++;

                                // Update address of chip-level processor (temporal buffer)
                                address_multi_chip = ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::OUTPUT] + 
                                                                                  m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()] + 
                                                                                  b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                  multi_chip->mask_bits[data_type_t::OUTPUT]) << multi_chip->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }


            /*
            // Update overlapped cost between PE array and global buffer
            unsigned ratio = ceil((double)(line_size[data_type_t::WEIGHT])/(double)(pe_array->line_size[data_type_t::WEIGHT]));

            // At the 1, 2, before last, and last stages
            unsigned first_stage = u_read_cycle[data_type_t::WEIGHT];
            unsigned second_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                             u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT]*8*sizeof(data_t))/(double)(bitwidth)));
            unsigned last_before_stage = std::max(ratio*pe_array->u_write_cycle[data_type_t::WEIGHT],
                                                  u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT]*8*sizeof(data_t))/(double)(bitwidth)));
            unsigned last_stage = ratio*pe_array->u_write_cycle[data_type_t::WEIGHT];

            // Remainder stages
            unsigned other_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                   std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT]*8*sizeof(data_t))/(double)(bitwidth)),
                                            ratio*pe_array->u_write_cycle[data_type_t::WEIGHT]));  
            */

            /*
            // Update Global buffer read cycle and energy.
            write_back_cycle += tile_size[data_type_t::OUTPUT]*u_read_cycle[data_type_t::OUTPUT]/(line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            // Update DRAM write cycle and energy.
            multi_chip->write_back_cycle += tile_size[data_type_t::OUTPUT]*multi_chip->u_write_cycle[data_type_t::OUTPUT]/(multi_chip->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            */

            // Increase flush counter of weight and output data
            weight_flush_counter++;
            output_flush_counter++;

            // Waiting for weight and output data.
            wait_data();
        }
        // Case 2) request all data type
        else {
            // Input data, weight, and output data do not exist in Global buffer.
            exist_data[data_type_t::INPUT] = false, exist_data[data_type_t::WEIGHT] = false, exist_data[data_type_t::OUTPUT] = false;

#ifdef FUNCTIONAL
            // Write back Output data 
            m_scheduler->transfer_data(multi_chip->data, data, m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()], offsets[data_type_t::OUTPUT],
                                       component_type_t::CHIPS_Y, component_type_t::GLOBAL_BUFFER, 
                                       data_type_t::OUTPUT, get_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
            parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);

            uint64_t address_global_buffer = 0, address_multi_chip = 0;
            unsigned num_access_global_buffer = 0, num_access_multi_chip = 0;

            for(unsigned b = 0; b < parameters_global_buffer[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            // Check global buffer address
                            if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                                   b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                   mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT]) {

                                // Update global buffer cost
                                access_energy[data_type_t::OUTPUT] += u_read_energy[data_type_t::OUTPUT];
                                access_cycle[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT];
                                num_access_global_buffer++;

                                // Update global buffer address
                                address_global_buffer = ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                         b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                         k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                         p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                         mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT];

                            }
                            // Check multi chip address (temporal buffer)
                            if(address_multi_chip != ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::OUTPUT] + 
                                                                                  m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()] + 
                                                                                  b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                  multi_chip->mask_bits[data_type_t::OUTPUT]) << multi_chip->mask_bits[data_type_t::OUTPUT]) {

                                // Update costs of chip-level processor (temporal buffer)
                                multi_chip->access_energy[data_type_t::OUTPUT] += multi_chip->u_write_energy[data_type_t::OUTPUT];
                                multi_chip->access_cycle[data_type_t::OUTPUT] += multi_chip->u_write_cycle[data_type_t::OUTPUT];
                                num_access_multi_chip++;

                                // Update address of chip-level processor (temporal buffer)
                                address_multi_chip = ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::OUTPUT] + 
                                                                                  m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()] + 
                                                                                  b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                  multi_chip->mask_bits[data_type_t::OUTPUT]) << multi_chip->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }

            /*
            // Update Global buffer read cycle and energy.
            write_back_cycle += tile_size[data_type_t::OUTPUT]*u_read_cycle[data_type_t::OUTPUT]/(line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            // Update DRAM write cycle and energy.
            multi_chip->write_back_cycle += tile_size[data_type_t::OUTPUT]*multi_chip->u_write_cycle[data_type_t::OUTPUT]/(multi_chip->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            */

            weight_flush_counter = 0; 
            output_flush_counter = 0;
            
            // Waiting for input data, weight, and output data.
            wait_data();
        }
    }
    // Case 2. Weight stationary
    else if(stationary_type == stationary_type_t::WEIGHT_STATIONARY) {
        // When not all input data and output data are transferred from DRAM to Global buffer.
        // Request input data and output data to DRAM.
        if(input_flush_counter < m_scheduler->offset_size_dram[data_type_t::INPUT].front() - 1 &&
           output_flush_counter < m_scheduler->offset_size_dram[data_type_t::OUTPUT].front() - 1) {
            // Input data and output data do not exist in Global buffer.
            exist_data[data_type_t::INPUT] = false, exist_data[data_type_t::OUTPUT] = false;
#ifdef FUNCTIONAL
            // Write back Output data 
            m_scheduler->transfer_data(multi_chip->data, data, m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()], offsets[data_type_t::OUTPUT],
                                       component_type_t::CHIPS_Y, component_type_t::GLOBAL_BUFFER, 
                                       data_type_t::OUTPUT, get_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
            parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);

            uint64_t address_global_buffer = 0, address_multi_chip = 0;
            unsigned num_access_global_buffer = 0, num_access_multi_chip = 0;

            for(unsigned b = 0; b < parameters_global_buffer[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            // Check global buffer address
                            if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                                   b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                   mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT]) {

                                // Update global buffer cost
                                access_energy[data_type_t::OUTPUT] += u_read_energy[data_type_t::OUTPUT];
                                access_cycle[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT];
                                num_access_global_buffer++;

                                // Update global buffer address
                                address_global_buffer = ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                         b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                         k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                         p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                         mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT];

                            }
                            // Check multi chip address (temporal buffer)
                            if(address_multi_chip != ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::OUTPUT] + 
                                                                                  m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()] + 
                                                                                  b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                  multi_chip->mask_bits[data_type_t::OUTPUT]) << multi_chip->mask_bits[data_type_t::OUTPUT]) {

                                // Update costs of chip-level processor (temporal buffer)
                                multi_chip->access_energy[data_type_t::OUTPUT] += multi_chip->u_write_energy[data_type_t::OUTPUT];
                                multi_chip->access_cycle[data_type_t::OUTPUT] += multi_chip->u_write_cycle[data_type_t::OUTPUT];
                                num_access_multi_chip++;

                                // Update address of chip-level processor (temporal buffer)
                                address_multi_chip = ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::OUTPUT] + 
                                                                                  m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()] + 
                                                                                  b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                  multi_chip->mask_bits[data_type_t::OUTPUT]) << multi_chip->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }

            /*
            // Update Global buffer read cycle and energy.
            write_back_cycle += tile_size[data_type_t::OUTPUT]*u_read_cycle[data_type_t::OUTPUT]/(line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            // Update DRAM write cycle and energy.
            multi_chip->write_back_cycle += tile_size[data_type_t::OUTPUT]*multi_chip->u_write_cycle[data_type_t::OUTPUT]/(multi_chip->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            */
            
            input_flush_counter++;
            output_flush_counter++;

            // Waiting for input data and output data.
            // Global buffer is in idle state.
            wait_data();
        }
        // Case 2) Request all data type
        else {
            // Input data, weight, and output data do not exist in Global buffer.
            exist_data[data_type_t::INPUT] = false, exist_data[data_type_t::WEIGHT] = false, exist_data[data_type_t::OUTPUT] = false;

#ifdef FUNCTIONAL
            // Write back Output data 
            m_scheduler->transfer_data(multi_chip->data, data, m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()], offsets[data_type_t::OUTPUT],
                                       component_type_t::CHIPS_Y, component_type_t::GLOBAL_BUFFER, 
                                       data_type_t::OUTPUT, get_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
            parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);

            uint64_t address_global_buffer = 0, address_multi_chip = 0;
            unsigned num_access_global_buffer = 0, num_access_multi_chip = 0;

            for(unsigned b = 0; b < parameters_global_buffer[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            // Check global buffer address
                            if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                                   b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                   mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT]) {

                                // Update global buffer cost
                                access_energy[data_type_t::OUTPUT] += u_read_energy[data_type_t::OUTPUT];
                                access_cycle[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT];
                                num_access_global_buffer++;

                                // Update global buffer address
                                address_global_buffer = ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                         b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                         k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                         p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                         mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT];

                            }
                            // Check multi chip address (temporal buffer)
                            if(address_multi_chip != ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::OUTPUT] + 
                                                                                  m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()] + 
                                                                                  b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                  multi_chip->mask_bits[data_type_t::OUTPUT]) << multi_chip->mask_bits[data_type_t::OUTPUT]) {

                                // Update costs of chip-level processor (temporal buffer)
                                multi_chip->access_energy[data_type_t::OUTPUT] += multi_chip->u_write_energy[data_type_t::OUTPUT];
                                multi_chip->access_cycle[data_type_t::OUTPUT] += multi_chip->u_write_cycle[data_type_t::OUTPUT];
                                num_access_multi_chip++;

                                // Update address of chip-level processor (temporal buffer)
                                address_multi_chip = ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::OUTPUT] + 
                                                                                  m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()] + 
                                                                                  b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                  multi_chip->mask_bits[data_type_t::OUTPUT]) << multi_chip->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }

            /*
            // Update Global buffer read cycle and energy.
            write_back_cycle += tile_size[data_type_t::OUTPUT]*u_read_cycle[data_type_t::OUTPUT/(line_size[data_type_t::OUTPUT]/8/sizeof(data_t))];
            // Update DRAM write cycle and energy.
            multi_chip->write_back_cycle += tile_size[data_type_t::OUTPUT]*multi_chip->u_write_cycle[data_type_t::OUTPUT]/(multi_chip->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            */
         
            input_flush_counter = 0;
            output_flush_counter = 0;

            // Waiting for input data, weight, and output data.
            // Global buffer is in idle state.
            wait_data();
        }
    }
    // Case 3. Output stationary
    else if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
        // Case 1) Reuse output data && Request input data and weight
        if(input_flush_counter < m_scheduler->offset_size_dram[data_type_t::INPUT].front() - 1 &&
           weight_flush_counter < m_scheduler->offset_size_dram[data_type_t::WEIGHT].front() - 1) {
            // Input data and weight do not exist in Global buffer.
            exist_data[data_type_t::INPUT] = false, exist_data[data_type_t::WEIGHT] = false;

            input_flush_counter++;
            weight_flush_counter++;

            // Waiting for input data and weight.
            // Global buffer is in idle state.
            wait_data();
        }
        // Case 2) Request all data types
        else {
            // Input data, weight, and output data do not exist in Global buffer.
            exist_data[data_type_t::INPUT] = false, exist_data[data_type_t::WEIGHT] = false, exist_data[data_type_t::OUTPUT] = false;
#ifdef FUNCTIONAL
            // Write back Output data 
            m_scheduler->transfer_data(multi_chip->data, data, m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()], offsets[data_type_t::OUTPUT],
                                       component_type_t::CHIPS_Y, component_type_t::GLOBAL_BUFFER, 
                                       data_type_t::OUTPUT, get_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
            parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);

            uint64_t address_global_buffer = 0, address_multi_chip = 0;
            unsigned num_access_global_buffer = 0, num_access_multi_chip = 0;

            for(unsigned b = 0; b < parameters_global_buffer[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            // Check global buffer address
                            if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                                   b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                   mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT]) {

                                // Update global buffer cost
                                access_energy[data_type_t::OUTPUT] += u_read_energy[data_type_t::OUTPUT];
                                access_cycle[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT];
                                num_access_global_buffer++;

                                // Update global buffer address
                                address_global_buffer = ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                         b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                         k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                         p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                         mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT];

                            }
                            // Check multi chip address (temporal buffer)
                            if(address_multi_chip != ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::OUTPUT] + 
                                                                                  m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()] + 
                                                                                  b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                  multi_chip->mask_bits[data_type_t::OUTPUT]) << multi_chip->mask_bits[data_type_t::OUTPUT]) {

                                // Update costs of chip-level processor (temporal buffer)
                                multi_chip->access_energy[data_type_t::OUTPUT] += multi_chip->u_write_energy[data_type_t::OUTPUT];
                                multi_chip->access_cycle[data_type_t::OUTPUT] += multi_chip->u_write_cycle[data_type_t::OUTPUT];
                                num_access_multi_chip++;

                                // Update address of chip-level processor (temporal buffer)
                                address_multi_chip = ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::OUTPUT] + 
                                                                                  m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()] + 
                                                                                  b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                  multi_chip->mask_bits[data_type_t::OUTPUT]) << multi_chip->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }

            /*
            // Update Global buffer read cycle and energy.
            write_back_cycle += tile_size[data_type_t::OUTPUT]*u_read_cycle[data_type_t::OUTPUT]/(line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            // Update DRAM write cycle and energy.
            multi_chip->write_back_cycle += tile_size[data_type_t::OUTPUT]*multi_chip->u_write_cycle[data_type_t::OUTPUT]/(multi_chip->line_size[data_type_t::OUTPUT]/8/sizeof(data_t));
            */

            // Reset flush counter of input data and weight
            input_flush_counter = 0;
            weight_flush_counter = 0;

            // Waiting for input data, weight, and output data.
            // Global buffer is in idle state.
            wait_data();
        }
    }
    // The counter of Global buffer should be initialized as 0.
    input_index = 0, weight_index = 0, output_index = 0;
}

void global_buffer_t::reset() {
    memset(data, 0.0, size);
    
    idle = false;
    initial = true;

    write_back_cycle = 0.0;
    overlapped_transfer_cycle = 0.0;

    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);

    num_request.assign(data_type_t::NUM_DATA_TYPES, 0);
    num_data_transfer.assign(data_type_t::NUM_DATA_TYPES, 0);

    access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    cycle_pe_array_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    utilization.assign(data_type_t::NUM_DATA_TYPES, 0.0);
}

/* Separated Global buffer */
separate_buffer_t::separate_buffer_t(section_config_t m_section_config) :
    global_buffer_t(m_section_config),
    input_size(1.0),
    weight_size(1.0),
    output_size(1.0) {

    init(m_section_config);
}

separate_buffer_t::~separate_buffer_t() {
    delete [] data;
}

// Initialize the Global buffer
void separate_buffer_t::init(section_config_t m_section_config) {

    /* Initialize separate buffer's specifications */

    // Initialize buffer size of input, weight, and output in KB
    m_section_config.get_setting("input_size", &input_size);
    m_section_config.get_setting("weight_size", &weight_size);
    m_section_config.get_setting("output_size", &output_size);
    input_size *= 1024, weight_size *= 1024, output_size *= 1024;
    size = input_size + weight_size + output_size;

    unsigned num_entry = (unsigned)size/sizeof(data_t);
    data = new data_t[num_entry]();
    memset(data, 1.0, num_entry*sizeof(data_t));

    // Initialize the frequency and bandwidth of the separate buffer
    m_section_config.get_setting("frequency", &frequency);
    m_section_config.get_setting("bandwidth", &bandwidth);
    bitwidth = 8*bandwidth/frequency;
    m_section_config.get_setting("bitwidth", &bitwidth);

    // Initialize global buffer type (double buffer or single buffer)
    m_section_config.get_setting("double_buffer", &double_buffer);

    bypass.reserve(data_type_t::NUM_DATA_TYPES);
    bypass.assign(data_type_t::NUM_DATA_TYPES, 0);
    m_section_config.get_vector_setting("bypass", &bypass);


    // Initialize line size and mask bits of the global buffer
    line_size.reserve(data_type_t::NUM_DATA_TYPES);
    line_size.assign(data_type_t::NUM_DATA_TYPES, sizeof(data_t));

    mask_bits.reserve(data_type_t::NUM_DATA_TYPES);
    mask_bits.assign(data_type_t::NUM_DATA_TYPES, 0);

    m_section_config.get_vector_setting("line_size", &line_size);

    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        while(line_size[i] > 8) {
            line_size[i] /= 2;
            mask_bits[i]++;
        }
    }
    m_section_config.get_vector_setting("line_size", &line_size);

    // Initialize the stationary type of the global buffer
    std::string stationary_str;
    if(m_section_config.get_setting("stationary_type", &stationary_str)) {
        stationary_type = (stationary_type_t)get_type(stationary_type_str, stationary_str);
    }

    if(stationary_type == stationary_type_t::INPUT_STATIONARY) {
        parameter_order = "BCPQKRS";
    } else if(stationary_type == stationary_type_t::WEIGHT_STATIONARY) {
        parameter_order = "KCRSBPQ";
    } else if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
        parameter_order = "BKPQKRS";
    }
    m_section_config.get_setting("parameter_order", &parameter_order);

    // Initialize the memory type.
    memory_type = memory_type_t::SEPARATE;

    skip_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);
    
    // Initialize the tile size of Global buffer.
    tile_size.reserve(data_type_t::NUM_DATA_TYPES);
    tile_size.assign(data_type_t::NUM_DATA_TYPES, 1);

    offsets.reserve(data_type_t::NUM_DATA_TYPES);
    offsets.assign(data_type_t::NUM_DATA_TYPES, 0);

    /* Initialize signals of the global buffer */

    // Initialize data exist signal
    exist_data.reserve(data_type_t::NUM_DATA_TYPES);
    exist_data.assign(data_type_t::NUM_DATA_TYPES, false);

    // Initialize data request signal
    request_to_multi_chip.reserve(data_type_t::NUM_DATA_TYPES);
    request_to_multi_chip.assign(data_type_t::NUM_DATA_TYPES, false);

    /* Initialize unit cost of the global buffer */

    // Initialize the unit cycle and energy of Global buffer
    u_read_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("read_cycle", &u_read_cycle);

    u_read_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("read_energy", &u_read_energy);

    u_write_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("write_cycle", &u_write_cycle);

    u_write_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("write_energy", &u_write_energy);

    u_dynamic_power.reserve(data_type_t::NUM_DATA_TYPES);
    u_dynamic_power.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("dynamic_power", &u_dynamic_power);
    
    u_static_power.reserve(data_type_t::NUM_DATA_TYPES);
    u_static_power.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("static_power", &u_static_power);

    m_section_config.get_setting("transfer_cycle", &u_transfer_cycle);
    m_section_config.get_setting("transfer_energy", &u_transfer_energy);

    /* Initialize stats of the global buffer */

    // Initialize the number of request to Global buffer
    num_request.reserve(data_type_t::NUM_DATA_TYPES);
    num_request.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the number of data transfer to PE array
    num_data_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the total access cycle
    access_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize the total access energy
    access_energy.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    cycle_pe_array_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    cycle_pe_array_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize the total transfer cycle between PE array and Global buffer
    transfer_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize the total transfer energy between PE array and Global buffer
    transfer_energy.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize dynamic power of global buffer
    dynamic_power.reserve(data_type_t::NUM_DATA_TYPES);
    dynamic_power.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    utilization.reserve(data_type_t::NUM_DATA_TYPES);
    utilization.assign(data_type_t::NUM_DATA_TYPES, 0.0);

}

void separate_buffer_t::update_offset() {
    offsets[data_type_t::INPUT] = 0;
    offsets[data_type_t::WEIGHT] = input_size/sizeof(data_t);
    offsets[data_type_t::OUTPUT] = input_size/sizeof(data_t) + weight_size/sizeof(data_t);

}

void separate_buffer_t::check_tile_size() {
    
    if(tile_size[data_type_t::INPUT]*sizeof(data_t) > input_size && !bypass[data_type_t::INPUT]) {
        std::cerr << "Input data size : " << tile_size[data_type_t::INPUT]*sizeof(data_t)
                  << " is bigger than input buffer size : " << input_size << std::endl;
        exit(1);
    }
    if(tile_size[data_type_t::WEIGHT]*sizeof(data_t) > weight_size && !bypass[data_type_t::WEIGHT]) {
        std::cerr << "Weight data size : " << tile_size[data_type_t::WEIGHT]*sizeof(data_t)
                  << " is bigger than weight buffer size : " << weight_size << std::endl;
        exit(1);
    }
    if(tile_size[data_type_t::OUTPUT]*sizeof(data_t) > output_size && !bypass[data_type_t::OUTPUT]) {
        std::cerr << "Output data size : " << tile_size[data_type_t::OUTPUT]*sizeof(data_t)
                  << " is bigger than output buffer size : " << output_size << std::endl;
        exit(1);
    }
}

double separate_buffer_t::get_dynamic_power() {
    double t_dynamic_power = 0.0;
    t_dynamic_power += dynamic_power[data_type_t::INPUT]
                   + dynamic_power[data_type_t::WEIGHT]
                   + dynamic_power[data_type_t::OUTPUT];
    return t_dynamic_power;
}

double separate_buffer_t::get_static_power() {
    double static_power = 0.0;
    static_power = u_static_power[data_type_t::INPUT] 
                 + u_static_power[data_type_t::WEIGHT] 
                 + u_static_power[data_type_t::OUTPUT];
    return static_power;
}


// Print out the stat of separate Global buffer.
void separate_buffer_t::print_specification() {
    std::cout << "============== Global buffer ===============" << std::endl;

    // Print out the buffer size.
    std::cout << "Global buffer type :" << std::setw(24)
                                        << "Separate" << std::endl;

    std::cout << "Input size         :" << std::setw(21) << std::setprecision(1)
                                        << input_size/1024 << " KB" << std::endl;
    std::cout << "Weight size        :" << std::setw(21) << std::setprecision(1)
                                        << weight_size/1024 << " KB" << std::endl;
    std::cout << "Output size        :" << std::setw(21) << std::setprecision(1)
                                        << output_size/1024 << " KB" << std::endl;

    // Print line size.
    std::cout << "Line size" << std::endl;
    std::cout << " * Input data      :" << std::setw(19) << std::setprecision(0)
                                        << line_size[data_type_t::INPUT] << " bits" << std::endl;
    std::cout << " * Weight          :" << std::setw(19) << std::setprecision(0)
                                        << line_size[data_type_t::WEIGHT] << " bits" << std::endl;
    std::cout << " * Output data     :" << std::setw(19) << std::setprecision(0)
                                        << line_size[data_type_t::OUTPUT] << " bits" << std::endl;    
    // Print out the stationary type.
    if(stationary_type == stationary_type_t::INPUT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24)
                                            << "Input stationary" << std::endl;
    }
    else if(stationary_type == stationary_type_t::WEIGHT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24)
                                            << "Weight stationary" << std::endl;
    } 
    else if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24)
                                            << "Output stationary" << std::endl;
    }
    std::cout << "Bandwidth          :" << std::setw(19) << std::setprecision(0)
                                        << bandwidth << " GB/s" << std::endl;
    
    std::cout << "Access energy (read/write)" << std::endl;
    std::cout << " * Input buffer    :" << std::setw(16) << std::setprecision(2)
                                        << u_read_energy[data_type_t::INPUT] << "/" << u_write_energy[data_type_t::INPUT] << " pJ" << std::endl;
    std::cout << " * Weight buffer   :" << std::setw(16) << std::setprecision(2)
                                        << u_read_energy[data_type_t::WEIGHT] << "/" << u_write_energy[data_type_t::WEIGHT] << " pJ" << std::endl;
    std::cout << " * Output buffer   :" << std::setw(16) << std::setprecision(2)
                                        << u_read_energy[data_type_t::OUTPUT] << "/" << u_write_energy[data_type_t::OUTPUT] << " pJ" << std::endl;

    std::cout << "Access cycle (read/write)" << std::endl;
    std::cout << " * Input buffer    :" << std::setw(13) << std::setprecision(1)
                                        << u_read_cycle[data_type_t::INPUT] << "/" << u_write_cycle[data_type_t::INPUT] << " cycles" << std::endl;
    std::cout << " * Weight buffer   :" << std::setw(13) << std::setprecision(1)
                                        << u_read_cycle[data_type_t::WEIGHT] << "/" << u_write_cycle[data_type_t::WEIGHT] << " cycles" << std::endl;
    std::cout << " * Output buffer   :" << std::setw(13) << std::setprecision(1)
                                        << u_read_cycle[data_type_t::OUTPUT] << "/" << u_write_cycle[data_type_t::OUTPUT] << " cycles" << std::endl;
    std::cout << std::endl;
}

shared_buffer_t::shared_buffer_t(section_config_t m_section_config) :
    global_buffer_t(m_section_config) {

    init(m_section_config);
}

shared_buffer_t::~shared_buffer_t() {
    delete [] data;
}

// Initialize the Global buffer.
void shared_buffer_t::init(section_config_t m_section_config) {

    /* Initialize shared buffer's specifications */

    // Initialize size of shared buffer in KB
    m_section_config.get_setting("memory_size", &size);
    // KB -> Byte
    size *= 1024;

    unsigned num_entry = size/sizeof(data_t);
    data = new data_t[num_entry]();
    memset(data, 1.0, num_entry*sizeof(data_t));

    // Initialize frequency and bandwidth of the shared buffer
    m_section_config.get_setting("frequency", &frequency);
    m_section_config.get_setting("bandwidth", &bandwidth);
    bitwidth = 8*bandwidth/frequency;
    m_section_config.get_setting("bitwidth", &bitwidth);

    m_section_config.get_setting("double_buffer", &double_buffer);

    bypass.reserve(data_type_t::NUM_DATA_TYPES);
    bypass.assign(data_type_t::NUM_DATA_TYPES, 0);
    m_section_config.get_vector_setting("bypass", &bypass);
    
    // Initialize line size and mask bits of the global buffer
    line_size.reserve(data_type_t::NUM_DATA_TYPES);
    line_size.assign(data_type_t::NUM_DATA_TYPES, sizeof(data_t));

    mask_bits.reserve(data_type_t::NUM_DATA_TYPES);
    mask_bits.assign(data_type_t::NUM_DATA_TYPES, 0);

    m_section_config.get_vector_setting("line_size", &line_size);

    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        while(line_size[i] > 8) {
            line_size[i] /= 2;
            mask_bits[i]++;
        }
    }
    m_section_config.get_vector_setting("line_size", &line_size);

    // Initialize the stationary type.
    std::string stationary_str;
    if(m_section_config.get_setting("stationary_type", &stationary_str)) {
        stationary_type = (stationary_type_t)get_type(stationary_type_str, stationary_str);
    }
    if(stationary_type == stationary_type_t::INPUT_STATIONARY) {
        parameter_order = "BCPQKRS";
    } else if(stationary_type == stationary_type_t::WEIGHT_STATIONARY) {
        parameter_order = "KCRSBPQ";
    } else if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
        parameter_order = "BKPQKRS";
    }
    m_section_config.get_setting("parameter_order", &parameter_order);

    // Initialize the memory type.
    memory_type = memory_type_t::SHARED;


    skip_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);

    // Initialize the tile size of the global buffer
    tile_size.reserve(data_type_t::NUM_DATA_TYPES);
    tile_size.assign(data_type_t::NUM_DATA_TYPES, 1);

    offsets.reserve(data_type_t::NUM_DATA_TYPES);
    offsets.assign(data_type_t::NUM_DATA_TYPES, 0);

    /* Initialize signals of shared buffer */

    // Initialize the signal to indicate whether the data exists or not.
    exist_data.reserve(data_type_t::NUM_DATA_TYPES);
    exist_data.assign(data_type_t::NUM_DATA_TYPES, false);

    // Initialize the request signal
    request_to_multi_chip.reserve(data_type_t::NUM_DATA_TYPES);
    request_to_multi_chip.assign(data_type_t::NUM_DATA_TYPES, false);

    /* Initialize the unit cost of the shared buffer */

    // Initialize the unit cycle and energy of Global buffer
    u_read_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("read_cycle", &u_read_cycle);

    u_read_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("read_energy", &u_read_energy);

    u_write_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("write_cycle", &u_write_cycle);

    u_write_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("write_energy", &u_write_energy);

    u_dynamic_power.reserve(data_type_t::NUM_DATA_TYPES);
    u_dynamic_power.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("dynamic_power", &u_dynamic_power);

    u_static_power.reserve(data_type_t::NUM_DATA_TYPES);
    u_static_power.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("static_power", &u_static_power);

    m_section_config.get_setting("transfer_cycle", &u_transfer_cycle);
    m_section_config.get_setting("transfer_energy", &u_transfer_energy);

    /* Initialize stats of the shared buffer */
 
    // Initialize the number of request 
    num_request.reserve(data_type_t::NUM_DATA_TYPES);
    num_request.assign(data_type_t::NUM_DATA_TYPES, 0);

    num_data_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the access cycle 
    access_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize the access energy.
    access_energy.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);   

    cycle_pe_array_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    cycle_pe_array_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize the data transfer cycle.
    transfer_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    
    // Initialize the data transfer energy.
    transfer_energy.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize dynamic power of global buffer
    dynamic_power.reserve(data_type_t::NUM_DATA_TYPES);
    dynamic_power.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    utilization.reserve(data_type_t::NUM_DATA_TYPES);
    utilization.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    
}

void shared_buffer_t::update_offset() {
    offsets[data_type_t::INPUT] = 0;
    offsets[data_type_t::WEIGHT] = tile_size[data_type_t::INPUT];
    offsets[data_type_t::OUTPUT] = tile_size[data_type_t::INPUT] + tile_size[data_type_t::WEIGHT];
}


void shared_buffer_t::check_tile_size() {
    unsigned data_size = tile_size[data_type_t::INPUT] + tile_size[data_type_t::WEIGHT] + tile_size[data_type_t::OUTPUT];

    if(bypass[data_type_t::INPUT]) {
        data_size -= tile_size[data_type_t::INPUT];
    }
    if(bypass[data_type_t::WEIGHT]) {
        data_size -= tile_size[data_type_t::WEIGHT];
    }
    if(bypass[data_type_t::OUTPUT]) {
        data_size -= tile_size[data_type_t::OUTPUT];
    }

    if(data_size*sizeof(data_t) > size) {
        std::cout << "The data size is bigger than Global buffer size\n"
                  << "Data : " << data_size*sizeof(data_t)
                  << " Buffer : " << size << std::endl;
        exit(1);
    }
}


double shared_buffer_t::get_dynamic_power() {
    double dynamic_power = 0.0;
    dynamic_power = u_dynamic_power[data_type_t::INPUT];
    return dynamic_power;
}

double shared_buffer_t::get_static_power() {
    double static_power = 0.0;
    static_power = u_static_power[data_type_t::INPUT];
    return static_power;
}


// Print out the stat of shared Global buffer.
void shared_buffer_t::print_specification() {
    std::cout << "============== Global buffer ===============" << std::endl;

    std::cout << "Global buffer type :" << std::setw(24) 
                                        << "Shared" << std::endl;
    std::cout << "Buffer size        :" << std::setw(21)
                                        << size/1024 << " KB" << std::endl;
    std::cout << "Line size          :" << std::setw(19) 
                                        << line_size[data_type_t::INPUT] << " bits" << std::endl;

    // Print out the stationary type.
    if(stationary_type == stationary_type_t::INPUT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24)
                                            << "Input stationary" << std::endl;
    }
    else if(stationary_type == stationary_type_t::WEIGHT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24)                     
                                            << "Weight stationary" << std::endl;
    } 
    else {
        std::cout << "Stationary type    :" << std::setw(24)                     
                                            << "Output stationary" << std::endl;
    }
    std::cout << "Bandwidth          :" << std::setw(19) << std::setprecision(0)
                                        << bandwidth << " GB/s" << std::endl;
    
    std::cout << "Access energy (read/write)" << std::endl;
    std::cout << " * Global buffer   :" << std::setw(16) << std::setprecision(2) 
                                        << u_read_energy[data_type_t::INPUT] << "/" << u_write_energy[data_type_t::INPUT] << " pJ" << std::endl;
    
    std::cout << "Access cycle (read/write) " << std::endl;
    std::cout << " * Global buffer    :" << std::setw(12) << std::setprecision(1)
                                        << u_read_cycle[data_type_t::INPUT] << "/" << u_write_cycle[data_type_t::INPUT] << " cycles" << std::endl;
    std::cout << std::endl;
}

