#include <cstring>
#include <cmath>
#include "dram.h"


dram_t::dram_t(section_config_t m_section_config) :
    u_transfer_cycle(0.0),
    u_transfer_energy(0.0),
    transfer_output(false),
    data(NULL),
    size(0),
    frequency(0.0),
    bandwidth(0.0),
    bitwidth(0.0),
    input_index(0),
    weight_index(0),
    output_index(0),
    done(false),
    multi_chip(NULL) {
    
    init(m_section_config);
}


dram_t::~dram_t() {

#ifdef DRAMSIM3
    delete memory;
#endif
}


void dram_t::init(section_config_t m_section_config) {

    // Initialize the size of off-chip memory in GB.
    m_section_config.get_setting("size", &size);

    size *= 1024*1024*1024; 
    size_t num_entry = size/sizeof(data_t);
    data = new data_t[num_entry]();
    npu_mmu::npu_malloc(num_entry);

    // Initialize frequency, bandwidth, and bitwidth of the off-chip memory
    m_section_config.get_setting("frequency", &frequency);
    m_section_config.get_setting("bandwidth", &bandwidth);
    bitwidth = 8*bandwidth/frequency;
    m_section_config.get_setting("bitwidth", &bitwidth);

    // Initialize off-chip memory line size.
    line_size.reserve(data_type_t::NUM_DATA_TYPES);
    line_size.assign(data_type_t::NUM_DATA_TYPES, sizeof(data_t));
    m_section_config.get_vector_setting("line_size", &line_size);
    mask_bits.reserve(data_type_t::NUM_DATA_TYPES);
    mask_bits.assign(data_type_t::NUM_DATA_TYPES, 0);
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        while(line_size[i] > 8) {
            line_size[i] /= 2;
            mask_bits[i]++;
        }
    }
    m_section_config.get_vector_setting("line_size", &line_size);

    // Initialize the tile size
    tile_size.reserve(data_type_t::NUM_DATA_TYPES);
    tile_size.assign(data_type_t::NUM_DATA_TYPES, 1);

    offsets.reserve(data_type_t::NUM_DATA_TYPES);
    offsets.assign(data_type_t::NUM_DATA_TYPES, 1);

    skip_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);

    /* Off-chip memory stats */

    // Initialize the number of request at DRAM.
    num_request.reserve(data_type_t::NUM_DATA_TYPES);
    num_request.assign(data_type_t::NUM_DATA_TYPES, 0);
    
    num_data_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the access cycle of DRAM.
    access_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize the access energy of DRAM.
    access_energy.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize the transfer cycle between Global buffer and DRAM.
    transfer_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize the transfer energy between Global buffer and DRAM.
    transfer_energy.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize total cycles at the off-chip memory.
    cycle_chip_dram.reserve(data_type_t::NUM_DATA_TYPES);
    cycle_chip_dram.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    /* The unit stats */ 
    // Initialize DRAM transfer cycle and energy
    m_section_config.get_setting("transfer_cycle", &u_transfer_cycle);
    m_section_config.get_setting("transfer_energy", &u_transfer_energy);

    // Initialize DRAM unit read cycle
    u_read_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("read_cycle", &u_read_cycle);

    // Initialize DRAM unit read energy
    u_read_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("read_energy", &u_read_energy);

    // Initialize DRAM unit write cycle
    u_write_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("write_cycle", &u_write_cycle);

    // Initialize DRAM unit write energy
    u_write_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("write_energy", &u_write_energy);

#ifdef DRAMSIM3
    // Initialize DRAMsim3 configuration.
    m_section_config.get_setting("dram_config", &dram_config);
    m_section_config.get_setting("output_dir", &output_dir);

    // Initialize the DRAMsim3 wrapper.
    memory = new memory_controller_t(dram_config, output_dir);
#endif

}


// Connect DRAM to Global buffer.
void dram_t::connect(multi_chip_t *m_multi_chip) {
    multi_chip = m_multi_chip;
}

/* TODO */
// Copy neural data from PyTorch for Functional simulation
void dram_t::connect_layer(layer_t *m_layer) {
    layer = m_layer;
}

void dram_t::disconnect_layer() {
}

void dram_t::update_tile_size(scheduler_t *m_scheduler) {
    tile_size = m_scheduler->tile_size[component_type_t::DRAM];

    update_offset();
    check_tile_size();
}

void dram_t::update_offset() {
    offsets[data_type_t::INPUT] = 0;
    offsets[data_type_t::WEIGHT] = tile_size[data_type_t::INPUT];
    offsets[data_type_t::OUTPUT] = tile_size[data_type_t::INPUT] + tile_size[data_type_t::WEIGHT];

}

void dram_t::check_tile_size() {
    size_t data_size = tile_size[data_type_t::INPUT]
                     + tile_size[data_type_t::WEIGHT]
                     + tile_size[data_type_t::OUTPUT];
    if(data_size*sizeof(data_t) > size) {
        std::cerr << "The data size is bigger than the off-chip memory size\n" 
                  << "Data : " << data_size*sizeof(data_t) 
                  << ", Buffer : " << size << std::endl;
        exit(1);
    }
}

unsigned dram_t::get_bitwidth() {
    return bitwidth;
}

bool dram_t::is_idle() {
    return done;
}

void dram_t::data_transfer(scheduler_t *m_scheduler) {
    if(multi_chip->request_to_dram[data_type_t::INPUT]) {
#ifdef DRAMSIM3
        //send_request((data_t*)layer->input_data, m_scheduler->mapping_table, m_scheduler->input_offset_dram.front(), data_type_t::INPUT, action_type_t::LOAD);
        send_request((data_t*)data, m_scheduler->mapping_table, offsets[data_type_t::INPUT] + m_scheduler->input_offset_dram.front(), data_type_t::INPUT, action_type_t::LOAD);
#endif

#ifdef FUNCTIONAL
        // Input data transfer
        m_scheduler->transfer_data(multi_chip->data, (data_t*)layer->input_data, multi_chip->offsets[data_type_t::INPUT], m_scheduler->input_offset_dram.front(), 
                                   component_type_t::CHIPS_Y, component_type_t::DRAM, 
                                   data_type_t::INPUT, multi_chip->get_stationary_type(), action_type_t::LOAD);
        // Update for NPUsim ver2
        //m_scheduler->transfer_data_ver2(multi_chip->data, (data_t*)layer->input_data, 
        //                                component_type_t::CHIPS_Y, component_type_t::DRAM,
        //                                data_type_t::INPUT, multi_chip->get_stationary_type(), 
        //                                action_type_t::LOAD, true);

        // Case 1. Dense data format
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::INPUT]) {
                num_data_transfer[data_type_t::INPUT]++;
                std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_dram(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);
                parameters_dram = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::DRAM);

                uint64_t address_multi_chip = 0, address_dram = 0;
                unsigned num_access_multi_chip = 0, num_access_dram = 0;

                for(unsigned b = 0; b < parameters_multi_chip[parameter_type_t::BATCH_SIZE]; b++) {
                    for(unsigned c = 0; c < parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned h = 0; h < parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]; h++) {
                            for(unsigned w = 0; w < parameters_multi_chip[parameter_type_t::INPUT_WIDTH]; w++) {
                                // Check address of multi-chip processor
                                if(address_multi_chip != ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::INPUT] + 
                                                                                      b*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                                       *parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                                       *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] +
                                                                                      c*parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                                       *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                                      h*parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                      multi_chip->mask_bits[data_type_t::INPUT]) << multi_chip->mask_bits[data_type_t::INPUT]) {

                                    // Update Multi-chip processor cost (Temporal buffer)
                                    multi_chip->access_energy[data_type_t::INPUT] += multi_chip->u_write_energy[data_type_t::INPUT];
                                    if(!multi_chip->double_buffer) {
                                        multi_chip->access_cycle[data_type_t::INPUT] += multi_chip->u_write_cycle[data_type_t::INPUT];
                                    }
                                    num_access_multi_chip++;

                                    // Update Multi-chip processor address
                                    address_multi_chip = ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::INPUT] + 
                                                                                      b*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                                       *parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                                       *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                                      c*parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                                       *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                                      h*parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                      multi_chip->mask_bits[data_type_t::INPUT]) << multi_chip->mask_bits[data_type_t::INPUT];
                                }

                                // Check DRAM address
                                if(address_dram != ((uint64_t)&layer->input_data[m_scheduler->input_offset_dram.front() +
                                                                                 b*parameters_dram[parameter_type_t::INPUT_CHANNEL]
                                                                                  *parameters_dram[parameter_type_t::INPUT_HEIGHT]
                                                                                  *parameters_dram[parameter_type_t::INPUT_WIDTH] + 
                                                                                 c*parameters_dram[parameter_type_t::INPUT_HEIGHT]
                                                                                  *parameters_dram[parameter_type_t::INPUT_WIDTH] + 
                                                                                 h*parameters_dram[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                 mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT]) {
                                    // Update DRAM cost 
                                    access_energy[data_type_t::INPUT] += u_read_energy[data_type_t::INPUT];
                                    if(!multi_chip->double_buffer) {
                                        access_cycle[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT];
                                    }
                                    num_access_dram++;
                                    
                                    address_dram = ((uint64_t)&layer->input_data[m_scheduler->input_offset_dram.front() +
                                                                                 b*parameters_dram[parameter_type_t::INPUT_CHANNEL]
                                                                                  *parameters_dram[parameter_type_t::INPUT_HEIGHT]
                                                                                  *parameters_dram[parameter_type_t::INPUT_WIDTH] + 
                                                                                 c*parameters_dram[parameter_type_t::INPUT_HEIGHT]
                                                                                  *parameters_dram[parameter_type_t::INPUT_WIDTH] + 
                                                                                 h*parameters_dram[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                 mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT];
                                }
                            }
                        }
                    }
                }

                unsigned ratio = ceil((double)(line_size[data_type_t::INPUT])/(double)(multi_chip->line_size[data_type_t::INPUT]));

                if(multi_chip->exist_temporal_buffer) {
                    // At the 1, 2, before last, and last stages
                    unsigned first_stage = u_read_cycle[data_type_t::INPUT];
                    unsigned second_stage = std::max(u_read_cycle[data_type_t::INPUT], 
                                                     u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)));
                    unsigned last_before_stage = std::max(ratio*multi_chip->u_write_cycle[data_type_t::INPUT], 
                                                          u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)));
                    unsigned last_stage = ratio*multi_chip->u_write_cycle[data_type_t::INPUT];

                    // Remainder stages
                    unsigned other_stage = std::max(u_read_cycle[data_type_t::INPUT],
                                           std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)),
                                                    ratio*multi_chip->u_write_cycle[data_type_t::INPUT]));

                    // Update overlapped cycle between off-chip memory and multi-chip processor
                    if(num_access_dram == 1) {
                        cycle_chip_dram[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT] +
                                                               u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)) +
                                                               ratio*multi_chip->u_write_cycle[data_type_t::INPUT];
                    } else {
                        cycle_chip_dram[data_type_t::INPUT] += first_stage + second_stage +
                                                               (num_access_dram-2)*other_stage + 
                                                               last_before_stage + last_stage;
                    }
                }
                // Update data transfer cycle and energy between DRAM and Global buffer. 
                transfer_cycle[data_type_t::INPUT] += num_access_dram*u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth));
                transfer_energy[data_type_t::INPUT] += num_access_dram*u_transfer_energy*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth));
                multi_chip->skip_transfer[data_type_t::INPUT] = false;
            }
        }
        // Case 2. COO data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
            if(!skip_transfer[data_type_t::INPUT]) {
                std::cout << "Current version does not support COO format for sparse data" << std::endl;
                exit(1);
            }
        }
        // Case 3. CSC data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::INPUT]) {
                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::CHIPS_Y);
                unsigned row = parameters[parameter_type_t::INPUT_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer[data_type_t::INPUT]++;

                // Update off-chip memory access cost
                access_cycle[data_type_t::INPUT] += (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT]  + // Non-zero data
                                                    (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT]
                                                   /(sizeof(data_t)*8/row_bit) + // Row index
                                                    parameters[parameter_type_t::BATCH_SIZE]
                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                   *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                   *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/row_bit); // Column pointer
                access_energy[data_type_t::INPUT] += (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT] + // Non-zero data
                                                     (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT]
                                                    /(sizeof(data_t)*8/row_bit) + // Row index
                                                     parameters[parameter_type_t::BATCH_SIZE]
                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                    *u_read_energy[data_type_t::INPUT]/(sizeof(data_t)*8/row_bit); // Column pointer

                // Update on-chip processor access cost
                multi_chip->access_cycle[data_type_t::INPUT] += (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                                (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT]
                                                               /(sizeof(data_t)*8/row_bit) + // row index
                                                                parameters[parameter_type_t::BATCH_SIZE]
                                                               *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/row_bit); // Column pointer
                multi_chip->access_energy[data_type_t::INPUT] += (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                *multi_chip->u_write_energy[data_type_t::INPUT] + // Non-zero data
                                                                 (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                *multi_chip->u_write_energy[data_type_t::INPUT]
                                                                /(sizeof(data_t)*8/row_bit) + // Row index
                                                                 parameters[parameter_type_t::BATCH_SIZE]
                                                                *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                                *multi_chip->u_write_energy[data_type_t::INPUT]/(sizeof(data_t)*8/row_bit); // Column pointer

                // Update overlapped cycle between the off-chip memory and temporal buffer in on-chip processor
                cycle_chip_dram[data_type_t::INPUT] += std::max((multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                                (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *u_read_cycle[data_type_t::INPUT]
                                                               /(sizeof(data_t)*8/row_bit) + // Row index
                                                                parameters[parameter_type_t::BATCH_SIZE]
                                                               *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                               *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/row_bit), // Column pointer
                                                                (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                                (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT]
                                                               /(sizeof(data_t)*8/row_bit) + // Row index
                                                                parameters[parameter_type_t::BATCH_SIZE]
                                                               *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/row_bit)); // Column pointer

                // Update transfer cost between the off-chip memory and temporal buffer in on-chip processor
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((float)((multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                      u_transfer_cycle*ceil((float)((multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)/(float)bitwidth) + // Row index
                                                      u_transfer_cycle*ceil((float)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                   *(parameters[parameter_type_t::INPUT_WIDTH+1])*row_bit)/(float)bitwidth); // Column pointer
                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((float)((multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((float)((multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)/(float)bitwidth) + // Row index
                                                       u_transfer_energy*ceil((float)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                     *(parameters[parameter_type_t::INPUT_WIDTH+1])*row_bit)/(float)bitwidth); // Column pointer
                multi_chip->skip_transfer[data_type_t::INPUT] = false;

            }
        }
        // Case 4. CSR data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::INPUT]) {
                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::CHIPS_Y);
                unsigned column = parameters[parameter_type_t::INPUT_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer[data_type_t::INPUT]++;

                // Update off-chip memory access cost
                access_cycle[data_type_t::INPUT] += (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT]  + // Non-zero data
                                                    (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT]
                                                   /(sizeof(data_t)*8/column_bit) + // Column index
                                                    parameters[parameter_type_t::BATCH_SIZE]
                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                   *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                   *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/column_bit); // Row pointer
                access_energy[data_type_t::INPUT] += (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT] + // Non-zero data
                                                     (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT]
                                                    /(sizeof(data_t)*8/column_bit) + // Column index
                                                     parameters[parameter_type_t::BATCH_SIZE]
                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                    *u_read_energy[data_type_t::INPUT]/(sizeof(data_t)*8/column_bit); // Row pointer

                // Update on-chip processor access cost
                multi_chip->access_cycle[data_type_t::INPUT] += (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                                (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT]
                                                               /(sizeof(data_t)*8/column_bit) + // Column index
                                                                parameters[parameter_type_t::BATCH_SIZE]
                                                               *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/column_bit); // Row pointer
                multi_chip->access_energy[data_type_t::INPUT] += (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                *multi_chip->u_write_energy[data_type_t::INPUT] + // Non-zero data
                                                                 (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                *multi_chip->u_write_energy[data_type_t::INPUT]
                                                                /(sizeof(data_t)*8/column_bit) + // Column index
                                                                 parameters[parameter_type_t::BATCH_SIZE]
                                                                *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                                *multi_chip->u_write_energy[data_type_t::INPUT]/(sizeof(data_t)*8/column_bit); // Row pointer

                // Update overlapped cycle between the off-chip memory and temporal buffer in on-chip processor
                cycle_chip_dram[data_type_t::INPUT] += std::max((multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                                (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *u_read_cycle[data_type_t::INPUT]
                                                               /(sizeof(data_t)*8/column_bit) + // Column index
                                                                parameters[parameter_type_t::BATCH_SIZE]
                                                               *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                               *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/column_bit), // Row pointer
                                                                (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                                (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT]
                                                               /(sizeof(data_t)*8/column_bit) + // Column index
                                                                parameters[parameter_type_t::BATCH_SIZE]
                                                               *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/column_bit)); // Row pointer

                // Update transfer cost between the off-chip memory and temporal buffer in on-chip processor
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((float)((multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                      u_transfer_cycle*ceil((float)((multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)/(float)bitwidth) + // Column index
                                                      u_transfer_cycle*ceil((float)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                   *(parameters[parameter_type_t::INPUT_HEIGHT+1])*column_bit)/(float)bitwidth); // Row pointer
                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((float)((multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((float)((multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)/(float)bitwidth) + // Column index
                                                       u_transfer_energy*ceil((float)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                     *(parameters[parameter_type_t::INPUT_HEIGHT+1])*column_bit)/(float)bitwidth); // Row pointer
                multi_chip->skip_transfer[data_type_t::INPUT] = false;
            }
        }
        // Case 4. SparseMap
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::INPUT]) {
                num_data_transfer[data_type_t::INPUT]++;
    
                // Update off-chip memory access cost
                access_cycle[data_type_t::INPUT] += (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                    multi_chip->tile_size[data_type_t::INPUT]
                                                    *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
                access_energy[data_type_t::INPUT] += (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT] + // Non-zero data
                                                     multi_chip->tile_size[data_type_t::INPUT]
                                                    *u_read_energy[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
    
                // Update on-chip processor access cost (if temporal buffer exist)
                multi_chip->access_cycle[data_type_t::INPUT] += (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                                multi_chip->tile_size[data_type_t::INPUT]
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
                multi_chip->access_energy[data_type_t::INPUT] += (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                *multi_chip->u_write_energy[data_type_t::INPUT] + // Non-zero data
                                                                 multi_chip->tile_size[data_type_t::INPUT]
                                                                *multi_chip->u_write_energy[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
    
                // Update overlapped cycle between the off-chip memory and on-chip processor
                cycle_chip_dram[data_type_t::INPUT] += std::max((multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                                multi_chip->tile_size[data_type_t::INPUT]
                                                               *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8), // metadata
                                                                (multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                                multi_chip->tile_size[data_type_t::INPUT]
                                                               *multi_chip->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8)); // meta data
    
                // Update transfer cost between the off-chip memory and on-chip processor
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((float)((multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                     *8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                      u_transfer_cycle*ceil((float)(multi_chip->tile_size[data_type_t::INPUT])/(float)bitwidth); // Metadata
                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((float)((multi_chip->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                      *8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((float)(multi_chip->tile_size[data_type_t::INPUT])/(float)bitwidth); // Metadata

                multi_chip->skip_transfer[data_type_t::INPUT] = false;

            }
        }
        else {
            std::cerr << "Undefined data format" << std::endl;
            exit(1);
        }
        move_front(&m_scheduler->input_offset_dram);

#else
        if(!skip_transfer[data_type_t::INPUT]) {
            num_data_transfer[data_type_t::INPUT]++;

            std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_dram(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);
            parameters_dram = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::DRAM);

            uint64_t address_multi_chip = 0, address_dram = 0;
            unsigned num_access_multi_chip = 0, num_access_dram = 0;

            for(unsigned b = 0; b < parameters_multi_chip[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned c = 0; c < parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]; c++) {
                    for(unsigned h = 0; h < parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]; h++) {
                        for(unsigned w = 0; w < parameters_multi_chip[parameter_type_t::INPUT_WIDTH]; w++) {
                            // Check address of multi-chip processor
                            if(address_multi_chip != ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::INPUT] + 
                                                                                  b*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] +
                                                                                  c*parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                                  h*parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                  multi_chip->mask_bits[data_type_t::INPUT]) << multi_chip->mask_bits[data_type_t::INPUT]) {

                                // Update Multi-chip processor cost (Temporal buffer)
                                multi_chip->access_energy[data_type_t::INPUT] += multi_chip->u_write_energy[data_type_t::INPUT];
                                if(!multi_chip->double_buffer) {
                                    multi_chip->access_cycle[data_type_t::INPUT] += multi_chip->u_write_cycle[data_type_t::INPUT];
                                }
                                num_access_multi_chip++;

                                // Update Multi-chip processor address
                                address_multi_chip = ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::INPUT] + 
                                                                                  b*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                                  c*parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                                  h*parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                  multi_chip->mask_bits[data_type_t::INPUT]) << multi_chip->mask_bits[data_type_t::INPUT];
                            }

                            // Check DRAM address
                            if(address_dram != ((uint64_t)&layer->input_data[m_scheduler->input_offset_dram.front() +
                                                                             b*parameters_dram[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_dram[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_dram[parameter_type_t::INPUT_WIDTH] + 
                                                                             c*parameters_dram[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_dram[parameter_type_t::INPUT_WIDTH] + 
                                                                             h*parameters_dram[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                             mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT]) {
                                // Update DRAM cost 
                                access_energy[data_type_t::INPUT] += u_read_energy[data_type_t::INPUT];
                                if(!multi_chip->double_buffer) {
                                    access_cycle[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT];
                                }
                                num_access_dram++;
                                
                                address_dram = ((uint64_t)&layer->input_data[m_scheduler->input_offset_dram.front() +
                                                                             b*parameters_dram[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_dram[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_dram[parameter_type_t::INPUT_WIDTH] + 
                                                                             c*parameters_dram[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_dram[parameter_type_t::INPUT_WIDTH] + 
                                                                             h*parameters_dram[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                             mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT];
                            }
                        }
                    }
                }
            }

            unsigned ratio = ceil((double)(line_size[data_type_t::INPUT])/(double)(multi_chip->line_size[data_type_t::INPUT]));

            if(multi_chip->exist_temporal_buffer) {
                // At the 1, 2, before last, and last stages
                unsigned first_stage = u_read_cycle[data_type_t::INPUT];
                unsigned second_stage = std::max(u_read_cycle[data_type_t::INPUT], 
                                                 u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)));
                unsigned last_before_stage = std::max(ratio*multi_chip->u_write_cycle[data_type_t::INPUT], 
                                                      u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)));
                unsigned last_stage = ratio*multi_chip->u_write_cycle[data_type_t::INPUT];

                // Remainder stages
                unsigned other_stage = std::max(u_read_cycle[data_type_t::INPUT],
                                       std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)),
                                                ratio*multi_chip->u_write_cycle[data_type_t::INPUT]));

                // Update overlapped cycle between off-chip memory and multi-chip processor
                if(num_access_dram == 1) {
                    cycle_chip_dram[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT] +
                                                           u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)) +
                                                           ratio*multi_chip->u_write_cycle[data_type_t::INPUT];
                } else {
                    cycle_chip_dram[data_type_t::INPUT] += first_stage + second_stage +
                                                           (num_access_dram-2)*other_stage + 
                                                           last_before_stage + last_stage;
                }
            }

            // Update data transfer cycle and energy between DRAM and Global buffer. 
            transfer_cycle[data_type_t::INPUT] += num_access_dram*u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth));
            transfer_energy[data_type_t::INPUT] += num_access_dram*u_transfer_energy*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth));
            multi_chip->skip_transfer[data_type_t::INPUT] = false;
        /* Stats */
        }
#endif
        // Increment the input index
        input_index++;
        multi_chip->exist_data[data_type_t::INPUT] = true, multi_chip->request_to_dram[data_type_t::INPUT] = false;
        if(multi_chip->tile_size[data_type_t::INPUT] == tile_size[data_type_t::INPUT]) { skip_transfer[data_type_t::INPUT] = true;}
    }
    if(multi_chip->request_to_dram[data_type_t::WEIGHT]) {

#ifdef DRAMSIM3
        //send_request((data_t*)layer->weight, m_scheduler->mapping_table, m_scheduler->weight_offset_dram.front(), data_type_t::WEIGHT, action_type_t::LOAD);
        send_request((data_t*)data, m_scheduler->mapping_table, offsets[data_type_t::WEIGHT] + m_scheduler->weight_offset_dram.front(), data_type_t::WEIGHT, action_type_t::LOAD);
#endif

#ifdef FUNCTIONAL
        // Weight transfer
        m_scheduler->transfer_data(multi_chip->data, (data_t*)layer->weight, multi_chip->offsets[data_type_t::WEIGHT], m_scheduler->weight_offset_dram.front(), 
                                   component_type_t::CHIPS_Y, component_type_t::DRAM, 
                                   data_type_t::WEIGHT, multi_chip->get_stationary_type(), action_type_t::LOAD);
        // Update for NPUsim ver2
        //m_scheduler->transfer_data_ver2(multi_chip->data, (data_t*)layer->weight,
        //                                component_type_t::CHIPS_Y, component_type_t::DRAM,
        //                                data_type_t::WEIGHT, multi_chip->get_stationary_type(), 
        //                                action_type_t::LOAD, true);
        // Case 1. Dense data format
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                num_data_transfer[data_type_t::WEIGHT]++;

                std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_dram(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);
                parameters_dram = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::DRAM);

                uint64_t address_multi_chip = 0, address_dram = 0;
                unsigned num_access_multi_chip = 0, num_access_dram = 0;

                for(unsigned k = 0; k < parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned c = 0; c < parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned r = 0; r < parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]; r++) {
                            for(unsigned s = 0; s < parameters_multi_chip[parameter_type_t::FILTER_WIDTH]; s++) {
                                if(address_multi_chip != ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::WEIGHT] + 
                                                                                      k*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                                       *parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                                       *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] +
                                                                                      c*parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                                       *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                                      r*parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                      multi_chip->mask_bits[data_type_t::WEIGHT]) << multi_chip->mask_bits[data_type_t::WEIGHT]) {

                                    // Update multi-chip processor cost (temporal buffer)
                                    multi_chip->access_energy[data_type_t::WEIGHT] += multi_chip->u_write_energy[data_type_t::WEIGHT];
                                    if(!multi_chip->double_buffer) {
                                        multi_chip->access_cycle[data_type_t::WEIGHT] += multi_chip->u_write_cycle[data_type_t::WEIGHT];
                                    }
                                    num_access_multi_chip++;

                                    // Update multi-chip processor address
                                    address_multi_chip = ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::WEIGHT] + 
                                                                                      k*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                                       *parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                                       *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                                      c*parameters_multi_chip[parameter_type_t::FILTER_HEIGHT] 
                                                                                       *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                                      r*parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                      multi_chip->mask_bits[data_type_t::WEIGHT]) << multi_chip->mask_bits[data_type_t::WEIGHT];
                                }

                                if(address_dram != ((uint64_t)&layer->weight[m_scheduler->weight_offset_dram.front() +
                                                                             k*parameters_dram[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_dram[parameter_type_t::FILTER_HEIGHT]
                                                                              *parameters_dram[parameter_type_t::FILTER_WIDTH] + 
                                                                             c*parameters_dram[parameter_type_t::FILTER_HEIGHT]
                                                                              *parameters_dram[parameter_type_t::FILTER_WIDTH] + 
                                                                             r*parameters_dram[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                             mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT]) {
                                    // Update off-chip memory cost
                                    access_energy[data_type_t::WEIGHT] += u_read_energy[data_type_t::WEIGHT];
                                    if(!multi_chip->double_buffer) {
                                        access_cycle[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT];
                                    }
                                    num_access_dram++;

                                    // Update off-chip memory address
                                    address_dram = ((uint64_t)&layer->weight[m_scheduler->weight_offset_dram.front() +
                                                                             k*parameters_dram[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_dram[parameter_type_t::FILTER_HEIGHT]
                                                                              *parameters_dram[parameter_type_t::FILTER_WIDTH] + 
                                                                             c*parameters_dram[parameter_type_t::FILTER_HEIGHT]
                                                                              *parameters_dram[parameter_type_t::FILTER_WIDTH] + 
                                                                             r*parameters_dram[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                             mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT];
                                }
                            }
                        }
                    }
                }

                unsigned ratio = ceil((double)(line_size[data_type_t::WEIGHT])/(double)(multi_chip->line_size[data_type_t::WEIGHT]));

                // At the 1, 2, before last, and last stages
                double first_stage = u_read_cycle[data_type_t::WEIGHT];
                double second_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                                 u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
                double last_before_stage = std::max(ratio*multi_chip->u_write_cycle[data_type_t::WEIGHT],
                                                      u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
                double last_stage = ratio*multi_chip->u_write_cycle[data_type_t::WEIGHT];

                // Remainder stages
                double other_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                       std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)),
                                                ratio*multi_chip->u_write_cycle[data_type_t::WEIGHT]));

                // Update overlapped cycle between off-chip memory and multi-chip processor
                if(num_access_dram == 1) {
                    cycle_chip_dram[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT] +
                                                            u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)) +
                                                            ratio*multi_chip->u_write_cycle[data_type_t::WEIGHT];
                } else {
                    cycle_chip_dram[data_type_t::WEIGHT] += first_stage + second_stage +
                                                           (num_access_dram-2)*other_stage +
                                                           last_before_stage + last_stage;
                }

                // Update data transfer cycle and energy between DRAM and Global buffer.
                transfer_cycle[data_type_t::WEIGHT] += num_access_dram*u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth));
                transfer_energy[data_type_t::WEIGHT] += num_access_dram*u_transfer_energy*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth));

                multi_chip->skip_transfer[data_type_t::WEIGHT] = false;
            }
        }
        // Case 2. COO data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                std::cout << "Current version does not support COO format" << std::endl;
                exit(1);
            }
        }
        // Case 3. CSC data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                // Row bit calculation
                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::CHIPS_Y);
                unsigned row = parameters[parameter_type_t::FILTER_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer[data_type_t::WEIGHT]++;

                // Update off-chip memory access cost
                access_cycle[data_type_t::WEIGHT] += (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                     (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT]
                                                    /(sizeof(data_t)*8/row_bit) + // row index
                                                     parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                    *u_read_cycle[data_type_t::WEIGHT]
                                                    /(sizeof(data_t)*8/row_bit); // Column pointer
                access_energy[data_type_t::WEIGHT] += (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                      (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT]
                                                     /(sizeof(data_t)*8/row_bit) + // Row index
                                                      parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                     *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                     *u_read_energy[data_type_t::WEIGHT]
                                                     /(sizeof(data_t)*8/row_bit); // Column pointer

                // Update on-chip processor access cost
                multi_chip->access_cycle[data_type_t::WEIGHT] += (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                *multi_chip->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                 (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                *multi_chip->u_write_cycle[data_type_t::WEIGHT]
                                                                /(sizeof(data_t)*8/row_bit) + // row index
                                                                 parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                                *multi_chip->u_write_cycle[data_type_t::WEIGHT]
                                                                /(sizeof(data_t)*8/row_bit); // Column pointer
                multi_chip->access_energy[data_type_t::WEIGHT] += (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                 *multi_chip->u_write_energy[data_type_t::WEIGHT] + // Non-zero data
                                                                  (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                 *multi_chip->u_write_energy[data_type_t::WEIGHT]
                                                                 /(sizeof(data_t)*8/row_bit) + // Row index
                                                                  parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                 *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                 *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                                 *multi_chip->u_write_energy[data_type_t::WEIGHT]
                                                                 /(sizeof(data_t)*8/row_bit); // Column pointer

                // Update overlapped cycle between the off-chip memory and on-chip processor
                cycle_chip_dram[data_type_t::WEIGHT] += std::max((multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                        (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle[data_type_t::WEIGHT]
                                                       /(sizeof(data_t)*8/row_bit) + // Row index
                                                        parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                       *u_read_cycle[data_type_t::WEIGHT]
                                                       /(sizeof(data_t)*8/row_bit), // Column pointer
                                                        (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *multi_chip->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                        (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *multi_chip->u_write_cycle[data_type_t::WEIGHT]
                                                       /(sizeof(data_t)*8/row_bit) + // Row index
                                                        parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                       *multi_chip->u_write_cycle[data_type_t::WEIGHT]
                                                       /(sizeof(data_t)*8/row_bit)); // Column pointer

                // Update transfer cost between the off-chip memory and on-chip processor
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle
                                                      *ceil((float)((multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))
                                                      /(float)bitwidth) + // Non-zero data
                                                       u_transfer_cycle
                                                      *ceil((float)((multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*row_bit)
                                                      /(float)bitwidth) + // Column index
                                                       u_transfer_cycle
                                                      *ceil((float)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                      *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                      *(parameters[parameter_type_t::FILTER_WIDTH+1])*row_bit)
                                                      /(float)bitwidth); // Row pointer
                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy
                                                       *ceil((float)((multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))
                                                       /(float)bitwidth) + // Non-zero data
                                                        u_transfer_energy
                                                       *ceil((float)((multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*row_bit)
                                                       /(float)bitwidth) + // Column index
                                                        u_transfer_energy
                                                       *ceil((float)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::FILTER_WIDTH+1])*row_bit)
                                                       /(float)bitwidth); // Row pointer
                multi_chip->skip_transfer[data_type_t::WEIGHT] = false;

            }
        }
        // Case 4. CSR data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                // Column bit calculation
                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::CHIPS_Y);
                unsigned column = parameters[parameter_type_t::FILTER_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer[data_type_t::WEIGHT]++;

                // Update off-chip memory access cost
                access_cycle[data_type_t::WEIGHT] += (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                     (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT]
                                                    /(sizeof(data_t)*8/column_bit) + // Column index
                                                     parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                    *u_read_cycle[data_type_t::WEIGHT]
                                                    /(sizeof(data_t)*8/column_bit); // Row pointer
                access_energy[data_type_t::WEIGHT] += (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                      (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT]
                                                     /(sizeof(data_t)*8/column_bit) + // Column index
                                                      parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                     *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                     *u_read_energy[data_type_t::WEIGHT]
                                                     /(sizeof(data_t)*8/column_bit); // Row pointer

                // Update on-chip processor access cost
                multi_chip->access_cycle[data_type_t::WEIGHT] += (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                *multi_chip->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                 (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                *multi_chip->u_write_cycle[data_type_t::WEIGHT]
                                                                /(sizeof(data_t)*8/column_bit) + // Column index
                                                                 parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                                *multi_chip->u_write_cycle[data_type_t::WEIGHT]
                                                                /(sizeof(data_t)*8/column_bit); // row pointer
                multi_chip->access_energy[data_type_t::WEIGHT] += (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                 *multi_chip->u_write_energy[data_type_t::WEIGHT] + // Non-zero data
                                                                  (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                 *multi_chip->u_write_energy[data_type_t::WEIGHT]
                                                                 /(sizeof(data_t)*8/column_bit) + // Column index
                                                                  parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                 *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                 *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                                 *multi_chip->u_write_energy[data_type_t::WEIGHT]
                                                                 /(sizeof(data_t)*8/column_bit); // row pointer

                // Update overlapped cycle between the off-chip memory and on-chip processor
                cycle_chip_dram[data_type_t::WEIGHT] += std::max((multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                        (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle[data_type_t::WEIGHT]
                                                       /(sizeof(data_t)*8/column_bit) + // Column index
                                                        parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                       *u_read_cycle[data_type_t::WEIGHT]
                                                       /(sizeof(data_t)*8/column_bit), // Row pointer
                                                        (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *multi_chip->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                        (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *multi_chip->u_write_cycle[data_type_t::WEIGHT]
                                                       /(sizeof(data_t)*8/column_bit) + // Column index
                                                        parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                       *multi_chip->u_write_cycle[data_type_t::WEIGHT]/(sizeof(data_t)*8/column_bit)); // row pointer

                // Update transfer cost between the off-chip memory and on-chip processor
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle
                                                      *ceil((float)((multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))
                                                      /(float)bitwidth) + // Non-zero data
                                                       u_transfer_cycle
                                                      *ceil((float)((multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*column_bit)
                                                      /(float)bitwidth) + // Column index
                                                       u_transfer_cycle
                                                      *ceil((float)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                      *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                      *(parameters[parameter_type_t::FILTER_HEIGHT+1])*column_bit)
                                                      /(float)bitwidth); // Row pointer
                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy
                                                       *ceil((float)((multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))
                                                       /(float)bitwidth) + // Non-zero data
                                                        u_transfer_energy
                                                       *ceil((float)((multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*column_bit)
                                                       /(float)bitwidth) + // Column index
                                                        u_transfer_energy
                                                       *ceil((float)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::FILTER_HEIGHT+1])*column_bit)
                                                       /(float)bitwidth); // Row pointer
                multi_chip->skip_transfer[data_type_t::WEIGHT] = false;

            }
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                num_data_transfer[data_type_t::WEIGHT]++;

                // Update off-chip memory access cost
                access_cycle[data_type_t::WEIGHT] += (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                     multi_chip->tile_size[data_type_t::WEIGHT]
                                                    *u_read_cycle[data_type_t::WEIGHT]
                                                    /(sizeof(data_t)*8); // Metadata
                access_energy[data_type_t::WEIGHT] += (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                      multi_chip->tile_size[data_type_t::WEIGHT]
                                                     *u_read_energy[data_type_t::WEIGHT]
                                                     /(sizeof(data_t)*8); // Metadata

                // Update on-chip processor access cost
                multi_chip->access_cycle[data_type_t::WEIGHT] += (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                *multi_chip->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                 multi_chip->tile_size[data_type_t::WEIGHT]
                                                                *multi_chip->u_write_cycle[data_type_t::WEIGHT]
                                                                /(sizeof(data_t)*8); // Metadata
                multi_chip->access_energy[data_type_t::WEIGHT] += (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                 *multi_chip->u_write_energy[data_type_t::WEIGHT] + // Non-zero data
                                                                  multi_chip->tile_size[data_type_t::WEIGHT]
                                                                 *multi_chip->u_write_energy[data_type_t::WEIGHT]
                                                                 /(sizeof(data_t)*8); // Metadata

                // Update overlapped cycle between the off-chip memory and on-chip processor
                cycle_chip_dram[data_type_t::WEIGHT] += std::max((multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                 multi_chip->tile_size[data_type_t::WEIGHT]
                                                                *u_read_cycle[data_type_t::WEIGHT]
                                                                /(sizeof(data_t)*8), // Metadata
                                                                 (multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                *multi_chip->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                 multi_chip->tile_size[data_type_t::WEIGHT]
                                                                *multi_chip->u_write_cycle[data_type_t::WEIGHT]); // Metadata

                // Update transfer cost between the off-chip memory and on-chip processor
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle
                                                      *ceil((float)((multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))
                                                      /(float)bitwidth) + // Non-zero data
                                                       u_transfer_cycle
                                                      *ceil((float)(multi_chip->tile_size[data_type_t::WEIGHT])
                                                      /(float)bitwidth); // Metadata
                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy
                                                       *ceil((float)((multi_chip->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))
                                                       /(float)bitwidth) + // Non-zero data
                                                        u_transfer_energy
                                                       *ceil((float)(multi_chip->tile_size[data_type_t::WEIGHT])
                                                       /(float)bitwidth); // Metadata

                // Update data transfer cycle and energy between DRAM and Global buffer
                multi_chip->skip_transfer[data_type_t::WEIGHT] = false;
                
            }
        }
        else {
            std::cerr << "Undefined data format" << std::endl;
            exit(1);
        }

#else
        if(!skip_transfer[data_type_t::WEIGHT]) {
            num_data_transfer[data_type_t::WEIGHT]++;

            std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_dram(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);
            parameters_dram = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::DRAM);

            uint64_t address_multi_chip = 0, address_dram = 0;
            unsigned num_access_multi_chip = 0, num_access_dram = 0;

            for(unsigned k = 0; k < parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                for(unsigned c = 0; c < parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]; c++) {
                    for(unsigned r = 0; r < parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]; r++) {
                        for(unsigned s = 0; s < parameters_multi_chip[parameter_type_t::FILTER_WIDTH]; s++) {
                            if(address_multi_chip != ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::WEIGHT] + 
                                                                                  k*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] +
                                                                                  c*parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                                  r*parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                  multi_chip->mask_bits[data_type_t::WEIGHT]) << multi_chip->mask_bits[data_type_t::WEIGHT]) {

                                // Update multi-chip processor cost (temporal buffer)
                                multi_chip->access_energy[data_type_t::WEIGHT] += multi_chip->u_write_energy[data_type_t::WEIGHT];
                                if(!multi_chip->double_buffer) {
                                    multi_chip->access_cycle[data_type_t::WEIGHT] += multi_chip->u_write_cycle[data_type_t::WEIGHT];
                                }
                                num_access_multi_chip++;

                                // Update multi-chip processor address
                                address_multi_chip = ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::WEIGHT] + 
                                                                                  k*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                                   *parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                                   *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                                  c*parameters_multi_chip[parameter_type_t::FILTER_HEIGHT] 
                                                                                   *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                                  r*parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                  multi_chip->mask_bits[data_type_t::WEIGHT]) << multi_chip->mask_bits[data_type_t::WEIGHT];
                            }

                            if(address_dram != ((uint64_t)&layer->weight[m_scheduler->weight_offset_dram.front() +
                                                                         k*parameters_dram[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_dram[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_dram[parameter_type_t::FILTER_WIDTH] + 
                                                                         c*parameters_dram[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_dram[parameter_type_t::FILTER_WIDTH] + 
                                                                         r*parameters_dram[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                         mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT]) {
                                // Update off-chip memory cost
                                access_energy[data_type_t::WEIGHT] += u_read_energy[data_type_t::WEIGHT];
                                if(!multi_chip->double_buffer) {
                                    access_cycle[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT];
                                }
                                num_access_dram++;

                                // Update off-chip memory address
                                address_dram = ((uint64_t)&layer->weight[m_scheduler->weight_offset_dram.front() +
                                                                         k*parameters_dram[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_dram[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_dram[parameter_type_t::FILTER_WIDTH] + 
                                                                         c*parameters_dram[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_dram[parameter_type_t::FILTER_WIDTH] + 
                                                                         r*parameters_dram[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                         mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT];
                            }
                        }
                    }
                }
            }

            unsigned ratio = ceil((double)(line_size[data_type_t::WEIGHT])/(double)(multi_chip->line_size[data_type_t::WEIGHT]));

            // At the 1, 2, before last, and last stages
            double first_stage = u_read_cycle[data_type_t::WEIGHT];
            double second_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                             u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
            double last_before_stage = std::max(ratio*multi_chip->u_write_cycle[data_type_t::WEIGHT],
                                                  u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
            double last_stage = ratio*multi_chip->u_write_cycle[data_type_t::WEIGHT];

            // Remainder stages
            double other_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                   std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)),
                                            ratio*multi_chip->u_write_cycle[data_type_t::WEIGHT]));

            // Update overlapped cycle between off-chip memory and multi-chip processor
            if(num_access_dram == 1) {
                cycle_chip_dram[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT] +
                                                        u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)) +
                                                        ratio*multi_chip->u_write_cycle[data_type_t::WEIGHT];
            } else {
                cycle_chip_dram[data_type_t::WEIGHT] += first_stage + second_stage +
                                                       (num_access_dram-2)*other_stage +
                                                       last_before_stage + last_stage;
            }

            // Update data transfer cycle and energy between DRAM and Global buffer.
            transfer_cycle[data_type_t::WEIGHT] += num_access_dram*u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth));
            transfer_energy[data_type_t::WEIGHT] += num_access_dram*u_transfer_energy*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth));

            multi_chip->skip_transfer[data_type_t::WEIGHT] = false;

        }
#endif
        weight_index++;
        multi_chip->exist_data[data_type_t::WEIGHT] = true, multi_chip->request_to_dram[data_type_t::WEIGHT] = false;
        if(multi_chip->tile_size[data_type_t::WEIGHT] == tile_size[data_type_t::WEIGHT]) {skip_transfer[data_type_t::WEIGHT] =true;}
    }
    if(multi_chip->request_to_dram[data_type_t::OUTPUT]) {
        // Check whether the output data should be read 
        // or does not have to be transferred but just be initialized.
        if(m_scheduler->output_read_dram[m_scheduler->output_offset_dram.front()]) {
            if(!skip_transfer[data_type_t::OUTPUT]) {

#ifdef DRAMSIM3
                //send_request((data_t*)layer->output_data, m_scheduler->mapping_table, m_scheduler->output_offset_dram.front(), data_type_t::OUTPUT, action_type_t::LOAD);
                send_request((data_t*)data, m_scheduler->mapping_table, offsets[data_type_t::OUTPUT] + m_scheduler->output_offset_dram.front(), data_type_t::OUTPUT, action_type_t::LOAD);
#endif

#ifdef FUNCTIONAL
                // Weight transfer
                m_scheduler->transfer_data(multi_chip->data, (data_t*)layer->output_data, multi_chip->offsets[data_type_t::OUTPUT], m_scheduler->output_offset_dram.front(), 
                                           component_type_t::CHIPS_Y, component_type_t::DRAM, 
                                           data_type_t::OUTPUT, multi_chip->get_stationary_type(), action_type_t::LOAD);
                // Update for NPUsim ver2
                //m_scheduler->transfer_data_ver2(multi_chip->data, (data_t*)layer->output_data,
                //                                component_type_t::CHIPS_Y, component_type_t::DRAM,
                //                                data_type_t::OUTPUT, multi_chip->get_stationary_type(), 
                //                                action_type_t::LOAD, true);

#endif
                /* Stats */
                num_data_transfer[data_type_t::OUTPUT]++;
                std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_dram(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);
                parameters_dram = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::DRAM);

                uint64_t address_multi_chip = 0, address_dram = 0;
                unsigned num_access_multi_chip = 0, num_access_dram = 0;

                for(unsigned b = 0; b < parameters_multi_chip[parameter_type_t::BATCH_SIZE]; b++) {
                    for(unsigned k = 0; k < parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                        for(unsigned p = 0; p < parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                            for(unsigned q = 0; q < parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH]; q++) {
                                if(address_multi_chip != ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::OUTPUT] + 
                                                                                      b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                                       *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                       *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] +
                                                                                      k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                       *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                      p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                      multi_chip->mask_bits[data_type_t::OUTPUT]) << multi_chip->mask_bits[data_type_t::OUTPUT]) {

                                    // Update multi-chip processor cost
                                    multi_chip->access_energy[data_type_t::OUTPUT] += multi_chip->u_write_energy[data_type_t::OUTPUT];
                                    if(!multi_chip->double_buffer) {
                                        multi_chip->access_cycle[data_type_t::OUTPUT] += multi_chip->u_write_cycle[data_type_t::OUTPUT];
                                    }
                                    num_access_multi_chip++;

                                    // Update multi-chip processor address
                                    address_multi_chip = ((uint64_t)&multi_chip->data[multi_chip->offsets[data_type_t::OUTPUT] + 
                                                                                      b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                                       *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                       *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] +
                                                                                      k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                                       *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                      p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                      multi_chip->mask_bits[data_type_t::OUTPUT]) << multi_chip->mask_bits[data_type_t::OUTPUT];
                                }

                                // Check off-chip memory address
                                if(address_dram != ((uint64_t)&layer->output_data[m_scheduler->output_offset_dram.front() + 
                                                                                  b*parameters_dram[parameter_type_t::OUTPUT_CHANNEL] 
                                                                                   *parameters_dram[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_dram[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  k*parameters_dram[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_dram[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  p*parameters_dram[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                  mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT]) {
                                    // Update off-chip memory cost
                                    access_energy[data_type_t::OUTPUT] += u_read_energy[data_type_t::OUTPUT];
                                    if(!multi_chip->double_buffer) {
                                        access_cycle[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT];
                                    }
                                    num_access_dram++;

                                    // Update off-chip memory address
                                    address_dram = ((uint64_t)&layer->output_data[m_scheduler->output_offset_dram.front() + 
                                                                                  b*parameters_dram[parameter_type_t::OUTPUT_CHANNEL] 
                                                                                   *parameters_dram[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_dram[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  k*parameters_dram[parameter_type_t::OUTPUT_HEIGHT]
                                                                                   *parameters_dram[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                  p*parameters_dram[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                  mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT];
                                }
                            }
                        }
                    }
                }

                unsigned ratio = ceil((double)(line_size[data_type_t::OUTPUT])/(double)(multi_chip->line_size[data_type_t::OUTPUT]));
                                                                                                                            
                // At the 1, 2, before last, and last stages
                unsigned first_stage = u_read_cycle[data_type_t::OUTPUT];
                unsigned second_stage = std::max(u_read_cycle[data_type_t::OUTPUT],
                                                 u_transfer_cycle*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth)));
                unsigned last_before_stage = std::max(ratio*multi_chip->u_write_cycle[data_type_t::OUTPUT],
                                                      u_transfer_cycle*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth)));
                unsigned last_stage = ratio*multi_chip->u_write_cycle[data_type_t::OUTPUT];

                // Remainder stages
                unsigned other_stage = std::max(u_read_cycle[data_type_t::OUTPUT],
                                       std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth)),
                                                ratio*multi_chip->u_write_cycle[data_type_t::OUTPUT]));
                                                                                                                            
                // Update overlapped cycle between off-chip memory and multi-chip processor                                 
                if(num_access_dram == 1) {
                    cycle_chip_dram[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT] +
                                                            u_transfer_cycle*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth)) +
                                                            ratio*multi_chip->u_write_cycle[data_type_t::OUTPUT];
                } else {
                    cycle_chip_dram[data_type_t::OUTPUT] += first_stage + second_stage +
                                                            (num_access_dram-2)*other_stage +
                                                            last_before_stage + last_stage;
                }
                                                                                                                            
                // Update data transfer cycle and energy between DRAM and Global buffer.                                    
                transfer_cycle[data_type_t::OUTPUT] += num_access_dram*u_transfer_cycle*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth));
                transfer_energy[data_type_t::OUTPUT] += num_access_dram*u_transfer_energy*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth));

                multi_chip->skip_transfer[data_type_t::OUTPUT] = false;
                multi_chip->equal_output_tile = false;
                transfer_output = true;
            }
        }
        else {
            transfer_output = false;
            m_scheduler->output_read_dram[m_scheduler->output_offset_dram.front()] = true;
            for(auto it = m_scheduler->output_read_global_buffer.begin(); it != m_scheduler->output_read_global_buffer.end(); ++it) {
                it->second = false;
            }
            for(auto it = m_scheduler->output_read_pe.begin(); it != m_scheduler->output_read_pe.end(); ++it) {
                it->second = false;
            }
        }
        move_front(&m_scheduler->output_offset_dram);
        output_index++;
        multi_chip->exist_data[data_type_t::OUTPUT] = true, multi_chip->request_to_dram[data_type_t::OUTPUT] = false;
    }
    multi_chip->fill_data();

    if(multi_chip->get_stationary_type() == stationary_type_t::INPUT_STATIONARY) {
        if(weight_index == m_scheduler->offset_size_dram[data_type_t::WEIGHT].front() &&
           output_index == m_scheduler->offset_size_dram[data_type_t::OUTPUT].front() &&
           input_index < m_scheduler->offset_size_dram[data_type_t::INPUT].front()){

            weight_index = 0, output_index = 0;

            move_front(&m_scheduler->offset_size_dram[data_type_t::WEIGHT]);
            move_front(&m_scheduler->offset_size_dram[data_type_t::OUTPUT]);
        }
    }
    else if(multi_chip->get_stationary_type() == stationary_type_t::WEIGHT_STATIONARY) {
        if(input_index == m_scheduler->offset_size_dram[data_type_t::INPUT].front() &&
           output_index == m_scheduler->offset_size_dram[data_type_t::OUTPUT].front() &&
           weight_index < m_scheduler->offset_size_dram[data_type_t::WEIGHT].front()) {

            input_index = 0, output_index = 0;

            move_front(&m_scheduler->offset_size_dram[data_type_t::INPUT]);
            move_front(&m_scheduler->offset_size_dram[data_type_t::OUTPUT]);
        }
    }
    else if(multi_chip->get_stationary_type() == stationary_type_t::OUTPUT_STATIONARY) {
        if(input_index == m_scheduler->offset_size_dram[data_type_t::INPUT].front() &&
           weight_index == m_scheduler->offset_size_dram[data_type_t::WEIGHT].front() &&
           output_index < m_scheduler->offset_size_dram[data_type_t::OUTPUT].front()) {

            input_index = 0, weight_index = 0;

            move_front(&m_scheduler->offset_size_dram[data_type_t::INPUT]);
            move_front(&m_scheduler->offset_size_dram[data_type_t::WEIGHT]);
        }
    }
    // Check if all data of three data types are transferred from DRAM to Global buffer.
    if(input_index == m_scheduler->offset_size_dram[data_type_t::INPUT].front() && weight_index == m_scheduler->offset_size_dram[data_type_t::WEIGHT].front() && output_index == m_scheduler->offset_size_dram[data_type_t::OUTPUT].front()) {
        done = true;
        input_index = 0, weight_index = 0, output_index = 0;
    }
}

void dram_t::print_specification() {
    std::cout << "===== DRAM =====" << std::endl;
    std::cout << std::endl;

}

// Reset stats of the off-chip memory.
void dram_t::reset() {
    done = false;

    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);

    num_request.assign(data_type_t::NUM_DATA_TYPES, 0);
    num_data_transfer.assign(data_type_t::NUM_DATA_TYPES, 0);

    access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    cycle_chip_dram.assign(data_type_t::NUM_DATA_TYPES, 0.0);
}

#ifdef DRAMSIM3
void dram_t::print_result() {
    memory->print_stats();
}

void dram_t::send_request(data_t *m_data, mapping_table_t *m_mapping_table, unsigned m_offset, data_type_t m_data_type, action_type_t m_action_type) {
    memory->send_request(m_data, m_mapping_table, m_offset, m_data_type, m_action_type);
}
#endif
