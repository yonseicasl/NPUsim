#include "adder_tree.h"

adder_tree_t::adder_tree_t(section_config_t m_section_config) :
    pe_array_t(m_section_config) {

    init(m_section_config);

}

adder_tree_t::~adder_tree_t() {

	delete [] input_data;
	delete [] weight;
	delete [] output_data;

	for(unsigned i = 0; i < pes.size(); i++) {
		delete pes[i];
	}
}

// Initialize the PE array.
void adder_tree_t::init(section_config_t m_section_config) {

    /* Initialize adder tree specifications */

    // Initialize the number of PE array.
    m_section_config.get_setting("height", &height);
    m_section_config.get_setting("width", &width);
    num_pes = height * width;

    // Initialize frequency and bandwidth of adder tree
    m_section_config.get_setting("frequency", &frequency);
    m_section_config.get_setting("bandwidth", &bandwidth);
    bitwidth = 8*bandwidth/frequency;
    m_section_config.get_setting("bitwidth", &bitwidth);

    // Initialize line size and mask bits of temporal buffer in PE array.
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

    // Initialize stationary type at the PE array
    std::string array_stationary_str;
    if(m_section_config.get_setting("pe_stationary", &array_stationary_str)) {
        stationary_type = (stationary_type_t)get_type(stationary_type_str, array_stationary_str);
    }

    // Initialize the order of parameters
    if(stationary_type == stationary_type_t::INPUT_STATIONARY) {
        array_parameter_order = "BCPQKRS";
    } else if(stationary_type == stationary_type_t::WEIGHT_STATIONARY) {
        array_parameter_order = "KCRSBPQ";
    } else if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
        array_parameter_order = "BKPQKRS";
    }
    m_section_config.get_setting("pe_array_parameter_order", &array_parameter_order);

    m_section_config.get_setting("exist_temporal_buffer", &exist_temporal_buffer);
    std::string memory_type_str;
    m_section_config.get_setting("memory_type", &memory_type_str);
    if(memory_type_str == "shared") {
        memory_type = memory_type_t::SEPARATE;

        m_section_config.get_setting("input_buffer", &input_size);
        m_section_config.get_setting("weight_buffer", &weight_size);
        m_section_config.get_setting("output_buffer", &output_size);

        input_size *= num_pes, weight_size *= num_pes, output_size *= num_pes;

        unsigned num_input = input_size*num_pes/sizeof(data_t);
        unsigned num_weight = weight_size*num_pes/sizeof(data_t);
        unsigned num_output = output_size*num_pes/sizeof(data_t);

        input_data  = new data_t[num_input]();
        weight      = new data_t[num_weight]();
        output_data = new data_t[num_output]();

        memset(input_data, 1.0, num_input*sizeof(data_t));
        memset(weight, 1.0, num_weight*sizeof(data_t));
        memset(output_data, 1.0, num_output*sizeof(data_t));
    }
    else if(memory_type_str == "separate") {
        memory_type = memory_type_t::SHARED;

        m_section_config.get_setting("input_buffer", &input_size);
        m_section_config.get_setting("weight_buffer", &weight_size);
        m_section_config.get_setting("output_buffer", &output_size);

        input_size *= num_pes, weight_size *= num_pes, output_size *= num_pes;

        unsigned num_input = input_size*num_pes/sizeof(data_t);
        unsigned num_weight = weight_size*num_pes/sizeof(data_t);
        unsigned num_output = output_size*num_pes/sizeof(data_t);

        input_data  = new data_t[num_input]();
        weight      = new data_t[num_weight]();
        output_data = new data_t[num_output]();

        memset(input_data, 1.0, num_input*sizeof(data_t));
        memset(weight, 1.0, num_weight*sizeof(data_t));
        memset(output_data, 1.0, num_output*sizeof(data_t));
    }

    // Initialize the NoC type.
    std::string noc_str;
    if(m_section_config.get_setting("noc", &noc_str)) {
        noc_type = (noc_type_t)get_type(noc_type_str, noc_str);
    }

    // Define stationary of each PE (between MAC and Local buffer inside the PE).
	pe_t* pe;
    std::string pe_stationary_type_str;
    m_section_config.get_setting("mac_stationary", &pe_stationary_type_str);
    if(pe_stationary_type_str == "input_stationary") {
        for(unsigned i = 0; i < num_pes; i++) {
            pe = new input_stationary_t(m_section_config);
            pe->index = i;
            pes.emplace_back(pe);
        }
    }
    else if(pe_stationary_type_str == "weight_stationary") {
        for(unsigned i = 0; i < num_pes; i++) {
            pe = new weight_stationary_t(m_section_config);
            pe->index = i;
            pes.emplace_back(pe);
        }
    }
    else if(pe_stationary_type_str == "output_stationary") {
        for(unsigned i = 0; i < num_pes; i++) {
            pe = new output_stationary_t(m_section_config);
            pe->index = i;
            pes.emplace_back(pe);
        }
    }
	else {
		std::cerr << "Error: Wrong stationary type name : " << pe_stationary_type_str << std::endl;
		exit(1);
	}

    skip_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);

    // Initialize the tile size of PE array.
    tile_size.reserve(data_type_t::NUM_DATA_TYPES);
    tile_size.assign(data_type_t::NUM_DATA_TYPES, 1);

    offsets.reserve(data_type_t::NUM_DATA_TYPES);
    offsets.assign(data_type_t::NUM_DATA_TYPES, 0);

    /* Initialize signals of adder tree */

    exist_data.reserve(data_type_t::NUM_DATA_TYPES);
    exist_data.assign(data_type_t::NUM_DATA_TYPES, false);

    request_to_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    request_to_global_buffer.assign(data_type_t::NUM_DATA_TYPES, false);

    /* Initialize unit cost of adder tree */

    // Initialize the stats between PE array and Global buffer.
    u_read_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("pe_array_read_cycle", &u_read_cycle);

    u_read_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("pe_array_read_energy", &u_read_energy);

    u_write_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("pe_array_write_cycle", &u_write_cycle);

    u_write_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("pe_array_write_energy", &u_write_energy);

    // Initialize the NoC cycle of PE array.
    m_section_config.get_setting("noc_cycle", &noc_cycle);
    m_section_config.get_setting("noc_energy", &noc_energy);

    /* Initialize stats of adder tree */

    num_request.reserve(data_type_t::NUM_DATA_TYPES);
    num_request.assign(data_type_t::NUM_DATA_TYPES, 0);

    num_data_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer.assign(data_type_t::NUM_DATA_TYPES, 0);

    access_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    access_energy.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    transfer_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    transfer_energy.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    cycle_temporal_pe.reserve(data_type_t::NUM_DATA_TYPES);
    cycle_temporal_pe.assign(data_type_t::NUM_DATA_TYPES, 0.0);

}

void adder_tree_t::update_tile_size(scheduler_t *m_scheduler) {
    num_active_pe_x = m_scheduler->num_active_pe_x;
    num_active_pe_y = m_scheduler->num_active_pe_y;
    utilization = (double)(get_number_of_active_pes())/(double)(get_number_of_pes());
    
    // Update PE array tile size.
    tile_size = m_scheduler->tile_size[component_type_t::PE_Y];

    // Update PEs' tile size.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        pes[i]->update_tile_size(m_scheduler);
        pes[i]->check_tile_size();
    }
}

void adder_tree_t::data_transfer(scheduler_t *m_scheduler) {

    bool request_to_pe_array_input = false;
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::INPUT]) {
            request_to_pe_array_input = true;
            break;
        }
    }
    if(request_to_pe_array_input) {
#ifdef FUNCTIONAL
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::INPUT]) {
                num_data_transfer[data_type_t::INPUT]++;

                std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
                parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

                unsigned num_access_pe = 0, num_access_pe_array = 0;
                for(unsigned i = 0; i < get_number_of_active_pes(); i++) { 
                    uint64_t address_pe = 0, address_pe_array = 0;
                    m_scheduler->transfer_data(pes[i]->input_data_lb, input_data, 0, m_scheduler->input_offset_pe_array[i%m_scheduler->input_offset_pe_array.size()], 
                                               component_type_t::PE, component_type_t::PE_Y, 
                                               data_type_t::INPUT, pes[i]->get_local_buffer_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
                    //m_scheduler->transfer_data_network_on_chip(pes[i]->input_data_lb, input_data, 
                    //                                           component_type_t::PE, component_type_t::PE_Y,
                    //                                           data_type_t::INPUT, pes[i]->get_local_buffer_stationary_type(), 
                    //                                           action_type_t::LOAD, last_component);

                    for(unsigned b = 0; b < parameters_pe[parameter_type_t::BATCH_SIZE]; b++) {
                        for(unsigned c = 0; c < parameters_pe[parameter_type_t::INPUT_CHANNEL]; c++) {
                            for(unsigned h = 0; h < parameters_pe[parameter_type_t::INPUT_HEIGHT]; h++) {
                                for(unsigned w = 0; w < parameters_pe[parameter_type_t::INPUT_WIDTH]; w++) {

                                    if(address_pe != ((uint64_t)&pes[i]->input_data_lb[b*parameters_pe[parameter_type_t::INPUT_CHANNEL]
                                                                                        *parameters_pe[parameter_type_t::INPUT_HEIGHT]
                                                                                        *parameters_pe[parameter_type_t::INPUT_WIDTH] + 
                                                                                       c*parameters_pe[parameter_type_t::INPUT_HEIGHT]
                                                                                        *parameters_pe[parameter_type_t::INPUT_WIDTH] + 
                                                                                       h*parameters_pe[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                       pes[i]->mask_bits_lb[data_type_t::INPUT]) << pes[i]->mask_bits_lb[data_type_t::INPUT]) {

                                        // Update local buffer cost in PEs
                                        pes[i]->access_cycle_lb[data_type_t::INPUT] += pes[i]->u_write_cycle_lb[data_type_t::INPUT];
                                        pes[i]->access_energy_lb[data_type_t::INPUT] += pes[i]->u_write_energy_lb[data_type_t::INPUT];

                                        num_access_pe++;

                                        // Update local buffer address of PEs
                                        address_pe = ((uint64_t)&pes[i]->input_data_lb[b*parameters_pe[parameter_type_t::INPUT_CHANNEL]
                                                                                        *parameters_pe[parameter_type_t::INPUT_HEIGHT]
                                                                                        *parameters_pe[parameter_type_t::INPUT_WIDTH] + 
                                                                                       c*parameters_pe[parameter_type_t::INPUT_HEIGHT]
                                                                                        *parameters_pe[parameter_type_t::INPUT_WIDTH] + 
                                                                                       h*parameters_pe[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                       pes[i]->mask_bits_lb[data_type_t::INPUT]) << pes[i]->mask_bits_lb[data_type_t::INPUT];
                                    }

                                    // Check address of temporal buffer in PE array
                                    if(!m_scheduler->read_tile_granular_pe_input[i%m_scheduler->input_offset_pe_array.size()] && 
                                       address_pe_array != ((uint64_t)&input_data[m_scheduler->input_offset_pe_array[i%m_scheduler->input_offset_pe_array.size()] + 
                                                                                  b*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                   *parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                   *parameters_pe_array[parameter_type_t::INPUT_WIDTH] + 
                                                                                  c*parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                   *parameters_pe_array[parameter_type_t::INPUT_WIDTH] + 
                                                                                  h*parameters_pe_array[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                  mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT]) { 

                                        // Update temporal buffer cost in PE array
                                        access_cycle[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT];
                                        access_energy[data_type_t::INPUT] += u_read_energy[data_type_t::INPUT];
                                        num_access_pe_array++;

                                        // Update address of temporal buffer in PE array
                                        m_scheduler->read_tile_granular_pe_input[i%m_scheduler->input_offset_pe_array.size()] = true;
                                        address_pe_array = ((uint64_t)&input_data[m_scheduler->input_offset_pe_array[i%m_scheduler->input_offset_pe_array.size()] + 
                                                                                  b*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                   *parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                   *parameters_pe_array[parameter_type_t::INPUT_WIDTH] + 
                                                                                  c*parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                   *parameters_pe_array[parameter_type_t::INPUT_WIDTH] + 
                                                                                  h*parameters_pe_array[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                  mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT];
                                    }
                                }
                            }
                        }
                    }
                    pes[i]->skip_transfer[data_type_t::INPUT] = false;
                }

                // Update transfer cycle and energy
                transfer_energy[data_type_t::INPUT] += num_access_pe*noc_energy*ceil((double)(pes[0]->line_size_lb[data_type_t::INPUT])/(double)(bitwidth));
                transfer_cycle[data_type_t::INPUT] += num_access_pe_array*noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::INPUT])/(double)(bitwidth));

                // At the 1, 2, before last, and last stages
                double first_stage = u_read_cycle[data_type_t::INPUT];
                double second_stage = std::max(u_read_cycle[data_type_t::INPUT],
                                               noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::INPUT])/(double)(bitwidth)));

                double last_before_stage = std::max(pes[0]->u_write_cycle_lb[data_type_t::INPUT], 
                                                    noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::INPUT])/(double)(bitwidth)));
                double last_stage = pes[0]->u_write_cycle_lb[data_type_t::INPUT];

                double other_stage = std::max(u_read_cycle[data_type_t::INPUT],
                                     std::max(pes[0]->u_write_cycle_lb[data_type_t::INPUT],
                                              noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::INPUT])/(double)(bitwidth))));

                // if temporal buffer is exist at the PE array
                if(exist_temporal_buffer) {
                    // Update overlapped cycle of Scattering
                    if(num_access_pe_array == 1) {
                        cycle_temporal_pe[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT] + 
                                                                 noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::INPUT])/(double)(bitwidth)) + 
                                                                 pes[0]->u_write_cycle_lb[data_type_t::INPUT];
                   } else {
                       cycle_temporal_pe[data_type_t::INPUT] += first_stage + second_stage +
                                                                (num_access_pe_array-2)*other_stage + 
                                                                last_before_stage + last_stage; 
                   }
                }
            }
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::INPUT]) {

                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::PE);
                unsigned row = parameters[parameter_type_t::INPUT_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer[data_type_t::INPUT]++;
                
                for(unsigned i = 0; i < get_number_of_active_pes(); i++) {

                    m_scheduler->transfer_data(pes[i]->input_data_lb, input_data, 0, m_scheduler->input_offset_pe_array[i%m_scheduler->input_offset_pe_array.size()], 
                                               component_type_t::PE, component_type_t::PE_Y, 
                                               data_type_t::INPUT, pes[i]->get_local_buffer_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
                    //m_scheduler->transfer_data_network_on_chip(pes[i]->input_data_lb, input_data, 
                    //                                           component_type_t::PE, component_type_t::PE_Y,
                    //                                           data_type_t::INPUT, pes[i]->get_local_buffer_stationary_type(), 
                    //                                           action_type_t::LOAD, last_component);

                    if(!m_scheduler->read_tile_granular_pe_input[i%m_scheduler->input_offset_pe_array.size()]) {
                        // Update access cost to PE array
                        access_cycle[data_type_t::INPUT] += (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                           *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                            (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                           *u_read_cycle[data_type_t::INPUT]
                                                           /(sizeof(data_t)*8/row_bit) + // Row index
                                                            parameters[parameter_type_t::BATCH_SIZE]
                                                           *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                           *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                           *u_read_cycle[data_type_t::INPUT]
                                                           /(sizeof(data_t)*8/row_bit); // Column pointer

                        access_energy[data_type_t::INPUT] += (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_energy[data_type_t::INPUT] + // Non-zero data
                                                             (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_energy[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/row_bit) + // Row index
                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                            *u_read_energy[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/row_bit); // Column pointer

                        m_scheduler->read_tile_granular_pe_input[i%m_scheduler->input_offset_pe_array.size()] = true;

                        transfer_cycle[data_type_t::INPUT] += noc_cycle
                                                             *ceil((double)(pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *8*sizeof(data_t)
                                                             /(double)bitwidth) + // Non-zero data
                                                              noc_cycle
                                                             *ceil((double)((pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)
                                                             /(double)(bitwidth)) + // Row index
                                                              noc_cycle
                                                             *ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                           *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                           *(parameters[parameter_type_t::INPUT_WIDTH]+1)*row_bit)
                                                             /(double)(bitwidth)); // Column pointer
                    }

                    // Update local buffers access cost
                    pes[i]->access_cycle_lb[data_type_t::INPUT] += (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                  *pes[i]->u_write_cycle_lb[data_type_t::INPUT] + // Non-zero data
                                                                   (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                  *pes[i]->u_write_cycle_lb[data_type_t::INPUT]
                                                                  /(sizeof(data_t)*8/row_bit) + // Row index
                                                                   parameters[parameter_type_t::BATCH_SIZE]
                                                                  *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                  *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                                  *pes[i]->u_write_cycle_lb[data_type_t::INPUT]
                                                                  /(sizeof(data_t)*8/row_bit); // Column pointer
                    pes[i]->access_energy_lb[data_type_t::INPUT] += (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                   *pes[i]->u_write_energy_lb[data_type_t::INPUT] + // Non-zero data
                                                                    (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                   *pes[i]->u_write_energy_lb[data_type_t::INPUT]
                                                                   /(sizeof(data_t)*8/row_bit) + // Row index
                                                                    parameters[parameter_type_t::BATCH_SIZE]
                                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                   *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                                   *pes[i]->u_write_energy_lb[data_type_t::INPUT]
                                                                   /(sizeof(data_t)*8/row_bit); // Column pointer

                    // Update transfer energy between PE array and local buffers in PEs.
                    transfer_energy[data_type_t::INPUT] += noc_energy
                                                         *ceil((double)(pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                         *8*sizeof(data_t)
                                                         /(double)bitwidth) + // Non-zero data
                                                          noc_energy
                                                         *ceil((double)((pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)
                                                         /(double)(bitwidth)) + // Row index
                                                          noc_energy
                                                         *ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                       *(parameters[parameter_type_t::INPUT_WIDTH]+1)*row_bit)
                                                         /(double)bitwidth); // Column pointer

                    unsigned data_size = (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                        *sizeof(data_t) + // Non-zero data
                                         (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                        *row_bit/8 + // Row index
                                         parameters[parameter_type_t::BATCH_SIZE]
                                        *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                        *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                        *row_bit/8; // Column pointer

                    pes[i]->utilization_local_buffer[data_type_t::INPUT] = (double)data_size/(double)pes[i]->input_size;

                    pes[i]->skip_transfer[data_type_t::INPUT] = false;
                }
            }
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::INPUT]) {

                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::PE);
                unsigned column = parameters[parameter_type_t::INPUT_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer[data_type_t::INPUT]++;
                
                for(unsigned i = 0; i < get_number_of_active_pes(); i++) {

                    m_scheduler->transfer_data(pes[i]->input_data_lb, input_data, 0, m_scheduler->input_offset_pe_array[i%m_scheduler->input_offset_pe_array.size()], 
                                               component_type_t::PE, component_type_t::PE_Y, 
                                               data_type_t::INPUT, pes[i]->get_local_buffer_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
                    //m_scheduler->transfer_data_network_on_chip(pes[i]->input_data_lb, input_data, 
                    //                                           component_type_t::PE, component_type_t::PE_Y,
                    //                                           data_type_t::INPUT, pes[i]->get_local_buffer_stationary_type(), 
                    //                                           action_type_t::LOAD, last_component);

                    if(!m_scheduler->read_tile_granular_pe_input[i%m_scheduler->input_offset_pe_array.size()]) {
                        // Update access cost to PE array
                        access_cycle[data_type_t::INPUT] += (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                           *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                            (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                           *u_read_cycle[data_type_t::INPUT]
                                                           /(sizeof(data_t)*8/column_bit) + // Column index
                                                            parameters[parameter_type_t::BATCH_SIZE]
                                                           *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                           *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                           *u_read_cycle[data_type_t::INPUT]
                                                           /(sizeof(data_t)*8/column_bit); // Row pointer

                        access_energy[data_type_t::INPUT] += (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_energy[data_type_t::INPUT] + // Non-zero data
                                                             (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_energy[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/column_bit) + // Column index
                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                            *u_read_energy[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/column_bit); // Row pointer

                        m_scheduler->read_tile_granular_pe_input[i%m_scheduler->input_offset_pe_array.size()] = true;

                        transfer_cycle[data_type_t::INPUT] += noc_cycle
                                                             *ceil((double)(pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *8*sizeof(data_t)
                                                             /(double)bitwidth) + // Non-zero data
                                                              noc_cycle
                                                             *ceil((double)((pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)
                                                             /(double)(bitwidth)) + // Column index
                                                              noc_cycle
                                                             *ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                           *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                           *(parameters[parameter_type_t::INPUT_HEIGHT]+1)*column_bit)
                                                             /(double)(bitwidth)); // Row pointer
                    }

                    // Update local buffers access cost
                    pes[i]->access_cycle_lb[data_type_t::INPUT] += (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                  *pes[i]->u_write_cycle_lb[data_type_t::INPUT] + // Non-zero data
                                                                   (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                  *pes[i]->u_write_cycle_lb[data_type_t::INPUT]
                                                                  /(sizeof(data_t)*8/column_bit) + // Column index
                                                                   parameters[parameter_type_t::BATCH_SIZE]
                                                                  *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                  *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                                  *pes[i]->u_write_cycle_lb[data_type_t::INPUT]
                                                                  /(sizeof(data_t)*8/column_bit); // Row pointer
                    pes[i]->access_energy_lb[data_type_t::INPUT] += (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                   *pes[i]->u_write_energy_lb[data_type_t::INPUT] + // Non-zero data
                                                                    (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                   *pes[i]->u_write_energy_lb[data_type_t::INPUT]
                                                                   /(sizeof(data_t)*8/column_bit) + // Column index
                                                                    parameters[parameter_type_t::BATCH_SIZE]
                                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                   *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                                   *pes[i]->u_write_energy_lb[data_type_t::INPUT]
                                                                   /(sizeof(data_t)*8/column_bit); // Row pointer

                    // Update transfer energy between PE array and local buffers in PEs.
                    transfer_energy[data_type_t::INPUT] += noc_energy
                                                         *ceil((double)(pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                         *8*sizeof(data_t)
                                                         /(double)bitwidth) + // Non-zero data
                                                          noc_energy
                                                         *ceil((double)((pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)
                                                         /(double)(bitwidth)) + // Column index
                                                          noc_energy
                                                         *ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                       *(parameters[parameter_type_t::INPUT_WIDTH]+1)*column_bit)
                                                         /(double)bitwidth); // row pointer

                    unsigned data_size = (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                        *sizeof(data_t) + // Non-zero data
                                         (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                        *column_bit/8 + // Column index
                                         parameters[parameter_type_t::BATCH_SIZE]
                                        *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                        *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                        *column_bit/8; // Row pointer

                    pes[i]->utilization_local_buffer[data_type_t::INPUT] = (double)data_size/(double)pes[i]->input_size;

                    pes[i]->skip_transfer[data_type_t::INPUT] = false;
                }
            }
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::INPUT]) {

                num_data_transfer[data_type_t::INPUT]++;
                for(unsigned i = 0; i < get_number_of_active_pes(); i++) {

                    m_scheduler->transfer_data(pes[i]->input_data_lb, input_data, 0, m_scheduler->input_offset_pe_array[i%m_scheduler->input_offset_pe_array.size()], 
                                               component_type_t::PE, component_type_t::PE_Y, 
                                               data_type_t::INPUT, pes[i]->get_local_buffer_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
                    //m_scheduler->transfer_data_network_on_chip(pes[i]->input_data_lb, input_data, 
                    //                                           component_type_t::PE, component_type_t::PE_Y,
                    //                                           data_type_t::INPUT, pes[i]->get_local_buffer_stationary_type(), 
                    //                                           action_type_t::LOAD, last_component);

                    if(!m_scheduler->read_tile_granular_pe_input[i%m_scheduler->input_offset_pe_array.size()]) {
                        // Update access cost of PE array 
                        access_cycle[data_type_t::INPUT] += (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                           *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                            pes[i]->tile_size_lb[data_type_t::INPUT]
                                                           *u_read_cycle[data_type_t::INPUT]
                                                           /(sizeof(data_t)*8); // Metadata
                    
                        access_energy[data_type_t::INPUT] += (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_energy[data_type_t::INPUT] + // Non-zero data
                                                             pes[i]->tile_size_lb[data_type_t::INPUT]
                                                            *u_read_energy[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8); // Metadata
                        // Update transfer cycle between the chip-level processor and the global buffer
                        transfer_cycle[data_type_t::INPUT] += noc_cycle
                                                             *ceil((double)((pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *8*sizeof(data_t))
                                                             /(double)bitwidth) + // Non-zero data
                                                              noc_cycle 
                                                             *ceil((double)(pes[i]->tile_size_lb[data_type_t::INPUT])/(double)bitwidth); // Metadata
                    }

                    // Update local buffers access cost
                    pes[i]->access_cycle_lb[data_type_t::INPUT] += (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                  *pes[i]->u_write_cycle_lb[data_type_t::INPUT] + // Non-zero data
                                                                   pes[i]->tile_size_lb[data_type_t::INPUT]
                                                                  *pes[i]->u_write_cycle_lb[data_type_t::INPUT]
                                                                  /(sizeof(data_t)*8); // Metadata
                   
                    pes[i]->access_energy_lb[data_type_t::INPUT] += (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                   *pes[i]->u_write_energy_lb[data_type_t::INPUT] + // Non-zero data
                                                                    pes[i]->tile_size_lb[data_type_t::INPUT]
                                                                   *pes[i]->u_write_energy_lb[data_type_t::INPUT]
                                                                   /(sizeof(data_t)*8); // Metadata

                    transfer_energy[data_type_t::INPUT] += noc_energy
                                                          *ceil((double)((pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                          *8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                           noc_energy
                                                          *ceil((double)(pes[i]->tile_size_lb[data_type_t::INPUT])/(double)bitwidth); // Metadata

                    unsigned data_size = (pes[i]->tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                        *sizeof(data_t) + // Non-zero data
                                         pes[i]->tile_size_lb[data_type_t::INPUT]/8; // Metadata

                    pes[i]->utilization_local_buffer[data_type_t::INPUT] = (double)data_size/(double)pes[i]->input_size;

                    pes[i]->skip_transfer[data_type_t::INPUT] = false;
                }    
            }
        }
        else {
            std::cerr << "Undefined compression type" << std::endl;
            exit(1);
        }
#else
        if(!skip_transfer[data_type_t::INPUT]) {
            num_data_transfer[data_type_t::INPUT]++;

            std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

            unsigned num_access_pe = 0, num_access_pe_array = 0;
            for(unsigned i = 0; i < get_number_of_active_pes(); i++) { 
                uint64_t address_pe = 0, address_pe_array = 0;
                for(unsigned b = 0; b < parameters_pe[parameter_type_t::BATCH_SIZE]; b++) {
                    for(unsigned c = 0; c < parameters_pe[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned h = 0; h < parameters_pe[parameter_type_t::INPUT_HEIGHT]; h++) {
                            for(unsigned w = 0; w < parameters_pe[parameter_type_t::INPUT_WIDTH]; w++) {

                                if(address_pe != ((uint64_t)&pes[i]->input_data_lb[b*parameters_pe[parameter_type_t::INPUT_CHANNEL]
                                                                                    *parameters_pe[parameter_type_t::INPUT_HEIGHT]
                                                                                    *parameters_pe[parameter_type_t::INPUT_WIDTH] + 
                                                                                   c*parameters_pe[parameter_type_t::INPUT_HEIGHT]
                                                                                    *parameters_pe[parameter_type_t::INPUT_WIDTH] + 
                                                                                   h*parameters_pe[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                   pes[i]->mask_bits_lb[data_type_t::INPUT]) << pes[i]->mask_bits_lb[data_type_t::INPUT]) {

                                    // Update local buffer cost in PEs
                                    pes[i]->access_cycle_lb[data_type_t::INPUT] += pes[i]->u_write_cycle_lb[data_type_t::INPUT];
                                    pes[i]->access_energy_lb[data_type_t::INPUT] += pes[i]->u_write_energy_lb[data_type_t::INPUT];

                                    num_access_pe++;

                                    // Update local buffer address of PEs
                                    address_pe = ((uint64_t)&pes[i]->input_data_lb[b*parameters_pe[parameter_type_t::INPUT_CHANNEL]
                                                                                    *parameters_pe[parameter_type_t::INPUT_HEIGHT]
                                                                                    *parameters_pe[parameter_type_t::INPUT_WIDTH] + 
                                                                                   c*parameters_pe[parameter_type_t::INPUT_HEIGHT]
                                                                                    *parameters_pe[parameter_type_t::INPUT_WIDTH] + 
                                                                                   h*parameters_pe[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                   pes[i]->mask_bits_lb[data_type_t::INPUT]) << pes[i]->mask_bits_lb[data_type_t::INPUT];
                                }

                                // Check address of temporal buffer in PE array
                                if(!m_scheduler->read_tile_granular_pe_input[i%m_scheduler->input_offset_pe_array.size()] && 
                                   address_pe_array != ((uint64_t)&input_data[m_scheduler->input_offset_pe_array[i%m_scheduler->input_offset_pe_array.size()] + 
                                                                              b*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                               *parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                               *parameters_pe_array[parameter_type_t::INPUT_WIDTH] + 
                                                                              c*parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                               *parameters_pe_array[parameter_type_t::INPUT_WIDTH] + 
                                                                              h*parameters_pe_array[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                              mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT]) { 

                                    // Update temporal buffer cost in PE array
                                    access_cycle[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT];
                                    access_energy[data_type_t::INPUT] += u_read_energy[data_type_t::INPUT];
                                    num_access_pe_array++;

                                    // Update address of temporal buffer in PE array
                                    m_scheduler->read_tile_granular_pe_input[i%m_scheduler->input_offset_pe_array.size()] = true;
                                    address_pe_array = ((uint64_t)&input_data[m_scheduler->input_offset_pe_array[i%m_scheduler->input_offset_pe_array.size()] + 
                                                                              b*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                               *parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                               *parameters_pe_array[parameter_type_t::INPUT_WIDTH] + 
                                                                              c*parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                               *parameters_pe_array[parameter_type_t::INPUT_WIDTH] + 
                                                                              h*parameters_pe_array[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                              mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT];
                                }
                            }
                        }
                    }
                }
                pes[i]->skip_transfer[data_type_t::INPUT] = false;
            }

            // Update transfer cycle and energy
            transfer_energy[data_type_t::INPUT] += num_access_pe*noc_energy*ceil((double)(pes[0]->line_size_lb[data_type_t::INPUT])/(double)(bitwidth));
            transfer_cycle[data_type_t::INPUT] += num_access_pe_array*noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::INPUT])/(double)(bitwidth));

            // At the 1, 2, before last, and last stages
            double first_stage = u_read_cycle[data_type_t::INPUT];
            double second_stage = std::max(u_read_cycle[data_type_t::INPUT],
                                           noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::INPUT])/(double)(bitwidth)));

            double last_before_stage = std::max(pes[0]->u_write_cycle_lb[data_type_t::INPUT], 
                                                noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::INPUT])/(double)(bitwidth)));
            double last_stage = pes[0]->u_write_cycle_lb[data_type_t::INPUT];

            double other_stage = std::max(u_read_cycle[data_type_t::INPUT],
                                 std::max(pes[0]->u_write_cycle_lb[data_type_t::INPUT],
                                          noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::INPUT])/(double)(bitwidth))));

            // if temporal buffer is exist at the PE array
            if(exist_temporal_buffer) {
                // Update overlapped cycle of Scattering
                if(num_access_pe_array == 1) {
                    cycle_temporal_pe[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT] + 
                                                             noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::INPUT])/(double)(bitwidth)) + 
                                                             pes[0]->u_write_cycle_lb[data_type_t::INPUT];
               } else {
                   cycle_temporal_pe[data_type_t::INPUT] += first_stage + second_stage +
                                                            (num_access_pe_array-2)*other_stage + 
                                                            last_before_stage + last_stage; 
               }
            }
        }
#endif
        for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
            pes[i]->exist_data_lb[data_type_t::INPUT] = true, pes[i]->request_to_pe_array[data_type_t::INPUT] = false;
        }
    }

    bool request_to_pe_array_weight = false;
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::WEIGHT]) {
            request_to_pe_array_weight = true;
            break;
        }
    }

    if(request_to_pe_array_weight) {
#ifdef FUNCTIONAL
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                num_data_transfer[data_type_t::WEIGHT]++;

                std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
                parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

                unsigned num_access_pe = 0, num_access_pe_array = 0;

                for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
                    uint64_t address_pe = 0, address_pe_array = 0;
                    m_scheduler->transfer_data(pes[i]->weight_lb, weight, 0, m_scheduler->weight_offset_pe_array[i%m_scheduler->weight_offset_pe_array.size()], 
                                               component_type_t::PE, component_type_t::PE_Y, 
                                               data_type_t::WEIGHT, pes[i]->get_local_buffer_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
                    //m_scheduler->transfer_data_network_on_chip(pes[i]->weight_lb, weight, 
                    //                                           component_type_t::PE, component_type_t::PE_Y,
                    //                                           data_type_t::WEIGHT, pes[i]->get_local_buffer_stationary_type(), 
                    //                                           action_type_t::LOAD, last_component);
                    for(unsigned k = 0; k < parameters_pe[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                        for(unsigned c = 0; c < parameters_pe[parameter_type_t::INPUT_CHANNEL]; c++) {
                            for(unsigned r = 0; r < parameters_pe[parameter_type_t::FILTER_HEIGHT]; r++) {
                                for(unsigned s = 0; s < parameters_pe[parameter_type_t::FILTER_WIDTH]; s++) {

                                    if(address_pe != ((uint64_t)&pes[i]->weight_lb[k*parameters_pe[parameter_type_t::INPUT_CHANNEL]
                                                                                    *parameters_pe[parameter_type_t::FILTER_HEIGHT]
                                                                                    *parameters_pe[parameter_type_t::FILTER_WIDTH] + 
                                                                                   c*parameters_pe[parameter_type_t::FILTER_HEIGHT]
                                                                                    *parameters_pe[parameter_type_t::FILTER_WIDTH] + 
                                                                                   r*parameters_pe[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                   pes[i]->mask_bits_lb[data_type_t::WEIGHT]) << pes[i]->mask_bits_lb[data_type_t::WEIGHT]) {

                                        // Update cost of local buffer in PEs
                                        pes[i]->access_cycle_lb[data_type_t::WEIGHT] += pes[i]->u_write_cycle_lb[data_type_t::WEIGHT];
                                        pes[i]->access_energy_lb[data_type_t::WEIGHT] += pes[i]->u_write_energy_lb[data_type_t::WEIGHT];
                                        num_access_pe++;

                                        // Update address of local buffer in PEs
                                        address_pe = ((uint64_t)&pes[i]->weight_lb[k*parameters_pe[parameter_type_t::INPUT_CHANNEL]
                                                                                    *parameters_pe[parameter_type_t::FILTER_HEIGHT]
                                                                                    *parameters_pe[parameter_type_t::FILTER_WIDTH] + 
                                                                                   c*parameters_pe[parameter_type_t::FILTER_HEIGHT]
                                                                                    *parameters_pe[parameter_type_t::FILTER_WIDTH] + 
                                                                                   r*parameters_pe[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                   pes[i]->mask_bits_lb[data_type_t::WEIGHT]) << pes[i]->mask_bits_lb[data_type_t::WEIGHT];
                                    }

                                    // Check address of temporal buffer in PE array
                                    if(!m_scheduler->read_tile_granular_pe_weight[i%m_scheduler->weight_offset_pe_array.size()] &&
                                       address_pe_array != ((uint64_t)&weight[m_scheduler->weight_offset_pe_array[i%m_scheduler->weight_offset_pe_array.size()] + 
                                                                              k*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                               *parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                               *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                              c*parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                               *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                              r*parameters_pe_array[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                              mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT]) { 

                                        // Update cost of temporal buffer in PE array
                                        access_cycle[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT];
                                        access_energy[data_type_t::WEIGHT] += u_read_energy[data_type_t::WEIGHT];
                                        num_access_pe_array++;

                                        // Update address of temporal buffer in PE array
                                        m_scheduler->read_tile_granular_pe_weight[i%m_scheduler->weight_offset_pe_array.size()] = true;
                                        address_pe_array = ((uint64_t)&weight[m_scheduler->weight_offset_pe_array[i%m_scheduler->weight_offset_pe_array.size()] + 
                                                                              k*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                               *parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                               *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                              c*parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                               *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                              r*parameters_pe_array[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                              mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT];
                                    }
                                }
                            }
                        }
                    }
                    pes[i]->skip_transfer[data_type_t::WEIGHT] = false;
                }
                    
                    // Update transfer cycle and energy
                transfer_cycle[data_type_t::WEIGHT] += num_access_pe_array*noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth));
                transfer_energy[data_type_t::WEIGHT] += num_access_pe*noc_energy*ceil((double)(pes[0]->line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth));

                double first_stage = u_read_cycle[data_type_t::WEIGHT];
                double second_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                               noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)));

                double last_before_stage = std::max(pes[0]->u_write_cycle_lb[data_type_t::WEIGHT],
                                                    noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)));
                double last_stage = pes[0]->u_write_cycle_lb[data_type_t::WEIGHT];

                double other_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                     std::max(pes[0]->u_write_cycle_lb[data_type_t::WEIGHT],
                                              noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth))));

                // if temporal buffer is exist at the PE array
                if(exist_temporal_buffer) {
                    // Update overlapped cycle of Scattering
                    if(num_access_pe_array == 1) {
                        cycle_temporal_pe[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT] +
                                                                 noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)) +
                                                                 pes[0]->u_write_cycle_lb[data_type_t::WEIGHT];
                   } else {
                       cycle_temporal_pe[data_type_t::INPUT] += first_stage + second_stage +
                                                                (num_access_pe_array-2)*other_stage +
                                                                last_before_stage + last_stage;
                   }
                }
            }
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::PE);
                unsigned row = parameters[parameter_type_t::FILTER_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer[data_type_t::WEIGHT]++;

                for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
                    m_scheduler->transfer_data(pes[i]->weight_lb, weight, 0, m_scheduler->weight_offset_pe_array[i%m_scheduler->weight_offset_pe_array.size()], 
                                               component_type_t::PE, component_type_t::PE_Y, 
                                               data_type_t::WEIGHT, pes[i]->get_local_buffer_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
                    //m_scheduler->transfer_data_network_on_chip(pes[i]->weight_lb, weight, 
                    //                                           component_type_t::PE, component_type_t::PE_Y,
                    //                                           data_type_t::WEIGHT, pes[i]->get_local_buffer_stationary_type(), 
                    //                                           action_type_t::LOAD, last_component);

                    if(!m_scheduler->read_tile_granular_pe_weight[i%m_scheduler->weight_offset_pe_array.size()]) {
                        // Update access cost of PE array
                        access_cycle[data_type_t::WEIGHT] += (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                            *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                             (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                            *u_read_cycle[data_type_t::WEIGHT]
                                                            /(sizeof(data_t)*8/row_bit) + // Row index
                                                             parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                            *u_read_cycle[data_type_t::WEIGHT]
                                                            /(sizeof(data_t)*8/row_bit); // Column pointer

                        access_energy[data_type_t::WEIGHT] += (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                              (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_energy[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/row_bit) + // Row index
                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                            *u_read_energy[data_type_t::WEIGHT]
                                                            /(sizeof(data_t)*8/row_bit); // Column pointer


                        // Update transfer cycle between PE array and local buffer
                        transfer_cycle[data_type_t::WEIGHT] += noc_cycle
                                                             *ceil((double)(pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *8*sizeof(data_t)
                                                             /(double)bitwidth) + // Non-zero data
                                                              noc_cycle
                                                             *ceil((double)((pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*row_bit)
                                                             /(double)(bitwidth)) + // Row index
                                                              noc_cycle
                                                             *ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                           *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                           *(parameters[parameter_type_t::FILTER_WIDTH]+1)*row_bit)
                                                             /(double)(bitwidth)); // Column pointer

                        m_scheduler->read_tile_granular_pe_weight[i%m_scheduler->weight_offset_pe_array.size()] = true;
                        
                    }
                    
                    // Update local buffer access cost
                    pes[i]->access_cycle_lb[data_type_t::WEIGHT] += (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                   *pes[i]->u_write_cycle_lb[data_type_t::WEIGHT] + // Non-zero data
                                                                    (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                   *pes[i]->u_write_cycle_lb[data_type_t::WEIGHT]
                                                                   /(sizeof(data_t)*8/row_bit) + // Row index
                                                                    parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                   *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                                   *pes[i]->u_write_cycle_lb[data_type_t::WEIGHT]
                                                                   /(sizeof(data_t)*8/row_bit); // Column pointer

                    pes[i]->access_energy_lb[data_type_t::WEIGHT] += (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                    *pes[i]->u_write_energy_lb[data_type_t::WEIGHT] + // Non-zero data
                                                                     (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                    *pes[i]->u_write_energy_lb[data_type_t::WEIGHT]
                                                                    /(sizeof(data_t)*8/row_bit) + // Row index
                                                                     parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                    *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                                    *pes[i]->u_write_energy_lb[data_type_t::WEIGHT]
                                                                    /(sizeof(data_t)*8/row_bit); // Column pointer

                    // Update transfer energy between PE array and local buffer
                    transfer_energy[data_type_t::WEIGHT] += noc_energy
                                                           *ceil((double)(pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                           *8*sizeof(data_t)
                                                           /(double)bitwidth) + // Non-zero data
                                                            noc_energy
                                                           *ceil((double)((pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*row_bit)
                                                           /(double)(bitwidth)) + // Row index
                                                            noc_energy
                                                           *ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                         *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                         *(parameters[parameter_type_t::FILTER_WIDTH]+1)*row_bit)
                                                           /(double)(bitwidth)); // Column pointer

                    unsigned data_size = (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                        *sizeof(data_t) + // Non-zero data
                                         (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                        *row_bit/8 + // Row index
                                         parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                        *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                        *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                        *row_bit/8; // Column pointer

                    pes[i]->utilization_local_buffer[data_type_t::WEIGHT] = (double)data_size/(double)pes[i]->weight_size;

                    pes[i]->skip_transfer[data_type_t::WEIGHT] = false;
                }
            }
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::PE);
                unsigned column = parameters[parameter_type_t::FILTER_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer[data_type_t::WEIGHT]++;

                for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
                    m_scheduler->transfer_data(pes[i]->weight_lb, weight, 0, m_scheduler->weight_offset_pe_array[i%m_scheduler->weight_offset_pe_array.size()], 
                                               component_type_t::PE, component_type_t::PE_Y, 
                                               data_type_t::WEIGHT, pes[i]->get_local_buffer_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
                    //m_scheduler->transfer_data_network_on_chip(pes[i]->weight_lb, weight, 
                    //                                           component_type_t::PE, component_type_t::PE_Y,
                    //                                           data_type_t::WEIGHT, pes[i]->get_local_buffer_stationary_type(), 
                    //                                           action_type_t::LOAD, last_component);

                    if(!m_scheduler->read_tile_granular_pe_weight[i%m_scheduler->weight_offset_pe_array.size()]) {
                        // Update access cost of PE array
                        access_cycle[data_type_t::WEIGHT] += (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                            *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                             (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                            *u_read_cycle[data_type_t::WEIGHT]
                                                            /(sizeof(data_t)*8/column_bit) + // Column index
                                                             parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                            *u_read_cycle[data_type_t::WEIGHT]
                                                            /(sizeof(data_t)*8/column_bit); // Row pointer

                        access_energy[data_type_t::WEIGHT] += (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                              (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_energy[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/column_bit) + // Column index
                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                            *u_read_energy[data_type_t::WEIGHT]
                                                            /(sizeof(data_t)*8/column_bit); // Row pointer


                        // Update transfer cycle between PE array and local buffer
                        transfer_cycle[data_type_t::WEIGHT] += noc_cycle
                                                             *ceil((double)(pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *8*sizeof(data_t)
                                                             /(double)bitwidth) + // Non-zero data
                                                              noc_cycle
                                                             *ceil((double)((pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*column_bit)
                                                             /(double)(bitwidth)) + // Column index
                                                              noc_cycle
                                                             *ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                           *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                           *(parameters[parameter_type_t::FILTER_HEIGHT]+1)*column_bit)
                                                             /(double)(bitwidth)); // Row pointer

                        m_scheduler->read_tile_granular_pe_weight[i%m_scheduler->weight_offset_pe_array.size()] = true;
                        
                    }
                    
                    // Update local buffer access cost
                    pes[i]->access_cycle_lb[data_type_t::WEIGHT] += (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                   *pes[i]->u_write_cycle_lb[data_type_t::WEIGHT] + // Non-zero data
                                                                    (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                   *pes[i]->u_write_cycle_lb[data_type_t::WEIGHT]
                                                                   /(sizeof(data_t)*8/column_bit) + // Column index
                                                                    parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                   *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                                   *pes[i]->u_write_cycle_lb[data_type_t::WEIGHT]
                                                                   /(sizeof(data_t)*8/column_bit); // Row pointer

                    pes[i]->access_energy_lb[data_type_t::WEIGHT] += (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                    *pes[i]->u_write_energy_lb[data_type_t::WEIGHT] + // Non-zero data
                                                                     (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                    *pes[i]->u_write_energy_lb[data_type_t::WEIGHT]
                                                                    /(sizeof(data_t)*8/column_bit) + // Column index
                                                                     parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                    *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                                    *pes[i]->u_write_energy_lb[data_type_t::WEIGHT]
                                                                    /(sizeof(data_t)*8/column_bit); // Row pointer

                    // Update transfer energy between PE array and local buffer
                    transfer_energy[data_type_t::WEIGHT] += noc_energy
                                                         *ceil((double)(pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                         *8*sizeof(data_t)
                                                         /(double)bitwidth) + // Non-zero data
                                                          noc_energy
                                                         *ceil((double)((pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*column_bit)
                                                         /(double)(bitwidth)) + // Column index
                                                          noc_energy
                                                         *ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                       *(parameters[parameter_type_t::FILTER_HEIGHT]+1)*column_bit)
                                                         /(double)(bitwidth)); // Row pointer

                    unsigned data_size = (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                        *sizeof(data_t) + // Non-zero data
                                         (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                        *column_bit/8 + // Column index
                                         parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                        *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                        *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                        *column_bit/8; // Row pointer

                    pes[i]->utilization_local_buffer[data_type_t::WEIGHT] = (double)data_size/(double)pes[i]->weight_size;


                    pes[i]->skip_transfer[data_type_t::WEIGHT] = false;
                }
            }
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::WEIGHT]) {

                num_data_transfer[data_type_t::WEIGHT]++;

                for(unsigned i = 0; i < get_number_of_active_pes(); i++) {

                    m_scheduler->transfer_data(pes[i]->weight_lb, weight, 0, m_scheduler->weight_offset_pe_array[i%m_scheduler->weight_offset_pe_array.size()], 
                                               component_type_t::PE, component_type_t::PE_Y, 
                                               data_type_t::WEIGHT, pes[i]->get_local_buffer_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
                    //m_scheduler->transfer_data_network_on_chip(pes[i]->weight_lb, weight, 
                    //                                           component_type_t::PE, component_type_t::PE_Y,
                    //                                           data_type_t::WEIGHT, pes[i]->get_local_buffer_stationary_type(), 
                    //                                           action_type_t::LOAD, last_component);
                    
                    if(!m_scheduler->read_tile_granular_pe_weight[i%m_scheduler->weight_offset_pe_array.size()]) { 
                        // Update access cost of PE array
                        access_cycle[data_type_t::WEIGHT] += (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                            *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                             pes[i]->tile_size_lb[data_type_t::WEIGHT]
                                                            *u_read_cycle[data_type_t::WEIGHT]
                                                            /(sizeof(data_t)*8); // Metadata
                       
                        access_energy[data_type_t::WEIGHT] += (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                              pes[i]->tile_size_lb[data_type_t::WEIGHT]
                                                             *u_read_energy[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8); // Metadata

                        m_scheduler->read_tile_granular_pe_weight[i%m_scheduler->weight_offset_pe_array.size()] = true;

                        // Update transfer cycle between the chip-level processor and the global buffer
                        transfer_cycle[data_type_t::WEIGHT] += noc_cycle
                                                             *ceil((double)((pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *8*sizeof(data_t))
                                                             /(double)bitwidth) +  // Non-zero data
                                                              noc_cycle 
                                                             *ceil((double)(pes[i]->tile_size_lb[data_type_t::WEIGHT])/(double)bitwidth); // Metadata
                    }
                    // Update local buffer access cost
                    pes[i]->access_cycle_lb[data_type_t::WEIGHT] += (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                   *pes[i]->u_write_cycle_lb[data_type_t::WEIGHT] + // Non-zero data
                                                                    pes[i]->tile_size_lb[data_type_t::WEIGHT]
                                                                   *pes[i]->u_write_cycle_lb[data_type_t::WEIGHT]
                                                                   /(sizeof(data_t)*8); // Metadata
                   
                    pes[i]->access_energy_lb[data_type_t::WEIGHT] += (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                    *pes[i]->u_write_energy_lb[data_type_t::WEIGHT] + // Non-zero data
                                                                     pes[i]->tile_size_lb[data_type_t::WEIGHT]
                                                                    *pes[i]->u_write_energy_lb[data_type_t::WEIGHT]
                                                                    /(sizeof(data_t)*8); // Metadata

                    // Update transfer energy between PE array and local buffers 
                    transfer_energy[data_type_t::WEIGHT] += noc_energy
                                                           *ceil((double)((pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                           *8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                            noc_energy
                                                           *ceil((double)(pes[i]->tile_size_lb[data_type_t::WEIGHT])/(double)bitwidth); // Metadata

                    unsigned data_size = (pes[i]->tile_size_lb[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*sizeof(data_t) + //Non-zero data
                                         pes[i]->tile_size_lb[data_type_t::WEIGHT]/8; // Metadata

                    pes[i]->utilization_local_buffer[data_type_t::WEIGHT] = (double)data_size/(double)pes[i]->weight_size;

                    pes[i]->skip_transfer[data_type_t::WEIGHT] = false;
                }
            }
        }
        else {
            std::cerr << "Undefined compression type" << std::endl;
            exit(1);
        }

#else
        if(!skip_transfer[data_type_t::WEIGHT]) {
            num_data_transfer[data_type_t::WEIGHT]++;

            std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

            unsigned num_access_pe = 0, num_access_pe_array = 0;

            for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
                uint64_t address_pe = 0, address_pe_array = 0;
                for(unsigned k = 0; k < parameters_pe[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned c = 0; c < parameters_pe[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned r = 0; r < parameters_pe[parameter_type_t::FILTER_HEIGHT]; r++) {
                            for(unsigned s = 0; s < parameters_pe[parameter_type_t::FILTER_WIDTH]; s++) {

                                if(address_pe != ((uint64_t)&pes[i]->weight_lb[k*parameters_pe[parameter_type_t::INPUT_CHANNEL]
                                                                                *parameters_pe[parameter_type_t::FILTER_HEIGHT]
                                                                                *parameters_pe[parameter_type_t::FILTER_WIDTH] + 
                                                                               c*parameters_pe[parameter_type_t::FILTER_HEIGHT]
                                                                                *parameters_pe[parameter_type_t::FILTER_WIDTH] + 
                                                                               r*parameters_pe[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                               pes[i]->mask_bits_lb[data_type_t::WEIGHT]) << pes[i]->mask_bits_lb[data_type_t::WEIGHT]) {

                                    // Update cost of local buffer in PEs
                                    pes[i]->access_cycle_lb[data_type_t::WEIGHT] += pes[i]->u_write_cycle_lb[data_type_t::WEIGHT];
                                    pes[i]->access_energy_lb[data_type_t::WEIGHT] += pes[i]->u_write_energy_lb[data_type_t::WEIGHT];
                                    num_access_pe++;

                                    // Update address of local buffer in PEs
                                    address_pe = ((uint64_t)&pes[i]->weight_lb[k*parameters_pe[parameter_type_t::INPUT_CHANNEL]
                                                                                *parameters_pe[parameter_type_t::FILTER_HEIGHT]
                                                                                *parameters_pe[parameter_type_t::FILTER_WIDTH] + 
                                                                               c*parameters_pe[parameter_type_t::FILTER_HEIGHT]
                                                                                *parameters_pe[parameter_type_t::FILTER_WIDTH] + 
                                                                               r*parameters_pe[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                               pes[i]->mask_bits_lb[data_type_t::WEIGHT]) << pes[i]->mask_bits_lb[data_type_t::WEIGHT];
                                }

                                // Check address of temporal buffer in PE array
                                if(!m_scheduler->read_tile_granular_pe_weight[i%m_scheduler->weight_offset_pe_array.size()] &&
                                   address_pe_array != ((uint64_t)&weight[m_scheduler->weight_offset_pe_array[i%m_scheduler->weight_offset_pe_array.size()] + 
                                                                          k*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                           *parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                           *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                          c*parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                           *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                          r*parameters_pe_array[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                          mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT]) { 

                                    // Update cost of temporal buffer in PE array
                                    access_cycle[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT];
                                    access_energy[data_type_t::WEIGHT] += u_read_energy[data_type_t::WEIGHT];
                                    num_access_pe_array++;

                                    // Update address of temporal buffer in PE array
                                    m_scheduler->read_tile_granular_pe_weight[i%m_scheduler->weight_offset_pe_array.size()] = true;
                                    address_pe_array = ((uint64_t)&weight[m_scheduler->weight_offset_pe_array[i%m_scheduler->weight_offset_pe_array.size()] + 
                                                                          k*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                           *parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                           *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                          c*parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                           *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                          r*parameters_pe_array[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                          mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT];
                                }
                            }
                        }
                    }
                }
                pes[i]->skip_transfer[data_type_t::WEIGHT] = false;
            }
                
                // Update transfer cycle and energy
            transfer_cycle[data_type_t::WEIGHT] += num_access_pe_array*noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth));
            transfer_energy[data_type_t::WEIGHT] += num_access_pe*noc_energy*ceil((double)(pes[0]->line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth));

            double first_stage = u_read_cycle[data_type_t::WEIGHT];
            double second_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                           noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)));

            double last_before_stage = std::max(pes[0]->u_write_cycle_lb[data_type_t::WEIGHT],
                                                noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)));
            double last_stage = pes[0]->u_write_cycle_lb[data_type_t::WEIGHT];

            double other_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                 std::max(pes[0]->u_write_cycle_lb[data_type_t::WEIGHT],
                                          noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth))));

            // if temporal buffer is exist at the PE array
            if(exist_temporal_buffer) {
                // Update overlapped cycle of Scattering
                if(num_access_pe_array == 1) {
                    cycle_temporal_pe[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT] +
                                                             noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)) +
                                                             pes[0]->u_write_cycle_lb[data_type_t::WEIGHT];
               } else {
                   cycle_temporal_pe[data_type_t::INPUT] += first_stage + second_stage +
                                                            (num_access_pe_array-2)*other_stage +
                                                            last_before_stage + last_stage;
               }
            }
        }
#endif
        for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
            pes[i]->exist_data_lb[data_type_t::WEIGHT] = true, pes[i]->request_to_pe_array[data_type_t::WEIGHT] = false;
        }
    }

    bool request_to_pe_array_output = false;
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::OUTPUT]) {
            request_to_pe_array_output = true;
            break;
        }
    }


    if(request_to_pe_array_output) {
        if(!skip_transfer[data_type_t::OUTPUT]) {
            /* Stats */
            if(global_buffer->transfer_output) {
                num_data_transfer[data_type_t::OUTPUT]++;

                std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
                parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

                unsigned num_access_pe = 0, num_access_pe_array = 0;

                for(unsigned i = 0; i < get_number_of_active_pes(); i++) {

                    uint64_t address_pe = 0, address_pe_array = 0;
#ifdef FUNCTIONAL
                    // Transfer weight from temporal weight buffer of PE array to PE.
                    m_scheduler->transfer_data(pes[i]->output_data_lb, output_data, 0, m_scheduler->output_offset_pe_array[i%m_scheduler->output_offset_pe_array.size()], 
                                               component_type_t::PE, component_type_t::PE_Y, 
                                               data_type_t::OUTPUT, pes[i]->get_local_buffer_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
                    //m_scheduler->transfer_data_network_on_chip(pes[i]->output_data_lb, output_data, 
                    //                                           component_type_t::PE, component_type_t::PE_Y,
                    //                                           data_type_t::OUTPUT, pes[i]->get_local_buffer_stationary_type(), 
                    //                                           action_type_t::LOAD, last_component);
#endif
                    for(unsigned b = 0; b < parameters_pe[parameter_type_t::BATCH_SIZE]; b++) {
                        for(unsigned k = 0; k < parameters_pe[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                            for(unsigned p = 0; p < parameters_pe[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                                for(unsigned q = 0; q < parameters_pe[parameter_type_t::OUTPUT_WIDTH]; q++) {

                                    // Check address of local buffer in PEs
                                    if(address_pe != ((uint64_t)&pes[i]->output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                        pes[i]->mask_bits_lb[data_type_t::OUTPUT]) << pes[i]->mask_bits_lb[data_type_t::OUTPUT]) {

                                        // Update cost of local buffer in PEs
                                        pes[i]->access_cycle_lb[data_type_t::OUTPUT] += pes[i]->u_write_cycle_lb[data_type_t::OUTPUT];
                                        pes[i]->access_energy_lb[data_type_t::OUTPUT] += pes[i]->u_write_energy_lb[data_type_t::OUTPUT];
                                        num_access_pe++;

                                        // Update address of local buffer in PEs
                                        address_pe = ((uint64_t)&pes[i]->output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                        pes[i]->mask_bits_lb[data_type_t::OUTPUT]) << pes[i]->mask_bits_lb[data_type_t::OUTPUT];
                                    }

                                    // Check address of temporal buffer in PE array
                                    if(!m_scheduler->read_tile_granular_pe_output[i%m_scheduler->output_offset_pe_array.size()] &&
                                       address_pe_array != ((uint64_t)&output_data[m_scheduler->output_offset_pe_array[i%m_scheduler->output_offset_pe_array.size()] + 
                                                                                   b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                    *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                   mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT]) { 

                                        // Update cost of temporal buffer in PE array
                                        access_cycle[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT];
                                        access_energy[data_type_t::OUTPUT] += u_read_energy[data_type_t::OUTPUT];
                                        num_access_pe_array++;

                                        // Update address of temporal buffer in PE array
                                        m_scheduler->read_tile_granular_pe_output[i%m_scheduler->output_offset_pe_array.size()] = true;
                                        address_pe_array = ((uint64_t)&output_data[m_scheduler->output_offset_pe_array[i%m_scheduler->output_offset_pe_array.size()] + 
                                                                                   b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                    *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                    *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                   p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                   mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT];
                                    }
                                }
                            }
                        }
                    }
                    pes[i]->skip_transfer[data_type_t::OUTPUT] = false;
                }

                transfer_cycle[data_type_t::OUTPUT] += num_access_pe_array*noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth));
                transfer_energy[data_type_t::OUTPUT] += num_access_pe*noc_energy*ceil((double)(pes[0]->line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth));

                double first_stage = u_read_cycle[data_type_t::OUTPUT];
                double second_stage = std::max(u_read_cycle[data_type_t::OUTPUT],
                                               noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
                double last_before_stage = std::max(pes[0]->u_write_cycle_lb[data_type_t::OUTPUT],
                                                    noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
                double last_stage = pes[0]->u_write_cycle_lb[data_type_t::OUTPUT];


                double other_stage = std::max(u_read_cycle[data_type_t::OUTPUT],
                                     std::max(pes[0]->u_write_cycle_lb[data_type_t::OUTPUT],
                                              noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth))));
                                // if temporal buffer is exist at the PE array
                if(exist_temporal_buffer) {
                    // Update overlapped cycle of Scattering
                    if(num_access_pe_array == 1) {
                        cycle_temporal_pe[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT] +
                                                                 noc_cycle*ceil((double)(pes[0]->line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)) +
                                                                 pes[0]->u_write_cycle_lb[data_type_t::OUTPUT];
                    } else {
                        cycle_temporal_pe[data_type_t::OUTPUT] += first_stage + second_stage +
                                                                 (num_access_pe_array-2)*other_stage +
                                                                 last_before_stage + last_stage;
                    }
                }
            }
        }
        for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
            pes[i]->exist_data_lb[data_type_t::OUTPUT] = true, pes[i]->request_to_pe_array[data_type_t::OUTPUT] = false;
        }
    }
        
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        pes[i]->fill_data();
    }
    /*
    for(unsigned i = get_number_of_active_pes(); i < get_number_of_pes(); i++) {
        pes[i]->static_energy[data_type_t::INPUT] += pes[i]->u_static_energy[data_type_t::INPUT];
        pes[i]->static_energy[data_type_t::WEIGHT] += pes[i]->u_static_energy[data_type_t::WEIGHT];
        pes[i]->static_energy[data_type_t::OUTPUT] += pes[i]->u_static_energy[data_type_t::OUTPUT];
    }
    */

    if(tile_size[data_type_t::INPUT] == global_buffer->tile_size[data_type_t::INPUT]) {skip_transfer[data_type_t::INPUT] = true;}
    if(tile_size[data_type_t::WEIGHT] == global_buffer->tile_size[data_type_t::WEIGHT]) {skip_transfer[data_type_t::WEIGHT] = true;}
    if(tile_size[data_type_t::OUTPUT] == global_buffer->tile_size[data_type_t::OUTPUT]) {skip_transfer[data_type_t::OUTPUT] = true;}

#ifdef PRINT
    std::cout << "Transfer data from temporal buffer in PE array to each PE" << std::endl;

#endif

}

// Print out the specification of PE array.
void adder_tree_t::print_specification() {
	pes[0]->print_specification();
    
    std::cout << "========== PE array specification ==========" << std::endl;
    std::cout << "PE array type      :" << std::setw(24) 
                                        << "Adder Tree" << std::endl;
    std::cout << "Array size         :" << std::setw(24) 
                                        << num_pes << std::endl;
    // Print the stationary type between the PE array and Global buffer.
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

    std::cout << "NoC type           :" << std::setw(24) 
                                        << "Adder Tree" << std::endl;
    std::cout << "NoC cycle          :" << std::setw(17) 
                                        << noc_cycle << " cycles" << std::endl;
    std::cout << "NoC energy         :" << std::setw(21) 
                                        << noc_energy << " pJ" << std::endl;
	std::cout << std::endl;
}

