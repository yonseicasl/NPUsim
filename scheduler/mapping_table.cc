#include <iostream>
#include "mapping_table.h"

mapping_table_t::mapping_table_t(section_config_t m_section_config) {	
    init(m_section_config);
}

mapping_table_t::~mapping_table_t() {

}

// Initialize the mapping table.
void mapping_table_t::init(section_config_t m_section_config) {

	mapping_table.reserve(component_type_t::NUM_COMPONENT_TYPES);

	std::vector<unsigned> parameters(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    // Convert mapping values from string.
	for(unsigned i = 0; i < component_type_str.size() - 1; i++) {
		parameters.assign(parameters.size(), 1);
		if(m_section_config.exists(component_type_str[i])) {
			parameters = get_value(m_section_config, component_type_str[i]);
		}
		mapping_table.emplace_back(parameters);
	}
}

// Read parameters and fill values to mapping table.
std::vector<unsigned> mapping_table_t::get_value(section_config_t m_section_config, std::string m_name) {
	std::string line;
	std::vector<unsigned> parameters;
	parameters.reserve(parameter_type_t::NUM_PARAMETER_TYPES);
	m_section_config.get_setting(m_name, &line);

	line.erase(remove(line.begin(), line.end(), ' '),line.end());
	unsigned found = 0;
	for(unsigned i = 0; i < parameter_type_t::NUM_PARAMETER_TYPES; i++) {
		found = line.find(',');
		parameters.emplace_back(stoi(line.substr(0, found)));
		line.erase(0, found + 1);
	}

	return parameters;
}

// Calculate the parameter size.
std::vector<unsigned> mapping_table_t::calculate_parameter_size(component_type_t m_component_type) {
    std::vector<unsigned> parameters(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    for(unsigned i = 0; i <= m_component_type; i++) {
        for(unsigned j = 0; j < parameter_type_t::NUM_PARAMETER_TYPES; j++) {
            parameters[j] *= mapping_table[i][j];
        }
    }
    parameters[parameter_type_t::STRIDE]       = mapping_table[m_component_type][parameter_type_t::STRIDE];
    parameters[parameter_type_t::INPUT_HEIGHT] = (parameters[parameter_type_t::OUTPUT_HEIGHT] - 1)*parameters[parameter_type_t::STRIDE] + parameters[parameter_type_t::FILTER_HEIGHT];
    parameters[parameter_type_t::INPUT_WIDTH]  = (parameters[parameter_type_t::OUTPUT_WIDTH] - 1)*parameters[parameter_type_t::STRIDE] + parameters[parameter_type_t::FILTER_WIDTH];

    return parameters;
}

std::vector<unsigned> mapping_table_t::calculate_parameter_size(component_type_t m_component_type, std::string m_parameter_order) {
    std::vector<unsigned> parameters(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    char t_parameter_order[parameter_type_t::NUM_PARAMETER_TYPES];
    strcpy(t_parameter_order, m_parameter_order.c_str());

    for(unsigned i = 0; i < parameter_type_t::NUM_PARAMETER_TYPES; i++) {
        std::cout << t_parameter_order[i] << " " << std::endl;
    }

    for(unsigned i = 0; i <= m_component_type; i++) {
        for(unsigned j = 0; j < parameter_type_t::NUM_PARAMETER_TYPES; j++) {
            if(t_parameter_order[j] == 'K' || t_parameter_order[j] == 'k') {
                parameters[j] *= mapping_table[i][parameter_type_t::OUTPUT_CHANNEL];
            } else if(t_parameter_order[j] == 'B' || t_parameter_order[j] == 'b') {
                parameters[j] *= mapping_table[i][parameter_type_t::BATCH_SIZE];
            } else if(t_parameter_order[j] == 'P' || t_parameter_order[j] == 'P') {
                parameters[j] *= mapping_table[i][parameter_type_t::OUTPUT_HEIGHT];
            } else if(t_parameter_order[j] == 'Q' || t_parameter_order[j] == 'q') {
                parameters[j] *= mapping_table[i][parameter_type_t::OUTPUT_WIDTH];
            } else if(t_parameter_order[j] == 'C' || t_parameter_order[j] == 'c') {
                parameters[j] *= mapping_table[i][parameter_type_t::INPUT_CHANNEL];
            } else if(t_parameter_order[j] == 'R' || t_parameter_order[j] == 'r') {
                parameters[j] *= mapping_table[i][parameter_type_t::FILTER_HEIGHT];
            } else if(t_parameter_order[j] == 'S' || t_parameter_order[j] == 's') {
                parameters[j] *= mapping_table[i][parameter_type_t::FILTER_WIDTH];
            }
        }
    }
    parameters[parameter_type_t::STRIDE]       = mapping_table[m_component_type][parameter_type_t::STRIDE];
    parameters[parameter_type_t::INPUT_HEIGHT] = (parameters[parameter_type_t::OUTPUT_HEIGHT] - 1)*parameters[parameter_type_t::STRIDE] + parameters[parameter_type_t::FILTER_HEIGHT];
    parameters[parameter_type_t::INPUT_WIDTH]  = (parameters[parameter_type_t::OUTPUT_WIDTH] - 1)*parameters[parameter_type_t::STRIDE] + parameters[parameter_type_t::FILTER_WIDTH];

    return parameters;

}

void mapping_table_t::calculate_num_tile_granular_data(component_type_t m_component_type, std::vector<unsigned> *m_tile_granular_data) {


    std::vector<unsigned> parameters(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    for(unsigned j = 0; j < parameter_type_t::NUM_PARAMETER_TYPES; j++) {
        parameters[j] *= mapping_table[m_component_type][j];
    }
    parameters[parameter_type_t::STRIDE]       = mapping_table[m_component_type][parameter_type_t::STRIDE];
    parameters[parameter_type_t::INPUT_HEIGHT] = (parameters[parameter_type_t::OUTPUT_HEIGHT] - 1)*parameters[parameter_type_t::STRIDE] + parameters[parameter_type_t::FILTER_HEIGHT];
    parameters[parameter_type_t::INPUT_WIDTH]  = (parameters[parameter_type_t::OUTPUT_WIDTH] - 1)*parameters[parameter_type_t::STRIDE] + parameters[parameter_type_t::FILTER_WIDTH];

    m_tile_granular_data->at(data_type_t::INPUT) *= parameters[parameter_type_t::BATCH_SIZE]
                                                   *parameters[parameter_type_t::INPUT_CHANNEL]
                                                   *parameters[parameter_type_t::INPUT_HEIGHT]
                                                   *parameters[parameter_type_t::INPUT_WIDTH];

    m_tile_granular_data->at(data_type_t::WEIGHT) *= parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                    *parameters[parameter_type_t::INPUT_CHANNEL]
                                                    *parameters[parameter_type_t::FILTER_HEIGHT]
                                                    *parameters[parameter_type_t::FILTER_WIDTH];

    m_tile_granular_data->at(data_type_t::OUTPUT) *= parameters[parameter_type_t::BATCH_SIZE]
                                                    *parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                    *parameters[parameter_type_t::OUTPUT_HEIGHT]
                                                    *parameters[parameter_type_t::OUTPUT_WIDTH];

}

// Calculate the number of active components (i.e., MAC, PE, Chips)
unsigned mapping_table_t::calculate_active_component(component_type_t m_component_type) {
    unsigned num_comp = 1;
    for(unsigned i = 0; i < 7; i++) {
        num_comp *= mapping_table[m_component_type][i];
    }
    num_comp /= mapping_table[m_component_type][parameter_type_t::GROUP];
    return num_comp;
}

// Print out the value of mapping table.
// For Debugging.
void mapping_table_t::print() {

	for(unsigned i = 0; i < mapping_table.size(); i++) {
		for(unsigned j = 0; j < mapping_table[i].size(); j++) {
			std::cout << mapping_table[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

