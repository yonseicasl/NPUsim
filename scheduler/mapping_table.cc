#include <iostream>
#include <cstdint>
#include <limits>
#include <sstream>
#include "mapping_table.h"
namespace {

unsigned checked_multiply(unsigned lhs, unsigned rhs, const char *context) {
    const uint64_t product = static_cast<uint64_t>(lhs) * rhs;
    if(product > std::numeric_limits<unsigned>::max()) {
        std::cerr << "Error: unsigned overflow while calculating " << context << std::endl;
        exit(1);
    }
    return static_cast<unsigned>(product);
}

unsigned checked_extent(unsigned output, unsigned stride, unsigned filter, const char *context) {
    if(output == 0 || stride == 0 || filter == 0) {
        std::cerr << "Error: zero dimension while calculating " << context << std::endl;
        exit(1);
    }
    const uint64_t extent = static_cast<uint64_t>(output - 1) * stride + filter;
    if(extent > std::numeric_limits<unsigned>::max()) {
        std::cerr << "Error: unsigned overflow while calculating " << context << std::endl;
        exit(1);
    }
    return static_cast<unsigned>(extent);
}

} // namespace

mapping_table_t::mapping_table_t(section_config_t m_section_config) {	
    init(m_section_config);
}

mapping_table_t::~mapping_table_t() {

}

// Initialize the mapping table.
void mapping_table_t::init(section_config_t m_section_config) {

	mapping_table.reserve(component_type_t::NUM_COMPONENT_TYPES);

    std::vector<unsigned> parameters(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    legacy_global_buffer_mapping.assign(parameters.size(), 1);

    // Convert mapping values from string.
    for(unsigned i = 0; i < component_type_str.size() - 1; i++) {
        parameters.assign(parameters.size(), 1);
        if(m_section_config.exists(component_type_str[i])) {
            parameters = get_value(m_section_config, component_type_str[i]);
        }
        // Legacy mapping files use "GLB" for an outer temporal repetition
        // factor.  Retain it for accounting without changing the historical
        // lower-level tile scheduler, which expects "global_buffer".
        else if(i == component_type_t::GLOBAL_BUFFER && m_section_config.exists("glb")) {
            legacy_global_buffer_mapping = get_value(m_section_config, "glb");
        }
        mapping_table.emplace_back(parameters);
    }
}
std::vector<unsigned> mapping_table_t::get_value(section_config_t m_section_config, std::string m_name) {
    std::string line;
    if(!m_section_config.get_setting(m_name, &line)) {
        std::cerr << "Error: missing mapping entry '" << m_name << "'" << std::endl;
        exit(1);
    }
    std::stringstream stream(line);
    std::vector<unsigned> parameters;
    parameters.reserve(parameter_type_t::NUM_PARAMETER_TYPES);
    for(unsigned i = 0; i < parameter_type_t::NUM_PARAMETER_TYPES; i++) {
        std::string token;
        if(!std::getline(stream, token, ',') || token.empty()) {
            std::cerr << "Error: mapping entry '" << m_name << "' requires "
                      << parameter_type_t::NUM_PARAMETER_TYPES << " values" << std::endl;
            exit(1);
        }
        try {
            size_t parsed = 0;
            unsigned long value = std::stoul(token, &parsed);
            if(parsed != token.size() || value > std::numeric_limits<unsigned>::max()) {
                throw std::out_of_range("mapping value");
            }
            parameters.emplace_back((unsigned)value);
            if(value == 0 && (i <= parameter_type_t::FILTER_WIDTH ||
                              i == parameter_type_t::GROUP || i == parameter_type_t::STRIDE)) {
                std::cerr << "Error: zero mapping factor in '" << m_name << "'" << std::endl;
                exit(1);
            }
        } catch(const std::exception &) {
            std::cerr << "Error: invalid unsigned mapping value '" << token
                      << "' in '" << m_name << "'" << std::endl;
            exit(1);
        }
    }
    std::string remainder;
    std::getline(stream, remainder, '\0');
    if(!remainder.empty()) {
        std::cerr << "Error: too many values in mapping entry '" << m_name << "'" << std::endl;
        exit(1);
    }
    return parameters;
}

// Calculate the parameter size.
std::vector<unsigned> mapping_table_t::calculate_parameter_size(component_type_t m_component_type) {
    std::vector<unsigned> parameters(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    for(unsigned i = 0; i <= m_component_type; i++) {
        for(unsigned j = 0; j < parameter_type_t::NUM_PARAMETER_TYPES; j++) {
            parameters[j] = checked_multiply(parameters[j], mapping_table[i][j], "parameter size");
        }
    }
    parameters[parameter_type_t::STRIDE]       = mapping_table[m_component_type][parameter_type_t::STRIDE];
    parameters[parameter_type_t::INPUT_HEIGHT] =
        checked_extent(parameters[parameter_type_t::OUTPUT_HEIGHT], parameters[parameter_type_t::STRIDE],
                       parameters[parameter_type_t::FILTER_HEIGHT], "input height");
    parameters[parameter_type_t::INPUT_WIDTH] =
        checked_extent(parameters[parameter_type_t::OUTPUT_WIDTH], parameters[parameter_type_t::STRIDE],
                       parameters[parameter_type_t::FILTER_WIDTH], "input width");

    return parameters;
}


void mapping_table_t::calculate_num_tile_granular_data(component_type_t m_component_type, std::vector<unsigned> *m_tile_granular_data) {


    std::vector<unsigned> parameters(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    for(unsigned j = 0; j < parameter_type_t::NUM_PARAMETER_TYPES; j++) {
        parameters[j] = checked_multiply(parameters[j], mapping_table[m_component_type][j], "tile parameter size");
    }
    parameters[parameter_type_t::STRIDE]       = mapping_table[m_component_type][parameter_type_t::STRIDE];
    parameters[parameter_type_t::INPUT_HEIGHT] =
        checked_extent(parameters[parameter_type_t::OUTPUT_HEIGHT], parameters[parameter_type_t::STRIDE],
                       parameters[parameter_type_t::FILTER_HEIGHT], "tile input height");
    parameters[parameter_type_t::INPUT_WIDTH] =
        checked_extent(parameters[parameter_type_t::OUTPUT_WIDTH], parameters[parameter_type_t::STRIDE],
                       parameters[parameter_type_t::FILTER_WIDTH], "tile input width");

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
    const std::vector<unsigned> &component_mapping =
        (m_component_type == component_type_t::GLOBAL_BUFFER)
            ? legacy_global_buffer_mapping
            : mapping_table[m_component_type];
    uint64_t num_comp = 1;
    for(unsigned i = 0; i < 7; i++) {
        num_comp *= component_mapping[i];
        if(num_comp > std::numeric_limits<unsigned>::max()) {
            std::cerr << "Error: unsigned overflow while calculating active components" << std::endl;
            exit(1);
        }
    }
    const unsigned group = component_mapping[parameter_type_t::GROUP];
    if(group == 0) {
        std::cerr << "Error: mapping GROUP must be non-zero" << std::endl;
        exit(1);
    }
    return static_cast<unsigned>(num_comp / group);
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
