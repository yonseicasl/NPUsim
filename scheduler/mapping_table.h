#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <cstring>
#include "config.h"
#include "def.h"
//#include "convolutional_layer.h"

class mapping_table_t {

public:
	// Constructor.
	mapping_table_t(section_config_t m_section_config);
	// Destructor.
	~mapping_table_t();
	// Initialize the mapping table.
	void init(section_config_t m_section_config);
    // Read parameter values and fill values to mapping table.
	std::vector<unsigned> get_value(section_config_t m_section_config, std::string m_name);
    // Calculate the parameter size.
    std::vector<unsigned> calculate_parameter_size(component_type_t m_component_type);

    std::vector<unsigned> calculate_parameter_size(component_type_t m_component_type, std::string m_parameter_order);

    // Calculate the number of tile-granular data 
    void calculate_num_tile_granular_data(component_type_t m_component_type, std::vector<unsigned> *m_tile_granumlar_data);

    // Calculate the number of active component (e.g., MAC, PE, and Chips)
    unsigned calculate_active_component(component_type_t m_component_type);
    // Print out the component of mapping table.
	void print();

private:
	std::vector<std::vector<unsigned>> mapping_table;   // The mapping table of neural layer.

};

#endif
