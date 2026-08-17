#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <cstring>
#include <string>
#include <vector>

#include "config.h"
#include "def.h"

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
    // Full per-dimension coverage of the mapping: the DRAM-level cumulative product
    // times the legacy GLB temporal-repetition factors (K/B/P/Q/C/R/S only).
    std::vector<unsigned> calculate_total_parameter_size();
    // Per-datatype GLB temporal-repetition factor: the product of only those legacy
    // GLB loop factors that index the datatype (a repetition over a dimension a
    // tensor does not depend on revisits the SAME tile, which the on-chip buffers
    // retain -- off-chip traffic must not scale with it).
    std::vector<unsigned> datatype_repetitions();

    // Calculate the number of tile-granular data 
    void calculate_num_tile_granular_data(component_type_t m_component_type, std::vector<unsigned> *m_tile_granumlar_data);

    // Calculate the number of active component (e.g., MAC, PE, and Chips)
    unsigned calculate_active_component(component_type_t m_component_type);
    // Print out the component of mapping table.
	void print();

private:
    std::vector<std::vector<unsigned>> mapping_table;   // The mapping table of neural layer.
    // Legacy GLB factors are retained for temporal accounting only.
    std::vector<unsigned> legacy_global_buffer_mapping;

};
#endif
