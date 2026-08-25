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
size_t checked_multiply_size(size_t lhs, size_t rhs, const char *context) {
    if(lhs != 0 && rhs > std::numeric_limits<size_t>::max()/lhs) {
        std::cerr << "Error: size overflow while calculating " << context << std::endl;
        exit(1);
    }
    return lhs*rhs;
}
size_t checked_add_size(size_t lhs, size_t rhs, const char *context) {
    if(rhs > std::numeric_limits<size_t>::max() - lhs) {
        std::cerr << "Error: size overflow while calculating " << context << std::endl;
        exit(1);
    }
    return lhs + rhs;
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

std::vector<unsigned> mapping_table_t::calculate_total_parameter_size() {
    std::vector<unsigned> parameters = calculate_parameter_size(component_type_t::DRAM);
    for(unsigned j = 0; j < 7; j++) {
        parameters[j] = checked_multiply(parameters[j], legacy_global_buffer_mapping[j],
                                         "total parameter size");
    }
    // INPUT_HEIGHT/WIDTH are derived extents, not multiplicative loop factors. Recompute them
    // after applying the legacy P/Q/R/S coverage so diagnostics see the same full convolution
    // footprint as the halo model.
    parameters[parameter_type_t::INPUT_HEIGHT] =
        checked_extent(parameters[parameter_type_t::OUTPUT_HEIGHT],
                       parameters[parameter_type_t::STRIDE],
                       parameters[parameter_type_t::FILTER_HEIGHT], "total input height");
    parameters[parameter_type_t::INPUT_WIDTH] =
        checked_extent(parameters[parameter_type_t::OUTPUT_WIDTH],
                       parameters[parameter_type_t::STRIDE],
                       parameters[parameter_type_t::FILTER_WIDTH], "total input width");
    return parameters;
}

std::vector<unsigned> mapping_table_t::datatype_repetitions() {
    const std::vector<unsigned> &glb = legacy_global_buffer_mapping;
    std::vector<unsigned> repetitions(data_type_t::NUM_DATA_TYPES, 1);
    // DR5 group convention (mirrors calculate_active_component): a legacy K entry
    // spans ALL groups and the row product is divided by GROUP. Hence K-bearing
    // factor sets (weight, output) divide by GROUP -- distinct chunks are (g, K/G)
    // -- while INPUT (no K) iterates once per group and multiplies by it.
    const unsigned group = glb[parameter_type_t::GROUP] ? glb[parameter_type_t::GROUP] : 1;
    repetitions[data_type_t::INPUT] = checked_multiply(
        checked_multiply(glb[parameter_type_t::BATCH_SIZE], glb[parameter_type_t::INPUT_CHANNEL], "input repetitions"),
        checked_multiply(glb[parameter_type_t::OUTPUT_HEIGHT],
                         checked_multiply(glb[parameter_type_t::OUTPUT_WIDTH], group, "input repetitions"),
                         "input repetitions"),
        "input repetitions");
    repetitions[data_type_t::WEIGHT] = checked_multiply(
        checked_multiply(glb[parameter_type_t::OUTPUT_CHANNEL], glb[parameter_type_t::INPUT_CHANNEL], "weight repetitions"),
        checked_multiply(glb[parameter_type_t::FILTER_HEIGHT], glb[parameter_type_t::FILTER_WIDTH], "weight repetitions"),
        "weight repetitions");
    repetitions[data_type_t::OUTPUT] = checked_multiply(
        checked_multiply(glb[parameter_type_t::OUTPUT_CHANNEL], glb[parameter_type_t::BATCH_SIZE], "output repetitions"),
        checked_multiply(glb[parameter_type_t::OUTPUT_HEIGHT], glb[parameter_type_t::OUTPUT_WIDTH], "output repetitions"),
        "output repetitions");
    if(group > 1) {
        if(repetitions[data_type_t::WEIGHT] % group || repetitions[data_type_t::OUTPUT] % group) {
            std::cerr << "Error: grouped legacy GLB mapping requires K-bearing factor products divisible by GROUP" << std::endl;
            exit(1);
        }
        repetitions[data_type_t::WEIGHT] /= group;
        repetitions[data_type_t::OUTPUT] /= group;
    }
    return repetitions;
}
input_halo_reuse_t mapping_table_t::input_halo_reuse() const {
    input_halo_reuse_t reuse;
    const std::vector<unsigned> &glb = legacy_global_buffer_mapping;

    // A filter loop in the legacy row partitions R/S rather than moving an adjacent output
    // window. Its input-union and loop-order semantics are different; do not silently apply the
    // P/Q contract to it.
    if(glb[parameter_type_t::FILTER_HEIGHT] != 1 ||
       glb[parameter_type_t::FILTER_WIDTH] != 1) return reuse;

    // Rebuild the cumulative DRAM-level tile without mutating the table. The legacy GLB row is
    // deliberately absent from calculate_parameter_size(DRAM), which is exactly the one live
    // pass whose off-chip counters are later repetition-scaled.
    std::vector<unsigned> base(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    for(unsigned level = 0; level <= component_type_t::DRAM; ++level) {
        for(unsigned dim = 0; dim < parameter_type_t::NUM_PARAMETER_TYPES; ++dim) {
            base[dim] = checked_multiply(base[dim], mapping_table[level][dim],
                                         "input halo base parameter");
        }
    }
    base[parameter_type_t::STRIDE] =
        mapping_table[component_type_t::DRAM][parameter_type_t::STRIDE];
    const unsigned base_h = checked_extent(base[parameter_type_t::OUTPUT_HEIGHT],
                                           base[parameter_type_t::STRIDE],
                                           base[parameter_type_t::FILTER_HEIGHT],
                                           "input halo base height");
    const unsigned base_w = checked_extent(base[parameter_type_t::OUTPUT_WIDTH],
                                           base[parameter_type_t::STRIDE],
                                           base[parameter_type_t::FILTER_WIDTH],
                                           "input halo base width");
    const unsigned full_p = checked_multiply(base[parameter_type_t::OUTPUT_HEIGHT],
                                             glb[parameter_type_t::OUTPUT_HEIGHT],
                                             "input halo full output height");
    const unsigned full_q = checked_multiply(base[parameter_type_t::OUTPUT_WIDTH],
                                             glb[parameter_type_t::OUTPUT_WIDTH],
                                             "input halo full output width");
    const unsigned full_h = checked_extent(full_p, base[parameter_type_t::STRIDE],
                                           base[parameter_type_t::FILTER_HEIGHT],
                                           "input halo full height");
    const unsigned full_w = checked_extent(full_q, base[parameter_type_t::STRIDE],
                                           base[parameter_type_t::FILTER_WIDTH],
                                           "input halo full width");

    size_t nonspatial = checked_multiply_size(base[parameter_type_t::BATCH_SIZE],
                                              base[parameter_type_t::INPUT_CHANNEL],
                                              "input halo nonspatial tile");
    const size_t legacy_nonspatial = checked_multiply_size(
        checked_multiply_size(glb[parameter_type_t::BATCH_SIZE],
                              glb[parameter_type_t::INPUT_CHANNEL],
                              "input halo legacy nonspatial"),
        glb[parameter_type_t::GROUP], "input halo legacy groups");
    reuse.unique_elements = checked_multiply_size(
        checked_multiply_size(nonspatial, legacy_nonspatial, "input halo unique elements"),
        checked_multiply_size(full_h, full_w, "input halo unique spatial footprint"),
        "input halo unique elements");
    reuse.replicated_elements = checked_multiply_size(
        checked_multiply_size(nonspatial, legacy_nonspatial,
                              "input halo replicated elements"),
        checked_multiply_size(
            checked_multiply_size(base_h, base_w, "input halo base spatial footprint"),
            checked_multiply_size(glb[parameter_type_t::OUTPUT_HEIGHT],
                                  glb[parameter_type_t::OUTPUT_WIDTH],
                                  "input halo spatial repetitions"),
            "input halo replicated spatial footprint"),
        "input halo replicated elements");

    // Q sliding needs one base tile. Crossing a legacy-P boundary additionally retains the
    // overlapping rows across the full Q width; only the new rows of the current base tile add
    // storage. This is the ring-buffer working set, not the whole layer input tensor.
    const unsigned p_shift = checked_multiply(base[parameter_type_t::OUTPUT_HEIGHT],
                                              base[parameter_type_t::STRIDE],
                                              "input halo P shift");
    const unsigned overlap_h = glb[parameter_type_t::OUTPUT_HEIGHT] > 1 && base_h > p_shift
                             ? base_h - p_shift : 0;
    const size_t working_spatial = overlap_h
        ? checked_add_size(
              checked_multiply_size(overlap_h, full_w, "input halo retained rows"),
              checked_multiply_size(base_h - overlap_h, base_w, "input halo new tile rows"),
              "input halo working rows")
        : checked_multiply_size(base_h, base_w, "input halo working tile");
    reuse.working_set_elements = checked_multiply_size(nonspatial, working_spatial,
                                                        "input halo working set");
    reuse.active = reuse.unique_elements < reuse.replicated_elements;
    return reuse;
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

// E20-3: see the header. A reduction factor above the array means a psum is revisited after other
// tiles have passed through, so it cannot be retained.
bool mapping_table_t::psum_must_leave_array() const {
    static const parameter_type_t REDUCTION[] = { parameter_type_t::INPUT_CHANNEL,
                                                  parameter_type_t::FILTER_HEIGHT,
                                                  parameter_type_t::FILTER_WIDTH };
    static const parameter_type_t OUTPUT_BEARING[] = { parameter_type_t::OUTPUT_CHANNEL,
                                                       parameter_type_t::BATCH_SIZE,
                                                       parameter_type_t::OUTPUT_HEIGHT,
                                                       parameter_type_t::OUTPUT_WIDTH,
                                                       parameter_type_t::GROUP };
    // Levels above the PE array, INNERMOST FIRST. The nest runs outermost (DRAM) to innermost,
    // so walking it this way lets one pass answer "is there an output loop inside this one".
    static const component_type_t ABOVE_ARRAY[] = { component_type_t::GLOBAL_BUFFER,
                                                    component_type_t::CHIPS_X,
                                                    component_type_t::CHIPS_Y,
                                                    component_type_t::DRAM };
    bool output_loop_inside = false;
    for(unsigned level = 0; level < sizeof(ABOVE_ARRAY)/sizeof(ABOVE_ARRAY[0]); ++level) {
        // The GLB level keeps its factors in the legacy row, which is deliberately absent from
        // mapping_table[GLOBAL_BUFFER] (see calculate_total_parameter_size). Reading only the
        // table made every GLB-level factor invisible here: a reduction mapped AT the GLB --
        // every validated GEMM does exactly that -- read as "no reduction above the array".
        // Multiplying the two rows is correct whichever one carries the factor.
        const bool glb = ABOVE_ARRAY[level] == component_type_t::GLOBAL_BUFFER;
        bool has_reduction = false, has_output = false;
        for(unsigned dim = 0; dim < sizeof(REDUCTION)/sizeof(REDUCTION[0]); ++dim) {
            const unsigned factor = mapping_table[ABOVE_ARRAY[level]][REDUCTION[dim]]
                                  * (glb ? legacy_global_buffer_mapping[REDUCTION[dim]] : 1);
            if(factor > 1) has_reduction = true;
        }
        for(unsigned dim = 0; dim < sizeof(OUTPUT_BEARING)/sizeof(OUTPUT_BEARING[0]); ++dim) {
            const unsigned factor = mapping_table[ABOVE_ARRAY[level]][OUTPUT_BEARING[dim]]
                                  * (glb ? legacy_global_buffer_mapping[OUTPUT_BEARING[dim]] : 1);
            if(factor > 1) has_output = true;
        }
        // A reduction loop alone does NOT push the partial sum out. What pushes it out is the
        // array visiting a DIFFERENT output tile between two steps of that reduction, and only
        // an output-bearing loop INSIDE the reduction (or beside it at the same level, where the
        // within-level order is not modeled -- taken conservatively) can do that. With no output
        // loop under it the array simply keeps accumulating into the tile it already holds.
        if(has_reduction && (output_loop_inside || has_output)) return true;
        if(has_output) output_loop_inside = true;
    }
    return false;
}
