#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <cstring>
#include <cstddef>
#include <string>
#include <vector>

#include "config.h"
#include "def.h"

// E20-3: the dense input footprint implied by a legacy-GLB spatial loop nest.
// `unique_elements` is the union of all overlapping convolution windows. The live pass fetches
// `replicated_elements` when every P/Q tile is charged independently. `working_set_elements` is
// the maximum input state a sliding implementation must retain in the GLB to realize that union.
struct input_halo_reuse_t {
    bool active;
    size_t unique_elements;
    size_t replicated_elements;
    size_t working_set_elements;

    input_halo_reuse_t() : active(false), unique_elements(0), replicated_elements(0),
                           working_set_elements(0) {}
};

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
    // times the legacy GLB temporal-repetition factors (K/B/P/Q/C/R/S), followed by
    // recomputation of the derived input H/W extents.
    std::vector<unsigned> calculate_total_parameter_size();
    // Per-datatype GLB temporal-repetition factor: the product of only those legacy
    // GLB loop factors that index the datatype (a repetition over a dimension a
    // tensor does not depend on revisits the SAME tile, which the on-chip buffers
    // retain -- off-chip traffic must not scale with it).
    std::vector<unsigned> datatype_repetitions();

    // E20-3: is a REDUCTION dimension (input channel, filter height, filter width) tiled at a level
    // ABOVE the PE array?
    //
    // This decides whether a partial sum may stay resident in the array. If every reduction step
    // Dense convolution input-halo contract for the legacy GLB P/Q loops. This computes the
    // exact union footprint, not a fractional "effective repetition": conv1/conv2 tile both P
    // and Q, and descriptor-tail rounding makes repetition-only scaling inexact. Filter factors
    // in the legacy row are intentionally left unsupported until their loop-order semantics are
    // validated; in that case the returned contract is inactive.
    input_halo_reuse_t input_halo_reuse() const;

    // for an output completes before the array moves on, the psum never leaves and charging a
    // write-back would be wrong. If a reduction dimension is split above the array, the array walks
    // ALL of the outer level's output tiles and then comes back -- it cannot hold them all, so each
    // revisit is a physical read-back and write-out through the global buffer.
    //
    // The model used to suppress that traffic whenever the array's output TILE equalled the global
    // buffer's, comparing one tile while the loop above it iterates many. On Eyeriss conv3 (input
    // channel split 64 ways at DRAM, 312 output tiles per step) the suppressed traffic is what the
    // measured chip spends most of its GLB bandwidth on.
    // Must a partial sum physically leave the PE array between reduction steps? True when a
    // reduction loop above the array has an output-bearing loop inside it, so the array walks
    // other output tiles and has to come back. A reduction above the array is NOT sufficient on
    // its own -- see the definition for why that weaker test misread every validated GEMM.
    bool psum_must_leave_array() const;

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
