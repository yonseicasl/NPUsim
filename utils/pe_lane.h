#ifndef __PE_LANE_H__
#define __PE_LANE_H__

#include <cstddef>

struct mac_lane_state_t {
    bool valid;
    size_t scalar_capacity;
    unsigned lane_width;                  // structural lanes per accumulator (mac_width)
    unsigned active_accumulator_units;
    unsigned final_accumulator_lanes;     // lanes feeding the LAST active accumulator
    double utilization;
};

mac_lane_state_t calculate_mac_lane_state(unsigned number_of_macs, unsigned mac_width,
                                          unsigned active_scalar_lanes);
double accumulate_issue_cycles(size_t issue_steps, double unit_cycle);

// A3/PE1: latency of the lane->accumulator reduction. `mac_width` scalar FMA lanes
// feed one accumulator through a ceil(log2(width))-stage adder tree; the tree is
// pipelined, so it adds fill latency per issue step without changing throughput or
// per-FMA energy. Width 1 has no reduction.
double lane_reduction_fill_cycles(unsigned mac_width, double unit_cycle);

// P4-2/PE2: energy of the same lane->accumulator reduction. Unlike the cycle fill
// (pipeline DEPTH, ceil(log2(width)) stages), energy is proportional to the total
// WORK done -- a standard N-leaf binary adder tree performs N-1 additions
// regardless of depth (mirrors adder_tree_reduction_cost()'s num_additions
// convention). Defaults to a no-op unless a config calibrates unit_energy > 0.
// Width 1 has no reduction.
double lane_reduction_energy(unsigned mac_width, double unit_energy);

// L9/PE1-PE2: the same reduction charged against the ACTIVE lane state instead of the
// structural width. A mapping activates `active_accumulator_units` accumulators; all but
// the last are fed by a full `mac_width` lanes and the last by `final_accumulator_lanes`
// (see calculate_mac_lane_state()). Charging the structural width was wrong in both
// directions: it billed a full-width tree when only a few lanes were live, and it billed
// only ONE tree's work when several accumulators reduced concurrently.
//
//   LATENCY  is the depth of the DEEPEST live tree, because the accumulators reduce in
//            parallel -- ceil(log2(mac_width)) as soon as any unit is full, otherwise
//            ceil(log2(final_accumulator_lanes)).
//   ENERGY   is the total work over ALL live trees: an N-leaf tree performs N-1 additions,
//            so (units-1)*(mac_width-1) + (final_accumulator_lanes-1).
double lane_reduction_fill_cycles(const mac_lane_state_t &lanes, double unit_cycle);
double lane_reduction_energy(const mac_lane_state_t &lanes, double unit_energy);

double calculate_time_based_mac_utilization(double busy_scalar_mac_cycles,
                                            double available_scalar_mac_cycles);

#endif
