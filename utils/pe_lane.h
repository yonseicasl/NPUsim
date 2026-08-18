#ifndef __PE_LANE_H__
#define __PE_LANE_H__

#include <cstddef>

struct mac_lane_state_t {
    bool valid;
    size_t scalar_capacity;
    unsigned active_accumulator_units;
    unsigned final_accumulator_lanes;
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

double calculate_time_based_mac_utilization(double busy_scalar_mac_cycles,
                                            double available_scalar_mac_cycles);

#endif
