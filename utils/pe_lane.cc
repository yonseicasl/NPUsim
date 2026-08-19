#include "pe_lane.h"

#include <cmath>
#include <limits>

mac_lane_state_t calculate_mac_lane_state(unsigned number_of_macs, unsigned mac_width,
                                          unsigned active_scalar_lanes) {
    mac_lane_state_t state = {false, 0, mac_width, 0, 0, 0.0};
    if(number_of_macs == 0 || mac_width == 0 ||
       number_of_macs > std::numeric_limits<size_t>::max() / mac_width) {
        return state;
    }

    state.scalar_capacity = static_cast<size_t>(number_of_macs) * mac_width;
    if(active_scalar_lanes == 0 || active_scalar_lanes > state.scalar_capacity) return state;

    state.valid = true;
    state.active_accumulator_units = (active_scalar_lanes + mac_width - 1) / mac_width;
    state.final_accumulator_lanes = active_scalar_lanes -
        (state.active_accumulator_units - 1) * mac_width;
    state.utilization = static_cast<double>(active_scalar_lanes) /
                        static_cast<double>(state.scalar_capacity);
    return state;
}

double accumulate_issue_cycles(size_t issue_steps, double unit_cycle) {
    return static_cast<double>(issue_steps) * unit_cycle;
}

double lane_reduction_fill_cycles(unsigned mac_width, double unit_cycle) {
    if(mac_width <= 1) return 0.0;
    return std::ceil(std::log2(static_cast<double>(mac_width))) * unit_cycle;
}

double lane_reduction_energy(unsigned mac_width, double unit_energy) {
    if(mac_width <= 1) return 0.0;
    return static_cast<double>(mac_width - 1) * unit_energy;
}

double lane_reduction_fill_cycles(const mac_lane_state_t &lanes, double unit_cycle) {
    if(!lanes.valid || lanes.active_accumulator_units == 0) return 0.0;
    // The live accumulators reduce CONCURRENTLY, so the fill is the depth of the deepest
    // live tree. Only the last unit can be partial, so a full-width tree exists as soon as
    // more than one unit is active.
    const unsigned leaves = (lanes.active_accumulator_units > 1)
        ? lanes.lane_width : lanes.final_accumulator_lanes;
    if(leaves <= 1) return 0.0;
    return std::ceil(std::log2(static_cast<double>(leaves))) * unit_cycle;
}

double lane_reduction_energy(const mac_lane_state_t &lanes, double unit_energy) {
    if(!lanes.valid || lanes.active_accumulator_units == 0) return 0.0;
    // Energy is total WORK across every live tree, not the depth of one: an N-leaf binary
    // adder tree performs N-1 additions.
    const unsigned full_units = lanes.active_accumulator_units - 1;
    double additions = 0.0;
    if(lanes.lane_width > 1) {
        additions += static_cast<double>(full_units)*static_cast<double>(lanes.lane_width - 1);
    }
    if(lanes.final_accumulator_lanes > 1) {
        additions += static_cast<double>(lanes.final_accumulator_lanes - 1);
    }
    return additions * unit_energy;
}

double calculate_time_based_mac_utilization(double busy_scalar_mac_cycles,
                                            double available_scalar_mac_cycles) {
    if(busy_scalar_mac_cycles < 0.0 || available_scalar_mac_cycles <= 0.0) return 0.0;
    return busy_scalar_mac_cycles / available_scalar_mac_cycles;
}
