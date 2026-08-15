#include "pe_lane.h"

#include <cmath>
#include <limits>

mac_lane_state_t calculate_mac_lane_state(unsigned number_of_macs, unsigned mac_width,
                                          unsigned active_scalar_lanes) {
    mac_lane_state_t state = {false, 0, 0, 0, 0.0};
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

double calculate_time_based_mac_utilization(double busy_scalar_mac_cycles,
                                            double available_scalar_mac_cycles) {
    if(busy_scalar_mac_cycles < 0.0 || available_scalar_mac_cycles <= 0.0) return 0.0;
    return busy_scalar_mac_cycles / available_scalar_mac_cycles;
}
