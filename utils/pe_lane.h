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

double calculate_time_based_mac_utilization(double busy_scalar_mac_cycles,
                                            double available_scalar_mac_cycles);

#endif
