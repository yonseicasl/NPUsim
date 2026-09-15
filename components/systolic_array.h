#ifndef __SYSTOLIC_ARRAY_H__
#define __SYSTOLIC_ARRAY_H__

#include "pe_array.h"
#include "interconnect_timing.h"
#include "scheduler.h"

class systolic_array_t : public pe_array_t {

public:
    systolic_array_t(section_config_t m_section_config);
    ~systolic_array_t();
    
    // Initialize the PU.
    void init(section_config_t m_section_config);

    // Update tile size of PE array
    void update_tile_size(scheduler_t *m_scheduler);

    // Transfer data to the local buffers
    void data_transfer(scheduler_t *m_scheduler);

    // Print out the configuration of PU.
    void print_specification();

protected:
    // P4-3/SY2: the array is a 2D grid on the drain-out (output write-back) direction
    // just as on the load direction (see data_transfer()'s override), regardless of
    // the configured noc label.
    noc_type_t writeback_noc_type() const { return noc_type_t::MESH; }
    // SY2/L9: the accumulation pipeline drains down the active columns before a different
    // weight residency can take effect. This is a WEIGHT-STATIONARY phenomenon: stationary
    // weights are swapped, and the psums must flush before the new weights compute.
    //
    // OS1 (2026-09-15): under OUTPUT_STATIONARY the psum stays resident in the PE and the
    // reduction STREAMS through it -- weights (and inputs) flow in one operand per step with
    // no per-step flush. The pipeline fills once when an output tile begins accumulating and
    // drains once when it is evicted (a per-output-tile event, negligible next to the
    // reduction it overlaps), NOT once per reduction-step weight distribution. Charging this
    // drain per weight distribution (the caller's event) is the same class of per-step
    // over-count the WS mapping fix removed, and validation/phase2/os_size_sweep measured it
    // directly: comp + this analytical fold overshoots the 32x32 OS RTL golden by ~675% MAPE
    // (documented there as "a loose upper bound"), while the Computation cycle alone is exact.
    // So do not charge it in OS. The genuine OS overhead (fixed per-layer setup + exposed DRAM,
    // os_size_sweep fit: setup~=2270, ~72 B/cycle) is a separate, still-uncalibrated term and
    // is deliberately NOT synthesized here.
    double weight_fold_bubble_cycles() const {
        if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) return 0.0;
        return systolic_pipeline_cost(num_active_pe_y, num_active_pe_x).drain_hops*noc_cycle;
    }
};

#endif
