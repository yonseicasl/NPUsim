#ifndef __SYSTOLIC_ARRAY_H__
#define __SYSTOLIC_ARRAY_H__

#include "pe_array.h"
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
};

#endif
