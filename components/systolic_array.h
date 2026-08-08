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

};

#endif
