#ifndef __SPATIAL_ARCH_H__
#define __SPATIAL_ARCH_H__

#include "pe_array.h"
#include "scheduler.h"

class spatial_arch_t : public pe_array_t {

public:
    spatial_arch_t(section_config_t m_section_config);
    ~spatial_arch_t();
    
    // Initialize the PU.
    void init(section_config_t m_section_config);

    // Update tile size of PE array
    void update_tile_size(scheduler_t *m_scheduler);

    // Transfer data to the local buffers in PE
    void data_transfer(scheduler_t *m_scheduler);

    // Print out the configuration of PU.
    void print_specification();
};

#endif
