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

	// Flatten DNN data.
	void flatten(unsigned m_channel, unsigned m_height, unsigned m_width,
				unsigned m_weight_height, unsigned m_weight_width, unsigned m_stride, 
				unsigned m_padding_height, unsigned m_padding_width,
                data_t *m_data, data_t *m_workspace);

    // Print out the configuration of PU.
    void print_specification();

    data_t *workspace_input;                        // Input data for GEMM used in systolic array
    data_t *workspace_weight;                       // Weight for GEMM used in systolic array
    data_t *workspace_output;                       // Output data for GEMM used in systolic array
};

#endif
