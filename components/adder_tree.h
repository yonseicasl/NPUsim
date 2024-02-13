#ifndef __ADDER_TREE_H__
#define __ADDER_TREE_H__

#include "pe_array.h"
#include "scheduler.h"

class adder_tree_t : public pe_array_t {

public:
	adder_tree_t(section_config_t m_section_config);
	~adder_tree_t();

	// Initialize the PE array.
	void init(section_config_t m_section_config);

    // Update tile size of PE array
    void update_tile_size(scheduler_t *m_scheduler);

    // Transfer data to the local buffers in PE.
    void data_transfer(scheduler_t *m_scheduler);

	// print out the configuration of PE array.
	void print_specification();
};

#endif
