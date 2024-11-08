#include <thread>

#include "npu.h"
#include "config.h"
#include "python_interface.h"

npu_t::npu_t() :
 num_processors(1),
 total_cycle(0),
 power_consumption(0.0),
 data_format(data_format_t::CONVOLUTION),
 compression_type(compression_type_t::DENSE),
 multi_chip(NULL),
 dram(NULL) {

}

npu_t::~npu_t() {

    // Free the memory for accelerator components
    for(unsigned i = 0; i < num_processors; i++) {
        delete pe_arrays[i];
        delete global_buffers[i];
    }

    delete multi_chip;
    delete dram;

	// Free the memory for the network.
	delete network;

	// Free the memory for mapping table.
    for(auto vector : mapping_tables) { delete vector;}
    mapping_tables.clear();
	//for(unsigned i = 0; i < mapping_tables.size(); i++) { delete mapping_tables[i]; }

    // Free the memory for scheduler.
    schedulers.clear();
    for(unsigned i = 0; i < schedulers.size(); i++) {
        delete schedulers[i];
    }

    for(unsigned i = 0; i < layer_stats.size(); i++) {
        delete layer_stats[i];
    }
    delete network_stats;

}

// Initialize the simulation environment.
void npu_t::init(const std::string m_accelerator_config, const std::string m_network_config, const std::string m_mapping_config) {

    /* Initialize DNN Accelerator */
    config_t accelerator_config;
    accelerator_config.parse(m_accelerator_config);

    // Initialize the components.
    for(unsigned i = 0 ; i < accelerator_config.sections.size(); i++) {
        section_config_t section_config = accelerator_config.sections[i];

        if(section_config.name == "accelerator") {
            section_config.get_setting("num_processors", &num_processors);
    
            // Initialize data format : Convolution or GEMM.
            std::string format_str;
            if(section_config.get_setting("data_format", &format_str)) {
                data_format = (data_format_t)get_type(data_format_str, format_str);
            }
            // Initialize compression type : Dense, CSR, CSC, SparseMap.
            std::string compression_str;
            if(section_config.get_setting("compression_type", &compression_str)) {
                compression_type = (compression_type_t)get_type(compression_type_str, compression_str);
            }
        }
        // Initialize PE array.
        else if(section_config.name == "adder_tree" || section_config.name == "ADDER_TREE") {
            pe_array_t *pe_array;
            for(unsigned i = 0; i < num_processors; i++) {
                pe_array = new adder_tree_t(section_config);
                pe_array->index = i;
                pe_arrays.emplace_back(pe_array);
            }
        }
        else if(section_config.name == "spatial_arch" || section_config.name == "SPATIAL_ARCH") {
            pe_array_t *pe_array;
            for(unsigned i = 0; i < num_processors; i++) {
                pe_array = new spatial_arch_t(section_config);
                pe_array->index = i;
                pe_arrays.emplace_back(pe_array);
            }
        }
        else if(section_config.name == "systolic_array" || section_config.name == "SYSTOLIC_ARRAY") {
            pe_array_t *pe_array;
            for(unsigned i = 0; i < num_processors; i++) {
                pe_array = new systolic_array_t(section_config);
                pe_array->index = i;
                pe_arrays.emplace_back(pe_array);
            }
        }
        // Initialize Global buffer.
        else if(section_config.name == "separate" || section_config.name == "SEPARATE") {
            global_buffer_t *global_buffer;
            for(unsigned i = 0; i < num_processors; i++) {
                global_buffer = new separate_buffer_t(section_config);
                global_buffer->index = i;
                global_buffers.emplace_back(global_buffer);
            }
        }
        else if(section_config.name == "shared" || section_config.name == " SHARED") {
            global_buffer_t *global_buffer;
            for(unsigned i = 0; i < num_processors; i++) {
                global_buffer = new shared_buffer_t(section_config);
                global_buffer->index = i;
                global_buffers.emplace_back(global_buffer);
            }
        }
        // Initialize Processors.
        else if(section_config.name == "multi_chip" || section_config.name == "MULTI_CHIP") {
            multi_chip = new multi_chip_t(section_config);
        }
        // Initialize off-chip memory 
        else if(section_config.name == "dram") {
            dram = new dram_t(section_config);
        }
        else {
            std::cerr << "Error: unknown accelerator component " << section_config.name << std::endl;
            exit(1);
        }
    }

    // Connect components
    connect();

    // Print out the stats of component
    print_accelerator_specification();

    /* Initialize the Neural network */
    std::cout << "# Initialize neural network model ..." << std::endl;

#ifdef PyTorch
    Py_Initialize();
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");

    pModule = PyImport_ImportModule("main");

    network = new network_t();
    network->init(pModule, m_network_config);
    network->init_weight(pModule);
#endif

    std::cout << "  Done!" << std::endl;

	/* Initialize the mapping table. */
    std::cout << "# Initialize the mapping table ..." << std::endl;
	mapping_table_t* mapping_table;
	config_t mapping_config;
	mapping_config.parse(m_mapping_config);
	mapping_tables.reserve(mapping_config.sections.size());
	for(unsigned i = 0; i < mapping_config.sections.size(); i++) {
		section_config_t section_config = mapping_config.sections[i];
		mapping_table = new mapping_table_t(section_config);
		mapping_tables.emplace_back(mapping_table);
	}
    std::cout << "  Done!" << std::endl;

    /* Initialize the scheduler */
    std::cout << "# Initialize the scheduler ..." << std::endl;
    scheduler_t* scheduler_;
    for(unsigned i = 0; i < mapping_tables.size(); i++) {
        scheduler_ = new scheduler_t(mapping_tables[i], pe_arrays[0]->pes[0]->get_mac_stationary_type(), pe_arrays[0]->pes[0]->get_parameter_order(),
                                                        pe_arrays[0]->get_stationary_type(), pe_arrays[0]->get_parameter_order(),
                                                        multi_chip->get_stationary_type(), multi_chip->get_parameter_order());
        scheduler_->compression_type = compression_type;
        schedulers.emplace_back(scheduler_);
    }
    std::cout << "  Done!" << std::endl;

    /* Initialize stats. */
    std::cout << "# Initialize the stat ..." << std::endl;
    stats_t* stats_;
    for(unsigned i = 0; i < mapping_tables.size(); i++) {
        stats_ = new stats_t();
        layer_stats.emplace_back(stats_);
    }
    network_stats = new stats_t();
    std::cout << "  Done!" << std::endl;
}

// Connect accelerator components.
void npu_t::connect() {
    for(unsigned i = 0; i < num_processors; i++) {
        // Connect PE array to Global buffer.
        pe_arrays[i]->connect(global_buffers[i]);

        // Connect global buffer to PE array and Multi Chip
        global_buffers[i]->connect(pe_arrays[i]);
        global_buffers[i]->connect(multi_chip);
    }

    // Connect Multi Chip to Global buffer and DRAM.
    multi_chip->connect(global_buffers);
    multi_chip->connect(dram);

    // Connect DRAM to Multi Chip
    dram->connect(multi_chip);
}

void npu_t::run(const std::string m_accelerator_config, const std::string m_network_config) {
	// Initialize the mapping table.
	std::cout << "# Run the network" << std::endl;

    // the number of iterations.
	unsigned num_iteration = 1;
    // Run the network
	for(unsigned iteration = 0; iteration < num_iteration; iteration++) {
		network->load_data(pModule, m_network_config, iteration);

        // Layer-wise simulation.
		for(unsigned index = 0; index < network->num_layers; index++) {
            network->layers[index]->input_data = index > 0 ? network->layers[index-1]->output_data : network->input_data;
            // Run on accelerator if the layer type is
            // Convolutional, Fully-connected, Max pooling, or Average pooling.
			if(network->layers[index]->layer_type == layer_name_t::CONVOLUTIONAL_LAYER ||
			   network->layers[index]->layer_type == layer_name_t::CONNECTED_LAYER) {
               //network->layers[index]->layer_type == layer_name_t::AVGPOOL_LAYER ||
               //network->layers[index]->layer_type == layer_name_t::MAXPOOL_LAYER) {

                layer = network->layers[index];
                dram->connect_layer(layer);

                scheduler = schedulers[index];
                
                print_network_configuration(index);
                reset();
                update_tile_size();

                std::string output_file_name = m_accelerator_config + "_" 
                                             + m_network_config + "_layer_" 
                                             + std::to_string(index) + "_power_trace.txt";
                std::ofstream output_file;
                output_file.open(output_file_name, std::ios::out);
                while(!is_idle()) {
                    power_consumption = 0.0;
                    // Process simulation in backward path.
                    execute();
                    // Transfer data from PE array to PE
                    transfer_data_to_pe();
                    // Transfer data from GLB to PE array
                    transfer_data_to_pe_array();
                    // Transfer data from Multi Chip to GLB
                    transfer_data_to_global_buffer();
                    // Transfer data from DRAM to Multi Chip
                    transfer_data_to_multi_chip();

                    // Request data from Multi Chip to DRAM.
                    request_to_dram();
                    // Request data from GLB to Multi Chip.
                    request_to_multi_chip();
                    // Request data from PE array to PE.
                    request_to_global_buffer();
                    // Request data from PEs to PE array.
                    request_to_pe_array();
                    power_measurement(output_file);
                }
                output_file.close();
                print_layerwise_results(m_accelerator_config, m_network_config, index);
			}
#ifdef PyTorch
            network->forward(pModule, iteration, index);
#endif
#ifdef FUNCTIONAL
            //network->layers[index]->forward();
#endif
		}
#ifdef PyTorch
        network->print_result(pModule);
#endif
        print_total_result(m_accelerator_config, m_network_config);
	}
#ifdef PyTorch
    Py_Finalize();
#endif
}

bool npu_t::is_idle() {
    // Check whether all PEs and PE array are idle or not
    for(unsigned i = 0; i < multi_chip->get_number_of_active_chips(); i++) {
        if(!pe_arrays[i]->is_idle()) {
            return false;
        }
    }
    // Check whether all processors are idle or not
    for(unsigned i = 0; i < multi_chip->get_number_of_active_chips(); i++) {
        if(!global_buffers[i]->is_idle()) {
            return false;
        }
    }
    // Check whether multi chip is idle or not
    if(!multi_chip->is_idle()) {
        return false;
    }
    // Check whether the off-chip is idle or not
    if(!dram->is_idle()) {
        return false;
    }
    return true;
}

// DNN execution (e.g., MAC and pooling) at PEs
void npu_t::execute() {
    for(unsigned i = 0; i < multi_chip->get_number_of_active_chips(); i++) {
        for(unsigned j = 0; j < pe_arrays[i]->get_number_of_active_pes(); j++) {
            if(pe_arrays[i]->pes[j]->is_exist_data()) {
                pe_arrays[i]->pes[j]->data_transfer_to_mac(scheduler);
            }
            else {
                update_power_consumption(pe_arrays[i]->pes[j]->get_static_power());
            }
        }
    }
}

// Transfer data from temporal buffer in PE array to PE (NoC)
void npu_t::transfer_data_to_pe() {
    for(unsigned i = 0; i < multi_chip->get_number_of_active_chips(); i++) {
        if(pe_arrays[i]->is_exist_request_at_pe() && pe_arrays[i]->is_exist_data()) {
            pe_arrays[i]->data_transfer(scheduler);
        }
    }
}


// Transfer data from global buffer to temporal buffer in PE array
void npu_t::transfer_data_to_pe_array() {
    for(unsigned i = 0; i < multi_chip->get_number_of_active_chips(); i++) {
        if(global_buffers[i]->is_exist_data() && pe_arrays[i]->is_exist_request_at_buffer()) {
            global_buffers[i]->data_transfer(scheduler);
        }
        else {
            update_power_consumption(global_buffers[i]->get_static_power());
        }
    }
}

// Transfer data from temporal buffer in the chip-level processors to the global buffer
void npu_t::transfer_data_to_global_buffer() {
    if(multi_chip->is_exist_request_at_global_buffer() && multi_chip->is_exist_data()) {
        multi_chip->data_transfer(scheduler);
    }
}

// Transfer data from the off-chip memory to temporal buffer in the chip-level processors
void npu_t::transfer_data_to_multi_chip() {
    if(multi_chip->is_exist_request_at_buffer()) {
        dram->data_transfer(scheduler);
    }
}

// Send a request signal from the chip-level processors to the off-chip memory
void npu_t::request_to_dram() {
    if(!multi_chip->is_exist_request_at_buffer() && multi_chip->is_exist_request_at_global_buffer() && !multi_chip->wait_data()) {
        multi_chip->request_data();
    }
}

// Send a request signal from the global buffer to chip-level processors
void npu_t::request_to_multi_chip() {
    for(unsigned i = 0; i < multi_chip->get_number_of_active_chips(); i++) {
        if(!global_buffers[i]->is_exist_data() && !global_buffers[i]->is_exist_request() && pe_arrays[i]->is_exist_request_at_buffer()) {
            global_buffers[i]->request_data();
        }
    }
}

// Send a request signal from PE array to the global buffer
void npu_t::request_to_global_buffer() {
    for(unsigned i = 0; i < multi_chip->get_number_of_active_chips(); i++) {
        if(!pe_arrays[i]->is_exist_request_at_buffer() && pe_arrays[i]->is_exist_request_at_pe() && !pe_arrays[i]->wait_data()) {
            pe_arrays[i]->request_data();
        }
    }
}

// Send a request signal from PEs to PE array
void npu_t::request_to_pe_array() {
    for(unsigned i = 0; i < multi_chip->get_number_of_active_chips(); i++) {
        for(unsigned j = 0; j < pe_arrays[i]->get_number_of_active_pes(); j++) {
            if(!pe_arrays[i]->pes[j]->is_exist_data() && ! pe_arrays[i]->pes[j]->is_exist_request()) {
                pe_arrays[i]->pes[j]->request_data();
            }
        }
    }
}

void npu_t::update_power_consumption(double m_power) {
    power_consumption += m_power;
}

// Print out the accelerator specification.
void npu_t::print_accelerator_specification() {
    pe_arrays[0]->print_specification();
    global_buffers[0]->print_specification();
    multi_chip->print_specification();
    dram->print_specification();
}

// Print out the network stats (e.g., tile size)
void npu_t::print_network_configuration(unsigned m_index) {
    std::cout << "The network configuration of #" << m_index << " layer" << std::endl;
    schedulers[m_index]->print_stats();
}

// Print out the simulation result
void npu_t::print_layerwise_results(const std::string m_accelerator_config, const std::string m_network_config, unsigned m_index) {
    std::cout << "The simulation result of #" << m_index << " layer" << std::endl;
    std::cout << std::endl;
    
    // Concatenate the name of output file.
    std::string output_file_name = m_accelerator_config + "_" + m_network_config + "_layer_" + std::to_string(m_index) + ".txt";
    std::cout << output_file_name << std::endl;

    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out);

    //layer_stats[m_index]->print_stats(output_file);
    schedulers[m_index]->print_stats(output_file);

    //schedulers
    layer_stats[m_index]->update_stats(pe_arrays, global_buffers, multi_chip, dram);
    layer_stats[m_index]->print_results(output_file);

#ifdef DRAMSIM3
    dram->print_result();
#endif

    output_file.close();
    network_stats->update_network_stats(layer_stats[m_index]);
}

void npu_t::print_total_result(const std::string m_accelerator_config, const std::string m_network_config) {
    std::cout << "The simulator result of " << m_network_config << std::endl;
    std::cout << std::endl;

    std::string output_file_name = m_accelerator_config + "_" + m_network_config + ".txt";
    std::cout << output_file_name << std::endl;

    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out);

    network_stats->print_results(output_file);

    output_file.close();
}

void npu_t::power_measurement(std::ofstream &m_output_file) {
    
    m_output_file << total_cycle << " cycles, " << power_consumption << " mW" << std::endl;
}

// Reset the simulation result and stats
void npu_t::reset() {
    dram->reset();
    multi_chip->reset();
    for(unsigned i = 0; i < multi_chip->get_number_of_active_chips(); i++) {
        global_buffers[i]->reset();
        pe_arrays[i]->reset();
    }
}

void npu_t::update_tile_size() {
    dram->update_tile_size(scheduler);
    multi_chip->update_tile_size(scheduler);

    for(unsigned i = 0; i < multi_chip->get_number_of_active_chips(); i++) {
        global_buffers[i]->update_tile_size(scheduler);
        pe_arrays[i]->update_tile_size(scheduler);
    }
}
