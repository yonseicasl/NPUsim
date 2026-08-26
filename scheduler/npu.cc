#include <thread>

#include "npu.h"
#include "energy_units.h"
#include "config.h"
#include "datatype.h"

npu_t::npu_t() :
 num_processors(1),
 num_pes(1),
 compression_type(compression_type_t::DENSE),
 num_skipped_timing_layers(0),
 multi_chip(NULL),
 dram(NULL),
 network(NULL),
 layer(NULL),
 scheduler(NULL),
 network_stats(NULL) {

}

npu_t::~npu_t() {

    // Free the memory for accelerator components
    for(auto pe_array : pe_arrays) { delete pe_array; }
    for(auto global_buffer : global_buffers) { delete global_buffer; }
    for(auto sfu : sfus) { delete sfu; }

    delete multi_chip;
    delete dram;

	// Free the memory for the network.
	delete network;

	// Free the memory for mapping table.
    for(auto mapping_table : mapping_tables) { delete mapping_table; }

    // Free the memory for scheduler.
    for(auto scheduler_ : schedulers) { delete scheduler_; }
    for(auto stats : layer_stats) { delete stats; }
    for(auto stats : sfu_layer_stats) { delete stats; }
    delete network_stats;

}

// Initialize the simulation environment.
void npu_t::init(const std::string m_accelerator_config, const std::string m_network_config, const std::string m_mapping_config) {

    /* Initialize DNN Accelerator */
    config_t accelerator_config;
    accelerator_config.parse(m_accelerator_config);

    // Read the accelerator-wide chip count before creating per-chip components.
    unsigned accelerator_sections = 0;
    for(unsigned i = 0; i < accelerator_config.sections.size(); i++) {
        section_config_t section_config = accelerator_config.sections[i];
        std::string section_name = section_config.name;
        lowercase(section_name);
        if(section_name != "accelerator") continue;

        accelerator_sections++;
        runtime_datatypes().configure(section_config);
        // E7: what the energy numbers mean (absolute pJ vs normalized) and where they came
        // from. Printed with the energy summary so a relative fixture cannot be read as an
        // absolute one.
        energy_units().configure(section_config);
        std::cout << "# Energy unit: " << energy_units().describe() << std::endl;
        std::cout << "# Runtime formats: input=" << runtime_datatypes().describe(data_type_t::INPUT)
                  << " weight=" << runtime_datatypes().describe(data_type_t::WEIGHT)
                  << " output=" << runtime_datatypes().describe(data_type_t::OUTPUT)
                  << " accumulator=" << runtime_datatypes().accumulator_format().name << std::endl;
        if(!section_config.get_setting("num_chips", &num_processors) &&
           !section_config.get_setting("num_processors", &num_processors)) {
            std::cerr << "Error: [accelerator] requires num_chips" << std::endl;
            exit(1);
        }
    }
    // E8: reject unusable energy unit costs before any component reads them, so a bad config
    // fails fast instead of producing negative or NaN energy that looks like a result.
    const std::string energy_error = validate_energy_settings(accelerator_config);
    if(!energy_error.empty()) {
        std::cerr << "Error: invalid energy unit cost: " << energy_error << std::endl;
        exit(1);
    }
    // RE5: derive each component's DECLARATION state now, so the energy breakdown can tell a
    // modeled zero from a missing cost from a component the layer simply never touched.
    energy_cost_schema().configure(accelerator_config);
    if(accelerator_sections != 1 || num_processors == 0) {
        std::cerr << "Error: accelerator config requires exactly one non-zero [accelerator] section"
                  << std::endl;
        exit(1);
    }

    // Initialize the components.
    for(unsigned i = 0 ; i < accelerator_config.sections.size(); i++) {
        section_config_t section_config = accelerator_config.sections[i];

        if(section_config.name == "accelerator") {
            section_config.get_setting("num_pes", &num_pes);
    
            // Initialize compression type : Dense, CSR, CSC, SparseMap.
            std::string compression_str;
            if(section_config.get_setting("compression_type", &compression_str)) {
                compression_type = (compression_type_t)get_type(compression_type_str, compression_str);
            }
            if(compression_type != compression_type_t::DENSE) {
                std::cerr << "Error: sparse PE execution is not implemented; use compression_type=dense" << std::endl;
                exit(1);
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
        else if(section_config.name == "shared" || section_config.name == "SHARED") {
            global_buffer_t *global_buffer;
            for(unsigned i = 0; i < num_processors; i++) {
                global_buffer = new shared_buffer_t(section_config);
                global_buffer->index = i;
                global_buffers.emplace_back(global_buffer);
            }
        }
        // Initialize Processors.
        else if(section_config.name == "multi_chip" || section_config.name == "MULTI_CHIP") {
            if(multi_chip != NULL) {
                std::cerr << "Error: duplicate [multi_chip] section" << std::endl;
                exit(1);
            }
            multi_chip = new multi_chip_t(section_config);
        }
        // Initialize off-chip memory
        else if(section_config.name == "dram") {
            if(dram != NULL) {
                std::cerr << "Error: duplicate [dram] section" << std::endl;
                exit(1);
            }
            dram = new dram_t(section_config);
        }
        // Initialize the per-chip Special Function Unit (opt-in; plan/plan_sfu.md).
        else if(section_config.name == "sfu" || section_config.name == "SFU") {
            if(!sfus.empty()) {
                std::cerr << "Error: duplicate [sfu] section" << std::endl;
                exit(1);
            }
            for(unsigned i = 0; i < num_processors; i++) {
                sfu_t *sfu = new sfu_t(section_config);
                sfu->index = i;
                sfus.emplace_back(sfu);
            }
        }
        else {
            std::cerr << "Error: unknown accelerator component " << section_config.name << std::endl;
            exit(1);
        }
    }
    validate_accelerator_components();


    // Connect components
    connect();

    // Print out the stats of component
    print_accelerator_specification();

    /* Initialize the Neural network */
    std::cout << "# Initialize neural network model ..." << std::endl;
	network = new nebula::convolutional_t();
	network->init(m_network_config);
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
    if(mapping_tables.empty()) {
        std::cerr << "Error: mapping config contains no layers" << std::endl;
        exit(1);
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
        layer_stats[i]->update_tile_size(schedulers[i]);
    }
    network_stats = new stats_t();
    std::cout << "  Done!" << std::endl;
}
void npu_t::validate_accelerator_components() {
    if(pe_arrays.size() != num_processors || global_buffers.size() != num_processors) {
        std::cerr << "Error: expected " << num_processors
                  << " PE arrays and global buffers, but found "
                  << pe_arrays.size() << " and " << global_buffers.size() << std::endl;
        exit(1);
    }
    if(multi_chip == NULL || dram == NULL) {
        std::cerr << "Error: accelerator config requires one [multi_chip] and one [dram] section" << std::endl;
        exit(1);
    }
    if(multi_chip->get_number_of_chips() != num_processors) {
        std::cerr << "Error: num_chips=" << num_processors
                  << " does not match [multi_chip] height*width="
                  << multi_chip->get_number_of_chips() << std::endl;
        exit(1);
    }
    for(unsigned i = 0; i < pe_arrays.size(); i++) {
        if(pe_arrays[i] == NULL || pe_arrays[i]->pes.empty() || global_buffers[i] == NULL) {
            std::cerr << "Error: chip " << i << " has an incomplete PE-array/global-buffer hierarchy" << std::endl;
            exit(1);
        }
    }
}

void npu_t::validate_active_components() {
    const size_t active_chips = static_cast<size_t>(scheduler->num_active_chips_x) * scheduler->num_active_chips_y;
    const size_t active_pes = static_cast<size_t>(scheduler->num_active_pe_x) * scheduler->num_active_pe_y;
    if(active_chips == 0 || active_chips > pe_arrays.size() || active_chips > global_buffers.size()) {
        std::cerr << "Error: mapping activates " << active_chips
                  << " chips, but the accelerator provides " << num_processors << std::endl;
        exit(1);
    }
    for(size_t i = 0; i < active_chips; i++) {
        if(active_pes == 0 || active_pes > pe_arrays[i]->get_number_of_pes()) {
            std::cerr << "Error: mapping activates " << active_pes
                      << " PEs on chip " << i << ", but only "
                      << pe_arrays[i]->get_number_of_pes() << " are available" << std::endl;
            exit(1);
        }
    }
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

    // SFU: validate every layer's activation BEFORE simulating anything. An unsupported
    // operation must fail fast at layer-start/config time (plan_sfu.md) -- aborting at
    // the END of that layer's simulation would discard every completed layer first.
    if(!sfus.empty()) {
        for(unsigned index = 0; index < network->num_layers; index++) {
            if(network->layers[index]->layer_type != nebula::CONVOLUTIONAL_LAYER &&
               network->layers[index]->layer_type != nebula::CONNECTED_LAYER) continue;
            const unsigned type = static_cast<unsigned>(network->layers[index]->activation_type);
            if(type == nebula::UNDEFINED_ACTIVATION ||
               type >= nebula::NUM_ACTIVATION_TYPES) continue;   // identity: linear bypass
            sfu_op_t op;
            if(!sfu_t::op_from_name(nebula::activation_type_str[type], &op)) {
                std::cerr << "Error: network layer " << index << " activation '"
                          << nebula::activation_type_str[type]
                          << "' is not supported by the SFU; unsupported operations fail"
                          << " fast before simulation instead of silently falling back to"
                          << " ReLU" << std::endl;
                exit(1);
            }
        }
    }

    // the number of iterations.
	unsigned num_iteration = 1;
    // Run the network
	for(unsigned iteration = 0; iteration < num_iteration; iteration++) {
		network->load_data(iteration);
        num_skipped_timing_layers = 0;

        // Layer-wise simulation.
		for(unsigned index = 0; index < network->num_layers; index++) {
            network->layers[index]->input_data = index > 0 ? network->layers[index-1]->output_data : network->input_data;
            // The current accelerator timing path supports convolution and fully-connected layers.
			if(network->layers[index]->layer_type == nebula::CONVOLUTIONAL_LAYER ||
			   network->layers[index]->layer_type == nebula::CONNECTED_LAYER) {
                if(index >= schedulers.size()) {
                    std::cerr << "Error: no mapping section for executable network layer "
                              << index << std::endl;
                    exit(1);
                }
                if(network->layers[index]->layer_type == nebula::CONVOLUTIONAL_LAYER) {
                    schedulers[index]->layer_name = layer_name_t::CONVOLUTIONAL_LAYER;
                } else if(network->layers[index]->layer_type == nebula::CONNECTED_LAYER) {
                    schedulers[index]->layer_name = layer_name_t::CONNECTED_LAYER;
                }

                layer = network->layers[index];
                dram->connect_layer(layer);

                scheduler = schedulers[index];
                const unsigned global_buffer_repetitions = scheduler->mapping_table->calculate_active_component(component_type_t::GLOBAL_BUFFER);

                // Spatial-padding guard: mapping factor products may exceed the layer
                // dimensions (silicon-style padding, e.g. Eyeriss computes 55 output
                // columns as 7x8 = 56), but that padding -- or an accidental
                // under-coverage -- must be visible, not silent.
                {
                    const std::vector<unsigned> mapped =
                        scheduler->mapping_table->calculate_total_parameter_size();
                    // Connected layers fold the input volume into a flat C dimension.
                    const unsigned layer_c =
                        (scheduler->layer_name == layer_name_t::CONNECTED_LAYER)
                            ? layer->input_size : layer->input_channel;
                    const struct { const char *name; unsigned mapped_size; unsigned layer_size; } dims[] = {
                        {"K", mapped[parameter_type_t::OUTPUT_CHANNEL], layer->output_channel},
                        {"P", mapped[parameter_type_t::OUTPUT_HEIGHT], layer->output_height},
                        {"Q", mapped[parameter_type_t::OUTPUT_WIDTH], layer->output_width},
                        {"C", mapped[parameter_type_t::INPUT_CHANNEL], layer_c},
                    };
                    for(const auto &dim : dims) {
                        if(dim.layer_size == 0) continue;
                        if(dim.mapped_size > dim.layer_size) {
                            std::cerr << "Warning: layer " << index << " mapping pads " << dim.name
                                      << " from " << dim.layer_size << " to " << dim.mapped_size
                                      << " (padded work is charged as compute)" << std::endl;
                        } else if(dim.mapped_size < dim.layer_size) {
                            std::cerr << "Warning: layer " << index << " mapping covers only "
                                      << dim.mapped_size << " of " << dim.layer_size << " in " << dim.name
                                      << " (layer is partially simulated)" << std::endl;
                        }
                    }
                }

                print_network_configuration(index);
                reset();
                update_tile_size();

                while(!is_idle()) {
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
                }
                // E20-3b: store the array's final resident output tile, then DR6's final
                // multi-chip -> DRAM store. Innermost boundary first, so the tile is written out
                // of the array before the level above flushes it off-chip.
                for(unsigned i = 0; i < pe_arrays.size(); i++) pe_arrays[i]->flush_psum_writeback();
                // DR6: store the layer's final resident output tile.
                multi_chip->flush_output_writeback();
                layer_stats[index]->update_stats(pe_arrays, global_buffers, multi_chip, dram);
                const input_halo_reuse_t input_halo =
                    scheduler->mapping_table->input_halo_reuse();
                const bool halo_capacity_sufficient = !global_buffers.empty() &&
                    global_buffers[0]->can_retain_input_halo(input_halo.working_set_elements);
                layer_stats[index]->scale_serial_repetitions(global_buffer_repetitions,
                    scheduler->mapping_table->datatype_repetitions(), input_halo,
                    halo_capacity_sufficient);
                // SFU: the fused activation fires once per valid output element, on the
                // FINAL (repetition-scaled) timeline -- reduction tiling and GLB
                // repetitions must not multiply the activation count.
                apply_fused_sfu_activation(index);
                print_layerwise_results(m_accelerator_config, m_network_config, index);
                //dram->disconnect_layer();
			}
            // Phase-7 (plan_sfu.md): a standalone softmax layer runs on the SFU's
            // multi-pass microprogram when an [sfu] section is configured. It needs no
            // mapping section -- the SFU is the only modeled resource.
            else if(network->layers[index]->layer_type == nebula::SOFTMAX_LAYER && !sfus.empty()) {
                run_standalone_softmax(index, m_accelerator_config, m_network_config);
            }
            else {
                num_skipped_timing_layers++;
                std::cerr << "Warning: network layer " << index
                          << " is excluded from accelerator timing (only convolution/connected are supported)" << std::endl;
            }
#ifdef FUNCTIONAL
            network->layers[index]->forward();
#endif
		}
        print_total_result(m_accelerator_config, m_network_config);
	}
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

// Print out the accelerator specification.
void npu_t::print_accelerator_specification() {
    pe_arrays[0]->print_specification();
    global_buffers[0]->print_specification();
    multi_chip->print_specification();
    dram->print_specification();
    if(!sfus.empty()) { sfus[0]->print_specification(); }
}

// Print out the network stats (e.g., tile size)
void npu_t::print_network_configuration(unsigned m_index) {
    std::cout << "The network configuration of #" << m_index << " layer" << std::endl;
    layer_stats[m_index]->print_stats();
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
    layer_stats[m_index]->print_stats(output_file);

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

    // L11/P4-14: hand the rollup its timing scope before printing, so the scope line appears
    // next to the latency it qualifies rather than only as a trailing warning.
    network_stats->excluded_timing_layers = num_skipped_timing_layers;
    network_stats->print_results(output_file);
    if(num_skipped_timing_layers > 0) {
        const std::string warning = "WARNING: partial timing result; " +
            std::to_string(num_skipped_timing_layers) +
            " non-convolution/connected layers were excluded.";
        std::cerr << warning << std::endl;
        output_file << warning << std::endl << std::endl;
    }

    output_file.close();
}

// Reset the simulation result and stats
void npu_t::reset() {
    dram->reset();
    multi_chip->reset();
    for(unsigned i = 0; i < multi_chip->get_number_of_active_chips(); i++) {
        global_buffers[i]->reset();
        pe_arrays[i]->reset();
    }
    for(auto sfu : sfus) { sfu->reset(); }
}

void npu_t::update_tile_size() {
    validate_active_components();
    dram->update_tile_size(scheduler);
    multi_chip->update_tile_size(scheduler);

    for(unsigned i = 0; i < multi_chip->get_number_of_active_chips(); i++) {
        global_buffers[i]->update_tile_size(scheduler);
        pe_arrays[i]->update_tile_size(scheduler);
        pe_arrays[i]->set_psum_retention_scope(scheduler);
    }
}

// SFU (plan/plan_sfu.md): fused-activation cost event for one finished layer.
//
// The event contract: the activation fires exactly once per NETWORK-valid output element
// (B x K x P x Q; mapping padding excluded), only after every C/R/S reduction completed --
// which is guaranteed here because the whole layer, including the final psum flush and
// repetition scaling, is already accounted. Nebula's forward() remains the functional
// owner; this only generates the hardware cost event.
void npu_t::apply_fused_sfu_activation(unsigned m_index) {
    nebula::layer_t *current_layer = network->layers[m_index];
    const unsigned activation_index = static_cast<unsigned>(current_layer->activation_type);
    // An undeclared activation behaves as identity in Nebula, so it maps to the linear
    // bypass; every other name must map to a supported SFU operation or fail fast.
    const std::string activation_name =
        (activation_index == nebula::UNDEFINED_ACTIVATION ||
         activation_index >= nebula::NUM_ACTIVATION_TYPES)
        ? "linear" : nebula::activation_type_str[activation_index];

    // Valid output elements: the NETWORK dimensions clamped by what the mapping actually
    // covered. min(mapped, layer) excludes mapping padding (padded outputs are never
    // committed) AND partial coverage (a mapping that simulates 48 of 64 output channels
    // must not charge activation for outputs the timing model never produced).
    size_t layer_k, layer_p, layer_q;
    if(scheduler->layer_name == layer_name_t::CONNECTED_LAYER) {
        layer_k = current_layer->output_size ? current_layer->output_size
                                             : current_layer->output_channel;
        layer_p = 1;
        layer_q = 1;
    } else {
        layer_k = current_layer->output_channel;
        layer_p = current_layer->output_height;
        layer_q = current_layer->output_width;
        if(layer_k*layer_p*layer_q == 0) {
            layer_k = current_layer->output_size;
            layer_p = 1;
            layer_q = 1;
        }
    }
    const std::vector<unsigned> mapped =
        scheduler->mapping_table->calculate_total_parameter_size();
    auto covered = [](unsigned m_mapped, size_t m_layer) -> size_t {
        if(m_layer == 0) return m_mapped ? m_mapped : 1;
        if(m_mapped == 0) return m_layer;
        return std::min(static_cast<size_t>(m_mapped), m_layer);
    };
    const size_t valid_elements =
        covered(mapped[parameter_type_t::BATCH_SIZE], network->batch_size)*
        covered(mapped[parameter_type_t::OUTPUT_CHANNEL], layer_k)*
        covered(mapped[parameter_type_t::OUTPUT_HEIGHT], layer_p)*
        covered(mapped[parameter_type_t::OUTPUT_WIDTH], layer_q);

    if(sfus.empty()) {
        // Compatibility policy: no [sfu] section keeps every legacy number bit-identical,
        // but a nonlinear activation that DID execute must be stated as out of scope.
        if(activation_name != "linear" && valid_elements > 0) {
            layer_stats[m_index]->mark_unmodeled_activation(
                "activation '" + activation_name + "' over " +
                std::to_string(valid_elements) + " output element(s) executed with no [sfu]"
                " section; its cycles/traffic/energy are ABSENT from this report");
        }
        return;
    }

    sfu_op_t op;
    if(!sfu_t::op_from_name(activation_name, &op)) {
        // Normally unreachable: run() validates every layer's activation up front.
        std::cerr << "Error: layer " << m_index << " activation '" << activation_name
                  << "' is not supported by the SFU; unsupported operations fail fast"
                  << " instead of silently falling back to ReLU" << std::endl;
        exit(1);
    }

    // Distribute outputs over the chips that own DISTINCT output elements only. A chip
    // factor on a reduction dimension (C/R/S) replicates the SAME outputs as partial
    // sums, which merge before the activation fires -- splitting the element count over
    // those replicas would understate every output-owning chip's SFU window.
    const unsigned active_chips = scheduler->num_active_chips_x*scheduler->num_active_chips_y;
    unsigned output_chips = scheduler->mapping_table->calculate_output_partition_chips();
    output_chips = std::max(1u, std::min(output_chips, active_chips));
    sfu_invocation_t combined;
    const size_t base = valid_elements/output_chips;
    const size_t remainder = valid_elements % output_chips;
    for(unsigned c = 0; c < output_chips; ++c) {
        const size_t share = base + (c < remainder ? 1 : 0);
        combined.merge_parallel(sfus[c]->elementwise_invocation(op, share));
    }
    layer_stats[m_index]->sfu_contract_note =
        "fused post-accumulator activation; input/weight/output DRAM traffic is unchanged"
        " by the SFU (fused invariant)";
    layer_stats[m_index]->apply_sfu_activation(combined,
        num_processors*sfus[0]->get_num_units(), sfus[0]->get_lanes(),
        sfus[0]->get_static_energy_per_cycle());
}

// Phase-7 (plan_sfu.md): standalone softmax on the SFU multi-pass microprogram
// (max -> subtract -> exp -> sum -> reciprocal -> normalize). The layer has no mapping
// section; the SFU of chip 0 is the modeled execution resource and the layer latency is
// its serial multi-pass window.
void npu_t::run_standalone_softmax(unsigned m_index, const std::string &m_accelerator_config,
                                   const std::string &m_network_config) {
    nebula::layer_t *softmax_layer = network->layers[m_index];
    const size_t rows = network->batch_size;
    const size_t row_length = softmax_layer->output_size;

    std::cout << "The network configuration of #" << m_index
              << " layer (standalone softmax on SFU)" << std::endl;
    std::cout << " * softmax rows x length : " << rows << " x " << row_length << std::endl;

    sfu_invocation_t invocation = sfus[0]->softmax_invocation(rows, row_length);

    stats_t *stats = new stats_t();
    sfu_layer_stats.push_back(stats);

    // Clock contract: the SFU runs on the chip clock, so the same single-domain check
    // update_stats() applies to mapped layers holds here.
    const double frequencies[5] = {
        pe_arrays[0]->pes[0]->clock_mhz(),
        pe_arrays[0]->clock_mhz(),
        global_buffers[0]->clock_mhz(),
        multi_chip->clock_mhz(),
        dram->clock_mhz()};
    bool single_clock = frequencies[0] > 0.0;
    for(unsigned f = 1; f < 5; ++f) {
        if(frequencies[f] != frequencies[0]) single_clock = false;
    }
    const std::string clock_note = single_clock
        ? "all modeled components share one clock"
        : "mixed or undeclared clock domains; cycles are not convertible to time";

    // Scope: the operand tensor is materialized by the neighboring layers (the previous
    // layer already stored its output off-chip), and streaming it to/from the SFU is not
    // yet part of this model -- stated instead of silently read as included.
    const size_t stream_bytes = runtime_datatypes().storage_bytes(
        data_type_t::OUTPUT, rows*row_length);
    stats->sfu_contract_note =
        "standalone softmax executes on chip 0's SFU; operand streaming (~" +
        std::to_string(stream_bytes) + " bytes each way at output precision) between the"
        " memory hierarchy and the SFU is NOT included (materialized-tensor accounting is"
        " a later phase)";

    stats->record_sfu_only_layer(invocation,
        num_processors*sfus[0]->get_num_units(), sfus[0]->get_lanes(),
        sfus[0]->get_static_energy_per_cycle(),
        single_clock ? frequencies[0] : 0.0, single_clock, clock_note);

    std::cout << "The simulation result of #" << m_index << " layer" << std::endl;
    std::cout << std::endl;
    const std::string output_file_name = m_accelerator_config + "_" + m_network_config +
                                         "_layer_" + std::to_string(m_index) + ".txt";
    std::cout << output_file_name << std::endl;
    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out);
    output_file << "Standalone softmax layer executed on the SFU (multi-pass microprogram)"
                << std::endl << std::endl;
    stats->print_results(output_file);
    output_file.close();
    network_stats->update_network_stats(stats);
}
