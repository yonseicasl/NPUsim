#include <thread>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <cstdio>
#include <map>
#include <set>
#include <stdexcept>
#include <unistd.h>


#include "npu.h"
#include "energy_units.h"
#include "config.h"
#include "datatype.h"
#include "interconnect_timing.h"

namespace {

size_t pool_reduction_operations(const workload_geometry_t &geometry,
                                 size_t output_begin, size_t output_count) {
    const size_t spatial = static_cast<size_t>(geometry.output_height)*
                           geometry.output_width;
    const size_t window = static_cast<size_t>(geometry.kernel_height)*
                          geometry.kernel_width;
    size_t reductions = 0;
    for(size_t flat = output_begin; flat < output_begin + output_count; ++flat) {
        const size_t location = flat % spatial;
        const size_t output_h = location/geometry.output_width;
        const size_t output_w = location%geometry.output_width;
        size_t valid = 0;
        for(unsigned kernel_h = 0; kernel_h < geometry.kernel_height; ++kernel_h) {
            const int64_t input_h =
                static_cast<int64_t>(output_h*geometry.stride_height) -
                geometry.padding_height +
                static_cast<int64_t>(kernel_h)*geometry.dilation_height;
            if(input_h < 0 || input_h >= static_cast<int64_t>(geometry.input_height)) continue;
            for(unsigned kernel_w = 0; kernel_w < geometry.kernel_width; ++kernel_w) {
                const int64_t input_w =
                    static_cast<int64_t>(output_w*geometry.stride_width) -
                    geometry.padding_width +
                    static_cast<int64_t>(kernel_w)*geometry.dilation_width;
                if(input_w >= 0 && input_w < static_cast<int64_t>(geometry.input_width)) ++valid;
            }
        }
        const size_t samples = geometry.mode == "average" && geometry.count_include_pad
            ? window : valid;
        if(samples > 0) {
            if(reductions > std::numeric_limits<size_t>::max() - (samples - 1)) {
                throw std::runtime_error("pool reduction operation count overflows");
            }
            reductions += samples - 1;
        }
    }
    return reductions;
}

} // namespace

npu_t::npu_t() :
 num_processors(1),
 num_pes(1),
 compression_type(compression_type_t::DENSE),
 num_skipped_timing_layers(0),
 multi_chip(NULL),
 dram(NULL),
 decomp(NULL),
 kvcache(NULL),
 workload(NULL),
 workload_lifetime(NULL),
 executable_ir_mode(false),
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
    delete decomp;
    delete kvcache;

    delete multi_chip;
    delete dram;

	// Free the memory for the network.
	delete network;
    delete workload;
    delete workload_lifetime;

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
        // Initialize the weight-decompression engine (opt-in; evaluation.md Sec 4).
        else if(section_config.name == "decomp" || section_config.name == "DECOMP") {
            if(decomp != NULL) {
                std::cerr << "Error: duplicate [decomp] section" << std::endl;
                exit(1);
            }
            decomp = new decomp_t(section_config);
        }
        else if(section_config.name == "kvcache" || section_config.name == "KVCACHE") {
            if(kvcache != NULL) {
                std::cerr << "Error: duplicate [kvcache] section" << std::endl;
                exit(1);
            }
            kvcache = new kvcache_t(section_config);
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
    executable_ir_mode = m_network_config.size() >= 5 &&
        m_network_config.substr(m_network_config.size() - 5) == ".json";
    if(executable_ir_mode) {
        workload = new workload_graph_t();
        try {
            workload->load(m_network_config);
            const std::string generated = workload->legacy_network_config();
            char temporary_path[] = "/tmp/npusim-executable-XXXXXX";
            const int temporary_fd = mkstemp(temporary_path);
            if(temporary_fd == -1) throw std::runtime_error("cannot create transitional network config");
            close(temporary_fd);
            {
                std::ofstream output(temporary_path, std::ios::out | std::ios::trunc);
                if(!output.good()) {
                    std::remove(temporary_path);
                    throw std::runtime_error("cannot write transitional network config");
                }
                output << generated;
            }
            network->init(temporary_path);
            std::remove(temporary_path);
            override_executable_layer_geometry();
            std::map<std::string, size_t> runtime_bytes;
            for(const workload_tensor_t &tensor_value : workload->tensors) {
                // Alias views take the classification and bytes of their storage
                // tensor so lifetime accounting never double-books a reshape.
                const workload_tensor_t &storage = workload->storage_tensor(tensor_value.id);
                data_type_t type = data_type_t::OUTPUT;
                if(storage.kind == "parameter" || storage.kind == "buffer" ||
                   storage.kind == "constant") {
                    type = data_type_t::WEIGHT;
                } else if(std::find(workload->inputs.begin(), workload->inputs.end(),
                                    storage.id) != workload->inputs.end()) {
                    type = data_type_t::INPUT;
                }
                runtime_bytes[tensor_value.id] = runtime_datatypes().storage_bytes(
                    type, storage.elements());
            }
            const size_t per_chip_capacity = global_buffers[0]->tensor_residency_capacity();
            if(num_processors != 0 && per_chip_capacity > std::numeric_limits<size_t>::max()/num_processors) {
                throw std::runtime_error("aggregate GLB residency capacity overflows");
            }
            workload_lifetime = new workload_lifetime_t(
                *workload, per_chip_capacity*num_processors, runtime_bytes);
            std::cout << "# Frontend IR: " << workload->schema_version
                      << " model=" << workload->model_name
                      << " source_sha256=" << workload->source_sha256 << std::endl;
        } catch(const std::exception &error) {
            std::cerr << "Error: " << error.what() << std::endl;
            exit(1);
        }
    } else {
        network->init(m_network_config);
    }
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
    bind_executable_mappings();
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
void npu_t::override_executable_layer_geometry() {
    if(workload == NULL) return;
    if(network->layers.size() != workload->operations.size()) {
        throw std::runtime_error("transitional layer count disagrees with executable operations");
    }
    for(size_t index = 0; index < workload->operations.size(); ++index) {
        const workload_operation_t &operation = workload->operations[index];
        nebula::layer_t *current = network->layers[index];
        const workload_tensor_t &input = workload->tensor(operation.inputs.front());
        const workload_tensor_t &output = workload->tensor(operation.outputs.front());
        const size_t batch = input.shape.empty() ? 1 : input.shape.front();
        if(batch == 0 || input.elements() % batch != 0 || output.elements() % batch != 0 ||
           input.elements()/batch > std::numeric_limits<unsigned>::max() ||
           output.elements()/batch > std::numeric_limits<unsigned>::max()) {
            throw std::runtime_error("operation " + operation.id + " tensor volume exceeds Nebula bridge range");
        }
        current->input_size = static_cast<unsigned>(input.elements()/batch);
        current->output_size = static_cast<unsigned>(output.elements()/batch);
        if(input.shape.size() == 4) {
            current->input_channel = static_cast<unsigned>(input.shape[1]);
            current->input_height = static_cast<unsigned>(input.shape[2]);
            current->input_width = static_cast<unsigned>(input.shape[3]);
        } else {
            current->input_channel = 1;
            current->input_height = 1;
            current->input_width = current->input_size;
        }
        if(output.shape.size() == 4) {
            current->output_channel = static_cast<unsigned>(output.shape[1]);
            current->output_height = static_cast<unsigned>(output.shape[2]);
            current->output_width = static_cast<unsigned>(output.shape[3]);
        } else {
            current->output_channel = current->output_size;
            current->output_height = 1;
            current->output_width = 1;
        }
        if(operation.kind == WORKLOAD_LINEAR) {
            current->input_channel = operation.geometry.input_features;
            current->input_height = current->input_width = 1;
            current->output_channel = operation.geometry.output_features;
            current->output_height = current->output_width = 1;
            current->weight_size = operation.geometry.input_features*operation.geometry.output_features;
        } else if(operation.kind == WORKLOAD_CONV2D) {
            if(operation.geometry.stride_height != operation.geometry.stride_width ||
               operation.geometry.dilation_height != 1 || operation.geometry.dilation_width != 1) {
                throw std::runtime_error("operation " + operation.id +
                    " uses asymmetric stride or dilation, which the mapped MAC core cannot represent yet");
            }
            current->filter_height = operation.geometry.filter_height;
            current->filter_width = operation.geometry.filter_width;
            current->filter_size = operation.geometry.filter_height*operation.geometry.filter_width;
            current->stride = operation.geometry.stride_height;
            current->padding_h = operation.geometry.padding_height;
            current->padding_w = operation.geometry.padding_width;
            current->group = operation.geometry.groups;
            current->num_filters = operation.geometry.output_channels;
            current->weight_size = operation.geometry.output_channels*
                (operation.geometry.input_channels/operation.geometry.groups)*
                operation.geometry.filter_height*operation.geometry.filter_width;
        }
    }
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
    std::cout << "# Run the network" << std::endl;

    if(executable_ir_mode && sfus.empty()) {
        for(const workload_operation_t &operation : workload->operations) {
            if(!operation.mapping_required) {
                std::cerr << "Error: executable operation " << operation.id
                          << " requires an [sfu] section for non-MAC timing" << std::endl;
                exit(1);
            }
        }
    }
    // Validate all SFU capabilities before the first operation, so a missing primitive
    // never leaves a partially simulated DAG behind.
    if(!sfus.empty()) {
        auto require_sfu = [&](unsigned index, sfu_op_t op) {
            if(!sfus[0]->op_supported(op)) {
                std::cerr << "Error: network operation " << index << " needs SFU primitive '"
                          << sfu_t::op_name(op) << "', outside this architecture's"
                          << " [sfu] supported_ops contract" << std::endl;
                exit(1);
            }
        };
        if(executable_ir_mode) {
            for(unsigned index = 0; index < workload->operations.size(); ++index) {
                const workload_operation_t &operation = workload->operations[index];
                if(operation.kind == WORKLOAD_SOFTMAX) {
                    const sfu_op_t ops[5] = {SFU_OP_VMAX, SFU_OP_VADD, SFU_OP_EXP,
                                             SFU_OP_RECIP, SFU_OP_VMUL};
                    for(unsigned op = 0; op < 5; ++op) require_sfu(index, ops[op]);
                } else if(operation.kind == WORKLOAD_POOL2D) {
                    require_sfu(index, operation.geometry.mode == "max" ? SFU_OP_VMAX : SFU_OP_VADD);
                    if(operation.geometry.mode == "average") require_sfu(index, SFU_OP_VMUL);
                } else if(operation.kind == WORKLOAD_ELEMENTWISE) {
                    require_sfu(index, operation.geometry.elementwise_operator == "add"
                        ? SFU_OP_VADD : SFU_OP_VMUL);
                } else if(operation.kind == WORKLOAD_BATCH_NORM) {
                    require_sfu(index, SFU_OP_VMUL);
                    require_sfu(index, SFU_OP_VADD);
                }
            }
        } else {
            for(unsigned index = 0; index < network->num_layers; index++) {
                if(network->layers[index]->layer_type == nebula::SOFTMAX_LAYER) {
                    const sfu_op_t ops[5] = {SFU_OP_VMAX, SFU_OP_VADD, SFU_OP_EXP,
                                             SFU_OP_RECIP, SFU_OP_VMUL};
                    for(unsigned op = 0; op < 5; ++op) require_sfu(index, ops[op]);
                    continue;
                }
                if(network->layers[index]->layer_type != nebula::CONVOLUTIONAL_LAYER &&
                   network->layers[index]->layer_type != nebula::CONNECTED_LAYER) continue;
                const unsigned type = static_cast<unsigned>(network->layers[index]->activation_type);
                if(type == nebula::UNDEFINED_ACTIVATION || type >= nebula::NUM_ACTIVATION_TYPES) continue;
                sfu_op_t op;
                if(!sfu_t::op_from_name(nebula::activation_type_str[type], &op)) {
                    std::cerr << "Error: network layer " << index << " activation '"
                              << nebula::activation_type_str[type] << "' is unsupported" << std::endl;
                    exit(1);
                }
                require_sfu(index, op);
            }
        }
    }

    const unsigned num_iteration = 1;
    for(unsigned iteration = 0; iteration < num_iteration; iteration++) {
#ifdef FUNCTIONAL
        if(executable_ir_mode) {
            std::cerr << "Error: executable IR currently supports timing-only builds; "
                      << "functional mode requires a tensor artifact" << std::endl;
            exit(1);
        }
#endif
        if(!executable_ir_mode) network->load_data(iteration);
        num_skipped_timing_layers = 0;
        unsigned mapping_index = 0;

        for(unsigned index = 0; index < network->num_layers; index++) {
            const workload_operation_t *operation = executable_ir_mode
                ? &workload->operations[index] : NULL;
            workload_residency_plan_t residency;
            if(executable_ir_mode) residency = workload_lifetime->plan(index);
            if(!executable_ir_mode) {
                network->layers[index]->input_data = index > 0
                    ? network->layers[index-1]->output_data : network->input_data;
            }
            const bool mapped = executable_ir_mode
                ? operation->mapping_required
                : (network->layers[index]->layer_type == nebula::CONVOLUTIONAL_LAYER ||
                   network->layers[index]->layer_type == nebula::CONNECTED_LAYER);
            if(mapped) {
                const unsigned stats_index = executable_ir_mode ? mapping_index : index;
                if(stats_index >= schedulers.size() || stats_index >= layer_stats.size()) {
                    std::cerr << "Error: no mapping section for network operation " << index << std::endl;
                    exit(1);
                }
                const bool convolution = executable_ir_mode
                    ? operation->kind == WORKLOAD_CONV2D
                    : network->layers[index]->layer_type == nebula::CONVOLUTIONAL_LAYER;
                schedulers[stats_index]->layer_name = convolution
                    ? layer_name_t::CONVOLUTIONAL_LAYER : layer_name_t::CONNECTED_LAYER;
                layer = network->layers[index];
                dram->connect_layer(layer);
                scheduler = schedulers[stats_index];
                const unsigned global_buffer_repetitions =
                    scheduler->mapping_table->calculate_active_component(component_type_t::GLOBAL_BUFFER);
                {
                    const std::vector<unsigned> mapped_size =
                        scheduler->mapping_table->calculate_total_parameter_size();
                    const unsigned layer_c = scheduler->layer_name == layer_name_t::CONNECTED_LAYER
                        ? layer->input_size : layer->input_channel;
                    const struct { const char *name; unsigned mapped_size; unsigned layer_size; } dims[] = {
                        {"K", mapped_size[parameter_type_t::OUTPUT_CHANNEL], layer->output_channel},
                        {"P", mapped_size[parameter_type_t::OUTPUT_HEIGHT], layer->output_height},
                        {"Q", mapped_size[parameter_type_t::OUTPUT_WIDTH], layer->output_width},
                        {"C", mapped_size[parameter_type_t::INPUT_CHANNEL], layer_c},
                    };
                    for(const auto &dim : dims) {
                        if(dim.layer_size == 0) continue;
                        if(dim.mapped_size > dim.layer_size) {
                            std::cerr << "Warning: layer " << index << " mapping pads " << dim.name
                                      << " from " << dim.layer_size << " to " << dim.mapped_size
                                      << " (padded work is charged as compute)" << std::endl;
                        } else if(dim.mapped_size < dim.layer_size) {
                            if(executable_ir_mode) {
                                std::cerr << "Error: executable operation " << operation->id
                                          << " mapping covers only " << dim.mapped_size << " of "
                                          << dim.layer_size << " in " << dim.name << std::endl;
                                exit(1);
                            }
                            std::cerr << "Warning: layer " << index << " mapping covers only "
                                      << dim.mapped_size << " of " << dim.layer_size << " in " << dim.name
                                      << " (layer is partially simulated)" << std::endl;
                        }
                    }
                }

                print_network_configuration(index, stats_index);
                reset();
                update_tile_size();
                while(!is_idle()) {
                    execute();
                    transfer_data_to_pe();
                    transfer_data_to_pe_array();
                    transfer_data_to_global_buffer();
                    transfer_data_to_multi_chip();
                    request_to_dram();
                    request_to_multi_chip();
                    request_to_global_buffer();
                    request_to_pe_array();
                }
                for(unsigned chip = 0; chip < pe_arrays.size(); chip++) {
                    pe_arrays[chip]->flush_psum_writeback();
                }
                multi_chip->flush_output_writeback();
                for(unsigned chip = 0; chip < pe_arrays.size(); chip++) {
                    for(unsigned pe = 0; pe < pe_arrays[chip]->get_number_of_pes(); pe++) {
                        pe_arrays[chip]->pes[pe]->suppress_streaming_cycles();
                    }
                }
                layer_stats[stats_index]->update_stats(pe_arrays, global_buffers, multi_chip, dram);
                if(executable_ir_mode) {
                    workload_lifetime->commit(index, &residency);
                    const bool input_in_glb = !residency.inputs.empty() &&
                        residency.inputs.front() == WORKLOAD_RESIDENCY_GLB;
                    const size_t retained_inputs = static_cast<size_t>(std::count(
                        residency.retain_inputs.begin(), residency.retain_inputs.end(), true));
                    const std::string note = std::string("input ") +
                        (input_in_glb ? "GLB-resident" : "DRAM-backed") + ", output " +
                        (residency.retain_output ? "retained in GLB" : "materialized to DRAM") +
                        (retained_inputs ? "; " + std::to_string(retained_inputs) +
                         " future-use input(s) pinned" : "") +
                        "; GLB occupancy " + std::to_string(residency.occupied_before) + " -> " +
                        std::to_string(residency.occupied_after) + " / " +
                        std::to_string(residency.capacity) + " bytes";
                    layer_stats[stats_index]->apply_graph_residency(
                        input_in_glb, residency.retain_output, note);
                }
                const input_halo_reuse_t input_halo = scheduler->mapping_table->input_halo_reuse();
                const bool halo_capacity_sufficient = !global_buffers.empty() &&
                    global_buffers[0]->can_retain_input_halo(input_halo.working_set_elements);
                apply_fused_sfu_activation(index, stats_index);
                apply_kv_cache_read(stats_index);
                apply_weight_decompression(stats_index);
                layer_stats[stats_index]->scale_serial_repetitions(
                    global_buffer_repetitions, scheduler->mapping_table->datatype_repetitions(),
                    input_halo, halo_capacity_sufficient);
                print_layerwise_results(m_accelerator_config, m_network_config, index, stats_index);
                ++mapping_index;
            } else if(executable_ir_mode) {
                run_standalone_graph_operation(index, residency,
                    m_accelerator_config, m_network_config);
            } else if(network->layers[index]->layer_type == nebula::SOFTMAX_LAYER && !sfus.empty()) {
                run_standalone_softmax(index, m_accelerator_config, m_network_config);
            } else {
                ++num_skipped_timing_layers;
                std::cerr << "Warning: network layer " << index
                          << " is excluded from accelerator timing (only convolution/connected are supported)"
                          << std::endl;
            }
#ifdef FUNCTIONAL
            network->layers[index]->forward();
#endif
        }
        if(executable_ir_mode && mapping_index != schedulers.size()) {
            std::cerr << "Error: executable DAG did not consume every mapping" << std::endl;
            exit(1);
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
    if(decomp != NULL) { decomp->print_specification(); }
    if(kvcache != NULL) { kvcache->print_specification(); }
}

// Print out the network stats (e.g., tile size)
void npu_t::print_network_configuration(unsigned m_layer_index, unsigned m_stats_index) {
    std::cout << "The network configuration of #" << m_layer_index << " layer" << std::endl;
    layer_stats[m_stats_index]->print_stats();
}

// Print out the simulation result
void npu_t::print_layerwise_results(const std::string m_accelerator_config,
                                    const std::string m_network_config,
                                    unsigned m_layer_index, unsigned m_stats_index) {
    std::cout << "The simulation result of #" << m_layer_index << " layer" << std::endl;
    std::cout << std::endl;

    // Concatenate the name of output file.
    std::string output_file_name = m_accelerator_config + "_" + m_network_config + "_layer_" + std::to_string(m_layer_index) + ".txt";
    std::cout << output_file_name << std::endl;

    std::ofstream output_file;
    output_file.open(output_file_name, std::ios::out);
    print_workload_provenance(output_file);
    layer_stats[m_stats_index]->print_stats(output_file);

    layer_stats[m_stats_index]->print_results(output_file);

#ifdef DRAMSIM3
    dram->print_result();
#endif

    output_file.close();
    network_stats->update_network_stats(layer_stats[m_stats_index]);
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
    print_workload_provenance(output_file);
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
    if(decomp != NULL) { decomp->reset(); }
    if(kvcache != NULL) { kvcache->reset(); }
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
void npu_t::apply_fused_sfu_activation(unsigned m_layer_index, unsigned m_stats_index) {
    nebula::layer_t *current_layer = network->layers[m_layer_index];
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
            layer_stats[m_stats_index]->mark_unmodeled_activation(
                "activation '" + activation_name + "' over " +
                std::to_string(valid_elements) + " output element(s) executed with no [sfu]"
                " section; its cycles/traffic/energy are ABSENT from this report");
        }
        return;
    }

    sfu_op_t op;
    if(!sfu_t::op_from_name(activation_name, &op)) {
        // Normally unreachable: run() validates every layer's activation up front.
        std::cerr << "Error: layer " << m_layer_index << " activation '" << activation_name
                  << "' is not supported by the SFU; unsupported operations fail fast"
                  << " instead of silently falling back to ReLU" << std::endl;
        exit(1);
    }

    // Phase-2: final_output_tile events. The live pass counts one event per output
    // commit at the multi-chip -> DRAM boundary (the reduction-complete, once-per-element
    // point RE1/DR6 establish); output-datatype repetitions replay it with DISTINCT
    // output tiles, while reduction repetitions revisit the same identity and commit
    // once. The identity gate: the committed (mapped, padding-included) elements must
    // reproduce the mapping's output volume exactly -- otherwise the commit stream is
    // not a trustworthy event source and the model falls back to the layer-granular
    // single invocation, saying so (plan gate 1: no timing hookup on an uncertain
    // reduction-completion event).
    const std::vector<unsigned> repetitions = scheduler->mapping_table->datatype_repetitions();
    const size_t output_repetitions = std::max(1u, repetitions[data_type_t::OUTPUT]);
    const size_t commit_events = multi_chip->final_output_tile_events*output_repetitions;
    const size_t committed_elements = multi_chip->final_output_tile_elements*output_repetitions;
    const size_t mapped_output_volume =
        static_cast<size_t>(mapped[parameter_type_t::BATCH_SIZE])*
        mapped[parameter_type_t::OUTPUT_CHANNEL]*
        mapped[parameter_type_t::OUTPUT_HEIGHT]*
        mapped[parameter_type_t::OUTPUT_WIDTH];
    const bool commit_identity_ok = commit_events > 0 &&
                                    committed_elements == mapped_output_volume;
    const size_t model_events = commit_identity_ok ? commit_events : 1;
    layer_stats[m_stats_index]->sfu_commit_events = model_events;
    layer_stats[m_stats_index]->sfu_commit_note = commit_identity_ok
        ? std::to_string(multi_chip->final_output_tile_events) + " commit(s) x " +
          std::to_string(output_repetitions) + " output repetition(s); identity OK: " +
          std::to_string(committed_elements) + " committed = mapped output volume"
        : "identity MISMATCH: " + std::to_string(committed_elements) +
          " committed vs " + std::to_string(mapped_output_volume) +
          " mapped output elements -- layer-granular single-invocation fallback";

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
        combined.merge_parallel(sfus[c]->elementwise_invocation(op, share, model_events));
    }
    layer_stats[m_stats_index]->sfu_contract_note =
        "fused post-accumulator activation; input/weight/output DRAM traffic is unchanged"
        " by the SFU (fused invariant)";
    if(!sfus[0]->get_precision_note().empty()) {
        layer_stats[m_stats_index]->sfu_contract_note += "; " + sfus[0]->get_precision_note();
    }
    layer_stats[m_stats_index]->sfu_profile_reference = sfus[0]->get_profile_reference();
    layer_stats[m_stats_index]->set_sfu_activation(combined,
        num_processors*sfus[0]->get_num_units(), sfus[0]->get_lanes(),
        sfus[0]->get_static_energy_per_cycle(), sfus[0]->get_queue_depth());
}

void npu_t::apply_weight_decompression(unsigned m_stats_index) {
    if(decomp == NULL) return;
    // Dense weight footprint from the mapping: K x C x R x S (grouped conv folds the
    // group into C already via the mapping factors), at the runtime weight precision.
    const std::vector<unsigned> mapped =
        scheduler->mapping_table->calculate_total_parameter_size();
    const size_t weight_elements =
        static_cast<size_t>(mapped[parameter_type_t::OUTPUT_CHANNEL])*
        mapped[parameter_type_t::INPUT_CHANNEL]*
        mapped[parameter_type_t::FILTER_HEIGHT]*
        mapped[parameter_type_t::FILTER_WIDTH];
    const size_t dense_weight_bytes =
        runtime_datatypes().storage_bytes(data_type_t::WEIGHT, weight_elements);
    // DRAM link rate (bytes/cycle) for the reported decoder ratio; the weight DRAM saving
    // itself is applied as the compression ratio on the measured weight DRAM cycles.
    const double dram_bytes_per_cycle = static_cast<double>(dram->get_bitwidth())/8.0;
    // Decoder/queue granularity: the weight tile the multi-chip stages off-chip.
    const size_t tile_elements = std::max<size_t>(1,
        multi_chip->tile_size[data_type_t::WEIGHT]);
    const size_t tile_bytes =
        runtime_datatypes().storage_bytes(data_type_t::WEIGHT, tile_elements);
    layer_stats[m_stats_index]->apply_decompression(decomp, dense_weight_bytes,
                                                    dram_bytes_per_cycle, tile_bytes);
}

void npu_t::apply_kv_cache_read(unsigned m_stats_index) {
    if(kvcache == NULL) return;
    if(kvcache->attention_enabled()) {
        // Attention consumer mode: cost the KV read/write as a DEDICATED stream on the
        // live DRAM component -- device accesses per line, off-chip link beats, and the
        // open-page row activations of two contiguous sequential streams (K then V; that
        // sequential layout is the declared address model) -- the same machinery the
        // standalone-softmax operand stream uses. KV is activation-like, so it is priced
        // at the OUTPUT datatype's declared unit costs, exactly like that stream.
        const data_type_t T = data_type_t::OUTPUT;
        const double per_activation_cycle =
            (dram->t_ras_cycle > 0.0 && dram->t_rp_cycle > 0.0)
                ? dram->t_ras_cycle + dram->t_rp_cycle : dram->u_row_miss_cycle;
        kv_stream_cost_t cost;
        // One sequential stream of m_bytes: device/link makespans overlap packet-level
        // (max), the busiest bank's row activations serialize on top.
        auto read_stream = [&](size_t m_bytes, double *m_energy) {
            const size_t bits = m_bytes*8;
            const size_t accesses = (bits + std::max(1u, dram->line_size[T]) - 1)/
                                    std::max(1u, dram->line_size[T]);
            const size_t beats = (bits + std::max(1u, dram->get_bitwidth()) - 1)/
                                 std::max(1u, dram->get_bitwidth());
            cost.link_transactions += beats;
            double cycle = std::max(static_cast<double>(accesses)*dram->u_read_cycle[T],
                                    static_cast<double>(beats)*dram->u_transfer_cycle);
            *m_energy += static_cast<double>(accesses)*dram->u_read_energy[T] +
                         static_cast<double>(beats)*dram->u_transfer_energy;
            if(dram->row_buffer_bytes > 0) {
                const dram_row_activation_cost_t rows = dram_row_activations(
                    m_bytes, dram->row_buffer_bytes, dram->num_banks);
                cost.row_activations += rows.activations;
                cycle += static_cast<double>(rows.busiest_bank)*per_activation_cycle;
                *m_energy += static_cast<double>(rows.activations)*dram->u_row_miss_energy;
            }
            return cycle;
        };
        // K and V passes fetch the COMPRESSED halves; the decoder reconstitutes dense at
        // its declared throughput, half per pass.
        const size_t comp_half = kvcache->compressed_read_bytes()/2;
        const double decoder_half =
            kvcache->decoder_cycles(kvcache->dense_read_bytes())/2.0;
        cost.k_supply_cycle = read_stream(comp_half, &cost.dram_energy) + decoder_half;
        cost.v_supply_cycle = read_stream(comp_half, &cost.dram_energy) + decoder_half;
        // Cache append: the current token's K/V, stored at the cache's compression. One
        // token is far below a DRAM row, so no row-activation charge for the append.
        const size_t dense_w = kvcache->kv_write_bytes();
        const size_t write_b = kvcache->bypassed() ? dense_w
            : std::max<size_t>(1, static_cast<size_t>(
                  static_cast<double>(dense_w)*kvcache->compressed_read_bytes()/
                  std::max<size_t>(1, kvcache->dense_read_bytes())));
        const size_t wbits = write_b*8;
        const size_t waccesses = (wbits + std::max(1u, dram->line_size[T]) - 1)/
                                 std::max(1u, dram->line_size[T]);
        const size_t wbeats = (wbits + std::max(1u, dram->get_bitwidth()) - 1)/
                              std::max(1u, dram->get_bitwidth());
        cost.link_transactions += wbeats;
        cost.write_cycle = std::max(static_cast<double>(waccesses)*dram->u_write_cycle[T],
                                    static_cast<double>(wbeats)*dram->u_transfer_cycle);
        cost.dram_energy += static_cast<double>(waccesses)*dram->u_write_energy[T] +
                            static_cast<double>(wbeats)*dram->u_transfer_energy;
        layer_stats[m_stats_index]->apply_attention_step(kvcache, cost);
        return;
    }
    // The dense weight footprint of this layer is the per-byte DRAM cost reference used by
    // stats_t (measured weight-DRAM cycles / dense weight bytes). Same K x C x R x S the
    // decompression path computes.
    const std::vector<unsigned> mapped =
        scheduler->mapping_table->calculate_total_parameter_size();
    const size_t weight_elements =
        static_cast<size_t>(mapped[parameter_type_t::OUTPUT_CHANNEL])*
        mapped[parameter_type_t::INPUT_CHANNEL]*
        mapped[parameter_type_t::FILTER_HEIGHT]*
        mapped[parameter_type_t::FILTER_WIDTH];
    const size_t dense_weight_bytes =
        runtime_datatypes().storage_bytes(data_type_t::WEIGHT, weight_elements);
    layer_stats[m_stats_index]->apply_kv_cache_read(kvcache, dense_weight_bytes);
}

sfu_operand_stream_t npu_t::graph_operand_stream(
    const workload_operation_t &m_operation, const workload_residency_plan_t &m_plan) {
    if(m_plan.inputs.size() != m_operation.inputs.size()) {
        throw std::runtime_error("graph residency/input count mismatch");
    }
    sfu_operand_stream_t stream;
    stream.active = true;
    global_buffer_t *glb = global_buffers[0];
    size_t resident_inputs = 0;
    for(size_t index = 0; index < m_operation.inputs.size(); ++index) {
        // Datatype classification follows the storage tensor: reading through an
        // elided reshape view streams the aliased buffer, not a new tensor.
        const workload_tensor_t &tensor_value = workload->storage_tensor(m_operation.inputs[index]);
        data_type_t type = data_type_t::OUTPUT;
        if(tensor_value.kind == "parameter" || tensor_value.kind == "buffer" ||
           tensor_value.kind == "constant") type = data_type_t::WEIGHT;
        else if(std::find(workload->inputs.begin(), workload->inputs.end(), tensor_value.id) !=
                workload->inputs.end()) type = data_type_t::INPUT;
        const size_t elements = tensor_value.elements();
        stream.ingress_bytes += runtime_datatypes().storage_bytes(type, elements);
        const size_t bits = runtime_datatypes().storage_bits(type, elements);
        const size_t glb_accesses = (bits + glb->line_size[type] - 1)/glb->line_size[type];
        const double feed_cycle = static_cast<double>(glb_accesses)*glb->u_read_cycle[type];
        stream.glb_access_cycle += feed_cycle;
        stream.glb_access_energy += static_cast<double>(glb_accesses)*glb->u_read_energy[type];
        if(m_plan.inputs[index] == WORKLOAD_RESIDENCY_GLB) {
            ++resident_inputs;
            stream.ingress_cycle += feed_cycle;
            continue;
        }
        const datatype_transfer_timing_t timing = datatype_transfer_timing(
            type, elements, dram->line_size[type], glb->line_size[type], dram->get_bitwidth());
        stream.dram_access_cycle += static_cast<double>(timing.source_accesses)*dram->u_read_cycle[type];
        stream.dram_access_energy += static_cast<double>(timing.source_accesses)*dram->u_read_energy[type];
        stream.dram_link_cycle += static_cast<double>(timing.link_transactions)*dram->u_transfer_cycle;
        stream.dram_link_energy += static_cast<double>(timing.link_transactions)*dram->u_transfer_energy;
        stream.dram_link_transactions += timing.link_transactions;
        stream.glb_access_cycle += static_cast<double>(timing.destination_accesses)*glb->u_write_cycle[type];
        stream.glb_access_energy += static_cast<double>(timing.destination_accesses)*glb->u_write_energy[type];
        double row_cycle = 0.0;
        if(dram->row_buffer_bytes > 0) {
            const dram_row_activation_cost_t rows = dram_row_activations(
                runtime_datatypes().storage_bytes(type, elements), dram->row_buffer_bytes,
                dram->num_banks);
            const double unit_cycle = (dram->t_ras_cycle > 0.0 && dram->t_rp_cycle > 0.0)
                ? dram->t_ras_cycle + dram->t_rp_cycle : dram->u_row_miss_cycle;
            row_cycle = static_cast<double>(rows.busiest_bank)*unit_cycle;
            stream.dram_row_activations += rows.activations;
            stream.dram_row_activation_cycle += row_cycle;
            stream.dram_row_activation_energy += static_cast<double>(rows.activations)*dram->u_row_miss_energy;
        }
        stream.ingress_cycle += pipelined_transfer_cycles(
            timing.groups, dram->u_read_cycle[type], dram->u_transfer_cycle,
            glb->u_write_cycle[type]) + feed_cycle + row_cycle;
    }

    const workload_tensor_t &output = workload->tensor(m_operation.outputs.front());
    const size_t output_elements = output.elements();
    stream.egress_bytes = runtime_datatypes().storage_bytes(data_type_t::OUTPUT, output_elements);
    const size_t output_bits = runtime_datatypes().storage_bits(data_type_t::OUTPUT, output_elements);
    const size_t output_glb_accesses =
        (output_bits + glb->line_size[data_type_t::OUTPUT] - 1)/glb->line_size[data_type_t::OUTPUT];
    const double result_cycle = static_cast<double>(output_glb_accesses)*
                                glb->u_write_cycle[data_type_t::OUTPUT];
    stream.glb_access_cycle += result_cycle;
    stream.glb_access_energy += static_cast<double>(output_glb_accesses)*
                                glb->u_write_energy[data_type_t::OUTPUT];
    stream.egress_cycle = result_cycle;
    if(!m_plan.retain_output) {
        const datatype_transfer_timing_t timing = datatype_transfer_timing(
            data_type_t::OUTPUT, output_elements, glb->line_size[data_type_t::OUTPUT],
            dram->line_size[data_type_t::OUTPUT], dram->get_bitwidth());
        stream.dram_access_cycle += static_cast<double>(timing.destination_accesses)*
                                    dram->u_write_cycle[data_type_t::OUTPUT];
        stream.dram_access_energy += static_cast<double>(timing.destination_accesses)*
                                     dram->u_write_energy[data_type_t::OUTPUT];
        stream.dram_link_cycle += static_cast<double>(timing.link_transactions)*dram->u_transfer_cycle;
        stream.dram_link_energy += static_cast<double>(timing.link_transactions)*dram->u_transfer_energy;
        stream.dram_link_transactions += timing.link_transactions;
        stream.glb_access_cycle += static_cast<double>(timing.source_accesses)*
                                   glb->u_read_cycle[data_type_t::OUTPUT];
        stream.glb_access_energy += static_cast<double>(timing.source_accesses)*
                                    glb->u_read_energy[data_type_t::OUTPUT];
        double row_cycle = 0.0;
        if(dram->row_buffer_bytes > 0) {
            const dram_row_activation_cost_t rows = dram_row_activations(
                stream.egress_bytes, dram->row_buffer_bytes, dram->num_banks);
            const double unit_cycle = (dram->t_ras_cycle > 0.0 && dram->t_rp_cycle > 0.0)
                ? dram->t_ras_cycle + dram->t_rp_cycle : dram->u_row_miss_cycle;
            row_cycle = static_cast<double>(rows.busiest_bank)*unit_cycle;
            stream.dram_row_activations += rows.activations;
            stream.dram_row_activation_cycle += row_cycle;
            stream.dram_row_activation_energy += static_cast<double>(rows.activations)*dram->u_row_miss_energy;
        }
        stream.egress_cycle += pipelined_transfer_cycles(
            timing.groups, glb->u_read_cycle[data_type_t::OUTPUT], dram->u_transfer_cycle,
            dram->u_write_cycle[data_type_t::OUTPUT]) + row_cycle;
    }
    const size_t retained_inputs = static_cast<size_t>(std::count(
        m_plan.retain_inputs.begin(), m_plan.retain_inputs.end(), true));
    stream.residency = std::to_string(resident_inputs) + "/" +
        std::to_string(m_operation.inputs.size()) + " input tensor(s) in GLB; output " +
        (m_plan.retain_output ? "retained in GLB" : "materialized to DRAM") +
        (retained_inputs ? "; " + std::to_string(retained_inputs) +
         " future-use input(s) pinned" : "");
    return stream;
}

void npu_t::run_standalone_graph_operation(
    unsigned m_index, workload_residency_plan_t m_plan,
    const std::string &m_accelerator_config, const std::string &m_network_config) {
    if(sfus.empty()) {
        std::cerr << "Error: executable operation " << workload->operations[m_index].id
                  << " requires an [sfu] section for non-MAC timing" << std::endl;
        exit(1);
    }
    const workload_operation_t &operation = workload->operations[m_index];
    const size_t output_elements = workload->tensor(operation.outputs.front()).elements();
    sfu_invocation_t combined;
    const size_t base = output_elements/num_processors;
    const size_t remainder = output_elements % num_processors;
    size_t output_begin = 0;
    for(unsigned chip = 0; chip < num_processors; ++chip) {
        const size_t share = base + (chip < remainder ? 1 : 0);
        sfu_invocation_t local;
        if(operation.kind == WORKLOAD_SOFTMAX) {
            const size_t row_base = operation.geometry.rows/num_processors;
            const size_t row_remainder = operation.geometry.rows % num_processors;
            local = sfus[chip]->softmax_invocation(
                row_base + (chip < row_remainder ? 1 : 0), operation.geometry.row_length);
        } else if(operation.kind == WORKLOAD_POOL2D) {
            const size_t reductions =
                pool_reduction_operations(operation.geometry, output_begin, share);
            if(operation.geometry.mode == "max") {
                local = sfus[chip]->elementwise_invocation(SFU_OP_VMAX, reductions);
            } else {
                local = sfus[chip]->elementwise_invocation(SFU_OP_VADD, reductions);
                local.merge_serial(sfus[chip]->elementwise_invocation(SFU_OP_VMUL, share));
            }
        } else if(operation.kind == WORKLOAD_ELEMENTWISE) {
            const sfu_op_t op = operation.geometry.elementwise_operator == "add"
                ? SFU_OP_VADD : SFU_OP_VMUL;
            local = sfus[chip]->elementwise_invocation(op, share);
        } else if(operation.kind == WORKLOAD_BATCH_NORM) {
            local = sfus[chip]->elementwise_invocation(SFU_OP_VMUL, share);
            local.merge_serial(sfus[chip]->elementwise_invocation(SFU_OP_VADD, share));
        } else if(operation.kind == WORKLOAD_CONCAT) {
            local.operation = "concat (copy-only)";
        } else {
            std::cerr << "Error: unsupported standalone graph operation " << operation.id << std::endl;
            exit(1);
        }
        combined.merge_parallel(local);
        output_begin += share;
    }
    combined.valid_elements = output_elements;
    if(operation.kind == WORKLOAD_POOL2D) {
        combined.operation = operation.geometry.mode + " pool2d";
    } else if(operation.kind == WORKLOAD_ELEMENTWISE) {
        combined.operation = "elementwise " + operation.geometry.elementwise_operator;
    } else if(operation.kind == WORKLOAD_BATCH_NORM) {
        combined.operation = "batch_norm inference (mul-add)";
    } else if(operation.kind == WORKLOAD_CONCAT) {
        combined.operation = "concat (copy-only)";
    }

    const sfu_operand_stream_t stream = graph_operand_stream(operation, m_plan);
    workload_lifetime->commit(m_index, &m_plan);
    stats_t *stats = new stats_t();
    sfu_layer_stats.push_back(stats);
    const double frequencies[5] = {
        pe_arrays[0]->pes[0]->clock_mhz(), pe_arrays[0]->clock_mhz(),
        global_buffers[0]->clock_mhz(), multi_chip->clock_mhz(), dram->clock_mhz()};
    bool single_clock = frequencies[0] > 0.0;
    for(unsigned f = 1; f < 5; ++f) if(frequencies[f] != frequencies[0]) single_clock = false;
    const std::string clock_note = single_clock
        ? "all modeled components share one clock"
        : "mixed or undeclared clock domains; cycles are not convertible to time";
    stats->sfu_contract_note = "framework-neutral DAG operation " + operation.id +
        "; memory stream uses tensor liveness/residency; parameter traffic is aggregated"
        " on the standalone stream rows";
    stats->sfu_profile_reference = sfus[0]->get_profile_reference();
    stats->record_sfu_only_layer(combined,
        num_processors*sfus[0]->get_num_units(), sfus[0]->get_lanes(),
        sfus[0]->get_static_energy_per_cycle(), single_clock ? frequencies[0] : 0.0,
        single_clock, clock_note, stream);
    const std::string residency_note = stream.residency + "; GLB occupancy " +
        std::to_string(m_plan.occupied_before) + " -> " +
        std::to_string(m_plan.occupied_after) + " / " +
        std::to_string(m_plan.capacity) + " bytes";
    stats->apply_graph_residency(false, false, residency_note);

    std::cout << "The simulation result of #" << m_index << " operation "
              << operation.id << std::endl << std::endl;
    const std::string output_file_name = m_accelerator_config + "_" + m_network_config +
        "_layer_" + std::to_string(m_index) + ".txt";
    std::ofstream output_file(output_file_name.c_str(), std::ios::out);
    print_workload_provenance(output_file);
    output_file << "Executable DAG operation: " << workload->operation_kind_name(operation.kind)
                << " (" << operation.id << ")" << std::endl << std::endl;
    stats->print_results(output_file);
    output_file.close();
    network_stats->update_network_stats(stats);
}

// Phase-7 (plan_sfu.md): standalone softmax on the SFU multi-pass microprogram
// (max -> subtract -> exp -> sum -> reciprocal -> normalize). The layer has no mapping
// section; rows distribute across every chip's SFU (independent rows, latency = busiest
// chip), and the operand tensor's streaming between the memory hierarchy and the SFU is
// charged per [sfu] softmax_operand_residency.
void npu_t::run_standalone_softmax(unsigned m_index, const std::string &m_accelerator_config,
                                   const std::string &m_network_config) {
    nebula::layer_t *softmax_layer = network->layers[m_index];
    // Softmax groups (Darknet/Nebula semantics): the vector splits into `groups`
    // independent normalization spans, so the SFU sees batch x groups rows of
    // output_size/groups elements. Nebula's init fail-fasts on non-dividing groups;
    // the guard here keeps the cost model safe against a frontend regression.
    const size_t groups = std::max(1u, softmax_layer->group);
    if(softmax_layer->output_size % groups != 0) {
        std::cerr << "Error: softmax layer " << m_index << " groups = " << groups
                  << " does not divide output size " << softmax_layer->output_size
                  << std::endl;
        exit(1);
    }
    const size_t rows = static_cast<size_t>(network->batch_size)*groups;
    const size_t row_length = softmax_layer->output_size/groups;

    std::cout << "The network configuration of #" << m_index
              << " layer (standalone softmax on SFU)" << std::endl;
    std::cout << " * softmax rows x length : " << rows << " x " << row_length
              << "  (batch " << network->batch_size << " x groups " << groups << ")"
              << std::endl;

    // Rows are independent, so they partition across every physical chip's SFU; the
    // window follows the busiest chip and work/energy/traffic sum (merge_parallel).
    sfu_invocation_t invocation;
    const size_t base_rows = rows/num_processors;
    const size_t remainder_rows = rows % num_processors;
    for(unsigned c = 0; c < num_processors; ++c) {
        const size_t share = base_rows + (c < remainder_rows ? 1 : 0);
        invocation.merge_parallel(sfus[c]->softmax_invocation(share, row_length));
    }

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

    // Phase-7: operand streaming. The producing layer committed the tensor at OUTPUT
    // precision; the softmax result returns the same way. Residency decides the walk:
    //   dram : DRAM -> GLB staging -> SFU, and the mirror on egress (matches this
    //          simulator's layer flow, which materializes every layer's output off-chip);
    //   glb  : the tensor is retained on-chip by a fused schedule -- only the GLB feed
    //          and result ports are exercised (the tensor must fit).
    const sfu_operand_stream_t stream = softmax_operand_stream(rows*row_length);
    stats->sfu_contract_note =
        "standalone softmax distributes " + std::to_string(rows) + " row(s) (batch x"
        " groups) across " + std::to_string(num_processors) + " chip SFU(s); operand"
        " residency: " + stream.residency +
        (stream.dram_row_activations > 0
             ? "; " + std::to_string(stream.dram_row_activations) +
               " DRAM row activation(s) charged on the transfer axis"
             : "") +
        " (the NoP hop cost is not modeled for this layer)";
    if(!sfus[0]->get_precision_note().empty()) {
        stats->sfu_contract_note += "; " + sfus[0]->get_precision_note();
    }

    stats->sfu_profile_reference = sfus[0]->get_profile_reference();
    stats->record_sfu_only_layer(invocation,
        num_processors*sfus[0]->get_num_units(), sfus[0]->get_lanes(),
        sfus[0]->get_static_energy_per_cycle(),
        single_clock ? frequencies[0] : 0.0, single_clock, clock_note, stream);

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

// Phase-7 (plan_sfu.md): operand-tensor streaming cost for a standalone softmax, at
// OUTPUT precision, from the live components' declared unit costs. Nothing here mutates
// a component counter -- the costs land in the softmax layer's own stats object.
sfu_operand_stream_t npu_t::softmax_operand_stream(size_t m_elements) {
    sfu_operand_stream_t stream;
    stream.active = true;
    const std::string &residency = sfus[0]->get_softmax_operand_residency();
    global_buffer_t *glb = global_buffers[0];
    const size_t tensor_bytes = runtime_datatypes().storage_bytes(data_type_t::OUTPUT,
                                                                  m_elements);
    stream.ingress_bytes = tensor_bytes;
    stream.egress_bytes = tensor_bytes;

    // GLB ports common to both residencies: feed the SFU (read) and take back the result
    // (write), one access per GLB output line.
    const size_t storage_bits = runtime_datatypes().storage_bits(data_type_t::OUTPUT,
                                                                 m_elements);
    const size_t glb_line_bits = std::max(1u, glb->line_size[data_type_t::OUTPUT]);
    const size_t glb_port_accesses = (storage_bits + glb_line_bits - 1)/glb_line_bits;
    const double feed_read_cycle = static_cast<double>(glb_port_accesses)*
                                   glb->u_read_cycle[data_type_t::OUTPUT];
    const double result_write_cycle = static_cast<double>(glb_port_accesses)*
                                      glb->u_write_cycle[data_type_t::OUTPUT];
    stream.glb_access_cycle = feed_read_cycle + result_write_cycle;
    stream.glb_access_energy =
        static_cast<double>(glb_port_accesses)*glb->u_read_energy[data_type_t::OUTPUT] +
        static_cast<double>(glb_port_accesses)*glb->u_write_energy[data_type_t::OUTPUT];

    if(residency == "glb") {
        // On-chip retained (fused-schedule scenario): the tensor never leaves the chip,
        // so it must FIT -- input operand and result at once.
        const double required_bytes = 2.0*static_cast<double>(tensor_bytes);
        if(required_bytes > glb->get_buffer_size()) {
            std::cerr << "Error: [sfu] softmax_operand_residency = glb, but the softmax"
                      << " tensor needs " << required_bytes << " bytes (operand + result)"
                      << " and the global buffer holds " << glb->get_buffer_size()
                      << "; use softmax_operand_residency = dram or a larger buffer"
                      << std::endl;
            exit(1);
        }
        stream.residency = "glb (on-chip retained; fits the global buffer)";
        stream.ingress_cycle = feed_read_cycle;
        stream.egress_cycle = result_write_cycle;
        return stream;
    }

    // Materialized round trip (default; matches the simulator's layer flow -- the
    // producing layer committed the tensor off-chip): DRAM -> GLB staging -> SFU on
    // ingress, and the mirror on egress. Each hop is the standard three-resource
    // transfer: source access, link crossing, destination access.
    stream.residency = "dram (materialized round trip via GLB staging)";
    const datatype_transfer_timing_t ingress = datatype_transfer_timing(
        data_type_t::OUTPUT, m_elements, dram->line_size[data_type_t::OUTPUT],
        glb->line_size[data_type_t::OUTPUT], dram->get_bitwidth());
    const datatype_transfer_timing_t egress = datatype_transfer_timing(
        data_type_t::OUTPUT, m_elements, glb->line_size[data_type_t::OUTPUT],
        dram->line_size[data_type_t::OUTPUT], dram->get_bitwidth());

    stream.dram_access_cycle =
        static_cast<double>(ingress.source_accesses)*dram->u_read_cycle[data_type_t::OUTPUT] +
        static_cast<double>(egress.destination_accesses)*dram->u_write_cycle[data_type_t::OUTPUT];
    stream.dram_access_energy =
        static_cast<double>(ingress.source_accesses)*dram->u_read_energy[data_type_t::OUTPUT] +
        static_cast<double>(egress.destination_accesses)*dram->u_write_energy[data_type_t::OUTPUT];
    stream.dram_link_cycle = dram->u_transfer_cycle*
        static_cast<double>(ingress.link_transactions + egress.link_transactions);
    stream.dram_link_energy = dram->u_transfer_energy*
        static_cast<double>(ingress.link_transactions + egress.link_transactions);
    stream.dram_link_transactions = ingress.link_transactions + egress.link_transactions;
    // Open-page row activations of the two sequential streams -- the SAME model and cost
    // resolution dram_t applies to its own streams (tRC when calibrated, else the flat
    // row_miss cost; bank parallelism helps latency, never energy). Disabled, exactly
    // like dram_t, when no row-buffer geometry is declared.
    if(dram->row_buffer_bytes > 0) {
        const double per_activation_cycle =
            (dram->t_ras_cycle > 0.0 && dram->t_rp_cycle > 0.0)
                ? dram->t_ras_cycle + dram->t_rp_cycle : dram->u_row_miss_cycle;
        const dram_row_activation_cost_t ingress_rows = dram_row_activations(
            tensor_bytes, dram->row_buffer_bytes, dram->num_banks);
        const dram_row_activation_cost_t egress_rows = dram_row_activations(
            tensor_bytes, dram->row_buffer_bytes, dram->num_banks);
        stream.dram_row_activations = ingress_rows.activations + egress_rows.activations;
        stream.dram_row_activation_cycle =
            static_cast<double>(ingress_rows.busiest_bank + egress_rows.busiest_bank)*
            per_activation_cycle;
        stream.dram_row_activation_energy =
            static_cast<double>(stream.dram_row_activations)*dram->u_row_miss_energy;
    }
    // GLB staging ports: the off-chip transfers land in (and drain from) the GLB.
    stream.glb_access_cycle +=
        static_cast<double>(ingress.destination_accesses)*glb->u_write_cycle[data_type_t::OUTPUT] +
        static_cast<double>(egress.source_accesses)*glb->u_read_cycle[data_type_t::OUTPUT];
    stream.glb_access_energy +=
        static_cast<double>(ingress.destination_accesses)*glb->u_write_energy[data_type_t::OUTPUT] +
        static_cast<double>(egress.source_accesses)*glb->u_read_energy[data_type_t::OUTPUT];

    // Serial makespans on the critical path: the off-chip hop pipelines internally
    // (packet-level source/link/destination overlap), then the GLB->SFU feed runs. Each
    // direction additionally pays its own stream's busiest-bank row-activation
    // serialization (split evenly: the two directions activate the same row count).
    stream.ingress_cycle = pipelined_transfer_cycles(ingress.groups,
        dram->u_read_cycle[data_type_t::OUTPUT], dram->u_transfer_cycle,
        glb->u_write_cycle[data_type_t::OUTPUT]) + feed_read_cycle +
        stream.dram_row_activation_cycle/2.0;
    stream.egress_cycle = result_write_cycle + pipelined_transfer_cycles(egress.groups,
        glb->u_read_cycle[data_type_t::OUTPUT], dram->u_transfer_cycle,
        dram->u_write_cycle[data_type_t::OUTPUT]) +
        stream.dram_row_activation_cycle/2.0;
    return stream;
}
