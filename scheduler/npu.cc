#include <thread>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <cstdio>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
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

#ifdef FUNCTIONAL
// Round an fp32 value to the bf16 grid (round-to-nearest-even): keep the top 16 bits with a
// rounding bias. Matches the Python golden's identical bit manipulation.
inline float round_bf16(float f) {
    uint32_t x; std::memcpy(&x, &f, sizeof(x));
    if((x & 0x7fffffffu) > 0x7f800000u) return f;                 // NaN: leave as-is
    x += 0x7fffu + ((x >> 16) & 1u);                              // round-to-nearest-even
    x &= 0xffff0000u;
    float r; std::memcpy(&r, &x, sizeof(r)); return r;
}
// Round to the IEEE fp16 grid (round-to-nearest-even) via the compiler half type; matches
// Python struct 'e'. _Float16 arithmetic type is available on the x86-64 GCC used here.
inline float round_fp16(float f) { return static_cast<float>(static_cast<_Float16>(f)); }
inline float round_lowp(const std::string &fmt, float f) {
    if(fmt == "bf16") return round_bf16(f);
    if(fmt == "fp16") return round_fp16(f);
    return f;
}
#endif

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
#ifdef FUNCTIONAL
    functional_layers_checked = 0;
    functional_layers_failed = 0;
    functional_external_golden = false;
    functional_last_mapped = -1;
    functional_requant_shift = 0;
    functional_requant_min = -127;
    functional_requant_max = 127;
    functional_input_zero_point = 0;
    functional_weight_zero_point = 0;
#endif
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
#ifdef FUNCTIONAL
        // G1: a FUNCTIONAL executable-IR run needs the npusim.tensor.v1 value artifact.
        // Load + validate it (hash-bound to this executable), seed the tensor store with
        // the graph inputs and parameters, and arm the acceptance gate.
        if(!functional_artifact_path.empty()) {
            // G2 for executable runs: the artifact carries float32 values, so the
            // accelerator must declare fp32 tensors (no override channel here -- pick an
            // fp32 accelerator config for functional executable runs).
            const bool fp32_formats =
                runtime_datatypes().format(data_type_t::INPUT).kind  == data_format_kind_t::FP32 &&
                runtime_datatypes().format(data_type_t::WEIGHT).kind == data_format_kind_t::FP32 &&
                runtime_datatypes().format(data_type_t::OUTPUT).kind == data_format_kind_t::FP32;
            if(!fp32_formats) {
                std::cerr << "Error: functional executable runs carry float32 values; the"
                          << " accelerator config must declare input/weight/output_format ="
                          << " fp32" << std::endl;
                exit(1);
            }
            functional_artifact.load(functional_artifact_path, *workload);
            for(const auto &entry : functional_artifact.tensors) {
                functional_store(entry.first) = entry.second;
            }
            functional_semantics = "fp32";
            functional_external_golden = true;
        }
#endif
    } else {
        network->init(m_network_config);
    }
#ifdef FUNCTIONAL
    // Functional simulation needs REAL parameters: nebula's network init leaves
    // init_weight()/init_data() commented out, so layer->weight/bias stay zero and both the
    // accelerator datapath (which reads layer->weight) and the reference forward() would
    // compute all-zero. Load them here for the legacy path. init_weight() takes the WEIGHT
    // FILE path (not the network config), so parse [data] weight from the config and resolve
    // it relative to the config directory (nebula opens it via the same relative convention).
    // (The executable-IR path will instead receive parameters through the tensor artifact --
    // see the correctness plan.)
    if(!executable_ir_mode) {
        config_t functional_config;
        functional_config.parse(m_network_config);
        std::string weight_path;
        for(unsigned i = 0; i < functional_config.sections.size(); i++) {
            if(functional_config.sections[i].get_setting("weight", &weight_path) &&
               !weight_path.empty()) {
                break;
            }
        }
        if(weight_path.empty()) {
            std::cerr << "Error: FUNCTIONAL build requires a [data] weight file in "
                      << m_network_config << std::endl;
            exit(1);
        }
        network->init_weight(weight_path);

        // Fixture-driven input + external golden (correctness plan Phase 1, lite). When
        // [data] functional_input / functional_golden name raw fp32 files, load them and
        // bypass the nebula image loader + forward() oracle. A helper reads a whole file of
        // little-endian float32.
        std::string finput_path, fgolden_path, rmult_path;
        for(unsigned i = 0; i < functional_config.sections.size(); i++) {
            functional_config.sections[i].get_setting("functional_input", &finput_path);
            functional_config.sections[i].get_setting("functional_golden", &fgolden_path);
            // INT8 requantization parameters (optional, [data] section).
            functional_config.sections[i].get_setting("requant_shift", &functional_requant_shift);
            functional_config.sections[i].get_setting("requant_min",   &functional_requant_min);
            functional_config.sections[i].get_setting("requant_max",   &functional_requant_max);
            functional_config.sections[i].get_setting("input_zero_point", &functional_input_zero_point);
            functional_config.sections[i].get_setting("weight_zero_point", &functional_weight_zero_point);
            functional_config.sections[i].get_setting("requant_mult", &rmult_path);   // per-channel
            functional_config.sections[i].get_setting("output_format", &functional_output_format); // fp16/bf16
            functional_config.sections[i].get_setting("weight_layout", &functional_weight_layout); // ""|ktile
        }
        auto read_fp32 = [](const std::string &m_path, std::vector<float> *m_out) {
            std::ifstream in(m_path.c_str(), std::ios::binary | std::ios::ate);
            if(!in.is_open()) {
                std::cerr << "Error: cannot open functional fixture " << m_path << std::endl;
                exit(1);
            }
            const std::streamsize bytes = in.tellg();
            if(bytes % static_cast<std::streamsize>(sizeof(float)) != 0) {
                std::cerr << "Error: fixture " << m_path << " size is not a float multiple"
                          << std::endl;
                exit(1);
            }
            in.seekg(0);
            m_out->resize(static_cast<size_t>(bytes)/sizeof(float));
            in.read(reinterpret_cast<char*>(m_out->data()), bytes);
        };
        if(!finput_path.empty() && !fgolden_path.empty()) {
            read_fp32(finput_path, &functional_input_buffer);
            read_fp32(fgolden_path, &functional_golden);
            functional_external_golden = true;
        }
        if(!rmult_path.empty()) read_fp32(rmult_path, &functional_requant_mult);  // per-channel

        // G4 (gaps plan Step 3): optional per-layer goldens -- [data] functional_golden<i>
        // (post-finalize) and functional_golden_raw<i> (raw accumulator, pre-finalize) --
        // localize the first mismatching operation of a multi-layer DAG.
        if(functional_external_golden) {
            for(unsigned l = 0; l < network->num_layers; ++l) {
                std::string layer_path, raw_path;
                for(unsigned i = 0; i < functional_config.sections.size(); i++) {
                    functional_config.sections[i].get_setting(
                        "functional_golden" + std::to_string(l), &layer_path);
                    functional_config.sections[i].get_setting(
                        "functional_golden_raw" + std::to_string(l), &raw_path);
                }
                if(!layer_path.empty()) read_fp32(layer_path, &functional_layer_golden[l]);
                if(!raw_path.empty())   read_fp32(raw_path, &functional_layer_golden_raw[l]);
            }
        }

        // G2 (gaps plan Step 2): bind the accelerator config's declared tensor formats to
        // the functional arithmetic. The fixture states its semantics ([data]
        // functional_semantics, or inferred from its quantization/rounding knobs); the
        // accelerator states its formats (input/weight/output_format, parsed into
        // runtime_datatypes() during init). A disagreement -- e.g. an int8-declared config
        // computing fp32 values -- is refused before simulation, unless the fixture
        // explicitly opts into an fp32 reference run with functional_semantics =
        // fp32_reference (recorded as such in the report).
        {
            std::string declared;
            for(unsigned i = 0; i < functional_config.sections.size(); i++) {
                functional_config.sections[i].get_setting("functional_semantics", &declared);
            }
            if(declared.empty()) {
                if(functional_output_format == "fp16" || functional_output_format == "bf16")
                    declared = functional_output_format;
                else if(functional_requant_shift > 0 || functional_input_zero_point != 0 ||
                        functional_weight_zero_point != 0 || !functional_requant_mult.empty())
                    declared = "int8";
                else
                    declared = "fp32";
            }
            auto semantics_of = [](const tensor_format_t &f) -> std::string {
                switch(f.kind) {
                    case data_format_kind_t::INT:
                    case data_format_kind_t::UINT: return f.payload_bits == 8 ? "int8" : "unsupported";
                    case data_format_kind_t::FP16: return "fp16";
                    case data_format_kind_t::BF16: return "bf16";
                    case data_format_kind_t::FP32: return "fp32";
                    default:                       return "unsupported";
                }
            };
            const std::string in_sem  = semantics_of(runtime_datatypes().format(data_type_t::INPUT));
            const std::string wt_sem  = semantics_of(runtime_datatypes().format(data_type_t::WEIGHT));
            const std::string out_sem = semantics_of(runtime_datatypes().format(data_type_t::OUTPUT));
            if(declared == "fp32_reference") {
                functional_semantics = "fp32_reference(" + out_sem + ")";
            } else {
                if(in_sem != out_sem || wt_sem != out_sem) {
                    std::cerr << "Error: functional simulation needs matching input/weight/"
                              << "output_format declarations (got " << in_sem << "/" << wt_sem
                              << "/" << out_sem << "); declare [data] functional_semantics = "
                              << "fp32_reference to run an fp32 reference anyway" << std::endl;
                    exit(1);
                }
                if(out_sem == "unsupported") {
                    std::cerr << "Error: the accelerator's declared tensor format has no "
                              << "functional arithmetic yet; declare [data] "
                              << "functional_semantics = fp32_reference for an fp32 reference run"
                              << std::endl;
                    exit(1);
                }
                if(declared != out_sem) {
                    std::cerr << "Error: accelerator declares " << out_sem << " tensors "
                              << "(input/weight/output_format) but the fixture's functional "
                              << "semantics is '" << declared << "'; run it on a matching "
                              << "accelerator config or declare [data] functional_semantics = "
                              << "fp32_reference for an explicit fp32 reference run" << std::endl;
                    exit(1);
                }
                functional_semantics = declared;
                // Config-driven binding: an fp16/bf16 accelerator rounds the finalized
                // output to its grid even when the fixture omits output_format.
                if((out_sem == "fp16" || out_sem == "bf16") && functional_output_format.empty())
                    functional_output_format = out_sem;
            }
        }
    }
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
        if(executable_ir_mode && functional_artifact_path.empty()) {
            std::cerr << "Error: a FUNCTIONAL executable-IR run requires the npusim.tensor.v1"
                      << " value artifact; use run-ir-functional <accelerator> <executable>"
                      << " <mapping> <tensors.json>" << std::endl;
            exit(1);
        }
#endif
#ifdef FUNCTIONAL
        // With an external-golden fixture the image loader's output is fully replaced by
        // functional_input_buffer, so skip it: it would only constrain fixtures to image-
        // shaped inputs (e.g. it exits on channel counts other than 1/3) and fill accuracy
        // bookkeeping this run never reads.
        const bool skip_image_loader = functional_external_golden;
#else
        const bool skip_image_loader = false;
#endif
        if(!executable_ir_mode && !skip_image_loader) network->load_data(iteration);
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
#ifdef FUNCTIONAL
                // Fixture input: point layer 0 at our OWNED input buffer, fully bypassing
                // nebula's image loader. The connected/conv layer's input_data is a borrowed
                // pointer (never freed by the layer), so re-pointing it is safe. Size must
                // match the layer's expected input footprint.
                if(index == 0 && functional_external_golden) {
                    const size_t n = static_cast<size_t>(network->layers[0]->input_size)*
                                     network->batch_size;
                    if(functional_input_buffer.size() != n) {
                        std::cerr << "Error: functional_input has " << functional_input_buffer.size()
                                  << " floats, layer 0 expects " << n
                                  << " (input_size " << network->layers[0]->input_size
                                  << " x batch " << network->batch_size << ")" << std::endl;
                        exit(1);
                    }
                    network->layers[0]->input_data = functional_input_buffer.data();
                }
#endif
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
                multi_chip->connect_layer(layer);
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

#ifdef FUNCTIONAL
                if(executable_ir_mode) functional_bind_executable_operation(index, *operation);
                functional_reject_unsupported_mapping(index);
#endif
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
                    pe_arrays[chip]->flush_psum_writeback(scheduler);
                }
                // Functional output write-back chain, in order: PE array -> GLB (above) ->
                // multi-chip (here) -> DRAM/layer tensor (below). Each level flushes the
                // last retained output tile the eviction loop never wrote back.
                for(unsigned i = 0; i < global_buffers.size(); i++) {
                    global_buffers[i]->flush_output_writeback(scheduler);
                }
                multi_chip->flush_output_writeback(scheduler);
                for(unsigned chip = 0; chip < pe_arrays.size(); chip++) {
                    for(unsigned pe = 0; pe < pe_arrays[chip]->get_number_of_pes(); pe++) {
                        pe_arrays[chip]->pes[pe]->suppress_streaming_cycles();
                    }
                }
                layer_stats[stats_index]->update_stats(pe_arrays, global_buffers, multi_chip, dram);
#ifdef FUNCTIONAL
                // G3: a convolution with P/Q > 1 gets its VALUES from the in-simulator
                // im2col kernel (the datapath above already produced this layer's timing);
                // the temporal-fold replay below is GEMM machinery and is skipped for it.
                const bool conv_value_kernel =
                    scheduler->layer_name == layer_name_t::CONVOLUTIONAL_LAYER &&
                    (layer->output_height > 1 || layer->output_width > 1);
                if(conv_value_kernel) functional_conv_im2col(index);
                // TEMPORAL-FOLD FUNCTIONAL REPLAY (batch): the analytical engine simulates ONE
                // representative tile and scales timing by repetitions, so only the first
                // batch's VALUES were just computed. Timing is already captured (update_stats
                // above) and later multiplied by scale_serial_repetitions, so re-running the
                // datapath here changes ONLY layer->output_data. We slide the layer's input and
                // output windows by one sample per batch and replay the exact same datapath +
                // write-back chain, so every batch's outputs land in their own tensor region.
                // Restricted to a pure-batch temporal fold (weight reused; no other GLB/DRAM
                // output fold); other folds fall through unreplayed (see functional-sim plan).
                if(!conv_value_kernel) {
                    // Per-dimension GLB TEMPORAL-REPETITION fold count. The GLB row is stored as
                    // "legacy GLB temporal-repetition factors" SEPARATE from the mapping-table
                    // cumulative product, so calculate_parameter_size() misses it. The full
                    // per-dimension coverage (calculate_total_parameter_size) equals the DRAM
                    // cumulative product TIMES those GLB factors, hence fold = full / spatial.
                    // (DRAM-queue folds live in the spatial product here but are excluded by the
                    // size-1 DRAM-queue guard below, so `fold` is exactly the GLB repetition.)
                    mapping_table_t *mt = scheduler->mapping_table;
                    const std::vector<unsigned> full    = mt->calculate_total_parameter_size();
                    const std::vector<unsigned> spatial = mt->calculate_parameter_size(component_type_t::DRAM);
                    auto fold = [&](parameter_type_t d) -> unsigned {
                        return spatial[d] ? full[d]/spatial[d] : 1;
                    };
                    const unsigned Bf = fold(parameter_type_t::BATCH_SIZE);
                    const unsigned Nf = fold(parameter_type_t::OUTPUT_CHANNEL);
                    const unsigned Pf = fold(parameter_type_t::OUTPUT_HEIGHT);
                    const unsigned Qf = fold(parameter_type_t::OUTPUT_WIDTH);
                    const unsigned Kf = fold(parameter_type_t::INPUT_CHANNEL);
                    // One functional datapath pass at the current layer input/weight/output
                    // windows, followed by the output write-back chain (PE psum -> GLB -> MC ->
                    // DRAM/layer tensor). Timing was already captured at update_stats, so extra
                    // passes only move values.
                    auto run_pass = [&]() {
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
                        for(unsigned chip = 0; chip < pe_arrays.size(); chip++)
                            pe_arrays[chip]->flush_psum_writeback(scheduler);
                        for(unsigned i = 0; i < global_buffers.size(); i++)
                            global_buffers[i]->flush_output_writeback(scheduler);
                        multi_chip->flush_output_writeback(scheduler);
                    };
                    // UNIFIED temporal-fold replay over the distinct-output tiles (batch b,
                    // conv output positions p/q, output-channel block nf) with an inner REDUCTION
                    // fold (kt) that ACCUMULATES. Timing was captured at update_stats, so this only
                    // moves values.
                    //  - distinct outputs (b,p,q,nf): each names its own output region; scatter.
                    //  - reduction (kt): every K-tile adds a PARTIAL to the SAME output; accumulate
                    //    (kt=0 writes, kt>0 adds a scratch tile).
                    // Layouts: input [M][K] (row b, conv window p/q, K-tile kt at kt*sK -- all
                    // contiguous); weight [Kf][N][sK] reduction-tile-major so the array's spatial_N
                    // channels stay contiguous at kt*full_N*sK + nf*spatial_N*sK (for Kf==1,
                    // sK==K_full and this is the plain [N][K] layout -- existing fixtures unchanged);
                    // output POSITION-MAJOR [B][P][Q][N] (P/Q folded => the per-pass tile has
                    // OH=OW=1, so spatial_N channels are written contiguously). Guards: no filter
                    // fold, P/Q fully folded (spatial P=Q=1), a reduction fold only for GEMM
                    // (P=Q=1), and size-1 DRAM queues (the DRAM-queue mechanism is out of scope).
                    const unsigned full_N = full[parameter_type_t::OUTPUT_CHANNEL];
                    const unsigned full_P = full[parameter_type_t::OUTPUT_HEIGHT];
                    const unsigned full_Q = full[parameter_type_t::OUTPUT_WIDTH];
                    const bool filter_fold =
                        fold(parameter_type_t::FILTER_HEIGHT) != 1 ||
                        fold(parameter_type_t::FILTER_WIDTH)  != 1;
                    const bool pq_partly_spatial =
                        spatial[parameter_type_t::OUTPUT_HEIGHT] != 1 ||
                        spatial[parameter_type_t::OUTPUT_WIDTH]  != 1;
                    const bool reduction_needs_gemm = (Kf > 1) && (full_P != 1 || full_Q != 1);
                    const bool dram_queue =
                        scheduler->input_offset_dram.size()  != 1 ||
                        scheduler->weight_offset_dram.size()  != 1 ||
                        scheduler->output_offset_dram.size() != 1;
                    const bool replayable =
                        (static_cast<size_t>(Bf)*Nf*Pf*Qf*Kf > 1) &&
                        !filter_fold && !pq_partly_spatial && !reduction_needs_gemm &&
                        !dram_queue && Nf > 0 && full_N % Nf == 0;
                    if(replayable) {
                        float *in_base  = layer->input_data;
                        float *wt_base  = layer->weight;
                        float *out_base = layer->output_data;
                        const unsigned spatial_N = full_N / Nf;                        // channels per array tile
                        const unsigned sK        = spatial[parameter_type_t::INPUT_CHANNEL]; // reduction on PE_Y
                        const size_t   in_W      = layer->input_width;
                        const size_t   cstride   = layer->stride ? layer->stride : 1;
                        std::vector<float> scratch(spatial_N, 0.0f);
                        for(unsigned b = 0; b < Bf; ++b) {
                        for(unsigned p = 0; p < Pf; ++p) {
                        for(unsigned q = 0; q < Qf; ++q) {
                        for(unsigned n = 0; n < Nf; ++n) {
                            float *out_tile = out_base
                                + static_cast<size_t>(b)*layer->output_size          // one output sample
                                + (static_cast<size_t>(p)*full_Q + q)*full_N          // position (p,q)
                                + static_cast<size_t>(n)*spatial_N;                   // N-block
                            for(unsigned kt = 0; kt < Kf; ++kt) {
                                if(b==0 && p==0 && q==0 && n==0 && kt==0) continue;   // representative did (0..)
                                layer->input_data  = in_base
                                    + static_cast<size_t>(b)*layer->input_size        // input sample
                                    + static_cast<size_t>(p)*cstride*in_W             // conv row window
                                    + static_cast<size_t>(q)*cstride                  // conv col window
                                    + static_cast<size_t>(kt)*sK;                     // K-tile slice
                                layer->weight = wt_base
                                    + static_cast<size_t>(kt)*full_N*sK               // reduction-tile-major
                                    + static_cast<size_t>(n)*spatial_N*sK;            // N-block
                                const bool first_tile = (kt == 0);
                                if(first_tile) {
                                    layer->output_data = out_tile;
                                    run_pass();
                                } else {
                                    std::fill(scratch.begin(), scratch.end(), 0.0f);
                                    layer->output_data = scratch.data();
                                    run_pass();
                                    for(unsigned c = 0; c < spatial_N; ++c) out_tile[c] += scratch[c];
                                }
                            }
                        }}}}
                        layer->input_data  = in_base;
                        layer->weight      = wt_base;
                        layer->output_data = out_base;
                    }
                }
                // G4 raw-stage checkpoint (plan §5.3): the datapath + replay just produced the
                // RAW reduction accumulators and the finalize has not run yet -- compare them
                // here when the fixture supplies functional_golden_raw<i>.
                if(functional_external_golden && functional_layer_golden_raw.count(index))
                    verify_buffer_against(index, functional_layer_golden_raw[index], "raw");
                // FINALIZE (plan §5.3): the datapath produced the RAW accumulators; now apply
                // bias then activation EXACTLY ONCE per output element, after its reduction is
                // complete. bias is per-output-channel (layer->output_channel); the channel is
                // the innermost axis of the functional POSITION-MAJOR output ([B][P][Q][N] for
                // conv, [B][N] for connected), so channel = i % output_channel. Activation
                // formulas mirror Nebula's so the accelerator agrees with the external golden,
                // and the finalized values feed the next layer (activation-correct chaining).
                // (Executable-IR mapped ops finalize here too: their transitional layer
                // carries the fused activation and the artifact-copied bias; BN/requant/
                // zero-point knobs are legacy-fixture settings and stay at their defaults.)
                if(mapped && layer->output_channel > 0) {
                    const size_t elems = static_cast<size_t>(layer->output_size)*network->batch_size;
                    const unsigned Nch = layer->output_channel;
                    // Channel index per element: the GEMM datapath writes POSITION-MAJOR
                    // ([B][P][Q][N] / [B][N]) so channel = i % N; the conv im2col kernel
                    // writes nebula's CHANNEL-MAJOR [B][N][P][Q] so channel = (i/(P*Q)) % N.
                    // For P=Q=1 the two coincide.
                    const size_t spatial = conv_value_kernel
                        ? static_cast<size_t>(layer->output_height)*layer->output_width : 1;
                    float *od = layer->output_data;
                    const float *bias = layer->get_bias();   // NULL if this layer has no bias
                    // Inference BatchNorm (fused, applied ONCE before bias/activation, Nebula
                    // semantics): v = scale*(raw - rolling_mean)/(sqrt(rolling_variance)+1e-5).
                    const bool   bn  = layer->has_batchnorm();
                    const float *bsc = bn ? layer->get_bn_scale()    : NULL;
                    const float *bmu = bn ? layer->get_bn_mean()     : NULL;
                    const float *bvar= bn ? layer->get_bn_variance() : NULL;
                    // Asymmetric int8 zero-point correction. Full expansion:
                    //   Sum((qi-zp_i)(qw-zp_w)) = Sum(qi*qw) - zp_i*colsum_w[n]
                    //                             - zp_w*rowsum_i[m] + Kdim*zp_i*zp_w.
                    // colsum_w[n] (per output channel) folds like bias; rowsum_i[m] (per output
                    // ROW) does NOT, so the simulator computes it from the injected input rows.
                    const int zp_i = functional_input_zero_point;
                    const int zp_w = functional_weight_zero_point;
                    const size_t Kdim = Nch ? static_cast<size_t>(layer->weight_size)/Nch : 0;   // reduction len
                    std::vector<double> colsum;   // per output channel n (weight column sum)
                    if(zp_i != 0 && layer->weight != NULL && layer->weight_size > 0) {
                        colsum.assign(Nch, 0.0);
                        for(unsigned n = 0; n < Nch; ++n) {
                            double s = 0.0;
                            for(size_t k = 0; k < Kdim; ++k) s += layer->weight[static_cast<size_t>(n)*Kdim + k];
                            colsum[n] = s;
                        }
                    }
                    std::vector<double> rowsum;   // per output row m (input row sum)
                    const size_t rows = Nch ? elems / Nch : 0;
                    if(zp_w != 0 && layer->input_data != NULL && layer->input_size > 0) {
                        rowsum.assign(rows, 0.0);
                        for(size_t m = 0; m < rows; ++m) {
                            double s = 0.0;
                            for(size_t k = 0; k < static_cast<size_t>(layer->input_size); ++k)
                                s += layer->input_data[m*static_cast<size_t>(layer->input_size) + k];
                            rowsum[m] = s;
                        }
                    }
                    const double zp_const = static_cast<double>(Kdim)*zp_i*zp_w;
                    for(size_t i = 0; i < elems; ++i) {
                        const unsigned c = spatial > 1 ? (i/spatial) % Nch : i % Nch;
                        const size_t   m = Nch ? i / Nch : 0;
                        float v = od[i];
                        if(bn) v = bsc[c]*(v - bmu[c])/(std::sqrt(bvar[c]) + 0.00001f);
                        if(!colsum.empty()) v += static_cast<float>(-static_cast<double>(zp_i)*colsum[c]);
                        if(!rowsum.empty()) v += static_cast<float>(-static_cast<double>(zp_w)*rowsum[m]);
                        if(zp_i != 0 && zp_w != 0) v += static_cast<float>(zp_const);
                        v += (bias ? bias[c] : 0.0f);
                        switch(layer->activation_type) {
                            case nebula::RELU_ACTIVATION:  v = v > 0.0f ? v : 0.0f;        break;
                            case nebula::LEAKY_ACTIVATION: v = v > 0.0f ? v : 0.1f*v;      break;
                            case nebula::LINEAR_ACTIVATION: default:                       break;
                        }
                        if(functional_requant_shift > 0) {
                            // INT8 requantization: (acc*mult + round) >> shift, then clamp. Done
                            // in int64 (v is an exact integer; per-channel mult can push the
                            // product beyond float's 2^24). mult[c] per output channel, or 1.
                            long long a = llroundf(v);
                            const long long mult = functional_requant_mult.empty()
                                ? 1LL : static_cast<long long>(llroundf(functional_requant_mult[c]));
                            a = (a*mult + (1LL << (functional_requant_shift - 1))) >> functional_requant_shift;
                            if(a < functional_requant_min) a = functional_requant_min;
                            if(a > functional_requant_max) a = functional_requant_max;
                            v = static_cast<float>(a);
                        }
                        // Low-precision OUTPUT rounding (fp16/bf16); no-op when format is fp32.
                        od[i] = round_lowp(functional_output_format, v);
                    }
                }
                // G1: publish the finalized executable-op values into the tensor store
                // (the DAG's value medium) and compare against the artifact golden.
                if(executable_ir_mode) functional_commit_executable_operation(index, *operation);
#endif
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
#ifdef FUNCTIONAL
                functional_execute_graph_operation(index, *operation);
#endif
            } else if(network->layers[index]->layer_type == nebula::SOFTMAX_LAYER && !sfus.empty()) {
                run_standalone_softmax(index, m_accelerator_config, m_network_config);
            } else {
                ++num_skipped_timing_layers;
                std::cerr << "Warning: network layer " << index
                          << " is excluded from accelerator timing (only convolution/connected are supported)"
                          << std::endl;
            }
#ifdef FUNCTIONAL
            // Non-MAC functional kernels (plan §5, "elementwise add"): layers with no MAC
            // datapath still carry values in a functional run. Compute them with the layer's own
            // kernel, which READS the accelerator's prior-layer output_data (already finalized
            // with bias+activation) and writes this layer's output_data for the next layer.
            // Residual/shortcut add: output = activation(prev_layer_output + skip_source_output).
            // Pooling (max/avg): output = pool(prev_layer_output over the filter window).
            // MAC layers are never forwarded here -- that would overwrite the accelerator result.
            if(!mapped && !executable_ir_mode) {
                nebula::layer_t *fl = network->layers[index];
                const nebula::layer_type_t lt = fl->layer_type;
                if(lt == nebula::SHORTCUT_LAYER ||
                   lt == nebula::MAXPOOL_LAYER  ||
                   lt == nebula::AVGPOOL_LAYER  ||
                   lt == nebula::SOFTMAX_LAYER) {
                    fl->forward();                        // softmax: per-(batch,group) exp/normalize
                } else if(lt == nebula::CONCAT_LAYER) {
                    // Concat copies each source layer's output into this layer's buffer along the
                    // channel axis. Nebula's forward() ADVANCES output_data as it copies, leaving
                    // the member pointer past the end; save/restore so verify reads the buffer base.
                    float *saved = fl->output_data;
                    fl->forward();
                    fl->output_data = saved;
                }
            }
            verify_functional_layer(index, mapped && !executable_ir_mode);
#endif
        }
        if(executable_ir_mode && mapping_index != schedulers.size()) {
            std::cerr << "Error: executable DAG did not consume every mapping" << std::endl;
            exit(1);
        }
        print_total_result(m_accelerator_config, m_network_config);
#ifdef FUNCTIONAL
        std::cout << "[FUNCTIONAL] summary: " << functional_layers_checked
                  << " layer(s) verified, " << functional_layers_failed
                  << " failed"
                  << (functional_external_golden && functional_layers_checked == 0
                      ? "  (NO LAYER COMPARED -- gate fails)" : "")
                  << std::endl;
        if(functional_external_golden) write_functional_report(m_network_config);
#endif
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
        multi_chip->request_data(scheduler);
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
            pe_arrays[i]->request_data(scheduler);
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
    // Scratchpad absorb rate for the decoder's dense output: the GLB weight partition's
    // write port (line bits -> bytes, per its unit write cycle). This is the sink stage of
    // the supply -> decode -> scratchpad pipeline.
    double sink_bytes_per_cycle = 0.0;
    if(!global_buffers.empty()) {
        global_buffer_t *glb = global_buffers[0];
        const double write_cycle = glb->u_write_cycle[data_type_t::WEIGHT];
        if(write_cycle > 0.0) {
            sink_bytes_per_cycle =
                static_cast<double>(glb->line_size[data_type_t::WEIGHT])/8.0/write_cycle;
        }
    }
    layer_stats[m_stats_index]->apply_decompression(decomp, dense_weight_bytes,
                                                    dram_bytes_per_cycle, tile_bytes,
                                                    sink_bytes_per_cycle);
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
#ifdef FUNCTIONAL
// G5 (gaps plan Step 1): refuse mappings whose values the functional path cannot compute
// correctly. Each rejected class would otherwise surface as a confusing value-mismatch FAIL
// (or a silently-unreplayed fold): the DRAM-queue mechanism iterates tiles without
// re-initializing accumulators, chip-boundary reduction still moves psums with data_copy
// (accumulate exists only at the PE->PE_Y adder-tree boundary, scheduler.cc), filter folds
// have no replay ownership, and the native-conv on-chip offset network derives the
// per-channel input stride from the PE tile's extent instead of the full input whenever
// P/Q > 1. The supported envelope is documented in validation/functional/README.md.
void npu_t::functional_reject_unsupported_mapping(unsigned m_index) {
    auto reject = [&](const std::string &m_reason) {
        std::cerr << "Error: unsupported functional mapping for layer " << m_index
                  << ": " << m_reason
                  << " (see validation/functional/README.md for the supported envelope)"
                  << std::endl;
        exit(1);
    };
    // G3: a convolution with P/Q > 1 is computed by the in-simulator im2col value kernel
    // (functional_conv_im2col), which is mapping-independent -- the mapping then only
    // shapes timing, so the value-path mapping checks below do not apply. The kernel's
    // finalize supports bias/BN/activation/requant but not the GEMM-row zero-point
    // corrections, which assume a [M][K] input layout.
    if(scheduler->layer_name == layer_name_t::CONVOLUTIONAL_LAYER &&
       (network->layers[m_index]->output_height > 1 ||
        network->layers[m_index]->output_width  > 1)) {
        if(functional_input_zero_point != 0 || functional_weight_zero_point != 0) {
            reject("zero-point correction on a convolution computed by the im2col value "
                   "kernel; asymmetric quantization is GEMM-only for now");
        }
        return;
    }
    if(scheduler->input_offset_dram.size()  != 1 ||
       scheduler->weight_offset_dram.size() != 1 ||
       scheduler->output_offset_dram.size() != 1) {
        reject("DRAM-level temporal fold (offset queue size > 1); "
               "express the fold as GLB repetitions instead");
    }
    mapping_table_t *mt = scheduler->mapping_table;
    const std::vector<unsigned> full    = mt->calculate_total_parameter_size();
    const std::vector<unsigned> spatial = mt->calculate_parameter_size(component_type_t::DRAM);
    // G7: a CHIPS_Y reduction split is supported (GLB->multi-chip accumulate); CHIPS_X
    // has no accumulate convention (the chip index keys outputs by its X part).
    if(scheduler->chip_reduction_x) {
        reject("reduction dimension (C/R/S) split across CHIPS_X; "
               "map chip-level reduction onto CHIPS_Y");
    }
    auto fold = [&](parameter_type_t d) -> unsigned {
        return spatial[d] ? full[d]/spatial[d] : 1;
    };
    if(fold(parameter_type_t::FILTER_HEIGHT) != 1 ||
       fold(parameter_type_t::FILTER_WIDTH)  != 1) {
        reject("filter (R/S) temporal fold; no replay ownership for filter tiles");
    }
    // An INPUT_CHANNEL temporal fold is replayed assuming the reduction-tile-major
    // [Kf][N][sK] weight layout; a standard [N][K] fixture would be read with the wrong
    // strides and fail as a value mismatch. The fixture must declare the layout.
    if(fold(parameter_type_t::INPUT_CHANNEL) > 1 && functional_weight_layout != "ktile") {
        reject("INPUT_CHANNEL temporal fold with a standard [N][K] weight; regenerate the "
               "fixture with the reduction-tile-major layout and declare "
               "[data] weight_layout = ktile, or keep the reduction spatial (PE_Y)");
    }
    // The zero-point finalize corrections read layer->weight as [N][K]; the ktile layout
    // would produce wrong per-channel column sums.
    if((functional_input_zero_point != 0 || functional_weight_zero_point != 0) &&
       functional_weight_layout == "ktile") {
        reject("zero-point correction with the reduction-tile-major weight layout; "
               "use a spatial-reduction mapping for asymmetric quantization");
    }
}

// G1 (gaps plan Step 4): executable-IR functional execution over a tensor store keyed by
// executable tensor id. Aliases resolve to their storage tensor, so a consumer reading
// through an elided reshape sees the producer's values.
std::vector<float> &npu_t::functional_store(const std::string &m_tensor_id) {
    return functional_tensor_store[workload->storage_tensor(m_tensor_id).id];
}

namespace {
// Does this operation produce a graph output? (labels its golden comparison "final")
bool produces_graph_output(const workload_graph_t *m_graph,
                           const workload_operation_t &m_operation) {
    for(const std::string &out : m_operation.outputs) {
        for(const std::string &graph_out : m_graph->outputs) {
            if(out == graph_out || out == m_graph->storage_tensor(graph_out).id) return true;
        }
    }
    return false;
}
} // namespace

// Bind one mapped executable operation (linear/conv2d) to its transitional nebula layer
// BEFORE the datapath runs: point input_data at the store, copy the artifact weight into
// the layer's own weight buffer (re-laid reduction-tile-major when the mapping folds
// INPUT_CHANNEL, which also satisfies the G5 layout contract), zero-then-copy the bias,
// and allocate the operation's output in the store. Pointer ownership note: only
// input_data is re-pointed (it is a borrowed pointer by nebula's contract); weight, bias
// and output_data keep nebula's own allocations.
void npu_t::functional_bind_executable_operation(unsigned m_index,
                                                 const workload_operation_t &m_operation) {
    nebula::layer_t *current = network->layers[m_index];
    auto bind_fail = [&](const std::string &m_message) {
        std::cerr << "Error: executable operation " << m_operation.id << ": " << m_message
                  << std::endl;
        exit(1);
    };
    // Data input from the store (a graph input from the artifact, or a prior op's output).
    const workload_tensor_t &in_decl = workload->tensor(m_operation.inputs.front());
    std::vector<float> &in = functional_store(m_operation.inputs.front());
    if(in.size() != in_decl.elements())
        bind_fail("input tensor " + in_decl.id + " has not been produced yet");
    current->input_data = in.data();

    // Output store allocation; the nebula buffer must agree on the element count.
    const workload_tensor_t &out_decl = workload->tensor(m_operation.outputs.front());
    functional_store(m_operation.outputs.front()).assign(out_decl.elements(), 0.0f);
    const size_t layer_elements =
        static_cast<size_t>(current->output_size)*network->batch_size;
    if(layer_elements != out_decl.elements())
        bind_fail("transitional layer holds " + std::to_string(layer_elements) +
                  " output elements, executable declares " +
                  std::to_string(out_decl.elements()));
    std::memset(current->output_data, 0, layer_elements*sizeof(float));

    if(m_operation.inputs.size() < 2) bind_fail("missing weight tensor");
    std::vector<float> &wt = functional_store(m_operation.inputs[1]);
    if(current->weight == NULL || current->weight_size != wt.size())
        bind_fail("weight buffer holds " + std::to_string(current->weight_size) +
                  " elements, artifact supplies " + std::to_string(wt.size()));
    functional_weight_layout.clear();
    bool relaid = false;
    if(m_operation.kind == WORKLOAD_LINEAR) {
        // The GEMM replay consumes an INPUT_CHANNEL temporal fold in reduction-tile-major
        // [Kf][N][sK] order; the artifact weight is the model's [N][K]. Re-lay here when
        // the mapping declares such a fold, which also satisfies the G5 layout contract.
        mapping_table_t *mt = scheduler->mapping_table;
        const std::vector<unsigned> full    = mt->calculate_total_parameter_size();
        const std::vector<unsigned> spatial = mt->calculate_parameter_size(component_type_t::DRAM);
        const unsigned sK = spatial[parameter_type_t::INPUT_CHANNEL];
        const unsigned Kf = sK ? full[parameter_type_t::INPUT_CHANNEL]/sK : 1;
        const unsigned N  = current->output_channel;
        const size_t   K  = N ? wt.size()/N : 0;
        if(Kf > 1 && K == static_cast<size_t>(Kf)*sK) {
            for(unsigned kt = 0; kt < Kf; ++kt)
                for(unsigned n = 0; n < N; ++n)
                    for(unsigned kk = 0; kk < sK; ++kk)
                        current->weight[(static_cast<size_t>(kt)*N + n)*sK + kk] =
                            wt[static_cast<size_t>(n)*K + kt*sK + kk];
            functional_weight_layout = "ktile";
            relaid = true;
        }
    }
    if(!relaid) std::memcpy(current->weight, wt.data(), wt.size()*sizeof(float));

    float *bias = current->get_bias();
    if(bias == NULL) {
        if(m_operation.inputs.size() >= 3) bind_fail("layer has no bias storage");
    } else {
        std::memset(bias, 0, sizeof(float)*current->output_channel);
        if(m_operation.inputs.size() >= 3) {
            std::vector<float> &bv = functional_store(m_operation.inputs[2]);
            if(bv.size() != current->output_channel)
                bind_fail("bias tensor holds " + std::to_string(bv.size()) +
                          " elements, layer has " + std::to_string(current->output_channel) +
                          " channels");
            std::memcpy(bias, bv.data(), bv.size()*sizeof(float));
        }
    }
}

// After the datapath + replay/kernel + finalize: publish the finalized values into the
// tensor store (the DAG's value medium) and compare against the artifact golden.
void npu_t::functional_commit_executable_operation(unsigned m_index,
                                                   const workload_operation_t &m_operation) {
    nebula::layer_t *current = network->layers[m_index];
    std::vector<float> &out = functional_store(m_operation.outputs.front());
    std::memcpy(out.data(), current->output_data, out.size()*sizeof(float));
    const auto golden = functional_artifact.golden.find(m_operation.id);
    if(golden != functional_artifact.golden.end()) {
        verify_buffer_against(m_index, out.data(), out.size(), golden->second,
                              produces_graph_output(workload, m_operation) ? "final" : "op");
    }
}

// Non-MAC executable operations carry values through dedicated kernels over the tensor
// store: softmax, max/average pool, elementwise add/multiply, concat, inference
// BatchNorm. Branch inputs are fetched by tensor ID, so DAG fan-in reads the correct
// producer regardless of operation order.
void npu_t::functional_execute_graph_operation(unsigned m_index,
                                               const workload_operation_t &m_operation) {
    auto op_fail = [&](const std::string &m_message) {
        std::cerr << "Error: executable operation " << m_operation.id << ": " << m_message
                  << std::endl;
        exit(1);
    };
    auto fetch = [&](const std::string &m_id) -> std::vector<float>& {
        std::vector<float> &values = functional_store(m_id);
        if(values.size() != workload->tensor(m_id).elements())
            op_fail("input tensor " + m_id + " has not been produced yet");
        return values;
    };
    const workload_tensor_t &out_decl = workload->tensor(m_operation.outputs.front());
    std::vector<float> &out = functional_store(m_operation.outputs.front());
    out.assign(out_decl.elements(), 0.0f);
    const workload_geometry_t &g = m_operation.geometry;

    switch(m_operation.kind) {
        case WORKLOAD_SOFTMAX: {
            std::vector<float> &in = fetch(m_operation.inputs.front());
            const size_t rows = g.rows, len = g.row_length;
            if(rows*len != in.size()) op_fail("softmax geometry disagrees with input");
            for(size_t r = 0; r < rows; ++r) {
                const float *x = in.data() + r*len;
                float *y = out.data() + r*len;
                float peak = x[0];
                for(size_t i = 1; i < len; ++i) peak = std::max(peak, x[i]);
                float sum = 0.0f;
                for(size_t i = 0; i < len; ++i) { y[i] = std::exp(x[i] - peak); sum += y[i]; }
                for(size_t i = 0; i < len; ++i) y[i] /= sum;
            }
            break;
        }
        case WORKLOAD_POOL2D: {
            std::vector<float> &in = fetch(m_operation.inputs.front());
            const unsigned B = g.batch, C = g.input_channels;
            const unsigned H = g.input_height, W = g.input_width;
            const unsigned P = g.output_height, Q = g.output_width;
            const bool max_mode = g.mode == "max";
            for(unsigned b = 0; b < B; ++b)
            for(unsigned c = 0; c < C; ++c)
            for(unsigned p = 0; p < P; ++p)
            for(unsigned q = 0; q < Q; ++q) {
                float best = -std::numeric_limits<float>::infinity();
                double sum = 0.0;
                unsigned valid = 0;
                for(unsigned kh = 0; kh < g.kernel_height; ++kh) {
                    const long ih = static_cast<long>(p)*g.stride_height - g.padding_height +
                                    static_cast<long>(kh)*g.dilation_height;
                    if(ih < 0 || ih >= static_cast<long>(H)) continue;
                    for(unsigned kw = 0; kw < g.kernel_width; ++kw) {
                        const long iw = static_cast<long>(q)*g.stride_width - g.padding_width +
                                        static_cast<long>(kw)*g.dilation_width;
                        if(iw < 0 || iw >= static_cast<long>(W)) continue;
                        const float v = in[((static_cast<size_t>(b)*C + c)*H + ih)*W + iw];
                        best = std::max(best, v);
                        sum += v;
                        ++valid;
                    }
                }
                const unsigned window = g.kernel_height*g.kernel_width;
                const unsigned samples = (!max_mode && g.count_include_pad) ? window : valid;
                out[((static_cast<size_t>(b)*C + c)*P + p)*Q + q] = max_mode
                    ? best : (samples ? static_cast<float>(sum/samples) : 0.0f);
            }
            break;
        }
        case WORKLOAD_ELEMENTWISE: {
            if(m_operation.inputs.size() != 2) op_fail("elementwise needs two inputs");
            std::vector<float> &a = fetch(m_operation.inputs[0]);
            std::vector<float> &b = fetch(m_operation.inputs[1]);
            if(a.size() != out.size() || b.size() != out.size())
                op_fail("elementwise shapes disagree");
            const bool multiply = g.elementwise_operator == "multiply";
            for(size_t i = 0; i < out.size(); ++i)
                out[i] = multiply ? a[i]*b[i] : a[i] + b[i];
            break;
        }
        case WORKLOAD_CONCAT: {
            // Copy each input's [axis:] block per outer row, in declared input order.
            const std::vector<size_t> &out_shape = out_decl.shape;
            if(g.axis >= out_shape.size()) op_fail("concat axis out of range");
            size_t outer = 1;
            for(size_t d = 0; d < g.axis; ++d) outer *= out_shape[d];
            std::vector<size_t> inner(m_operation.inputs.size());
            size_t inner_total = 0;
            for(size_t i = 0; i < m_operation.inputs.size(); ++i) {
                const workload_tensor_t &decl = workload->tensor(m_operation.inputs[i]);
                inner[i] = 1;
                for(size_t d = g.axis; d < decl.shape.size(); ++d) inner[i] *= decl.shape[d];
                inner_total += inner[i];
            }
            if(outer*inner_total != out.size()) op_fail("concat geometry disagrees");
            for(size_t o = 0; o < outer; ++o) {
                size_t cursor = o*inner_total;
                for(size_t i = 0; i < m_operation.inputs.size(); ++i) {
                    std::vector<float> &in = fetch(m_operation.inputs[i]);
                    std::memcpy(out.data() + cursor, in.data() + o*inner[i],
                                inner[i]*sizeof(float));
                    cursor += inner[i];
                }
            }
            break;
        }
        case WORKLOAD_BATCH_NORM: {
            // Torch aten order: (input, weight, bias, running_mean, running_var).
            if(m_operation.inputs.size() != 5)
                op_fail("batch_norm functional execution requires the full affine form "
                        "(input, weight, bias, running_mean, running_var)");
            std::vector<float> &in    = fetch(m_operation.inputs[0]);
            std::vector<float> &scale = fetch(m_operation.inputs[1]);
            std::vector<float> &shift = fetch(m_operation.inputs[2]);
            std::vector<float> &mean  = fetch(m_operation.inputs[3]);
            std::vector<float> &var   = fetch(m_operation.inputs[4]);
            const unsigned C = g.output_channels;
            if(in.size() != out.size() || in.size() % C != 0)
                op_fail("batch_norm geometry disagrees");
            const size_t spatial = in.size()/workload->tensor(m_operation.inputs[0]).shape[0]/C;
            for(size_t i = 0; i < in.size(); ++i) {
                const unsigned c = (i/spatial) % C;
                out[i] = scale[c]*(in[i] - mean[c])/
                         std::sqrt(var[c] + static_cast<float>(g.epsilon)) + shift[c];
            }
            break;
        }
        default:
            op_fail("no functional kernel for this operation kind");
    }

    // Fused activation (exec IR: linear/relu/leaky).
    if(m_operation.activation == "relu") {
        for(size_t i = 0; i < out.size(); ++i) out[i] = out[i] > 0.0f ? out[i] : 0.0f;
    } else if(m_operation.activation == "leaky") {
        for(size_t i = 0; i < out.size(); ++i) out[i] = out[i] > 0.0f ? out[i] : 0.1f*out[i];
    }

    const auto golden = functional_artifact.golden.find(m_operation.id);
    if(golden != functional_artifact.golden.end()) {
        verify_buffer_against(m_index, out.data(), out.size(), golden->second,
                              produces_graph_output(workload, m_operation) ? "final" : "op");
    }
}

// G3 (gaps plan Step 5): in-simulator conv value kernel. Lowers the mapped convolution to
// the im2col GEMM a GEMM accelerator actually executes -- zero padding (no OOB source
// reads), stride, and grouped/depthwise -- and writes the RAW accumulators (no bias, no
// activation; the shared FINALIZE applies those once). Output is CHANNEL-MAJOR
// [B][N][P][Q], nebula's tensor layout, so the result chains into nebula's non-MAC
// kernels (pool/shortcut/...) and matches a channel-major golden. Accumulation order is
// (c,r,s), the same as the PE mac_operation loop, so fp32 rounding matches the scalar
// CPU reference. Timing is untouched: the mapped datapath already ran for this layer.
void npu_t::functional_conv_im2col(unsigned m_index) {
    nebula::layer_t *l = network->layers[m_index];
    const unsigned batch  = network->batch_size;
    const unsigned C = l->input_channel,  H = l->input_height,  W = l->input_width;
    const unsigned N = l->output_channel, P = l->output_height, Q = l->output_width;
    const unsigned R = l->filter_height,  S = l->filter_width;
    const unsigned stride = l->stride ? l->stride : 1;
    const unsigned groups = l->group ? l->group : 1;
    if(N % groups != 0 || C % groups != 0) {
        std::cerr << "Error: layer " << m_index << " groups=" << groups
                  << " does not divide channels C=" << C << " N=" << N << std::endl;
        exit(1);
    }
    const unsigned Cg = C/groups, Ng = N/groups;
    const float *in = l->input_data;
    const float *wt = l->weight;                     // [N][Cg][R][S]
    float *out = l->output_data;
    for(unsigned b = 0; b < batch; ++b) {
        for(unsigned n = 0; n < N; ++n) {
            const unsigned g = n/Ng;
            for(unsigned p = 0; p < P; ++p) {
                for(unsigned q = 0; q < Q; ++q) {
                    float acc = 0.0f;
                    for(unsigned c = 0; c < Cg; ++c) {
                        for(unsigned r = 0; r < R; ++r) {
                            const long ih = static_cast<long>(p)*stride + r - l->padding_h;
                            if(ih < 0 || ih >= static_cast<long>(H)) continue;   // zero pad
                            for(unsigned s = 0; s < S; ++s) {
                                const long iw = static_cast<long>(q)*stride + s - l->padding_w;
                                if(iw < 0 || iw >= static_cast<long>(W)) continue;
                                acc += in[((static_cast<size_t>(b)*C + g*Cg + c)*H + ih)*W + iw]*
                                       wt[((static_cast<size_t>(n)*Cg + c)*R + r)*S + s];
                            }
                        }
                    }
                    out[((static_cast<size_t>(b)*N + n)*P + p)*Q + q] = acc;
                }
            }
        }
    }
    functional_kernel_layers.insert(m_index);
}

// Functional verification of one layer. The accelerator datapath has already moved the
// tensors through DRAM -> multi-chip -> GLB -> PE array -> MAC and written its computed
// output back into layer->output_data (raw MAC accumulation -- no bias, no activation).
// Nebula's forward() is the single functional owner of bias + activation (see
// components/sfu.h policy), so the comparison applies the layer's activation to the
// accelerator snapshot and adds nothing else: any residual difference is a datapath defect.
// The reference-nonzero count is printed so a zero-weight fixture cannot masquerade as a
// meaningful PASS.
//
// KNOWN LIMITATION (surfaced by this check): the functional data path MOVES output tiles
// with data_copy, so it does not SUM partial sums across a reduction dimension -- a
// reduction mapped spatially (PE_Y) or temporally (a GLB/DRAM fold) yields only a partial
// output, so this check reports FAIL for any mapping with output reduction. Making it PASS
// needs load-accumulate-store semantics on the output write-back, which is a separate piece
// of work from the output write-back CHAIN (offsets + per-level flush) fixed here.
// Compare the mapped layer's accelerator output against the external CPU golden. The
// accelerator writes the RAW reduction result (no bias, no activation) into output_data, so
// the fixture golden is the raw A x W^T (bias 0, linear) for this single-op milestone;
// bias/activation stages join once Phase 5 wires post-op ownership. Any element outside
// tolerance, or an all-zero golden, is a failure.
void npu_t::verify_against_golden(unsigned m_index) {
    verify_buffer_against(m_index, functional_golden, "final");
}

void npu_t::verify_buffer_against(unsigned m_index, const std::vector<float> &m_golden,
                                  const char *m_stage) {
    nebula::layer_t *current = network->layers[m_index];
    const size_t elements = static_cast<size_t>(current->output_size)*network->batch_size;
    verify_buffer_against(m_index, current->output_data, elements, m_golden, m_stage);
}

void npu_t::verify_buffer_against(unsigned m_index, const float *m_actual, size_t m_elements,
                                  const std::vector<float> &m_golden, const char *m_stage) {
    nebula::layer_t *current = network->layers[m_index];
    const size_t elements = m_elements;
    ++functional_layers_checked;
    if(m_golden.size() != elements) {
        std::cout << "[FUNCTIONAL] layer " << m_index << " (stage " << m_stage
                  << "): FAIL -- golden has " << m_golden.size()
                  << " floats, output has " << elements << std::endl;
        ++functional_layers_failed;
        return;
    }
    if(getenv("FVERIFY")) {
        for(size_t i = 0; i < elements; ++i)
            fprintf(stderr, "[V] L%u %s i=%zu golden=%.6f accel=%.6f\n",
                    m_index, m_stage, i, (double)m_golden[i], (double)m_actual[i]);
    }
    size_t mismatches = 0, golden_nonzeros = 0, accel_nonzeros = 0;
    double max_abs_diff = 0.0;
    int first_bad = -1;
    for(size_t i = 0; i < elements; ++i) {
        const double ref = m_golden[i];
        const double got = m_actual[i];
        if(ref != 0.0) ++golden_nonzeros;
        if(got != 0.0) ++accel_nonzeros;
        const double diff = std::fabs(got - ref);
        if(diff > max_abs_diff) max_abs_diff = diff;
        if(diff > 1.0e-4 + 1.0e-4*std::fabs(ref)) {
            ++mismatches;
            if(first_bad < 0) first_bad = static_cast<int>(i);
        }
    }
    const bool pass = mismatches == 0 && golden_nonzeros > 0;
    if(!pass) ++functional_layers_failed;
    std::cout << std::setprecision(6) << std::defaultfloat
              << "[FUNCTIONAL] layer " << m_index << " (external golden, stage "
              << m_stage << "): "
              << (pass ? "PASS" : "FAIL") << " -- " << elements << " elements, "
              << mismatches << " mismatch(es), max |diff| " << max_abs_diff
              << ", nonzeros golden/accel " << golden_nonzeros << "/" << accel_nonzeros
              << (golden_nonzeros == 0 ? "  (VACUOUS: all-zero golden)" : "");
    if(first_bad >= 0) {
        std::cout << "; first mismatch @" << first_bad << " golden="
                  << m_golden[first_bad] << " accel="
                  << m_actual[first_bad];
    }
    std::cout << std::endl;

    // G6 (pruning masked-dense): input/weight zero statistics, so a masked run documents
    // the sparsity it actually consumed (value-impact reporting only -- no speedup claim).
    // In executable-IR mode only MAPPED operations have live layer operand pointers (the
    // store-bound input and the artifact-copied weight); a non-MAC operation's transitional
    // layer keeps stale placeholder pointers, which must not be scanned.
    const bool layer_operands_valid = !executable_ir_mode ||
        (workload != NULL && m_index < workload->operations.size() &&
         workload->operations[m_index].mapping_required);
    size_t weight_zeros = 0, weight_elems = 0, input_zeros = 0, input_elems = 0;
    if(layer_operands_valid && current->weight != NULL && current->weight_size > 0) {
        weight_elems = current->weight_size;
        for(size_t i = 0; i < weight_elems; ++i)
            if(current->weight[i] == 0.0f) ++weight_zeros;
    }
    if(layer_operands_valid && current->input_data != NULL && current->input_size > 0) {
        input_elems = static_cast<size_t>(current->input_size)*network->batch_size;
        for(size_t i = 0; i < input_elems; ++i)
            if(current->input_data[i] == 0.0f) ++input_zeros;
    }

    // Machine-readable record for the per-op report (plan §5.2).
    std::ostringstream row;
    row << std::setprecision(9) << std::defaultfloat
        << "{\"layer\":" << m_index
        << ",\"stage\":\"" << m_stage << "\""
        << (functional_kernel_layers.count(m_index) ? ",\"kernel\":\"im2col\"" : "")
        << ",\"pass\":" << (pass ? "true" : "false")
        << ",\"elements\":" << elements
        << ",\"mismatches\":" << mismatches
        << ",\"max_abs_diff\":" << max_abs_diff
        << ",\"golden_nonzeros\":" << golden_nonzeros
        << ",\"accel_nonzeros\":" << accel_nonzeros
        << ",\"weight_zeros\":" << weight_zeros
        << ",\"weight_elements\":" << weight_elems
        << ",\"input_zeros\":" << input_zeros
        << ",\"input_elements\":" << input_elems
        << ",\"vacuous\":" << (golden_nonzeros == 0 ? "true" : "false")
        << ",\"first_mismatch\":" << first_bad << "}";
    functional_report.push_back(row.str());
}

void npu_t::write_functional_report(const std::string &m_network_label) const {
    const std::string dir = "result/functional/" + m_network_label;
    // result/ and result/functional/ are created by the build; make the leaf best-effort.
    std::string cmd = "mkdir -p '" + dir + "'";
    if(system(cmd.c_str()) != 0) { /* fall through; ofstream failure is reported below */ }
    const std::string path = dir + "/report.json";
    std::ofstream out(path.c_str());
    if(!out.good()) {
        std::cerr << "Warning: could not write functional report to " << path << std::endl;
        return;
    }
    out << "{\"network\":\"" << m_network_label << "\""
        << ",\"semantics\":\"" << functional_semantics << "\""
        << ",\"layers_checked\":" << functional_layers_checked
        << ",\"layers_failed\":" << functional_layers_failed
        << ",\"gate_pass\":" << (functional_failed() ? "false" : "true")
        << ",\"per_layer\":[";
    for(size_t i = 0; i < functional_report.size(); ++i)
        out << (i ? "," : "") << functional_report[i];
    out << "]}" << std::endl;
    std::cout << "[FUNCTIONAL] report -> " << path << std::endl;
}

void npu_t::verify_functional_layer(unsigned m_index, bool m_mapped) {
    // Executable-IR verification happens per operation against the artifact goldens
    // (functional_commit/execute); the legacy per-layer machinery does not apply.
    if(executable_ir_mode) return;
    nebula::layer_t *current = network->layers[m_index];
    // Fixture path: compare the accelerator's raw output against an EXTERNAL golden and do
    // NOT run nebula forward() (no shared-oracle, no output_data overwrite). Only the mapped
    // op carries a golden in this single-op milestone.
    if(functional_external_golden) {
        // The external golden is the network's FINAL output. Compare only the last
        // VALUE-BEARING layer (conv/connected MAC layers plus non-MAC functional kernels such
        // as the residual/shortcut add). Earlier layers were already executed and their
        // output_data feeds the next layer (chaining), so they need no per-layer golden.
        if(functional_last_mapped < 0) {
            for(unsigned i = 0; i < network->layers.size(); ++i) {
                const nebula::layer_type_t t = network->layers[i]->layer_type;
                if(t == nebula::CONNECTED_LAYER || t == nebula::CONVOLUTIONAL_LAYER ||
                   t == nebula::SHORTCUT_LAYER  || t == nebula::MAXPOOL_LAYER ||
                   t == nebula::AVGPOOL_LAYER   || t == nebula::SOFTMAX_LAYER ||
                   t == nebula::CONCAT_LAYER)
                    functional_last_mapped = static_cast<int>(i);
            }
        }
        // G4: a per-layer golden (functional_golden<i>) verifies this layer's finalized
        // output wherever it sits in the DAG, localizing the first mismatching operation.
        if(functional_layer_golden.count(m_index))
            verify_buffer_against(m_index, functional_layer_golden[m_index], "layer");
        if(static_cast<int>(m_index) == functional_last_mapped)
            verify_against_golden(m_index);
        return;
    }
    if(!m_mapped) {
        current->forward();
        return;
    }
    const size_t elements = static_cast<size_t>(current->output_size)*network->batch_size;
    std::vector<float> accelerator(current->output_data, current->output_data + elements);
    current->forward();

    bool activation_comparable = true;
    switch(current->activation_type) {
        case nebula::LINEAR_ACTIVATION:
            break;
        case nebula::RELU_ACTIVATION:
            for(size_t i = 0; i < elements; ++i) {
                if(accelerator[i] < 0.0f) accelerator[i] = 0.0f;
            }
            break;
        default:
            activation_comparable = false;
            break;
    }

    size_t mismatches = 0, reference_nonzeros = 0, accelerator_nonzeros = 0;
    double max_abs_diff = 0.0;
    for(size_t i = 0; i < elements; ++i) {
        const double reference = current->output_data[i];
        if(reference != 0.0) ++reference_nonzeros;
        if(accelerator[i] != 0.0f) ++accelerator_nonzeros;
        const double diff = std::fabs(static_cast<double>(accelerator[i]) - reference);
        if(diff > max_abs_diff) max_abs_diff = diff;
        if(diff > 1.0e-4 + 1.0e-4*std::fabs(reference)) ++mismatches;
    }

    ++functional_layers_checked;
    const bool pass = activation_comparable && mismatches == 0;
    if(!pass) ++functional_layers_failed;
    std::cout << "[FUNCTIONAL] layer " << m_index << ": "
              << (activation_comparable
                  ? (pass ? "PASS" : "FAIL")
                  : "NOT COMPARED (unsupported activation for the functional check)")
              << " -- " << elements << " elements, " << mismatches << " mismatch(es), "
              << "max |diff| " << max_abs_diff << ", nonzeros ref/accel "
              << reference_nonzeros << "/" << accelerator_nonzeros
              << (reference_nonzeros == 0 ? "  (VACUOUS: all-zero reference)" : "")
              << std::endl;
}
#endif

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
