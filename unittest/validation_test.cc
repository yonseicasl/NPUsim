#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "datatype.h"
#include "mapping_table.h"
#include "pe_lane.h"
#include "interconnect_timing.h"

namespace {

void fail(const std::string &message) {
    std::cerr << "validation test failed: " << message << std::endl;
    std::exit(1);
}

void validate_parser_contract() {
    section_config_t section("test");
    if(!section.add_setting("vector", "1:2:3") ||
       !section.add_setting("broadcast", "7") ||
       !section.add_setting("negative", "-1") ||
       section.add_setting("vector", "4:5:6")) {
        fail("duplicate-key contract");
    }

    std::vector<unsigned> values(3, 0);
    if(!section.get_vector_setting("vector", &values) ||
       values[0] != 1 || values[1] != 2 || values[2] != 3) {
        fail("exact vector parsing");
    }
    if(!section.get_vector_setting("broadcast", &values) ||
       values[0] != 7 || values[1] != 7 || values[2] != 7) {
        fail("scalar vector broadcast");
    }

    unsigned unsigned_value = 0;
    if(section.get_setting("negative", &unsigned_value)) {
        fail("negative unsigned value was accepted");
    }

    section_config_t short_vector("test");
    short_vector.add_setting("value", "1:2");
    values.assign(3, 0);
    if(short_vector.get_vector_setting("value", &values)) {
        fail("short vector was accepted");
    }

    section_config_t long_vector("test");
    long_vector.add_setting("value", "1:2:3:4");
    if(long_vector.get_vector_setting("value", &values)) {
        fail("long vector was accepted");
    }
}

void validate_legacy_glb_mapping_contract() {
    section_config_t section("mapping");
    section.add_setting("glb", "2,3,1,1,1,1,1,0,0,2,1");
    mapping_table_t table(section);
    if(table.calculate_active_component(component_type_t::GLOBAL_BUFFER) != 3) {
        fail("legacy GLB temporal repetition mapping");
    }
}

void validate_runtime_datatypes() {
    struct dense_format_case_t {
        const char *name;
        size_t bits;
    };
    const dense_format_case_t dense_formats[] = {
        {"fp32", 32}, {"fp16", 16}, {"bf16", 16}, {"int8", 8},
        {"int4", 4}, {"int2", 2}, {"uint8", 8}
    };
    for(const dense_format_case_t &format : dense_formats) {
        section_config_t format_section("accelerator");
        format_section.add_setting("data_format", format.name);
        runtime_datatypes().configure(format_section);
        if(runtime_datatypes().describe(data_type_t::INPUT) != format.name ||
           runtime_datatypes().describe(data_type_t::WEIGHT) != format.name ||
           runtime_datatypes().describe(data_type_t::OUTPUT) != format.name ||
           runtime_datatypes().payload_bits(data_type_t::INPUT, 3) != 3*format.bits ||
           runtime_datatypes().metadata_bits(data_type_t::INPUT, 3) != 0 ||
           runtime_datatypes().storage_transactions(data_type_t::INPUT, 3, 8) !=
               (3*format.bits + 7)/8) {
            fail(std::string("dense datatype format contract: ") + format.name);
        }
    }

    section_config_t section("accelerator");
    section.add_setting("data_format", "bf16");
    section.add_setting("weight_format", "mxfp4");
    section.add_setting("output_format", "int4");
    section.add_setting("accumulator_format", "fp32");
    runtime_datatypes().configure(section);

    if(runtime_datatypes().describe(data_type_t::INPUT) != "bf16" ||
       runtime_datatypes().storage_bytes(data_type_t::INPUT, 32) != 64 ||
       runtime_datatypes().describe(data_type_t::WEIGHT) != "mxfp4_b32_e8m0" ||
       runtime_datatypes().payload_bits(data_type_t::WEIGHT, 32) != 128 ||
       runtime_datatypes().metadata_bits(data_type_t::WEIGHT, 32) != 8 ||
       runtime_datatypes().storage_bytes(data_type_t::WEIGHT, 32) != 17 ||
       runtime_datatypes().storage_bytes(data_type_t::WEIGHT, 3) != 3 ||
       runtime_datatypes().payload_transactions(data_type_t::WEIGHT, 33, 16) != 9 ||
       runtime_datatypes().metadata_transactions(data_type_t::WEIGHT, 33, 16) != 1 ||
       runtime_datatypes().storage_transactions(data_type_t::WEIGHT, 33, 32) != 6 ||
       runtime_datatypes().storage_transactions(data_type_t::WEIGHT, 33, 16) != 10 ||
       runtime_datatypes().accumulator_format().name != "fp32" ||
       runtime_datatypes().describe(data_type_t::OUTPUT) != "int4" ||
       runtime_datatypes().storage_bytes(data_type_t::OUTPUT, 3) != 2 ||
       runtime_datatypes().storage_transactions(data_type_t::OUTPUT, 3, 8) != 2) {
        fail("runtime datatype storage accounting");
    }

    section_config_t interleaved("accelerator");
    interleaved.add_setting("weight_format", "mxfp4");
    interleaved.add_setting("mxfp_metadata_layout", "interleaved");
    runtime_datatypes().configure(interleaved);
    if(runtime_datatypes().metadata_layout() != mxfp_metadata_layout_t::INTERLEAVED ||
       runtime_datatypes().storage_bytes(data_type_t::WEIGHT, 3) != 3 ||
       runtime_datatypes().storage_transactions(data_type_t::WEIGHT, 33, 32) != 5) {
        fail("interleaved MXFP metadata layout accounting");
    }

    runtime_datatypes().configure(section);
    const datatype_transfer_timing_t mxfp_transfer = datatype_transfer_timing(
        data_type_t::WEIGHT, 33, 16, 32, 64);
    const datatype_transfer_timing_t int4_transfer = datatype_transfer_timing(
        data_type_t::OUTPUT, 3, 8, 16, 32);
    if(mxfp_transfer.source_accesses != 10 ||
       mxfp_transfer.destination_accesses != 6 ||
       mxfp_transfer.link_transactions != 4 ||
       mxfp_transfer.pipeline_transactions != 10 ||
       mxfp_transfer.payload_link_transactions != 3 ||
       mxfp_transfer.metadata_link_transactions != 1 ||
       int4_transfer.source_accesses != 2 ||
       int4_transfer.destination_accesses != 1 ||
       int4_transfer.link_transactions != 1 ||
       int4_transfer.pipeline_transactions != 2 ||
       int4_transfer.payload_link_transactions != 1 ||
       int4_transfer.metadata_link_transactions != 0) {
        fail("global buffer bit-width transaction contract");
    }

}
void validate_hierarchy_datatype_transactions() {
    section_config_t section("accelerator");
    section.add_setting("input_format", "int2");
    section.add_setting("weight_format", "mxfp8");
    section.add_setting("output_format", "fp32");
    section.add_setting("mxfp_metadata_layout", "separate");
    runtime_datatypes().configure(section);

    const datatype_transfer_timing_t int2 = datatype_transfer_timing(
        data_type_t::INPUT, 17, 8, 32, 16);
    const datatype_transfer_timing_t mxfp8 = datatype_transfer_timing(
        data_type_t::WEIGHT, 33, 64, 32, 16);
    const datatype_transfer_timing_t fp32 = datatype_transfer_timing(
        data_type_t::OUTPUT, 3, 64, 128, 32);

    if(int2.source_accesses != 5 || int2.destination_accesses != 2 ||
       int2.link_transactions != 3 || int2.pipeline_transactions != 5 ||
       int2.payload_link_transactions != 3 || int2.metadata_link_transactions != 0 ||
       mxfp8.source_accesses != 6 || mxfp8.destination_accesses != 10 ||
       mxfp8.link_transactions != 18 || mxfp8.pipeline_transactions != 18 ||
       mxfp8.payload_link_transactions != 17 || mxfp8.metadata_link_transactions != 1 ||
       fp32.source_accesses != 2 || fp32.destination_accesses != 1 ||
       fp32.link_transactions != 3 || fp32.pipeline_transactions != 3 ||
       !is_valid_memory_line_bits(8) || !is_valid_memory_line_bits(256) ||
       is_valid_memory_line_bits(0) || is_valid_memory_line_bits(12)) {
        fail("hierarchy-wide datatype transaction contract");
    }
}

void validate_pe_lane_contract() {
    const mac_lane_state_t full = calculate_mac_lane_state(8, 8, 64);
    if(!full.valid || full.scalar_capacity != 64 ||
       full.active_accumulator_units != 8 || full.final_accumulator_lanes != 8 ||
       full.utilization != 1.0) {
        fail("full PE lane contract");
    }

    const mac_lane_state_t partial = calculate_mac_lane_state(8, 8, 20);
    if(!partial.valid || partial.scalar_capacity != 64 ||
       partial.active_accumulator_units != 3 || partial.final_accumulator_lanes != 4 ||
       partial.utilization != 0.3125) {
        fail("partial PE lane contract");
    }

    if(calculate_mac_lane_state(8, 8, 0).valid ||
       calculate_mac_lane_state(8, 8, 65).valid ||
       calculate_mac_lane_state(0, 8, 1).valid ||
       accumulate_issue_cycles(3, 0.5) != 1.5 ||
       accumulate_issue_cycles(0, 0.5) != 0.0 ||
       calculate_time_based_mac_utilization(84.0, 168.0) != 0.5 ||
       calculate_time_based_mac_utilization(0.0, 168.0) != 0.0 ||
       calculate_time_based_mac_utilization(1.0, 0.0) != 0.0) {
        fail("PE lane capacity or fractional issue cycle");
    }
}

void validate_spatial_interconnect_contract() {
    if(!is_supported_spatial_noc(noc_type_t::BUS) ||
       !is_supported_spatial_noc(noc_type_t::STORE_AND_FORWARD) ||
       !is_supported_spatial_noc(noc_type_t::MESH) ||
       !is_supported_multi_chip_nop(noc_type_t::BUS) ||
       spatial_noc_cost(noc_type_t::MESH, 2, 3).latency_multiplier != 3.0 ||
       spatial_noc_cost(noc_type_t::MESH, 2, 3).energy_multiplier != 1.5 ||
       spatial_noc_cost(noc_type_t::CROSSBAR, 2, 3).latency_multiplier != 1.0 ||
       is_supported_multi_chip_nop(noc_type_t::MESH) ||
       !is_supported_spatial_noc(noc_type_t::CROSSBAR) ||
       !has_valid_active_shape(4, 8, 4, 8) ||
       has_valid_active_shape(4, 8, 5, 8) ||
       has_valid_active_shape(4, 8, 4, 9) ||
       // (transactions, source, link, dest). One transaction has no overlap: source+link+dest
       // (=8, not the pre-fix 5+max(5,1)+2=12). N>=2 pipelines the stages.
       pipelined_transfer_cycles(0, 5.0, 1.0, 2.0) != 0.0 ||
       pipelined_transfer_cycles(1, 5.0, 1.0, 2.0) != 8.0 ||
       pipelined_transfer_cycles(4, 5.0, 1.0, 2.0) != 24.0 ||
       static_energy_for_cycles(2.5, 4.0) != 10.0) {
        fail("spatial interconnect timing contract");
    }
}

void validate_accelerator(config_t &config, const std::string &path) {
    unsigned accelerator_count = 0;
    unsigned pe_array_count = 0;
    unsigned global_buffer_count = 0;
    unsigned multi_chip_count = 0;
    unsigned dram_count = 0;
    unsigned num_chips = 0;
    unsigned height = 0;
    unsigned width = 0;

    for(section_config_t &section : config.sections) {
        if(section.name == "accelerator") {
            accelerator_count++;
            if(!section.get_setting("num_chips", &num_chips) || num_chips == 0) {
                fail(path + ": invalid num_chips");
            }
        } else if(section.name == "spatial_arch" ||
                  section.name == "adder_tree" ||
                  section.name == "systolic_array") {
            pe_array_count++;
        } else if(section.name == "shared" || section.name == "separate") {
            global_buffer_count++;
        } else if(section.name == "multi_chip") {
            multi_chip_count++;
            if(!section.get_setting("height", &height) ||
               !section.get_setting("width", &width) ||
               height == 0 || width == 0) {
                fail(path + ": invalid multi-chip dimensions");
            }
        } else if(section.name == "dram") {
            dram_count++;
        } else {
            fail(path + ": unknown accelerator section " + section.name);
        }
    }

    if(accelerator_count != 1 || pe_array_count != 1 ||
       global_buffer_count != 1 || multi_chip_count != 1 || dram_count != 1) {
        fail(path + ": component cardinality");
    }
    if(height > 0 && width > 0 && height * width != num_chips) {
        fail(path + ": num_chips does not match height*width");
    }
}
void validate_mapping(config_t &config, const std::string &path) {
    if(config.sections.empty()) {
        fail(path + ": empty mapping");
    }
    for(section_config_t &section : config.sections) {
        mapping_table_t table(section);
        (void)table.calculate_parameter_size(component_type_t::DRAM);
    }
}

bool contains(const std::string &value, const std::string &needle) {
    return value.find(needle) != std::string::npos;
}

} // namespace

int main(int argc, char **argv) {
    if(argc < 2) {
        std::cerr << "Usage: " << argv[0] << " CONFIG..." << std::endl;
        return 2;
    }

    validate_parser_contract();
    validate_runtime_datatypes();
    validate_hierarchy_datatype_transactions();
    validate_legacy_glb_mapping_contract();

    validate_pe_lane_contract();
    validate_spatial_interconnect_contract();
    for(int i = 1; i < argc; i++) {
        const std::string path = argv[i];
        config_t config;
        config.parse(path);
        if(contains(path, "/accelerators/")) {
            validate_accelerator(config, path);
        } else if(contains(path, ".map")) {
            validate_mapping(config, path);
        }
    }

    std::cout << "validated " << argc - 1 << " configuration files" << std::endl;
    return 0;
}
