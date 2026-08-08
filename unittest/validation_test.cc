#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "datatype.h"
#include "mapping_table.h"

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

void validate_runtime_datatypes() {
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
       runtime_datatypes().accumulator_format().name != "fp32" ||
       runtime_datatypes().describe(data_type_t::OUTPUT) != "int4" ||
       runtime_datatypes().storage_bytes(data_type_t::OUTPUT, 3) != 2) {
        fail("runtime datatype storage accounting");
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
