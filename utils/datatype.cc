#include "datatype.h"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <limits>

#include "config.h"

namespace {

std::string normalized(std::string value) {
    value.erase(std::remove_if(value.begin(), value.end(),
                               [](unsigned char c) { return c == '-' || c == '_' || std::isspace(c); }),
                value.end());
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return value;
}

tensor_format_t scalar(data_format_kind_t kind, unsigned bits, bool is_signed, const char *name) {
    tensor_format_t format = {kind, bits, is_signed, 0, 0, name};
    return format;
}

size_t ceil_div(size_t numerator, size_t denominator) {
    return numerator / denominator + (numerator % denominator != 0 ? 1 : 0);
}

const char *config_key(data_type_t type) {
    switch(type) {
    case data_type_t::INPUT: return "input_format";
    case data_type_t::WEIGHT: return "weight_format";
    case data_type_t::OUTPUT: return "output_format";
    default: return "";
    }
}

} // namespace

runtime_datatypes_t::runtime_datatypes_t() : formats(data_type_t::NUM_DATA_TYPES,
                                                       scalar(data_format_kind_t::UINT, 8, false, "uint8")),
                                             accumulator(scalar(data_format_kind_t::FP32, 32, true, "fp32")) {
}

tensor_format_t parse_data_format(const std::string &value) {
    const std::string name = normalized(value);
    if(name == "fp32" || name == "float" || name == "float32") return scalar(data_format_kind_t::FP32, 32, true, "fp32");
    if(name == "fp16" || name == "float16" || name == "half") return scalar(data_format_kind_t::FP16, 16, true, "fp16");
    if(name == "bf16" || name == "bfloat16") return scalar(data_format_kind_t::BF16, 16, true, "bf16");
    if(name == "int8" || name == "i8") return scalar(data_format_kind_t::INT, 8, true, "int8");
    if(name == "int4" || name == "i4") return scalar(data_format_kind_t::INT, 4, true, "int4");
    if(name == "int2" || name == "i2") return scalar(data_format_kind_t::INT, 2, true, "int2");
    if(name == "uint8" || name == "u8") return scalar(data_format_kind_t::UINT, 8, false, "uint8");
    if(name == "mxfp8" || name == "mxfp8e4m3" || name == "mxfp8e5m2") {
        tensor_format_t format = {data_format_kind_t::MXFP, 8, true, 32, 8, "mxfp8_b32_e8m0"};
        return format;
    }
    if(name == "mxfp4" || name == "mxfp4e2m1") {
        tensor_format_t format = {data_format_kind_t::MXFP, 4, true, 32, 8, "mxfp4_b32_e8m0"};
        return format;
    }
    std::cerr << "Error: unsupported runtime data format '" << value
              << "' (supported: fp32, fp16, bf16, int8, int4, int2, uint8, mxfp8, mxfp4)" << std::endl;
    exit(1);
}

void runtime_datatypes_t::configure(section_config_t &section) {
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        const data_type_t type = static_cast<data_type_t>(i);
        std::string value;
        if(section.get_setting(config_key(type), &value)) formats[i] = parse_data_format(value);
    }

    std::string accumulator_value;
    if(section.get_setting("accumulator_format", &accumulator_value)) accumulator = parse_data_format(accumulator_value);

    std::string common;
    if(section.get_setting("data_format", &common)) {
        const tensor_format_t format = parse_data_format(common);
        for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
            if(!section.exists(config_key(static_cast<data_type_t>(i)))) formats[i] = format;
        }
    }
}

const tensor_format_t &runtime_datatypes_t::format(data_type_t type) const {
    if(static_cast<unsigned>(type) >= formats.size()) {
        std::cerr << "Error: invalid data type format lookup" << std::endl;
        exit(1);
    }
    return formats[static_cast<unsigned>(type)];
}

const tensor_format_t &runtime_datatypes_t::accumulator_format() const {
    return accumulator;
}

size_t runtime_datatypes_t::payload_bits(data_type_t type, size_t elements) const {
    const tensor_format_t &value = format(type);
    if(elements != 0 && value.payload_bits > std::numeric_limits<size_t>::max() / elements) {
        std::cerr << "Error: datatype payload size overflow" << std::endl;
        exit(1);
    }
    return elements * value.payload_bits;
}

size_t runtime_datatypes_t::metadata_bits(data_type_t type, size_t elements) const {
    const tensor_format_t &value = format(type);
    if(!value.is_block_scaled() || elements == 0) return 0;
    return ceil_div(elements, value.block_elements) * value.scale_bits;
}

size_t runtime_datatypes_t::storage_bits(data_type_t type, size_t elements) const {
    const size_t payload = payload_bits(type, elements);
    const size_t metadata = metadata_bits(type, elements);
    if(payload > std::numeric_limits<size_t>::max() - metadata) {
        std::cerr << "Error: datatype storage size overflow" << std::endl;
        exit(1);
    }
    return payload + metadata;
}

size_t runtime_datatypes_t::storage_bytes(data_type_t type, size_t elements) const {
    return ceil_div(storage_bits(type, elements), 8);
}

std::string runtime_datatypes_t::describe(data_type_t type) const {
    const tensor_format_t &value = format(type);
    return value.name;
}

runtime_datatypes_t &runtime_datatypes() {
    static runtime_datatypes_t instance;
    return instance;
}

size_t runtime_datatypes_t::payload_transactions(data_type_t type, size_t elements,
                                                  size_t transaction_bits) const {
    if(transaction_bits == 0) {
        std::cerr << "Error: datatype transaction width must be non-zero" << std::endl;
        exit(1);
    }
    return ceil_div(payload_bits(type, elements), transaction_bits);
}

size_t runtime_datatypes_t::metadata_transactions(data_type_t type, size_t elements,
                                                   size_t transaction_bits) const {
    if(transaction_bits == 0) {
        std::cerr << "Error: datatype transaction width must be non-zero" << std::endl;
        exit(1);
    }
    return ceil_div(metadata_bits(type, elements), transaction_bits);
}

size_t runtime_datatypes_t::storage_transactions(data_type_t type, size_t elements,
                                                  size_t transaction_bits) const {
    if(transaction_bits == 0) {
        std::cerr << "Error: datatype transaction width must be non-zero" << std::endl;
        exit(1);
    }
    return ceil_div(storage_bits(type, elements), transaction_bits);
}
