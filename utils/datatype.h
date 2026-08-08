#ifndef __DATATYPE_H__
#define __DATATYPE_H__

#include <cstddef>
#include <string>
#include <vector>

#include "def.h"

class section_config_t;

enum class data_format_kind_t {
    UINT,
    INT,
    FP32,
    FP16,
    BF16,
    MXFP
};

struct tensor_format_t {
    data_format_kind_t kind;
    unsigned payload_bits;
    bool is_signed;
    unsigned block_elements;
    unsigned scale_bits;
    std::string name;

    bool is_block_scaled() const { return kind == data_format_kind_t::MXFP; }
};

class runtime_datatypes_t {
public:
    runtime_datatypes_t();

    void configure(section_config_t &section);
    const tensor_format_t &format(data_type_t type) const;
    const tensor_format_t &accumulator_format() const;
    size_t payload_bits(data_type_t type, size_t elements) const;
    size_t metadata_bits(data_type_t type, size_t elements) const;
    size_t storage_bits(data_type_t type, size_t elements) const;
    size_t storage_bytes(data_type_t type, size_t elements) const;
    size_t payload_transactions(data_type_t type, size_t elements, size_t transaction_bits) const;
    size_t metadata_transactions(data_type_t type, size_t elements, size_t transaction_bits) const;
    size_t storage_transactions(data_type_t type, size_t elements, size_t transaction_bits) const;
    std::string describe(data_type_t type) const;

private:
    std::vector<tensor_format_t> formats;
    tensor_format_t accumulator;
};

runtime_datatypes_t &runtime_datatypes();

tensor_format_t parse_data_format(const std::string &name);

#endif
