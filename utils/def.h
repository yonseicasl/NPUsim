#ifndef __DEF_H__
#define __DEF_H__

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

enum data_format_t {
	CONVOLUTION, 
	GEMM,
	NUM_DATA_FORMATS,
};

static std::vector<std::string> data_format_str __attribute((unused)) = {
    "convolution",
    "gemm",
    "num_data_formats",
};

enum data_type_t {
    INPUT, //input, 
    WEIGHT, //weight, 
    OUTPUT, //output, 
    NUM_DATA_TYPES, //num_data_types,
};

// The parameter types of neural network.
// It includes the parameters of input data, weight, and output data.
enum parameter_type_t {
	OUTPUT_CHANNEL,         // output_channel,
	BATCH_SIZE,             // batch_size,
	OUTPUT_HEIGHT,          // output_height,
	OUTPUT_WIDTH,           // output_width,
	INPUT_CHANNEL,          // input_channel, 
	FILTER_HEIGHT,          // filter_height,
	FILTER_WIDTH,           // filter_width, 
	INPUT_HEIGHT,           // input_height, 
	INPUT_WIDTH,            // input_width,
    GROUP, 
    STRIDE,                 // stride,
	NUM_PARAMETER_TYPES,    // num_parameter_types,
};

enum action_type_t {
    LOAD, 
    STORE,
    NUM_ACTION_TYPES,
};

// The component type of neural network accelerator.
enum component_type_t {
    MAC, 
	PE, 
    PE_X,
    PE_Y,
	GLOBAL_BUFFER, 
    CHIPS_X,
    CHIPS_Y,
    DRAM, 
	NUM_COMPONENT_TYPES,
};

static std::vector<std::string> component_type_str __attribute((unused)) = {
    "mac",
	"pe",
    "pe_x",
    "pe_y",
	"global_buffer",
    "chips_x",
    "chips_y",
    "dram",
	"num_component_types",
};

// The MAC type of PE.
enum mac_type_t {
    UNDEFINED_MAC = 0,
    SDSO,
    SVSO, 
    MDMO,
    MVMO,
};

static std::vector<std::string> mac_type_str __attribute((unused)) = {
    "undefined_mac",
    "sdso",
    "svso",
    "mdmo",
    "mvmo",
};

// The memory type of Local buffer and Global buffer.
enum memory_type_t {
    UNDEFINED_MEMORY = 0,
    SEPARATE,
    SHARED,
};

static std::vector<std::string> memory_type_str __attribute((unused)) = {
    "undefined_memory",
    "separate",
    "shared",
};

// The stationary type of dataflow.
// It includes input stationary (IS), weight stationary (WS), and output stationary (OS).
enum stationary_type_t {
    UNDEFINED_STATIONARY = 0,
    INPUT_STATIONARY,
    WEIGHT_STATIONARY, 
    OUTPUT_STATIONARY,
};

static std::vector<std::string> stationary_type_str __attribute((unused)) = {
    "undefined_stationary",
    "input_stationary",
    "weight_stationary",
    "output_stationary",
};

enum compression_type_t {
    DENSE,
    SPARSE_COO, 
    SPARSE_CSC,
    SPARSE_CSR, 
    SPARSEMAP, // Should be renamed
};

static std::vector<std::string> compression_type_str __attribute((unused)) = {
    "dense",
    "sparse_coo",
    "sparse_csc",
    "sparse_csr",
    "sparse_sparten",
};


enum noc_type_t {
    UNDEFINED_NOC = 0,
    MESH,
    BUS, 
    STORE_AND_FORWARD,
    ADDER_TREE,
    CROSSBAR,
};

static std::vector<std::string> noc_type_str __attribute((unused)) = {
    "undefined_noc",
    "mesh",
    "bus",
    "store_and_forward",
    "adder_tree",
    "crossbar",
};


enum layer_name_t {
    UNDEFINED_LAYER = 0,
    AVGPOOL_LAYER,
    CONVOLUTIONAL_LAYER,
    CONNECTED_LAYER,
    MAXPOOL_LAYER,
};


#define is_valid_type(m_vector, m_string) \
    (find(m_vector.begin(), m_vector.end(), m_string.c_str()) != m_vector.end())

#define get_type(m_vector, m_string) \
    distance(m_vector.begin(), find(m_vector.begin(), m_vector.end(), m_string.c_str()))

#endif
