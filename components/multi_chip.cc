#include <iomanip>
#include <cmath>
#include <cstring>
#include "multi_chip.h"
#include "datatype.h"
#include <limits>
#include "interconnect_timing.h"

multi_chip_t::multi_chip_t(section_config_t m_section_config) :
    data(NULL),
    double_buffer(false),
    equal_output_tile(false),
    duplicated_input(0),
    exist_temporal_buffer(true),
    u_transfer_cycle(0.0),
    u_transfer_energy(0.0),
    utilization(0.0),
    write_back_cycle(0.0),
    overlapped_transfer_cycle(0.0),
    stationary_type(stationary_type_t::UNDEFINED_STATIONARY),
    parameter_order("kbpqcrs"),
    nop_type(noc_type_t::UNDEFINED_NOC), 
    memory_type(memory_type_t::UNDEFINED_MEMORY),
    height(0),
    width(0),
    num_chips(1),
    num_active_chips_x(1),
    num_active_chips_y(1),
    frequency(0.0),
    bandwidth(0.0),
    bitwidth(0),
    memory_size(0),
    input_size(0),
    weight_size(0),
    output_size(0),
    initial(true),
    nop_cycle(0.0),
    nop_energy(0.0) {

    init(m_section_config);
}

multi_chip_t::~multi_chip_t() {
    delete [] data;

}

// Initialize the chip-level processor
void multi_chip_t::init(section_config_t m_section_config) {
    
    /* Initialize chip-level processor's specifications */

    // Initialize X and Y dimension of chip-level processors
    m_section_config.get_setting("height", &height);
    m_section_config.get_setting("width", &width);
    if(height == 0 || width == 0 || height > std::numeric_limits<unsigned>::max()/width) {
        std::cerr << "Error: multi_chip height and width must define a non-zero non-overflowing grid" << std::endl;
        exit(1);
    }
    num_chips = height * width;

    // Initialize frequency and bandwidth of chip-level processors.
    m_section_config.get_setting("frequency", &frequency);
    m_section_config.get_setting("bandwidth", &bandwidth);
    // B12: an explicit 'bitwidth' wins silently; a derived width warns when the
    // bandwidth/frequency ratio truncates fractionally.
    unsigned explicit_bitwidth = 0;
    if(m_section_config.get_setting("bitwidth", &explicit_bitwidth)) {
        bitwidth = explicit_bitwidth;
    } else {
        bitwidth = derived_link_bitwidth("multi_chip", bandwidth, frequency);
    }
    if(frequency <= 0.0f || bandwidth < 0.0f || bitwidth == 0) {
        std::cerr << "Error: multi_chip frequency and bitwidth must be positive, and bandwidth non-negative" << std::endl;
        exit(1);
    }

    // Initialize global buffer type (double buffer or single buffer)
    m_section_config.get_setting("double_buffer", &double_buffer);

    // Initialize line size and mask bits of temporal buffer
    line_size.reserve(data_type_t::NUM_DATA_TYPES);
    line_size.assign(data_type_t::NUM_DATA_TYPES, 8);

    mask_bits.reserve(data_type_t::NUM_DATA_TYPES);
    mask_bits.assign(data_type_t::NUM_DATA_TYPES, 0);

    m_section_config.get_vector_setting("line_size", &line_size);

    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        while(line_size[i] > 8) {
            line_size[i] /= 2;
            mask_bits[i]++;
        }
    }
    m_section_config.get_vector_setting("line_size", &line_size);
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        if(!is_valid_memory_line_bits(line_size[i])) {
            std::cerr << "Error: multi_chip line_size must be a power-of-two bit width of at least 8" << std::endl;
            exit(1);
        }
    }

    // Initialize stationary type at the chip-level processor.
    std::string stationary_str;
    if(m_section_config.get_setting("stationary_type", &stationary_str)) {
        stationary_type = (stationary_type_t)get_type(stationary_type_str, stationary_str);
    }

    // Initialize the order of parameters referring to the stationary type
    if(stationary_type == stationary_type_t::INPUT_STATIONARY) {
        parameter_order = "BCPQKRS";
    } else if(stationary_type == stationary_type_t::WEIGHT_STATIONARY) {
        parameter_order = "KCRSBPQ";
    } else if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
        parameter_order = "BKPQKRS";
    }
    m_section_config.get_setting("parameter_order", &parameter_order);

    // Initialize temporal buffer of the chip-level processor.
    m_section_config.get_setting("exist_temporal_buffer", &exist_temporal_buffer);
    std::string memory_type_str;
    m_section_config.get_setting("memory_type", &memory_type_str);
    if(memory_type_str == "shared") {
        memory_type = memory_type_t::SHARED;
        m_section_config.get_setting("memory_size", &memory_size);
        memory_size *= 1024*num_chips;
    }
    else if(memory_type_str == "separate") {
        memory_type = memory_type_t::SEPARATE;
        m_section_config.get_setting("input_size", &input_size);
        m_section_config.get_setting("weight_size", &weight_size);
        m_section_config.get_setting("output_size", &output_size);
        
        input_size *= 1024*num_chips, weight_size *= 1024*num_chips, output_size *= 1024*num_chips;
        memory_size = input_size + weight_size + output_size;
    }
    else {
        std::cerr << "Error: Wrong memory type name : " << memory_type << std::endl;
        exit(1);
    }
    unsigned num_entry = (memory_size + sizeof(data_t) - 1)/sizeof(data_t);

    data = new data_t[num_entry](); 

    std::string nop_str;
    // Existing accelerator configurations name the package-level interconnect
    // "noc". Keep "nop" as the explicit NoP spelling for newer configs.
    if(!m_section_config.get_setting("nop", &nop_str)) {
        m_section_config.get_setting("noc", &nop_str);
    }
    if(!nop_str.empty()) {
        nop_type = (noc_type_t)get_type(noc_type_str, nop_str);
    }
    if(!is_supported_multi_chip_nop(nop_type)) {
        std::cerr << "Error: multi_chip supports only bus or mesh NoP timing models" << std::endl;
        exit(1);
    }

    // Initialize skip transfer signal
    skip_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);

    // Initialize tile size
    tile_size.reserve(data_type_t::NUM_DATA_TYPES);
    tile_size.assign(data_type_t::NUM_DATA_TYPES, 1);

    // Initialize data offsets at the buffer
    offsets.reserve(data_type_t::NUM_DATA_TYPES);
    offsets.assign(data_type_t::NUM_DATA_TYPES, 0);

    /* Initialize signals of chip-level processor */

    // Initialize data exist signal
    exist_data.reserve(data_type_t::NUM_DATA_TYPES);
    exist_data.assign(data_type_t::NUM_DATA_TYPES, false);

    // Initialize request to the off-chip memory signal
    request_to_dram.reserve(data_type_t::NUM_DATA_TYPES);
    request_to_dram.assign(data_type_t::NUM_DATA_TYPES, false);

    // Initialize waiting data from the off-chip memory signal
    is_waiting_data.reserve(data_type_t::NUM_DATA_TYPES);
    is_waiting_data.assign(data_type_t::NUM_DATA_TYPES, false);

    /* Initialize unit cost of chip-level processor */

    // Initial unit read cost
    u_read_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("read_cycle", &u_read_cycle);

    u_read_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("read_energy", &u_read_energy);

    // Initialize unit write cost
    u_write_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("write_cycle", &u_write_cycle);

    u_write_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("write_energy", &u_write_energy);

    // Initialize unit static (leakage) energy of the temporal buffer, pJ/cycle.
    u_static_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("static_energy", &u_static_energy);

    // Initialize NoP cost. Existing configs spell the unit cost "noc_cycle"/"noc_energy"
    // (matching the topology-key fallback above); accept both, "nop_*" taking priority.
    // (Previously only "nop_*" was read, so every shipped config silently ran with a
    // zero NoP unit cost.)
    if(!m_section_config.get_setting("nop_cycle", &nop_cycle)) {
        m_section_config.get_setting("noc_cycle", &nop_cycle);
    }
    if(!m_section_config.get_setting("nop_energy", &nop_energy)) {
        m_section_config.get_setting("noc_energy", &nop_energy);
    }
    // The GLB->multi-chip write-back path charges the same NoP link unit cost
    // (u_transfer_* was previously never assigned).
    u_transfer_cycle = nop_cycle;
    u_transfer_energy = nop_energy;

    /* Initialize stats of chip-level processor */

    // Initialize the number of request the chip-level processor
    num_request.reserve(data_type_t::NUM_DATA_TYPES);
    num_request.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the number of data transfer of chip-level processor
    num_data_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize access cost of chip-level processor
    access_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    access_energy.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize transfer cost between the chip-level processor and the global buffer
    transfer_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    transfer_energy.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize overlapped cycle between the chip-level processor and the global buffer
    cycle_temporal_chips.reserve(data_type_t::NUM_DATA_TYPES);
    cycle_temporal_chips.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    payload_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    metadata_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    storage_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);

    static_energy.reserve(data_type_t::NUM_DATA_TYPES);
    static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
}

// Connect the chip-level processor and the global buffer
void multi_chip_t::connect(std::vector<global_buffer_t*> m_global_buffers) {
    chips = m_global_buffers;
}

// Connect the chip-level processor and the off-chip memory
void multi_chip_t::connect(dram_t *m_dram) {
    dram = m_dram;
}

// Update tile-granular data size
void multi_chip_t::update_tile_size(scheduler_t *m_scheduler) {
    // Update active number of chip-level processors
    num_active_chips_x = m_scheduler->num_active_chips_x;
    num_active_chips_y = m_scheduler->num_active_chips_y;
    if(!has_valid_active_shape(height, width, num_active_chips_y, num_active_chips_x)) {
        std::cerr << "Error: active multi_chip shape exceeds physical chip grid" << std::endl;
        exit(1);
    }
    utilization = (double)(get_number_of_active_chips())/(double)(get_number_of_chips());

    // Update tile size of chip-level processor
    tile_size = m_scheduler->tile_size[component_type_t::CHIPS_Y];

    const size_t required_input = runtime_datatypes().storage_bytes(data_type_t::INPUT, tile_size[data_type_t::INPUT]);
    const size_t required_weight = runtime_datatypes().storage_bytes(data_type_t::WEIGHT, tile_size[data_type_t::WEIGHT]);
    const size_t required_output = runtime_datatypes().storage_bytes(data_type_t::OUTPUT, tile_size[data_type_t::OUTPUT]);
    bool capacity_overflow = false;
    if(memory_type == memory_type_t::SHARED) {
        capacity_overflow = required_input > memory_size || required_weight > memory_size - required_input ||
                            required_output > memory_size - required_input - required_weight;
    } else {
        capacity_overflow = required_input > input_size || required_weight > weight_size || required_output > output_size;
    }
    if(capacity_overflow) {
        std::cerr << "Error: runtime datatype multi-chip tile exceeds temporal-buffer capacity" << std::endl;
        exit(1);
    }

    // Update data offset stored in chip-level processor 
    update_offset();

}

// Update offsets of each data
void multi_chip_t::update_offset() {
    // Update offsets in the case of shared buffer
    if(memory_type == memory_type_t::SHARED) {
        offsets[data_type_t::INPUT] = 0;
        offsets[data_type_t::WEIGHT] = tile_size[data_type_t::INPUT];
        offsets[data_type_t::OUTPUT] = tile_size[data_type_t::INPUT] + tile_size[data_type_t::WEIGHT];
    }
    // Update offsets in the case of separate buffer
    if(memory_type == memory_type_t::SEPARATE) {
        offsets[data_type_t::INPUT] = 0;
        offsets[data_type_t::WEIGHT] = input_size/sizeof(data_t);
        offsets[data_type_t::OUTPUT] = input_size/sizeof(data_t) + weight_size/sizeof(data_t);
    }
}

// Get the stationary type of chip-level processor
stationary_type_t multi_chip_t::get_stationary_type() { return stationary_type; }

// Get parameter order
std::string multi_chip_t::get_parameter_order() { return parameter_order; }

// Get the number of the chip-level processors
unsigned multi_chip_t::get_number_of_chips() { return num_chips; }

// Get the number of active chip-level processors
unsigned multi_chip_t::get_number_of_active_chips() { return num_active_chips_x*num_active_chips_y; }

// A signal that checks whether the PE array is idle state or not.
bool multi_chip_t::is_idle() {
    bool idle = true;
    for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
        idle = chips[i]->is_idle();
        if(idle == false) { break; }
    }
    return idle;
}

// Check whether the request at temporal buffer exists.
bool multi_chip_t::is_exist_request_at_buffer() {
    if(request_to_dram[data_type_t::INPUT] || 
       request_to_dram[data_type_t::WEIGHT] || 
       request_to_dram[data_type_t::OUTPUT]) { 
        return true; 
    }
    else { return false; }
}

// Check whether the request at global buffers.
bool multi_chip_t::is_exist_request_at_global_buffer() {
    bool request_input_ = false, request_weight_ = false, request_output_ = false;

    // Request input data
    for(unsigned i = 0; i < num_active_chips_x*num_active_chips_y; i++) {
        if(chips[i]->request_to_multi_chip[data_type_t::INPUT]) {
            request_input_ = true;
            break;
        }
    }
    // Request weight
    for(unsigned i = 0; i < num_active_chips_x*num_active_chips_y; i++) {
        if(chips[i]->request_to_multi_chip[data_type_t::WEIGHT]) {
            request_weight_ = true;
            break;
        }
    }
    // Request output data
    for(unsigned i = 0; i < num_active_chips_x*num_active_chips_y; i++) {
        if(chips[i]->request_to_multi_chip[data_type_t::OUTPUT]) {
            request_output_ = true;
            break;
        }
    }

    if(request_input_ || request_weight_ || request_output_) { return true; }
    else { return false; }
}

bool multi_chip_t::is_exist_data() {
    if(exist_data[data_type_t::INPUT] && 
       exist_data[data_type_t::WEIGHT] && 
       exist_data[data_type_t::OUTPUT] ) {
        return true;
    }
    else {
        return false;
    }
}

// Wait for the data comes from Global buffer.
bool multi_chip_t::wait_data() {
    if(is_waiting_data[data_type_t::INPUT] || 
       is_waiting_data[data_type_t::WEIGHT] || 
       is_waiting_data[data_type_t::OUTPUT]) {
        return true;
    }
    else {
        return false;
    }
}

// A signal that the data exist in the PE array.
void multi_chip_t::fill_data() {
    for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
        chips[i]->fill_data();
    }
}

void multi_chip_t::account_descriptor_dense_distribution(data_type_t type) {
    num_data_transfer[type]++;
    for(unsigned i = 0; i < get_number_of_active_chips(); ++i) {
        const datatype_transfer_timing_t timing = datatype_transfer_timing(
            type, chips[i]->tile_size[type], line_size[type], chips[i]->line_size[type], bitwidth);
        payload_link_transactions[type] += timing.payload_link_transactions;
        metadata_link_transactions[type] += timing.metadata_link_transactions;
        storage_link_transactions[type] += timing.link_transactions;
        access_energy[type] += timing.source_accesses*u_read_energy[type];
        // GLB fill (write) side: separate accumulator so it can scale with the
        // per-datatype off-chip supply it mirrors (not the full GLB repetition).
        chips[i]->fill_access_energy[type] += timing.destination_accesses*chips[i]->u_write_energy[type];
        if(!chips[i]->double_buffer) {
            access_cycle[type] += timing.source_accesses*u_read_cycle[type];
            chips[i]->fill_access_cycle[type] += timing.destination_accesses*chips[i]->u_write_cycle[type];
        }
        // R2: routed-unicast mesh -- each chip's stream burns per-hop energy per
        // transaction and pays its Manhattan route depth as a one-time pipeline fill
        // (SP1 contract). On the serialized bus hops == 1, so nothing changes.
        const unsigned hops = nop_hops_for_chip(i);
        const double hop_fill = static_cast<double>(hops - 1)*nop_cycle;
        transfer_cycle[type] += timing.link_transactions*nop_cycle + hop_fill;
        transfer_energy[type] += timing.link_transactions*nop_energy*hops;
        if(timing.pipeline_transactions > std::numeric_limits<unsigned>::max()) {
            std::cerr << "Error: multi_chip pipeline transaction count overflow" << std::endl;
            exit(1);
        }
        if(exist_temporal_buffer) {
            cycle_temporal_chips[type] += pipelined_transfer_cycles(
                timing.source_accesses, u_read_cycle[type],
                timing.link_transactions, nop_cycle,
                timing.destination_accesses, chips[i]->u_write_cycle[type]) + hop_fill;
        }
        chips[i]->skip_transfer[type] = false;
    }
}

unsigned multi_chip_t::get_bitwidth() {
    return bitwidth;
}

unsigned multi_chip_t::nop_hops_for_chip(unsigned chip_index) const {
    if(nop_type != noc_type_t::MESH) return 1;
    return nop_route_hops(chip_index, width);
}

void multi_chip_t::account_output_writeback_to_dram() {
    dram->num_request[data_type_t::OUTPUT]++;

    const datatype_transfer_timing_t timing = datatype_transfer_timing(
        data_type_t::OUTPUT, tile_size[data_type_t::OUTPUT], line_size[data_type_t::OUTPUT],
        dram->line_size[data_type_t::OUTPUT], dram->get_bitwidth());
    access_cycle[data_type_t::OUTPUT] += timing.source_accesses*u_read_cycle[data_type_t::OUTPUT];
    access_energy[data_type_t::OUTPUT] += timing.source_accesses*u_read_energy[data_type_t::OUTPUT];
    dram->access_cycle[data_type_t::OUTPUT] += timing.destination_accesses*dram->u_write_cycle[data_type_t::OUTPUT];
    dram->access_energy[data_type_t::OUTPUT] += timing.destination_accesses*dram->u_write_energy[data_type_t::OUTPUT];
    dram->cycle_chip_dram[data_type_t::OUTPUT] += timing.destination_accesses*dram->u_write_cycle[data_type_t::OUTPUT];
    dram->transfer_cycle[data_type_t::OUTPUT] += dram->u_transfer_cycle*timing.link_transactions;
    dram->transfer_energy[data_type_t::OUTPUT] += dram->u_transfer_energy*timing.link_transactions;
    dram->account_row_activations(data_type_t::OUTPUT, tile_size[data_type_t::OUTPUT]);
    dram->payload_link_transactions[data_type_t::OUTPUT] += timing.payload_link_transactions;
    dram->metadata_link_transactions[data_type_t::OUTPUT] += timing.metadata_link_transactions;
    dram->storage_link_transactions[data_type_t::OUTPUT] += timing.link_transactions;
}

// DR6: the in-loop write-back fires on output-tile changes, so the LAST resident
// tile of the live pass never stores. Charged once at layer end, before repetition
// scaling, so each repetition's final tile is covered as well.
void multi_chip_t::flush_output_writeback() {
    account_output_writeback_to_dram();
}

void multi_chip_t::request_data() {
    // Request input data to the off-chip memory
    for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
        if(chips[i]->request_to_multi_chip[data_type_t::INPUT]) {
            request_to_dram[data_type_t::INPUT] = true;
            is_waiting_data[data_type_t::INPUT] = true;

            dram->num_request[data_type_t::INPUT]++;

            break;
        }
    }
    // Request weight to DRAM.
    for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
        if(chips[i]->request_to_multi_chip[data_type_t::WEIGHT]) {
            request_to_dram[data_type_t::WEIGHT] = true;
            is_waiting_data[data_type_t::WEIGHT] = true;

            dram->num_request[data_type_t::WEIGHT]++;

            break;
        }
    }
    // Request output data to DRAM.
    for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
        if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
            if(tile_size[data_type_t::OUTPUT] == dram->tile_size[data_type_t::OUTPUT]) {equal_output_tile = true;}
        }
        else {
            if(tile_size[data_type_t::OUTPUT] == dram->tile_size[data_type_t::OUTPUT] &&
                    tile_size[data_type_t::WEIGHT] != dram->tile_size[data_type_t::WEIGHT]) {equal_output_tile = true;}
        }
        if(chips[i]->request_to_multi_chip[data_type_t::OUTPUT]) {
            if(!initial && !equal_output_tile) {
#ifdef PRINT
                std::cout << "Write back output data from temporal buffer in Multi Chip to DRAM" << std::endl;
#endif
#ifdef FUNCTIONAL
                //m_scheduler->transfer_data((data_t*)layer->output_data, data, m_scheduler->output_offset_dram.front(), multi_chip->offsets[data_type_t::OUTPUT],
                //                           component_type_t::DRAM, component_type_t::CHIPS_Y,
                //                           data_type_t::OUTPUT, multi_chip->get_stationary_type(), action_type_t::STORE);
#endif
                // MC3: charge the DRAM write-back only when an actual write-back occurs.
                // On the initial pass, or when the output tile already matches DRAM's (no
                // eviction), no data is stored, so its cost must not be counted.
                account_output_writeback_to_dram();
            }
            else {
                initial = 0;
            }

            request_to_dram[data_type_t::OUTPUT] = true;
            is_waiting_data[data_type_t::OUTPUT] = true;

            break;
        }
    }
}

void multi_chip_t::data_transfer(scheduler_t *m_scheduler) {
#ifndef FUNCTIONAL
    if(m_scheduler->compression_type != compression_type_t::DENSE) {
        std::cerr << "Error: timing multi_chip supports dense descriptor traffic only" << std::endl;
        exit(1);
    }
#endif

    bool request_to_multi_chip_input = false;
    for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
        if(chips[i]->request_to_multi_chip[data_type_t::INPUT]) {
            request_to_multi_chip_input = true;
            break;
        }
    }

    if(request_to_multi_chip_input) {
#ifndef FUNCTIONAL
        if(!skip_transfer[data_type_t::INPUT]) {
            account_descriptor_dense_distribution(data_type_t::INPUT);
        }
#endif
#ifdef FUNCTIONAL
        // Case 1. Dense data format
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::INPUT]) {
                num_data_transfer[data_type_t::INPUT]++;

                std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
                parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);

                unsigned num_access_global_buffer = 0, num_access_multi_chip = 0;
                for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
                    uint64_t address_global_buffer = 0, address_multi_chip = 0;
                    m_scheduler->transfer_data(chips[i]->data, data, 
                                               chips[i]->offsets[data_type_t::INPUT], offsets[data_type_t::INPUT] + m_scheduler->input_offset_multi_chip[i%m_scheduler->input_offset_multi_chip.size()], 
                                               component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                                               data_type_t::INPUT, chips[i]->get_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //m_scheduler->transfer_data_network_on_chip(chips[i]->data, data, 
                    //                                           component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                    //                                           data_type_t::INPUT, chips[i]->get_stationary_type(), 
                    //                                           action_type_t::LOAD, false);

                    for(unsigned b = 0; b < parameters_global_buffer[parameter_type_t::BATCH_SIZE]; b++) {
                        for(unsigned c = 0; c < parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]; c++) {
                            for(unsigned h = 0; h < parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]; h++) {
                                for(unsigned w = 0; w < parameters_global_buffer[parameter_type_t::INPUT_WIDTH]; w++) {

                                    if(address_global_buffer != ((uint64_t)&chips[i]->data[offsets[data_type_t::INPUT] + 
                                                                                           b*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                                            *parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                                            *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                                           c*parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                                            *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                                           h*parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                           chips[i]->mask_bits[data_type_t::INPUT]) << chips[i]->mask_bits[data_type_t::INPUT]) {

                                        chips[i]->access_energy[data_type_t::INPUT] += chips[i]->u_write_energy[data_type_t::INPUT];
                                        if(!chips[i]->double_buffer) chips[i]->access_cycle[data_type_t::INPUT] += chips[i]->u_write_cycle[data_type_t::INPUT];

                                        num_access_global_buffer++;

                                        address_global_buffer = ((uint64_t)&chips[i]->data[offsets[data_type_t::INPUT] + 
                                                                                           b*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                                            *parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                                            *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                                           c*parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                                            *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                                           h*parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                           chips[i]->mask_bits[data_type_t::INPUT]) << chips[i]->mask_bits[data_type_t::INPUT];
                                    }
                                    if(!m_scheduler->read_tile_granular_chip_input[i%m_scheduler->input_offset_multi_chip.size()] &&
                                       address_multi_chip != ((uint64_t)&data[offsets[data_type_t::INPUT] + 
                                                                              m_scheduler->input_offset_multi_chip[i%m_scheduler->input_offset_multi_chip.size()] + 
                                                                              b*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                               *parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                               *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                              c*parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                               *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                              h*parameters_multi_chip[parameter_type_t::INPUT_HEIGHT] + w] >> 
                                                                              mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT]) {

                                        access_energy[data_type_t::INPUT] += u_read_energy[data_type_t::INPUT];
                                        if(!chips[i]->double_buffer) access_cycle[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT];
                                        num_access_multi_chip++;
                            
                                        // Update address of temporal buffer in on-chip processors
                                        m_scheduler->read_tile_granular_chip_input[i%m_scheduler->input_offset_multi_chip.size()] = true;
                                        address_multi_chip = ((uint64_t)&data[offsets[data_type_t::INPUT] + 
                                                                              m_scheduler->input_offset_multi_chip[i%m_scheduler->input_offset_multi_chip.size()] + 
                                                                              b*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                               *parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                               *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                              c*parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                               *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                              h*parameters_multi_chip[parameter_type_t::INPUT_HEIGHT] + w] >> 
                                                                              mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT];
                                    }
                                }
                            }
                        }
                    }
                    chips[i]->skip_transfer[data_type_t::INPUT] = false;
                }

                // Update transfer cycle and energy
                transfer_cycle[data_type_t::INPUT] += num_access_multi_chip*nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::INPUT])/(double)(bitwidth));
                transfer_energy[data_type_t::INPUT] += num_access_multi_chip*nop_energy*ceil((double)(chips[0]->line_size[data_type_t::INPUT])/(double)(bitwidth));

                // At the 1, 2, before last, and last stages
                double first_stage = u_read_cycle[data_type_t::INPUT];
                double second_stage = std::max(u_read_cycle[data_type_t::INPUT], 
                                               nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::INPUT])/(double)(bitwidth)));

                double last_before_stage = std::max(chips[get_number_of_active_chips()-1]->u_write_cycle[data_type_t::INPUT], 
                                                    nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::INPUT])/(double)(bitwidth)));
                double last_stage = chips[get_number_of_active_chips()-1]->u_write_cycle[data_type_t::INPUT];

                // Other stages
                double other_stage = std::max(u_read_cycle[data_type_t::INPUT],
                                     std::max(chips[0]->u_write_cycle[data_type_t::INPUT],
                                              nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::INPUT])/(double)(bitwidth))));

                // If temporal buffer is exist at the chip-level processors.
                if(exist_temporal_buffer) {
                    if(num_access_multi_chip == 0) { } else if(num_access_multi_chip == 1) {
                        cycle_temporal_chips[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT] +
                                                                    nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::INPUT])/(double)(bitwidth)) + 
                                                                    chips[0]->u_write_cycle[data_type_t::INPUT];
                    } else {
                        cycle_temporal_chips[data_type_t::INPUT] += first_stage + second_stage +
                                                                    (num_access_multi_chip - 2)*other_stage + 
                                                                    last_before_stage + last_stage;
                    }
                }
            }
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
            std::cerr << "Current version does not support COO format for sparse data" << std::endl;
            exit(1);
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::INPUT]) {
                // Row bit calculation
                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
                unsigned row = parameters[parameter_type_t::INPUT_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer[data_type_t::INPUT]++;

                for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
                    m_scheduler->transfer_data(chips[i]->data, data, 
                                               chips[i]->offsets[data_type_t::INPUT], offsets[data_type_t::INPUT] + m_scheduler->input_offset_multi_chip[i%m_scheduler->input_offset_multi_chip.size()], 
                                               component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                                               data_type_t::INPUT, chips[i]->get_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //m_scheduler->transfer_data_network_on_chip(chips[i]->data, data, 
                    //                                           component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                    //                                           data_type_t::INPUT, chips[i]->get_stationary_type(), 
                    //                                           action_type_t::LOAD, false);
                    if(!m_scheduler->read_tile_granular_chip_input[i%m_scheduler->input_offset_multi_chip.size()]) {
                        // Update chip-level processor access cost
                        access_cycle[data_type_t::INPUT] += (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                           *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                            (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                           *u_read_cycle[data_type_t::INPUT]
                                                           /(sizeof(data_t)*8/row_bit) + // Row index
                                                            parameters[parameter_type_t::BATCH_SIZE]
                                                           *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                           *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                           *u_read_cycle[data_type_t::INPUT]
                                                           /(sizeof(data_t)*8/row_bit); // Column pointer
                        access_energy[data_type_t::INPUT] += (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_energy[data_type_t::INPUT] + // Non-zero data
                                                             (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_energy[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/row_bit) + // Row index
                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                            *u_read_energy[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/row_bit); // Column pointer

                        m_scheduler->read_tile_granular_chip_input[i%m_scheduler->input_offset_multi_chip.size()] = true;

                        transfer_cycle[data_type_t::INPUT] += nop_cycle
                                                             *ceil((double)(chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *8*sizeof(data_t)
                                                             /(double)bitwidth) + // Non-zero data
                                                              nop_cycle
                                                             *ceil((double)((chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)
                                                             /(double)(bitwidth)) + // Row index
                                                              nop_cycle
                                                             *ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                     *(parameters[parameter_type_t::INPUT_WIDTH]+1)*row_bit)
                                                             /(double)(bitwidth)); // Column pointer
                    }

                    // Update global buffer access cost
                    chips[i]->access_cycle[data_type_t::INPUT] += (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                 *chips[i]->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                                  (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                 *chips[i]->u_write_cycle[data_type_t::INPUT]
                                                                 /(sizeof(data_t)*8/row_bit) + // Row index
                                                                  parameters[parameter_type_t::BATCH_SIZE]
                                                                 *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                 *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                                 *chips[i]->u_write_cycle[data_type_t::INPUT]
                                                                 /(sizeof(data_t)*8/row_bit); // Column pointer
                    chips[i]->access_energy[data_type_t::INPUT] += (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                  *chips[i]->u_write_energy[data_type_t::INPUT] + // Non-zero data
                                                                   (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                  *chips[i]->u_write_energy[data_type_t::INPUT]
                                                                  /(sizeof(data_t)*8/row_bit) + // Row index
                                                                   parameters[parameter_type_t::BATCH_SIZE]
                                                                  *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                  *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                                  *chips[i]->u_write_energy[data_type_t::INPUT]
                                                                  /(sizeof(data_t)*8/row_bit); // Column pointer


                    transfer_energy[data_type_t::INPUT] += nop_energy
                                                         *ceil((double)(chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                         *8*sizeof(data_t)
                                                         /(double)bitwidth) + // Non-zero data
                                                          nop_energy
                                                         *ceil((double)((chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)
                                                         /(double)(bitwidth)) + // Row index
                                                          nop_energy
                                                         *ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                 *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                 *(parameters[parameter_type_t::INPUT_WIDTH]+1)*row_bit)
                                                         /(double)(bitwidth)); // Column pointer

                    double data_size = (double)((chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*sizeof(data_t))
                                     + (double)((chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)/8.0
                                     + (double)(parameters[parameter_type_t::BATCH_SIZE]*parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]*(parameters[parameter_type_t::INPUT_WIDTH]+1)*row_bit)/8.0;
                    chips[i]->utilization[data_type_t::INPUT] = std::max(chips[i]->utilization[data_type_t::INPUT], data_size/(double)chips[i]->get_buffer_size());

                    chips[i]->skip_transfer[data_type_t::INPUT] = false;
                }
            }
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::INPUT]) {
                // Column bit calculation
                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
                unsigned column = parameters[parameter_type_t::INPUT_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer[data_type_t::INPUT]++;

                for(unsigned i = 0; i < get_number_of_active_chips(); i++) {

                    // Transfer input data from the chip-level processor to the global buffer
                    m_scheduler->transfer_data(chips[i]->data, data, 
                                               chips[i]->offsets[data_type_t::INPUT], offsets[data_type_t::INPUT] + m_scheduler->input_offset_multi_chip[i%m_scheduler->input_offset_multi_chip.size()], 
                                               component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                                               data_type_t::INPUT, chips[i]->get_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //m_scheduler->transfer_data_network_on_chip(chips[i]->data, data, 
                    //                                           component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                    //                                           data_type_t::INPUT, chips[i]->get_stationary_type(), 
                    //                                           action_type_t::LOAD, false);
                    
                    if(!m_scheduler->read_tile_granular_chip_input[i%m_scheduler->input_offset_multi_chip.size()]) {
                        // Update chip-level processor access cost
                        access_cycle[data_type_t::INPUT] += (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                           *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                            (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                           *u_read_cycle[data_type_t::INPUT]
                                                           /(sizeof(data_t)*8/column_bit) + // Column index
                                                            parameters[parameter_type_t::BATCH_SIZE]
                                                           *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                           *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                           *u_read_cycle[data_type_t::INPUT]
                                                           /(sizeof(data_t)*8/column_bit); // Row pointer
                        access_energy[data_type_t::INPUT] += (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_energy[data_type_t::INPUT] + // Non-zero data
                                                             (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_energy[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/column_bit) + // Column index
                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                            *u_read_energy[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/column_bit); // Row pointer

                        m_scheduler->read_tile_granular_chip_input[i%m_scheduler->input_offset_multi_chip.size()] = true;

                        transfer_cycle[data_type_t::INPUT] += nop_cycle
                                                             *ceil((double)(chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *8*sizeof(data_t)
                                                             /(double)bitwidth) + // Non-zero data
                                                              nop_cycle
                                                             *ceil((double)((chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)
                                                             /(double)(bitwidth)) + // Column index
                                                              nop_cycle
                                                             *ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                     *(parameters[parameter_type_t::INPUT_HEIGHT]+1)*column_bit)
                                                             /(double)(bitwidth)); // Row pointer
                    }

                    // Update global buffer access cost
                    chips[i]->access_cycle[data_type_t::INPUT] += (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                 *chips[i]->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                                  (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                 *chips[i]->u_write_cycle[data_type_t::INPUT]
                                                                 /(sizeof(data_t)*8/column_bit) + // Column index
                                                                  parameters[parameter_type_t::BATCH_SIZE]
                                                                 *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                 *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                                 *chips[i]->u_write_cycle[data_type_t::INPUT]
                                                                 /(sizeof(data_t)*8/column_bit); // Row pointer
                    chips[i]->access_energy[data_type_t::INPUT] += (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                  *chips[i]->u_write_energy[data_type_t::INPUT] + // Non-zero data
                                                                   (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                  *chips[i]->u_write_energy[data_type_t::INPUT]
                                                                  /(sizeof(data_t)*8/column_bit) + // Column index
                                                                   parameters[parameter_type_t::BATCH_SIZE]
                                                                  *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                  *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                                  *chips[i]->u_write_energy[data_type_t::INPUT]
                                                                  /(sizeof(data_t)*8/column_bit); // Row pointer

                    transfer_energy[data_type_t::INPUT] += nop_energy
                                                         *ceil((double)(chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                         *8*sizeof(data_t)
                                                         /(double)bitwidth) + // Non-zero data
                                                          nop_energy
                                                         *ceil((double)((chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)
                                                         /(double)(bitwidth)) + // Column index
                                                          nop_energy
                                                         *ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                 *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                 *(parameters[parameter_type_t::INPUT_HEIGHT]+1)*column_bit)
                                                         /(double)(bitwidth)); // Row pointer

                    double data_size = (double)((chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*sizeof(data_t))
                                     + (double)((chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)/8.0
                                     + (double)(parameters[parameter_type_t::BATCH_SIZE]*parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]*(parameters[parameter_type_t::INPUT_HEIGHT]+1)*column_bit)/8.0;
                    chips[i]->utilization[data_type_t::INPUT] = std::max(chips[i]->utilization[data_type_t::INPUT], data_size/(double)chips[i]->get_buffer_size());

                    chips[i]->skip_transfer[data_type_t::INPUT] = false;
                }
            }
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::INPUT]) {
                num_data_transfer[data_type_t::INPUT]++;
                for(unsigned i = 0; i < get_number_of_active_chips(); i++) {

                    // Transfer input data from the chip-level processor to the global buffer
                    m_scheduler->transfer_data(chips[i]->data, data, 
                                               chips[i]->offsets[data_type_t::INPUT], offsets[data_type_t::INPUT] + m_scheduler->input_offset_multi_chip[i%m_scheduler->input_offset_multi_chip.size()], 
                                               component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                                               data_type_t::INPUT, chips[i]->get_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //m_scheduler->transfer_data_network_on_chip(chips[i]->data, data, 
                    //                                           component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                    //                                           data_type_t::INPUT, chips[i]->get_stationary_type(), 
                    //                                           action_type_t::LOAD, false);

                    if(!m_scheduler->read_tile_granular_chip_input[i%m_scheduler->input_offset_multi_chip.size()]) {
                        // Update access cost of chip-level processor
                        access_cycle[data_type_t::INPUT] += (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                           *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                            chips[i]->tile_size[data_type_t::INPUT]
                                                           *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
                   
                        access_energy[data_type_t::INPUT] += (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_energy[data_type_t::INPUT] + // Non-zero data
                                                             chips[i]->tile_size[data_type_t::INPUT]
                                                            *u_read_energy[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8); // Metadata

                        m_scheduler->read_tile_granular_chip_input[i%m_scheduler->input_offset_multi_chip.size()] = true;

                        // Update transfer cycle between the chip-level processor and the global buffer
                        transfer_cycle[data_type_t::INPUT] += nop_cycle
                                                             *ceil((double)((chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *8*sizeof(data_t))
                                                             /(double)bitwidth) +  // Non-zero data
                                                              nop_cycle 
                                                             *ceil((double)(chips[i]->tile_size[data_type_t::INPUT])/(double)bitwidth);
                    }

                    // Update global buffer access cost
                    chips[i]->access_cycle[data_type_t::INPUT] += (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                 *chips[i]->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                                  chips[i]->tile_size[data_type_t::INPUT]
                                                                 *chips[i]->u_write_cycle[data_type_t::INPUT]
                                                                 /(sizeof(data_t)*8); // Metadata
                   
                    chips[i]->access_energy[data_type_t::INPUT] += (chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                  *chips[i]->u_write_energy[data_type_t::INPUT] +  // Non-zero data
                                                                   chips[i]->tile_size[data_type_t::INPUT]
                                                                  *chips[i]->u_write_energy[data_type_t::INPUT]
                                                                  /(sizeof(data_t)*8); // Metadata

                    transfer_energy[data_type_t::INPUT] += nop_energy
                                                          *ceil((double)((chips[i]->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                          *8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                           nop_energy
                                                          *ceil((double)(chips[i]->tile_size[data_type_t::INPUT])/(double)bitwidth); // Metadata

                    chips[i]->skip_transfer[data_type_t::INPUT] = false;
                }
            }
        }
#elif defined(NPUSIM_LEGACY_ADDRESS_TIMING)
        if(!skip_transfer[data_type_t::INPUT]) {
            num_data_transfer[data_type_t::INPUT]++;

            std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
            parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);

            unsigned num_access_global_buffer = 0, num_access_multi_chip = 0;
            for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
                uint64_t address_global_buffer = 0, address_multi_chip = 0;

                for(unsigned b = 0; b < parameters_global_buffer[parameter_type_t::BATCH_SIZE]; b++) {
                    for(unsigned c = 0; c < parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned h = 0; h < parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]; h++) {
                            for(unsigned w = 0; w < parameters_global_buffer[parameter_type_t::INPUT_WIDTH]; w++) {

                                if(address_global_buffer != ((uint64_t)&chips[i]->data[offsets[data_type_t::INPUT] + 
                                                                                       b*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                                        *parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                                        *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                                       c*parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                                        *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                                       h*parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                       chips[i]->mask_bits[data_type_t::INPUT]) << chips[i]->mask_bits[data_type_t::INPUT]) {

                                    chips[i]->access_energy[data_type_t::INPUT] += chips[i]->u_write_energy[data_type_t::INPUT];
                                    if(!chips[i]->double_buffer) chips[i]->access_cycle[data_type_t::INPUT] += chips[i]->u_write_cycle[data_type_t::INPUT];

                                    num_access_global_buffer++;

                                    address_global_buffer = ((uint64_t)&chips[i]->data[offsets[data_type_t::INPUT] + 
                                                                                       b*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                                        *parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                                        *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                                       c*parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                                        *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                                       h*parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                       chips[i]->mask_bits[data_type_t::INPUT]) << chips[i]->mask_bits[data_type_t::INPUT];
                                }
                                if(!m_scheduler->read_tile_granular_chip_input[i%m_scheduler->input_offset_multi_chip.size()] &&
                                   address_multi_chip != ((uint64_t)&data[offsets[data_type_t::INPUT] + 
                                                                          m_scheduler->input_offset_multi_chip[i%m_scheduler->input_offset_multi_chip.size()] + 
                                                                          b*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                           *parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                           *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                          c*parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                           *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                          h*parameters_multi_chip[parameter_type_t::INPUT_HEIGHT] + w] >> 
                                                                          mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT]) {

                                    access_energy[data_type_t::INPUT] += u_read_energy[data_type_t::INPUT];
                                    if(!chips[i]->double_buffer) access_cycle[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT];
                                    num_access_multi_chip++;
                        
                                    // Update address of temporal buffer in on-chip processors
                                    m_scheduler->read_tile_granular_chip_input[i%m_scheduler->input_offset_multi_chip.size()] = true;
                                    address_multi_chip = ((uint64_t)&data[offsets[data_type_t::INPUT] + 
                                                                          m_scheduler->input_offset_multi_chip[i%m_scheduler->input_offset_multi_chip.size()] + 
                                                                          b*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                           *parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                           *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                          c*parameters_multi_chip[parameter_type_t::INPUT_HEIGHT]
                                                                           *parameters_multi_chip[parameter_type_t::INPUT_WIDTH] + 
                                                                          h*parameters_multi_chip[parameter_type_t::INPUT_HEIGHT] + w] >> 
                                                                          mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT];
                                }
                            }
                        }
                    }
                }
                chips[i]->skip_transfer[data_type_t::INPUT] = false;
            }

            // Update transfer cycle and energy
            transfer_cycle[data_type_t::INPUT] += num_access_multi_chip*nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::INPUT])/(double)(bitwidth));
            transfer_energy[data_type_t::INPUT] += num_access_multi_chip*nop_energy*ceil((double)(chips[0]->line_size[data_type_t::INPUT])/(double)(bitwidth));

            // At the 1, 2, before last, and last stages
            double first_stage = u_read_cycle[data_type_t::INPUT];
            double second_stage = std::max(u_read_cycle[data_type_t::INPUT], 
                                           nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::INPUT])/(double)(bitwidth)));

            double last_before_stage = std::max(chips[get_number_of_active_chips()-1]->u_write_cycle[data_type_t::INPUT], 
                                                nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::INPUT])/(double)(bitwidth)));
            double last_stage = chips[get_number_of_active_chips()-1]->u_write_cycle[data_type_t::INPUT];

            // Other stages
            double other_stage = std::max(u_read_cycle[data_type_t::INPUT],
                                 std::max(chips[0]->u_write_cycle[data_type_t::INPUT],
                                          nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::INPUT])/(double)(bitwidth))));

            // If temporal buffer is exist at the chip-level processors.
            if(exist_temporal_buffer) {
                if(num_access_multi_chip == 0) { } else if(num_access_multi_chip == 1) {
                    cycle_temporal_chips[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT] +
                                                                nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::INPUT])/(double)(bitwidth)) + 
                                                                chips[0]->u_write_cycle[data_type_t::INPUT];
                } else {
                    cycle_temporal_chips[data_type_t::INPUT] += first_stage + second_stage +
                                                                (num_access_multi_chip - 2)*other_stage + 
                                                                last_before_stage + last_stage;
                }
            }
        }
#endif
        for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
            chips[i]->exist_data[data_type_t::INPUT] = true, chips[i]->request_to_multi_chip[data_type_t::INPUT] = false;
        }
        is_waiting_data[data_type_t::INPUT] = false;
    }

    bool request_to_multi_chip_weight = false;
    for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
        if(chips[i]->request_to_multi_chip[data_type_t::WEIGHT]) {
            request_to_multi_chip_weight = true;
            break;
        }
    }

    if(request_to_multi_chip_weight) {
#ifndef FUNCTIONAL
        if(!skip_transfer[data_type_t::WEIGHT]) {
            account_descriptor_dense_distribution(data_type_t::WEIGHT);
        }
#endif
#ifdef FUNCTIONAL
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                num_data_transfer[data_type_t::WEIGHT]++;

                std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
                parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);

                unsigned num_access_global_buffer = 0, num_access_multi_chip = 0;
                for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
                    uint64_t address_global_buffer = 0, address_multi_chip = 0;
                    
                    m_scheduler->transfer_data(chips[i]->data, data, 
                                               chips[i]->offsets[data_type_t::WEIGHT], offsets[data_type_t::WEIGHT] + m_scheduler->weight_offset_multi_chip[i%m_scheduler->weight_offset_multi_chip.size()], 
                                               component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                                               data_type_t::WEIGHT, chips[i]->get_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //m_scheduler->transfer_data_network_on_chip(chips[i]->data, data, 
                    //                                           component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                    //                                           data_type_t::WEIGHT, chips[i]->get_stationary_type(), 
                    //                                           action_type_t::LOAD, false);

                    for(unsigned k = 0; k < parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                        for(unsigned c = 0; c < parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]; c++) {
                            for(unsigned r = 0; r < parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]; r++) {
                                for(unsigned s = 0; s < parameters_global_buffer[parameter_type_t::FILTER_WIDTH]; s++) {

                                    if(address_global_buffer != ((uint64_t)&chips[i]->data[offsets[data_type_t::WEIGHT] + 
                                                                                           k*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                                            *parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                                            *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                                           c*parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                                            *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                                           r*parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                           chips[i]->mask_bits[data_type_t::WEIGHT]) << chips[i]->mask_bits[data_type_t::WEIGHT]) {

                                        // Update costs of global buffer
                                        chips[i]->access_energy[data_type_t::WEIGHT] += chips[i]->u_write_energy[data_type_t::WEIGHT];
                                        if(!chips[i]->double_buffer) chips[i]->access_cycle[data_type_t::WEIGHT] += chips[i]->u_write_cycle[data_type_t::WEIGHT];
                                        num_access_global_buffer++;

                                        address_global_buffer = ((uint64_t)&chips[i]->data[offsets[data_type_t::INPUT] + 
                                                                                           k*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                                            *parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                                            *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                                           c*parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                                            *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                                           r*parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                           chips[i]->mask_bits[data_type_t::WEIGHT]) << chips[i]->mask_bits[data_type_t::WEIGHT];
                                    }
                                    if(!m_scheduler->read_tile_granular_chip_weight[i%m_scheduler->weight_offset_multi_chip.size()] &&
                                       address_multi_chip != ((uint64_t)&data[offsets[data_type_t::WEIGHT] + 
                                                                              m_scheduler->weight_offset_multi_chip[i%m_scheduler->weight_offset_multi_chip.size()] + 
                                                                              k*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                               *parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                               *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                              c*parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                               *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                              r*parameters_multi_chip[parameter_type_t::FILTER_HEIGHT] + s] >> 
                                                                              mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT]) {

                                        // Update cost of temporal buffer in chip-level processors
                                        access_energy[data_type_t::WEIGHT] += u_read_energy[data_type_t::WEIGHT];
                                        if(!chips[i]->double_buffer) access_cycle[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT];
                                        num_access_multi_chip++;

                                        // Update address of temporal buffer in chip-level processors
                                        m_scheduler->read_tile_granular_chip_weight[i%m_scheduler->weight_offset_multi_chip.size()] = true;
                                        address_multi_chip = ((uint64_t)&data[offsets[data_type_t::WEIGHT] + 
                                                                              m_scheduler->weight_offset_multi_chip[i%m_scheduler->weight_offset_multi_chip.size()] + 
                                                                              k*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                               *parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                               *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                              c*parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                               *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                              r*parameters_multi_chip[parameter_type_t::FILTER_HEIGHT] + s] >> 
                                                                              mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT];

                                    }
                                }
                            }
                        }
                    }
                    chips[i]->skip_transfer[data_type_t::WEIGHT] = false;
                }

                transfer_cycle[data_type_t::WEIGHT] += num_access_multi_chip*nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::WEIGHT])/(double)(bitwidth));
                transfer_energy[data_type_t::WEIGHT] += num_access_multi_chip*nop_energy*ceil((double)(chips[0]->line_size[data_type_t::WEIGHT])/(double)(bitwidth));

                double first_stage = u_read_cycle[data_type_t::WEIGHT];
                double second_stage = std::max(u_read_cycle[data_type_t::WEIGHT], 
                                               nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::WEIGHT])/(double)(bitwidth)));

                double last_before_stage = std::max(chips[get_number_of_active_chips()-1]->u_write_cycle[data_type_t::WEIGHT], 
                                                    nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
                double last_stage = chips[get_number_of_active_chips()-1]->u_write_cycle[data_type_t::WEIGHT];

                double other_stage = std::max(u_read_cycle[data_type_t::WEIGHT], 
                                     std::max(chips[0]->u_write_cycle[data_type_t::WEIGHT], 
                                              nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::WEIGHT])/(double)(bitwidth))));
                
                if(exist_temporal_buffer) {
                    if(num_access_multi_chip == 0) { } else if(num_access_multi_chip == 1) {
                        cycle_temporal_chips[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT] + 
                                                                     nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::WEIGHT])/(double)(bitwidth)) + 
                                                                     chips[0]->u_write_cycle[data_type_t::WEIGHT];
                    } else {
                        cycle_temporal_chips[data_type_t::WEIGHT] += first_stage + second_stage + 
                                                                     (num_access_multi_chip - 2)*other_stage + 
                                                                     last_before_stage + last_stage;
                    }
                }
            }
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
            std::cerr << "Current version does not support COO format for sparse data" << std::endl;
            exit(1);
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                // Row bit calculation
                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
                unsigned row = parameters[parameter_type_t::FILTER_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer[data_type_t::WEIGHT]++;

                for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
                    // Transfer weight from chip-level processor to the global buffer
                    m_scheduler->transfer_data(chips[i]->data, data, 
                                               chips[i]->offsets[data_type_t::WEIGHT], offsets[data_type_t::WEIGHT] + m_scheduler->weight_offset_multi_chip[i%m_scheduler->weight_offset_multi_chip.size()], 
                                               component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                                               data_type_t::WEIGHT, chips[i]->get_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //m_scheduler->transfer_data_network_on_chip(chips[i]->data, data, 
                    //                                           component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                    //                                           data_type_t::WEIGHT, chips[i]->get_stationary_type(), 
                    //                                           action_type_t::LOAD, false);
                    if(!m_scheduler->read_tile_granular_chip_weight[i%m_scheduler->weight_offset_multi_chip.size()]) {
                        // Update chip-level processor access cost
                        access_cycle[data_type_t::WEIGHT] += (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                            *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                             (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                            *u_read_cycle[data_type_t::WEIGHT]
                                                            /(sizeof(data_t)*8/row_bit) + // Row index
                                                             parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                            *u_read_cycle[data_type_t::WEIGHT]
                                                            /(sizeof(data_t)*8/row_bit); // Column pointer
                        access_energy[data_type_t::WEIGHT] += (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                              (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_energy[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/row_bit) + // Row index
                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                             *u_read_energy[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/row_bit); // Column pointer

                        m_scheduler->read_tile_granular_chip_weight[i%m_scheduler->weight_offset_multi_chip.size()] = true;

                        transfer_cycle[data_type_t::WEIGHT] += nop_cycle
                                                              *ceil((double)(chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *8*sizeof(data_t)/(double)bitwidth) + // Non-zero data
                                                               nop_cycle
                                                              *ceil((double)((chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *row_bit)
                                                              /(double)bitwidth) + // Row index
                                                               nop_cycle
                                                              *ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                            *(parameters[parameter_type_t::FILTER_WIDTH]+1)*row_bit)
                                                              /(double)bitwidth); // Column pointer
                                                              
                    }
                    // Update global buffer access cost
                    chips[i]->access_cycle[data_type_t::WEIGHT] += (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                  *chips[i]->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                   (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                  *chips[i]->u_write_cycle[data_type_t::WEIGHT]
                                                                  /(sizeof(data_t)*8/row_bit) + // Row index
                                                                   parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                  *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                  *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                                  *chips[i]->u_write_cycle[data_type_t::WEIGHT]
                                                                  /(sizeof(data_t)*8/row_bit); // Column pointer
                    chips[i]->access_energy[data_type_t::WEIGHT] += (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                   *chips[i]->u_write_energy[data_type_t::WEIGHT] + // Non-zero data
                                                                    (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                   *chips[i]->u_write_energy[data_type_t::WEIGHT]
                                                                   /(sizeof(data_t)*8/row_bit) + // Row index
                                                                    parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                   *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                                   *chips[i]->u_write_energy[data_type_t::WEIGHT]
                                                                   /(sizeof(data_t)*8/row_bit); // Column pointer
                    
                    // Update transfer between the chip-level processor and the global buffer 
                    transfer_energy[data_type_t::WEIGHT] += nop_energy
                                                           *ceil((double)(chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                           *8*sizeof(data_t)
                                                           /(double)bitwidth) + // Non-zero data
                                                            nop_energy
                                                           *ceil((double)((chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*row_bit)
                                                           /(double)bitwidth) + // Row index
                                                            nop_energy
                                                           *ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *(parameters[parameter_type_t::FILTER_WIDTH]+1)*row_bit)
                                                           /(double)bitwidth); // Column pointer

                    chips[i]->skip_transfer[data_type_t::WEIGHT] = false;
                }
            }
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                // Column bit calculation
                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
                unsigned column = parameters[parameter_type_t::FILTER_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer[data_type_t::WEIGHT]++;

                for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
                    // Transfer weight from chip-level processor to the global buffer
                    m_scheduler->transfer_data(chips[i]->data, data, 
                                               chips[i]->offsets[data_type_t::WEIGHT], offsets[data_type_t::WEIGHT] + m_scheduler->weight_offset_multi_chip[i%m_scheduler->weight_offset_multi_chip.size()], 
                                               component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                                               data_type_t::WEIGHT, chips[i]->get_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //m_scheduler->transfer_data_network_on_chip(chips[i]->data, data, 
                    //                                           component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                    //                                           data_type_t::WEIGHT, chips[i]->get_stationary_type(), 
                    //                                           action_type_t::LOAD, false);
                    if(!m_scheduler->read_tile_granular_chip_weight[i%m_scheduler->weight_offset_multi_chip.size()]) {
                        // Update chip-level processor access cost
                        access_cycle[data_type_t::WEIGHT] += (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                            *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                             (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                            *u_read_cycle[data_type_t::WEIGHT]
                                                            /(sizeof(data_t)*8/column_bit) + // Column index
                                                             parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                            *u_read_cycle[data_type_t::WEIGHT]
                                                            /(sizeof(data_t)*8/column_bit); // Row pointer
                        access_energy[data_type_t::WEIGHT] += (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                              (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_energy[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/column_bit) + // Column index
                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                             *u_read_energy[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/column_bit); // Row pointer

                        m_scheduler->read_tile_granular_chip_weight[i%m_scheduler->weight_offset_multi_chip.size()] = true;

                        transfer_cycle[data_type_t::WEIGHT] += nop_cycle
                                                              *ceil((double)(chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *8*sizeof(data_t)/(double)bitwidth) + // Non-zero data
                                                               nop_cycle
                                                              *ceil((double)((chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *column_bit)
                                                              /(double)bitwidth) + // Column index
                                                               nop_cycle
                                                              *ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                            *(parameters[parameter_type_t::FILTER_HEIGHT]+1)*column_bit)
                                                              /(double)bitwidth); // Row pointer
                                                              
                    }
                    // Update global buffer access cost
                    chips[i]->access_cycle[data_type_t::WEIGHT] += (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                  *chips[i]->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                   (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                  *chips[i]->u_write_cycle[data_type_t::WEIGHT]
                                                                  /(sizeof(data_t)*8/column_bit) + // Column index
                                                                   parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                  *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                  *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                                  *chips[i]->u_write_cycle[data_type_t::WEIGHT]
                                                                  /(sizeof(data_t)*8/column_bit); // Row pointer
                    chips[i]->access_energy[data_type_t::WEIGHT] += (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                   *chips[i]->u_write_energy[data_type_t::WEIGHT] + // Non-zero data
                                                                    (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                   *chips[i]->u_write_energy[data_type_t::WEIGHT]
                                                                   /(sizeof(data_t)*8/column_bit) + // Column index
                                                                    parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                   *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                                   *chips[i]->u_write_energy[data_type_t::WEIGHT]
                                                                   /(sizeof(data_t)*8/column_bit); // Row pointer
                    
                    // Update transfer between the chip-level processor and the global buffer 
                    transfer_energy[data_type_t::WEIGHT] += nop_energy
                                                           *ceil((double)(chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                           *8*sizeof(data_t)
                                                           /(double)bitwidth) + // Non-zero data
                                                            nop_energy
                                                           *ceil((double)((chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*column_bit)
                                                           /(double)bitwidth) + // Column index
                                                            nop_energy
                                                           *ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *(parameters[parameter_type_t::FILTER_HEIGHT]+1)*column_bit)
                                                           /(double)bitwidth); // Row pointer

                    chips[i]->skip_transfer[data_type_t::WEIGHT] = false;
                }
            }
        }
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
            
                num_data_transfer[data_type_t::WEIGHT]++;
                for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
                    // Transfer weight from chip-level processor to the global buffer
                    m_scheduler->transfer_data(chips[i]->data, data, 
                                               chips[i]->offsets[data_type_t::WEIGHT], offsets[data_type_t::WEIGHT] + m_scheduler->weight_offset_multi_chip[i%m_scheduler->weight_offset_multi_chip.size()], 
                                               component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                                               data_type_t::WEIGHT, chips[i]->get_stationary_type(), action_type_t::LOAD);
                    // Update for NPUsim ver2
                    //m_scheduler->transfer_data_network_on_chip(chips[i]->data, data, 
                    //                                           component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                    //                                           data_type_t::WEIGHT, chips[i]->get_stationary_type(), 
                    //                                           action_type_t::LOAD, false);
                   
                    if(!m_scheduler->read_tile_granular_chip_weight[i%m_scheduler->weight_offset_multi_chip.size()]) { 
                        // Update access cost of chip-level processor
                        access_cycle[data_type_t::WEIGHT] += (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                            *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                             chips[i]->tile_size[data_type_t::WEIGHT]
                                                            *u_read_cycle[data_type_t::WEIGHT]
                                                            /(sizeof(data_t)*8); // Metadata
                        
                        access_energy[data_type_t::WEIGHT] += (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                              chips[i]->tile_size[data_type_t::WEIGHT]
                                                             *u_read_energy[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8); // Metadata

                        m_scheduler->read_tile_granular_chip_weight[i%m_scheduler->weight_offset_multi_chip.size()] = true;

                        // Update transfer cycle between the chip-level processor and the global buffer
                        transfer_cycle[data_type_t::WEIGHT] += nop_cycle
                                                              *ceil((double)((chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *8*sizeof(data_t))
                                                              /(double)bitwidth) +  // Non-zero data
                                                               nop_cycle 
                                                              *ceil((double)(chips[i]->tile_size[data_type_t::WEIGHT])/(double)bitwidth);
                    }
                   
                    // Update the global buffer access cost
                    chips[i]->access_cycle[data_type_t::WEIGHT] += (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                  *chips[i]->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                   chips[i]->tile_size[data_type_t::WEIGHT]
                                                                  *chips[i]->u_write_cycle[data_type_t::WEIGHT]
                                                                  /(sizeof(data_t)*8); // Metadata
                   
                    chips[i]->access_energy[data_type_t::WEIGHT] += (chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                   *chips[i]->u_write_energy[data_type_t::WEIGHT] + // Non-zero data
                                                                    chips[i]->tile_size[data_type_t::WEIGHT]
                                                                   *chips[i]->u_write_energy[data_type_t::WEIGHT]
                                                                   /(sizeof(data_t)*8); // Metadata

                    // Update transfer energy between the chip-level processor and the global buffer
                    transfer_energy[data_type_t::WEIGHT] += nop_energy
                                                           *ceil((double)((chips[i]->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                           *8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                            nop_energy
                                                           *ceil((double)(chips[i]->tile_size[data_type_t::WEIGHT])/(double)bitwidth); // Metadata

                    chips[i]->skip_transfer[data_type_t::WEIGHT] = false;
                }
            }
        }
        else {
            std::cerr << "Undefined data format" << std::endl;
            exit(1);
        }
#elif defined(NPUSIM_LEGACY_ADDRESS_TIMING)
        if(!skip_transfer[data_type_t::WEIGHT]) {
            num_data_transfer[data_type_t::WEIGHT]++;

            std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
            parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);

            unsigned num_access_global_buffer = 0, num_access_multi_chip = 0;
            for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
                uint64_t address_global_buffer = 0, address_multi_chip = 0;

                for(unsigned k = 0; k < parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned c = 0; c < parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned r = 0; r < parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]; r++) {
                            for(unsigned s = 0; s < parameters_global_buffer[parameter_type_t::FILTER_WIDTH]; s++) {

                                if(address_global_buffer != ((uint64_t)&chips[i]->data[offsets[data_type_t::WEIGHT] + 
                                                                                       k*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                                        *parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                                        *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                                       c*parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                                        *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                                       r*parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                       chips[i]->mask_bits[data_type_t::WEIGHT]) << chips[i]->mask_bits[data_type_t::WEIGHT]) {

                                    // Update costs of global buffer
                                    chips[i]->access_energy[data_type_t::WEIGHT] += chips[i]->u_write_energy[data_type_t::WEIGHT];
                                    if(!chips[i]->double_buffer) chips[i]->access_cycle[data_type_t::WEIGHT] += chips[i]->u_write_cycle[data_type_t::WEIGHT];
                                    num_access_global_buffer++;

                                    address_global_buffer = ((uint64_t)&chips[i]->data[offsets[data_type_t::INPUT] + 
                                                                                       k*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                                        *parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                                        *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                                       c*parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                                        *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                                       r*parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                       chips[i]->mask_bits[data_type_t::WEIGHT]) << chips[i]->mask_bits[data_type_t::WEIGHT];
                                }
                                if(!m_scheduler->read_tile_granular_chip_weight[i%m_scheduler->weight_offset_multi_chip.size()] &&
                                   address_multi_chip != ((uint64_t)&data[offsets[data_type_t::WEIGHT] + 
                                                                          m_scheduler->weight_offset_multi_chip[i%m_scheduler->weight_offset_multi_chip.size()] + 
                                                                          k*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                           *parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                           *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                          c*parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                           *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                          r*parameters_multi_chip[parameter_type_t::FILTER_HEIGHT] + s] >> 
                                                                          mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT]) {

                                    // Update cost of temporal buffer in chip-level processors
                                    access_energy[data_type_t::WEIGHT] += u_read_energy[data_type_t::WEIGHT];
                                    if(!chips[i]->double_buffer) access_cycle[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT];
                                    num_access_multi_chip++;

                                    // Update address of temporal buffer in chip-level processors
                                    m_scheduler->read_tile_granular_chip_weight[i%m_scheduler->weight_offset_multi_chip.size()] = true;
                                    address_multi_chip = ((uint64_t)&data[offsets[data_type_t::WEIGHT] + 
                                                                          m_scheduler->weight_offset_multi_chip[i%m_scheduler->weight_offset_multi_chip.size()] + 
                                                                          k*parameters_multi_chip[parameter_type_t::INPUT_CHANNEL]
                                                                           *parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                           *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                          c*parameters_multi_chip[parameter_type_t::FILTER_HEIGHT]
                                                                           *parameters_multi_chip[parameter_type_t::FILTER_WIDTH] + 
                                                                          r*parameters_multi_chip[parameter_type_t::FILTER_HEIGHT] + s] >> 
                                                                          mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT];

                                }
                            }
                        }
                    }
                }
                chips[i]->skip_transfer[data_type_t::WEIGHT] = false;
            }

            transfer_cycle[data_type_t::WEIGHT] += num_access_multi_chip*nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::WEIGHT])/(double)(bitwidth));
            transfer_energy[data_type_t::WEIGHT] += num_access_multi_chip*nop_energy*ceil((double)(chips[0]->line_size[data_type_t::WEIGHT])/(double)(bitwidth));

            double first_stage = u_read_cycle[data_type_t::WEIGHT];
            double second_stage = std::max(u_read_cycle[data_type_t::WEIGHT], 
                                           nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::WEIGHT])/(double)(bitwidth)));

            double last_before_stage = std::max(chips[get_number_of_active_chips()-1]->u_write_cycle[data_type_t::WEIGHT], 
                                                nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
            double last_stage = chips[get_number_of_active_chips()-1]->u_write_cycle[data_type_t::WEIGHT];

            double other_stage = std::max(u_read_cycle[data_type_t::WEIGHT], 
                                 std::max(chips[0]->u_write_cycle[data_type_t::WEIGHT], 
                                          nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::WEIGHT])/(double)(bitwidth))));
            
            if(exist_temporal_buffer) {
                if(num_access_multi_chip == 0) { } else if(num_access_multi_chip == 1) {
                    cycle_temporal_chips[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT] + 
                                                                 nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::WEIGHT])/(double)(bitwidth)) + 
                                                                 chips[0]->u_write_cycle[data_type_t::WEIGHT];
                } else {
                    cycle_temporal_chips[data_type_t::WEIGHT] += first_stage + second_stage + 
                                                                 (num_access_multi_chip - 2)*other_stage + 
                                                                 last_before_stage + last_stage;
                }
            }
        }
#endif
        for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
            chips[i]->exist_data[data_type_t::WEIGHT] = true, chips[i]->request_to_multi_chip[data_type_t::WEIGHT] = false;
        }
        is_waiting_data[data_type_t::WEIGHT] = false;
    }

    bool request_to_multi_chip_output = false;
    for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
        if(chips[i]->request_to_multi_chip[data_type_t::OUTPUT]) {
            request_to_multi_chip_output = true;
            break;
        }
    }

    if(request_to_multi_chip_output) {
        if(!skip_transfer[data_type_t::OUTPUT]) {
#ifdef FUNCTIONAL
            for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
                // Transfer output data from temporal output buffer of PE array to local buffer of PE.
                m_scheduler->transfer_data(chips[i]->data, data, 
                                           chips[i]->offsets[data_type_t::OUTPUT], offsets[data_type_t::OUTPUT] + m_scheduler->output_offset_multi_chip[i%m_scheduler->output_offset_multi_chip.size()],
                                           component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                                           data_type_t::OUTPUT, chips[i]->get_stationary_type(), action_type_t::LOAD);
                // Update for NPUsim ver2
                //m_scheduler->transfer_data_network_on_chip(chips[i]->data, data, 
                //                                           component_type_t::GLOBAL_BUFFER, component_type_t::CHIPS_Y, 
                //                                           data_type_t::OUTPUT, chips[i]->get_stationary_type(), 
                //                                           action_type_t::LOAD, false);
            }
#endif

#ifndef FUNCTIONAL
            if(dram->transfer_output) {
                account_descriptor_dense_distribution(data_type_t::OUTPUT);
            }
#else
            if(dram->transfer_output) {
                num_data_transfer[data_type_t::OUTPUT]++;

                std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_multi_chip(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
                parameters_multi_chip = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::CHIPS_Y);

                unsigned num_access_global_buffer = 0, num_access_multi_chip = 0;
                for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
                    uint64_t address_global_buffer = 0, address_multi_chip = 0;

                    for(unsigned b = 0; b < parameters_global_buffer[parameter_type_t::BATCH_SIZE]; b++) {
                        for(unsigned k = 0; k < parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                            for(unsigned p = 0; p < parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                                for(unsigned q = 0; q < parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH]; q++) {

                                    if(address_global_buffer != ((uint64_t)&chips[i]->data[offsets[data_type_t::OUTPUT] + 
                                                                                           b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                                            *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                            *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                           k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                            *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                           p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                           chips[i]->mask_bits[data_type_t::OUTPUT]) << chips[i]->mask_bits[data_type_t::OUTPUT]) {

                                        // Update global buffer cost
                                        chips[i]->access_energy[data_type_t::OUTPUT] += chips[i]->u_write_energy[data_type_t::OUTPUT];
                                        if(!chips[i]->double_buffer) chips[i]->access_cycle[data_type_t::OUTPUT] += chips[i]->u_write_cycle[data_type_t::OUTPUT];
                                        num_access_global_buffer++;

                                        // Update global buffer address
                                        address_global_buffer = ((uint64_t)&chips[i]->data[offsets[data_type_t::OUTPUT] + 
                                                                                           b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                                            *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                            *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                           k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                                            *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                           p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                           chips[i]->mask_bits[data_type_t::OUTPUT]) << chips[i]->mask_bits[data_type_t::OUTPUT];
                                    }
                                    if(!m_scheduler->read_tile_granular_chip_output[i%m_scheduler->output_offset_multi_chip.size()] &&
                                       address_multi_chip != ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                              m_scheduler->output_offset_multi_chip[i%m_scheduler->output_offset_multi_chip.size()] + 
                                                                              b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                               *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                               *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                              k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                               *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                              p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                              mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT]) {

                                        // Update cost of temporal buffer in on-chip processors
                                        access_energy[data_type_t::OUTPUT] += u_read_energy[data_type_t::OUTPUT];
                                        if(!chips[i]->double_buffer) access_cycle[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT];
                                        num_access_multi_chip++;

                                        //Update address of temporal buffer in on-chip processors
                                        m_scheduler->read_tile_granular_chip_output[i%m_scheduler->output_offset_multi_chip.size()] = true;
                                        address_multi_chip = ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                              m_scheduler->output_offset_multi_chip[i%m_scheduler->output_offset_multi_chip.size()] + 
                                                                              b*parameters_multi_chip[parameter_type_t::OUTPUT_CHANNEL]
                                                                               *parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                               *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                              k*parameters_multi_chip[parameter_type_t::OUTPUT_HEIGHT]
                                                                               *parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + 
                                                                              p*parameters_multi_chip[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                              mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT];
                                    }
                                }
                            }
                        }
                    }
                    chips[i]->skip_transfer[data_type_t::OUTPUT] = false;
                }
   
                // Update transfer cycle and energy
                transfer_cycle[data_type_t::OUTPUT] += num_access_multi_chip*nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::OUTPUT])/(double)(bitwidth));
                transfer_energy[data_type_t::OUTPUT] += num_access_multi_chip*nop_energy*ceil((double)(chips[0]->line_size[data_type_t::OUTPUT])/(double)(bitwidth));

                // Overlapped cycle at 1, 2, before last, and last stages
                double first_stage = u_read_cycle[data_type_t::OUTPUT];
                double second_stage = std::max(u_read_cycle[data_type_t::OUTPUT], 
                                               nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::OUTPUT])/(double)(bitwidth)));

                double last_before_stage = std::max(chips[get_number_of_active_chips()-1]->u_write_cycle[data_type_t::OUTPUT], 
                                                    nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::OUTPUT])/(double)(bitwidth)));
                double last_stage = chips[get_number_of_active_chips()-1]->u_write_cycle[data_type_t::OUTPUT];

                // Overlapped cycle at remainder stages
                double other_stage = std::max(u_read_cycle[data_type_t::OUTPUT], 
                                     std::max(chips[0]->u_write_cycle[data_type_t::OUTPUT],
                                              nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::OUTPUT])/(double)(bitwidth))));

                if(exist_temporal_buffer) {
                    if(num_access_multi_chip == 0) { } else if(num_access_multi_chip == 1) {
                        cycle_temporal_chips[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT] +
                                                                     nop_cycle*ceil((double)(chips[0]->line_size[data_type_t::OUTPUT])/(double)(bitwidth)) + 
                                                                     chips[0]->u_write_cycle[data_type_t::OUTPUT];
                    } else {
                        cycle_temporal_chips[data_type_t::OUTPUT] += first_stage + second_stage + 
                                                                     (num_access_multi_chip - 2)*other_stage + 
                                                                     last_before_stage + last_stage;
                    }
                }
            }
#endif
        }
        for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
            chips[i]->exist_data[data_type_t::OUTPUT] = true, chips[i]->request_to_multi_chip[data_type_t::OUTPUT] = false;
        }
        is_waiting_data[data_type_t::OUTPUT] = false;
    }
    for(unsigned i = 0; i < get_number_of_active_chips(); i++) {
        chips[i]->fill_data();
    }

    if(tile_size[data_type_t::INPUT] == dram->tile_size[data_type_t::INPUT]) {skip_transfer[data_type_t::INPUT] = true;}
    if(tile_size[data_type_t::WEIGHT] == dram->tile_size[data_type_t::WEIGHT]) {skip_transfer[data_type_t::WEIGHT] = true;}
    if(tile_size[data_type_t::OUTPUT] == dram->tile_size[data_type_t::OUTPUT]) {skip_transfer[data_type_t::OUTPUT] = true;}

#ifdef PRINT
    std::cout << "Transfer data from Multi Chip to Global buffer" << std::endl;
#endif
}

// Flush data at temporal buffer in PE array.
void multi_chip_t::flush_data() {
    // At least one PE sends an input request to PE array.
    for(unsigned i = 0; i < num_active_chips_x*num_active_chips_y; i++) {
        if(chips[i]->request_to_multi_chip[data_type_t::INPUT]) {
            exist_data[data_type_t::INPUT] = false;
            break;
        }
    }
    // At least one PE sends a weight request to PE array.
    for(unsigned i = 0; i < num_active_chips_x*num_active_chips_y; i++) {
        if(chips[i]->request_to_multi_chip[data_type_t::WEIGHT]) {
            exist_data[data_type_t::WEIGHT] = false;
            break;
        }
    }
    // At least one PE sends an output request to PE array.
    for(unsigned i = 0; i < num_active_chips_x*num_active_chips_y; i++) {
        if(chips[i]->request_to_multi_chip[data_type_t::OUTPUT]) {
            exist_data[data_type_t::OUTPUT] = false;
            break;
        }
    }
}


void multi_chip_t::print_specification() {
    std::cout << "================ Multi Chip ================" << std::endl;

    std::cout << "Type               :" << std::setw(24) 
                                        << "Mesh" << std::endl;
    std::cout << "Height             :" << std::setw(24) 
                                        << height << std::endl;
    std::cout << "Width              :" << std::setw(24) 
                                        << width << std::endl;

    if(stationary_type == stationary_type_t::INPUT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24) 
                                            << "Input stationary" << std::endl;
    }
    else if(stationary_type == stationary_type_t::WEIGHT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24)
                                            << "Weight stationary" << std::endl;
    }
    else if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24) 
                                            << "Output stationary" << std::endl;
    }
    std::cout << "Bandwidth          :" << std::setw(19) << std::setprecision(0) 
                                        << bandwidth << " GB/s" << std::endl;

    std::cout << "NoP cycle          :" << std::setw(17) 
                                        << nop_cycle << " cycles" << std::endl;
    std::cout << "NoP energy         :" << std::setw(21)
                                        << nop_energy << " pJ" << std::endl;
    std::cout << std::endl;
}

// Accumulate leakage energy of the temporal buffer once for the modeled layer duration.
// static_energy in config is leakage in pJ/cycle; the temporal buffer is always-on so it
// leaks for the whole layer elapsed window regardless of how many chips are active.
void multi_chip_t::update_static_energy(double elapsed_cycles) {
    if(elapsed_cycles < 0.0) {
        std::cerr << "Error: multi_chip elapsed cycles must be non-negative" << std::endl;
        exit(1);
    }
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        if(u_static_energy[i] < 0.0) {
            std::cerr << "Error: multi_chip static_energy must be a non-negative pJ/cycle value" << std::endl;
            exit(1);
        }
        static_energy[i] = static_energy_for_cycles(u_static_energy[i], elapsed_cycles);
    }
}

// Modeled busy duration of the multi-chip fabric for the current layer: the max of its
// access, NoP transfer, and overlapped cost axes.
double multi_chip_t::modeled_elapsed_cycles() const {
    double elapsed = std::max(write_back_cycle, overlapped_transfer_cycle);
    for(unsigned type = 0; type < data_type_t::NUM_DATA_TYPES; ++type) {
        elapsed = std::max(elapsed, access_cycle[type]);
        elapsed = std::max(elapsed, transfer_cycle[type]);
        elapsed = std::max(elapsed, cycle_temporal_chips[type]);
    }
    return elapsed;
}

void multi_chip_t::reset() {
    std::fill_n(data, (memory_size + sizeof(data_t) - 1)/sizeof(data_t), data_t{});

    initial = true;
    equal_output_tile = false;

    // Reset utilization of chip-level processor
    utilization = 0.0;
    write_back_cycle = 0.0, overlapped_transfer_cycle = 0.0;

    num_request.assign(data_type_t::NUM_DATA_TYPES, 0);
    num_data_transfer.assign(data_type_t::NUM_DATA_TYPES, 0);

    access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    payload_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    metadata_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    storage_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);

    static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
}

