#include <iomanip>
#include <string>
#include <cassert>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <limits>
#include "pe.h"
#include "datatype.h"
#include "interconnect_timing.h"
#include "pe_lane.h"
#include "energy_units.h"

namespace {

size_t checked_product(size_t lhs, size_t rhs, const char *context) {
    if(lhs != 0 && rhs > std::numeric_limits<size_t>::max()/lhs) {
        std::cerr << "Error: size overflow while calculating " << context << std::endl;
        exit(1);
    }
    return lhs * rhs;
}

bool is_power_of_two(unsigned value) {
    return value != 0 && (value & (value - 1)) == 0;
}

struct mac_geometry_t {
    unsigned batch;
    unsigned groups;
    unsigned output_channels_per_group;
    unsigned input_channels_per_group;
    unsigned output_height;
    unsigned output_width;
    unsigned input_height;
    unsigned input_width;
    unsigned filter_height;
    unsigned filter_width;
    unsigned stride;
    size_t input_elements;
    size_t weight_elements;
    size_t output_elements;
    size_t operations;
};

mac_geometry_t get_mac_geometry(scheduler_t *m_scheduler) {
    const std::vector<unsigned> parameters =
        m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
    const unsigned groups = parameters[parameter_type_t::GROUP];
    const unsigned output_channels = parameters[parameter_type_t::OUTPUT_CHANNEL];
    const unsigned input_channels = parameters[parameter_type_t::INPUT_CHANNEL];

    if(groups == 0 || output_channels == 0 || input_channels == 0 ||
       output_channels % groups != 0 || input_channels % groups != 0 ||
       parameters[parameter_type_t::BATCH_SIZE] == 0 ||
       parameters[parameter_type_t::OUTPUT_HEIGHT] == 0 ||
       parameters[parameter_type_t::OUTPUT_WIDTH] == 0 ||
       parameters[parameter_type_t::INPUT_HEIGHT] == 0 ||
       parameters[parameter_type_t::INPUT_WIDTH] == 0 ||
       parameters[parameter_type_t::FILTER_HEIGHT] == 0 ||
       parameters[parameter_type_t::FILTER_WIDTH] == 0 ||
       parameters[parameter_type_t::STRIDE] == 0) {
        std::cerr << "Error: invalid MAC mapping geometry" << std::endl;
        exit(1);
    }

    mac_geometry_t geometry;
    geometry.batch = parameters[parameter_type_t::BATCH_SIZE];
    geometry.groups = groups;
    geometry.output_channels_per_group = output_channels/groups;
    geometry.input_channels_per_group = input_channels/groups;
    geometry.output_height = parameters[parameter_type_t::OUTPUT_HEIGHT];
    geometry.output_width = parameters[parameter_type_t::OUTPUT_WIDTH];
    geometry.input_height = parameters[parameter_type_t::INPUT_HEIGHT];
    geometry.input_width = parameters[parameter_type_t::INPUT_WIDTH];
    geometry.filter_height = parameters[parameter_type_t::FILTER_HEIGHT];
    geometry.filter_width = parameters[parameter_type_t::FILTER_WIDTH];
    geometry.stride = parameters[parameter_type_t::STRIDE];

    geometry.input_elements = checked_product(
        checked_product(checked_product(geometry.batch, input_channels, "MAC input elements"),
                        geometry.input_height, "MAC input elements"),
        geometry.input_width, "MAC input elements");
    geometry.weight_elements = checked_product(
        checked_product(checked_product(output_channels, geometry.input_channels_per_group, "MAC weight elements"),
                        geometry.filter_height, "MAC weight elements"),
        geometry.filter_width, "MAC weight elements");
    geometry.output_elements = checked_product(
        checked_product(checked_product(geometry.batch, output_channels, "MAC output elements"),
                        geometry.output_height, "MAC output elements"),
        geometry.output_width, "MAC output elements");
    geometry.operations = checked_product(
        checked_product(geometry.output_elements, geometry.input_channels_per_group, "MAC operations"),
        checked_product(geometry.filter_height, geometry.filter_width, "MAC operations"), "MAC operations");
    return geometry;
}

void validate_mac_geometry(const mac_geometry_t &geometry, size_t register_capacity,
                           unsigned active_macs) {
    if(geometry.input_elements > register_capacity ||
       geometry.weight_elements > register_capacity ||
       geometry.output_elements > register_capacity) {
        std::cerr << "Error: MAC tile exceeds PE register capacity" << std::endl;
        exit(1);
    }
    if(geometry.operations != active_macs) {
        std::cerr << "Error: MAC mapping operation count does not match active MAC count" << std::endl;
        exit(1);
    }
}

} // namespace

pe_t::pe_t(section_config_t m_section_config) :
    input_data_mac(NULL),
    weight_mac(NULL),
    output_data_mac(NULL),
    input_data_lb(NULL),
    weight_lb(NULL),
    output_data_lb(NULL),
    input_size(0),
    weight_size(0),
    output_size(0),
    edge_accumulation(false),
    double_buffer(true),
    index(0),
    duplicated_input_mac(0),
    duplicated_input_lb(0),
    /* Unit stats */
    u_computation_cycle(0.0),
    u_computation_energy(0.0),
    accumulator_reload_bytes(0),
    accumulator_spill_bytes(0),
    accumulator_create_events(0),
    accumulator_retained_events(0),
    accumulator_energy(0.0),
    reduction_energy(0.0),
    compute_energy_basis("computation_energy"),
    compute_energy_precision_calibrated(false),
    u_mac_reduction_energy(0.0),
    u_transfer_cycle(0.0),
    u_transfer_energy(0.0),
    /* Accelerator stats */
    num_computation(0),
    computation_cycle(0.0),
    computation_energy(0.0),
    utilization_mac(0.0),
    write_back_cycle_mac(0.0),
    write_back_cycle_lb(0.0),
    overlapped_transfer_cycle(0.0),
    pe_array(NULL),
    stationary_type_mac(stationary_type_t::UNDEFINED_STATIONARY),
    stationary_type_local_buffer(stationary_type_t::UNDEFINED_STATIONARY),
    parameter_order("kbpqcrs"),
    memory_type(memory_type_t::SEPARATE),
    mac_type(mac_type_t::UNDEFINED_MAC),
    num_macs(1),
    mac_width(1),
    num_active_macs(1),
    active_mac_width(1),
    lane_state(),
    active_mac_units(1),
    mac_register_capacity(1),
    frequency(0.0),
    bandwidth(0.0),
    bitwidth(0),
    input_index(0),
    weight_index(0),
    output_index(0),
    input_flush_counter(0),
    weight_flush_counter(0),
    output_flush_counter(0),
    idle(false) {

    init(m_section_config);

}

pe_t::~pe_t() {

    // Free memory space at MAC unit.
    delete [] input_data_mac;
    delete [] weight_mac;
    delete [] output_data_mac;

    // Free memory space at local buffer.
    delete [] input_data_lb;
    delete [] weight_lb;
    delete [] output_data_lb;

}

// Initialize the PE.
void pe_t::init(section_config_t m_section_config) {

    /* Initialize PE specifications */

    // Initialize the size of input, weight, and output buffer in Byte.
    const bool have_size_keys =
        m_section_config.get_setting("input_size", &input_size) &&
        m_section_config.get_setting("weight_size", &weight_size) &&
        m_section_config.get_setting("output_size", &output_size);
    if(!have_size_keys &&
       (!m_section_config.get_setting("input_buffer", &input_size) ||
        !m_section_config.get_setting("weight_buffer", &weight_size) ||
        !m_section_config.get_setting("output_buffer", &output_size))) {
        std::cerr << "Error: PE requires input/weight/output buffer sizes" << std::endl;
        exit(1);
    }
    if(input_size == 0 || weight_size == 0 || output_size == 0) {
        std::cerr << "Error: PE buffer sizes must be non-zero" << std::endl;
        exit(1);
    }

    // V3: outputs accumulate at the array edge (Gemmini-style accumulator); the OUTPUT
    // tile then streams through rather than residing in the local buffer.
    m_section_config.get_setting("edge_accumulation", &edge_accumulation);

    // LB7/P1-A: a single-buffered local buffer can't overlap the next tile's load with
    // the current tile's compute, so its LB<->MAC transfer makespan serializes against
    // the compute schedule (see stats_t::finalize_layer_timeline()). Defaults to true
    // (the model's original always-overlap behavior) so existing configs are unaffected
    // unless they opt out; every shipped accelerator config now states it explicitly so
    // the assumption is not silent.
    //
    // CAPACITY CONTRACT (decided, P1-A item 4): this flag is a TIMING property of the
    // LB->MAC-register fill, NOT a second resident local-buffer tile, so check_tile_size()
    // requires ONE live copy in both modes -- unlike global_buffer_t, whose double_buffer
    // genuinely allocates two tile halves. Rationale: what a fill-behind-compute overlap
    // needs twice is the MAC register file, and this model represents that file as exactly
    // the physical lane count (number_of_macs*mac_width) holding ONE tile, with no 2-deep
    // register model to charge; the LB tile itself is read, not held twice. Charging two
    // copies of either buffer would reject every validated config (Eyeriss conv3 already
    // fills 100% of its input LB partition, and every config's MAC tile exactly fills its
    // lane capacity), i.e. declare the validated silicon un-modelable. The active mode is
    // printed in the layer-timeline result section so the contract is visible per run.
    m_section_config.get_setting("double_buffer", &double_buffer);

    // Initialize the number of elements in each buffer.
    size_t num_input = (static_cast<size_t>(input_size) + sizeof(data_t) - 1)/sizeof(data_t);
    size_t num_weight = (static_cast<size_t>(weight_size) + sizeof(data_t) - 1)/sizeof(data_t);
    size_t num_output = (static_cast<size_t>(output_size) + sizeof(data_t) - 1)/sizeof(data_t);

    // Allocate the memory space for local buffer.
    input_data_lb  = new data_t[num_input]();
    weight_lb      = new data_t[num_weight]();
    output_data_lb = new data_t[num_output]();

    // Initialize the size of MAC register.
    m_section_config.get_setting("number_of_macs", &num_macs);
    m_section_config.get_setting("mac_width", &mac_width);
    if(num_macs == 0 || mac_width == 0) {
        std::cerr << "Error: number_of_macs and mac_width must be non-zero" << std::endl;
        exit(1);
    }
    mac_register_capacity = checked_product(num_macs, mac_width, "MAC register capacity");

    input_data_mac  = new data_t[mac_register_capacity]();
    weight_mac      = new data_t[mac_register_capacity]();
    // A mapping may expose up to num_macs * mac_width independent outputs.
    // Keep scalar storage for that maximum even when lanes reduce together.
    output_data_mac = new data_t[mac_register_capacity]();

    // Initialize the frequency (MHz) and bandwidth (GB/sec).
    const bool have_frequency = m_section_config.get_setting("frequency", &frequency);
    const bool have_bandwidth = m_section_config.get_setting("bandwidth", &bandwidth);
    const bool have_bitwidth = m_section_config.get_setting("bitwidth", &bitwidth);
    if(!have_bitwidth) {
        if(!have_frequency || !have_bandwidth || frequency <= 0.0f || bandwidth <= 0.0f) {
            std::cerr << "Error: PE requires a positive bitwidth or positive frequency and bandwidth" << std::endl;
            exit(1);
        }
        bitwidth = static_cast<unsigned>(8.0f*bandwidth/frequency);
    }
    if(bitwidth == 0) {
        std::cerr << "Error: PE bitwidth must be non-zero" << std::endl;
        exit(1);
    }

    bypass.reserve(data_type_t::NUM_DATA_TYPES);
    bypass.assign(data_type_t::NUM_DATA_TYPES, 0);
    m_section_config.get_vector_setting("bypass", &bypass);
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        if(bypass[i]) {
            std::cerr << "Error: PE local-buffer bypass is not implemented" << std::endl;
            exit(1);
        }
    }

    // Initialize line size and mask bits of MAC register and local buffer
    line_size_mac.reserve(data_type_t::NUM_DATA_TYPES);
    line_size_mac.assign(data_type_t::NUM_DATA_TYPES, 8);

    mask_bits_mac.reserve(data_type_t::NUM_DATA_TYPES);
    mask_bits_mac.assign(data_type_t::NUM_DATA_TYPES, 0);

    line_size_lb.reserve(data_type_t::NUM_DATA_TYPES);
    line_size_lb.assign(data_type_t::NUM_DATA_TYPES, 8);

    mask_bits_lb.reserve(data_type_t::NUM_DATA_TYPES);
    mask_bits_lb.assign(data_type_t::NUM_DATA_TYPES, 0);

    m_section_config.get_vector_setting("mac_line_size", &line_size_mac);
    m_section_config.get_vector_setting("lb_line_size", &line_size_lb);

    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        unsigned mac_line_size = line_size_mac[i];
        unsigned lb_line_size = line_size_lb[i];
        if(mac_line_size < 8 || lb_line_size < 8 ||
           !is_power_of_two(mac_line_size) || !is_power_of_two(lb_line_size)) {
            std::cerr << "Error: PE line sizes must be power-of-two values of at least 8 bits" << std::endl;
            exit(1);
        }
        while(mac_line_size > 8) {
            mac_line_size /= 2;
            mask_bits_mac[i]++;
        }

        while(lb_line_size > 8) {
            lb_line_size /= 2;
            mask_bits_lb[i]++;
        }
    }

    // Initialize the stationary type.
    // Initialize the stationary type between MAC and local buffer.
    std::string stationary_str;
    if(!m_section_config.get_setting("mac_stationary", &stationary_str) ||
       !is_valid_type(stationary_type_str, stationary_str) ||
       stationary_str == "undefined_stationary") {
        std::cerr << "Error: PE requires a valid mac_stationary" << std::endl;
        exit(1);
    }
    stationary_type_mac = (stationary_type_t)get_type(stationary_type_str, stationary_str);

    // Initialize the order of parameters.
    if(stationary_type_mac == stationary_type_t::INPUT_STATIONARY) {
        parameter_order = "BCPQKRS";
    } else if(stationary_type_mac == stationary_type_t::WEIGHT_STATIONARY) {
        parameter_order = "KCRSBPQ";
    } else if(stationary_type_mac == stationary_type_t::OUTPUT_STATIONARY) {
        parameter_order = "BKPQCRS";
    }
    m_section_config.get_setting("pe_parameter_order", &parameter_order);

    // Initialize the stationary type between PE array and Global buffer.
    if(!m_section_config.get_setting("pe_stationary", &stationary_str) ||
       !is_valid_type(stationary_type_str, stationary_str) ||
       stationary_str == "undefined_stationary") {
        std::cerr << "Error: PE requires a valid pe_stationary" << std::endl;
        exit(1);
    }
    stationary_type_local_buffer = (stationary_type_t)get_type(stationary_type_str, stationary_str);

    std::string memory_str;
    if(!m_section_config.get_setting("memory_type", &memory_str)) {
        std::cerr << "Error: PE requires memory_type" << std::endl;
        exit(1);
    }
    lowercase(memory_str);
    if(memory_str == "separated") memory_str = "separate";
    if(memory_str != "separate") {
        std::cerr << "Error: shared PE local buffers are not implemented" << std::endl;
        exit(1);
    }
    memory_type = memory_type_t::SEPARATE;

    skip_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);

    // Initialize the tile size of MAC unit.
    tile_size_mac.reserve(data_type_t::NUM_DATA_TYPES);
    tile_size_mac.assign(data_type_t::NUM_DATA_TYPES, 1);

    // Initialize the tile size of Local buffer.
    tile_size_lb.reserve(data_type_t::NUM_DATA_TYPES);
    tile_size_lb.assign(data_type_t::NUM_DATA_TYPES, 1);

    offsets_mac.reserve(data_type_t::NUM_DATA_TYPES);
    offsets_mac.assign(data_type_t::NUM_DATA_TYPES, 0);

    offsets_lb.reserve(data_type_t::NUM_DATA_TYPES);
    offsets_lb.assign(data_type_t::NUM_DATA_TYPES, 0);

    /* Initialize PE signals */

    // Exist data in MAC unit.
    exist_data_mac.reserve(data_type_t::NUM_DATA_TYPES);
    exist_data_mac.assign(data_type_t::NUM_DATA_TYPES, false);

    // Exist data in Local buffer.
    exist_data_lb.reserve(data_type_t::NUM_DATA_TYPES);
    exist_data_lb.assign(data_type_t::NUM_DATA_TYPES, false);

    // Exist request to Local buffer (from MAC unit).
    request_to_lb.reserve(data_type_t::NUM_DATA_TYPES);
    request_to_lb.assign(data_type_t::NUM_DATA_TYPES, true);

    // Exist request to PE array (from Local buffer).
    request_to_pe_array.reserve(data_type_t::NUM_DATA_TYPES);
    request_to_pe_array.assign(data_type_t::NUM_DATA_TYPES, false);

    /* Initialize PE unit costs */

    // Initialize the MAC unit cycle and energy.
    m_section_config.get_setting("computation_cycle", &u_computation_cycle);
    m_section_config.get_setting("computation_energy", &u_computation_energy);
    // E3/RE4: prefer a per-precision MAC energy for the precisions actually in use. The key names
    // all THREE datapath widths that set a MAC's energy:
    //
    //     mac_energy_<input>_<weight>_<accumulator>
    //
    // The first version of this used only the two operand formats, which made an INT8 x INT8 MAC
    // accumulating into FP32 cost exactly the same as one accumulating into FP16 -- the adder and
    // the accumulator register are a large part of a MAC's energy, and that difference was
    // invisible. The two-operand key is still honoured as a PARTIAL calibration: it prices the
    // multiplier for these operands but says nothing about the accumulator, so it does not qualify
    // the config for absolute power (RE2). Only the three-part key does.
    const std::string operand_key = "mac_energy_" +
        runtime_datatypes().format(data_type_t::INPUT).name + "_" +
        runtime_datatypes().format(data_type_t::WEIGHT).name;
    const std::string precision_key =
        operand_key + "_" + runtime_datatypes().accumulator_format().name;
    double precision_energy = 0.0;
    if(m_section_config.get_setting(precision_key, &precision_energy)) {
        u_computation_energy = precision_energy;
        compute_energy_basis = precision_key;
        compute_energy_precision_calibrated = true;
    } else if(m_section_config.get_setting(operand_key, &precision_energy)) {
        // Operand widths priced, accumulator width not. Use it -- it is closer than the scalar --
        // but do NOT call the compute cost calibrated: every accumulator format would share it.
        u_computation_energy = precision_energy;
        compute_energy_basis = operand_key + " (operands only; accumulator " +
                               runtime_datatypes().accumulator_format().name + " UNCALIBRATED)";
        compute_energy_precision_calibrated = false;
    } else {
        // No entry for this precision: keep the scalar, but say so instead of implying the
        // number is precision-aware.
        compute_energy_basis = "computation_energy";
        compute_energy_precision_calibrated = false;
    }
    // P4-2/PE2: unit energy of one lane->accumulator adder-tree addition.
    m_section_config.get_setting("mac_reduction_energy", &u_mac_reduction_energy);

    // Initialize the transfer cycle/energy between PE and MAC unit and local buffer.
    m_section_config.get_setting("transfer_cycle_pe", &u_transfer_cycle);
    m_section_config.get_setting("transfer_energy_pe", &u_transfer_energy);

    // Initialize the unit cycle and energy of MAC unit
    u_read_cycle_mac.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_cycle_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("mac_read_cycle", &u_read_cycle_mac);

    u_read_energy_mac.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_energy_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("mac_read_energy", &u_read_energy_mac);

    u_write_cycle_mac.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_cycle_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("mac_write_cycle", &u_write_cycle_mac);

    u_write_energy_mac.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_energy_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("mac_write_energy", &u_write_energy_mac);

    // Initialize the unit cycle and energy of local buffer
    u_read_cycle_lb.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_cycle_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("lb_read_cycle", &u_read_cycle_lb);

    u_read_energy_lb.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_energy_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("lb_read_energy", &u_read_energy_lb);

    u_write_cycle_lb.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_cycle_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("lb_write_cycle", &u_write_cycle_lb);

    u_write_energy_lb.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_energy_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("lb_write_energy", &u_write_energy_lb);

    u_static_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("static_energy", &u_static_energy);

    // LB4: optional local-buffer leakage, separately sweepable from the MAC-side
    // `static_energy`. Both accrue over the same layer elapsed window.
    u_lb_static_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_lb_static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("lb_static_energy", &u_lb_static_energy);

    // Optional format-IP unit costs. Each value is per packed transaction.
    format_payload_events.assign(data_type_t::NUM_DATA_TYPES, 0);
    format_metadata_events.assign(data_type_t::NUM_DATA_TYPES, 0);
    u_format_payload_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    u_format_payload_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    u_format_metadata_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    u_format_metadata_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    u_accumulator_spill_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    u_accumulator_spill_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("format_payload_cycle", &u_format_payload_cycle);
    m_section_config.get_vector_setting("format_payload_energy", &u_format_payload_energy);
    m_section_config.get_vector_setting("format_metadata_cycle", &u_format_metadata_cycle);
    m_section_config.get_vector_setting("format_metadata_energy", &u_format_metadata_energy);
    m_section_config.get_vector_setting("accumulator_spill_cycle", &u_accumulator_spill_cycle);
    m_section_config.get_vector_setting("accumulator_spill_energy", &u_accumulator_spill_energy);
    // E4: the final output cast/pack is a separate event from an accumulator spill, on a
    // different precision, so it gets its own unit cost instead of reusing the spill's.

    /* Initialize PE stats */

    // Total number of data request to local buffer
    num_request_to_lb.reserve(data_type_t::NUM_DATA_TYPES);
    num_request_to_lb.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Total number of data transfer to MAC unit
    num_data_transfer_to_mac.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer_to_mac.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Total access cycles to MAC unit
    access_cycle_mac.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Total access energies to MAC unit
    access_energy_mac.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Total access cycles to local buffer
    access_cycle_lb.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Total access energies to local buffer
    access_energy_lb.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Total transfer cycle between MAC unit and local buffer
    transfer_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Total transfer energies between MAC unit and local buffer
    transfer_energy.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    payload_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    metadata_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    storage_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Overlapped cycle between MAC units and local buffer
    cycle_mac_lb.reserve(data_type_t::NUM_DATA_TYPES);
    cycle_mac_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Total static energy of PE.
    static_energy.reserve(data_type_t::NUM_DATA_TYPES);
    static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    format_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    format_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    utilization_local_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    utilization_local_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

}

// Connect PE to PE array.
void pe_t::connect(pe_array_t *m_pe_array) {
    pe_array = m_pe_array;
}

// Update tile size of local buffer and MAC unit.
void pe_t::update_tile_size(scheduler_t *m_scheduler) {

    num_active_macs = m_scheduler->num_active_mac;
    const mac_lane_state_t lane_state = calculate_mac_lane_state(num_macs, mac_width, num_active_macs);
    if(!lane_state.valid || lane_state.scalar_capacity != mac_register_capacity) {
        std::cerr << "The number of Active macs : " << num_active_macs
                  << " is bigger than the scalar FMA lane capacity : " << mac_register_capacity << std::endl;
        exit(1);
    }
    active_mac_units = lane_state.active_accumulator_units;
    active_mac_width = lane_state.final_accumulator_lanes;
    // L9: keep the whole resolved state; the reduction cost is charged against it.
    this->lane_state = lane_state;
    // PE3: peak across the layer rather than last-call overwrite (identical today;
    // hardening for per-tile variation).
    utilization_mac = std::max(utilization_mac, lane_state.utilization);
    // Update tile sizes.
    tile_size_mac = m_scheduler->tile_size[component_type_t::MAC];
    tile_size_lb = m_scheduler->tile_size[component_type_t::PE];
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        if(tile_size_mac[i] == 0 || tile_size_mac[i] > mac_register_capacity) {
            std::cerr << "Error: MAC tile exceeds PE register capacity" << std::endl;
            exit(1);
        }
    }

}

void pe_t::update_offset() {
    // Update the offset of memory.
    if(memory_type == memory_type_t::SEPARATE) {
        offsets_lb[data_type_t::INPUT] = 0, offsets_lb[data_type_t::WEIGHT] = 0, offsets_lb[data_type_t::OUTPUT] = 0;
    }
    else if(memory_type == memory_type_t::SHARED) {
        offsets_lb[data_type_t::INPUT] = 0;
        offsets_lb[data_type_t::WEIGHT] = tile_size_lb[data_type_t::INPUT];
        offsets_lb[data_type_t::OUTPUT] = tile_size_lb[data_type_t::INPUT] + tile_size_lb[data_type_t::WEIGHT];
    }
}

void pe_t::check_tile_size() {
    // Check the validity.
    if(memory_type == memory_type_t::SEPARATE) {
        if(runtime_datatypes().storage_bytes(data_type_t::INPUT, tile_size_lb[data_type_t::INPUT]) > input_size ||
           runtime_datatypes().storage_bytes(data_type_t::WEIGHT, tile_size_lb[data_type_t::WEIGHT]) > weight_size ||
           (!edge_accumulation &&
            runtime_datatypes().storage_bytes(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT]) > output_size)) {
            std::cerr << "The data size is bigger than local buffer size\n"
                      << "Input data : " << runtime_datatypes().storage_bytes(data_type_t::INPUT, tile_size_lb[data_type_t::INPUT])
                      << " Input buffer : " << input_size << "\t"
                      << "Weight : " << runtime_datatypes().storage_bytes(data_type_t::WEIGHT, tile_size_lb[data_type_t::WEIGHT])
                      << " Weight buffer : " << weight_size << "\t"
                      << "Output data : " << runtime_datatypes().storage_bytes(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT])
                      << " Output buffer : " << output_size << std::endl;
            exit(1);
        }
    }
    else if(memory_type == memory_type_t::SHARED) {
        size_t data_size = 0;
        for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
            const data_type_t type = static_cast<data_type_t>(i);
            // V3: edge-accumulated outputs do not occupy the local buffer.
            if(edge_accumulation && type == data_type_t::OUTPUT) continue;
            if(!bypass[type]) {
                data_size += runtime_datatypes().storage_bytes(type, tile_size_lb[type]);
            }
        }

        const size_t buffer_size = static_cast<size_t>(input_size) + weight_size + output_size;
        if(data_size > buffer_size) {
            std::cerr << "The data size is bigger than local buffer size\n"
                      << "Data : " << data_size
                      << " Buffer : " << buffer_size << std::endl;
            exit(1);
        }
    }
}

// Get stationary_type of MAC register
stationary_type_t pe_t::get_mac_stationary_type() {
    return stationary_type_mac;
}

stationary_type_t pe_t::get_local_buffer_stationary_type() {
    return stationary_type_local_buffer;
}

// Return parameter order.
std::string pe_t::get_parameter_order() {
    return parameter_order;
}

memory_type_t pe_t::get_memory_type() {
    return memory_type;
}

// A signal to check whether the PE is idle state or not.
bool pe_t::is_idle() {
    return idle;
}

// A signal to check whether all data exist in the local buffer or not.
bool pe_t::is_exist_data() {
    if(exist_data_lb[data_type_t::INPUT] && exist_data_lb[data_type_t::WEIGHT] && exist_data_lb[data_type_t::OUTPUT]) {
        return true;
    }
    else {
        return false;
    }
}

// A signal to check whether at least one request exists from local buffer to PE array or not.
bool pe_t::is_exist_request() {
    if(request_to_pe_array[data_type_t::INPUT] || request_to_pe_array[data_type_t::WEIGHT] || request_to_pe_array[data_type_t::OUTPUT]) {
        return true;
    }
    else {
        return false;
    }
}


// Wait for the data comes from Global buffer.
void pe_t::wait_data() {
    idle = true;
}

// A signal that the data exist in the PE.
void pe_t::fill_data() {
    idle = false;
}

// Request data to PE array.
void pe_t::request_data() {
    // If input data does not exist in Local buffer.
    if(!exist_data_lb[data_type_t::INPUT]) {
        // Request input data to PE array.
        request_to_pe_array[data_type_t::INPUT] = true;
        pe_array->exist_data[data_type_t::INPUT] = false;

        // Update Stats.
        pe_array->num_request[data_type_t::INPUT]++;
    }
    // If weight does not exist in Local buffer.
    if(!exist_data_lb[data_type_t::WEIGHT]) {
        // Request weight to PE array.
        request_to_pe_array[data_type_t::WEIGHT] = true;
        pe_array->exist_data[data_type_t::WEIGHT] = false;

        // Update Stats.
        pe_array->num_request[data_type_t::WEIGHT]++;
    }
    // If output data does not exist in Local buffer.
    if(!exist_data_lb[data_type_t::OUTPUT]) {
        // Request output data to PE array.
        request_to_pe_array[data_type_t::OUTPUT] = true;
        pe_array->exist_data[data_type_t::OUTPUT] = false;

        // Update stats.
        pe_array->num_request[data_type_t::OUTPUT]++;
    }
}

void pe_t::account_descriptor_dense_mac_transfer(data_type_t type, size_t elements, bool to_mac) {
    const unsigned source_line = to_mac ? line_size_lb[type] : line_size_mac[type];
    const unsigned destination_line = to_mac ? line_size_mac[type] : line_size_lb[type];
    const datatype_transfer_timing_t timing = datatype_transfer_timing(
        type, elements, source_line, destination_line, bitwidth);

    payload_link_transactions[type] += timing.payload_link_transactions;
    metadata_link_transactions[type] += timing.metadata_link_transactions;
    storage_link_transactions[type] += timing.link_transactions;

    const double source_cycle = to_mac ? u_read_cycle_lb[type] : u_read_cycle_mac[type];
    const double source_energy = to_mac ? u_read_energy_lb[type] : u_read_energy_mac[type];
    const double destination_cycle = to_mac ? u_write_cycle_mac[type] : u_write_cycle_lb[type];
    const double destination_energy = to_mac ? u_write_energy_mac[type] : u_write_energy_lb[type];
    if(to_mac) {
        access_cycle_lb[type] += timing.source_accesses*source_cycle;
        access_energy_lb[type] += timing.source_accesses*source_energy;
        access_cycle_mac[type] += timing.destination_accesses*destination_cycle;
        access_energy_mac[type] += timing.destination_accesses*destination_energy;
    } else {
        access_cycle_mac[type] += timing.source_accesses*source_cycle;
        access_energy_mac[type] += timing.source_accesses*source_energy;
        access_cycle_lb[type] += timing.destination_accesses*destination_cycle;
        access_energy_lb[type] += timing.destination_accesses*destination_energy;
        write_back_cycle_mac += timing.source_accesses*source_cycle;
        write_back_cycle_lb += timing.destination_accesses*destination_cycle;
        overlapped_transfer_cycle += timing.link_transactions*u_transfer_cycle;
    }
    transfer_cycle[type] += timing.link_transactions*u_transfer_cycle;
    transfer_energy[type] += timing.link_transactions*u_transfer_energy;
    if(timing.pipeline_transactions > std::numeric_limits<unsigned>::max()) {
        std::cerr << "Error: PE pipeline transaction count overflow" << std::endl;
        exit(1);
    }
    // L2: the packet decomposition carries where the LB/MAC line widths pack and unpack.
    cycle_mac_lb[type] += pipelined_transfer_cycles(
        timing.groups, source_cycle, u_transfer_cycle, destination_cycle);
}

// Transfer data from local buffer to MAC unit.
// And execute MAC operation.
void pe_t::data_transfer_to_mac(scheduler_t *m_scheduler) {


    // LB2: track the peak occupancy across the layer instead of overwriting per call
    // (identical today since the LB tile is fixed per layer; hardens against future
    // per-tile variation), and keep double precision.
    if(!bypass[data_type_t::INPUT]) {
        utilization_local_buffer[data_type_t::INPUT] = std::max(utilization_local_buffer[data_type_t::INPUT],
            static_cast<double>(runtime_datatypes().storage_bytes(data_type_t::INPUT, tile_size_lb[data_type_t::INPUT]))/static_cast<double>(input_size));
    }
    if(!bypass[data_type_t::WEIGHT]) {
        utilization_local_buffer[data_type_t::WEIGHT] = std::max(utilization_local_buffer[data_type_t::WEIGHT],
            static_cast<double>(runtime_datatypes().storage_bytes(data_type_t::WEIGHT, tile_size_lb[data_type_t::WEIGHT]))/static_cast<double>(weight_size));
    }
    if(!bypass[data_type_t::OUTPUT]) {
        utilization_local_buffer[data_type_t::OUTPUT] = std::max(utilization_local_buffer[data_type_t::OUTPUT],
            static_cast<double>(runtime_datatypes().storage_bytes(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT]))/static_cast<double>(output_size));
    }
#ifndef FUNCTIONAL
    if(m_scheduler->compression_type != compression_type_t::DENSE) {
        std::cerr << "Error: timing PE supports dense descriptor traffic only" << std::endl;
        exit(1);
    }
    const bool last_pe = pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 &&
                         index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1;
    if(request_to_lb[data_type_t::INPUT]) {
        account_format_events(data_type_t::INPUT, tile_size_mac[data_type_t::INPUT]);
        if(!skip_transfer[data_type_t::INPUT]) {
            num_data_transfer_to_mac[data_type_t::INPUT]++;
            account_descriptor_dense_mac_transfer(data_type_t::INPUT, tile_size_mac[data_type_t::INPUT], true);
            if(last_pe) move_front(&m_scheduler->input_offset_pe);
        }
        input_index++;
        exist_data_mac[data_type_t::INPUT] = true;
        request_to_lb[data_type_t::INPUT] = false;
        if(tile_size_mac[data_type_t::INPUT] == tile_size_lb[data_type_t::INPUT]) skip_transfer[data_type_t::INPUT] = true;
    }
    if(request_to_lb[data_type_t::WEIGHT]) {
        account_format_events(data_type_t::WEIGHT, tile_size_mac[data_type_t::WEIGHT]);
        if(!skip_transfer[data_type_t::WEIGHT]) {
            num_data_transfer_to_mac[data_type_t::WEIGHT]++;
            account_descriptor_dense_mac_transfer(data_type_t::WEIGHT, tile_size_mac[data_type_t::WEIGHT], true);
            if(last_pe) move_front(&m_scheduler->weight_offset_pe);
        }
        weight_index++;
        exist_data_mac[data_type_t::WEIGHT] = true;
        request_to_lb[data_type_t::WEIGHT] = false;
        if(tile_size_mac[data_type_t::WEIGHT] == tile_size_lb[data_type_t::WEIGHT]) skip_transfer[data_type_t::WEIGHT] = true;
    }
    if(request_to_lb[data_type_t::OUTPUT]) {
        account_format_events(data_type_t::OUTPUT, tile_size_mac[data_type_t::OUTPUT]);
        if(m_scheduler->output_read_pe[m_scheduler->output_offset_pe.front()]) {
            // RE1/E20-4: a prior partial sum exists -- but that alone is not a reload. A reload is
            // the physical READ-BACK over the LB->MAC path, and that path is skipped when the
            // output tile is already resident in the MAC. Charging before the guard billed a
            // read-back that never happened, on every retained pass.
            if(!skip_transfer[data_type_t::OUTPUT]) {
                account_accumulator_reload(tile_size_mac[data_type_t::OUTPUT]);
                num_data_transfer_to_mac[data_type_t::OUTPUT]++;
                account_descriptor_dense_mac_transfer(data_type_t::OUTPUT, tile_size_mac[data_type_t::OUTPUT], true);
            } else {
                // E20-4: the accumulator stayed in the MAC. Counted so a retained pass is visible
                // as a retention, not as an absent event.
                ++accumulator_retained_events;
            }
            if(tile_size_mac[data_type_t::OUTPUT] == tile_size_lb[data_type_t::OUTPUT]) skip_transfer[data_type_t::OUTPUT] = true;
        } else {
            // RE1: a fresh accumulator is zero-initialized in place -- nothing is read back, so no
            // reload energy is charged here. This branch used to pay for a reload that never
            // happened.
            account_accumulator_create(tile_size_mac[data_type_t::OUTPUT]);
            clear_output_accumulators();
            m_scheduler->output_read_pe[m_scheduler->output_offset_pe.front()] = true;
        }
        output_index++;
        exist_data_mac[data_type_t::OUTPUT] = true;
        request_to_lb[data_type_t::OUTPUT] = false;
    }

    computation(m_scheduler);
    if(input_index == m_scheduler->offset_size_pe[data_type_t::INPUT].front() &&
       weight_index == m_scheduler->offset_size_pe[data_type_t::WEIGHT].front() &&
       output_index == m_scheduler->offset_size_pe[data_type_t::OUTPUT].front()) {
        flush_data(m_scheduler);
    }
    return;
#endif
    if(request_to_lb[data_type_t::INPUT]) {
        account_format_events(data_type_t::INPUT, tile_size_mac[data_type_t::INPUT]);
#ifdef FUNCTIONAL
        // Input data transfer
        m_scheduler->transfer_data(input_data_mac, input_data_lb, 0, m_scheduler->input_offset_pe.front(),
                                   component_type_t::MAC, component_type_t::PE,
                                   data_type_t::INPUT, get_mac_stationary_type(), action_type_t::LOAD);
        // Update for NPUsim ver2
        //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 &&
        //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
        //m_scheduler->transfer_data_ver2(input_data_mac, input_data_lb,
        //                                component_type_t::MAC, component_type_t::PE,
        //                                data_type_t::INPUT, get_mac_stationary_type(), action_type_t::LOAD, last_component);

        // Case 1. Dense format
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::INPUT]) {
                num_data_transfer_to_mac[data_type_t::INPUT]++;

                std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
                parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

                uint64_t address_mac = 0, address_lb = 0;
                unsigned num_access_mac = 0, num_access_lb = 0;
                for(unsigned b = 0; b < parameters_mac[parameter_type_t::BATCH_SIZE]; b++) {
                    for(unsigned c = 0; c < parameters_mac[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned h = 0; h < parameters_mac[parameter_type_t::INPUT_HEIGHT]; h++) {
                            for(unsigned w = 0; w < parameters_mac[parameter_type_t::INPUT_WIDTH]; w++) {
                                if(address_mac != ((uint64_t)&input_data_mac[b*parameters_mac[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_mac[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_mac[parameter_type_t::INPUT_WIDTH] +
                                                                             c*parameters_mac[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_mac[parameter_type_t::INPUT_WIDTH] +
                                                                             h*parameters_mac[parameter_type_t::INPUT_WIDTH] + w] >>
                                                                             mask_bits_mac[data_type_t::INPUT]) << mask_bits_mac[data_type_t::INPUT]) {
                                    // Update write cost to MAC unit
                                    access_energy_mac[data_type_t::INPUT] += u_write_energy_mac[data_type_t::INPUT];
                                    access_cycle_mac[data_type_t::INPUT] += u_write_cycle_mac[data_type_t::INPUT];
                                    num_access_mac++;

                                    // Update MAC address
                                    address_mac = ((uint64_t)&input_data_mac[b*parameters_mac[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_mac[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_mac[parameter_type_t::INPUT_WIDTH] +
                                                                             c*parameters_mac[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_mac[parameter_type_t::INPUT_WIDTH] +
                                                                             h*parameters_mac[parameter_type_t::INPUT_WIDTH] + w] >>
                                                                             mask_bits_mac[data_type_t::INPUT]) << mask_bits_mac[data_type_t::INPUT];
                                }
                                if(address_lb != ((uint64_t)&input_data_lb[m_scheduler->input_offset_pe.front() +
                                                                           b*parameters_lb[parameter_type_t::INPUT_CHANNEL]
                                                                            *parameters_lb[parameter_type_t::INPUT_HEIGHT]
                                                                            *parameters_lb[parameter_type_t::INPUT_WIDTH] +
                                                                           c*parameters_lb[parameter_type_t::INPUT_HEIGHT]
                                                                            *parameters_lb[parameter_type_t::INPUT_WIDTH] +
                                                                           h*parameters_lb[parameter_type_t::INPUT_WIDTH] + w] >>
                                                                           mask_bits_lb[data_type_t::INPUT]) << mask_bits_lb[data_type_t::INPUT]) {

                                    // Update read cost to local buffer.
                                    access_energy_lb[data_type_t::INPUT] += u_read_energy_lb[data_type_t::INPUT];
                                    access_cycle_lb[data_type_t::INPUT] += u_read_cycle_lb[data_type_t::INPUT];
                                    num_access_lb++;

                                    // Update local buffer address
                                    address_lb = ((uint64_t)&input_data_lb[m_scheduler->input_offset_pe.front() +
                                                                           b*parameters_lb[parameter_type_t::INPUT_CHANNEL]
                                                                            *parameters_lb[parameter_type_t::INPUT_HEIGHT]
                                                                            *parameters_lb[parameter_type_t::INPUT_WIDTH] +
                                                                           c*parameters_lb[parameter_type_t::INPUT_HEIGHT]
                                                                            *parameters_lb[parameter_type_t::INPUT_WIDTH] +
                                                                           h*parameters_lb[parameter_type_t::INPUT_WIDTH] + w] >>
                                                                           mask_bits_lb[data_type_t::INPUT]) << mask_bits_lb[data_type_t::INPUT];
                                }
                            }
                        }
                    }
                }

                // Update overlapped cycle at local buffer and MAC units.

                unsigned ratio = ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(line_size_mac[data_type_t::INPUT]));

                // At the 1, 2, before last, last stages
                double first_stage = u_read_cycle_lb[data_type_t::INPUT];
                double second_stage = std::max(u_read_cycle_lb[data_type_t::INPUT],
                                                 u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth)));
                double last_before_stage = std::max(ratio*u_write_cycle_mac[data_type_t::INPUT],
                                                      u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth)));
                double last_stage = ratio*u_write_cycle_mac[data_type_t::INPUT];

                // Remainder stages.
                double other_stage = std::max(u_read_cycle_lb[data_type_t::INPUT],
                                       std::max(ratio*u_write_cycle_mac[data_type_t::INPUT],
                                                u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth))));

                if(num_access_lb == 0) {
                    // No transaction: leave the overlapped stage unchanged.
                } else if(num_access_lb == 1) {
                    cycle_mac_lb[data_type_t::INPUT] += u_read_cycle_lb[data_type_t::INPUT] +
                                                        u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth)) +
                                                        ratio*u_write_cycle_mac[data_type_t::INPUT];
                } else {
                    cycle_mac_lb[data_type_t::INPUT] += first_stage + second_stage +
                                                        (num_access_lb-2)*other_stage +
                                                        last_before_stage + last_stage;
                }
                // Update per data stats between MAC unit and local buffer
                transfer_cycle[data_type_t::INPUT] += num_access_lb*u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth));
                transfer_energy[data_type_t::INPUT] += num_access_lb*u_transfer_energy*ceil((double)(line_size_lb[data_type_t::INPUT])/(double)(bitwidth));
            }
            if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                move_front(&m_scheduler->input_offset_pe);
            }
        }
        // Case 2. COO format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
            std::cerr << "Sparse format COO is not supported in this version" << std::endl;
            exit(1);
        }
        // Case 3. CSC format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::INPUT]) {
                // Row index calculation
                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::MAC);
                unsigned row = parameters[parameter_type_t::INPUT_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer_to_mac[data_type_t::INPUT]++;

                // Update local buffer access cost
                access_cycle_lb[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                      *u_read_cycle_lb[data_type_t::INPUT] + // Value
                                                       (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                      *u_read_cycle_lb[data_type_t::INPUT]
                                                      /(sizeof(data_t)*8/row_bit) + // Row index
                                                       parameters[parameter_type_t::BATCH_SIZE]
                                                      *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                      *(parameters[parameter_type_t::INPUT_WIDTH]+1)*u_read_cycle_lb[data_type_t::INPUT]
                                                      /(sizeof(data_t)*8/row_bit); // Column pointer

                access_energy_lb[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_read_energy_lb[data_type_t::INPUT] + // Non-zero data
                                                        (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_read_energy_lb[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/row_bit) + // Row index
                                                        parameters[parameter_type_t::BATCH_SIZE]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                       *u_read_energy_lb[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/row_bit); // Column pointer

                // Update MAC access cost
                access_cycle_mac[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_write_cycle_mac[data_type_t::INPUT] + // Non-zero data
                                                        (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_write_cycle_mac[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/row_bit) + // Row index
                                                        parameters[parameter_type_t::BATCH_SIZE]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                       *u_write_cycle_mac[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/row_bit); // Column pointer

                access_energy_mac[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                        *u_write_energy_mac[data_type_t::INPUT] + // DATA
                                                         (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                        *u_write_energy_mac[data_type_t::INPUT]
                                                        /(sizeof(data_t)*8/row_bit) + // Row index
                                                         parameters[parameter_type_t::BATCH_SIZE]
                                                        *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                        *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                        *u_write_energy_mac[data_type_t::INPUT]
                                                        /(sizeof(data_t)*8/row_bit); // Column pointer

                /* TODO */
                cycle_mac_lb[data_type_t::INPUT] += std::max((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_cycle_lb[data_type_t::INPUT] +  // Non-zero data
                                                             (tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_cycle_lb[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/row_bit)  + // Row index
                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                            *u_read_cycle_lb[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/row_bit), // Column pointer
                                                             (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_write_cycle_mac[data_type_t::INPUT] + // Non-zero data
                                                             (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_write_cycle_mac[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/row_bit) + // Row index
                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                            *u_write_cycle_mac[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/row_bit)); // Column pointer

                // Update transfer cost between MAC units and local buffers
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                      u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)/(double)bitwidth) + // Row index
                                                      u_transfer_cycle*ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *(parameters[parameter_type_t::INPUT_WIDTH+1])*row_bit)/(double)bitwidth);  // Column pointer
                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)/(double)bitwidth) + // Row index
                                                       u_transfer_energy*ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                      *parameters[parameter_type_t::INPUT_CHANNEL]
                                                                                      /parameters[parameter_type_t::GROUP]
                                                                                      *(parameters[parameter_type_t::INPUT_WIDTH+1])*row_bit)/(double)bitwidth); // Column pointer
                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->input_offset_pe);
                }
            }
        }
        // Case 4. CSR format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::INPUT]) {
                // Column index calculation
                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::MAC);
                unsigned column = parameters[parameter_type_t::INPUT_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer_to_mac[data_type_t::INPUT]++;

                // Update local buffer access cost
                access_cycle_lb[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                      *u_read_cycle_lb[data_type_t::INPUT] + // Non-zero data
                                                       (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                      *u_read_cycle_lb[data_type_t::INPUT]
                                                      /(sizeof(data_t)*8/column_bit) + // Column index
                                                       parameters[parameter_type_t::BATCH_SIZE]
                                                      *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                      *(parameters[parameter_type_t::INPUT_HEIGHT]+1)*u_read_cycle_lb[data_type_t::INPUT]
                                                      /(sizeof(data_t)*8/column_bit); // Row pointer

                access_energy_lb[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_read_energy_lb[data_type_t::INPUT] + // Non-zero data
                                                        (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_read_energy_lb[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/column_bit) + // Column index
                                                        parameters[parameter_type_t::BATCH_SIZE]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                       *u_read_energy_lb[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/column_bit); // Row pointer

                // Update MAC access cost
                access_cycle_mac[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_write_cycle_mac[data_type_t::INPUT] + // Non-zero data
                                                        (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_write_cycle_mac[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/column_bit) + // Column index
                                                        parameters[parameter_type_t::BATCH_SIZE]
                                                       *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                       *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                       *u_write_cycle_mac[data_type_t::INPUT]
                                                       /(sizeof(data_t)*8/column_bit); // Row pointer

                access_energy_mac[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                        *u_write_energy_mac[data_type_t::INPUT] + // DATA
                                                         (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                        *u_write_energy_mac[data_type_t::INPUT]
                                                        /(sizeof(data_t)*8/column_bit) + // Column index
                                                         parameters[parameter_type_t::BATCH_SIZE]
                                                        *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                        *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                        *u_write_energy_mac[data_type_t::INPUT]
                                                        /(sizeof(data_t)*8/column_bit);

                /* TODO */
                cycle_mac_lb[data_type_t::INPUT] += std::max((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_cycle_lb[data_type_t::INPUT] +  // Non-zero data
                                                             (tile_size_lb[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_read_cycle_lb[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/column_bit)  + // Column index
                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                            *u_read_cycle_lb[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/column_bit), // Row pointer
                                                             (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_write_cycle_mac[data_type_t::INPUT] + // Non-zero data
                                                             (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                            *u_write_cycle_mac[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/column_bit) + // Column index
                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                            *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                            *u_write_cycle_mac[data_type_t::INPUT]
                                                            /(sizeof(data_t)*8/column_bit)); // Row pointer

                // Update transfer cost between MAC units and local buffers
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                      u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)/(double)bitwidth) + // Column index
                                                      u_transfer_cycle*ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *(parameters[parameter_type_t::INPUT_HEIGHT+1])*column_bit)/(double)bitwidth);  // Row pointer
                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)/(double)bitwidth) + // Column index
                                                       u_transfer_energy*ceil((double)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                      *parameters[parameter_type_t::INPUT_CHANNEL]
                                                                                      /parameters[parameter_type_t::GROUP]
                                                                                      *(parameters[parameter_type_t::INPUT_HEIGHT+1])*column_bit)/(double)bitwidth); // Row pointer

                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->input_offset_pe);
                }
            }
        }
        // Case 5. SparseMap
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::INPUT]) {
                num_data_transfer_to_mac[data_type_t::INPUT]++;

                // Update local buffer access cost
                access_cycle_lb[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                      *u_read_cycle_lb[data_type_t::INPUT] + // Non-zero value
                                                       tile_size_mac[data_type_t::INPUT]
                                                      *u_read_cycle_lb[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
                access_energy_lb[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_read_energy_lb[data_type_t::INPUT] + // Non-zero data
                                                        tile_size_mac[data_type_t::INPUT]
                                                       *u_read_energy_lb[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata

                // Update MAC access cost
                access_cycle_mac[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                       *u_write_cycle_mac[data_type_t::INPUT] + // Non-zero data
                                                        tile_size_mac[data_type_t::INPUT]
                                                       *u_write_cycle_mac[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
                access_energy_mac[data_type_t::INPUT] += (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*u_write_energy_mac[data_type_t::INPUT]
                                                       + tile_size_mac[data_type_t::INPUT]*u_write_energy_mac[data_type_t::INPUT]/(sizeof(data_t)*8);

                // Update overlapped cycle between MAC and local buffers
                cycle_mac_lb[data_type_t::INPUT] += std::max((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *u_read_cycle_lb[data_type_t::INPUT] + // Non-zero data
                                                              tile_size_mac[data_type_t::INPUT]
                                                             *u_read_cycle_lb[data_type_t::INPUT]/(sizeof(data_t)*8), // Metadata
                                                              (tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *u_write_cycle_mac[data_type_t::INPUT] +  // Non-zero data
                                                              tile_size_mac[data_type_t::INPUT]
                                                             *u_write_cycle_mac[data_type_t::INPUT]/(sizeof(data_t)*8)); // Metadata

                // Transfer cycle between MAC unit and local buffer = transfer cycle for data + transfer cycle for index vector.
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(double)bitwidth) +  // Non-zero data
                                                      u_transfer_cycle*ceil((double)(tile_size_mac[data_type_t::INPUT])/(double)bitwidth); // Metadata

                // Transfer energy between MAC unit and local buffer = transfer energy for data + transfer energy for index vector
                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((double)(tile_size_mac[data_type_t::INPUT])/(double)bitwidth); // Metadata
                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->input_offset_pe);
                }
            }
        }
        else {
            std::cerr << "Undefined compression type" << std::endl;
            exit(1);
        }

#endif

        // Increment the input index that finds the offset of input data in local buffer.
        input_index++;
        // Input data exists in MAC unit and request input data is not required.
        exist_data_mac[data_type_t::INPUT] = true, request_to_lb[data_type_t::INPUT] = false;

        if(tile_size_mac[data_type_t::INPUT] == tile_size_lb[data_type_t::INPUT]) { skip_transfer[data_type_t::INPUT] = true; }
    }
    // Transfer weight from local buffer to MAC unit.
    if(request_to_lb[data_type_t::WEIGHT]) {
        account_format_events(data_type_t::WEIGHT, tile_size_mac[data_type_t::WEIGHT]);
#ifdef FUNCTIONAL

        // Weight data transfer
        m_scheduler->transfer_data(weight_mac, weight_lb, 0, m_scheduler->weight_offset_pe.front(),
                                   component_type_t::MAC, component_type_t::PE,
                                   data_type_t::WEIGHT, get_mac_stationary_type(), action_type_t::LOAD);
        // Update for NPUsim ver2
        //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 &&
        //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
        //m_scheduler->transfer_data_ver2(weight_mac, weight_lb,
        //                                component_type_t::MAC, component_type_t::PE,
        //                                data_type_t::WEIGHT, get_mac_stationary_type(), action_type_t::LOAD, last_component);

        // Case 1. Dense data format
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                num_data_transfer_to_mac[data_type_t::WEIGHT]++;

                std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
                parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

                uint64_t address_mac = 0, address_lb = 0;
                unsigned num_access_mac = 0, num_access_lb = 0;
                for(unsigned k = 0; k < parameters_mac[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned c = 0; c < parameters_mac[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned r = 0; r < parameters_mac[parameter_type_t::FILTER_HEIGHT]; r++) {
                            for(unsigned s = 0; s < parameters_mac[parameter_type_t::FILTER_WIDTH]; s++) {
                                if(address_mac != ((uint64_t)&weight_mac[k*parameters_mac[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_mac[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_mac[parameter_type_t::FILTER_WIDTH] +
                                                                         c*parameters_mac[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_mac[parameter_type_t::FILTER_WIDTH] +
                                                                         r*parameters_mac[parameter_type_t::FILTER_WIDTH] + s] >>
                                                                         mask_bits_mac[data_type_t::WEIGHT]) << mask_bits_mac[data_type_t::WEIGHT]) {

                                    // Update MAC cost.
                                    access_energy_mac[data_type_t::WEIGHT] += u_write_energy_mac[data_type_t::WEIGHT];
                                    access_cycle_mac[data_type_t::WEIGHT] += u_write_cycle_mac[data_type_t::WEIGHT];
                                    num_access_mac++;

                                    // Update MAC address
                                    address_mac = ((uint64_t)&weight_mac[k*parameters_mac[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_mac[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_mac[parameter_type_t::FILTER_WIDTH] +
                                                                         c*parameters_mac[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_mac[parameter_type_t::FILTER_WIDTH] +
                                                                         r*parameters_mac[parameter_type_t::FILTER_WIDTH] + s] >>
                                                                         mask_bits_mac[data_type_t::WEIGHT]) << mask_bits_mac[data_type_t::WEIGHT];
                                }

                                if(address_lb != ((uint64_t)&weight_lb[m_scheduler->weight_offset_pe.front() +
                                                                       k*parameters_lb[parameter_type_t::INPUT_CHANNEL]
                                                                        *parameters_lb[parameter_type_t::FILTER_HEIGHT]
                                                                        *parameters_lb[parameter_type_t::FILTER_WIDTH] +
                                                                       c*parameters_lb[parameter_type_t::FILTER_HEIGHT]
                                                                        *parameters_lb[parameter_type_t::FILTER_WIDTH] +
                                                                       r*parameters_lb[parameter_type_t::FILTER_WIDTH] + s] >>
                                                                       mask_bits_lb[data_type_t::WEIGHT]) << mask_bits_lb[data_type_t::WEIGHT]) {

                                    access_energy_lb[data_type_t::WEIGHT] += u_read_energy_lb[data_type_t::WEIGHT];
                                    access_cycle_lb[data_type_t::WEIGHT] += u_read_cycle_lb[data_type_t::WEIGHT];
                                    num_access_lb++;

                                    // Update Local buffer address
                                    address_lb = ((uint64_t)&weight_lb[m_scheduler->weight_offset_pe.front() +
                                                                       k*parameters_lb[parameter_type_t::INPUT_CHANNEL]
                                                                        *parameters_lb[parameter_type_t::FILTER_HEIGHT]
                                                                        *parameters_lb[parameter_type_t::FILTER_WIDTH] +
                                                                       c*parameters_lb[parameter_type_t::FILTER_HEIGHT]
                                                                        *parameters_lb[parameter_type_t::FILTER_WIDTH] +
                                                                       r*parameters_lb[parameter_type_t::FILTER_WIDTH] + s] >>
                                                                       mask_bits_lb[data_type_t::WEIGHT]) << mask_bits_lb[data_type_t::WEIGHT];
                                }
                            }
                        }
                    }
                }

                // Update overlapped cycle at local buffer and MAC units.
                unsigned ratio = ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(line_size_mac[data_type_t::WEIGHT]));

                // At the 1, 2, before last, last stages
                double first_stage = u_read_cycle_lb[data_type_t::WEIGHT];
                double second_stage = std::max(u_read_cycle_lb[data_type_t::WEIGHT],
                                                 u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)));
                double last_before_stage = std::max(ratio*u_write_cycle_mac[data_type_t::WEIGHT],
                                                      u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)));
                double last_stage = ratio*u_write_cycle_mac[data_type_t::WEIGHT];

                // Remainder stages
                double other_stage = std::max(u_read_cycle_lb[data_type_t::WEIGHT],
                                       std::max(ratio*u_write_cycle_mac[data_type_t::WEIGHT],
                                                u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth))));

                if(num_access_lb == 0) {
                    // No transaction: leave the overlapped stage unchanged.
                } else if(num_access_lb == 1) {
                    cycle_mac_lb[data_type_t::WEIGHT] += u_read_cycle_lb[data_type_t::WEIGHT] +
                                                         u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)(bitwidth)) +
                                                         ratio*u_write_cycle_mac[data_type_t::WEIGHT];
                } else {
                    cycle_mac_lb[data_type_t::WEIGHT] += first_stage + second_stage +
                                                         (num_access_lb-2)*other_stage +
                                                         last_before_stage + last_stage;
                }

                // Update per data stats between MAC unit and local buffer
                transfer_cycle[data_type_t::WEIGHT] += num_access_lb*u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)bitwidth);
                transfer_energy[data_type_t::WEIGHT] += num_access_lb*u_transfer_energy*ceil((double)(line_size_lb[data_type_t::WEIGHT])/(double)bitwidth);

                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->weight_offset_pe);
                }
            }
        }
        // Case 2. COO data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
            std::cerr << "Data format COO is not supported in this version." << std::endl;
            exit(1);
        }
        // Case 3. CSC format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                // Row bit calculation
                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::MAC);
                unsigned row = parameters[parameter_type_t::FILTER_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer_to_mac[data_type_t::WEIGHT]++;

                // Update local buffer access cost
                access_cycle_lb[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle_lb[data_type_t::WEIGHT] + // Non-zero data
                                                        (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle_lb[data_type_t::WEIGHT]
                                                       /(sizeof(data_t)*8/row_bit) +  // Row index
                                                       parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                      *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                      *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                      *u_read_cycle_lb[data_type_t::WEIGHT]
                                                      /(sizeof(data_t)*8/row_bit); // Column pointer
                access_energy_lb[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_read_energy_lb[data_type_t::WEIGHT] + // Non-zero data
                                                         (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_read_energy_lb[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/row_bit) + // Row index
                                                         parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                        *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                        *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                        *u_read_energy_lb[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/row_bit); // Column pointer

                // Update MAC access cost
                access_cycle_mac[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_write_cycle_mac[data_type_t::WEIGHT] + // Non-zero data
                                                         (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_write_cycle_mac[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/row_bit) + // Row index
                                                         parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                        *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                        *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                        *u_write_cycle_mac[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/row_bit); // Column pointer
                access_energy_mac[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                         *u_write_energy_mac[data_type_t::WEIGHT] + // Non-zero data
                                                          (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                         *u_write_energy_mac[data_type_t::WEIGHT]
                                                         /(sizeof(data_t)*8/row_bit) + // Row index
                                                          parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                         /parameters[parameter_type_t::GROUP]
                                                         *parameters[parameter_type_t::INPUT_CHANNEL]
                                                         /parameters[parameter_type_t::GROUP]
                                                         *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                         *u_write_energy_mac[data_type_t::WEIGHT]
                                                         /(sizeof(data_t)*8/row_bit); // Column pointer

                /* TODO */
                // Update overlapped cycle between MAC units and local buffers
                cycle_mac_lb[data_type_t::WEIGHT] += std::max((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_cycle_lb[data_type_t::WEIGHT] + // Non-zero data
                                                              (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_cycle_lb[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/row_bit) + // Row index
                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                             *u_read_cycle_lb[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/row_bit), // Column pointer
                                                              (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_write_cycle_mac[data_type_t::WEIGHT] + // Non-zero data
                                                              (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_write_cycle_mac[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/row_bit) + // Row index
                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                             *u_write_cycle_mac[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/row_bit)); // Column pointer

                // Update transfer cost between MAC units and local buffers
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                       u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*row_bit)/(double)bitwidth) +  // Row index
                                                       u_transfer_cycle*ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                     *(parameters[parameter_type_t::FILTER_WIDTH+1])
                                                                                     *row_bit)/(double)bitwidth); // Column pointer
                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                        u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*row_bit)/(double)bitwidth) + // Row index
                                                        u_transfer_energy*ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                                                       /parameters[parameter_type_t::GROUP]
                                                                                       *parameters[parameter_type_t::INPUT_CHANNEL]
                                                                                       /parameters[parameter_type_t::GROUP]
                                                                                       *(parameters[parameter_type_t::FILTER_WIDTH+1])
                                                                                       *row_bit)/(double)bitwidth); // Column pointer
                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->weight_offset_pe);
                }
            }
        }
        // Case 4. CSR format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                // Column bit calculation
                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::MAC);
                unsigned column = parameters[parameter_type_t::FILTER_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer_to_mac[data_type_t::WEIGHT]++;

                // Update local buffer access cost
                access_cycle_lb[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle_lb[data_type_t::WEIGHT] + // Non-zero data
                                                        (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle_lb[data_type_t::WEIGHT]
                                                       /(sizeof(data_t)*8/column_bit) +  // Column index
                                                       parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                      /parameters[parameter_type_t::GROUP]
                                                      *parameters[parameter_type_t::INPUT_CHANNEL]
                                                      /parameters[parameter_type_t::GROUP]
                                                      *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                      *u_read_cycle_lb[data_type_t::WEIGHT]
                                                      /(sizeof(data_t)*8/column_bit); // Row pointer
                access_energy_lb[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_read_energy_lb[data_type_t::WEIGHT] + // Non-zero data
                                                         (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_read_energy_lb[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/column_bit) + // Column index
                                                         parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                        /parameters[parameter_type_t::GROUP]
                                                        *parameters[parameter_type_t::INPUT_CHANNEL]
                                                        /parameters[parameter_type_t::GROUP]
                                                        *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                        *u_read_energy_lb[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/column_bit); // Row pointer

                // Update MAC access cost
                access_cycle_mac[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_write_cycle_mac[data_type_t::WEIGHT] + // Non-zero data
                                                         (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_write_cycle_mac[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/column_bit) + // Column index
                                                         parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                        /parameters[parameter_type_t::GROUP]
                                                        *parameters[parameter_type_t::INPUT_CHANNEL]
                                                        /parameters[parameter_type_t::GROUP]
                                                        *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                        *u_write_cycle_mac[data_type_t::WEIGHT]
                                                        /(sizeof(data_t)*8/column_bit); // Row pointer
                access_energy_mac[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                         *u_write_energy_mac[data_type_t::WEIGHT] + // Non-zero data
                                                          (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                         *u_write_energy_mac[data_type_t::WEIGHT]
                                                         /(sizeof(data_t)*8/column_bit) + // Column index
                                                          parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                         /parameters[parameter_type_t::GROUP]
                                                         *parameters[parameter_type_t::INPUT_CHANNEL]
                                                         /parameters[parameter_type_t::GROUP]
                                                         *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                         *u_write_energy_mac[data_type_t::WEIGHT]
                                                         /(sizeof(data_t)*8/column_bit); // Row pointer

                /* TODO */
                // Update overlapped cycle between MAC units and local buffers
                cycle_mac_lb[data_type_t::WEIGHT] += std::max((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_cycle_lb[data_type_t::WEIGHT] + // Non-zero data
                                                              (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_read_cycle_lb[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/column_bit) + // Column index
                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                             /parameters[parameter_type_t::GROUP]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]
                                                             /parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                             *u_read_cycle_lb[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/column_bit), // Row pointer
                                                              (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_write_cycle_mac[data_type_t::WEIGHT] + // Non-zero data
                                                              (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                             *u_write_cycle_mac[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/column_bit) + // Column index
                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                             /parameters[parameter_type_t::GROUP]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]
                                                             /parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                             *u_write_cycle_mac[data_type_t::WEIGHT]
                                                             /(sizeof(data_t)*8/column_bit)); // Row pointer

                // Update transfer cost between MAC units and local buffers
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                       u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*column_bit)/(double)bitwidth) +  // Column index
                                                       u_transfer_cycle*ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                                                     /parameters[parameter_type_t::GROUP]
                                                                                     *parameters[parameter_type_t::INPUT_CHANNEL]
                                                                                     /parameters[parameter_type_t::GROUP]
                                                                                     *(parameters[parameter_type_t::FILTER_HEIGHT+1])
                                                                                     *column_bit)/(double)bitwidth); // Row pointer
                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                        u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*column_bit)/(double)bitwidth) + // Column index
                                                        u_transfer_energy*ceil((double)(parameters[parameter_type_t::OUTPUT_CHANNEL]
                                                                                       /parameters[parameter_type_t::GROUP]
                                                                                       *parameters[parameter_type_t::INPUT_CHANNEL]
                                                                                       /parameters[parameter_type_t::GROUP]
                                                                                       *(parameters[parameter_type_t::FILTER_HEIGHT+1])
                                                                                       *column_bit)/(double)bitwidth); // Row pointer

                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->weight_offset_pe);
                }
            }
        }
        // Case 4. SparseMap
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                num_data_transfer_to_mac[data_type_t::WEIGHT]++;

                // Update local buffer access cost
                access_cycle_lb[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                       *u_read_cycle_lb[data_type_t::WEIGHT] + // Non-zero data
                                                        tile_size_mac[data_type_t::WEIGHT]
                                                       *u_read_cycle_lb[data_type_t::WEIGHT]/(sizeof(data_t)*8); // Metadata
                access_energy_lb[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_read_energy_lb[data_type_t::WEIGHT] + // Non-zero data
                                                         tile_size_mac[data_type_t::WEIGHT]
                                                        *u_read_energy_lb[data_type_t::WEIGHT]/(sizeof(data_t)*8); // Metadata
                // Update MAC access cost
                access_cycle_mac[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                        *u_write_cycle_mac[data_type_t::WEIGHT] + // Non-zero data
                                                         tile_size_mac[data_type_t::WEIGHT]
                                                        *u_write_cycle_mac[data_type_t::WEIGHT]/(sizeof(data_t)*8); // Metadata

                access_energy_mac[data_type_t::WEIGHT] += (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*u_write_energy_mac[data_type_t::WEIGHT]
                                                        + tile_size_mac[data_type_t::WEIGHT]*u_write_energy_mac[data_type_t::WEIGHT]/(sizeof(data_t)*8);

                cycle_mac_lb[data_type_t::WEIGHT] += std::max((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*u_read_cycle_lb[data_type_t::WEIGHT]
                                                                         + tile_size_mac[data_type_t::WEIGHT]*u_read_cycle_lb[data_type_t::WEIGHT]/(sizeof(data_t)*8),
                                                                         (tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*u_write_cycle_mac[data_type_t::WEIGHT]
                                                                         + tile_size_mac[data_type_t::WEIGHT]*u_write_cycle_mac[data_type_t::WEIGHT]/(sizeof(data_t)*8));

                // Transfer cycle between MAC unit and Local buffer = transfer cycle for data + transfer cycle for index vector
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(double)bitwidth)
                                                     + u_transfer_cycle*ceil((double)(tile_size_mac[data_type_t::WEIGHT])/(double)bitwidth);
                // Transfer energy between MAC unit and local buffer = transfer energy for data + transfer energy for index vector
                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy*ceil((double)((tile_size_mac[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(double)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((double)(tile_size_mac[data_type_t::WEIGHT])/(double)bitwidth); // Metadata

                if(pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1 && index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1) {
                    move_front(&m_scheduler->weight_offset_pe);
                }
            }
        }
        else {
            std::cerr << "Undefined compression type" << std::endl;
            exit(1);
        }
#endif
        /* Stats */

        // Increase the weight index that finds the offset of weight in local buffer.
        weight_index++;
        // Weight exists in MAC unit and request weight is not required.
        exist_data_mac[data_type_t::WEIGHT] = true, request_to_lb[data_type_t::WEIGHT] = false;

        if(tile_size_mac[data_type_t::WEIGHT] == tile_size_lb[data_type_t::WEIGHT]) { skip_transfer[data_type_t::WEIGHT] = true;}
    }
    // Transfer output data from local buffer to MAC unit.
    if(request_to_lb[data_type_t::OUTPUT]) {
        account_format_events(data_type_t::OUTPUT, tile_size_mac[data_type_t::OUTPUT]);
        // Load output data from local buffer to MAC unit.
        if(m_scheduler->output_read_pe[m_scheduler->output_offset_pe.front()]) {
            // RE1/E20-4: a prior partial sum exists -- but that alone is not a reload. A reload is
            // the physical READ-BACK over the LB->MAC path, and that path is skipped when the
            // output tile is already resident in the MAC. Charging before the guard billed a
            // read-back that never happened, on every retained pass.
            if(!skip_transfer[data_type_t::OUTPUT]) {
                account_accumulator_reload(tile_size_mac[data_type_t::OUTPUT]);
#ifdef FUNCTIONAL

                // Output data transfer
                m_scheduler->transfer_data(output_data_mac, output_data_lb, 0, m_scheduler->output_offset_pe.front(),
                                           component_type_t::MAC, component_type_t::PE,
                                           data_type_t::OUTPUT, get_mac_stationary_type(), action_type_t::LOAD);
                // Update for NPUsim ver2
                //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 &&
                //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
                //m_scheduler->transfer_data_ver2(output_data_mac, output_data_lb,
                //                                component_type_t::MAC, component_type_t::PE,
                //                                data_type_t::OUTPUT, get_mac_stationary_type(), action_type_t::LOAD, last_component);
#endif

                // Update the number of output data transfer from local buffer to MAC unit
                num_data_transfer_to_mac[data_type_t::OUTPUT]++;

                std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
                parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

                uint64_t address_mac = 0, address_lb = 0;
                unsigned num_access_mac = 0, num_access_lb = 0;
                for(unsigned b = 0; b < parameters_mac[parameter_type_t::BATCH_SIZE]; b++) {
                    for(unsigned k = 0; k < parameters_mac[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                        for(unsigned p = 0; p < parameters_mac[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                            for(unsigned q = 0; q < parameters_mac[parameter_type_t::OUTPUT_WIDTH]; q++) {
                                if(address_mac != ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                               *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                               *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                              k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                               *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                              p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                              mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT]) {

                                    // Update write cost of MAC unit.
                                    access_energy_mac[data_type_t::OUTPUT] += u_write_energy_mac[data_type_t::OUTPUT];
                                    access_cycle_mac[data_type_t::OUTPUT] += u_write_cycle_mac[data_type_t::OUTPUT];
                                    num_access_mac++;

                                    // Update MAC address
                                    address_mac = ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                               *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                               *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                              k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                               *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                              p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                              mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT];
                                }

                                if(address_lb != ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                            b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                             *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                             *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                            k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                             *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                            p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                            mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                    // Update read cost of Local buffer
                                    access_energy_lb[data_type_t::OUTPUT] += u_read_energy_lb[data_type_t::OUTPUT];
                                    access_cycle_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT];
                                    num_access_lb++;

                                    // Update local buffer address
                                    address_lb = ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                            b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                             *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                             *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                            k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                             *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                            p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                            mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];

                                }
                            }
                        }
                    }
                }

                // Update overlapped cycle at local buffer and MAC units.
                unsigned ratio = ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(line_size_mac[data_type_t::OUTPUT]));

                // at the 1, 2, before last, and last stages
                double first_stage = u_read_cycle_lb[data_type_t::OUTPUT];
                double second_stage = std::max(u_read_cycle_lb[data_type_t::OUTPUT],
                                                 u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
                double last_before_stage = std::max(ratio*u_write_cycle_mac[data_type_t::OUTPUT],
                                                      u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
                double last_stage = ratio*u_write_cycle_mac[data_type_t::OUTPUT];

                // Remainder stages
                double other_stage = std::max(u_read_cycle_lb[data_type_t::OUTPUT],
                                       std::max(ratio*u_write_cycle_mac[data_type_t::OUTPUT],
                                                u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth))));

                if(num_access_lb == 0) {
                    // No transaction: leave the overlapped stage unchanged.
                } else if(num_access_lb == 1) {
                    cycle_mac_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT] +
                                                         u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)) +
                                                         ratio*u_write_cycle_mac[data_type_t::OUTPUT];
                } else {
                    cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage +
                                                         (num_access_lb-2)*other_stage +
                                                         last_before_stage + last_stage;
                }

                // Update data transfer cycle and energy between MAC and local buffer.
                transfer_cycle[data_type_t::OUTPUT] += num_access_lb*u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth));
                transfer_energy[data_type_t::OUTPUT] += num_access_lb*u_transfer_energy*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth));
            } else {
                // E20-4: keep the FUNCTIONAL and analytical paths on the same event contract.
                // The partial sum stayed in the MAC, so this is retention rather than a reload.
                ++accumulator_retained_events;
            }

            if(tile_size_mac[data_type_t::OUTPUT] == tile_size_lb[data_type_t::OUTPUT]) { skip_transfer[data_type_t::OUTPUT] = true;}
            /* Stats */
        }
        // This output has no prior partial sum: start from zero rather than
        // retaining the accumulator used for the preceding output coordinate.
        else {
            // RE1: a fresh accumulator is zero-initialized in place -- no reload energy.
            account_accumulator_create(tile_size_mac[data_type_t::OUTPUT]);
            clear_output_accumulators();
            m_scheduler->output_read_pe[m_scheduler->output_offset_pe.front()] = true;
        }

        // Increment the output data index that finds the offset of output data in local buffer.
        output_index++;
        // Output data exists in MAC unit and request output data is not required.
        exist_data_mac[data_type_t::OUTPUT]  = true, request_to_lb[data_type_t::OUTPUT] = false;
    }

    // Execute MAC operation.
    computation(m_scheduler);

    // Flush the data in local buffer.
    // When all the data in local buffer are used.
    if(input_index == m_scheduler->offset_size_pe[data_type_t::INPUT].front() &&
       weight_index == m_scheduler->offset_size_pe[data_type_t::WEIGHT].front() &&
       output_index == m_scheduler->offset_size_pe[data_type_t::OUTPUT].front()) {
        flush_data(m_scheduler);
    }
}

// Flush the data in local buffer.
void pe_t::flush_data(scheduler_t *m_scheduler) {
#ifndef FUNCTIONAL
    if(m_scheduler->compression_type != compression_type_t::DENSE) {
        std::cerr << "Error: timing PE flush supports dense descriptor traffic only" << std::endl;
        exit(1);
    }
    bool write_output = false;
    if(stationary_type_local_buffer == stationary_type_t::INPUT_STATIONARY) {
        write_output = true;
        if(weight_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front() - 1 &&
           output_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front() - 1) {
            exist_data_lb[data_type_t::WEIGHT] = false;
            exist_data_lb[data_type_t::OUTPUT] = false;
            weight_flush_counter++;
            output_flush_counter++;
        } else {
            exist_data_lb[data_type_t::INPUT] = false;
            exist_data_lb[data_type_t::WEIGHT] = false;
            exist_data_lb[data_type_t::OUTPUT] = false;
            weight_flush_counter = 0;
            output_flush_counter = 0;
        }
    } else if(stationary_type_local_buffer == stationary_type_t::WEIGHT_STATIONARY) {
        write_output = true;
        if(input_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front() - 1 &&
           output_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front() - 1) {
            exist_data_lb[data_type_t::INPUT] = false;
            exist_data_lb[data_type_t::OUTPUT] = false;
            input_flush_counter++;
            output_flush_counter++;
        } else {
            exist_data_lb[data_type_t::INPUT] = false;
            exist_data_lb[data_type_t::WEIGHT] = false;
            exist_data_lb[data_type_t::OUTPUT] = false;
            input_flush_counter = 0;
            output_flush_counter = 0;
        }
    } else if(stationary_type_local_buffer == stationary_type_t::OUTPUT_STATIONARY) {
        if(input_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front() - 1 &&
           weight_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front() - 1) {
            exist_data_lb[data_type_t::INPUT] = false;
            exist_data_lb[data_type_t::WEIGHT] = false;
            input_flush_counter++;
            weight_flush_counter++;
        } else {
            exist_data_lb[data_type_t::INPUT] = false;
            exist_data_lb[data_type_t::WEIGHT] = false;
            exist_data_lb[data_type_t::OUTPUT] = false;
            input_flush_counter = 0;
            weight_flush_counter = 0;
            write_output = true;
        }
    } else {
        exist_data_lb[data_type_t::INPUT] = false;
        exist_data_lb[data_type_t::WEIGHT] = false;
        exist_data_lb[data_type_t::OUTPUT] = false;
        write_output = true;
    }
    // RE1: the final cast is NOT charged here. A PE writes the same output tile back once per
    // reduction pass, so a cast charged at this boundary would follow the MAC count rather than the
    // final output element count. It is charged where the completed accumulation is read out of the
    // chip instead (global_buffer_t::account_output_writeback_link()).
    if(write_output) pe_array->account_descriptor_dense_writeback(this, tile_size_lb[data_type_t::OUTPUT]);
    wait_data();
    input_index = 0;
    weight_index = 0;
    output_index = 0;
    return;
#endif
    // Case 1. Input stationary
    if(stationary_type_local_buffer == stationary_type_t::INPUT_STATIONARY) {
        // Case 1) Reuse input data && Request weight and output data
        if(weight_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front() - 1 &&
           output_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front() - 1) {
            // Weight and output data do not exist in local buffer.
            exist_data_lb[data_type_t::WEIGHT] = false, exist_data_lb[data_type_t::OUTPUT] = false;

#ifdef FUNCTIONAL
            // Output data transfer
            m_scheduler->transfer_data(pe_array->output_data, output_data_lb, m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()], 0,
                                       component_type_t::PE_Y, component_type_t::PE,
                                       data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
            // Update for NPUsim ver2
            //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 &&
            //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
            //m_scheduler->transfer_data_ver2(pe_array->output_data, output_data_lb,
            //                                component_type_t::PE_Y, component_type_t::PE,
            //                                data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

            uint64_t address_pe = 0, address_pe_array = 0;
            unsigned num_access_pe = 0, num_access_pe_array = 0;
            for(unsigned b = 0; b < parameters_pe[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_pe[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_pe[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_pe[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            if(address_pe != ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                // Update read cost of Local buffer
                                access_energy_lb[data_type_t::OUTPUT] += u_read_energy_lb[data_type_t::OUTPUT];
                                access_cycle_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT];
                                num_access_pe++;

                                // Update local buffer address
                                address_pe = ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];
                            }

                            if(address_pe_array != ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT]) {

                                // Update write cost of accumulator in PE array.
                                pe_array->access_energy[data_type_t::OUTPUT] += pe_array->u_write_energy[data_type_t::OUTPUT];
                                pe_array->access_cycle[data_type_t::OUTPUT] += pe_array->u_write_cycle[data_type_t::OUTPUT];
                                num_access_pe_array++;

                                // Update accumulator address
                                address_pe_array = ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }

            // Increase flush counter of weight and output data.
            weight_flush_counter++;
            output_flush_counter++;

            // Waiting for weight and output data.
            // PE is in idle state.
            wait_data();
        }
        // Case 2) Request all data types
        else {
            // Input data, weight, and output data do not exist in local buffer.
            exist_data_lb[data_type_t::INPUT] = false; exist_data_lb[data_type_t::WEIGHT] = false; exist_data_lb[data_type_t::OUTPUT] = false;

#ifdef FUNCTIONAL
            // Write back output data
            m_scheduler->transfer_data(pe_array->output_data, output_data_lb, m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()], 0,
                                       component_type_t::PE_Y, component_type_t::PE,
                                       data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
            // Update for NPUsim ver2
            //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 &&
            //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
            //m_scheduler->transfer_data_ver2(pe_array->output_data, output_data_lb,
            //                                component_type_t::PE_Y, component_type_t::PE,
            //                                data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

            uint64_t address_pe = 0, address_pe_array = 0;
            unsigned num_access_pe = 0, num_access_pe_array = 0;
            for(unsigned b = 0; b < parameters_pe[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_pe[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_pe[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_pe[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            if(address_pe != ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                // Update read cost of Local buffer
                                access_energy_lb[data_type_t::OUTPUT] += u_read_energy_lb[data_type_t::OUTPUT];
                                access_cycle_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT];
                                num_access_pe++;

                                // Update local buffer address
                                address_pe = ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];
                            }

                            if(address_pe_array != ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT]) {

                                // Update write cost of accumulator in PE array.
                                pe_array->access_energy[data_type_t::OUTPUT] += pe_array->u_write_energy[data_type_t::OUTPUT];
                                pe_array->access_cycle[data_type_t::OUTPUT] += pe_array->u_write_cycle[data_type_t::OUTPUT];
                                num_access_pe_array++;

                                // Update accumulator address
                                address_pe_array = ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }

            // Set flush counter of weight and output data as zero.
            weight_flush_counter = 0;
            output_flush_counter = 0;

            // Waiting for input data, weight, and output data.
            // PE is in idle state.
            wait_data();
        }
    }
    // Case 2. Weight stationary
    else if(stationary_type_local_buffer == stationary_type_t::WEIGHT_STATIONARY) {
        // Case 1) Reuse weight && Request new input and output data.
        if(input_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front() - 1 &&
           output_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front() - 1) {
            // Input data and output data do not exist in local buffer.
            exist_data_lb[data_type_t::INPUT] = false, exist_data_lb[data_type_t::OUTPUT] = false;

#ifdef FUNCTIONAL
            // Write back output data
            m_scheduler->transfer_data(pe_array->output_data, output_data_lb, m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()], 0,
                                       component_type_t::PE_Y, component_type_t::PE,
                                       data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
            // Update for NPUsim ver2
            //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 &&
            //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
            //m_scheduler->transfer_data_ver2(pe_array->output_data, output_data_lb,
            //                                component_type_t::PE_Y, component_type_t::PE,
            //                                data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

            uint64_t address_pe = 0, address_pe_array = 0;
            unsigned num_access_pe = 0, num_access_pe_array = 0;
            for(unsigned b = 0; b < parameters_pe[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_pe[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_pe[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_pe[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            if(address_pe != ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                // Update read cost of Local buffer
                                access_energy_lb[data_type_t::OUTPUT] += u_read_energy_lb[data_type_t::OUTPUT];
                                access_cycle_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT];
                                num_access_pe++;

                                // Update local buffer address
                                address_pe = ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];
                            }

                            if(address_pe_array != ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT]) {

                                // Update write cost of accumulator in PE array.
                                pe_array->access_energy[data_type_t::OUTPUT] += pe_array->u_write_energy[data_type_t::OUTPUT];
                                pe_array->access_cycle[data_type_t::OUTPUT] += pe_array->u_write_cycle[data_type_t::OUTPUT];
                                num_access_pe_array++;

                                // Update accumulator address
                                address_pe_array = ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }

            /* TODO */
            // Update overlapped cycle between the local buffer and PE array.
            // Overlap read cycle of local buffer, write_back_cycle, and write cycle at the accumulator.
            const size_t output_lines_lb = runtime_datatypes().storage_transactions(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT], line_size_lb[data_type_t::OUTPUT]);
            access_cycle_lb[data_type_t::OUTPUT] += output_lines_lb*u_read_cycle_lb[data_type_t::OUTPUT];
            write_back_cycle_lb += output_lines_lb*u_read_cycle_lb[data_type_t::OUTPUT];
            access_energy_lb[data_type_t::OUTPUT] += output_lines_lb*u_read_energy_lb[data_type_t::OUTPUT];

            // Update PE array write cycle and energy.
            pe_array->access_cycle[data_type_t::OUTPUT] += runtime_datatypes().storage_transactions(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT], pe_array->line_size[data_type_t::OUTPUT])*pe_array->u_write_cycle[data_type_t::OUTPUT];
            pe_array->write_back_cycle += runtime_datatypes().storage_transactions(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT], pe_array->line_size[data_type_t::OUTPUT])*pe_array->u_write_cycle[data_type_t::OUTPUT];
            pe_array->access_energy[data_type_t::OUTPUT] += runtime_datatypes().storage_transactions(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT], pe_array->line_size[data_type_t::OUTPUT])*pe_array->u_write_energy[data_type_t::OUTPUT];
            // Update data transfer cycle between PE and PE array.

            // Increase flush counter of input data and output data.
            input_flush_counter++;
            output_flush_counter++;

            // Waiting for input and output data.
            // PE is in idle state.
            wait_data();
        }
        // Case 2) Request all data types
        else {
            // Input data, weight, and output data do not exist in local buffer.
            exist_data_lb[data_type_t::INPUT] = false, exist_data_lb[data_type_t::WEIGHT] = false, exist_data_lb[data_type_t::OUTPUT] = false;
#ifdef FUNCTIONAL
            // Write back output data
            m_scheduler->transfer_data(pe_array->output_data, output_data_lb, m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()], 0,
                                       component_type_t::PE_Y, component_type_t::PE,
                                       data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
            // Update for NPUsim ver2
            //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 &&
            //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
            //m_scheduler->transfer_data_ver2(pe_array->output_data, output_data_lb,
            //                                component_type_t::PE_Y, component_type_t::PE,
            //                                data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

            uint64_t address_pe = 0, address_pe_array = 0;
            unsigned num_access_pe = 0, num_access_pe_array = 0;
            for(unsigned b = 0; b < parameters_pe[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_pe[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_pe[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_pe[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            if(address_pe != ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                // Update read cost of Local buffer
                                access_energy_lb[data_type_t::OUTPUT] += u_read_energy_lb[data_type_t::OUTPUT];
                                access_cycle_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT];
                                num_access_pe++;

                                // Update local buffer address
                                address_pe = ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];
                            }

                            if(address_pe_array != ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT]) {

                                // Update write cost of accumulator in PE array.
                                pe_array->access_energy[data_type_t::OUTPUT] += pe_array->u_write_energy[data_type_t::OUTPUT];
                                pe_array->access_cycle[data_type_t::OUTPUT] += pe_array->u_write_cycle[data_type_t::OUTPUT];
                                num_access_pe_array++;

                                // Update accumulator address
                                address_pe_array = ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }

            /*
            // Update stats.
            // Update local buffer read cycle and energy.
            const size_t output_lines_lb = runtime_datatypes().storage_transactions(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT], line_size_lb[data_type_t::OUTPUT]);
            access_cycle_lb[data_type_t::OUTPUT] += output_lines_lb*u_read_cycle_lb[data_type_t::OUTPUT];
            write_back_cycle_lb += output_lines_lb*u_read_cycle_lb[data_type_t::OUTPUT];
            access_energy_lb[data_type_t::OUTPUT] += output_lines_lb*u_read_energy_lb[data_type_t::OUTPUT];

            // Update PE array write cycle and energy.
            pe_array->access_cycle[data_type_t::OUTPUT] += runtime_datatypes().storage_transactions(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT], pe_array->line_size[data_type_t::OUTPUT])*pe_array->u_write_cycle[data_type_t::OUTPUT];
            pe_array->write_back_cycle += runtime_datatypes().storage_transactions(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT], pe_array->line_size[data_type_t::OUTPUT])*pe_array->u_write_cycle[data_type_t::OUTPUT];
            pe_array->access_energy[data_type_t::OUTPUT] += runtime_datatypes().storage_transactions(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT], pe_array->line_size[data_type_t::OUTPUT])*pe_array->u_write_energy[data_type_t::OUTPUT];
            */

            // Set flush counter of input data and output data as zero.
            input_flush_counter = 0;
            output_flush_counter = 0;

            // Waiting for input data, weight, and output data.
            // PE is in idle state.
            wait_data();
        }
    }
    // Case 3. Output stationary
    else if(stationary_type_local_buffer == stationary_type_t::OUTPUT_STATIONARY) {
        // Case 1) Reuse output data && Request input data and weight
        if(input_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front() - 1 &&
           weight_flush_counter < m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front() - 1) {
            // input data and weight do not exist in local buffer.
            exist_data_lb[data_type_t::INPUT] = false, exist_data_lb[data_type_t::WEIGHT] = false;

            // Increase flush counter of input data and weight.
            input_flush_counter++;
            weight_flush_counter++;

            // Waiting for input data and weight.
            // PE is in idle state
            wait_data();
        }
        // Case 2) Request all data types
        else {
            // Input data, weight, and output data are not in local buffer.
            exist_data_lb[data_type_t::INPUT] = false, exist_data_lb[data_type_t::WEIGHT] = false, exist_data_lb[data_type_t::OUTPUT] = false;
#ifdef FUNCTIONAL
            // Write back output data
            m_scheduler->transfer_data(pe_array->output_data, output_data_lb, m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()], 0,
                                       component_type_t::PE_Y, component_type_t::PE,
                                       data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
            // Update for NPUsim ver2
            //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 &&
            //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
            //m_scheduler->transfer_data_ver2(pe_array->output_data, output_data_lb,
            //                                component_type_t::PE_Y, component_type_t::PE,
            //                                data_type_t::OUTPUT, get_local_buffer_stationary_type(), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_pe(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);
            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);

            uint64_t address_pe = 0, address_pe_array = 0;
            unsigned num_access_pe = 0, num_access_pe_array = 0;
            for(unsigned b = 0; b < parameters_pe[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_pe[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_pe[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_pe[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            if(address_pe != ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                // Update read cost of Local buffer
                                access_energy_lb[data_type_t::OUTPUT] += u_read_energy_lb[data_type_t::OUTPUT];
                                access_cycle_lb[data_type_t::OUTPUT] += u_read_cycle_lb[data_type_t::OUTPUT];
                                num_access_pe++;

                                // Update local buffer address
                                address_pe = ((uint64_t)&output_data_lb[b*parameters_pe[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_pe[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_pe[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_pe[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];
                            }

                            if(address_pe_array != ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT]) {

                                // Update write cost of accumulator in PE array.
                                pe_array->access_energy[data_type_t::OUTPUT] += pe_array->u_write_energy[data_type_t::OUTPUT];
                                pe_array->access_cycle[data_type_t::OUTPUT] += pe_array->u_write_cycle[data_type_t::OUTPUT];
                                num_access_pe_array++;

                                // Update accumulator address
                                address_pe_array = ((uint64_t)&pe_array->output_data[m_scheduler->output_offset_pe_array[index%m_scheduler->output_offset_pe_array.size()] +
                                                                                     b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                      *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] +
                                                                                     p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                                     pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT];
                            }
                        }
                    }
                }
            }

            /*
            // Update local buffer read cycle and energy.
            const size_t output_lines_lb = runtime_datatypes().storage_transactions(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT], line_size_lb[data_type_t::OUTPUT]);
            access_cycle_lb[data_type_t::OUTPUT] += output_lines_lb*u_read_cycle_lb[data_type_t::OUTPUT];
            write_back_cycle_lb += output_lines_lb*u_read_cycle_lb[data_type_t::OUTPUT];
            access_energy_lb[data_type_t::OUTPUT] += output_lines_lb*u_read_energy_lb[data_type_t::OUTPUT];

            // Update PE array write cycle and energy.
            pe_array->access_cycle[data_type_t::OUTPUT] += runtime_datatypes().storage_transactions(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT], pe_array->line_size[data_type_t::OUTPUT])*pe_array->u_write_cycle[data_type_t::OUTPUT];
            pe_array->write_back_cycle += runtime_datatypes().storage_transactions(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT], pe_array->line_size[data_type_t::OUTPUT])*pe_array->u_write_cycle[data_type_t::OUTPUT];
            pe_array->access_energy[data_type_t::OUTPUT] += runtime_datatypes().storage_transactions(data_type_t::OUTPUT, tile_size_lb[data_type_t::OUTPUT], pe_array->line_size[data_type_t::OUTPUT])*pe_array->u_write_energy[data_type_t::OUTPUT];
            */

            // Set flush counter of input data and weight to zero.
            input_flush_counter = 0;
            weight_flush_counter = 0;

            // Waiting for input data, weight, and output data.
            // PE is in idle state.
            wait_data();
        }
    }
    // The counter of local buffer should be initialized to 0.
    input_index = 0, weight_index = 0, output_index = 0;
    //if(index == pe_array->get_number_of_active_pes()-1) {
    //    pe_array->update_pe_stats();
    //}
}

void pe_t::clear_output_accumulators() {
    std::fill_n(output_data_mac, mac_register_capacity, data_t{});
}

unsigned pe_t::count_nonzero_mac_operations(scheduler_t *m_scheduler) const {
    const mac_geometry_t geometry = get_mac_geometry(m_scheduler);
    validate_mac_geometry(geometry, mac_register_capacity, num_active_macs);

    unsigned nonzero_operations = 0;
    for(unsigned b = 0; b < geometry.batch; b++) {
        for(unsigned g = 0; g < geometry.groups; g++) {
            for(unsigned k = 0; k < geometry.output_channels_per_group; k++) {
                for(unsigned p = 0; p < geometry.output_height; p++) {
                    for(unsigned q = 0; q < geometry.output_width; q++) {
                        for(unsigned c = 0; c < geometry.input_channels_per_group; c++) {
                            for(unsigned r = 0; r < geometry.filter_height; r++) {
                                for(unsigned s = 0; s < geometry.filter_width; s++) {
                                    const size_t input_index = ((((static_cast<size_t>(b)*geometry.groups + g)*geometry.input_channels_per_group + c)*geometry.input_height + p*geometry.stride + r)*geometry.input_width + q*geometry.stride + s);
                                    const size_t weight_index = ((((static_cast<size_t>(g)*geometry.output_channels_per_group + k)*geometry.input_channels_per_group + c)*geometry.filter_height + r)*geometry.filter_width + s);
                                    if(data_is_nonzero(input_data_mac[input_index]) &&
                                       data_is_nonzero(weight_mac[weight_index])) {
                                        nonzero_operations++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return nonzero_operations;
}

void pe_t::mac_operation(scheduler_t *m_scheduler) {
    const mac_geometry_t geometry = get_mac_geometry(m_scheduler);
    validate_mac_geometry(geometry, mac_register_capacity, num_active_macs);

    for(unsigned b = 0; b < geometry.batch; b++) {
        for(unsigned g = 0; g < geometry.groups; g++) {
            for(unsigned k = 0; k < geometry.output_channels_per_group; k++) {
                for(unsigned p = 0; p < geometry.output_height; p++) {
                    for(unsigned q = 0; q < geometry.output_width; q++) {
                        const size_t output_index = ((((static_cast<size_t>(b)*geometry.groups + g)*geometry.output_channels_per_group + k)*geometry.output_height + p)*geometry.output_width + q);
                        for(unsigned c = 0; c < geometry.input_channels_per_group; c++) {
                            for(unsigned r = 0; r < geometry.filter_height; r++) {
                                for(unsigned s = 0; s < geometry.filter_width; s++) {
                                    const size_t input_index = ((((static_cast<size_t>(b)*geometry.groups + g)*geometry.input_channels_per_group + c)*geometry.input_height + p*geometry.stride + r)*geometry.input_width + q*geometry.stride + s);
                                    const size_t weight_index = ((((static_cast<size_t>(g)*geometry.output_channels_per_group + k)*geometry.input_channels_per_group + c)*geometry.filter_height + r)*geometry.filter_width + s);
                                    data_accumulate_product(output_data_mac[output_index],
                                                            input_data_mac[input_index],
                                                            weight_mac[weight_index]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void pe_t::activation() {
    for(size_t i = 0; i < mac_register_capacity; i++) {
        data_relu(output_data_mac[i]);
    }
}


// Print out the stats of the PE.
void pe_t::print_specification() {

    std::cout << std::fixed;
    std::cout << "============ MAC specification =============" << std::endl;
    std::cout << "# of MACs          :" << std::setw(24)
                                        << num_macs*mac_width << std::endl;
    std::cout << "Computation energy :" << std::setw(21) << std::setprecision(2)
                                        << u_computation_energy << " " << energy_units().label() << std::endl;
    std::cout << "Computation cycle  :" << std::setw(17) << std::setprecision(1)
                                        << u_computation_cycle << " cycles" << std::endl;

    std::cout << "Line size" << std::endl;
    std::cout << " * Input           :" << std::setw(19) << std::setprecision(0)
                                        << line_size_mac[data_type_t::INPUT] << " bits" << std::endl;
    std::cout << " * Weight          :" << std::setw(19) << std::setprecision(0)
                                        << line_size_mac[data_type_t::WEIGHT] << " bits" << std::endl;
    std::cout << " * Output          :" << std::setw(19) << std::setprecision(0)
                                        << line_size_mac[data_type_t::OUTPUT] << " bits" << std::endl;

    std::cout << "Access energy (read/write)" << std::endl;
    std::cout << " * MAC units       :" << std::setw(16) << std::setprecision(2)
                                        << u_read_energy_mac[data_type_t::INPUT] << "/" << u_write_energy_mac[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    std::cout << "Access cycle (read/write)" << std::endl;
    std::cout << " * MAC units       :" << std::setw(13) << std::setprecision(1)
                                        << u_read_cycle_mac[data_type_t::INPUT] << "/" << u_write_cycle_mac[data_type_t::INPUT] << " cycles" << std::endl;
    std::cout << std::endl;

    std::cout << "============= PE specification =============" << std::endl;
    // Print the size of local buffer.
    if(memory_type == memory_type_t::SEPARATE) {
        std::cout << "Local buffer type  :" << std::setw(24)
                                            << "Separated" << std::endl;
        std::cout << " * Input           :" << std::setw(19)
                                            << input_size << " Byte" << std::endl;
        std::cout << " * Weight          :" << std::setw(19)
                                            << weight_size << " Byte" << std::endl;
        std::cout << " * Output          :" << std::setw(19)
                                            << output_size << " Byte" << std::endl;

        std::cout << "Line size" << std::endl;
        std::cout << " * Input           :" << std::setw(19) << std::setprecision(0)
                                            << line_size_lb[data_type_t::INPUT] << " bits" << std::endl;
        std::cout << " * Weight          :" << std::setw(19) << std::setprecision(0)
                                            << line_size_lb[data_type_t::WEIGHT] << " bits" << std::endl;
        std::cout << " * Output          :" << std::setw(19) << std::setprecision(0)
                                            << line_size_lb[data_type_t::OUTPUT] << " bits" << std::endl;
    }
    else if(memory_type == memory_type_t::SHARED) {
        std::cout << "Local buffer type  :" << std::setw(24)
                                            << "Shared" << std::endl;
        std::cout << "Buffer size        :" << std::setw(19)
                                            << input_size + weight_size + output_size << " Byte" << std::endl;
        std::cout << "Line size          :" << std::setw(19) << std::setprecision(0)
                                            << line_size_lb[data_type_t::INPUT] << " bits" << std::endl;
    }
    // Print the stationary type.
    if(stationary_type_local_buffer == stationary_type_t::INPUT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24)
                                            << "Input stationary" << std::endl;
    }
    else if(stationary_type_local_buffer == stationary_type_t::WEIGHT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24)
                                            << "Weight stationary" << std::endl;
    }
    else if(stationary_type_local_buffer == stationary_type_t::OUTPUT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24)
                                            << "Output stationary" << std::endl;
    }
    // Print out the bandwidth
    std::cout << "Bandwidth          :" << std::setw(19) << std::setprecision(0)
                                        << bandwidth << " GB/s" << std::endl;
    // Print out unit access cost of the local buffer.
    if(memory_type == memory_type_t::SEPARATE) {
        std::cout << "Access energy (read/write)" << std::endl;
        std::cout << " * Input buffer    :" << std::setw(16) << std::setprecision(2)
                                            << u_read_energy_lb[data_type_t::INPUT]  << "/" << u_write_energy_lb[data_type_t::INPUT]  << " " << energy_units().label() << std::endl;
        std::cout << " * Weight buffer   :" << std::setw(16) << std::setprecision(2)
                                            << u_read_energy_lb[data_type_t::WEIGHT] << "/" << u_write_energy_lb[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
        std::cout << " * Output buffer   :" << std::setw(16) << std::setprecision(2)
                                            << u_read_energy_lb[data_type_t::OUTPUT] << "/" << u_write_energy_lb[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;

        std::cout << "Access cycle (read/write)" << std::endl;
        std::cout << " * Input buffer    :" << std::setw(13) << std::setprecision(1)
                                            << u_read_cycle_lb[data_type_t::INPUT] << "/" << u_write_cycle_lb[data_type_t::INPUT] << " cycles" << std::endl;
        std::cout << " * Weight buffer   :" << std::setw(13) << std::setprecision(1)
                                            << u_read_cycle_lb[data_type_t::WEIGHT] << "/" << u_write_cycle_lb[data_type_t::WEIGHT] << " cycles" << std::endl;
        std::cout << " * Output buffer   :" << std::setw(13) << std::setprecision(1)
                                            << u_read_cycle_lb[data_type_t::OUTPUT] << "/" << u_write_cycle_lb[data_type_t::OUTPUT] << " cycles" << std::endl;
    }
    else if(memory_type == memory_type_t::SHARED) {
        std::cout << "Access energy (read/write)" << std::endl;
        std::cout << " * Local buffer    :" << std::setw(16) << std::setprecision(2)
                                            << u_read_energy_lb[data_type_t::INPUT]  << "/" << u_write_energy_lb[data_type_t::INPUT]  << " " << energy_units().label() << std::endl;

        std::cout << "Access cycle (read/write)" << std::endl;
        std::cout << " * Local buffer    :" << std::setw(13) << std::setprecision(1)
                                            << u_read_cycle_lb[data_type_t::INPUT] << "/" << u_write_cycle_lb[data_type_t::INPUT] << " cycles" << std::endl;
    }
	std::cout << std::endl;
}

void pe_t::account_format_events(data_type_t type, size_t elements) {
    const size_t payload = runtime_datatypes().payload_transactions(type, elements, bitwidth);
    const size_t metadata = runtime_datatypes().metadata_transactions(type, elements, bitwidth);

    format_payload_events[type] += payload;
    format_metadata_events[type] += metadata;
    format_cycle[type] += payload*u_format_payload_cycle[type] +
                          metadata*u_format_metadata_cycle[type];
    format_energy[type] += payload*u_format_payload_energy[type] +
                           metadata*u_format_metadata_energy[type];
}

// RE1: the four output-accumulator events, each generated at its own boundary.
//
// All of them are charged per BYTE MOVED rather than per link transaction: a link-transaction
// count cannot express precision at this granularity (the MAC tile is a handful of elements, so
// ceil(elements * width / link_bits) rounds to the same 1 for fp32 and fp16 alike), and these are
// precision conversions on values, not link crossings.
//
// OWNERSHIP: with edge_accumulation the accumulator physically lives at the array edge, so its
// energy is handed to the PE array. Charging it to the PE would put the energy in a component that
// the config says does not hold the accumulator.
void pe_t::charge_accumulator_bytes(size_t bytes) {
    const double energy = bytes*u_accumulator_spill_energy[data_type_t::OUTPUT];
    format_cycle[data_type_t::OUTPUT] += bytes*u_accumulator_spill_cycle[data_type_t::OUTPUT];
    // RE1 ownership: with edge_accumulation the accumulator physically sits at the array edge, so
    // its energy belongs to the PE array. Charging it to the PE would attribute energy to a
    // component the config says does not hold the accumulator.
    if(edge_accumulation) {
        pe_array->accumulator_energy += energy;
    } else {
        accumulator_energy += energy;
    }
}

void pe_t::account_accumulator_create(size_t elements) {
    // A fresh accumulator is zero-initialized in place: nothing is read back and nothing is
    // spilled, so there is no reload/spill energy to charge. Counted so the report can show that
    // these requests exist and are deliberately free.
    (void)elements;
    ++accumulator_create_events;
}

void pe_t::account_accumulator_reload(size_t elements) {
    const size_t bytes = runtime_datatypes().accumulator_storage_bytes(elements);
    accumulator_reload_bytes += bytes;
    charge_accumulator_bytes(bytes);
}

void pe_t::account_accumulator_spill(size_t elements) {
    const size_t bytes = runtime_datatypes().accumulator_storage_bytes(elements);
    accumulator_spill_bytes += bytes;
    charge_accumulator_bytes(bytes);
}

double pe_t::modeled_elapsed_cycles() const {
    double elapsed = std::max(computation_cycle,
                              std::max(write_back_cycle_mac, write_back_cycle_lb));
    elapsed = std::max(elapsed, overlapped_transfer_cycle);
    for(unsigned type = 0; type < data_type_t::NUM_DATA_TYPES; ++type) {
        elapsed = std::max(elapsed, access_cycle_mac[type]);
        elapsed = std::max(elapsed, access_cycle_lb[type]);
        elapsed = std::max(elapsed, transfer_cycle[type]);
        elapsed = std::max(elapsed, cycle_mac_lb[type]);
        elapsed = std::max(elapsed, format_cycle[type]);
    }
    return elapsed;
}

size_t pe_t::scalar_mac_capacity() const {
    return mac_register_capacity;
}

void pe_t::update_static_energy(double elapsed_cycles) {
    if(elapsed_cycles < 0.0) {
        std::cerr << "Error: PE elapsed cycles cannot be negative" << std::endl;
        exit(1);
    }
    for(unsigned type = 0; type < data_type_t::NUM_DATA_TYPES; ++type) {
        if(u_static_energy[type] < 0.0 || u_lb_static_energy[type] < 0.0) {
            std::cerr << "Error: PE static_energy/lb_static_energy must be non-negative pJ/cycle values" << std::endl;
            exit(1);
        }
        // LB4: MAC-side and local-buffer leakage accrue over the same elapsed window.
        static_energy[type] = (u_static_energy[type] + u_lb_static_energy[type]) * elapsed_cycles;
    }
}

// Reset the stats
void pe_t::reset() {
    std::fill_n(input_data_mac, mac_register_capacity, data_t{});
    std::fill_n(weight_mac, mac_register_capacity, data_t{});
    clear_output_accumulators();

    std::fill_n(input_data_lb, (static_cast<size_t>(input_size) + sizeof(data_t) - 1)/sizeof(data_t), data_t{});
    std::fill_n(weight_lb, (static_cast<size_t>(weight_size) + sizeof(data_t) - 1)/sizeof(data_t), data_t{});
    std::fill_n(output_data_lb, (static_cast<size_t>(output_size) + sizeof(data_t) - 1)/sizeof(data_t), data_t{});

    input_index = 0;
    weight_index = 0;
    output_index = 0;
    input_flush_counter = 0;
    weight_flush_counter = 0;
    output_flush_counter = 0;
    exist_data_mac.assign(data_type_t::NUM_DATA_TYPES, false);
    exist_data_lb.assign(data_type_t::NUM_DATA_TYPES, false);
    request_to_lb.assign(data_type_t::NUM_DATA_TYPES, true);
    request_to_pe_array.assign(data_type_t::NUM_DATA_TYPES, false);
    offsets_mac.assign(data_type_t::NUM_DATA_TYPES, 0);
    offsets_lb.assign(data_type_t::NUM_DATA_TYPES, 0);
    idle = false;

    num_computation = 0;

    computation_cycle = 0;
    computation_energy = 0;
    // E4: the reduction-energy axis and the two format event counters reset with the layer.
    reduction_energy = 0.0;
    accumulator_reload_bytes = 0;
    accumulator_spill_bytes = 0;
    accumulator_create_events = 0;
    accumulator_retained_events = 0;
    accumulator_energy = 0.0;
    format_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    format_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    write_back_cycle_mac  = 0.0;
    write_back_cycle_lb = 0.0;
    overlapped_transfer_cycle = 0.0;

    utilization_mac = 0.0;

    num_request_to_lb.assign(data_type_t::NUM_DATA_TYPES, 0);
    num_data_transfer_to_mac.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Reset access cycle and energy of MAC unit.
    access_cycle_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    access_energy_mac.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Reset access cycle and energy of local buffer.
    access_cycle_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    access_energy_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    payload_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    metadata_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    storage_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Reset overlapped cycle between MAC unit and local buffer
    cycle_mac_lb.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Reset static energy of PE.
    static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    utilization_local_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);

}

// All functions of undefined stationary
undefined_stationary_t::undefined_stationary_t(section_config_t m_section_config) :
    pe_t(m_section_config) {

}

undefined_stationary_t::~undefined_stationary_t() {

}

void undefined_stationary_t::computation(scheduler_t *m_scheduler) {

    if(exist_data_mac[data_type_t::INPUT] && exist_data_mac[data_type_t::WEIGHT] && exist_data_mac[data_type_t::OUTPUT]) {
#ifdef FUNCTIONAL

#endif

        // PE7: scale computation count/energy by the active scalar-MAC lanes, matching the
        // input/weight/output_stationary paths (energy per active MAC; one issue step/tile).
        for(unsigned i = 0; i < num_active_macs; i++) {
            num_computation++;
            computation_energy += u_computation_energy;
        }
        computation_cycle += accumulate_issue_cycles(1, u_computation_cycle) +
                     lane_reduction_fill_cycles(lane_state, u_computation_cycle);
                     reduction_energy += lane_reduction_energy(lane_state, u_mac_reduction_energy);

        exist_data_mac[data_type_t::INPUT] = false, request_to_lb[data_type_t::INPUT] = true;
        exist_data_mac[data_type_t::WEIGHT] = false, request_to_lb[data_type_t::WEIGHT] = true;
        exist_data_mac[data_type_t::OUTPUT] = false, request_to_lb[data_type_t::OUTPUT] = true;

        // Update stats (Memory operation)
        /* Stats */
        num_request_to_lb[data_type_t::INPUT]++;
        num_request_to_lb[data_type_t::WEIGHT]++;
        num_request_to_lb[data_type_t::OUTPUT]++;
        /* Stats */

#ifdef FUNCTIONAL
        // Write back output data
        m_scheduler->transfer_data(output_data_lb, output_data_mac, m_scheduler->output_offset_pe.front(), 0,
                                   component_type_t::PE, component_type_t::MAC,
                                   data_type_t::OUTPUT, get_mac_stationary_type(), action_type_t::STORE);
        // Update for NPUsim ver2
        //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 &&
        //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
        //m_scheduler->transfer_data_ver2(output_data_lb, output_data_mac,
        //                                component_type_t::PE, component_type_t::MAC,
        //                                data_type_t::OUTPUT, get_mac_stationary_type()), action_type_t::STORE);
#endif
        /* Stats */
        account_descriptor_dense_mac_transfer(
            data_type_t::OUTPUT, tile_size_mac[data_type_t::OUTPUT], false);
        /* Stats */


        if(input_index == m_scheduler->offset_size_pe[data_type_t::INPUT].front() &&
           weight_index == m_scheduler->offset_size_pe[data_type_t::WEIGHT].front()) {

            input_index = 0, weight_index = 0;
        }
    }
}


input_stationary_t::input_stationary_t(section_config_t m_section_config) :
    pe_t(m_section_config) {

}

input_stationary_t::~input_stationary_t() {

}

// Execute MAC operation
void input_stationary_t::computation(scheduler_t *m_scheduler) {

    // When all the data exist, execute the MAC operation
    if(exist_data_mac[data_type_t::INPUT] && exist_data_mac[data_type_t::WEIGHT] && exist_data_mac[data_type_t::OUTPUT]) {
        if(m_scheduler->layer_name == layer_name_t::CONVOLUTIONAL_LAYER ||
           m_scheduler->layer_name == layer_name_t::CONNECTED_LAYER) {
#ifdef FUNCTIONAL
            mac_operation(m_scheduler);
            // Activation is applied only after a completed output reduction.
            // Split active_macs and mac_width
            if(m_scheduler->compression_type == compression_type_t::DENSE) {
                for(unsigned i = 0; i < num_active_macs; i++) {
                    num_computation++;
                    computation_energy += u_computation_energy;
                }
                computation_cycle += accumulate_issue_cycles(1, u_computation_cycle) +
                     lane_reduction_fill_cycles(lane_state, u_computation_cycle);
                     reduction_energy += lane_reduction_energy(lane_state, u_mac_reduction_energy);
            }
            else {
                const unsigned nonzero_operations = count_nonzero_mac_operations(m_scheduler);
                num_computation += nonzero_operations;
                computation_energy += nonzero_operations*u_computation_energy;
                if(nonzero_operations > 0) {
                    computation_cycle += accumulate_issue_cycles(1, u_computation_cycle) +
                     lane_reduction_fill_cycles(lane_state, u_computation_cycle);
                     reduction_energy += lane_reduction_energy(lane_state, u_mac_reduction_energy);
                }
            }
#else
            for(unsigned i = 0; i < num_active_macs; i++) {
                num_computation++;
                computation_energy += u_computation_energy;
            }
            computation_cycle += accumulate_issue_cycles(1, u_computation_cycle) +
                     lane_reduction_fill_cycles(lane_state, u_computation_cycle);
                     reduction_energy += lane_reduction_energy(lane_state, u_mac_reduction_energy);
#endif
        } else {
            std::cerr << "Error: PE computation supports only convolution/connected layers" << std::endl;
            exit(1);
        }
#ifndef FUNCTIONAL
        // RE1: a MAC -> LB write-back of a partial sum IS the accumulator spill.
        account_accumulator_spill(tile_size_mac[data_type_t::OUTPUT]);
        account_descriptor_dense_mac_transfer(data_type_t::OUTPUT, tile_size_mac[data_type_t::OUTPUT], false);
        exist_data_mac[data_type_t::WEIGHT] = false;
        request_to_lb[data_type_t::WEIGHT] = true;
        exist_data_mac[data_type_t::OUTPUT] = false;
        request_to_lb[data_type_t::OUTPUT] = true;
        num_request_to_lb[data_type_t::WEIGHT]++;
        num_request_to_lb[data_type_t::OUTPUT]++;
        if(weight_index == m_scheduler->offset_size_pe[data_type_t::WEIGHT].front() &&
           output_index == m_scheduler->offset_size_pe[data_type_t::OUTPUT].front()) {
            exist_data_mac[data_type_t::INPUT] = false;
            request_to_lb[data_type_t::INPUT] = true;
            num_request_to_lb[data_type_t::INPUT]++;
            if(input_index < m_scheduler->offset_size_pe[data_type_t::INPUT].front()) {
                weight_index = 0;
                output_index = 0;
            }
        }
        return;
#endif
#ifdef FUNCTIONAL
        // Write back output data
        m_scheduler->transfer_data(output_data_lb, output_data_mac, m_scheduler->output_offset_pe.front(), 0,
                                   component_type_t::PE, component_type_t::MAC,
                                   data_type_t::OUTPUT, get_mac_stationary_type(), action_type_t::STORE);
        // Update for NPUsim ver2
        //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 &&
        //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
        //m_scheduler->transfer_data_ver2(output_data_lb, output_data_mac,
        //                                component_type_t::PE, component_type_t::MAC,
        //                                data_type_t::OUTPUT, get_mac_stationary_type()), action_type_t::STORE);
#endif

        std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
        std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);

        parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
        parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

        uint64_t address_mac = 0, address_lb = 0;
        unsigned num_access_mac = 0, num_access_lb = 0;
        for(unsigned b = 0; b < parameters_mac[parameter_type_t::BATCH_SIZE]; b++) {
            for(unsigned k = 0; k < parameters_mac[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                for(unsigned p = 0; p < parameters_mac[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                    for(unsigned q = 0; q < parameters_mac[parameter_type_t::OUTPUT_WIDTH]; q++) {
                        if(address_mac != ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                      mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT]) {

                            // Update write cost of MAC unit.
                            access_energy_mac[data_type_t::OUTPUT] += u_read_energy_mac[data_type_t::OUTPUT];
                            access_cycle_mac[data_type_t::OUTPUT] += u_read_cycle_mac[data_type_t::OUTPUT];
                            num_access_mac++;

                            // Update MAC address
                            address_mac = ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                      mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT];
                        }

                        if(address_lb != ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                    b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                    mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                            // Update read cost of Local buffer
                            access_energy_lb[data_type_t::OUTPUT] += u_write_energy_lb[data_type_t::OUTPUT];
                            access_cycle_lb[data_type_t::OUTPUT] += u_write_cycle_lb[data_type_t::OUTPUT];
                            num_access_lb++;

                            // Update local buffer address
                            address_lb = ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                    b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                    mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];

                        }
                    }
                }
            }
        }

        // Update overlapped cycle at MAC units and local buffer
        unsigned ratio = ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(line_size_mac[data_type_t::OUTPUT]));

        double first_stage = ratio*u_read_cycle_mac[data_type_t::OUTPUT];
        double second_stage = std::max(ratio*u_read_cycle_mac[data_type_t::OUTPUT],
                                         u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
        double last_before_stage = std::max(u_write_cycle_lb[data_type_t::OUTPUT],
                                              u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
        double last_stage = u_write_cycle_lb[data_type_t::OUTPUT];

        double other_stage = std::max(ratio*u_read_cycle_mac[data_type_t::OUTPUT],
                               std::max(u_write_cycle_lb[data_type_t::OUTPUT],
                                        u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth))));

        if(num_access_lb == 0) {
            // No transaction: leave the overlapped stage unchanged.
        } else if(num_access_lb == 1) {
            cycle_mac_lb[data_type_t::OUTPUT] += ratio*u_read_cycle_mac[data_type_t::OUTPUT] +
                                                 u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)) +
                                                 u_write_cycle_lb[data_type_t::OUTPUT];
        } else if(num_access_lb == 2) {
            cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage + last_before_stage + last_stage;
        } else {
            cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage +
                                                 (num_access_lb - 2)*other_stage +
                                                 last_before_stage + last_stage;
        }

        transfer_energy[data_type_t::OUTPUT] += num_access_mac*u_transfer_energy*ceil((double)(line_size_mac[data_type_t::OUTPUT])/(double)(bitwidth));
        transfer_cycle[data_type_t::OUTPUT] += num_access_mac*u_transfer_cycle*ceil((double)(line_size_mac[data_type_t::OUTPUT])/(double)(bitwidth));

        // weight and output data should be requested from MAC to local buffer.
        exist_data_mac[data_type_t::WEIGHT] = false, request_to_lb[data_type_t::WEIGHT] = true;
        exist_data_mac[data_type_t::OUTPUT] = false, request_to_lb[data_type_t::OUTPUT] = true;

        // Update stats (Memory operation)
        num_request_to_lb[data_type_t::WEIGHT]++;
        num_request_to_lb[data_type_t::OUTPUT]++;

        // When all the weight and output data at local buffer are used.
        // A new input tile should be requested from MAC unit to local buffer.
        if(weight_index == m_scheduler->offset_size_pe[data_type_t::WEIGHT].front() &&
           output_index == m_scheduler->offset_size_pe[data_type_t::OUTPUT].front()) {

            // Request new input data.
            exist_data_mac[data_type_t::INPUT] = false, request_to_lb[data_type_t::INPUT] = true;

            // Update stats.
            num_request_to_lb[data_type_t::INPUT]++;

            // When not all input data in local buffer is transferred to MAC unit.
            if(input_index < m_scheduler->offset_size_pe[data_type_t::INPUT].front()) {
                weight_index = 0, output_index = 0;
            }
        }
    }
}

weight_stationary_t::weight_stationary_t(section_config_t m_section_config) :
    pe_t(m_section_config) {

}

weight_stationary_t::~weight_stationary_t() {

}

void weight_stationary_t::computation(scheduler_t *m_scheduler) {

    // Execute MAC operation when all DNN data exist in MAC register.
    if(exist_data_mac[data_type_t::INPUT] && exist_data_mac[data_type_t::WEIGHT] && exist_data_mac[data_type_t::OUTPUT]) {
        if(m_scheduler->layer_name == layer_name_t::CONVOLUTIONAL_LAYER ||
           m_scheduler->layer_name == layer_name_t::CONNECTED_LAYER) {
#ifdef FUNCTIONAL
            mac_operation(m_scheduler);
            // Activation is applied only after a completed output reduction.
            // Split active_macs and mac_width
            if(m_scheduler->compression_type == compression_type_t::DENSE) {
                for(unsigned i = 0; i < num_active_macs; i++) {
                    num_computation++;
                    computation_energy += u_computation_energy;
                }
                computation_cycle += accumulate_issue_cycles(1, u_computation_cycle) +
                     lane_reduction_fill_cycles(lane_state, u_computation_cycle);
                     reduction_energy += lane_reduction_energy(lane_state, u_mac_reduction_energy);
            }
            else {
                const unsigned nonzero_operations = count_nonzero_mac_operations(m_scheduler);
                num_computation += nonzero_operations;
                computation_energy += nonzero_operations*u_computation_energy;
                if(nonzero_operations > 0) {
                    computation_cycle += accumulate_issue_cycles(1, u_computation_cycle) +
                     lane_reduction_fill_cycles(lane_state, u_computation_cycle);
                     reduction_energy += lane_reduction_energy(lane_state, u_mac_reduction_energy);
                }
            }
#else
            for(unsigned i = 0; i < num_active_macs; i++) {
                num_computation++;
                computation_energy += u_computation_energy;
            }
            computation_cycle += accumulate_issue_cycles(1, u_computation_cycle) +
                     lane_reduction_fill_cycles(lane_state, u_computation_cycle);
                     reduction_energy += lane_reduction_energy(lane_state, u_mac_reduction_energy);
#endif
        }
        else {
            std::cerr << "Error: PE computation supports only convolution/connected layers" << std::endl;
            exit(1);
        }
#ifndef FUNCTIONAL
        // RE1: a MAC -> LB write-back of a partial sum IS the accumulator spill.
        account_accumulator_spill(tile_size_mac[data_type_t::OUTPUT]);
        account_descriptor_dense_mac_transfer(data_type_t::OUTPUT, tile_size_mac[data_type_t::OUTPUT], false);
        exist_data_mac[data_type_t::INPUT] = false;
        request_to_lb[data_type_t::INPUT] = true;
        exist_data_mac[data_type_t::OUTPUT] = false;
        request_to_lb[data_type_t::OUTPUT] = true;
        num_request_to_lb[data_type_t::INPUT]++;
        num_request_to_lb[data_type_t::OUTPUT]++;
        if(input_index == m_scheduler->offset_size_pe[data_type_t::INPUT].front() &&
           output_index == m_scheduler->offset_size_pe[data_type_t::OUTPUT].front()) {
            exist_data_mac[data_type_t::WEIGHT] = false;
            request_to_lb[data_type_t::WEIGHT] = true;
            num_request_to_lb[data_type_t::WEIGHT]++;
            if(weight_index < m_scheduler->offset_size_pe[data_type_t::WEIGHT].front()) {
                input_index = 0;
                output_index = 0;
            }
        }
        return;
#endif

#ifdef FUNCTIONAL
        // Write back output data
        m_scheduler->transfer_data(output_data_lb, output_data_mac, m_scheduler->output_offset_pe.front(), 0,
                                   component_type_t::PE, component_type_t::MAC,
                                   data_type_t::OUTPUT, get_mac_stationary_type(), action_type_t::STORE);
        // Update for NPUsim ver2
        //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 &&
        //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
        //m_scheduler->transfer_data_ver2(output_data_lb, output_data_mac,
        //                                component_type_t::PE, component_type_t::MAC,
        //                                data_type_t::OUTPUT, get_mac_stationary_type()), action_type_t::STORE);
#endif
        std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
        std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);

        parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
        parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

        uint64_t address_mac = 0, address_lb = 0;
        unsigned num_access_mac = 0, num_access_lb = 0;
        for(unsigned b = 0; b < parameters_mac[parameter_type_t::BATCH_SIZE]; b++) {
            for(unsigned k = 0; k < parameters_mac[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                for(unsigned p = 0; p < parameters_mac[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                    for(unsigned q = 0; q < parameters_mac[parameter_type_t::OUTPUT_WIDTH]; q++) {
                        if(address_mac != ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                      mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT]) {

                            // Update write cost of MAC unit.
                            access_energy_mac[data_type_t::OUTPUT] += u_read_energy_mac[data_type_t::OUTPUT];
                            access_cycle_mac[data_type_t::OUTPUT] += u_read_cycle_mac[data_type_t::OUTPUT];
                            num_access_mac++;

                            // Update MAC address
                            address_mac = ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                       *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                      p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                      mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT];
                        }

                        if(address_lb != ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                    b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                    mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                            // Update read cost of Local buffer
                            access_energy_lb[data_type_t::OUTPUT] += u_write_energy_lb[data_type_t::OUTPUT];
                            access_cycle_lb[data_type_t::OUTPUT] += u_write_cycle_lb[data_type_t::OUTPUT];
                            num_access_lb++;

                            // Update local buffer address
                            address_lb = ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                    b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                     *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                    p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                    mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];

                        }
                    }
                }
            }
        }

        // Update overlapped cycle at MAC units and local buffer
        unsigned ratio = ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(line_size_mac[data_type_t::OUTPUT]));

        double first_stage = ratio*u_read_cycle_mac[data_type_t::OUTPUT];
        double second_stage = std::max(ratio*u_read_cycle_mac[data_type_t::OUTPUT],
                                         u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
        double last_before_stage = std::max(u_write_cycle_lb[data_type_t::OUTPUT],
                                              u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
        double last_stage = u_write_cycle_lb[data_type_t::OUTPUT];

        double other_stage = std::max(ratio*u_read_cycle_mac[data_type_t::OUTPUT],
                               std::max(u_write_cycle_lb[data_type_t::OUTPUT],
                                        u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth))));

        if(num_access_lb == 0) {
            // No transaction: leave the overlapped stage unchanged.
        } else if(num_access_lb == 1) {
            cycle_mac_lb[data_type_t::OUTPUT] += ratio*u_read_cycle_mac[data_type_t::OUTPUT] +
                                                 u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)) +
                                                 u_write_cycle_lb[data_type_t::OUTPUT];
        } else if(num_access_lb == 2) {
            cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage + last_before_stage + last_stage;
        } else {
            cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage +
                                                 (num_access_lb - 2)*other_stage +
                                                 last_before_stage + last_stage;
        }

        transfer_energy[data_type_t::OUTPUT] += num_access_mac*u_transfer_energy*ceil((double)(line_size_mac[data_type_t::OUTPUT])/(double)(bitwidth));
        transfer_cycle[data_type_t::OUTPUT] += num_access_mac*u_transfer_cycle*ceil((double)(line_size_mac[data_type_t::OUTPUT])/(double)(bitwidth));

        // Input data and output data should be requested from MAC.
        exist_data_mac[data_type_t::INPUT] = false, request_to_lb[data_type_t::INPUT] = true;
        exist_data_mac[data_type_t::OUTPUT] = false, request_to_lb[data_type_t::OUTPUT] =true;

        // Update stats (Memory operation)
        num_request_to_lb[data_type_t::INPUT]++;
        num_request_to_lb[data_type_t::OUTPUT]++;

        // When all the input and output data at local buffer are used.
        // A new weight tile should be request from MAC to local buffer.
        if(input_index == m_scheduler->offset_size_pe[data_type_t::INPUT].front() &&
           output_index == m_scheduler->offset_size_pe[data_type_t::OUTPUT].front()) {

            // Request new weight.
            exist_data_mac[data_type_t::WEIGHT] = false, request_to_lb[data_type_t::WEIGHT] = true;

            // Update stats of memory operation.
            /* Stats */
            num_request_to_lb[data_type_t::WEIGHT]++;
            /* Stats */

            // The case when not all weight in local buffer is transferred to MAC.
            if(weight_index < m_scheduler->offset_size_pe[data_type_t::WEIGHT].front()) {
                input_index = 0, output_index = 0;
            }
        }
    }
}

output_stationary_t::output_stationary_t(section_config_t m_section_config) :
    pe_t(m_section_config) {

}

output_stationary_t::~output_stationary_t() {

}

void output_stationary_t::computation(scheduler_t *m_scheduler) {

    // Execute MAC operation
    // When all data exist.
    if(exist_data_mac[data_type_t::INPUT] && exist_data_mac[data_type_t::WEIGHT] && exist_data_mac[data_type_t::OUTPUT]) {
        if(m_scheduler->layer_name == layer_name_t::CONVOLUTIONAL_LAYER ||
           m_scheduler->layer_name == layer_name_t::CONNECTED_LAYER) {
#ifdef FUNCTIONAL
            mac_operation(m_scheduler);
            // Activation is applied only after a completed output reduction.
            if(m_scheduler->compression_type == compression_type_t::DENSE) {
                for(unsigned i = 0; i < num_active_macs; i++) {
                    num_computation++;
                    computation_energy += u_computation_energy;
                }
                computation_cycle += accumulate_issue_cycles(1, u_computation_cycle) +
                     lane_reduction_fill_cycles(lane_state, u_computation_cycle);
                     reduction_energy += lane_reduction_energy(lane_state, u_mac_reduction_energy);
            }
            else {
                const unsigned nonzero_operations = count_nonzero_mac_operations(m_scheduler);
                num_computation += nonzero_operations;
                computation_energy += nonzero_operations*u_computation_energy;
                if(nonzero_operations > 0) {
                    computation_cycle += accumulate_issue_cycles(1, u_computation_cycle) +
                     lane_reduction_fill_cycles(lane_state, u_computation_cycle);
                     reduction_energy += lane_reduction_energy(lane_state, u_mac_reduction_energy);
                }
            }
#else
            for(unsigned i = 0; i < num_active_macs; i++) {
                num_computation++;
                computation_energy += u_computation_energy;
            }
            computation_cycle += accumulate_issue_cycles(1, u_computation_cycle) +
                     lane_reduction_fill_cycles(lane_state, u_computation_cycle);
                     reduction_energy += lane_reduction_energy(lane_state, u_mac_reduction_energy);

#endif
        }

        else {
            std::cerr << "Error: PE computation supports only convolution/connected layers" << std::endl;
            exit(1);
        }
#ifndef FUNCTIONAL
        exist_data_mac[data_type_t::INPUT] = false;
        request_to_lb[data_type_t::INPUT] = true;
        exist_data_mac[data_type_t::WEIGHT] = false;
        request_to_lb[data_type_t::WEIGHT] = true;
        num_request_to_lb[data_type_t::INPUT]++;
        num_request_to_lb[data_type_t::WEIGHT]++;
        if(input_index == m_scheduler->offset_size_pe[data_type_t::INPUT].front() &&
           weight_index == m_scheduler->offset_size_pe[data_type_t::WEIGHT].front()) {
            exist_data_mac[data_type_t::OUTPUT] = false;
            request_to_lb[data_type_t::OUTPUT] = true;
            num_request_to_lb[data_type_t::OUTPUT]++;
            // RE1: a MAC -> LB write-back of a partial sum IS the accumulator spill.
        account_accumulator_spill(tile_size_mac[data_type_t::OUTPUT]);
        account_descriptor_dense_mac_transfer(data_type_t::OUTPUT, tile_size_mac[data_type_t::OUTPUT], false);
            if(output_index < m_scheduler->offset_size_pe[data_type_t::OUTPUT].front()) {
                input_index = 0;
                weight_index = 0;
            }
        }
        return;
#endif


        // After computation.
        // Input data and weight should be request from MAC to local buffer.
        exist_data_mac[data_type_t::INPUT] = false, request_to_lb[data_type_t::INPUT] = true;
        exist_data_mac[data_type_t::WEIGHT] = false, request_to_lb[data_type_t::WEIGHT] = true;

        // Update stats (memory operation).
        /* Stats */
        num_request_to_lb[data_type_t::INPUT]++;
        num_request_to_lb[data_type_t::WEIGHT]++;
        /* Stats */

        // When all the input and weight data at local buffer are used.
        // A new output data tile should be request from MAC to local buffer.
        if(input_index == m_scheduler->offset_size_pe[data_type_t::INPUT].front() &&
           weight_index == m_scheduler->offset_size_pe[data_type_t::WEIGHT].front()) {

            // Request new output data.
            exist_data_mac[data_type_t::OUTPUT] = false, request_to_lb[data_type_t::OUTPUT] = true;

            num_request_to_lb[data_type_t::OUTPUT]++;


#ifdef FUNCTIONAL
            // Write back output data
            m_scheduler->transfer_data(output_data_lb, output_data_mac, m_scheduler->output_offset_pe.front(), 0,
                                       component_type_t::PE, component_type_t::MAC,
                                       data_type_t::OUTPUT, get_mac_stationary_type(), action_type_t::STORE);
        // Update for NPUsim ver2
        //bool last_component = (index == m_scheduler->num_active_pe_x*m_scheduler->num_active_pe_y - 1 &&
        //                       pe_array->index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y);
        //m_scheduler->transfer_data_ver2(output_data_lb, output_data_mac,
        //                                component_type_t::PE, component_type_t::MAC,
        //                                data_type_t::OUTPUT, get_mac_stationary_type()), action_type_t::STORE);
#endif
            std::vector<unsigned> parameters_mac(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_lb(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_mac = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::MAC);
            parameters_lb = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE);

            uint64_t address_mac = 0, address_lb = 0;
            unsigned num_access_mac = 0, num_access_lb = 0;
            for(unsigned b = 0; b < parameters_mac[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned k = 0; k < parameters_mac[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned p = 0; p < parameters_mac[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                        for(unsigned q = 0; q < parameters_mac[parameter_type_t::OUTPUT_WIDTH]; q++) {
                            if(address_mac != ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                           *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                           *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                          k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                           *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                          p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                          mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT]) {

                                // Update write cost of MAC unit.
                                access_energy_mac[data_type_t::OUTPUT] += u_read_energy_mac[data_type_t::OUTPUT];
                                access_cycle_mac[data_type_t::OUTPUT] += u_read_cycle_mac[data_type_t::OUTPUT];
                                num_access_mac++;

                                // Update MAC address
                                address_mac = ((uint64_t)&output_data_mac[b*parameters_mac[parameter_type_t::OUTPUT_CHANNEL]
                                                                           *parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                           *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                          k*parameters_mac[parameter_type_t::OUTPUT_HEIGHT]
                                                                           *parameters_mac[parameter_type_t::OUTPUT_WIDTH] +
                                                                          p*parameters_mac[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                          mask_bits_mac[data_type_t::OUTPUT]) << mask_bits_mac[data_type_t::OUTPUT];
                            }

                            if(address_lb != ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                        b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT]) {

                                // Update read cost of Local buffer
                                access_energy_lb[data_type_t::OUTPUT] += u_write_energy_lb[data_type_t::OUTPUT];
                                access_cycle_lb[data_type_t::OUTPUT] += u_write_cycle_lb[data_type_t::OUTPUT];
                                num_access_lb++;

                                // Update local buffer address
                                address_lb = ((uint64_t)&output_data_lb[m_scheduler->output_offset_pe.front() +
                                                                        b*parameters_lb[parameter_type_t::OUTPUT_CHANNEL]
                                                                         *parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                        k*parameters_lb[parameter_type_t::OUTPUT_HEIGHT]
                                                                         *parameters_lb[parameter_type_t::OUTPUT_WIDTH] +
                                                                        p*parameters_lb[parameter_type_t::OUTPUT_WIDTH] + q] >>
                                                                        mask_bits_lb[data_type_t::OUTPUT]) << mask_bits_lb[data_type_t::OUTPUT];

                            }
                        }
                    }
                }
            }

            // Update overlapped cycle at MAC units and local buffer
            unsigned ratio = ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(line_size_mac[data_type_t::OUTPUT]));

            double first_stage = ratio*u_read_cycle_mac[data_type_t::OUTPUT];
            double second_stage = std::max(ratio*u_read_cycle_mac[data_type_t::OUTPUT],
                                             u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
            double last_before_stage = std::max(u_write_cycle_lb[data_type_t::OUTPUT],
                                                  u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)));
            double last_stage = u_write_cycle_lb[data_type_t::OUTPUT];

            double other_stage = std::max(ratio*u_read_cycle_mac[data_type_t::OUTPUT],
                                   std::max(u_write_cycle_lb[data_type_t::OUTPUT],
                                            u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth))));

            if(num_access_lb == 0) {
                // No transaction: leave the overlapped stage unchanged.
            } else if(num_access_lb == 1) {
                cycle_mac_lb[data_type_t::OUTPUT] += ratio*u_read_cycle_mac[data_type_t::OUTPUT] +
                                                     u_transfer_cycle*ceil((double)(line_size_lb[data_type_t::OUTPUT])/(double)(bitwidth)) +
                                                     u_write_cycle_lb[data_type_t::OUTPUT];
            } else if(num_access_lb == 2) {
                cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage + last_before_stage + last_stage;
            } else {
                cycle_mac_lb[data_type_t::OUTPUT] += first_stage + second_stage +
                                                     (num_access_lb - 2)*other_stage +
                                                     last_before_stage + last_stage;
            }

            transfer_energy[data_type_t::OUTPUT] += num_access_mac*u_transfer_energy*ceil((double)(line_size_mac[data_type_t::OUTPUT])/(double)(bitwidth));
            transfer_cycle[data_type_t::OUTPUT] += num_access_mac*u_transfer_cycle*ceil((double)(line_size_mac[data_type_t::OUTPUT])/(double)(bitwidth));

            // When not all output in local buffer is transferred to MAC.
            if(output_index < m_scheduler->offset_size_pe[data_type_t::OUTPUT].front()) {
                input_index = 0, weight_index = 0;
            }
        }
    }
}
