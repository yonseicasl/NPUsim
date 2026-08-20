#include <iostream>
#include <limits>
#include <cmath>
#include <cstring>
#include "pe_array.h"
#include "datatype.h"
#include "interconnect_timing.h"

pe_array_t::pe_array_t(section_config_t /*m_section_config*/) :
    input_data(NULL),
    weight(NULL),
    output_data(NULL),
    equal_output_tile(false),
    duplicated_input(0),
    index(0),
    reduction_additions(0.0),
    psum_retention_valid(true),
    exist_temporal_buffer(false),
    utilization(0.0),
    write_back_cycle(0.0),
    overlapped_transfer_cycle(0.0),
    u_weight_fold_fill_cycle(0.0),
    u_layer_setup_cycle(0.0),
    u_weight_fold_fill_energy(0.0),
    u_layer_setup_energy(0.0),
    weight_fold_events(0.0),
    weight_fold_energy(0.0),
    accumulator_energy(0.0),
    fold_fill_cycle(0.0),
    edge_accumulation(false),
    global_buffer(NULL),
    stationary_type(stationary_type_t::UNDEFINED_STATIONARY),
    array_parameter_order("kbpqcrs"),
    noc_type(noc_type_t::UNDEFINED_NOC), 
    memory_type(memory_type_t::SEPARATE),
    height(1),
    width(1),
    num_pes(1),
    num_active_pe_x(1),
    num_active_pe_y(1),
    input_size(0),
    weight_size(0),
    output_size(0),
    frequency(0.0),
    bandwidth(0.0),
    bitwidth(0),
    initial(true),
    noc_cycle(0),
    noc_energy(0.0) {
}

pe_array_t::~pe_array_t() {

}
void pe_array_t::initialize_temporal_buffer(section_config_t m_section_config) {
    m_section_config.get_setting("exist_temporal_buffer", &exist_temporal_buffer);

    std::string memory_type_str;
    if(!m_section_config.get_setting("memory_type", &memory_type_str)) {
        std::cerr << "Error: PE array requires memory_type" << std::endl;
        exit(1);
    }
    lowercase(memory_type_str);
    if(memory_type_str == "separated") memory_type_str = "separate";
    if(memory_type_str == "separate") memory_type = memory_type_t::SEPARATE;
    else if(memory_type_str == "shared") memory_type = memory_type_t::SHARED;
    else {
        std::cerr << "Error: invalid PE-array memory_type: " << memory_type_str << std::endl;
        exit(1);
    }

    unsigned per_pe_input = 0;
    unsigned per_pe_weight = 0;
    unsigned per_pe_output = 0;
    const bool have_explicit_buffers =
        m_section_config.get_setting("input_buffer", &per_pe_input) &&
        m_section_config.get_setting("weight_buffer", &per_pe_weight) &&
        m_section_config.get_setting("output_buffer", &per_pe_output);
    if(!have_explicit_buffers) {
        if(!m_section_config.get_setting("input_size", &per_pe_input) ||
           !m_section_config.get_setting("weight_size", &per_pe_weight) ||
           !m_section_config.get_setting("output_size", &per_pe_output)) {
            std::cerr << "Error: PE array requires input/weight/output buffer sizes" << std::endl;
            exit(1);
        }
    }
    if(per_pe_input == 0 || per_pe_weight == 0 || per_pe_output == 0) {
        std::cerr << "Error: PE-array buffer sizes must be non-zero" << std::endl;
        exit(1);
    }
    if(num_pes == 0 || num_pes > std::numeric_limits<unsigned>::max()/per_pe_input ||
       num_pes > std::numeric_limits<unsigned>::max()/per_pe_weight ||
       num_pes > std::numeric_limits<unsigned>::max()/per_pe_output) {
        std::cerr << "Error: PE-array buffer size overflow" << std::endl;
        exit(1);
    }

    input_size = per_pe_input * num_pes;
    weight_size = per_pe_weight * num_pes;
    output_size = per_pe_output * num_pes;

    const size_t num_input = (static_cast<size_t>(input_size) + sizeof(data_t) - 1) / sizeof(data_t);
    const size_t num_weight = (static_cast<size_t>(weight_size) + sizeof(data_t) - 1) / sizeof(data_t);
    const size_t num_output = (static_cast<size_t>(output_size) + sizeof(data_t) - 1) / sizeof(data_t);

    input_data = new data_t[num_input]();
    weight = new data_t[num_weight]();
    output_data = new data_t[num_output]();

    // Unit static (leakage) energy of the PE-array temporal buffer, pJ/cycle. Defaults to 0.
    u_static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("pe_array_static_energy", &u_static_energy);
    static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    buffer_utilization.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // V2: RTL-calibrated weight-residency fold fill and per-layer setup (default off).
    m_section_config.get_setting("weight_fold_fill_cycle", &u_weight_fold_fill_cycle);
    // E5: the energy counterparts of the fold/setup latency bubbles.
    m_section_config.get_setting("weight_fold_fill_energy", &u_weight_fold_fill_energy);
    m_section_config.get_setting("layer_setup_energy", &u_layer_setup_energy);
    m_section_config.get_setting("layer_setup_cycle", &u_layer_setup_cycle);
    // V3: array-edge output accumulation (Gemmini-style accumulator).
    m_section_config.get_setting("edge_accumulation", &edge_accumulation);
}

// Accumulate leakage energy of the PE-array temporal buffer once for the layer duration.
// The temporal buffer is always-on, so it leaks over the full layer elapsed window.
void pe_array_t::update_static_energy(double elapsed_cycles) {
    if(elapsed_cycles < 0.0) {
        std::cerr << "Error: PE-array elapsed cycles must be non-negative" << std::endl;
        exit(1);
    }
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        if(u_static_energy[i] < 0.0) {
            std::cerr << "Error: PE-array static_energy must be a non-negative pJ/cycle value" << std::endl;
            exit(1);
        }
        static_energy[i] = static_energy_for_cycles(u_static_energy[i], elapsed_cycles);
    }
}

// Modeled busy duration of the PE-array temporal buffer: the max of its access, transfer,
// and overlapped cost axes. Feeds the layer leakage window and consumes the otherwise
// write-only overlap/write-back counters.
double pe_array_t::modeled_elapsed_cycles() const {
    double elapsed = std::max(write_back_cycle, overlapped_transfer_cycle);
    for(unsigned type = 0; type < data_type_t::NUM_DATA_TYPES; ++type) {
        elapsed = std::max(elapsed, access_cycle[type]);
        elapsed = std::max(elapsed, transfer_cycle[type]);
        elapsed = std::max(elapsed, cycle_temporal_pe[type]);
    }
    return elapsed;
}

void pe_array_t::account_descriptor_dense_distribution(scheduler_t *m_scheduler,
                                                        double link_cycle, double link_energy,
                                                        double link_fill_cycles) {
    if(m_scheduler->compression_type != compression_type_t::DENSE) {
        std::cerr << "Error: timing PE array supports dense descriptor traffic only" << std::endl;
        exit(1);
    }
    const size_t required_input = runtime_datatypes().storage_bytes(data_type_t::INPUT, tile_size[data_type_t::INPUT]);
    const size_t required_weight = runtime_datatypes().storage_bytes(data_type_t::WEIGHT, tile_size[data_type_t::WEIGHT]);
    const size_t required_output = runtime_datatypes().storage_bytes(data_type_t::OUTPUT, tile_size[data_type_t::OUTPUT]);
    // V3: with edge accumulation the output tile lives in the array-edge accumulator
    // (modeled by the global-buffer output partition), not in the array/PE buffers.
    if(required_input > input_size || required_weight > weight_size ||
       (!edge_accumulation && required_output > output_size)) {
        std::cerr << "Error: runtime datatype PE-array tile exceeds temporal-buffer capacity" << std::endl;
        exit(1);
    }
    // PA6: record the temporal-buffer occupancy per data type (peak across the layer),
    // complementing the PE-count `utilization` which only reflects spatial mapping.
    buffer_utilization[data_type_t::INPUT] = std::max(buffer_utilization[data_type_t::INPUT],
        static_cast<double>(required_input)/static_cast<double>(input_size));
    buffer_utilization[data_type_t::WEIGHT] = std::max(buffer_utilization[data_type_t::WEIGHT],
        static_cast<double>(required_weight)/static_cast<double>(weight_size));
    if(!edge_accumulation) {
        buffer_utilization[data_type_t::OUTPUT] = std::max(buffer_utilization[data_type_t::OUTPUT],
            static_cast<double>(required_output)/static_cast<double>(output_size));
    }
    for(unsigned type_index = 0; type_index < data_type_t::NUM_DATA_TYPES; ++type_index) {
        const data_type_t type = static_cast<data_type_t>(type_index);
        bool requested = false;
        for(unsigned i = 0; i < get_number_of_active_pes(); ++i) {
            requested = requested || pes[i]->request_to_pe_array[type];
        }
        if(!requested) continue;

        if(!skip_transfer[type] && (type != data_type_t::OUTPUT || global_buffer->transfer_output)) {
            num_data_transfer[type]++;
            // V2: every array-wide weight residency stalls the compute schedule while the
            // new weights load down the columns and the previous accumulation drains
            // (RTL-measured ~DIM+1 cycles/fold on Gemmini). Each per-PE temporal weight
            // element is one sequential residency, so a distribution event carries
            // per-PE-weight-tile-size folds.
            //
            // SY2/L9: a config that has NOT calibrated that bubble from RTL now gets an
            // analytical one instead of nothing -- a systolic array's accumulation pipeline
            // must drain down the active columns before different weights take effect, which
            // an injection-wavefront-only model misses entirely.
            //
            // The two paths deliberately use DIFFERENT multipliers, because they measure
            // different things:
            //   * a calibrated weight_fold_fill_cycle is an RTL measurement of the bubble per
            //     per-PE weight ELEMENT (Gemmini reloads the array once per element, so each
            //     element is its own sequential residency), and it wins outright -- adding the
            //     analytical drain on top would double count the same stall;
            //   * the analytical drain is charged ONCE PER ARRAY-WIDE WEIGHT TILE, i.e. at the
            //     tile boundary. Charging it per element instead would assume the pipeline
            //     fully drains between consecutive weights of the SAME tile, contradicting the
            //     pipelining this model asserts everywhere else (on a 12x14 array with 352
            //     weights per PE it inflated the compute schedule 10x over the computation
            //     itself). It is therefore a LOWER bound: a design that really does drain per
            //     element must calibrate weight_fold_fill_cycle.
            if(type == data_type_t::WEIGHT) {
                // E5: the ENERGY of a residency follows the same event as its latency, and the
                // same selection rule -- a calibrated per-element fold cost is charged per
                // per-PE weight element, while the analytical drain is a per-tile boundary
                // bubble. Whichever supplies the latency supplies the event count, so the two
                // can never be double charged.
                if(u_weight_fold_fill_cycle > 0.0) {
                    const double folds =
                        static_cast<double>(pes[0]->tile_size_lb[data_type_t::WEIGHT]);
                    fold_fill_cycle += u_weight_fold_fill_cycle*folds;
                    weight_fold_events += folds;
                } else {
                    fold_fill_cycle += weight_fold_bubble_cycles();
                    if(weight_fold_bubble_cycles() > 0.0) weight_fold_events += 1.0;
                }
                weight_fold_energy += u_weight_fold_fill_energy*
                    (u_weight_fold_fill_cycle > 0.0
                        ? static_cast<double>(pes[0]->tile_size_lb[data_type_t::WEIGHT])
                        : (weight_fold_bubble_cycles() > 0.0 ? 1.0 : 0.0));
            }
            // PA1: operands shared by several PEs are read from the temporal buffer and
            // streamed over the fabric ONCE -- the array-level tile is the distinct data
            // set, and the fabric multicasts it. Per-destination replication is only the
            // local-buffer write below. (The old model charged the source read and the
            // link stream once per destination PE, over-counting reuse-heavy dataflows
            // by the sharing factor.)
            const datatype_transfer_timing_t shared = datatype_transfer_timing(
                type, tile_size[type], line_size[type], line_size[type], bitwidth);
            payload_link_transactions[type] += shared.payload_link_transactions;
            metadata_link_transactions[type] += shared.metadata_link_transactions;
            storage_link_transactions[type] += shared.link_transactions;
            // Reading the distributed tile out of the array temporal buffer only
            // costs cycles/energy when that buffer exists (pass-through otherwise).
            if(exist_temporal_buffer) {
                access_cycle[type] += shared.source_accesses*u_read_cycle[type];
                access_energy[type] += shared.source_accesses*u_read_energy[type];
            }
            // SP1/SY2: hops pipeline, so the route depth is a one-time fill per stream.
            transfer_cycle[type] += shared.link_transactions*link_cycle + link_fill_cycles;
            transfer_energy[type] += shared.link_transactions*link_energy;

            size_t max_destination_accesses = 0;
            double destination_write_cycle = 0.0;
            for(unsigned i = 0; i < get_number_of_active_pes(); ++i) {
                const size_t destination_accesses = runtime_datatypes().storage_transactions(
                    type, pes[i]->tile_size_lb[type], pes[i]->line_size_lb[type]);
                pes[i]->access_cycle_lb[type] += destination_accesses*pes[i]->u_write_cycle_lb[type];
                pes[i]->access_energy_lb[type] += destination_accesses*pes[i]->u_write_energy_lb[type];
                if(destination_accesses > max_destination_accesses) {
                    max_destination_accesses = destination_accesses;
                    destination_write_cycle = pes[i]->u_write_cycle_lb[type];
                }
                size_t capacity = pes[i]->output_size;
                if(type == data_type_t::INPUT) capacity = pes[i]->input_size;
                else if(type == data_type_t::WEIGHT) capacity = pes[i]->weight_size;
                // V3: edge-accumulated outputs stream through, so LB occupancy is undefined.
                if(!(edge_accumulation && type == data_type_t::OUTPUT)) {
                    pes[i]->utilization_local_buffer[type] = std::max(
                        pes[i]->utilization_local_buffer[type],
                        static_cast<double>(runtime_datatypes().storage_bytes(type, pes[i]->tile_size_lb[type])) /
                        static_cast<double>(capacity));
                }
                pes[i]->skip_transfer[type] = false;
            }

            const size_t pipeline_transactions = std::max(
                std::max(shared.source_accesses, shared.link_transactions), max_destination_accesses);
            if(pipeline_transactions > std::numeric_limits<unsigned>::max()) {
                std::cerr << "Error: PE-array pipeline transaction count overflow" << std::endl;
                exit(1);
            }
            if(exist_temporal_buffer) {
                // L2: the source/link packet boundaries come from the shared array tile;
                // the destination is one PE's SLICE of that tile written into its own local
                // buffer, so its total is spread over the same packets (see
                // with_destination_total()).
                cycle_temporal_pe[type] += pipelined_transfer_cycles(
                    shared.groups.with_destination_total(max_destination_accesses),
                    u_read_cycle[type], link_cycle, destination_write_cycle) + link_fill_cycles;
            }
        }
        for(unsigned i = 0; i < get_number_of_active_pes(); ++i) {
            pes[i]->exist_data_lb[type] = true;
            pes[i]->request_to_pe_array[type] = false;
        }
        is_waiting_data[type] = false;
    }
    for(unsigned i = 0; i < get_number_of_active_pes(); ++i) pes[i]->fill_data();
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        if(tile_size[i] == global_buffer->tile_size[i]) skip_transfer[i] = true;
    }
}

void pe_array_t::account_descriptor_dense_writeback(pe_t *source_pe, size_t elements) {
    const datatype_transfer_timing_t timing = datatype_transfer_timing(
        data_type_t::OUTPUT, elements, source_pe->line_size_lb[data_type_t::OUTPUT],
        line_size[data_type_t::OUTPUT], bitwidth);
    payload_link_transactions[data_type_t::OUTPUT] += timing.payload_link_transactions;
    metadata_link_transactions[data_type_t::OUTPUT] += timing.metadata_link_transactions;
    storage_link_transactions[data_type_t::OUTPUT] += timing.link_transactions;
    source_pe->access_cycle_lb[data_type_t::OUTPUT] += timing.source_accesses*source_pe->u_read_cycle_lb[data_type_t::OUTPUT];
    source_pe->access_energy_lb[data_type_t::OUTPUT] += timing.source_accesses*source_pe->u_read_energy_lb[data_type_t::OUTPUT];
    source_pe->write_back_cycle_lb += timing.source_accesses*source_pe->u_read_cycle_lb[data_type_t::OUTPUT];
    // A pass-through array (no temporal buffer) adds no buffer write on the
    // write-back path; outputs stream from the PE LBs straight toward the GLB.
    if(exist_temporal_buffer) {
        access_cycle[data_type_t::OUTPUT] += timing.destination_accesses*u_write_cycle[data_type_t::OUTPUT];
        access_energy[data_type_t::OUTPUT] += timing.destination_accesses*u_write_energy[data_type_t::OUTPUT];
        write_back_cycle += timing.destination_accesses*u_write_cycle[data_type_t::OUTPUT];
    }
    // PA5/SP1/SY2: apply the NoC topology cost to the output write-back too, mirroring
    // the distribution path -- per-hop energy per transaction, and the route depth as
    // a one-time pipeline fill (0 for single-hop fabrics, so BUS/S&F/crossbar
    // unchanged). writeback_noc_type() lets a systolic array force the same MESH
    // diameter treatment its load path already uses (see pe_array.h).
    const spatial_noc_cost_t topology_cost = spatial_noc_cost(writeback_noc_type(), num_active_pe_y, num_active_pe_x,
                                                       /* multicast = */ false);
    const double topology_energy = noc_energy*topology_cost.link_traversals;
    const double link_fill_cycles = topology_cost.latency_fill_hops*noc_cycle;
    transfer_cycle[data_type_t::OUTPUT] += timing.link_transactions*noc_cycle + link_fill_cycles;
    transfer_energy[data_type_t::OUTPUT] += timing.link_transactions*topology_energy;
    if(timing.pipeline_transactions > std::numeric_limits<unsigned>::max()) {
        std::cerr << "Error: PE-array write-back transaction count overflow" << std::endl;
        exit(1);
    }
    if(exist_temporal_buffer) {
        // L2: pipeline the physical packet decomposition of the write-back stream.
        cycle_temporal_pe[data_type_t::OUTPUT] += pipelined_transfer_cycles(
            timing.groups, source_pe->u_read_cycle_lb[data_type_t::OUTPUT], noc_cycle,
            u_write_cycle[data_type_t::OUTPUT]) + link_fill_cycles;
    }
}

void pe_array_t::connect(global_buffer_t *m_global_buffer) {
    global_buffer = m_global_buffer;
    index = global_buffer->index;

    for(unsigned i = 0; i < num_pes; i++) {
        pes[i]->connect(this);
    }
}
// Get the stationary type between the PE array and Global buffer.
stationary_type_t pe_array_t::get_stationary_type() { return stationary_type; }

// Return parameter order.
std::string pe_array_t::get_parameter_order() { return array_parameter_order; }

// Get the number of PEs.
unsigned pe_array_t::get_number_of_pes() { return num_pes; }

unsigned pe_array_t::get_number_of_active_pes() { return num_active_pe_x*num_active_pe_y; }

// A signal that checks whether the PE array is idle state or not.
bool pe_array_t::is_idle() {
    bool idle = true;
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        idle = pes[i]->is_idle();
        if(idle == false) {break;}
    }
    return idle;
}

// Check whether the request at PE exists.
bool pe_array_t::is_exist_request_at_pe() {
    bool request_input_ = false, request_weight_ = false, request_output_ = false;

    // At least one PE sends a request to PE array
    // Input request becomes true.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::INPUT]) {
            request_input_ = true;
            break;
        }
    }
    // At least one PE sends a request to PE array. 
    // Weight request becomes true.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::WEIGHT]) {
            request_weight_ = true;
            break;
        }
    }
    // At least one PE sends a request to PE array.
    // Output request becomes true.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::OUTPUT]) {
            request_output_ = true;
            break;
        }
    }

    if(request_input_ || request_weight_ || request_output_) { return true; }
    else { return false; }
}


// Check whether the request at temporal buffer exists.
bool pe_array_t::is_exist_request_at_buffer() {
    if(request_to_global_buffer[data_type_t::INPUT] || 
       request_to_global_buffer[data_type_t::WEIGHT] || 
       request_to_global_buffer[data_type_t::OUTPUT]) { 
        return true; 
    }
    else { return false; }
}

bool pe_array_t::is_exist_data() {
    if(exist_data[data_type_t::INPUT] && 
       exist_data[data_type_t::WEIGHT] && 
       exist_data[data_type_t::OUTPUT]) { 
        return true; 
    }
    else { return false; }
}


// Wait for the data comes from Global buffer.
bool pe_array_t::wait_data() {
    if(is_waiting_data[data_type_t::INPUT] || is_waiting_data[data_type_t::WEIGHT] || is_waiting_data[data_type_t::OUTPUT]) {
        return true;
    }
    else {
        return false;
    }
}

// A signal that the data exist in the PE array.
void pe_array_t::fill_data() {
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        pes[i]->fill_data();
    }
}

void pe_array_t::request_data() {
    // The case when at least on PE in PE array request data.
    // Request data from temporal buffer in PE array to Global buffer.
    
    // Request input data to Global buffer.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::INPUT]) {

            request_to_global_buffer[data_type_t::INPUT] = true;
            is_waiting_data[data_type_t::INPUT] = true;
            global_buffer->num_request[data_type_t::INPUT]++;
            break;
        }
    }
    // Request weight to Global buffer.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::WEIGHT]) {

            request_to_global_buffer[data_type_t::WEIGHT] = true;
            is_waiting_data[data_type_t::WEIGHT] = true;
            global_buffer->num_request[data_type_t::WEIGHT]++;
            break;
        }
    }
    // Request output data to Global buffer.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::OUTPUT]) {
            if(!initial && !equal_output_tile) {
#ifdef PRINT
                std::cout << "Write back output data from temporal buffer in PE array to Global buffer" << std::endl;
#endif
#ifdef FUNCTIONAL
#endif
                is_waiting_data[data_type_t::OUTPUT] = true;
                global_buffer->num_request[data_type_t::OUTPUT]++;
                const datatype_transfer_timing_t timing = datatype_transfer_timing(
                    data_type_t::OUTPUT, tile_size[data_type_t::OUTPUT], line_size[data_type_t::OUTPUT],
                    global_buffer->line_size[data_type_t::OUTPUT], global_buffer->get_bitwidth());
                const size_t pe_array_lines = timing.source_accesses;
                const size_t global_buffer_lines = timing.destination_accesses;
                const size_t link_transactions = timing.link_transactions;
                global_buffer->payload_link_transactions[data_type_t::OUTPUT] += timing.payload_link_transactions;
                global_buffer->metadata_link_transactions[data_type_t::OUTPUT] += timing.metadata_link_transactions;
                global_buffer->storage_link_transactions[data_type_t::OUTPUT] += timing.link_transactions;

                // GB3: hide the buffer access cycles when the destination global buffer is
                // double-buffered, matching the load-path convention (multi_chip.cc
                // distribution gates on the destination's flag). Energy is always charged.
                if(!global_buffer->double_buffer) {
                    if(exist_temporal_buffer) {
                        access_cycle[data_type_t::OUTPUT] += pe_array_lines*u_read_cycle[data_type_t::OUTPUT];
                        write_back_cycle += pe_array_lines*u_read_cycle[data_type_t::OUTPUT];
                    }
                    global_buffer->access_cycle[data_type_t::OUTPUT] += global_buffer_lines*global_buffer->u_write_cycle[data_type_t::OUTPUT];
                    global_buffer->write_back_cycle += global_buffer_lines*global_buffer->u_write_cycle[data_type_t::OUTPUT];
                }
                if(exist_temporal_buffer)
                    access_energy[data_type_t::OUTPUT] += pe_array_lines*u_read_energy[data_type_t::OUTPUT];
                global_buffer->access_energy[data_type_t::OUTPUT] += global_buffer_lines*global_buffer->u_write_energy[data_type_t::OUTPUT];

                global_buffer->transfer_cycle[data_type_t::OUTPUT] += global_buffer->u_transfer_cycle*link_transactions;
                global_buffer->overlapped_transfer_cycle += global_buffer->u_transfer_cycle*link_transactions;
                global_buffer->transfer_energy[data_type_t::OUTPUT] += global_buffer->u_transfer_energy*link_transactions;
                // Does the TILING alone suggest the array's output tile need not move?
                bool tiles_suggest_retention = false;
                if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
                    tiles_suggest_retention =
                        tile_size[data_type_t::OUTPUT] == global_buffer->tile_size[data_type_t::OUTPUT];
                }
                else {
                    if(global_buffer->multi_chip->tile_size[data_type_t::OUTPUT] == global_buffer->multi_chip->dram->tile_size[data_type_t::OUTPUT]) {
                        tiles_suggest_retention =
                            tile_size[data_type_t::OUTPUT] == global_buffer->tile_size[data_type_t::OUTPUT];
                    }
                    else {
                        tiles_suggest_retention =
                            tile_size[data_type_t::OUTPUT] == global_buffer->tile_size[data_type_t::OUTPUT] &&
                            tile_size[data_type_t::WEIGHT] != global_buffer->tile_size[data_type_t::WEIGHT];
                    }
                }
                // E20-3: equal tiles alone do not justify retention, whatever the stationary type.
                // The array holds ONE output tile; if a reduction dimension is split ABOVE it, the
                // loop walks every other tile of that level and comes back, so the psum is
                // physically read back and written out each time. Retention is real only when the
                // reduction for a tile finishes before the array moves on. The three branches above
                // compare one tile while the loop above them iterates many, which is what suppressed
                // the traffic the measured Eyeriss chip spends most of its GLB bandwidth on.
                if(tiles_suggest_retention && psum_retention_valid) {equal_output_tile = true;}
            }
            else {
                initial = 0;
            }
            request_to_global_buffer[data_type_t::OUTPUT] = true;
            
            break;
        }
    }
}

// Flush data at temporal buffer in PE array.
void pe_array_t::flush_data() {
    // At least one PE sends an input request to PE array.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::INPUT]) {
            exist_data[data_type_t::INPUT] = false;
            break;
        }
    }
    // At least one PE sends a weight request to PE array.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::WEIGHT]) {
            exist_data[data_type_t::WEIGHT] = false;
            break;
        }
    }
    // At least one PE sends an output request to PE array.
    for(unsigned i = 0; i < get_number_of_active_pes(); i++) {
        if(pes[i]->request_to_pe_array[data_type_t::OUTPUT]) {
            exist_data[data_type_t::OUTPUT] = false;
            break;
        }
    }
}

void pe_array_t::reset() {
    for(unsigned i = 0; i < num_pes; i++) {
        pes[i]->reset();
    }

    std::fill_n(input_data, (static_cast<size_t>(input_size) + sizeof(data_t) - 1)/sizeof(data_t), data_t{});
    std::fill_n(weight, (static_cast<size_t>(weight_size) + sizeof(data_t) - 1)/sizeof(data_t), data_t{});
    std::fill_n(output_data, (static_cast<size_t>(output_size) + sizeof(data_t) - 1)/sizeof(data_t), data_t{});

    initial = true;
    equal_output_tile = false;

    utilization = 0.0;
    write_back_cycle = 0.0, overlapped_transfer_cycle = 0.0;
    // V2: fold fills accrue per weight residency; the one-time per-layer setup
    // (u_layer_setup_cycle) is added at the stats level AFTER repetition scaling.
    fold_fill_cycle = 0.0;
    // E5: the fold energy and its event count reset with the layer.
    weight_fold_energy = 0.0;
    accumulator_energy = 0.0;
    weight_fold_events = 0.0;

    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);

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
    buffer_utilization.assign(data_type_t::NUM_DATA_TYPES, 0.0);
}

std::string pe_array_t::describe_noc_link_contract() const {
    return std::string("distribution -- ") + spatial_noc_link_contract(noc_type, true) +
           "; write-back -- " + spatial_noc_link_contract(writeback_noc_type(), false);
}

// E20-3: record whether a psum may stay in the array for this layer. Retention is legitimate only
// when no reduction dimension is split above the array; otherwise the array walks every other
// output tile of the outer level and comes back, which is a physical read-back and write-out.
void pe_array_t::set_psum_retention_scope(scheduler_t *m_scheduler) {
    psum_retention_valid = !m_scheduler->mapping_table->reduction_tiled_above_array();
}
