#include <iomanip>
#include <limits>
#include <cmath>
#include <cstring>
#include "global_buffer.h"
#include "datatype.h"
#include "interconnect_timing.h"
#include "energy_units.h"


namespace {

void validate_global_buffer_widths(unsigned link_bits, const std::vector<unsigned> &line_bits) {
    if(link_bits == 0) {
        std::cerr << "Error: global buffer bitwidth must be non-zero" << std::endl;
        exit(1);
    }
    for(unsigned i = 0; i < line_bits.size(); ++i) {
        const unsigned width = line_bits[i];
        if(width < 8 || (width & (width - 1)) != 0) {
            std::cerr << "Error: global buffer line_size must be a power-of-two bit width of at least 8" << std::endl;
            exit(1);
        }
    }
}

} // namespace
global_buffer_t::global_buffer_t(section_config_t /*m_section_config*/) :
    multi_chip(NULL),
    data(NULL),
    double_buffer(false),
    output_cast_bytes(0),
    output_cast_energy(0.0),
    index(0),
    duplicated_input(0),
    u_transfer_cycle(0.0),
    u_transfer_energy(0.0),
    write_back_cycle(0.0),
    overlapped_transfer_cycle(0.0),
    transfer_output(false),
    pe_array(NULL),
    memory_type(memory_type_t::UNDEFINED_MEMORY),
    stationary_type(stationary_type_t::UNDEFINED_STATIONARY),
    parameter_order("KBPQCRS"),
    size(0),
    frequency(0.0),
    bandwidth(0.0),
    bitwidth(0.0),
    input_index(0),
    weight_index(0),
    output_index(0),
    input_flush_counter(0),
    weight_flush_counter(0),
    output_flush_counter(0),
    idle(false),
    initial(true) {
}

global_buffer_t::~global_buffer_t() {

}

void global_buffer_t::connect(pe_array_t *m_pe_array) {
    pe_array = m_pe_array;
}

void global_buffer_t::connect(multi_chip_t *m_multi_chip) {
    multi_chip = m_multi_chip;
}

void global_buffer_t::update_tile_size(scheduler_t *m_scheduler) {
    tile_size = m_scheduler->tile_size[component_type_t::GLOBAL_BUFFER];
    
    update_offset();
    check_tile_size();
}

// Get global buffer's stationary type
stationary_type_t global_buffer_t::get_stationary_type() { return stationary_type; }

// Get global buffer's memory type
memory_type_t global_buffer_t::get_memory_type() { return memory_type;}

// Get global buffer's size
double global_buffer_t::get_buffer_size() { return size; }
bool global_buffer_t::can_retain_input_halo(size_t input_working_set_elements) const {
    if(bypass[data_type_t::INPUT] || input_working_set_elements == 0) return false;
    const size_t copies = double_buffer ? 2 : 1;
    const size_t input_bytes = runtime_datatypes().storage_bytes(data_type_t::INPUT,
                                                                 input_working_set_elements);
    if(input_bytes > std::numeric_limits<size_t>::max()/copies) return false;
    if(memory_type == memory_type_t::SEPARATE) {
        return static_cast<double>(copies*input_bytes) <= capacity_per_type[data_type_t::INPUT];
    }
    size_t working_bytes = input_bytes;
    for(unsigned i = data_type_t::WEIGHT; i < data_type_t::NUM_DATA_TYPES; ++i) {
        const data_type_t type = static_cast<data_type_t>(i);
        if(!bypass[type]) {
            const size_t resident_bytes = runtime_datatypes().storage_bytes(type, tile_size[type]);
            if(resident_bytes > std::numeric_limits<size_t>::max() - working_bytes) return false;
            working_bytes += resident_bytes;
        }
    }
    if(working_bytes > std::numeric_limits<size_t>::max()/copies) return false;
    return static_cast<double>(copies*working_bytes) <= size;
}


// Get global buffer's bitwidth
unsigned global_buffer_t::get_bitwidth() { return bitwidth; }

// A signal to check whether the Global buffer is idle state or not.
bool global_buffer_t::is_idle() { return idle; }

// A signal to check whether all data exist in the Global buffer or not.
bool global_buffer_t::is_exist_data() {
    if(exist_data[data_type_t::INPUT] && exist_data[data_type_t::WEIGHT] && exist_data[data_type_t::OUTPUT]) {
        return true;
    }
    else { return false; }
}

// A signal to check whether at least one request exists from Global buffer to DRAM or not.
bool global_buffer_t::is_exist_request() {
    if(request_to_multi_chip[data_type_t::INPUT] || request_to_multi_chip[data_type_t::WEIGHT] || request_to_multi_chip[data_type_t::OUTPUT]) {
        return true;
    }
    else { return false; }
}

// Wait for the data comes from DRAM.
void global_buffer_t::wait_data() { idle = true; }

// A signal indicates that the data exists in the Global buffer.
void global_buffer_t::fill_data() { idle = false; }

// Send a request signal to DRAM.
void global_buffer_t::request_data() {
    // When the data does not exist in Global buffer
    // Send a request signal to DRAM.
    // Case 1. Input data
    if(!exist_data[data_type_t::INPUT]) {
        // Send a input data request signal to DRAM.
        request_to_multi_chip[data_type_t::INPUT] = true;
        multi_chip->exist_data[data_type_t::INPUT] = false;

        multi_chip->num_request[data_type_t::INPUT]++;
    }
    // Case 2. Weight data
    if(!exist_data[data_type_t::WEIGHT]) {
        // Send a weight request signal to DRAM 
        request_to_multi_chip[data_type_t::WEIGHT] = true;
        multi_chip->exist_data[data_type_t::WEIGHT] = false;

        multi_chip->num_request[data_type_t::WEIGHT]++;
    }
    // Case 3. Output data
    if(!exist_data[data_type_t::OUTPUT]) {
        if(initial) {
            initial = false;
        }
        request_to_multi_chip[data_type_t::OUTPUT] = true;
        multi_chip->exist_data[data_type_t::OUTPUT] = false;

        multi_chip->num_request[data_type_t::OUTPUT]++;
    }
}

// P1-B GLB bypass contract (end to end): a bypassed datatype bypasses the GLB's
// STORAGE, not the chip's physical interconnect. The tile still arrives over the NoP
// (charged in multi_chip_t::account_descriptor_dense_distribution()) and still has to
// be delivered across the on-chip GLB<->PE-array fabric to reach the array, but it is
// never written into nor read back out of the GLB SRAM. So a bypassed datatype:
//   * pays NO GLB write (fill) access -- dropped in multi_chip_t::account_descriptor_dense_distribution()
//   * pays NO GLB read access        -- dropped below
//   * DOES pay the GLB<->PE-array link transactions and the PE-array destination write
//     (the wires and the temporal buffer are physically traversed either way)
// The previous behavior was the exact inverse -- it kept the mc->GLB fill and dropped
// the GLB->PE-array link/destination -- which charged the storage it claimed to bypass
// while making the delivery that actually happens free, and left the bypassed stream
// with zero on-chip traffic for the T4/T5 checks to validate.
void global_buffer_t::account_descriptor_dense_transfer(data_type_t type) {
    // A bypassed stream is sourced by the chip's NoP ingress, not by the GLB SRAM, so
    // its transaction counts come off the upstream line width and no GLB read is charged.
    const bool bypassed = bypass[type];
    const size_t source_line = bypassed ? multi_chip->line_size[type] : line_size[type];
    const datatype_transfer_timing_t timing = datatype_transfer_timing(
        type, pe_array->tile_size[type], source_line, pe_array->line_size[type], bitwidth);

    num_data_transfer[type]++;
    payload_link_transactions[type] += timing.payload_link_transactions;
    metadata_link_transactions[type] += timing.metadata_link_transactions;
    storage_link_transactions[type] += timing.link_transactions;
    if(!bypassed) {
        access_cycle[type] += timing.source_accesses*u_read_cycle[type];
        access_energy[type] += timing.source_accesses*u_read_energy[type];
    }
    // A pass-through PE array (no temporal buffer) imposes no buffer write cost;
    // the stream lands directly in the PE local buffers.
    if(pe_array->exist_temporal_buffer) {
        pe_array->access_cycle[type] += timing.destination_accesses*pe_array->u_write_cycle[type];
        pe_array->access_energy[type] += timing.destination_accesses*pe_array->u_write_energy[type];
    }
    transfer_cycle[type] += timing.link_transactions*u_transfer_cycle;
    transfer_energy[type] += timing.link_transactions*u_transfer_energy;
    if(timing.pipeline_transactions > std::numeric_limits<unsigned>::max()) {
        std::cerr << "Error: global buffer pipeline transaction count overflow" << std::endl;
        exit(1);
    }
    // A pass-through PE array has no temporal-buffer write stage (see above): drop the
    // destination leg from the overlap pipeline too, so cycle_pe_array_global_buffer
    // doesn't charge for a stage that doesn't exist.
    // L2: pipeline the physical packet decomposition. P1-B: a bypassed stream has no GLB
    // read stage here (its source stage is the NoP delivery already charged upstream), and
    // a pass-through PE array has no temporal-buffer write stage -- drop absent stages
    // rather than charging them at zero cost, which would still make them pipeline stages.
    transfer_packet_groups_t groups = timing.groups;
    if(bypassed) groups = groups.without_source();
    if(!pe_array->exist_temporal_buffer) groups = groups.without_destination();
    if(streams_pipelined) {
        // Streamed boundary: only the marginal beats add wall time; the stage fill is a
        // one-time cost already represented in the layer's pipeline composition.
        cycle_pe_array_global_buffer[type] += static_cast<double>(timing.link_transactions)*u_transfer_cycle;
    } else {
        cycle_pe_array_global_buffer[type] += pipelined_transfer_cycles(
            groups, u_read_cycle[type], u_transfer_cycle, pe_array->u_write_cycle[type]);
    }
    pe_array->skip_transfer[type] = false;
}

// Account a GLB->multi-chip OUTPUT write-back using runtime datatype transactions,
// mirroring the multi-chip->DRAM write-back path. This charges the GLB read access,
// the multi-chip temporal-buffer write access, the serialized NoP link transfer, and
// the write-back overlap counters -- all from packing-aware descriptor line/link
// counts, so the off-chip output store no longer reports zero timing/energy and no
// longer depends on host-pointer address granularity.
void global_buffer_t::account_output_writeback_link() {
    // P1-B: a bypassed OUTPUT bypasses the GLB SRAM, not the write-back path itself --
    // the tile still leaves the chip over the NoP and still lands in the multi-chip
    // temporal buffer. Only this chip's GLB read is dropped. (Returning early, as this
    // used to, silently deleted the whole off-chip output store for a bypassed OUTPUT.)
    const bool bypassed = bypass[data_type_t::OUTPUT];
    const size_t source_line = bypassed ? pe_array->line_size[data_type_t::OUTPUT]
                                        : line_size[data_type_t::OUTPUT];
    const datatype_transfer_timing_t timing = datatype_transfer_timing(
        data_type_t::OUTPUT, tile_size[data_type_t::OUTPUT], source_line,
        multi_chip->line_size[data_type_t::OUTPUT], multi_chip->get_bitwidth());

    // GB3: buffer access cycles are hidden when the destination (multi-chip temporal
    // buffer) is double-buffered, matching the load-path convention. Energy always charged.
    if(!multi_chip->double_buffer) {
        if(!bypassed) {
            access_cycle[data_type_t::OUTPUT] += timing.source_accesses*u_read_cycle[data_type_t::OUTPUT];
            write_back_cycle += timing.source_accesses*u_read_cycle[data_type_t::OUTPUT];
        }
        multi_chip->access_cycle[data_type_t::OUTPUT] += timing.destination_accesses*multi_chip->u_write_cycle[data_type_t::OUTPUT];
        multi_chip->write_back_cycle += timing.destination_accesses*multi_chip->u_write_cycle[data_type_t::OUTPUT];
    }
    if(!bypassed) {
        access_energy[data_type_t::OUTPUT] += timing.source_accesses*u_read_energy[data_type_t::OUTPUT];
    }
    multi_chip->access_energy[data_type_t::OUTPUT] += timing.destination_accesses*multi_chip->u_write_energy[data_type_t::OUTPUT];
    // RE1: the final cast is NOT charged here either. The GLB reads the output out once per
    // reduction pass (4 passes in the GEMM fixture's GLB-level C split), so a cast charged at this
    // boundary counts 4x the final output elements. It is charged at the off-chip output store,
    // the one boundary DR6/T1 establish fires exactly once per output element.

    // Serialized NoP link transfer over the GLB<->multi-chip fabric. On a mesh NoP the
    // stream burns per-hop energy and pays this chip's route depth as a one-time fill.
    const unsigned nop_hops = multi_chip->nop_hops_for_chip(index);
    const double nop_fill = static_cast<double>(nop_hops - 1)*multi_chip->u_transfer_cycle;
    multi_chip->transfer_cycle[data_type_t::OUTPUT] += multi_chip->u_transfer_cycle*timing.link_transactions + nop_fill;
    multi_chip->transfer_energy[data_type_t::OUTPUT] += multi_chip->u_transfer_energy*timing.link_transactions*nop_hops;
    multi_chip->overlapped_transfer_cycle += multi_chip->u_transfer_cycle*timing.link_transactions + nop_fill;
    multi_chip->payload_link_transactions[data_type_t::OUTPUT] += timing.payload_link_transactions;
    multi_chip->metadata_link_transactions[data_type_t::OUTPUT] += timing.metadata_link_transactions;
    multi_chip->storage_link_transactions[data_type_t::OUTPUT] += timing.link_transactions;
}

// Transfer the data to temporal buffer of PE array.
void global_buffer_t::data_transfer(scheduler_t *m_scheduler) {
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        const data_type_t type = static_cast<data_type_t>(i);
        if(!bypass[type]) {
            // GB4: a separate buffer's per-type occupancy is over its own partition,
            // not the summed capacity of all three partitions (capacity_per_type is the
            // partition for separate buffers and the total size for shared buffers).
            const double capacity = capacity_per_type[type];
            utilization[type] = (capacity > 0.0)
                ? static_cast<double>(runtime_datatypes().storage_bytes(type, tile_size[type]))/capacity
                : 0.0;
        }
    }
#ifndef FUNCTIONAL
    if(m_scheduler->compression_type != compression_type_t::DENSE) {
        std::cerr << "Error: timing Global Buffer supports dense descriptor traffic only; sparse metadata timing is not modeled" << std::endl;
        exit(1);
    }
#endif
    //std::cout << tile_size[data_type_t::INPUT] << " " << tile_size[data_type_t::WEIGHT] << " " << tile_size[data_type_t::OUTPUT] << " " << size << std::endl;
    //std::cout << utilization[data_type_t::INPUT] << utilization[data_type_t::WEIGHT] << utilization[data_type_t::OUTPUT] << std::endl;


    // Transfer input data from Global buffer to temporal buffer of PE array.
    if(pe_array->request_to_global_buffer[data_type_t::INPUT]) {
#ifndef FUNCTIONAL
        if(m_scheduler->compression_type == compression_type_t::DENSE && !skip_transfer[data_type_t::INPUT]) {
            account_descriptor_dense_transfer(data_type_t::INPUT);
        }
#endif
#ifdef FUNCTIONAL
        // Input data transfer
        m_scheduler->transfer_data(pe_array->input_data, data, 0, offsets[data_type_t::INPUT] + m_scheduler->input_offset_global_buffer.front(), 
                                   component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, 
                                   data_type_t::INPUT, pe_array->get_stationary_type(), action_type_t::LOAD);

        // Update for NPUsim ver2
        //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
        //m_scheduler->transfer_data_ver2(pe_array->input_data, data, 
        //                                component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, 
        //                                data_type_t::INPUT, pe_array->get_stationary_type(), action_type_t::LOAD, true);
                                   
        // Case 1. Dense data format
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::INPUT]) {
                // Update PE array write cycle and energy.
                num_data_transfer[data_type_t::INPUT]++;

                std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);
                parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);

                uint64_t address_pe_array = 0, address_global_buffer = 0;
                unsigned num_access_array = 0, num_access_global_buffer = 0;
                for(unsigned b = 0; b < parameters_pe_array[parameter_type_t::BATCH_SIZE]; b++) {
                    for(unsigned c = 0; c < parameters_pe_array[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned h = 0; h < parameters_pe_array[parameter_type_t::INPUT_HEIGHT]; h++) {
                            for(unsigned w = 0; w < parameters_pe_array[parameter_type_t::INPUT_WIDTH]; w++) {
                                // Check temporal buffer address of PE array.
                                if(address_pe_array != ((uint64_t)&pe_array->input_data[b*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                         *parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                         *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                        c*parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                         *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                        h*parameters_pe_array[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                        pe_array->mask_bits[data_type_t::INPUT]) << pe_array->mask_bits[data_type_t::INPUT]) {

                                    // Update temporal buffer's costs of PE array
                                    pe_array->access_energy[data_type_t::INPUT] += pe_array->u_write_energy[data_type_t::INPUT];
                                    pe_array->access_cycle[data_type_t::INPUT] += pe_array->u_write_cycle[data_type_t::INPUT];
                                    num_access_array++;

                                    // Update temporal buffer address
                                    address_pe_array = ((uint64_t)&pe_array->input_data[b*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                         *parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                         *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                        c*parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                         *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                        h*parameters_pe_array[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                        pe_array->mask_bits[data_type_t::INPUT]) << pe_array->mask_bits[data_type_t::INPUT];
                                }
                                // Check global buffer address
                                if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::INPUT] + 
                                                                             m_scheduler->input_offset_global_buffer.front() +
                                                                             b*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] +
                                                                             c*parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                             h*parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                             mask_bits[data_type_t::INPUT] ) << mask_bits[data_type_t::INPUT]) {

                                    // Update global buffer cost
                                    access_energy[data_type_t::INPUT] += u_read_energy[data_type_t::INPUT];
                                    access_cycle[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT];
                                    num_access_global_buffer++;

                                    // Update global buffer address
                                    address_global_buffer = ((uint64_t)&data[offsets[data_type_t::INPUT] + 
                                                                             m_scheduler->input_offset_global_buffer.front() +
                                                                             b*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] +
                                                                             c*parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                             h*parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                             mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT];
                                }
                            }
                        }
                    }
                }

                // Update overlapped cost between PE array and global buffer
                unsigned ratio = ceil((double)(line_size[data_type_t::INPUT])/(double)(pe_array->line_size[data_type_t::INPUT]));

                // At the 1, 2, before last, and last stages
                unsigned first_stage = u_read_cycle[data_type_t::INPUT];
                unsigned second_stage = std::max(u_read_cycle[data_type_t::INPUT], 
                                                 u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)));

                unsigned last_before_stage = std::max(ratio*pe_array->u_write_cycle[data_type_t::INPUT],
                                                      u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)));
                unsigned last_stage = ratio*pe_array->u_write_cycle[data_type_t::INPUT];

                // Remainder stage
                unsigned other_stage = std::max(u_read_cycle[data_type_t::INPUT], 
                                       std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)),
                                                ratio*pe_array->u_write_cycle[data_type_t::INPUT]));

                if(num_access_global_buffer == 0) { } else if(num_access_global_buffer == 1) {
                    cycle_pe_array_global_buffer[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT] + 
                                                                        u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)) +
                                                                        pe_array->u_write_cycle[data_type_t::INPUT];
                } else {
                    cycle_pe_array_global_buffer[data_type_t::INPUT] += first_stage + second_stage +
                                                                        (num_access_global_buffer - 2)*other_stage +
                                                                        last_before_stage + last_stage;
                }

                // Update data transfer cycle and energy between Global buffer and PE array.
                transfer_cycle[data_type_t::INPUT] += num_access_global_buffer*u_transfer_cycle
                                                      *ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth));
                transfer_energy[data_type_t::INPUT] += num_access_global_buffer*u_transfer_energy
                                                      *ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth));

                pe_array->skip_transfer[data_type_t::INPUT] = false;

                //if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y -1) {
                //    move_front(&m_scheduler->input_offset_dram);
                //}
            }
        }
        // Case 2. COO data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
        }
        // Case 3. CSC data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::INPUT]) {
                // Column index calculation
                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::PE_Y);
                unsigned row = parameters[parameter_type_t::INPUT_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer[data_type_t::INPUT]++;

                // Update global buffer access cost
                access_cycle[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                    (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT]
                                                   /(sizeof(data_t)*8/row_bit) + // Row index
                                                    parameters[parameter_type_t::BATCH_SIZE]
                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                   *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                   *u_read_cycle[data_type_t::INPUT]
                                                   /(sizeof(data_t)*8/row_bit); // Column pointer
                access_energy[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT]  + // Non-zero data
                                                     (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT]/(sizeof(data_t)*8/row_bit) + // Row index
                                                     parameters[parameter_type_t::BATCH_SIZE]
                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                    *u_read_energy[data_type_t::INPUT]
                                                    /(sizeof(data_t)*8/row_bit); // Column pointer

                // Update PE array access cost (if temporal buffer exist)
                pe_array->access_cycle[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *pe_array->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                              (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *pe_array->u_write_cycle[data_type_t::INPUT]
                                                             /(sizeof(data_t)*8/row_bit) + // Row index
                                                              parameters[parameter_type_t::BATCH_SIZE]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                             *pe_array->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/row_bit); // Column pointer

                pe_array->access_energy[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                              *pe_array->u_write_energy[data_type_t::INPUT] + // Non-zero data
                                                               (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                              *pe_array->u_write_energy[data_type_t::INPUT]
                                                              /(sizeof(data_t)*8/row_bit) + // Row index
                                                               parameters[parameter_type_t::BATCH_SIZE]
                                                              *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                              *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                              *pe_array->u_write_energy[data_type_t::INPUT]
                                                              /(sizeof(data_t)*8/row_bit); // Column pointer

                // Update overlapped cycle between the Global buffer and temporal buffer in PE array
                cycle_pe_array_global_buffer[data_type_t::INPUT] += std::max((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *u_read_cycle[data_type_t::INPUT]
                                                                            /(sizeof(data_t)*8/row_bit) + // Row index
                                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                            *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                                            *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/row_bit), // Column pointer
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT]
                                                                            /(sizeof(data_t)*8/row_bit) + // Row index
                                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                            *(parameters[parameter_type_t::INPUT_WIDTH]+1)
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT]
                                                                            /(sizeof(data_t)*8/row_bit)); // Column pointer

                // Update transfer cost between the global buffer and temporal buffer in PE array
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                      u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)/(float)bitwidth) + // Row index
                                                      u_transfer_cycle*ceil((float)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                   *(parameters[parameter_type_t::INPUT_WIDTH+1])*row_bit)/(float)bitwidth); // Column pointer
                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*row_bit)/(float)bitwidth) + // Row index
                                                       u_transfer_energy*ceil((float)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                     *(parameters[parameter_type_t::INPUT_WIDTH+1])*row_bit)/(float)bitwidth); // Column pointer
                pe_array->skip_transfer[data_type_t::INPUT] = false;

            }
        }
        // Case 4. CSR data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::INPUT]) {
                // Column index calculation
                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::PE_Y);
                unsigned column = parameters[parameter_type_t::INPUT_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer[data_type_t::INPUT]++;

                // Update global buffer access cost
                access_cycle[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                    (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT]
                                                   /(sizeof(data_t)*8/column_bit) + // Column index
                                                    parameters[parameter_type_t::BATCH_SIZE]
                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                   *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                   *u_read_cycle[data_type_t::INPUT]
                                                   /(sizeof(data_t)*8/column_bit); // Row pointer
                access_energy[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT]  + // Non-zero data
                                                     (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT]/(sizeof(data_t)*8/column_bit) + // Column index
                                                     parameters[parameter_type_t::BATCH_SIZE]
                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                    *u_read_energy[data_type_t::INPUT]
                                                    /(sizeof(data_t)*8/column_bit); // Row pointer

                // Update PE array access cost (if temporal buffer exist)
                pe_array->access_cycle[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *pe_array->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                              (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *pe_array->u_write_cycle[data_type_t::INPUT]
                                                             /(sizeof(data_t)*8/column_bit) + // Column index
                                                              parameters[parameter_type_t::BATCH_SIZE]
                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                             *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                             *pe_array->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/column_bit); // Row pointer

                pe_array->access_energy[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                              *pe_array->u_write_energy[data_type_t::INPUT] + // Non-zero data
                                                               (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                              *pe_array->u_write_energy[data_type_t::INPUT]
                                                              /(sizeof(data_t)*8/column_bit) + // Column index
                                                               parameters[parameter_type_t::BATCH_SIZE]
                                                              *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                              *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                              *pe_array->u_write_energy[data_type_t::INPUT]
                                                              /(sizeof(data_t)*8/column_bit); // Row pointer

                // Update overlapped cycle between the Global buffer and temporal buffer in PE array
                cycle_pe_array_global_buffer[data_type_t::INPUT] += std::max((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *u_read_cycle[data_type_t::INPUT]
                                                                            /(sizeof(data_t)*8/column_bit) + // Column index
                                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                            *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                                            *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8/column_bit), // Row pointer
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT] + // Non-zero data
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT]
                                                                            /(sizeof(data_t)*8/column_bit) + // Column index
                                                                             parameters[parameter_type_t::BATCH_SIZE]
                                                                            *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                            *(parameters[parameter_type_t::INPUT_HEIGHT]+1)
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT]
                                                                            /(sizeof(data_t)*8/column_bit)); // Row pointer

                // Update transfer cost between the global buffer and temporal buffer in PE array
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                      u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)/(float)bitwidth) + // Column index
                                                      u_transfer_cycle*ceil((float)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                   *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                   *(parameters[parameter_type_t::INPUT_HEIGHT+1])*column_bit)/(float)bitwidth); // Row pointer
                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                       u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*column_bit)/(float)bitwidth) + // Column index
                                                       u_transfer_energy*ceil((float)(parameters[parameter_type_t::BATCH_SIZE]
                                                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                     *(parameters[parameter_type_t::INPUT_HEIGHT+1])*column_bit)/(float)bitwidth); // Row pointer
                pe_array->skip_transfer[data_type_t::INPUT] = false;

                //if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y -1) {
                //    move_front(&m_scheduler->input_offset_dram);
                //}
            }
        }
        // Case 4. SparseMap
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::INPUT]) {

                num_data_transfer[data_type_t::INPUT]++;
                // Update Global buffer read cycle and energy.

                // Update global buffer cost
                access_cycle[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                   *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                    pe_array->tile_size[data_type_t::INPUT]
                                                   *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
                access_energy[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                    *u_read_energy[data_type_t::INPUT] + // Non-zero data
                                                     pe_array->tile_size[data_type_t::INPUT]
                                                     *u_read_energy[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata

                // Update PE array cost
                pe_array->access_cycle[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                             *pe_array->u_write_cycle[data_type_t::INPUT] + // Non-zero value
                                                              pe_array->tile_size[data_type_t::INPUT]
                                                             *pe_array->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
                pe_array->access_energy[data_type_t::INPUT] += (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                              *pe_array->u_write_energy[data_type_t::INPUT] + // Non-zero data
                                                               pe_array->tile_size[data_type_t::INPUT]
                                                              *pe_array->u_write_energy[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata

                // Update overlapped cycle between the global buffer and temporal buffer in PE array
                cycle_pe_array_global_buffer[data_type_t::INPUT] += std::max((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *u_read_cycle[data_type_t::INPUT] + // Non-zero data
                                                                             pe_array->tile_size[data_type_t::INPUT]
                                                                            *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8), // Metadata
                                                                             (pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT] + // Non-zero data 
                                                                             pe_array->tile_size[data_type_t::INPUT]
                                                                            *pe_array->u_write_cycle[data_type_t::INPUT]/(sizeof(data_t)*8)); // Metadata
               
                // Update transfer cost between the global buffer and temporal buffer in PE array
                transfer_cycle[data_type_t::INPUT] += u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                      u_transfer_cycle*ceil((float)(pe_array->tile_size[data_type_t::INPUT])/(float)bitwidth); // Meta data

                transfer_energy[data_type_t::INPUT] += u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::INPUT] - m_scheduler->num_zeros[data_type_t::INPUT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                     + u_transfer_energy*ceil((float)(pe_array->tile_size[data_type_t::INPUT])/(float)bitwidth); // Meta data

                pe_array->skip_transfer[data_type_t::INPUT] = false;

                //if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                //    move_front(&m_scheduler->input_offset_global_buffer);
                //}
            }
        }
        if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y -1) {
            move_front(&m_scheduler->input_offset_dram);
        }
#else
#if 0
        /* Legacy host-address timing path: replaced by descriptor accounting above. */
        if(!skip_transfer[data_type_t::INPUT]) {
            // Update PE array write cycle and energy.
            num_data_transfer[data_type_t::INPUT]++;

            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);
            parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);

            uint64_t address_pe_array = 0, address_global_buffer = 0;
            unsigned num_access_array = 0, num_access_global_buffer = 0;
            for(unsigned b = 0; b < parameters_pe_array[parameter_type_t::BATCH_SIZE]; b++) {
                for(unsigned c = 0; c < parameters_pe_array[parameter_type_t::INPUT_CHANNEL]; c++) {
                    for(unsigned h = 0; h < parameters_pe_array[parameter_type_t::INPUT_HEIGHT]; h++) {
                        for(unsigned w = 0; w < parameters_pe_array[parameter_type_t::INPUT_WIDTH]; w++) {
                            // Check temporal buffer address of PE array.
                            if(address_pe_array != ((uint64_t)&pe_array->input_data[b*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                     *parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                    c*parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                    h*parameters_pe_array[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                    pe_array->mask_bits[data_type_t::INPUT]) << pe_array->mask_bits[data_type_t::INPUT]) {

                                // Update temporal buffer's costs of PE array
                                pe_array->access_energy[data_type_t::INPUT] += pe_array->u_write_energy[data_type_t::INPUT];
                                pe_array->access_cycle[data_type_t::INPUT] += pe_array->u_write_cycle[data_type_t::INPUT];
                                num_access_array++;

                                // Update temporal buffer address
                                address_pe_array = ((uint64_t)&pe_array->input_data[b*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                     *parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                    c*parameters_pe_array[parameter_type_t::INPUT_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::INPUT_WIDTH] +
                                                                                    h*parameters_pe_array[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                                    pe_array->mask_bits[data_type_t::INPUT]) << pe_array->mask_bits[data_type_t::INPUT];
                            }
                            // Check global buffer address
                            if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::INPUT] + 
                                                                         m_scheduler->input_offset_global_buffer.front() +
                                                                         b*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] +
                                                                         c*parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                         h*parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                         mask_bits[data_type_t::INPUT] ) << mask_bits[data_type_t::INPUT]) {

                                // Update global buffer cost
                                access_energy[data_type_t::INPUT] += u_read_energy[data_type_t::INPUT];
                                access_cycle[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT];
                                num_access_global_buffer++;

                                // Update global buffer address
                                address_global_buffer = ((uint64_t)&data[offsets[data_type_t::INPUT] + 
                                                                         m_scheduler->input_offset_global_buffer.front() +
                                                                         b*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] +
                                                                         c*parameters_global_buffer[parameter_type_t::INPUT_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + 
                                                                         h*parameters_global_buffer[parameter_type_t::INPUT_WIDTH] + w] >> 
                                                                         mask_bits[data_type_t::INPUT]) << mask_bits[data_type_t::INPUT];
                            }
                        }
                    }
                }
            }

            // Update overlapped cost between PE array and global buffer
            unsigned ratio = ceil((double)(line_size[data_type_t::INPUT])/(double)(pe_array->line_size[data_type_t::INPUT]));

            // At the 1, 2, before last, and last stages
            unsigned first_stage = u_read_cycle[data_type_t::INPUT];
            unsigned second_stage = std::max(u_read_cycle[data_type_t::INPUT], 
                                             u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)));

            unsigned last_before_stage = std::max(ratio*pe_array->u_write_cycle[data_type_t::INPUT],
                                                  u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)));
            unsigned last_stage = ratio*pe_array->u_write_cycle[data_type_t::INPUT];

            // Remainder stage
            unsigned other_stage = std::max(u_read_cycle[data_type_t::INPUT], 
                                   std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)),
                                            ratio*pe_array->u_write_cycle[data_type_t::INPUT]));

            if(num_access_global_buffer == 0) { } else if(num_access_global_buffer == 1) {
                cycle_pe_array_global_buffer[data_type_t::INPUT] += u_read_cycle[data_type_t::INPUT] + 
                                                                    u_transfer_cycle*ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth)) +
                                                                    pe_array->u_write_cycle[data_type_t::INPUT];
            } else {
                cycle_pe_array_global_buffer[data_type_t::INPUT] += first_stage + second_stage +
                                                                    (num_access_global_buffer - 2)*other_stage +
                                                                    last_before_stage + last_stage;
            }

            // Update data transfer cycle and energy between Global buffer and PE array.
            transfer_cycle[data_type_t::INPUT] += num_access_global_buffer*u_transfer_cycle
                                                  *ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth));
            transfer_energy[data_type_t::INPUT] += num_access_global_buffer*u_transfer_energy
                                                  *ceil((double)(line_size[data_type_t::INPUT])/(double)(bitwidth));

            pe_array->skip_transfer[data_type_t::INPUT] = false;
            if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                move_front(&m_scheduler->input_offset_global_buffer);
            }
        }
#endif
#endif

        // Increment the input index that finds the offset of input data in Global buffer.
        input_index++;
        // Input data exists in PE array.
        // Request input data to Global buffer is not required.
        pe_array->exist_data[data_type_t::INPUT] = true, pe_array->request_to_global_buffer[data_type_t::INPUT] = false;
        if(pe_array->tile_size[data_type_t::INPUT] == tile_size[data_type_t::INPUT]) { skip_transfer[data_type_t::INPUT] = true; }
    }
    // Transfer weight from Global buffer to temporal buffer of PE array.
    if(pe_array->request_to_global_buffer[data_type_t::WEIGHT]) {
#ifndef FUNCTIONAL
        if(m_scheduler->compression_type == compression_type_t::DENSE && !skip_transfer[data_type_t::WEIGHT]) {
            account_descriptor_dense_transfer(data_type_t::WEIGHT);
        }
#endif
#ifdef FUNCTIONAL
        // Weight transfer
        m_scheduler->transfer_data(pe_array->weight, data, 0, offsets[data_type_t::WEIGHT] + m_scheduler->weight_offset_global_buffer.front(),
                                   component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, 
                                   data_type_t::WEIGHT, pe_array->get_stationary_type(), action_type_t::LOAD);
        // Update for NPUsim ver2
        //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
        //m_scheduler->transfer_data_ver2(pe_array->weight, data, 
        //                                component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, 
        //                                data_type_t::WEIGHT, pe_array->get_stationary_type(), action_type_t::LOAD, true);
        
        // Case 1. Dense
        if(m_scheduler->compression_type == compression_type_t::DENSE) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                num_data_transfer[data_type_t::WEIGHT]++;

                std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);
                parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);

                uint64_t address_pe_array = 0, address_global_buffer = 0;
                unsigned num_access_pe_array = 0, num_access_global_buffer = 0;

                for(unsigned k = 0; k < parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                    for(unsigned c = 0; c < parameters_pe_array[parameter_type_t::INPUT_CHANNEL]; c++) {
                        for(unsigned r = 0; r < parameters_pe_array[parameter_type_t::FILTER_HEIGHT]; r++) {
                            for(unsigned s = 0; s < parameters_pe_array[parameter_type_t::FILTER_WIDTH]; s++) {
                                // Check temporal buffer address at PE array
                                if(address_pe_array != ((uint64_t)&pe_array->weight[k*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                     *parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                                    c*parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::FILTER_WIDTH] +
                                                                                    r*parameters_pe_array[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                    pe_array->mask_bits[data_type_t::WEIGHT]) << pe_array->mask_bits[data_type_t::WEIGHT]) {

                                    // Update PE array cost
                                    pe_array->access_energy[data_type_t::WEIGHT] += pe_array->u_write_energy[data_type_t::WEIGHT];
                                    pe_array->access_cycle[data_type_t::WEIGHT] += pe_array->u_write_cycle[data_type_t::WEIGHT];
                                    num_access_pe_array++;

                                    // Update PE array address (temporal buffer)
                                    address_pe_array = ((uint64_t)&pe_array->weight[k*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                     *parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                                    c*parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                     *parameters_pe_array[parameter_type_t::FILTER_WIDTH] +
                                                                                    r*parameters_pe_array[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                    pe_array->mask_bits[data_type_t::WEIGHT]) << pe_array->mask_bits[data_type_t::WEIGHT];
                                }

                                // Check global buffer address
                                if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::WEIGHT] + 
                                                                             m_scheduler->weight_offset_global_buffer.front() +
                                                                             k*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] +
                                                                             c*parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                             r*parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                             mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT]) {

                                    // Update global buffer cost
                                    access_energy[data_type_t::WEIGHT] += u_read_energy[data_type_t::WEIGHT];
                                    access_cycle[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT];
                                    num_access_global_buffer++;

                                    // Update Global buffer address
                                    address_global_buffer = ((uint64_t)&data[offsets[data_type_t::WEIGHT] + 
                                                                             m_scheduler->weight_offset_global_buffer.front() +
                                                                             k*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                              *parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] +
                                                                             c*parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                             r*parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                             mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT];
                                }
                            }
                        }
                    }
                }

                // Update overlapped cost between PE array and global buffer
                unsigned ratio = ceil((double)(line_size[data_type_t::WEIGHT])/(double)(pe_array->line_size[data_type_t::WEIGHT]));

                // At the 1, 2, before last, and last stages
                unsigned first_stage = u_read_cycle[data_type_t::WEIGHT];
                unsigned second_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                                 u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
                unsigned last_before_stage = std::max(ratio*pe_array->u_write_cycle[data_type_t::WEIGHT],
                                                      u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
                unsigned last_stage = ratio*pe_array->u_write_cycle[data_type_t::WEIGHT];

                // Remainder stages
                unsigned other_stage = std::max(u_read_cycle[data_type_t::WEIGHT], 
                                       std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)),
                                                ratio*pe_array->u_write_cycle[data_type_t::WEIGHT]));

                if(num_access_global_buffer == 0) { } else if(num_access_global_buffer == 1) {
                    cycle_pe_array_global_buffer[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT] +   
                                                                         u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)) +
                                                                         ratio*pe_array->u_write_cycle[data_type_t::WEIGHT];
                } else {
                    cycle_pe_array_global_buffer[data_type_t::WEIGHT] += first_stage + second_stage +
                                                                         (num_access_global_buffer - 2)*other_stage +
                                                                         last_before_stage + last_stage;
                }

                // Update data transfer cycle and energy between Global buffer and PE array.
                transfer_cycle[data_type_t::WEIGHT] += num_access_global_buffer*u_transfer_cycle
                                                       *ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth));
                transfer_energy[data_type_t::WEIGHT] += num_access_global_buffer*u_transfer_energy
                                                        *ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth));
                pe_array->skip_transfer[data_type_t::WEIGHT] = false;

                if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                    move_front(&m_scheduler->weight_offset_global_buffer);
                }
            }
        }
        // Case 2. COO data format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_COO) {
        }
        // Case 3. CSC format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSC) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                // Row bit calculation
                unsigned row_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::PE_Y);
                unsigned row = parameters[parameter_type_t::FILTER_HEIGHT];
                while(row > 1) {
                    row /= 2;
                    row_bit++;
                }

                num_data_transfer[data_type_t::WEIGHT]++;

                // Update global buffer cost
                access_cycle[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                     (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT]/(sizeof(data_t)*8/row_bit) + // Row index
                                                     parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                    *u_read_cycle[data_type_t::WEIGHT]
                                                    /(sizeof(data_t)*8/row_bit); // Column pointer
                access_energy[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                      (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT]
                                                     /(sizeof(data_t)*8/row_bit) + // Row index
                                                      parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                     *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                     *u_read_energy[data_type_t::WEIGHT]
                                                     /(sizeof(data_t)*8/row_bit); // Column pointer
                                                     
                // Update PE array cost (if temporal buffer exist)
                pe_array->access_cycle[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                               (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                              /(sizeof(data_t)*8/row_bit) + // Row index
                                                               parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                              *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                              *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                              /(sizeof(data_t)*8/row_bit); // Column pointer
                pe_array->access_energy[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT] + // Non-zero data
                                                                (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT]
                                                               /(sizeof(data_t)*8/row_bit) + // Row index
                                                                parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT]
                                                               /(sizeof(data_t)*8/row_bit); // Column pointer

                // Update overlapped cycle between the global buffer and PE array
                cycle_pe_array_global_buffer[data_type_t::WEIGHT] += std::max((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                              (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *u_read_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/row_bit) + // Row index
                                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                                             *u_read_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/row_bit),  // Column pointer
                                                                              (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *pe_array->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                              (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/row_bit) + // Row index
                                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *(parameters[parameter_type_t::FILTER_WIDTH]+1)
                                                                             *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/row_bit)); // Column pointer

                // Update transfer cost between the global buffer and PE array
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                       u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*row_bit)/(float)bitwidth) + // Row index
                                                       u_transfer_cycle*ceil((float)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *(parameters[parameter_type_t::FILTER_WIDTH+1])*row_bit)/(float)bitwidth); // Column pointer
                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                        u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*row_bit)/(float)bitwidth) + // Row index
                                                        u_transfer_energy*ceil((float)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                      *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                      *(parameters[parameter_type_t::FILTER_WIDTH+1])*row)/(float)bitwidth); // Column pointer
                pe_array->skip_transfer[data_type_t::WEIGHT] = false;
                if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                    move_front(&m_scheduler->weight_offset_global_buffer);
                }
            }
        }
        // Case 4. CSR format
        else if(m_scheduler->compression_type == compression_type_t::SPARSE_CSR) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                // Column bit calculation
                unsigned column_bit = 1;
                std::vector<unsigned> parameters = m_scheduler->calculate_parameter_size(component_type_t::PE_Y);
                unsigned column = parameters[parameter_type_t::FILTER_WIDTH];
                while(column > 1) {
                    column /= 2;
                    column_bit++;
                }

                num_data_transfer[data_type_t::WEIGHT]++;

                // Update global buffer cost
                access_cycle[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                     (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT]/(sizeof(data_t)*8/column_bit) + // Column index
                                                     parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                    *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                    *u_read_cycle[data_type_t::WEIGHT]
                                                    /(sizeof(data_t)*8/column_bit); // Row pointer
                access_energy[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                      (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT]
                                                     /(sizeof(data_t)*8/column_bit) + // Column index
                                                      parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                     *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                     *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                     *u_read_energy[data_type_t::WEIGHT]
                                                     /(sizeof(data_t)*8/column_bit); // Row pointer
                                                     
                // Update PE array cost (if temporal buffer exist)
                pe_array->access_cycle[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                               (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                              /(sizeof(data_t)*8/column_bit) + // Column index
                                                               parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                              *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                              *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                              /(sizeof(data_t)*8/column_bit); // Row pointer
                pe_array->access_energy[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT] + // Non-zero data
                                                                (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT]
                                                               /(sizeof(data_t)*8/column_bit) + // Column index
                                                                parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                               *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT]
                                                               /(sizeof(data_t)*8/column_bit); // Row pointer

                // Update overlapped cycle between the global buffer and PE array
                cycle_pe_array_global_buffer[data_type_t::WEIGHT] += std::max((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                              (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *u_read_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/column_bit) + // Column index
                                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                                             *u_read_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/column_bit),  // Row pointer
                                                                              (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *pe_array->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                              (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                             *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/column_bit) + // Column index
                                                                              parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                             *(parameters[parameter_type_t::FILTER_HEIGHT]+1)
                                                                             *pe_array->u_write_cycle[data_type_t::WEIGHT]
                                                                             /(sizeof(data_t)*8/column_bit)); // row pointer

                // Update transfer cost between the global buffer and PE array
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                       u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*column_bit)/(float)bitwidth) + // Column index
                                                       u_transfer_cycle*ceil((float)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                    *(parameters[parameter_type_t::FILTER_HEIGHT+1])*column_bit)/(float)bitwidth); // Row pointer
                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                        u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*column_bit)/(float)bitwidth) + // Column index
                                                        u_transfer_energy*ceil((float)(parameters[parameter_type_t::OUTPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                      *parameters[parameter_type_t::INPUT_CHANNEL]/parameters[parameter_type_t::GROUP]
                                                                                      *(parameters[parameter_type_t::FILTER_HEIGHT+1])*column_bit)/(float)bitwidth); // Row pointer
                pe_array->skip_transfer[data_type_t::WEIGHT] = false;
                if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                    move_front(&m_scheduler->weight_offset_global_buffer);
                }
            }
        }
        // Case 5. SparseMap
        else if(m_scheduler->compression_type == compression_type_t::SPARSEMAP) {
            if(!skip_transfer[data_type_t::WEIGHT]) {
                num_data_transfer[data_type_t::WEIGHT]++;

                // Update Global buffer access cost
                access_cycle[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                    *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                     pe_array->tile_size[data_type_t::WEIGHT]
                                                    *u_read_cycle[data_type_t::INPUT]/(sizeof(data_t)*8); // Metadata
                access_energy[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                     *u_read_energy[data_type_t::WEIGHT] + // Non-zero data
                                                      pe_array->tile_size[data_type_t::WEIGHT]
                                                     *u_read_energy[data_type_t::WEIGHT]/(sizeof(data_t)*8); // Metadata

                // Update PE array access cost (if temporal buffer exist)
                pe_array->access_cycle[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                               pe_array->tile_size[data_type_t::WEIGHT]
                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT]/(sizeof(data_t)*8); // Metadata
                pe_array->access_energy[data_type_t::WEIGHT] += (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT] + // Non-zero data
                                                                pe_array->tile_size[data_type_t::WEIGHT]
                                                               *pe_array->u_write_energy[data_type_t::WEIGHT]/(sizeof(data_t)*8); // Metadata

                // Update overlapped cycle between the global buffer and PE array
                cycle_pe_array_global_buffer[data_type_t::WEIGHT] += std::max((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                              *u_read_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                               pe_array->tile_size[data_type_t::WEIGHT]
                                                                              *u_read_cycle[data_type_t::WEIGHT]/(sizeof(data_t)*8), // Metadata
                                                                               (pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])
                                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT] + // Non-zero data
                                                                               pe_array->tile_size[data_type_t::WEIGHT]
                                                                              *pe_array->u_write_cycle[data_type_t::WEIGHT]/(sizeof(data_t)*8)); // Metadata


                // Update transfer cost between the global buffer and PE array
                transfer_cycle[data_type_t::WEIGHT] += u_transfer_cycle*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                       u_transfer_cycle*ceil((float)(pe_array->tile_size[data_type_t::WEIGHT])/(float)bitwidth); // Metadata

                transfer_energy[data_type_t::WEIGHT] += u_transfer_energy*ceil((float)((pe_array->tile_size[data_type_t::WEIGHT] - m_scheduler->num_zeros[data_type_t::WEIGHT])*8*sizeof(data_t))/(float)bitwidth) + // Non-zero data
                                                        u_transfer_energy*ceil((float)(pe_array->tile_size[data_type_t::WEIGHT])/(float)bitwidth); // Meta data

                pe_array->skip_transfer[data_type_t::WEIGHT] = false;
                if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                    move_front(&m_scheduler->weight_offset_global_buffer);
                }
            }
        }
        else {
            std::cerr << "Undefined compressed data type" << std::endl;
            exit(1);
        }

        /* Stats */

#else
#if 0
        /* Legacy host-address timing path: replaced by descriptor accounting above. */
        if(!skip_transfer[data_type_t::WEIGHT]) {
            num_data_transfer[data_type_t::WEIGHT]++;

            std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);
            std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);

            parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);
            parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);

            uint64_t address_pe_array = 0, address_global_buffer = 0;
            unsigned num_access_pe_array = 0, num_access_global_buffer = 0;

            for(unsigned k = 0; k < parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                for(unsigned c = 0; c < parameters_pe_array[parameter_type_t::INPUT_CHANNEL]; c++) {
                    for(unsigned r = 0; r < parameters_pe_array[parameter_type_t::FILTER_HEIGHT]; r++) {
                        for(unsigned s = 0; s < parameters_pe_array[parameter_type_t::FILTER_WIDTH]; s++) {
                            // Check temporal buffer address at PE array
                            if(address_pe_array != ((uint64_t)&pe_array->weight[k*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                 *parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                 *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                                c*parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                 *parameters_pe_array[parameter_type_t::FILTER_WIDTH] +
                                                                                r*parameters_pe_array[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                pe_array->mask_bits[data_type_t::WEIGHT]) << pe_array->mask_bits[data_type_t::WEIGHT]) {

                                // Update PE array cost
                                pe_array->access_energy[data_type_t::WEIGHT] += pe_array->u_write_energy[data_type_t::WEIGHT];
                                pe_array->access_cycle[data_type_t::WEIGHT] += pe_array->u_write_cycle[data_type_t::WEIGHT];
                                num_access_pe_array++;

                                // Update PE array address (temporal buffer)
                                address_pe_array = ((uint64_t)&pe_array->weight[k*parameters_pe_array[parameter_type_t::INPUT_CHANNEL]
                                                                                 *parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                 *parameters_pe_array[parameter_type_t::FILTER_WIDTH] + 
                                                                                c*parameters_pe_array[parameter_type_t::FILTER_HEIGHT]
                                                                                 *parameters_pe_array[parameter_type_t::FILTER_WIDTH] +
                                                                                r*parameters_pe_array[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                                pe_array->mask_bits[data_type_t::WEIGHT]) << pe_array->mask_bits[data_type_t::WEIGHT];
                            }

                            // Check global buffer address
                            if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::WEIGHT] + 
                                                                         m_scheduler->weight_offset_global_buffer.front() +
                                                                         k*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] +
                                                                         c*parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                         r*parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                         mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT]) {

                                // Update global buffer cost
                                access_energy[data_type_t::WEIGHT] += u_read_energy[data_type_t::WEIGHT];
                                access_cycle[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT];
                                num_access_global_buffer++;

                                // Update Global buffer address
                                address_global_buffer = ((uint64_t)&data[offsets[data_type_t::WEIGHT] + 
                                                                         m_scheduler->weight_offset_global_buffer.front() + 
                                                                         k*parameters_global_buffer[parameter_type_t::INPUT_CHANNEL]
                                                                          *parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] +
                                                                         c*parameters_global_buffer[parameter_type_t::FILTER_HEIGHT]
                                                                          *parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + 
                                                                         r*parameters_global_buffer[parameter_type_t::FILTER_WIDTH] + s] >> 
                                                                         mask_bits[data_type_t::WEIGHT]) << mask_bits[data_type_t::WEIGHT];
                            }
                        }
                    }
                }
            }

            // Update overlapped cost between PE array and global buffer
            unsigned ratio = ceil((double)(line_size[data_type_t::WEIGHT])/(double)(pe_array->line_size[data_type_t::WEIGHT]));

            // At the 1, 2, before last, and last stages
            unsigned first_stage = u_read_cycle[data_type_t::WEIGHT];
            unsigned second_stage = std::max(u_read_cycle[data_type_t::WEIGHT],
                                             u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
            unsigned last_before_stage = std::max(ratio*pe_array->u_write_cycle[data_type_t::WEIGHT],
                                                  u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)));
            unsigned last_stage = ratio*pe_array->u_write_cycle[data_type_t::WEIGHT];

            // Remainder stages
            unsigned other_stage = std::max(u_read_cycle[data_type_t::WEIGHT], 
                                   std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)),
                                            ratio*pe_array->u_write_cycle[data_type_t::WEIGHT]));

            if(num_access_global_buffer == 0) { } else if(num_access_global_buffer == 1) {
                cycle_pe_array_global_buffer[data_type_t::WEIGHT] += u_read_cycle[data_type_t::WEIGHT] +   
                                                                     u_transfer_cycle*ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth)) +
                                                                     ratio*pe_array->u_write_cycle[data_type_t::WEIGHT];
            } else {
                cycle_pe_array_global_buffer[data_type_t::WEIGHT] += first_stage + second_stage +
                                                                     (num_access_global_buffer - 2)*other_stage +
                                                                     last_before_stage + last_stage;
            }

            // Update data transfer cycle and energy between Global buffer and PE array.
            transfer_cycle[data_type_t::WEIGHT] += num_access_global_buffer*u_transfer_cycle
                                                   *ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth));
            transfer_energy[data_type_t::WEIGHT] += num_access_global_buffer*u_transfer_energy
                                                    *ceil((double)(line_size[data_type_t::WEIGHT])/(double)(bitwidth));

            pe_array->skip_transfer[data_type_t::WEIGHT] = false;

            if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y - 1) {
                move_front(&m_scheduler->weight_offset_global_buffer);
            }
        }
#endif
#endif
        /* Stats */

        // Increase the weight index by one.
        weight_index++;

        // Update PE array signals
        pe_array->exist_data[data_type_t::WEIGHT] = true, pe_array->request_to_global_buffer[data_type_t::WEIGHT] = false;
        if(pe_array->tile_size[data_type_t::WEIGHT] == tile_size[data_type_t::WEIGHT]) { skip_transfer[data_type_t::WEIGHT] = true; }
    }
    // Transfer output data from Global buffer to temporal buffer of PE array.
    if(pe_array->request_to_global_buffer[data_type_t::OUTPUT]) {
        if(m_scheduler->output_read_global_buffer[m_scheduler->output_offset_global_buffer.front()]) {
#ifdef FUNCTIONAL
            // Output data transfer
            m_scheduler->transfer_data(pe_array->output_data, data, 0, offsets[data_type_t::OUTPUT] + m_scheduler->output_offset_global_buffer.front(), 
                                       component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, 
                                       data_type_t::OUTPUT, pe_array->get_stationary_type(), action_type_t::LOAD);

            // Update for NPUsim ver2
            //bool last_component = index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y;
            //m_scheduler->transfer_data_ver2(pe_array->output_data, data, 
            //                                component_type_t::PE_Y, component_type_t::GLOBAL_BUFFER, 
            //                                data_type_t::OUTPUT, pe_array->get_stationary_type(), action_type_t::LOAD, true);
#endif
#ifndef FUNCTIONAL
            if(!skip_transfer[data_type_t::OUTPUT]) {
                account_descriptor_dense_transfer(data_type_t::OUTPUT);
                pe_array->equal_output_tile = false;
                transfer_output = true;
            }
#else
            if(!skip_transfer[data_type_t::OUTPUT]) {
                num_data_transfer[data_type_t::OUTPUT]++;

                std::vector<unsigned> parameters_pe_array(parameter_type_t::NUM_PARAMETER_TYPES, 1);
                std::vector<unsigned> parameters_global_buffer(parameter_type_t::NUM_PARAMETER_TYPES, 1);

                parameters_pe_array = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::PE_Y);
                parameters_global_buffer = m_scheduler->mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);

                uint64_t address_pe_array = 0, address_global_buffer = 0;
                unsigned num_access_pe_array = 0, num_access_global_buffer = 0;

                for(unsigned b = 0; b < parameters_pe_array[parameter_type_t::BATCH_SIZE]; b++) {
                    for(unsigned k = 0; k < parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                        for(unsigned p = 0; p < parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                            for(unsigned q = 0; q < parameters_pe_array[parameter_type_t::OUTPUT_WIDTH]; q++) {
                                // Check address at PE array (temporal buffer)
                                if(address_pe_array != ((uint64_t)&pe_array->output_data[b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                          *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                          *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                         k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                          *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                         p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                         pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT]) {

                                    // Update PE array cost
                                    pe_array->access_energy[data_type_t::OUTPUT] += pe_array->u_write_energy[data_type_t::OUTPUT];
                                    pe_array->access_cycle[data_type_t::OUTPUT] += pe_array->u_write_cycle[data_type_t::OUTPUT];
                                    num_access_pe_array++;

                                    // Update PE array address (temporal buffer)
                                    address_pe_array = ((uint64_t)&pe_array->output_data[b*parameters_pe_array[parameter_type_t::OUTPUT_CHANNEL]
                                                                                          *parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                          *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                         k*parameters_pe_array[parameter_type_t::OUTPUT_HEIGHT]
                                                                                          *parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + 
                                                                                         p*parameters_pe_array[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                                         pe_array->mask_bits[data_type_t::OUTPUT]) << pe_array->mask_bits[data_type_t::OUTPUT];
                                }

                                // Check global buffer address 
                                if(address_global_buffer != ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                             m_scheduler->output_offset_global_buffer.front() + 
                                                                             b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                              *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                             k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                             p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                             mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT]) {
                                    // Update global buffer cost
                                    access_energy[data_type_t::OUTPUT] += u_read_energy[data_type_t::OUTPUT];
                                    access_cycle[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT];
                                    num_access_global_buffer++;

                                    // Update global buffer address 
                                    address_global_buffer = ((uint64_t)&data[offsets[data_type_t::OUTPUT] + 
                                                                             m_scheduler->output_offset_global_buffer.front() + 
                                                                             b*parameters_global_buffer[parameter_type_t::OUTPUT_CHANNEL]
                                                                              *parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                             k*parameters_global_buffer[parameter_type_t::OUTPUT_HEIGHT]
                                                                              *parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + 
                                                                             p*parameters_global_buffer[parameter_type_t::OUTPUT_WIDTH] + q] >> 
                                                                             mask_bits[data_type_t::OUTPUT]) << mask_bits[data_type_t::OUTPUT];
                                }
                            }
                        }
                    }
                }
                // Update overlapped cost between PE array and global buffer
                unsigned ratio = ceil((double)(line_size[data_type_t::OUTPUT])/(double)(pe_array->line_size[data_type_t::OUTPUT]));

                // At the 1, 2, before last, and last stages
                unsigned first_stage = u_read_cycle[data_type_t::OUTPUT];
                unsigned second_stage = std::max(u_read_cycle[data_type_t::OUTPUT],
                                                 u_transfer_cycle*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth)));
                unsigned last_before_stage = std::max(ratio*pe_array->u_write_cycle[data_type_t::OUTPUT],
                                                      u_transfer_cycle*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth)));
                unsigned last_stage = ratio*pe_array->u_write_cycle[data_type_t::OUTPUT];

                // Remainder stages
                unsigned other_stage = std::max(u_read_cycle[data_type_t::OUTPUT],
                                       std::max(u_transfer_cycle*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth)),
                                                ratio*pe_array->u_write_cycle[data_type_t::OUTPUT]));

                if(num_access_global_buffer == 0) { } else if(num_access_global_buffer == 1) {
                    cycle_pe_array_global_buffer[data_type_t::OUTPUT] += u_read_cycle[data_type_t::OUTPUT] +
                                                                         u_transfer_cycle*ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth)) +
                                                                         ratio*pe_array->u_write_cycle[data_type_t::OUTPUT];
                } else {
                    cycle_pe_array_global_buffer[data_type_t::OUTPUT] += first_stage + second_stage +
                                                                         (num_access_global_buffer - 2)*other_stage +
                                                                         last_before_stage + last_stage;
                }

                // Update data transfer cycle and energy between Global buffer and PE array.
                transfer_cycle[data_type_t::OUTPUT] += num_access_global_buffer*u_transfer_cycle
                                                       *ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth));
                transfer_energy[data_type_t::OUTPUT] += num_access_global_buffer*u_transfer_energy
                                                        *ceil((double)(line_size[data_type_t::OUTPUT])/(double)(bitwidth));


                pe_array->skip_transfer[data_type_t::OUTPUT] = false;
                pe_array->equal_output_tile = false;
                transfer_output = true;
            }
#endif
            // E20-3b: the psum LOAD and the psum STORE are one decision -- does the array's
            // output tile leave and come back? -- so this mirrors the store-side verdict rather
            // than deriving a second one. It used to LATCH on a bare tile-size comparison and
            // never reset within a layer, so one early match suppressed every later reload while
            // the store side kept writing psums out (Eyeriss conv3: 19656 stores against 312
            // loads). A tile-size comparison is also weaker than pe_array_t's retention test,
            // which additionally consults the multi-chip/DRAM output tiles, the weight tile and
            // psum_retention_valid; the two disagreed outright on the 512x512x512 GEMM (96
            // stores, 0 loads). A psum written out and never read back cannot be accumulated.
            skip_transfer[data_type_t::OUTPUT] = pe_array->equal_output_tile;
        }
        else {
            // E20-3b: a first touch means the array is being pointed at a DIFFERENT output tile,
            // so whatever it was holding is finished and has to be written out. Release the
            // residency latch here, mirroring what dram_t does for the multi-chip buffer. Without
            // this the latch, once set, survived the rest of the layer: with no reload to clear it
            // (the load path is the only other place that does) a retained array wrote its FIRST
            // finished tile back and then silently kept every later one.
            transfer_output = false;
            pe_array->equal_output_tile = false;
            m_scheduler->output_read_global_buffer[m_scheduler->output_offset_global_buffer.front()] = true;
            for(auto it = m_scheduler->output_read_pe.begin(); it != m_scheduler->output_read_pe.end(); ++it) {
                it->second = false;
            }
        }

        if(index == m_scheduler->num_active_chips_x*m_scheduler->num_active_chips_y -1) {
            move_front(&m_scheduler->output_offset_global_buffer);
        }
        // Increment the output index.
        output_index++;

        // Update PE array signals
        pe_array->exist_data[data_type_t::OUTPUT] = true, pe_array->request_to_global_buffer[data_type_t::OUTPUT] = false;
    }
    pe_array->fill_data();

    // Check iterations of data
    // Case 1. Input stationary
    if(pe_array->get_stationary_type() == stationary_type_t::INPUT_STATIONARY) {
        // Reuse input data
        if(weight_index == m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front() &&
           output_index == m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front() &&
           input_index < m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front()) {

            weight_index = 0, output_index = 0;

            move_front(&m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT]);
            move_front(&m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT]);
        }
    }
    // Case 2. Weight stationary
    else if(pe_array->get_stationary_type() == stationary_type_t::WEIGHT_STATIONARY) {
        // Reuse weight
        if(input_index == m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front() &&
           output_index == m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front() &&
           weight_index < m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front()) {

            input_index = 0, output_index = 0;

            move_front(&m_scheduler->offset_size_global_buffer[data_type_t::INPUT]);
            move_front(&m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT]);
        }
    }
    // Case 3. Output stationary
    else if(pe_array->get_stationary_type() == stationary_type_t::OUTPUT_STATIONARY) {
        // Reuse output data
        if(input_index == m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front() &&
           weight_index == m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front() &&
           output_index < m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front()) {

            input_index = 0, weight_index = 0;

            move_front(&m_scheduler->offset_size_global_buffer[data_type_t::INPUT]);
            move_front(&m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT]);
        }
    }
    // Refresh all data
    if(input_index == m_scheduler->offset_size_global_buffer[data_type_t::INPUT].front() &&
       weight_index == m_scheduler->offset_size_global_buffer[data_type_t::WEIGHT].front() &&
       output_index == m_scheduler->offset_size_global_buffer[data_type_t::OUTPUT].front()) {
        flush_data(m_scheduler);
    }
}


void global_buffer_t::flush_data(scheduler_t *m_scheduler) {
    // Case 1. Input stationary
    if(stationary_type == stationary_type_t::INPUT_STATIONARY) {
        // Case 1) Can reuse input data && request weight and output data
        if(weight_flush_counter < m_scheduler->offset_size_dram[data_type_t::WEIGHT].front() - 1 &&
           output_flush_counter < m_scheduler->offset_size_dram[data_type_t::OUTPUT].front() - 1) {
            // Weight and output data do not exist in Global buffer.
            exist_data[data_type_t::WEIGHT] = false, exist_data[data_type_t::OUTPUT] = false;

#ifdef FUNCTIONAL
            // Write back Output data 
            m_scheduler->transfer_data(multi_chip->data, data, m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()], offsets[data_type_t::OUTPUT],
                                       component_type_t::CHIPS_Y, component_type_t::GLOBAL_BUFFER, 
                                       data_type_t::OUTPUT, get_stationary_type(), action_type_t::STORE);
#endif
            // GB1+GB2: descriptor-based GLB->multi-chip OUTPUT write-back accounting
            // (GLB read + multi-chip write access, serialized NoP link transfer, overlap).
            account_output_writeback_link();

            // Increase flush counter of weight and output data
            weight_flush_counter++;
            output_flush_counter++;

            // Waiting for weight and output data.
            wait_data();
        }
        // Case 2) request all data type
        else {
            // Input data, weight, and output data do not exist in Global buffer.
            exist_data[data_type_t::INPUT] = false, exist_data[data_type_t::WEIGHT] = false, exist_data[data_type_t::OUTPUT] = false;

#ifdef FUNCTIONAL
            // Write back Output data 
            m_scheduler->transfer_data(multi_chip->data, data, m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()], offsets[data_type_t::OUTPUT],
                                       component_type_t::CHIPS_Y, component_type_t::GLOBAL_BUFFER, 
                                       data_type_t::OUTPUT, get_stationary_type(), action_type_t::STORE);
#endif
            // GB1+GB2: descriptor-based GLB->multi-chip OUTPUT write-back accounting
            // (GLB read + multi-chip write access, serialized NoP link transfer, overlap).
            account_output_writeback_link();

            weight_flush_counter = 0; 
            output_flush_counter = 0;
            
            // Waiting for input data, weight, and output data.
            wait_data();
        }
    }
    // Case 2. Weight stationary
    else if(stationary_type == stationary_type_t::WEIGHT_STATIONARY) {
        // When not all input data and output data are transferred from DRAM to Global buffer.
        // Request input data and output data to DRAM.
        if(input_flush_counter < m_scheduler->offset_size_dram[data_type_t::INPUT].front() - 1 &&
           output_flush_counter < m_scheduler->offset_size_dram[data_type_t::OUTPUT].front() - 1) {
            // Input data and output data do not exist in Global buffer.
            exist_data[data_type_t::INPUT] = false, exist_data[data_type_t::OUTPUT] = false;
#ifdef FUNCTIONAL
            // Write back Output data 
            m_scheduler->transfer_data(multi_chip->data, data, m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()], offsets[data_type_t::OUTPUT],
                                       component_type_t::CHIPS_Y, component_type_t::GLOBAL_BUFFER, 
                                       data_type_t::OUTPUT, get_stationary_type(), action_type_t::STORE);
#endif
            // GB1+GB2: descriptor-based GLB->multi-chip OUTPUT write-back accounting
            // (GLB read + multi-chip write access, serialized NoP link transfer, overlap).
            account_output_writeback_link();
            
            input_flush_counter++;
            output_flush_counter++;

            // Waiting for input data and output data.
            // Global buffer is in idle state.
            wait_data();
        }
        // Case 2) Request all data type
        else {
            // Input data, weight, and output data do not exist in Global buffer.
            exist_data[data_type_t::INPUT] = false, exist_data[data_type_t::WEIGHT] = false, exist_data[data_type_t::OUTPUT] = false;

#ifdef FUNCTIONAL
            // Write back Output data 
            m_scheduler->transfer_data(multi_chip->data, data, m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()], offsets[data_type_t::OUTPUT],
                                       component_type_t::CHIPS_Y, component_type_t::GLOBAL_BUFFER, 
                                       data_type_t::OUTPUT, get_stationary_type(), action_type_t::STORE);
#endif
            // GB1+GB2: descriptor-based GLB->multi-chip OUTPUT write-back accounting
            // (GLB read + multi-chip write access, serialized NoP link transfer, overlap).
            account_output_writeback_link();
         
            input_flush_counter = 0;
            output_flush_counter = 0;

            // Waiting for input data, weight, and output data.
            // Global buffer is in idle state.
            wait_data();
        }
    }
    // Case 3. Output stationary
    else if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
        // Case 1) Reuse output data && Request input data and weight
        if(input_flush_counter < m_scheduler->offset_size_dram[data_type_t::INPUT].front() - 1 &&
           weight_flush_counter < m_scheduler->offset_size_dram[data_type_t::WEIGHT].front() - 1) {
            // Input data and weight do not exist in Global buffer.
            exist_data[data_type_t::INPUT] = false, exist_data[data_type_t::WEIGHT] = false;

            input_flush_counter++;
            weight_flush_counter++;

            // Waiting for input data and weight.
            // Global buffer is in idle state.
            wait_data();
        }
        // Case 2) Request all data types
        else {
            // Input data, weight, and output data do not exist in Global buffer.
            exist_data[data_type_t::INPUT] = false, exist_data[data_type_t::WEIGHT] = false, exist_data[data_type_t::OUTPUT] = false;
#ifdef FUNCTIONAL
            // Write back Output data 
            m_scheduler->transfer_data(multi_chip->data, data, m_scheduler->output_offset_multi_chip[index%m_scheduler->output_offset_multi_chip.size()], offsets[data_type_t::OUTPUT],
                                       component_type_t::CHIPS_Y, component_type_t::GLOBAL_BUFFER, 
                                       data_type_t::OUTPUT, get_stationary_type(), action_type_t::STORE);
#endif
            // GB1+GB2: descriptor-based GLB->multi-chip OUTPUT write-back accounting
            // (GLB read + multi-chip write access, serialized NoP link transfer, overlap).
            account_output_writeback_link();

            // Reset flush counter of input data and weight
            input_flush_counter = 0;
            weight_flush_counter = 0;

            // Waiting for input data, weight, and output data.
            // Global buffer is in idle state.
            wait_data();
        }
    }
    // The counter of Global buffer should be initialized as 0.
    input_index = 0, weight_index = 0, output_index = 0;
}

void global_buffer_t::update_static_energy(double elapsed_cycles) {
    if(elapsed_cycles < 0.0) {
        std::cerr << "Error: global buffer elapsed cycles must be non-negative" << std::endl;
        exit(1);
    }
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        if(u_static_energy[i] < 0.0) {
            std::cerr << "Error: global buffer static_energy must be a non-negative pJ/cycle value" << std::endl;
            exit(1);
        }
        static_energy[i] = static_energy_for_cycles(u_static_energy[i], elapsed_cycles);
    }
}

// Modeled busy duration of the global buffer for the current layer: the max of its
// access, transfer, and overlapped cost axes. Used so that leakage is charged over
// the true layer wall-clock rather than only the PE-array compute window.
double global_buffer_t::modeled_elapsed_cycles() const {
    double elapsed = std::max(write_back_cycle, overlapped_transfer_cycle);
    for(unsigned type = 0; type < data_type_t::NUM_DATA_TYPES; ++type) {
        // Reads (serving the array) and fill writes (from the multi-chip) may
        // serialize on the SRAM port, so the busy axis is their sum. Charging the
        // fill here pre-scaling also keeps busy >= final per-type access totals:
        // the fill later scales by a per-datatype factor <= the uniform factor.
        elapsed = std::max(elapsed, access_cycle[type] + fill_access_cycle[type]);
        elapsed = std::max(elapsed, transfer_cycle[type]);
        elapsed = std::max(elapsed, cycle_pe_array_global_buffer[type]);
    }
    return elapsed;
}

void global_buffer_t::reset() {
    // RE1: the final-cast counters reset with the layer.
    output_cast_bytes = 0;
    output_cast_energy = 0.0;
    std::fill_n(data, ((unsigned)size + sizeof(data_t) - 1)/sizeof(data_t), data_t{});
    
    idle = false;
    initial = true;

    write_back_cycle = 0.0;
    overlapped_transfer_cycle = 0.0;

    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);
    psum_writeback_events = 0;

    num_request.assign(data_type_t::NUM_DATA_TYPES, 0);
    num_data_transfer.assign(data_type_t::NUM_DATA_TYPES, 0);
    payload_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    metadata_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    storage_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);

    access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    fill_access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    fill_access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    cycle_pe_array_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    utilization.assign(data_type_t::NUM_DATA_TYPES, 0.0);
}

/* Separated Global buffer */
separate_buffer_t::separate_buffer_t(section_config_t m_section_config) :
    global_buffer_t(m_section_config),
    input_size(1.0),
    weight_size(1.0),
    output_size(1.0) {

    init(m_section_config);
}

separate_buffer_t::~separate_buffer_t() {
    delete [] data;
}

// Initialize the Global buffer
void separate_buffer_t::init(section_config_t m_section_config) {

    /* Initialize separate buffer's specifications */

    // Initialize buffer size of input, weight, and output in KB
    m_section_config.get_setting("input_size", &input_size);
    m_section_config.get_setting("weight_size", &weight_size);
    m_section_config.get_setting("output_size", &output_size);
    input_size *= 1024, weight_size *= 1024, output_size *= 1024;
    size = input_size + weight_size + output_size;

    // GB4: per-type utilization denominators are the partitions, not the summed size.
    capacity_per_type.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    capacity_per_type[data_type_t::INPUT] = input_size;
    capacity_per_type[data_type_t::WEIGHT] = weight_size;
    capacity_per_type[data_type_t::OUTPUT] = output_size;

    unsigned num_entry = ((unsigned)size + sizeof(data_t) - 1)/sizeof(data_t);
    data = new data_t[num_entry]();

    // Initialize the frequency and bandwidth of the separate buffer
    m_section_config.get_setting("frequency", &frequency);
    m_section_config.get_setting("bandwidth", &bandwidth);
    // B12: an explicit 'bitwidth' wins silently; a derived width warns when the
    // bandwidth/frequency ratio truncates fractionally.
    unsigned explicit_bitwidth = 0;
    if(m_section_config.get_setting("bitwidth", &explicit_bitwidth)) {
        bitwidth = explicit_bitwidth;
    } else {
        bitwidth = derived_link_bitwidth("global_buffer", bandwidth, frequency);
    }
    if(bitwidth == 0) {
        std::cerr << "Error: global_buffer requires a positive link bitwidth (set 'bitwidth' or a positive 'frequency')" << std::endl;
        exit(1);
    }

    // Initialize global buffer type (double buffer or single buffer)
    m_section_config.get_setting("double_buffer", &double_buffer);
    m_section_config.get_setting("fabric_separate", &fabric_separate);
    m_section_config.get_setting("streams_pipelined", &streams_pipelined);

    bypass.reserve(data_type_t::NUM_DATA_TYPES);
    bypass.assign(data_type_t::NUM_DATA_TYPES, 0);
    m_section_config.get_vector_setting("bypass", &bypass);


    // Initialize line size and mask bits of the global buffer
    line_size.reserve(data_type_t::NUM_DATA_TYPES);
    line_size.assign(data_type_t::NUM_DATA_TYPES, 8); // bits

    mask_bits.reserve(data_type_t::NUM_DATA_TYPES);
    mask_bits.assign(data_type_t::NUM_DATA_TYPES, 0);

    m_section_config.get_vector_setting("line_size", &line_size);
    validate_global_buffer_widths(bitwidth, line_size);

    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        while(line_size[i] > 8) {
            line_size[i] /= 2;
            mask_bits[i]++;
        }
    }
    m_section_config.get_vector_setting("line_size", &line_size);

    // Initialize the stationary type of the global buffer
    std::string stationary_str;
    if(m_section_config.get_setting("stationary_type", &stationary_str)) {
        stationary_type = (stationary_type_t)get_type(stationary_type_str, stationary_str);
    }

    if(stationary_type == stationary_type_t::INPUT_STATIONARY) {
        parameter_order = "BCPQKRS";
    } else if(stationary_type == stationary_type_t::WEIGHT_STATIONARY) {
        parameter_order = "KCRSBPQ";
    } else if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
        parameter_order = "BKPQKRS";
    }
    m_section_config.get_setting("parameter_order", &parameter_order);

    // Initialize the memory type.
    memory_type = memory_type_t::SEPARATE;

    skip_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);
    psum_writeback_events = 0;
    
    // Initialize the tile size of Global buffer.
    tile_size.reserve(data_type_t::NUM_DATA_TYPES);
    tile_size.assign(data_type_t::NUM_DATA_TYPES, 1);

    offsets.reserve(data_type_t::NUM_DATA_TYPES);
    offsets.assign(data_type_t::NUM_DATA_TYPES, 0);

    /* Initialize signals of the global buffer */

    // Initialize data exist signal
    exist_data.reserve(data_type_t::NUM_DATA_TYPES);
    exist_data.assign(data_type_t::NUM_DATA_TYPES, false);

    // Initialize data request signal
    request_to_multi_chip.reserve(data_type_t::NUM_DATA_TYPES);
    request_to_multi_chip.assign(data_type_t::NUM_DATA_TYPES, false);

    /* Initialize unit cost of the global buffer */

    // Initialize the unit cycle and energy of Global buffer
    u_read_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("read_cycle", &u_read_cycle);

    u_read_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("read_energy", &u_read_energy);

    u_write_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("write_cycle", &u_write_cycle);

    u_write_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("write_energy", &u_write_energy);

    // GLB bank parallelism (ideal line interleaving): with B banks, B accesses proceed
    // in parallel, so the effective per-access cycle is u/B across every access path.
    // Per-access energy is unchanged. Bank conflicts, arbitration, and read/write port
    // contention are NOT modeled -- this knob covers the conflict-free upper bound.
    unsigned num_banks = 1;
    m_section_config.get_setting("num_banks", &num_banks);
    if(num_banks == 0) {
        std::cerr << "Error: global_buffer num_banks must be non-zero" << std::endl;
        exit(1);
    }
    if(num_banks > 1) {
        for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
            u_read_cycle[i] /= num_banks;
            u_write_cycle[i] /= num_banks;
        }
    }
    
    u_static_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("static_energy", &u_static_energy);
    // RE1: the final output cast/pack unit cost, charged per output byte read out.

    m_section_config.get_setting("transfer_cycle", &u_transfer_cycle);
    m_section_config.get_setting("transfer_energy", &u_transfer_energy);

    /* Initialize stats of the global buffer */
    payload_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    metadata_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    storage_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the number of request to Global buffer
    num_request.reserve(data_type_t::NUM_DATA_TYPES);
    num_request.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the number of data transfer to PE array
    num_data_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the total access cycle
    access_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    fill_access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    fill_access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize the total access energy
    access_energy.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    cycle_pe_array_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    cycle_pe_array_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    static_energy.reserve(data_type_t::NUM_DATA_TYPES);
    static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize the total transfer cycle between PE array and Global buffer
    transfer_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize the total transfer energy between PE array and Global buffer
    transfer_energy.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    utilization.reserve(data_type_t::NUM_DATA_TYPES);
    utilization.assign(data_type_t::NUM_DATA_TYPES, 0.0);

}

void separate_buffer_t::update_offset() {
    offsets[data_type_t::INPUT] = 0;
    offsets[data_type_t::WEIGHT] = input_size/sizeof(data_t);
    offsets[data_type_t::OUTPUT] = input_size/sizeof(data_t) + weight_size/sizeof(data_t);

}

void separate_buffer_t::check_tile_size() {
    // GB3: a double-buffered SRAM must hold two live tile copies (fill one half while
    // the other is consumed), so the capacity check doubles the required bytes.
    const size_t copies = double_buffer ? 2 : 1;
    const size_t input_bytes = copies*runtime_datatypes().storage_bytes(data_type_t::INPUT, tile_size[data_type_t::INPUT]);
    const size_t weight_bytes = copies*runtime_datatypes().storage_bytes(data_type_t::WEIGHT, tile_size[data_type_t::WEIGHT]);
    const size_t output_bytes = copies*runtime_datatypes().storage_bytes(data_type_t::OUTPUT, tile_size[data_type_t::OUTPUT]);
    if(input_bytes > input_size && !bypass[data_type_t::INPUT]) {
        std::cerr << "Input data size : " << input_bytes
                  << " is bigger than input buffer size : " << input_size << std::endl;
        exit(1);
    }
    if(weight_bytes > weight_size && !bypass[data_type_t::WEIGHT]) {
        std::cerr << "Weight data size : " << weight_bytes
                  << " is bigger than weight buffer size : " << weight_size << std::endl;
        exit(1);
    }
    if(output_bytes > output_size && !bypass[data_type_t::OUTPUT]) {
        std::cerr << "Output data size : " << output_bytes
                  << " is bigger than output buffer size : " << output_size << std::endl;
        exit(1);
    }
}
// Print out the stat of separate Global buffer.
void separate_buffer_t::print_specification() {
    std::cout << "============== Global buffer ===============" << std::endl;

    // Print out the buffer size.
    std::cout << "Global buffer type :" << std::setw(24)
                                        << "Separate" << std::endl;

    std::cout << "Input size         :" << std::setw(21) << std::setprecision(1)
                                        << input_size/1024 << " KB" << std::endl;
    std::cout << "Weight size        :" << std::setw(21) << std::setprecision(1)
                                        << weight_size/1024 << " KB" << std::endl;
    std::cout << "Output size        :" << std::setw(21) << std::setprecision(1)
                                        << output_size/1024 << " KB" << std::endl;

    // Print line size.
    std::cout << "Line size" << std::endl;
    std::cout << " * Input data      :" << std::setw(19) << std::setprecision(0)
                                        << line_size[data_type_t::INPUT] << " bits" << std::endl;
    std::cout << " * Weight          :" << std::setw(19) << std::setprecision(0)
                                        << line_size[data_type_t::WEIGHT] << " bits" << std::endl;
    std::cout << " * Output data     :" << std::setw(19) << std::setprecision(0)
                                        << line_size[data_type_t::OUTPUT] << " bits" << std::endl;    
    // Print out the stationary type.
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
    
    std::cout << "Access energy (read/write)" << std::endl;
    std::cout << " * Input buffer    :" << std::setw(16) << std::setprecision(2)
                                        << u_read_energy[data_type_t::INPUT] << "/" << u_write_energy[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    std::cout << " * Weight buffer   :" << std::setw(16) << std::setprecision(2)
                                        << u_read_energy[data_type_t::WEIGHT] << "/" << u_write_energy[data_type_t::WEIGHT] << " " << energy_units().label() << std::endl;
    std::cout << " * Output buffer   :" << std::setw(16) << std::setprecision(2)
                                        << u_read_energy[data_type_t::OUTPUT] << "/" << u_write_energy[data_type_t::OUTPUT] << " " << energy_units().label() << std::endl;

    std::cout << "Access cycle (read/write)" << std::endl;
    std::cout << " * Input buffer    :" << std::setw(13) << std::setprecision(1)
                                        << u_read_cycle[data_type_t::INPUT] << "/" << u_write_cycle[data_type_t::INPUT] << " cycles" << std::endl;
    std::cout << " * Weight buffer   :" << std::setw(13) << std::setprecision(1)
                                        << u_read_cycle[data_type_t::WEIGHT] << "/" << u_write_cycle[data_type_t::WEIGHT] << " cycles" << std::endl;
    std::cout << " * Output buffer   :" << std::setw(13) << std::setprecision(1)
                                        << u_read_cycle[data_type_t::OUTPUT] << "/" << u_write_cycle[data_type_t::OUTPUT] << " cycles" << std::endl;
    std::cout << std::endl;
}

shared_buffer_t::shared_buffer_t(section_config_t m_section_config) :
    global_buffer_t(m_section_config) {

    init(m_section_config);
}

shared_buffer_t::~shared_buffer_t() {
    delete [] data;
}

// Initialize the Global buffer.
void shared_buffer_t::init(section_config_t m_section_config) {

    /* Initialize shared buffer's specifications */

    // Initialize size of shared buffer in KB
    m_section_config.get_setting("memory_size", &size);
    // KB -> Byte
    size *= 1024;

    // A shared buffer holds every type in the same SRAM, so each type's utilization
    // denominator is the full capacity.
    capacity_per_type.assign(data_type_t::NUM_DATA_TYPES, size);

    unsigned num_entry = ((unsigned)size + sizeof(data_t) - 1)/sizeof(data_t);
    data = new data_t[num_entry]();

    // Initialize frequency and bandwidth of the shared buffer
    m_section_config.get_setting("frequency", &frequency);
    m_section_config.get_setting("bandwidth", &bandwidth);
    // B12: an explicit 'bitwidth' wins silently; a derived width warns when the
    // bandwidth/frequency ratio truncates fractionally.
    unsigned explicit_bitwidth = 0;
    if(m_section_config.get_setting("bitwidth", &explicit_bitwidth)) {
        bitwidth = explicit_bitwidth;
    } else {
        bitwidth = derived_link_bitwidth("global_buffer", bandwidth, frequency);
    }
    if(bitwidth == 0) {
        std::cerr << "Error: global_buffer requires a positive link bitwidth (set 'bitwidth' or a positive 'frequency')" << std::endl;
        exit(1);
    }

    m_section_config.get_setting("double_buffer", &double_buffer);
    m_section_config.get_setting("fabric_separate", &fabric_separate);
    m_section_config.get_setting("streams_pipelined", &streams_pipelined);

    bypass.reserve(data_type_t::NUM_DATA_TYPES);
    bypass.assign(data_type_t::NUM_DATA_TYPES, 0);
    m_section_config.get_vector_setting("bypass", &bypass);
    
    // Initialize line size and mask bits of the global buffer
    line_size.reserve(data_type_t::NUM_DATA_TYPES);
    line_size.assign(data_type_t::NUM_DATA_TYPES, 8); // bits

    mask_bits.reserve(data_type_t::NUM_DATA_TYPES);
    mask_bits.assign(data_type_t::NUM_DATA_TYPES, 0);

    m_section_config.get_vector_setting("line_size", &line_size);
    validate_global_buffer_widths(bitwidth, line_size);

    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; i++) {
        while(line_size[i] > 8) {
            line_size[i] /= 2;
            mask_bits[i]++;
        }
    }
    m_section_config.get_vector_setting("line_size", &line_size);

    // Initialize the stationary type.
    std::string stationary_str;
    if(m_section_config.get_setting("stationary_type", &stationary_str)) {
        stationary_type = (stationary_type_t)get_type(stationary_type_str, stationary_str);
    }
    if(stationary_type == stationary_type_t::INPUT_STATIONARY) {
        parameter_order = "BCPQKRS";
    } else if(stationary_type == stationary_type_t::WEIGHT_STATIONARY) {
        parameter_order = "KCRSBPQ";
    } else if(stationary_type == stationary_type_t::OUTPUT_STATIONARY) {
        parameter_order = "BKPQKRS";
    }
    m_section_config.get_setting("parameter_order", &parameter_order);

    // Initialize the memory type.
    memory_type = memory_type_t::SHARED;


    skip_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    skip_transfer.assign(data_type_t::NUM_DATA_TYPES, false);
    psum_writeback_events = 0;

    // Initialize the tile size of the global buffer
    tile_size.reserve(data_type_t::NUM_DATA_TYPES);
    tile_size.assign(data_type_t::NUM_DATA_TYPES, 1);

    offsets.reserve(data_type_t::NUM_DATA_TYPES);
    offsets.assign(data_type_t::NUM_DATA_TYPES, 0);

    /* Initialize signals of shared buffer */

    // Initialize the signal to indicate whether the data exists or not.
    exist_data.reserve(data_type_t::NUM_DATA_TYPES);
    exist_data.assign(data_type_t::NUM_DATA_TYPES, false);

    // Initialize the request signal
    request_to_multi_chip.reserve(data_type_t::NUM_DATA_TYPES);
    request_to_multi_chip.assign(data_type_t::NUM_DATA_TYPES, false);

    /* Initialize the unit cost of the shared buffer */

    // Initialize the unit cycle and energy of Global buffer
    u_read_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("read_cycle", &u_read_cycle);

    u_read_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_read_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("read_energy", &u_read_energy);

    u_write_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("write_cycle", &u_write_cycle);

    u_write_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_write_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("write_energy", &u_write_energy);

    // GLB bank parallelism (ideal line interleaving): with B banks, B accesses proceed
    // in parallel, so the effective per-access cycle is u/B across every access path.
    // Per-access energy is unchanged. Bank conflicts, arbitration, and read/write port
    // contention are NOT modeled -- this knob covers the conflict-free upper bound.
    unsigned num_banks = 1;
    m_section_config.get_setting("num_banks", &num_banks);
    if(num_banks == 0) {
        std::cerr << "Error: global_buffer num_banks must be non-zero" << std::endl;
        exit(1);
    }
    if(num_banks > 1) {
        for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
            u_read_cycle[i] /= num_banks;
            u_write_cycle[i] /= num_banks;
        }
    }

    u_static_energy.reserve(data_type_t::NUM_DATA_TYPES);
    u_static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    m_section_config.get_vector_setting("static_energy", &u_static_energy);
    // RE1: the final output cast/pack unit cost, charged per output byte read out.

    payload_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    metadata_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    storage_link_transactions.assign(data_type_t::NUM_DATA_TYPES, 0);
    m_section_config.get_setting("transfer_cycle", &u_transfer_cycle);
    m_section_config.get_setting("transfer_energy", &u_transfer_energy);

    /* Initialize stats of the shared buffer */
 
    // Initialize the number of request 
    num_request.reserve(data_type_t::NUM_DATA_TYPES);
    num_request.assign(data_type_t::NUM_DATA_TYPES, 0);

    num_data_transfer.reserve(data_type_t::NUM_DATA_TYPES);
    num_data_transfer.assign(data_type_t::NUM_DATA_TYPES, 0);

    // Initialize the access cycle 
    access_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    fill_access_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    fill_access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize the access energy.
    access_energy.reserve(data_type_t::NUM_DATA_TYPES);
    access_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);   

    cycle_pe_array_global_buffer.reserve(data_type_t::NUM_DATA_TYPES);
    cycle_pe_array_global_buffer.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    static_energy.reserve(data_type_t::NUM_DATA_TYPES);
    static_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    // Initialize the data transfer cycle.
    transfer_cycle.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_cycle.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    
    // Initialize the data transfer energy.
    transfer_energy.reserve(data_type_t::NUM_DATA_TYPES);
    transfer_energy.assign(data_type_t::NUM_DATA_TYPES, 0.0);

    utilization.reserve(data_type_t::NUM_DATA_TYPES);
    utilization.assign(data_type_t::NUM_DATA_TYPES, 0.0);
    
}

void shared_buffer_t::update_offset() {
    offsets[data_type_t::INPUT] = 0;
    offsets[data_type_t::WEIGHT] = tile_size[data_type_t::INPUT];
    offsets[data_type_t::OUTPUT] = tile_size[data_type_t::INPUT] + tile_size[data_type_t::WEIGHT];
}

void shared_buffer_t::check_tile_size() {
    size_t data_size = 0;
    for(unsigned i = 0; i < data_type_t::NUM_DATA_TYPES; ++i) {
        const data_type_t type = static_cast<data_type_t>(i);
        if(!bypass[type]) data_size += runtime_datatypes().storage_bytes(type, tile_size[type]);
    }
    // GB3: double buffering needs two live copies of the working set.
    if(double_buffer) data_size *= 2;

    if(data_size > size) {
        std::cout << "The data size is bigger than Global buffer size\n"
                  << "Data : " << data_size
                  << " Buffer : " << size << std::endl;
        exit(1);
    }
}

// Print out the stat of shared Global buffer.
void shared_buffer_t::print_specification() {
    std::cout << "============== Global buffer ===============" << std::endl;

    std::cout << "Global buffer type :" << std::setw(24) 
                                        << "Shared" << std::endl;
    std::cout << "Buffer size        :" << std::setw(21)
                                        << size/1024 << " KB" << std::endl;
    std::cout << "Line size          :" << std::setw(19) 
                                        << line_size[data_type_t::INPUT] << " bits" << std::endl;

    // Print out the stationary type.
    if(stationary_type == stationary_type_t::INPUT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24)
                                            << "Input stationary" << std::endl;
    }
    else if(stationary_type == stationary_type_t::WEIGHT_STATIONARY) {
        std::cout << "Stationary type    :" << std::setw(24)                     
                                            << "Weight stationary" << std::endl;
    } 
    else {
        std::cout << "Stationary type    :" << std::setw(24)                     
                                            << "Output stationary" << std::endl;
    }
    std::cout << "Bandwidth          :" << std::setw(19) << std::setprecision(0)
                                        << bandwidth << " GB/s" << std::endl;
    
    std::cout << "Access energy (read/write)" << std::endl;
    std::cout << " * Global buffer   :" << std::setw(16) << std::setprecision(2) 
                                        << u_read_energy[data_type_t::INPUT] << "/" << u_write_energy[data_type_t::INPUT] << " " << energy_units().label() << std::endl;
    
    std::cout << "Access cycle (read/write) " << std::endl;
    std::cout << " * Global buffer    :" << std::setw(13) << std::setprecision(1)
                                        << u_read_cycle[data_type_t::INPUT] << "/" << u_write_cycle[data_type_t::INPUT] << " cycles" << std::endl;
    std::cout << std::endl;
}

