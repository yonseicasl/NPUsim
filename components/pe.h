#ifndef __PE_H__
#define __PE_H__

#include <cstddef>
#include <string>

#include "def.h"
#include "utils.h"
#include "data.h"

#include "scheduler.h"
#include "pe_array.h"
#include "pe_lane.h"

class pe_array_t;
class pe_t {

public:

    // Phase-5: the component's clock in MHz. Power needs a single authoritative clock
    // across the modeled components (the timeline is one shared cycle axis), so stats_t
    // compares these and reports power as unsupported when they disagree.
    double clock_mhz() const { return static_cast<double>(frequency); }
    pe_t(section_config_t m_section_config);
    virtual ~pe_t();

    // Initialize the PE.
    void init(section_config_t m_section_config);

    // Connect PE to PE array.
    void connect(pe_array_t *m_pe_array);

    // Update tile size of PE.
    void update_tile_size(scheduler_t *m_scheduler);
    // Update the data offset
    void update_offset();
    // Check tile size of PE
    void check_tile_size();

    /* Get PE specifications */

    // Get stationary type of MAC register
    stationary_type_t get_mac_stationary_type();
    // Get stationary type of local buffer.
    stationary_type_t get_local_buffer_stationary_type();
    // Return parameter order.
    std::string get_parameter_order();
    // Return memory type.
    memory_type_t get_memory_type();

    /* Check PE status */

    // Return true when the PE has no data.
    bool is_idle();
    // Check whether all data exist in the Local buffer.
    bool is_exist_data();
    // Check whether at least one request exist in the local buffer.
    bool is_exist_request();

    /* PE actions */

    // Wait for the data.
    void wait_data();
    // Exist data in PE.
    void fill_data();
    // Request data to PE array.
    void request_data();
    // Execute the MAC operation.
    void data_transfer_to_mac(scheduler_t *m_scheduler);
    // Flush the PE data.
    void flush_data(scheduler_t *m_scheduler);
    // Execute MAC operation
    virtual void computation(scheduler_t *m_scheduler) = 0;

    // MAC operation
    void mac_operation(scheduler_t *m_scheduler);
    // NOTE: the old no-arg pe_t "activation" member (hardcoded per-PE ReLU) is REMOVED,
    // not merely deprecated (plan/plan_sfu.md). Activation cost is modeled by the SFU
    // (components/sfu.{h,cc}) at the final output commit, and functional activation is
    // owned solely by Nebula's forward(); reintroducing the method would double-apply
    // and double-charge. Removing the declaration makes any reintroduced call a compile
    // error instead of a runtime abort, and unittest/run_validation.sh greps against it.

    // Print out the configuration of PE.
    void print_specification();

    // Reset component and stats.
    void reset();

    // Timing-only helpers. static_energy is pJ/cycle multiplied by the
    // elapsed layer time, not an event counter.
    double modeled_elapsed_cycles() const;
    void update_static_energy(double elapsed_cycles);
    size_t scalar_mac_capacity() const;

    /* PE specifications */

    // Data at MAC unit.
    data_t *input_data_mac;                                 // Input data in MAC unit
    data_t *weight_mac;                                     // Weight in MAC unit
    data_t *output_data_mac;                                // Output data in MAC unit

    // Data at local buffer.
    data_t *input_data_lb;                                  // input data in local buffer
    data_t *weight_lb;                                      // weight in local buffer
    data_t *output_data_lb;                                 // output data in local buffer

    unsigned input_size;                                    // Input local buffer size
    unsigned weight_size;                                   // Weight local buffer size
    unsigned output_size;                                   // Output local buffer size

    std::vector<bool> bypass;                               // Check if bypass is applied at the local buffer
    bool edge_accumulation;                                 // V3: outputs accumulate at the array edge; OUTPUT exempt from LB capacity
    bool double_buffer;                                     // LB7/P1-A: timing-only overlap of the LB<->MAC transfer with compute (single buffer serializes it); one resident tile copy in BOTH modes -- see the capacity contract in pe_t::init()
    unsigned index;                                         // Index of PE in PE array.
    // L9/PE1-PE2: the ACTIVE lane state this layer's mapping resolves to (how many
    // accumulators are live and how many lanes feed the last one). The lane->accumulator
    // reduction is charged against this, not against the structural mac_width -- see
    // lane_reduction_fill_cycles()/lane_reduction_energy() in utils/pe_lane.h.
    mac_lane_state_t lane_state;

    std::vector<bool> skip_transfer;                        // Check whether skip data transfer from local buffer to MAC unit

    // Data tile size.
    std::vector<unsigned> tile_size_mac;                    // Tile size of MAC unit
    std::vector<unsigned> tile_size_lb;                     // Tile size of local buffer.

    std::vector<unsigned> offsets_mac;                      //
    std::vector<unsigned> offsets_lb;                       //

    unsigned duplicated_input_mac;                          //
    unsigned duplicated_input_lb;                           //

    /* PE signals */
    std::vector<bool> exist_data_mac;                       // Data exist in MAC
    std::vector<bool> request_to_lb;                        // Request data from MAC unit to Local buffer.

    std::vector<bool> exist_data_lb;                        // Data exist in local buffer
    std::vector<bool> request_to_pe_array;                  // Request data from Local buffer to PE array.


    /* Unit costs of PE */

    double u_computation_cycle;                             // The unit MAC cycle
    double u_computation_energy;                            // The unit MAC energy
    // E3: MAC energy used to be ONE scalar for every operand precision, so switching
    // input/weight formats changed the memory traffic and left the compute energy identical --
    // an INT4 multiply cost exactly what an FP16 multiply cost, and nothing flagged it. A config
    // may now declare a per-precision cost, keyed by the operand format pair:
    //     mac_energy_<input format>_<weight format>    e.g. mac_energy_int8_int8
    // When the running pair has no declared entry the scalar is still used (so existing configs
    // are unaffected), but the result is reported as UNCALIBRATED for that precision rather than
    // presented as a precision-aware number. `compute_energy_basis` names whichever key supplied
    // the value.
    // RE1: the final output cast is charged at the OFF-CHIP STORE, so its unit cost lives in
    // multi_chip_t, not here. Two earlier placements -- the PE's write_output (262,144 casts, one
    // per MAC issue) and the GLB readout (16,384, one per GLB pass) -- were both wrong, and each
    // left a `u_output_cast_*` pair behind in the component it had been tried in. Those pairs were
    // still parsed from the config and consumed by nothing, which is the exact defect the
    // dead-knob scan (validation/knobs KN9) exists to catch. They are gone; only the multi-chip
    // pair is real.
    // E4/RE1: the accumulator/output event counts, kept apart because they are DIFFERENT EVENTS
    // on different precisions and at different boundaries.
    //
    // RE1: these used to be produced together, unconditionally, by account_format_events(OUTPUT)
    // on every output request -- before that function could even tell whether the request was a
    // fresh zero-initialized accumulator or a genuine partial-sum reload. The consequence was that
    // the "final cast" count followed the MAC count (262,144 on GEMM 64x64x64) instead of the
    // final output element count (4,096), and a zero-init paid reload energy it never spent. The
    // four events are now generated at their own boundaries:
    //   CREATE  -- request_to_lb[OUTPUT] with no prior value: clear_output_accumulators(), free
    //   RELOAD  -- request_to_lb[OUTPUT] with a prior value: accumulator-precision read
    //   SPILL   -- MAC -> LB write-back of a partial sum: accumulator-precision write
    //   CAST    -- the output leaving the PE for good: output-precision pack, once per element
    size_t accumulator_reload_bytes;
    size_t accumulator_spill_bytes;
    size_t accumulator_create_events;
    // E20-4: passes where a prior partial sum existed but the output tile was already resident in
    // the MAC, so no LB->MAC read-back occurred. These used to be charged as reloads.
    size_t accumulator_retained_events;
    // RE1: the accumulator's own energy, kept separate from the format-IP energy so it can be
    // attributed to the component that OWNS the accumulator. With edge_accumulation the
    // accumulator lives at the array edge, not in the PE, so this energy is handed to the
    // PE array instead of being charged here.
    double accumulator_energy;
    // E4: lane->accumulator reduction energy, split out of computation_energy so the adder-tree
    // work is visible on its own axis instead of being folded into the MAC total.
    double reduction_energy;
    std::string compute_energy_basis;
    bool compute_energy_precision_calibrated;
    // P4-2/PE2: unit energy of one lane->accumulator adder-tree addition (see
    // lane_reduction_energy()). Defaults to 0 (no-op) until a config calibrates it,
    // mirroring adder_tree_t's u_adder_energy convention.
    double u_mac_reduction_energy;

    double u_transfer_cycle;                                // The unit data transfer cycle between MAC and local buffer
    double u_transfer_energy;                               // The unit data transfer energy between MAC and local buffer

    std::vector<double> u_read_cycle_mac;                   // The unit MAC read cycle
    std::vector<double> u_read_energy_mac;                  // The unit MAC read energy

    std::vector<double> u_write_cycle_mac;                  // The unit MAC write cycle
    std::vector<double> u_write_energy_mac;                 // The unit MAC write energy

    std::vector<double> u_read_cycle_lb;                    // The unit local buffer read cycle
    std::vector<double> u_read_energy_lb;                   // The unit local buffer read energy

    std::vector<double> u_write_cycle_lb;                   // The unit local buffer write cycle
    std::vector<double> u_write_energy_lb;                  // The unit local buffer write energy

    std::vector<double> u_static_energy;                    // PE (MAC-side) leakage, pJ/cycle
    std::vector<double> u_lb_static_energy;                 // LB4: local-buffer leakage, pJ/cycle (separately sweepable)

    // Optional format-IP cost per packed payload/metadata transaction.
    std::vector<double> u_format_payload_cycle;
    std::vector<double> u_format_payload_energy;
    std::vector<double> u_format_metadata_cycle;
    std::vector<double> u_format_metadata_energy;
    std::vector<double> u_accumulator_spill_cycle;
    std::vector<double> u_accumulator_spill_energy;


    /* PE stats */

    unsigned num_computation;                               // The number of computations

    double    computation_cycle;                            // Total computation cycle
    double    computation_energy;                           // Total computation energy.

    std::vector<unsigned> num_request_to_lb;                // The number of data request from MAC unit to local buffer
    std::vector<unsigned> num_data_transfer_to_mac;         // The number of data transfer from local buffer to MAC unit

    std::vector<double> access_cycle_mac;                   // Total access cycles to MAC unit
    std::vector<double> access_energy_mac;                  // Total access energies to MAC unit

    std::vector<double> access_cycle_lb;                    // Total access cycles to Local buffer.

    // Dynamic conversion/pack/decode/scale/spill costs by tensor type.
    std::vector<double> format_cycle;
    std::vector<double> format_energy;
    std::vector<double> access_energy_lb;                   // Total access energies to Local buffer.

    std::vector<double> transfer_cycle;                     // Total transfer cycle between MAC unit and local buffer
    std::vector<double> transfer_energy;                    // Total transfer energy between MAC unit and local buffer

    std::vector<size_t> payload_link_transactions;          // Logical payload transactions on MAC-local-buffer link.
    std::vector<size_t> metadata_link_transactions;         // Logical datatype metadata transactions.
    std::vector<size_t> storage_link_transactions;          // Physical serialized transactions.

    std::vector<double> cycle_mac_lb;                       // Total cycle between MAC unit and local buffer.

    std::vector<double> static_energy;                      // Static energy.

    double utilization_mac;                                 // Utilization of MAC units.
    std::vector<double> utilization_local_buffer;           // Local buffer utilization.

    double write_back_cycle_mac;                            // Access cycle to MAC units when writing back output data.
    double write_back_cycle_lb;                             // Access cycle to local buffers when writing back output data.
    double overlapped_transfer_cycle;                       // Total transfer cycles between MAC unit and local buffer

    std::vector<unsigned> line_size_mac;                    // Data block size of MAC unit.
    std::vector<unsigned> line_size_lb;                     // Data block size of local buffer.

    std::vector<unsigned> mask_bits_mac;                    // Mask bits
    std::vector<unsigned> mask_bits_lb;                     // Mask bits

    // E20-2: the Format-IP transaction counts, per stream. The energy axis alone cannot say
    // whether a zero means "no format work" or "format work whose unit cost was never
    // declared".
    std::vector<size_t> format_payload_events;
    std::vector<size_t> format_metadata_events;

protected:

    pe_array_t *pe_array;

    stationary_type_t stationary_type_mac;                  // Stationary type of MAC register
    stationary_type_t stationary_type_local_buffer;         // stationary type of local buffer
    std::string parameter_order;                            //
    memory_type_t memory_type;                              // Memory type.
    mac_type_t mac_type;                                    // Mac type.

    // Number of MACS.
    unsigned num_macs;                                      // Number of independent accumulator units.
    unsigned mac_width;                                     // Scalar FMA lanes per accumulator unit.
    unsigned num_active_macs;                               // Active scalar FMA lanes from the mapping.
    unsigned active_mac_width;                              // Lanes occupied in the final active MAC unit.
    unsigned active_mac_units;                              // Active accumulator units for the current tile.
    size_t mac_register_capacity;                            // Scalar register entries allocated per operand.

    float frequency;                                        // The frequency (GHz)
    float bandwidth;                                        // The bandwidth between MAC and local buffer (GB/sec)
    unsigned bitwidth;                                      // The bitwidth between MAC and local buffer

    // Counters for the number of data transfer.
    unsigned input_index;                                   // Input index  in local buffer.
    unsigned weight_index;                                  // Weight index in local buffer.
    unsigned output_index;                                  // Output index in local buffer.

    // Flush counter
    unsigned input_flush_counter;                           // Input data flush counter.
    unsigned weight_flush_counter;                          // Weight flush counter.
    unsigned output_flush_counter;                          // Output data flush counter.

    // Whether the PE is in idle state or not.
    bool idle;

    void clear_output_accumulators();
    unsigned count_nonzero_mac_operations(scheduler_t *m_scheduler) const;

    void account_format_events(data_type_t type, size_t elements);

    // RE1: the four output-accumulator events, separated so each is generated at its own
    // boundary instead of all of them firing on every output request.
    void account_accumulator_create(size_t elements);
    void account_accumulator_reload(size_t elements);
    void account_accumulator_spill(size_t elements);
    void charge_accumulator_bytes(size_t bytes);
    void account_descriptor_dense_mac_transfer(data_type_t type, size_t elements, bool to_mac);
};

class undefined_stationary_t : public pe_t {
public:
    undefined_stationary_t(section_config_t m_section_config);
    ~undefined_stationary_t();

    void computation(scheduler_t *m_sheduler);
};

class input_stationary_t : public pe_t {

public:
    input_stationary_t(section_config_t m_section_config);
    ~input_stationary_t();

    // Execute MAC operation
    void computation(scheduler_t *m_scheduler);
};

class weight_stationary_t : public pe_t {

public:
    weight_stationary_t(section_config_t m_section_config);
    ~weight_stationary_t();

    // Execute the MAC operation
    void computation(scheduler_t *m_scheduler);
};

class output_stationary_t : public pe_t {

public:
    output_stationary_t(section_config_t m_section_config);
    ~output_stationary_t();

    //Execute the MAC operation
    void computation(scheduler_t *m_scheduler);
};

#endif
