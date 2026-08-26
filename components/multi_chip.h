#ifndef __MULTI_CHIP_H__
#define __MULTI_CHIP_H__

#include <vector>
#include "global_buffer.h"
#include "dram.h"
#include "scheduler.h"

class global_buffer_t;
class dram_t;

class multi_chip_t {

public:

    // RE1: the FINAL output cast/pack to output precision, charged at the off-chip output store.
    //
    // Every earlier boundary fires more than once per output element: a PE writes its tile back
    // once per reduction pass, and the GLB reads it out once per GLB-level reduction pass. Charging
    // the cast at either made its count follow the MAC count (262,144 on GEMM 64x64x64) or a
    // multiple of the output count (16,384) instead of the 4,096 final elements. The off-chip store
    // is the boundary DR6/T1 establish happens exactly once per output element, and it is where the
    // value is committed in OUTPUT precision -- so it is both the correct count and a defensible
    // owner.
    std::vector<double> u_output_cast_cycle;
    std::vector<double> u_output_cast_energy;
    // RE6: which links this NoP's noc_energy prices, for the report.
    std::string describe_nop_link_contract() const;

    // The multi-chip temporal buffer's occupancy per datatype (peak over the layer). Every other
    // level reports this -- the PE local buffer, the PE-array temporal buffer and the GLB all do --
    // and the multi-chip level was the only one that did not, which left `memory_size` consumed
    // ONLY by the capacity check: shrinking it 16x changed no reported number at all. Found by the
    // dead-knob sweep (validation/knobs).
    std::vector<double> buffer_utilization;

    size_t output_cast_bytes;
    double output_cast_energy;
    // Phase-2 (plan_sfu.md): final_output_tile events. account_output_writeback_to_dram()
    // is the single boundary where the output commits in OUTPUT precision exactly once
    // per element (RE1/DR6/T1), which is precisely the reduction-complete,
    // before-writeback point the SFU activation contract needs -- so each call here IS
    // one final_output_tile event. Spills/reloads (GLB <-> array psum traffic) never
    // pass through this boundary, so the count is structurally spill-independent.
    size_t final_output_tile_events;
    size_t final_output_tile_elements;

    // Phase-5: the component's clock in MHz. Power needs a single authoritative clock
    // across the modeled components (the timeline is one shared cycle axis), so stats_t
    // compares these and reports power as unsupported when they disagree.
    double clock_mhz() const { return static_cast<double>(frequency); }
    
    multi_chip_t(section_config_t m_section_config);
    ~multi_chip_t();

    // Initialize the chip-level processors
    void init(section_config_t m_section_config);

    // Connect chip-level processors and the global buffer
    void connect(std::vector<global_buffer_t*> m_global_buffers);
    // Connect to chip-level processors and the off-chip memory.
    void connect(dram_t *m_dram);

    // Update tile size of chip-level processors.
    void update_tile_size(scheduler_t *m_scheduler);
    // Update offsets of each data of neural network.
    void update_offset();

    /* Get the chip-level processor's specifications */

    // Get stationary type.
    stationary_type_t get_stationary_type();
    // Get parameter order.
    std::string get_parameter_order();
    // Get the number of chips at chip-level processors
    unsigned get_number_of_chips();
    // Get the number of active chips at chip-level processors
    unsigned get_number_of_active_chips();

    // Get the NoP link bit-width between the global buffer and chip-level processors.
    unsigned get_bitwidth();

    // NoP hops from the package ingress to a chip: Manhattan distance on the mesh,
    // 1 on the serialized bus.
    unsigned nop_hops_for_chip(unsigned chip_index) const;

    /* Check chip-level processor's status */

    // Check whether the chip-level processor is idle or not
    bool is_idle();
    // Check whether the chip-level processor receives request from the global buffer
    bool is_exist_request_at_global_buffer();
    // Check whether the request exist in the PE array.
    bool is_exist_request_at_buffer();
    // Check whether the data exist is the Multi Chip
    bool is_exist_data();

    /* Chip-level processor's actions */

    // Wait for the data.
    bool wait_data();
    // Exist data in PE array.
    void fill_data();
    // Request data to Global buffer.
    void request_data();
    // DR6: charge the layer's final output tile write-back to DRAM (the in-loop
    // write-back only fires on tile CHANGES, so the last resident tile of the
    // live pass is otherwise never stored).
    void flush_output_writeback();
    // Shared accounting for one multi-chip -> DRAM output tile write-back.
    void account_output_writeback_to_dram();
    // Transfer the data each PE.
    void data_transfer(scheduler_t *m_scheduler);
    // Flush data at temporal buffer in PE array.
    void flush_data();

    // Print out specification of the chip-level processors.
    void print_specification();

    // Modeled busy duration of the multi-chip fabric for the current layer (max of its cost axes).
    double modeled_elapsed_cycles() const;

    // Accumulate leakage energy once for the modeled layer duration.
    void update_static_energy(double elapsed_cycles);

    // Reset stats.
    void reset();

    dram_t *dram;

    data_t *data;                               // Temporal buffer of chip-level processor
    bool   double_buffer;                       //
    // MC2: DBBIF-style outstanding-request limit for the DRAM <-> multi-chip boundary, in
    // TILES IN FLIGHT. 0 means unset -- the boundary keeps the depth its double_buffer flag
    // implies (1 or 2), so an absent key changes nothing. A positive value overrides that
    // depth directly, which is the only way to express "up to N requests outstanding": the
    // double_buffer flag is a boolean and cannot reach past 2.
    unsigned max_outstanding_requests;
   
    std::vector<bool> skip_transfer;              
    bool              equal_output_tile;         

    std::vector<unsigned> tile_size;            // Tile size of chip-level processor
    std::vector<unsigned> offsets;              // Offset index of chip-level processor
    unsigned duplicated_input;

    /* Chip-level processor's specifications */

    /* Chip-level processors signals */
    bool exist_temporal_buffer;
    std::vector<bool> exist_data;               // Data exist in chip-level processor (in temporal buffer)
    std::vector<bool> request_to_dram;          // Request data from the chip-level processors to off-chip memory
    std::vector<bool> is_waiting_data;          // Chip-level processor is waiting data or not

    /* Unit stats */
    double u_transfer_cycle;                    // The unit data transfer cycle between global buffer and chip-level processor (unit NoP cycle)
    double u_transfer_energy;                   // The unit data transfer energy between global buffer and chip-level processor (unit NoP energy)

    std::vector<double> u_read_cycle;           // The unit chip-level processor read cycle (if temporal buffer exist)
    std::vector<double> u_read_energy;          // The unit chip-level processor read energy (if temporal buffer exist)

    std::vector<double> u_write_cycle;          // The unit chip-level processor write cycle (if temporal buffer exist)
    std::vector<double> u_write_energy;         // The unit chip-level processor write energy (if temporal buffer exist)

    std::vector<double> u_static_energy;        // The unit static (leakage) energy of the chip-level processor, pJ/cycle

    /* Chip-level processors stats */
    std::vector<unsigned> num_request;          // The number of data request from the global buffer to chip-level processor
    std::vector<unsigned> num_data_transfer;    // The number of data transfer from the chip-level processor to global buffer

    std::vector<double> access_cycle;           // Total access cycles to the chip-level processor (if temporal buffer exist)
    std::vector<double> access_energy;          // Total access energies to the chip-level processor (if temporal buffer exist)

    std::vector<double> transfer_cycle;         // Total data transfer cycle between the global buffer and chip-level processor (NoP)
    std::vector<double> transfer_energy;        // Total data transfer energy between the global buffer and chip-level processor (NoP)

    std::vector<double> cycle_temporal_chips;

    std::vector<size_t> payload_link_transactions;
    std::vector<size_t> metadata_link_transactions;
    std::vector<size_t> storage_link_transactions;

    double utilization;                         // Utilization of chip-level processor

    // RE1 latency half: the final cast's TIME, kept on its own counter. It also enters
    // write_back_cycle, which reaches the fabric's busy time only as a MAX against the access and
    // transfer axes -- so a cast cheaper than those axes is invisible in the critical path. Keeping
    // and reporting the raw figure is what makes `output_cast_cycle` checkable at all; before this
    // the knob was read from the config and observable nowhere (validation/knobs KN9).
    double output_cast_cycle;
    double write_back_cycle;                    //
    double overlapped_transfer_cycle;           //

    std::vector<double> static_energy;          // Static (leakage) energy of the chip-level processor temporal buffer

    std::vector<unsigned> line_size;            // Line size of temporal buffer
    std::vector<unsigned> mask_bits;            // Mask bits of temporal buffer

protected:
    // P2-3: distinct_chunks is the number of distinct source tiles across the active
    // chip grid for this datatype (1 == fully broadcast/multicast; ==
    // get_number_of_active_chips() == fully partitioned per chip); see the scheduler's
    // {input,weight,output}_offset_multi_chip vector sizes.
    void account_descriptor_dense_distribution(data_type_t type, size_t distinct_chunks);


    std::vector<global_buffer_t*> chips;        // Global buffers to connect

    stationary_type_t stationary_type;          // Stationary type at the chip-level processor (equal to global buffer stationary type)
    std::string parameter_order;                // Order of parameters used in data transfer
    noc_type_t nop_type;                        // NoP type (Mesh)
    // L7: does the NoP forward ONE copy of a broadcast tile along a multicast tree (1) or
    // send a separate copy to every chip (0)? Default 1: on a bus one transmission is
    // physically seen by every receiver, so charging it per chip was never right, and a
    // mesh with multicast routers behaves the same way at the shared package ingress. Set
    // to 0 to model a package whose routers only do unicast.
    bool nop_multicast;
    // RE1: the FINAL output cast/pack to output precision, charged at the off-chip output store.
    //
    // Every earlier boundary fires more than once per output element: a PE writes its tile back
    // once per reduction pass, and the GLB reads it out once per GLB-level reduction pass. Charging
    // the cast at either made its count follow the MAC count (262,144 on GEMM 64x64x64) or a
    // multiple of the output count (16,384) instead of the 4,096 final elements. The off-chip store
    // is the boundary DR6/T1 establish happens exactly once per output element, and it is where the
    // value is committed in OUTPUT precision -- so it is both the correct count and a defensible
    // owner.
    memory_type_t memory_type;                  // Memory type of global buffer

    unsigned height;                            // Y-dimension of chip-level processors
    unsigned width;                             // X-dimension of chip-level processors

    unsigned num_chips;                         // The number of processors in chip-level processors
    unsigned num_active_chips_x;                // The number of active processors in X dimension
    unsigned num_active_chips_y;                // The number of active processors in Y dimension

    float frequency;                            // Frequency at chip-level processor
    float bandwidth;                            // The bandwidth between the global buffer and chip-level processors (temporal buffer) in GB/s
    unsigned bitwidth;              

    // Buffer size
    unsigned memory_size;                       // Temporal buffer size
    // KB, so FRACTIONAL: Eyeriss v2 declares 4.5 KB / 7.5 KB buffers. global_buffer_t has always
    // stored these as double; multi_chip_t stored them as unsigned, so `input_size = 4.5` failed
    // the unsigned parse and SILENTLY left the buffer at 0 -- which made eyerissv2.cfg unrunnable
    // (every tile "exceeded" a zero buffer) and went unnoticed because no gate ran that config.
    // The two components now agree on the type. validation/knobs KN11 rejects the silent-drop
    // class outright.
    double input_size;                          // Temporal input buffer size
    double weight_size;                         // Temporal weight buffer size
    double output_size;                         // Temporal output buffer size

    bool     initial;                           // 
    double nop_cycle;                           // NoP cycle
    double nop_energy;                          // NoP energy

    double scatter_cycle;                       // Scatter network cycle
    double scatter_energy;                      // Scatter network energy

    double local_cycle;                         // Local network cycle
    double local_energy;                        // Local network energy

    double gather_cycle;                        // Gather network cycle
    double gather_energy;                       // Gather network energy
};

#endif

