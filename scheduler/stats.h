#ifndef __STATS_H__
#define __STATS_H__


#include <string>
#include <vector>
#include "mapping_table.h"

#include "pe_array.h"

class stats_t {
	
public:
    stats_t();
    ~stats_t();
    // Initialize the stats.
    void init();

    // Print out the stats of accelerator and neural network.
    void print_stats();
    void print_stats(std::ofstream &m_output_file);

    // Update tile size
    void update_tile_size(scheduler_t *m_scheduler);

    // Update layer-wise stats.
    void update_stats(std::vector<pe_array_t*> m_pe_array, std::vector<global_buffer_t*> m_global_buffer, multi_chip_t *m_multi_chip, dram_t *m_dram);

    // Serially account for identical temporal tiles at or above the global buffer.
    // m_datatype_repetitions: per-datatype GLB repetition factors (see
    // mapping_table_t::datatype_repetitions) -- off-chip (multi-chip/DRAM) traffic
    // scales with these instead of the full repetition count, because repetitions
    // over dimensions a tensor does not depend on revisit tiles the on-chip
    // buffers retain.
    void merge_global_buffer_fill();
    // CE1/CE2/CE4: recompute stage busy, critical path, bottleneck inputs, leakage
    // window, and MAC availability from the FINAL (repetition-scaled) cycle vectors.
    void finalize_layer_timeline();
    void scale_serial_repetitions(unsigned m_repetitions,
                                  const std::vector<unsigned> &m_datatype_repetitions,
                                  const input_halo_reuse_t &m_input_halo = input_halo_reuse_t(),
                                  bool m_halo_capacity_sufficient = false);

    // Update network stats.
    void update_network_stats(stats_t *m_source);

    // Print out the result of simulation.
    void print_results(std::ofstream &m_output_file);
    // E1: component energy subtotals and the layer total. See the no-double-counting boundary
    // documented at the definition.
    void print_energy_summary(std::ofstream &m_output_file);
    // Phase-5: average power, EDP/ED2P and the explicit scope of what power does NOT include.
    void print_power_summary(std::ofstream &m_output_file, double dynamic_total,
                             double static_total);
	
    /* Layer timeline (A1): per-level busy time on a shared analytical timeline and the
       segment-combined critical-path latency (overlap at double-buffered boundaries). */
    double layer_latency;
    double busy_cycle_pe;
    // CE7: per-stage busy composition (0=DRAM,1=multi-chip,2=GLB,3=PE array,4=PE).
    double stage_axis_access[5];
    double stage_axis_link[5];
    double stage_axis_overlap[5];
    double stage_axis_compute;                 // PE compute schedule incl. fold fill (CE2)
    double stage_axis_format;                  // P4-4: PE format-IP axis (payload/metadata/spill)
    bool timeline_boundary_overlap[4];         // stage-boundary overlap flags captured pre-scale
    double timeline_physical_macs;             // physical scalar MACs for availability
    // P3: temporal GLB-repetition ("tile") count, set from scale_serial_repetitions()'s
    // m_repetitions. Lets finalize_layer_timeline() combine an overlapping run of
    // stage-busy totals via a fill+bottleneck pipeline makespan instead of an
    // unrealistic flat max() (see temporal_pipeline_run_cycles()).
    unsigned temporal_repetition_tiles;
    // L5: how many of those tiles the OFF-CHIP stages (DRAM, multi-chip) actually serve.
    // Off-chip traffic scales per datatype, so a GLB repetition over a dimension a tensor
    // does not depend on revisits a tile the on-chip buffers still hold and does NOT
    // re-fetch it. Taking the max across datatypes gives the number of tiles on the shared
    // tile clock that carry off-chip work; the rest carry zero cost for those stages. This
    // is what lets a rate-mismatched stage take part in the timeline at all -- the previous
    // closed form had to fall back to a flat max() whenever a run touched stage 0 or 1.
    unsigned offchip_repetition_tiles;
    // E20-3: dense input-window reuse at the off-chip boundary. The logical tile requests remain
    // visible, but their overlapping payload is coalesced to the exact union footprint when the
    // GLB can retain the required sliding working set.
    bool input_halo_overlap;
    bool input_halo_capacity_sufficient;
    bool input_halo_reuse_applied;
    size_t input_halo_unique_elements;
    size_t input_halo_replicated_elements;
    size_t input_halo_working_set_bytes;
    size_t input_halo_pre_dram_transactions;
    unsigned dram_input_link_bits;
    // Exact source-access reconstruction for a coalesced input union. Scaling an already rounded
    // descriptor aggregate by a floating ratio makes the physical access count fractional and can
    // break the read-energy/read-cycle unit identity even when both print the same rounded value.
    unsigned dram_input_source_line_bits;
    double dram_input_read_cycle;
    double dram_input_read_energy;
    bool dram_input_access_hidden;

    // L5: per-stage time spent blocked specifically by a full downstream buffer, from the
    // final per-tile timeline. Reported so the reader can see WHERE the pipeline stalled,
    // not only how long the layer took.
    double timeline_stall[5];
    // L5/L6: the staging-buffer depth the timeline used at each stage boundary, in tiles
    // in flight (1 = single buffer / no tile-level decoupling, 2 = double buffer). Reported
    // because it IS the overlap contract: it decides whether the two sides of a boundary
    // pipeline or alternate, and for the multi-chip -> GLB boundary it also carries the
    // bypass rule (a GLB-bypassed datatype has no buffer there, so that boundary cannot
    // decouple no matter what the GLB's double_buffer flag says). At network scope the min
    // across layers is reported -- a boundary is only as decoupled as its worst layer.
    unsigned timeline_boundary_depth[4];
    double busy_cycle_pe_array;
    double busy_cycle_global_buffer;
    double busy_cycle_multi_chip;
    double busy_cycle_dram;

    // CE4: max_entity(sum/max_type(...)) combined directly while iterating entities
    // (chip, or chip+PE) in update_stats(), before repetition scaling and before the
    // entity dimension is collapsed into the flat per-datatype vectors below. Consumed
    // by finalize_layer_timeline() instead of re-deriving stage_axis_* from the flat
    // vectors, which would incorrectly sum/max each entity's own per-type PEAK instead
    // of taking the max across entities of each entity's own per-type TOTAL.
    double entity_combined_access_global_buffer;   // access_cycle_global_buffer axis, fill excluded (see stage_fill_access_global_buffer)
    double entity_combined_link_global_buffer;      // transfer_cycle_global_buffer axis
    double entity_combined_overlap_global_buffer;   // cycle_pe_array_global_buffer axis
    double entity_combined_access_pe_array;         // access_cycle_pe_array axis
    double entity_combined_link_pe_array;           // transfer_cycle_pe_array axis
    double entity_combined_overlap_pe_array;        // cycle_temporal_pe_array axis
    double entity_combined_access_lb;               // access_cycle_lb axis
    double entity_combined_link_pe;                 // transfer_cycle_pe axis
    double entity_combined_overlap_mac_lb;           // cycle_mac_lb axis
    // P4-4: format_cycle_pe axis (payload/metadata/spill format-IP). Previously only
    // fed the PE-local max() in pe_t::modeled_elapsed_cycles(), which only feeds an
    // inert pre-scale estimate -- never the authoritative busy_cycle_pe/layer_latency
    // computed in finalize_layer_timeline(), so a slow format-IP was invisible to the
    // critical path. Now combined into busy_cycle_pe like the other PE axes.
    double entity_combined_format;
    // The mc->GLB fill (write) side scales per-datatype (not uniformly), so its type
    // combination is captured in merge_global_buffer_fill() -- right before it is folded
    // into access_cycle_global_buffer and zeroed -- at its FINAL scaled value.
    // Reported/diagnostic only: max_chip(combine_type(fill)). The authoritative GLB
    // access axis is entity_combined_access_global_buffer, which merge_global_buffer_fill()
    // recomputes as max_chip(base[chip] + combine_type(fill[chip])).
    double stage_fill_access_global_buffer;
    // CE4/P1-D: the GLB entity dimension must survive until AFTER repetition scaling,
    // because the base (GLB->PE-array read) side scales uniformly while the fill
    // (multi-chip->GLB write) side scales per datatype. Collapsing chips first --
    // fill[type] = max_chip(fill[chip][type]) and only then combining types -- produces
    // sum_type(max_chip(...)), which on a shared GLB can invent elapsed time that never
    // existed (chip 0 input-bound at 100 and chip 1 weight-bound at 100 report 200 while
    // the two chips actually run those fills in parallel in 100). So keep each chip's own
    // values and combine in the correct max_chip(combine_type(...)) order at the end.
    // L1: BOTH sides stay per-chip AND per-datatype. Type-combining either side early
    // destroys the datatype correspondence the two sides need: on a SEPARATE buffer each
    // datatype has its own partition and the partitions run in parallel, so
    // max_type(base) + max_type(fill) can sum the peaks of two DIFFERENT partitions.
    // The axis must be max_chip(combine_type(base[chip][t] + fill[chip][t])).
    std::vector<std::vector<double>> chip_access_cycle_global_buffer;      // per chip/datatype: BASE (GLB->PE-array read) access, scales uniformly
    std::vector<std::vector<double>> chip_fill_access_cycle_global_buffer; // per chip/datatype: FILL (multi-chip->GLB write) access, scales per datatype

    /* Tile size */
    std::vector<std::vector<unsigned>> tile_size;

    /* PE */
    memory_type_t local_buffer_type;
    // LB7: true only when every active PE's local buffer is double-buffered; a single
    // single-buffered PE forces the whole PE stage's compute<->LB axes to serialize.
    bool pe_double_buffer;
    size_t num_computation;                                             // Number of computations.
    double computation_cycle;                                           // Total computation cycle.
    double max_computation_cycle;
    double min_computation_cycle;
    double avg_computation_cycle;
    double computation_energy;                                          // Total computation energy.
    double mac_busy_cycle;                                              // Sum of scalar-MAC busy cycles.
    double mac_available_cycle;                                         // Sum of scalar-MAC available cycles.
    std::vector<double> format_cycle_pe;                                // PE format-IP cycle by tensor type.
    std::vector<double> format_energy_pe;                               // PE format-IP energy by tensor type.
    // E4: the lane->accumulator reduction's own energy axis (split out of computation_energy so
    // the adder-tree work is visible), and the two format event counts kept apart because they
    // move DIFFERENT precisions -- accumulator-format spills vs output-format casts.
    double reduction_energy;
    // E5: fold and setup dynamic energy, kept APART because they scale differently -- exactly as
    // their latency counterparts do. A weight residency repeats with every GLB repetition; the
    // per-layer schedule setup fires once. The event counts are recorded alongside so the
    // identity (events x unit cost) is checkable from the report.
    double weight_fold_energy;
    double layer_setup_energy;
    double weight_fold_events;
    double layer_setup_events;
    // RE1: the four accumulator/output events, each from its own boundary. The create count is
    // reported so the deliberately-free zero-init path is visible rather than implied.
    size_t accumulator_reload_bytes;
    size_t accumulator_spill_bytes;
    size_t accumulator_create_events;
    size_t accumulator_retained_events;      // E20-4: retained, not reloaded
    size_t output_cast_bytes;
    double accumulator_energy;
    double output_cast_energy;
    double output_cast_cycle;
    // E20-2: activity counts for the events whose cost keys are optional. An energy of 0 cannot
    // say whether the event happened; these can.
    size_t row_activation_events;
    size_t format_payload_events;
    size_t format_metadata_events;
    double reduction_additions;

    // E20-1/E20-2: events that FIRED while the key pricing them was absent from the config.
    // Crossing an activity count with the declaration is the only way to see this: the energy axis
    // reads 0 whether the event was free, absent, or simply unpriced. Every entry is a charge the
    // total is missing, of an amount nobody stated -- so a total carrying any of them is not
    // absolute, and a wattage derived from it is not either.
    std::vector<std::string> unpriced_active_events() const;
    // RE1: accumulator energy the PEs handed to the PE array because edge_accumulation puts the
    // accumulator at the array edge. Reported on the array's row, not the PE's.
    double pe_array_accumulator_energy;

    std::vector<unsigned> num_request_pe;                               // Number of request to local buffer of PE.
    std::vector<unsigned> num_data_transfer_pe;                         // Number of data transfer to MAC unit of PE.

    std::vector<double> access_cycle_mac;                               // Total access cycle to MAC unit of PE.
    std::vector<double> max_access_cycle_mac;
    std::vector<double> min_access_cycle_mac;
    std::vector<double> avg_access_cycle_mac;
    std::vector<double> access_energy_mac;                              // Total access energy to MAC unit of PE.
    double utilization_mac;                                             // Time-based utilization of physical scalar MACs.

    std::vector<double> access_cycle_lb;                                // Total access cycle to local buffer of PE.
    std::vector<double> max_access_cycle_lb;
    std::vector<double> min_access_cycle_lb;
    std::vector<double> avg_access_cycle_lb;
    std::vector<double> access_energy_lb;                               // Total access energy to local buffer in PE.
    std::vector<double> utilization_local_buffer;                       // Utilization of the local buffer in PE
    double total_utilization_local_buffer;

    std::vector<double> transfer_cycle_pe;                              // Total data transfer cycle between MAC unit and local buffer of PE.
    std::vector<double> transfer_energy_pe;                             // Total data transfer energy between MAC unit and local buffer of PE.
    std::vector<size_t> payload_link_transactions_pe;
    std::vector<size_t> metadata_link_transactions_pe;
    std::vector<size_t> storage_link_transactions_pe;
    std::vector<double> cycle_mac_lb;                                   // Overlapped cycle between computing and local buffer

    std::vector<double> static_energy_pe;                               // Static energy of PE
    std::vector<double> static_energy_pe_array;                         // Static (leakage) energy of the PE-array temporal buffer

    /* PE array */
    std::vector<unsigned> num_request_pe_array;                         // Number of request to PE array (from PE).
    std::vector<unsigned> num_data_transfer_pe_array;                   // Number of data transfer from PE array (to PE).

    std::vector<double> access_cycle_pe_array;                          // Total access cycle to PE array.
    std::vector<double> access_energy_pe_array;                         // Total access energy to PE array.
    double utilization_pe_array;
    std::vector<double> utilization_pe_array_buffer;                    // PE-array temporal-buffer occupancy per data type.
    double fold_fill_cycle_pe_array;                                    // V2: weight-residency fold fill + per-layer setup on the compute schedule.
    double layer_setup_cycle_pe_array;                                  // V2: one-time per-layer setup (added once, after repetition scaling).

    std::vector<double> transfer_cycle_pe_array;                        // Total data transfer cycle between PE and PE array.
    std::vector<double> cycle_temporal_pe_array;                        // PA4: PE-array temporal-buffer pipelined-overlap cycle.
    std::vector<double> transfer_energy_pe_array;                       // Total data transfer energy between PE and PE array.
    std::vector<size_t> payload_link_transactions_pe_array;
    std::vector<size_t> metadata_link_transactions_pe_array;
    std::vector<size_t> storage_link_transactions_pe_array;

    /* Global buffer */
    memory_type_t global_buffer_type;
    // P1-B: which datatypes bypass the GLB SRAM. A bypassed datatype is streamed from the
    // chip's NoP ingress straight to the PE array: it costs no GLB fill (write) and no GLB
    // read, but it DOES traverse the GLB<->PE-array fabric and land in the PE-array
    // temporal buffer, so its link/destination cost is charged like any other stream and
    // shows up in the PE-array <-> GLB transaction counters. Reported so a consumer that
    // wants SRAM-mediated traffic only (e.g. the Eyeriss GLB-traffic reference, which
    // counts SRAM accesses) can exclude exactly those datatypes instead of guessing.
    std::vector<bool> global_buffer_bypass;
    // L8: the analytical DRAM contract in force and what it does NOT capture, captured from
    // dram_t in update_stats() and printed with the DRAM results. The audit's completion
    // condition for the analytical path is that its supported scope is stated in the output
    // rather than left to be inferred -- an idealized bank-interleaved row model must not be
    // read as a bank-conflict model.
    // E3: which config key supplied the MAC energy, and whether it was declared for the operand
    // precision actually in use. Reported with the energy summary so a compute-energy number is
    // never read as precision-aware when the same scalar would have been used for any precision.
    std::string compute_energy_basis;
    bool compute_energy_precision_calibrated;
    std::string operand_precision;
    // E10/Phase-5: the clock contract that turns cycles into seconds, and therefore energy into
    // power.
    //
    // The layer timeline places every stage on ONE shared cycle axis (see pipeline_timeline_cycles),
    // so a cycle count is only convertible to time if every modeled component runs on the same
    // clock. When they do, that frequency is authoritative and power is well defined. When they do
    // NOT, the honest answer is that power is unsupported for this config: dividing a mixed-domain
    // cycle count by any one component's frequency would invent a number. Reported either way, so
    // the reader is never left guessing which case they are in.
    double authoritative_frequency_mhz;
    bool single_clock_domain;
    std::string clock_domain_note;
    std::string dram_timing_model;
    // RE6: which links `noc_energy` prices, stated in the report. Two conventions give different
    // numbers for the same fabric (a spanning tree over N endpoints has N-1 edges, while counting
    // one link per receiver gives N), and nothing in the output said which one produced the
    // energy -- so the number could not be checked against a physical link count at all.
    std::string array_noc_link_contract;
    bool output_tile_array_resident;   // see pe_array_t::output_tile_is_array_resident()
    std::string nop_link_contract;
    std::string dram_timing_limits;
    // L3: true only on the network rollup object. Layer and network scope have DIFFERENT
    // axis contracts and the report must say which one it is printing:
    //   layer   : busy = max(its own axes) (PE also folds in compute + format) -- the
    //             printed axes RECONSTRUCT busy.
    //   network : busy[stage] = sum_layer(layer busy[stage]) and axis = sum_layer(layer axis),
    //             two independent sums. sum_layer(max(...)) != max(sum_layer(...)) in
    //             general, so the printed network axes are per-axis WORK TOTALS and do NOT
    //             reconstruct network busy. Replacing network busy with max(network axes)
    //             would be worse, not better: layers run serially, so it would undercount
    //             the actual busy time of every stage.
    bool network_rollup;
    // L11/P4-14: the timing SCOPE of a network rollup -- how many layers were folded in and how
    // many the accelerator timing model does not support (pooling/activation/normalization).
    // Printed inside the layer-timeline block, next to the latency it qualifies: the trailing
    // warning at the end of the file was too far from the number to stop anyone reading a
    // partial rollup as an end-to-end network latency.
    unsigned network_timing_layers;
    unsigned excluded_timing_layers;
    std::vector<unsigned> num_request_global_buffer;                    // Number of request to global buffer (from PE array).
    std::vector<unsigned> num_data_transfer_global_buffer;              // Number of data transfer from global buffer (to PE array).
    // Link transactions are reported as logical payload/metadata and physical storage streams.
    std::vector<size_t> payload_link_transactions_global_buffer;
    std::vector<size_t> metadata_link_transactions_global_buffer;
    std::vector<size_t> storage_link_transactions_global_buffer;

    std::vector<double> access_cycle_global_buffer;                     // Total access cycle to global buffer.
    std::vector<double> fill_access_cycle_global_buffer;               // mc->GLB fill (write) side; scaled per datatype.
    std::vector<double> fill_access_energy_global_buffer;
    std::vector<double> access_energy_global_buffer;                    // Total access energy to global buffer.
    std::vector<double> cycle_pe_array_global_buffer;                   // Overlapped cycle between the global buffer and the PE array
    std::vector<double> utilization_global_buffer;                      // Global buffer utilization
    double total_utilization_global_buffer;

    std::vector<double> transfer_cycle_global_buffer;                   // Total data transfer cycle between global buffer and PE array.
    std::vector<double> transfer_energy_global_buffer;                  // Total data transfer energy between global buffer and PE array.

    std::vector<double> static_energy_global_buffer;                    // Static energy the global buffer

    /* Multi chip */
    std::vector<unsigned> num_request_multi_chip;                       // Number of request to Multi-chip (from global buffer).
    std::vector<unsigned> num_data_transfer_multi_chip;                 // Number of data transfer from Multi-chip (to global buffer).
    std::vector<size_t> payload_link_transactions_multi_chip;
    std::vector<size_t> metadata_link_transactions_multi_chip;
    std::vector<size_t> storage_link_transactions_multi_chip;

    std::vector<double> access_cycle_multi_chip;                        // Total access cycle to Multi-chip.
    std::vector<double> access_energy_multi_chip;                       // Total access energy to Multi-chip.
    double utilization_multi_chip;
    std::vector<double> utilization_multi_chip_buffer;                   // NoP temporal-buffer occupancy per data type.

    std::vector<double> transfer_cycle_multi_chip;                      // Total data transfer cycle between Multi-chip and global buffer.
    std::vector<double> transfer_energy_multi_chip;                     // Total data transfer energy between Multi-chip and global buffer.
    std::vector<double> static_energy_multi_chip;                       // Static (leakage) energy of the Multi-chip temporal buffer.

    /* DRAM */
    std::vector<unsigned> num_request_dram;                             // Number of request to DRAM (from Multi-chip).
    std::vector<unsigned> num_data_transfer_dram;                       // Number of data transfer from DRAM (to Multi-chip).
    std::vector<size_t> payload_link_transactions_dram;
    std::vector<size_t> metadata_link_transactions_dram;
    std::vector<size_t> storage_link_transactions_dram;

    std::vector<double> access_cycle_dram;                              // Total access cycle to DRAM.
    std::vector<double> access_energy_dram;                             // Total access energy to DRAM.
    std::vector<double> cycle_chip_dram;

    std::vector<double> transfer_cycle_dram;                            // Total data transfer cycle between DRAM and Multi-chip.
    std::vector<double> transfer_energy_dram;                           // Total data transfer energy between DRAM and Multi-chip.

};

#endif
