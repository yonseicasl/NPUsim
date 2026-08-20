#ifndef __PE_ARRAY_H__
#define __PE_ARRAY_H__

#include <vector>
#include <cmath>
#include "pe.h"
#include "global_buffer.h"
#include "scheduler.h"

class global_buffer_t;
class pe_t;

class pe_array_t {

public:

    // Phase-5: the component's clock in MHz. Power needs a single authoritative clock
    // across the modeled components (the timeline is one shared cycle axis), so stats_t
    // compares these and reports power as unsupported when they disagree.
    double clock_mhz() const { return static_cast<double>(frequency); }
    
    pe_array_t(section_config_t m_section_config);
    virtual ~pe_array_t();

    // Initialize the PE array.
    virtual void init(section_config_t m_section_config) = 0;

    // Connect PE array and Global buffer
    void connect(global_buffer_t *m_global_buffer);

    // Update tile size of PE array
    virtual void update_tile_size(scheduler_t *m_scheduler) = 0;

    /* Get the PE array's specifications */

    // Get stationary type (same as PE stationary type)
    stationary_type_t get_stationary_type();
    // Get parameter order.
    std::string get_parameter_order();
    // Get the number of PEs in the PE array.
    unsigned get_number_of_pes();
    // Get the number of active PEs
    unsigned get_number_of_active_pes();

    /* Check PE array's status */

    // Return true when the PE array has no data.
    bool is_idle();
    // Check whether the request exist 
    bool is_exist_request_at_pe();
    // Check whether the request exist in the PE array.
    bool is_exist_request_at_buffer();

    bool is_exist_data();

    /* PE array's actions */

    // Wait for the data.
    bool wait_data();
    // Exist data in PE array.
    void fill_data();
    // Request data to Global buffer.
    void request_data();
    // Transfer the data each PE.
    virtual void data_transfer(scheduler_t *m_scheduler) = 0;
    // Flush data at temporal buffer in PE array.
    void flush_data();
    // Virtual so reduction-network variants (adder tree) can layer their tree cost on
    // top of the generic gather accounting.
    virtual void account_descriptor_dense_writeback(pe_t *source_pe, size_t elements);

    // Print out the configuration of PE array.
    virtual void print_specification() = 0;

    // Accumulate leakage energy of the PE-array temporal buffer once for the layer duration.
    void update_static_energy(double elapsed_cycles);
    // Modeled busy duration of the PE-array temporal buffer for the current layer.
    double modeled_elapsed_cycles() const;

    // Reset the component and stat
    void reset();

    std::vector<pe_t*> pes;                                 // A group of PEs.
    data_t *input_data;
    data_t *weight;
    data_t *output_data;

    std::vector<bool> skip_transfer;                        // Check whether the data should transferred or not.
    bool              equal_output_tile;
   
    std::vector<unsigned> tile_size;                        // Tile size of PE array
    unsigned duplicated_input;                              //

    /* PE array specifications */

    unsigned index;                                         // Index of PE array

    /* PE array signals */
    // A signal to indicate the existence of data.
    bool exist_temporal_buffer;
    std::vector<bool> exist_data;                           // Data exist in PE array (in temporal buffer).
    std::vector<bool> request_to_global_buffer;             // Request data from PE array to Global buffer
    std::vector<bool> is_waiting_data;


    /* PE array unit costs */
    
    std::vector<double> u_read_cycle;                       // The unit PE array read cycle
    std::vector<double> u_read_energy;                      // The unit PE array read energy

    std::vector<double> u_write_cycle;                      // The unit PE array write cycle
    std::vector<double> u_write_energy;                     // The unit PE array write energy

    std::vector<double> u_static_energy;                    // The unit static (leakage) energy of the PE-array temporal buffer, pJ/cycle


    /* PE array stats */

    std::vector<unsigned> num_request;                      // The number of data request from PE to PE array
    std::vector<unsigned> num_data_transfer;                // The number of data transfer from PE array to local buffer

    std::vector<double> access_cycle;                       // Total access cycles to PE array
    std::vector<double> access_energy;                      // Total access energies to PE array

    std::vector<double> transfer_cycle;                     // Total data transfer cycle between PE and PE array
    std::vector<double> transfer_energy;                    // Total data transfer energy between PE and PE array

    std::vector<double> cycle_temporal_pe;

    std::vector<size_t> payload_link_transactions;
    std::vector<size_t> metadata_link_transactions;
    std::vector<size_t> storage_link_transactions;

    double utilization;                                     // Utilization of PE array
    std::vector<double> buffer_utilization;                 // Temporal-buffer occupancy per data type (peak across the layer)

    double write_back_cycle;                                //
    double overlapped_transfer_cycle;                       //

    // V2 (RTL-calibrated): per weight-residency fold fill. Each per-PE temporal weight
    // element is one array-wide residency (fold) in the real WS schedule; every fold
    // costs u_weight_fold_fill_cycle of load/accumulate-drain bubble on the compute
    // schedule. u_layer_setup_cycle is a one-time per-layer schedule setup cost.
    // Both default to 0 (off) so legacy configurations are unchanged.
    double u_weight_fold_fill_cycle;                        // Cycles per weight-element residency (fold)
    double u_layer_setup_cycle;                             // One-time per-layer setup cycles
    // E5: the ENERGY counterparts of those two latency bubbles. Calibrating the latency of a
    // weight reload or a layer setup used to make the schedule more accurate while leaving the
    // same activity FREE on the energy side (leakage aside) -- so the more precise the latency
    // model got, the more asymmetric energy became. A config declares these per event:
    //     weight_fold_fill_energy   per weight-element residency (fold)
    //     layer_setup_energy        once per layer
    // Both default to 0, so existing configs are unaffected. The EVENT COUNTS are recorded
    // separately from the energy so the identity (events x unit cost) is checkable, and the
    // selection rule mirrors the latency side exactly: a calibrated fold cost wins over the
    // analytical drain, and the two are never charged together.
    double u_weight_fold_fill_energy;
    double u_layer_setup_energy;
    // Accumulated fold/setup events and their energy for this layer.
    double weight_fold_events;
    double weight_fold_energy;
    // RE1: accumulator energy handed over by the PEs when edge_accumulation puts the accumulator
    // at the array edge. It is the array's own storage, so it belongs to the array's subtotal.
    double accumulator_energy;
    double fold_fill_cycle;                                 // Accumulated fold-fill + setup cycles (per layer)
    // V3: outputs accumulate at the array edge (Gemmini-style accumulator) instead of
    // residing in per-PE local buffers -- exempt OUTPUT from per-PE/array capacity checks.
    bool edge_accumulation;

    std::vector<double> static_energy;                      // Static (leakage) energy of the PE-array temporal buffer

    std::vector<unsigned> line_size;
    std::vector<unsigned> mask_bits;
    

    // RE6: which links this array's noc_energy prices, for the report. Two conventions give
    // different numbers for the same fabric, so the output has to say which one it used.
    std::string describe_noc_link_contract() const;

    // Whether the array's output tile equals the global buffer's. Under output stationary that
    // makes the tile array-resident, so intermediate write-backs are skipped until the drain. It is
    // reported because it changes what the array->GLB output request count means -- though it is
    // NOT the only reason that count can be 0 (see validation/hierarchy, which records the open
    // question: simba and fsd report 0 with neither edge accumulation nor an equal tile).
    bool output_tile_is_array_resident() const { return equal_output_tile; }

    // E20-2: array-level reduction additions. Only an [adder_tree] array performs any, and
    // its `adder_energy` was unreachable for two reasons at once (validation/knobs KN9); a
    // count is what lets an unpriced-but-active tree be told from an absent one.
    double reduction_additions;

    // E20-3: may a partial sum stay resident in the array between passes? Only when the reduction
    // for a tile completes before the array moves on -- see mapping_table_t::reduction_tiled_above_array().
    // Recorded once per layer, where a scheduler is in hand, because request_data() has none.
    bool psum_retention_valid;
    void set_psum_retention_scope(scheduler_t *m_scheduler);

protected:

    void initialize_temporal_buffer(section_config_t m_section_config);
    // link_fill_cycles: one-time pipeline fill through the fabric's route depth
    // (0 for single-hop fabrics), charged once per transferred data type.
    void account_descriptor_dense_distribution(scheduler_t *m_scheduler,
                                               double link_cycle, double link_energy,
                                               double link_fill_cycles = 0.0);
    // P4-3/SY2: the NoC topology account_descriptor_dense_writeback() charges the
    // output write-back fabric with. Defaults to the configured `noc_type` (unchanged
    // behavior). systolic_array_t overrides this to force MESH, matching how its own
    // data_transfer() override already forces MESH for the load path regardless of
    // the configured noc label -- a systolic array is structurally a 2D grid on both
    // the operand-skew-in and partial-sum-drain-out directions, not just the load side.
    virtual noc_type_t writeback_noc_type() const { return noc_type; }
    // SY2/L9: per-weight-residency tile-boundary bubble a SYSTOLIC array pays for its
    // accumulation pipeline to drain before different weights take effect. 0 for
    // non-systolic arrays (a spatial/adder-tree array has no such pipeline), and unused when
    // the config calibrates weight_fold_fill_cycle from RTL -- see
    // account_descriptor_dense_distribution().
    virtual double weight_fold_bubble_cycles() const { return 0.0; }


    global_buffer_t *global_buffer;                         // Global buffer to connect

    stationary_type_t stationary_type;                      // Stationary type at the PE array (equal to PE stationary type)
    std::string array_parameter_order;                      // Order of DNN parameters
    noc_type_t noc_type;                                    // NoC type
    memory_type_t memory_type;

    unsigned height;                                        // Y-dimension of PE array
    unsigned width;                                         // X-dimension of PE array
    
    unsigned num_pes;                                       // The number of PEs in PE array
    unsigned num_active_pe_x;                               // The number of active PEs in PE array
    unsigned num_active_pe_y;                               // The number of active PEs in PE array

    // Buffer size
    unsigned input_size;                                    // Input buffer size
    unsigned weight_size;                                   // Weight buffer size
    unsigned output_size;                                   // Output buffer size

    float frequency;                                        // The frequency.
    float bandwidth;                                        // The bandwidth between the PE array and Global buffer in GB/s.
    unsigned bitwidth;


    bool     initial;

    double noc_cycle;                                       // NoC cycle
    double noc_energy;                                      // NoC energy

};


#endif

