#ifndef __NPU_H__
#define __NPU_H__

#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>


#include "convolutional.h"
#include "fully_connected.h"
#include "recurrent.h"

#include "adder_tree.h"
#include "spatial_arch.h"
#include "systolic_array.h"
#include "pe_array.h"

#include "global_buffer.h"

#include "multi_chip.h"
#include "workload_graph.h"
#include "dram.h"
#include "sfu.h"
#include "decomp.h"

#include "mapping_table.h"
#include "scheduler.h"
#include "stats.h"

class network_t;
class scheduler_t;

class npu_t {

public:
    npu_t();
    ~npu_t();

    // Initialize the simulation environment.
    void init(const std::string m_accelerator_config, const std::string m_network_config, const std::string m_mapping_config);
    // Connect components
    void connect();
    // Execute simulation 
    void run(const std::string m_accelerator_config, const std::string m_network_config);

    // Check if the accelerator is idle or not.
    bool is_idle();

    /* Operation at accelerator components. */

    // Computation at PEs.
    void execute();
    // Data transfer from PE array to local buffers in PEs.
    void transfer_data_to_pe();
    // Transfer tiled data from global buffer to PE array.
    void transfer_data_to_pe_array();
    // Data transfer from chip to global buffers.
    void transfer_data_to_global_buffer();
    // Transfer tiled data from DRAM to Multi Chip
    void transfer_data_to_multi_chip();

    // Send data request from Multi Chip to DRAM.
    void request_to_dram();
    // Request from the global buffer to Multi Chip
    void request_to_multi_chip();
    // Send data request from PE array to the global buffer.
    void request_to_global_buffer();
    // Send data request from PE to PE array.
    void request_to_pe_array();


    // Print out the Accelerator specification.
    void print_accelerator_specification();
    // Print out DNN configuration.
    void print_network_configuration(unsigned m_layer_index, unsigned m_stats_index);
    // Print out the stats.
    void print_stats(const std::string m_accelerator_config, const std::string m_network_config, unsigned m_index);

    /* Print result of the simulation */

    // Attach executable-IR provenance to every result artifact.
    void print_workload_provenance(std::ofstream &m_output) const;
    void bind_executable_mappings();
    // Print out the simulation result.
    void print_layerwise_results(const std::string m_accelerator_config,
                                 const std::string m_network_config,
                                 unsigned m_layer_index, unsigned m_stats_index);
    // Print out the simulation result.
    void print_total_result(const std::string m_accelerator_config, const std::string m_network_config);

    // Reset performance counters and stats.
    void reset();
    // Update tile size for executing next layer
    void update_tile_size();

protected:
    // SFU (plan/plan_sfu.md): fire the fused-activation cost event for a finished
    // convolution/connected layer -- once per valid output element, AFTER repetition
    // scaling. Without an [sfu] section it only marks nonlinear activations as
    // out-of-scope (legacy numbers unchanged).
    void apply_fused_sfu_activation(unsigned m_layer_index, unsigned m_stats_index);
    // Standalone softmax layer executed on the SFU's multi-pass microprogram (Phase 7).
    void run_standalone_softmax(unsigned m_index, const std::string &m_accelerator_config,
                                const std::string &m_network_config);
    // Weight decompression (evaluation.md Sec 4): compute the layer's dense weight
    // footprint from the mapping and hand it to the engine, BEFORE repetition scaling.
    // No-op without a [decomp] section.
    void apply_weight_decompression(unsigned m_stats_index);
    // KV-cache read (evaluation.md Sec 4): inject the decode step's KV-cache DRAM read on
    // this layer. No-op without a [kvcache] section.
    void apply_kv_cache_read(unsigned m_stats_index);
#ifdef FUNCTIONAL
    // Functional verification: snapshot the accelerator-computed output of a mapped layer,
    // run Nebula's forward() (the single functional owner of bias + activation), and
    // compare element-by-element. Non-mapped layers just forward(). Prints a per-layer
    // verdict and feeds the run-level summary.
    void verify_functional_layer(unsigned m_index, bool m_mapped);
    size_t functional_layers_checked;
    size_t functional_layers_failed;
#endif
    // Phase-7: cost of streaming the softmax operand tensor between the memory hierarchy
    // and the SFU, per [sfu] softmax_operand_residency, from the live components' unit
    // costs (dram: DRAM device + off-chip link + GLB staging/feed ports; glb: GLB
    // feed/result ports only, with a capacity fail-fast).
    sfu_operand_stream_t softmax_operand_stream(size_t m_elements);
    sfu_operand_stream_t graph_operand_stream(const workload_operation_t &m_operation,
                                              const workload_residency_plan_t &m_plan);
    void run_standalone_graph_operation(unsigned m_index,
                                        workload_residency_plan_t m_plan,
                                        const std::string &m_accelerator_config,
                                        const std::string &m_network_config);
    void override_executable_layer_geometry();
    unsigned num_processors;                        // The number of on-chip processors.
    unsigned num_pes;                               // The number of processing elements for each processors.
    compression_type_t compression_type;            // Compression type : Dense, CSR, CSC, SparseMap
    unsigned num_skipped_timing_layers;             // Layers excluded from accelerator timing.
    // Validate physical component counts before connecting the hierarchy.
    void validate_accelerator_components();
    // Validate mapping-selected active components before a layer starts.
    void validate_active_components();


    /* Accelerator components */
    std::vector<pe_array_t*> pe_arrays;             // PE array
    std::vector<global_buffer_t*> global_buffers;   // Global buffer
    multi_chip_t *multi_chip;                       // On-chip processors
    dram_t *dram;                                   // DRAM
    std::vector<sfu_t*> sfus;                       // Per-chip SFU (empty without [sfu])
    decomp_t *decomp;                               // Weight-decompression engine (NULL without [decomp])
    kvcache_t *kvcache;                              // KV-cache read-traffic model (NULL without [kvcache])
    workload_graph_t *workload;                    // Framework-neutral executable IR, if used.
    workload_lifetime_t *workload_lifetime;        // DAG tensor liveness/GLB residency state.
    bool executable_ir_mode;

	nebula::network_t *network;                     // DNN model obtained from the software framework (PyTorch and Nebula)
	std::vector<mapping_table_t*> mapping_tables;	// Mapping tables.

    nebula::layer_t *layer;                         // Neural layers obtained from the software framework (PyTorch and Nebula)
    std::vector<scheduler_t*> schedulers;           // A set of schedulers.
    scheduler_t *scheduler;
    std::vector<stats_t*> layer_stats;
    // Stats of standalone SFU layers (softmax has no mapping section, so its stats live
    // outside the mapping-indexed layer_stats vector).
    std::vector<stats_t*> sfu_layer_stats;
    stats_t *network_stats;


};

#endif
