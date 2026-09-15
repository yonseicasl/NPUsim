#ifndef __NPU_H__
#define __NPU_H__

#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <fstream>
#include <map>
#include <set>


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
#ifdef FUNCTIONAL
#include "functional_artifact.h"
#endif

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
    // Fixture-driven functional path (correctness plan Phase 1, lite): bypass nebula's
    // OpenCV image loader and nebula-forward oracle entirely. `functional_input_buffer` is
    // an OWNED input tensor injected as layer 0's input; `functional_golden` is the external
    // CPU-computed reference for the final mapped layer's output. Both are raw little-endian
    // fp32 files named by [data] functional_input / functional_golden. Empty => fall back to
    // the legacy nebula load_data + forward() differential oracle.
    std::vector<float> functional_input_buffer;
    std::vector<float> functional_golden;
    bool functional_external_golden;
    // G4 (gaps plan Step 3): OPTIONAL per-layer goldens so a multi-layer DAG localizes its
    // first mismatching operation instead of only failing at the final output. [data]
    // functional_golden<i> compares layer i AFTER its finalize (bias/BN/activation);
    // [data] functional_golden_raw<i> compares layer i's RAW accumulator BEFORE the
    // finalize (stage separation, plan §5.3). The final golden remains mandatory.
    std::map<unsigned, std::vector<float>> functional_layer_golden;
    std::map<unsigned, std::vector<float>> functional_layer_golden_raw;
    // Index of the LAST mapped (conv/connected) layer. With an external golden the golden is
    // the network's final output, so only this layer is compared; earlier mapped layers are
    // still executed (their output feeds the next layer) but not verified. -1 = uncomputed.
    int functional_last_mapped;
    // INT8 requantization (plan §7), applied in the finalize after bias/activation when
    // requant_shift>0: round-half-up arithmetic right shift then clamp to [min,max]. Values
    // are exact integers in the float datapath, so this is bit-exact. 0 = disabled.
    int functional_requant_shift;
    int functional_requant_min;
    int functional_requant_max;
    // PER-CHANNEL requant multipliers (plan §7 per-channel scale). When non-empty, requant is
    // out = clamp((acc*mult[c] + round) >> shift, min, max) with a per-output-channel int32
    // multiplier mult[c] (TFLite-style fixed point). Done in int64 so the acc*mult product is
    // exact beyond float's 2^24. Empty => scalar multiplier 1 (plain >>shift). fp32-stored ints.
    std::vector<float> functional_requant_mult;
    // Low-precision OUTPUT format (plan §7 fp16/bf16): round each finalized output to the
    // reduced-mantissa grid. "" = fp32 (no rounding), "bf16", or "fp16". Operands are rounded
    // fixture-side (stored exactly); accumulation stays fp32 (mixed-precision policy).
    std::string functional_output_format;
    // Asymmetric-input INT8 zero-point (plan §7). With a symmetric weight, the raw MAC
    // Sum(qi*qw) is corrected by -zp_i * Sum_k(weight[n][k]) per output channel n (a
    // per-channel constant, so the simulator folds it into the finalize like bias). 0 = off.
    // Requires the standard [N][K] weight layout (not the kfold reduction-tile-major one).
    int functional_input_zero_point;
    // Weight zero-point (full asymmetric int8). Adds the per-output-ROW term -zp_w*Sum_k(input
    // row m) plus the constant +Kdim*zp_i*zp_w to the finalize correction (the row term does NOT
    // fold into per-channel bias, so the simulator computes it from the injected input rows). 0=off.
    int functional_weight_zero_point;
    // G2 (gaps plan Step 2): the EFFECTIVE functional arithmetic profile of this run,
    // recorded in the report. Derived from the fixture's [data] functional_semantics (or
    // inferred from its quantization/rounding knobs) and cross-checked against the
    // accelerator config's input/weight/output_format declarations: a config that declares
    // int8 tensors refuses an fp32 fixture unless the fixture explicitly declares
    // functional_semantics = fp32_reference.
    std::string functional_semantics;
    // DATE2027 zero-gating energy correction: enabled by `functional_zero_gating = 1` in
    // the accelerator config's PE-array section. For every mapped layer the run measures
    // the REAL input-tensor zero fraction (the functional values feeding this layer) and
    // scales the layer's MAC + weight-spad dynamic energy by (1 - fraction) -- the
    // Eyeriss data-gating semantics (JSSC'17 Sec. V-C), energy-only. Per-layer fractions
    // are recorded for the JSON report.
    bool functional_zero_gating;
    std::map<unsigned, double> functional_gating_fraction;
    // B-6: inference accuracy of a real classification functional run (top-1/top-5 hits
    // over the batch), reported at the end of run() and available to callers.
    size_t functional_top1 = 0;
    size_t functional_top5 = 0;
    size_t functional_accuracy_samples = 0;
    // A-3: EXACT per-MAC ifmap-zero fraction, weighted by actual MAC participation.
    // The reference kernels (functional_conv_im2col / functional_gemm_fallback) walk every
    // valid MAC and record (zero-ifmap MACs / total MACs) here. For a dense GEMM this
    // equals the unweighted input-tensor fraction (every input element drives N MACs), but
    // for convolution it differs: border ifmap elements drive fewer MACs, and padding
    // positions are gated (no spad read) exactly as the chip handles them -- so this is the
    // fraction the data-gating energy model should use. -1 = kernel did not record it
    // (datapath-value layer), in which case the unweighted tensor fraction is used.
    std::map<unsigned, double> functional_exact_gating_fraction;
    // Declared weight layout of the fixture ([data] weight_layout): "" = standard [N][K],
    // "ktile" = reduction-tile-major [Kf][N][sK] (required by an INPUT_CHANNEL temporal fold,
    // incompatible with the zero-point corrections that read weight as [N][K]). The mapping
    // alone cannot reveal the layout, so the fixture declares it and G5 cross-checks.
    std::string functional_weight_layout;
    // Machine-readable per-verified-layer comparison records (plan §5.2), one JSON object each.
    std::vector<std::string> functional_report;
    // G5 (gaps plan Step 1, extended): classify the mapped layer's mapping against the
    // datapath VALUE envelope. Returns false when the datapath (+ replay) computes the
    // values; true when the mapping is outside the envelope (DRAM-queue fold, filter
    // fold, mixed chip axis, INPUT_CHANNEL fold without the ktile layout) and the values
    // must come from a mapping-independent reference kernel instead -- timing still comes
    // from the mapped datapath run. Only a combination no kernel can serve (the ktile
    // weight layout together with an out-of-envelope mapping) is rejected outright.
    bool functional_mapping_needs_kernel(unsigned m_index);
    // G3 (gaps plan Step 5): in-simulator im2col value kernel for mapped convolution
    // layers with P/Q > 1 (or an out-of-envelope conv mapping) -- the native conv offset
    // network cannot compute them, so the VALUE path lowers the conv to a deterministic
    // scalar im2col GEMM (zero padding, stride, groups, in-accumulation zero points),
    // exactly the lowering a GEMM accelerator performs. Timing still comes from the
    // mapped datapath run; layers computed here are tagged "im2col" in the report.
    void functional_conv_im2col(unsigned m_index);
    // Reference GEMM kernel for a mapped CONNECTED layer whose mapping is outside the
    // datapath value envelope: out[M][N] = in[M][K] @ W[N][K]^T raw accumulators
    // (position-major, standard weight layout); the shared finalize applies
    // BN/zero-point/bias/activation/requant. Tagged "gemm" in the report.
    void functional_gemm_fallback(unsigned m_index);
    std::map<unsigned, std::string> functional_kernel_layers;   // layer -> kernel tag
    void verify_against_golden(unsigned m_index);
    // Compare one layer's output_data against a specific golden buffer at a named stage
    // ("final" | "layer" | "raw"). Feeds the summary counters and the JSON report.
    void verify_buffer_against(unsigned m_index, const std::vector<float> &m_golden,
                               const char *m_stage);
    // Raw-pointer form for buffers that do not live in a nebula layer (the executable-IR
    // tensor store); zero statistics still read the layer's weight/input when present.
    void verify_buffer_against(unsigned m_index, const float *m_actual, size_t m_elements,
                               const std::vector<float> &m_golden, const char *m_stage);

    // G1 (gaps plan Step 4): executable-IR functional execution. Values live in a tensor
    // store keyed by executable tensor id (aliases resolve to their storage tensor); the
    // npusim.tensor.v1 artifact supplies graph inputs and parameters, MAC operations run
    // on the mapped datapath (weights copied into the transitional nebula layer, re-laid
    // reduction-tile-major when the mapping folds INPUT_CHANNEL), convolutions with
    // P/Q > 1 use the im2col kernel, and non-MAC operations run dedicated value kernels.
    // Every operation with an artifact golden is compared; graph outputs are mandatory.
    functional_artifact_t functional_artifact;
    std::map<std::string, std::vector<float>> functional_tensor_store;
    std::vector<float> &functional_store(const std::string &m_tensor_id);
    void functional_bind_executable_operation(unsigned m_index,
                                              const workload_operation_t &m_operation);
    void functional_commit_executable_operation(unsigned m_index,
                                                const workload_operation_t &m_operation);
    void functional_execute_graph_operation(unsigned m_index,
                                            const workload_operation_t &m_operation);
    // B-4: batched matmul of two activations (attention Q*K^T / score*V, general bmm).
    // A mapped MAC op whose per-batch operands vary, so the single-tile datapath cannot
    // compute it -- this kernel reads A and B from the tensor store and writes the batched
    // result; timing still comes from the mapped datapath run.
    void functional_matmul(unsigned m_index, const workload_operation_t &m_operation);
public:
    // Path to the npusim.tensor.v1 manifest (set by main.cc for run-ir-functional before
    // init; empty = timing-only executable run, which a FUNCTIONAL build refuses).
    std::string functional_artifact_path;
protected:
    // Write functional_report to result/functional/<network>/report.json.
    void write_functional_report(const std::string &m_network_label) const;
public:
    // Acceptance gate (plan §8 "Failure semantics"): the run FAILS if any verified layer
    // mismatched, or if an external golden was supplied but nothing was ever compared (a
    // silently-uncompared or vacuous run must not read as success). main() maps this to a
    // non-zero process exit.
    bool functional_failed() const {
        return functional_layers_failed > 0 ||
               (functional_external_golden && functional_layers_checked == 0);
    }
protected:
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
