#ifndef __SFU_H__
#define __SFU_H__

#include <cstddef>
#include <set>
#include <string>
#include <vector>

#include "config.h"

// Special Function Unit (plan/plan_sfu.md).
//
// The SFU is a per-chip post-processing component placed between the final accumulator
// value and the output-format cast:
//
//     final accumulator -> [SFU activation] -> output cast -> writeback
//
// It models the hardware COST (cycles, internal traffic, dynamic/static energy) of the
// activation applied to each valid output element exactly once, after all C/R/S reduction
// has completed. It never mutates functional data on the main path: Nebula's forward()
// remains the single functional owner of activation values (sfu_t only carries a pure
// evaluator for unit tests).
//
// Element-wise operations use the plan's invocation formula
//     chunks     = ceil(valid_elements / lanes)
//     SFU cycles = setup + pipeline_latency + (chunks - 1) * initiation_interval
// and softmax is a multi-pass microprogram (max -> subtract -> exp -> sum -> reciprocal
// -> normalize) built from the micro-operation profiles below.

// SFU operations. The first block maps 1:1 onto Nebula's element-wise activations; GELU
// and SiLU are SFU-level operations for LLM workloads (Nebula's frontend has no enum for
// them yet, so they are reachable through the SFU API and priced by config, but cannot be
// selected by a network config until the frontend grows the type). The last block is the
// softmax micro-operation set.
enum sfu_op_t {
    SFU_OP_LINEAR = 0,   // bypass: no dynamic event (elements/cycles/op energy); the
                         // physical SFU still accrues always-on leakage over the layer
                         // window, since no power gating is modeled (plan energy model)
    SFU_OP_RELU,
    SFU_OP_LEAKY,        // slope 0.1 (Nebula/Darknet convention)
    SFU_OP_ELU,
    SFU_OP_SIGMOID,      // Nebula names this "logistic"; both map here
    SFU_OP_TANH,
    SFU_OP_HSIGMOID,
    SFU_OP_HSWISH,
    SFU_OP_GELU,         // tanh approximation
    SFU_OP_SILU,         // x * sigmoid(x) (a.k.a. swish)
    SFU_OP_LOGGY,        // 2*sigmoid(x) - 1 (Darknet/Nebula loggy)
    /* softmax micro-operations */
    SFU_OP_EXP,
    SFU_OP_RECIP,
    SFU_OP_VMAX,         // vector max-reduction step
    SFU_OP_VADD,         // vector add/subtract (subtract pass and sum-reduction)
    SFU_OP_VMUL,         // vector multiply (normalize pass)
    NUM_SFU_OPS,
};

// Per-operation timing/energy profile, parsed from the [sfu] section as
//   <name>_pipeline_latency / <name>_initiation_interval / sfu_op_energy_<name>
//   (+ Phase-6: <name>_approximation).
// The declared flags keep "defaulted 1/1, unpriced" distinguishable from a calibrated
// profile: an op that fires on a defaulted profile is reported UNCALIBRATED, and an op
// that fires without a declared energy key becomes an unpriced active event.
struct sfu_op_profile_t {
    const char *name;
    bool bypass;               // linear: modeled as free passthrough
    bool micro_op;             // softmax building block, not a layer-level activation
    // Phase-6: is this a transcendental operation (LUT/polynomial hardware) or a
    // piecewise-linear ALU operation? Decides which approximation modes are legal.
    bool transcendental;
    double pipeline_latency;
    double initiation_interval;
    double op_energy;
    bool latency_declared;
    bool ii_declared;
    bool energy_declared;
    // Phase-6: declared implementation the profile describes. ALU-class operations
    // accept {exact, piecewise}; transcendental ones accept {lut, polynomial,
    // piecewise}. Anything else -- and any declaration on the linear bypass -- is a
    // config error. Empty = not declared (allowed outside strict_profiles).
    std::string approximation;
};

// One SFU invocation event: cycles, operation counts, internal traffic and dynamic
// energy from the SAME event source (plan principle 6). Per-chip results combine with
// merge_parallel (latency = max, energy/traffic = sum); softmax passes combine with
// merge_serial (everything sums).
struct sfu_invocation_t {
    sfu_invocation_t();

    std::string operation;          // human-readable operation name
    size_t valid_elements;          // network-valid output elements (padding excluded)
    size_t operations;              // scalar SFU operations executed
    size_t invocations;             // pipeline setups (per unit, per pass)
    size_t chunks;                  // lane-wide issue slots
    double busy_cycle;              // SFU busy window (max over parallel units)
    double tail_lane_utilization;   // valid lanes of the final chunk / lanes

    /* internal (scratchpad/port) traffic, accumulator precision */
    size_t ingress_elements;
    size_t egress_elements;
    size_t ingress_bytes;
    size_t egress_bytes;
    size_t ingress_transactions;    // lane-wide internal transactions
    size_t egress_transactions;

    /* dynamic energy, split by event family */
    double op_energy;
    double read_energy;
    double write_energy;
    double setup_energy;

    bool timing_calibrated;         // every op used had a declared latency AND II
    std::vector<std::string> unpriced;  // events that fired with no declared cost key

    // Chips run in parallel: busy = max, counts/energy/traffic = sum.
    void merge_parallel(const sfu_invocation_t &m_other);
    // Softmax passes run serially: everything sums.
    void merge_serial(const sfu_invocation_t &m_other);
};

class sfu_t {

public:
    sfu_t(section_config_t m_section_config);
    ~sfu_t();

    // Parse and fail-fast validate the [sfu] section.
    void init(section_config_t m_section_config);
    // Print the SFU specification (units, lanes, per-op profiles).
    void print_specification();
    // Per-layer reset (the SFU keeps no cross-layer state; kept for component symmetry).
    void reset();

    // Map a Nebula activation name (nebula::activation_type_str) to an SFU operation.
    // Returns false for activations the SFU does not support -- the caller must fail
    // fast rather than fall back to ReLU (plan: no silent fallback).
    static bool op_from_name(const std::string &m_name, sfu_op_t *m_op);
    static const char *op_name(sfu_op_t m_op);

    // Element-wise activation over m_valid_elements final output elements on this chip.
    // Phase-2: m_commit_events is the number of final_output_tile events the elements
    // arrive in (from the multi-chip output-commit boundary). Each commit is one SFU
    // invocation -- its own pipeline setup and fill -- and the commits drain serially,
    // so busy = sum over events of [setup + latency + (chunks_e - 1) x II]. One event
    // (the default) is the layer-granular aggregate model.
    sfu_invocation_t elementwise_invocation(sfu_op_t m_op, size_t m_valid_elements,
                                            size_t m_commit_events = 1) const;
    // Standalone softmax: m_rows independent softmax vectors of length m_row_length,
    // decomposed into max -> subtract -> exp -> sum -> reciprocal -> normalize.
    sfu_invocation_t softmax_invocation(size_t m_rows, size_t m_row_length) const;

    // Pure functional evaluators for unit tests ONLY. The main FUNCTIONAL path must not
    // call these: Nebula's forward() is the single functional owner, and a second
    // application would double-transform the output (plan, functional-model policy).
    static float evaluate(sfu_op_t m_op, float m_x);
    static void evaluate_softmax(const float *m_input, float *m_output, size_t m_size);

    unsigned get_num_units() const { return num_units; }
    unsigned get_lanes() const { return lanes; }
    unsigned get_queue_depth() const { return queue_depth; }
    double get_setup_cycle() const { return setup_cycle; }
    double get_static_energy_per_cycle() const { return u_static_energy; }
    bool static_energy_declared() const { return static_energy_is_declared; }
    // Phase-7: where the standalone-softmax operand tensor lives -- "dram" (materialized
    // round trip through the memory hierarchy; matches the simulator's layer flow, which
    // commits every layer's output off-chip) or "glb" (retained on-chip by a fused
    // schedule; requires the tensor to fit).
    const std::string &get_softmax_operand_residency() const { return softmax_operand_residency; }
    // Phase-8: where the timing profile's numbers came from (RTL commit, trace set,
    // tool versions). Empty = not declared = the profile is not calibration-grade.
    const std::string &get_profile_reference() const { return profile_reference; }
    // Phase-6: precision qualification. `profile_precision` declares the operand format
    // the timing profile was characterized for; a mismatch with the runtime accumulator
    // format (the precision the SFU actually processes) fails fast under
    // strict_profiles and otherwise marks every invocation UNCALIBRATED, with the note
    // below carried into the layer report.
    const std::string &get_precision_note() const { return precision_note; }
    // Phase-6/8: architecture allowlist. When `supported_ops` is declared, operations
    // OUTSIDE the list do not exist in this SFU's hardware (e.g. nv_small's SDP has no
    // LUT/EW engine at all -- spec SDP_LUT_DISABLE / SDP_EW_DISABLE), and any use fails
    // fast: a missing execution unit is an architecture fact, not an uncalibrated cost.
    // Undeclared = every modeled operation is available.
    bool op_supported(sfu_op_t m_op) const;
    bool has_supported_ops_contract() const { return !supported_ops.empty(); }
    const sfu_op_profile_t &profile(sfu_op_t m_op) const { return profiles[m_op]; }

    unsigned index;

private:
    // Cycle cost of one continuous invocation of `chunks` lane-wide chunks.
    double invocation_cycles(const sfu_op_profile_t &m_profile, size_t m_chunks) const;
    // Strict-profile gate: abort when a used operation has no declared profile.
    void require_profile(const sfu_op_profile_t &m_profile) const;
    // Track defaulted timing and undeclared energy for a fired operation.
    void note_profile_use(const sfu_op_profile_t &m_profile, size_t m_operations,
                          sfu_invocation_t *m_invocation) const;
    // Add read/write/setup traffic-and-energy for one pass.
    void charge_traffic(size_t m_reads, size_t m_writes, size_t m_setups,
                        sfu_invocation_t *m_invocation) const;

    unsigned num_units;             // parallel SFU pipelines per chip
    unsigned lanes;                 // elements one pipeline accepts per initiation
    // Producer->SFU output-tile queue depth (Phase 5). 1 = the serial contract (the SFU
    // window is appended to the critical path); >= 2 = the SFU joins the per-tile
    // pipeline timeline as a post-processing stage with this staging depth, so its work
    // overlaps the producer up to the queue's back-pressure limit.
    unsigned queue_depth;
    double setup_cycle;             // per-invocation setup cycles
    bool strict_profiles;           // fail instead of defaulting an undeclared profile
    std::string placement;          // post_accumulator only (initial contract)
    std::string fusion;             // fused only (initial contract)
    std::string softmax_operand_residency;   // "dram" (default) | "glb"
    std::string profile_reference;           // Phase-8 provenance (free text)
    std::string profile_precision;           // Phase-6: format the profile describes
    bool precision_mismatch;                 // declared precision != runtime accumulator
    std::string precision_note;              // report text for a mismatch
    std::set<std::string> supported_ops;     // architecture allowlist (empty = all)

    double u_read_energy;           // per ingress element
    double u_write_energy;          // per egress element
    double u_setup_energy;          // per invocation
    double u_static_energy;         // pJ/cycle per physical SFU unit
    bool read_energy_is_declared;
    bool write_energy_is_declared;
    bool setup_energy_is_declared;
    bool static_energy_is_declared;

    sfu_op_profile_t profiles[NUM_SFU_OPS];
};

#endif
