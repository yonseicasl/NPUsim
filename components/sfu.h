#ifndef __SFU_H__
#define __SFU_H__

#include <cstddef>
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
    /* softmax micro-operations */
    SFU_OP_EXP,
    SFU_OP_RECIP,
    SFU_OP_VMAX,         // vector max-reduction step
    SFU_OP_VADD,         // vector add/subtract (subtract pass and sum-reduction)
    SFU_OP_VMUL,         // vector multiply (normalize pass)
    NUM_SFU_OPS,
};

// Per-operation timing/energy profile, parsed from the [sfu] section as
//   <name>_pipeline_latency / <name>_initiation_interval / sfu_op_energy_<name>.
// The declared flags keep "defaulted 1/1, unpriced" distinguishable from a calibrated
// profile: an op that fires on a defaulted profile is reported UNCALIBRATED, and an op
// that fires without a declared energy key becomes an unpriced active event.
struct sfu_op_profile_t {
    const char *name;
    bool bypass;               // linear: modeled as free passthrough
    bool micro_op;             // softmax building block, not a layer-level activation
    double pipeline_latency;
    double initiation_interval;
    double op_energy;
    bool latency_declared;
    bool ii_declared;
    bool energy_declared;
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
    sfu_invocation_t elementwise_invocation(sfu_op_t m_op, size_t m_valid_elements) const;
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
    unsigned queue_depth;           // producer->SFU tile queue; 1 = serial contract
    double setup_cycle;             // per-invocation setup cycles
    bool strict_profiles;           // fail instead of defaulting an undeclared profile
    std::string placement;          // post_accumulator only (initial contract)
    std::string fusion;             // fused only (initial contract)

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
