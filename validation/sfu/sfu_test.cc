// SFU validation (plan/plan_sfu.md, test strategy).
//
// Formula/boundary, event-count, energy-linearity, parallelism, softmax multi-pass and
// functional-evaluator checks for components/sfu.{h,cc}. Invalid-config fail-fast cases
// exit(1) inside sfu_t, so they run as separate subprocess modes driven by
// run_sfu_validation.sh (see `invalid-*` modes at the bottom of main()).

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "config.h"
#include "sfu.h"
// Phase-6 item 3: the SFU's pure evaluators are compared against the ACTUAL Nebula
// activation implementations (ext/nebula/common/activations.cc is linked into this
// test), not against re-derived formulas.
#include "activations.h"

namespace {

void fail(const std::string &message) {
    std::cerr << "sfu validation failed: " << message << std::endl;
    std::exit(1);
}

void expect(bool condition, const std::string &message) {
    if(!condition) fail(message);
}

void expect_near(double value, double expected, const std::string &message) {
    const double tolerance = 1e-9*std::max(1.0, std::fabs(expected));
    if(std::fabs(value - expected) > tolerance) {
        fail(message + " (got " + std::to_string(value) + ", expected " +
             std::to_string(expected) + ")");
    }
}

// Functional comparisons: the evaluator computes in float, the reference in double, so
// the agreement bound is single-precision rounding, not exactness.
void expect_close(double value, double expected, const std::string &message) {
    const double tolerance = 1e-6*std::max(1.0, std::fabs(expected));
    if(std::fabs(value - expected) > tolerance) {
        fail(message + " (got " + std::to_string(value) + ", expected " +
             std::to_string(expected) + ")");
    }
}

// The gemmini_sfu.cfg profile, built programmatically so every expectation below is a
// hand calculation against known unit costs.
section_config_t reference_section() {
    section_config_t section("sfu");
    section.add_setting("num_units_per_chip", "1");
    section.add_setting("lanes", "16");
    section.add_setting("queue_depth", "1");
    section.add_setting("setup_cycle", "0");
    section.add_setting("relu_pipeline_latency", "1");
    section.add_setting("relu_initiation_interval", "1");
    section.add_setting("leaky_pipeline_latency", "2");
    section.add_setting("leaky_initiation_interval", "1");
    section.add_setting("sigmoid_pipeline_latency", "6");
    section.add_setting("sigmoid_initiation_interval", "1");
    section.add_setting("exp_pipeline_latency", "6");
    section.add_setting("exp_initiation_interval", "1");
    section.add_setting("recip_pipeline_latency", "8");
    section.add_setting("recip_initiation_interval", "2");
    section.add_setting("vmax_pipeline_latency", "1");
    section.add_setting("vmax_initiation_interval", "1");
    section.add_setting("vadd_pipeline_latency", "1");
    section.add_setting("vadd_initiation_interval", "1");
    section.add_setting("vmul_pipeline_latency", "2");
    section.add_setting("vmul_initiation_interval", "1");
    section.add_setting("sfu_op_energy_relu", "0.10");
    section.add_setting("sfu_op_energy_leaky", "0.18");
    section.add_setting("sfu_op_energy_sigmoid", "0.60");
    section.add_setting("sfu_op_energy_exp", "0.60");
    section.add_setting("sfu_op_energy_recip", "0.70");
    section.add_setting("sfu_op_energy_vmax", "0.08");
    section.add_setting("sfu_op_energy_vadd", "0.08");
    section.add_setting("sfu_op_energy_vmul", "0.15");
    section.add_setting("sfu_read_energy", "0.02");
    section.add_setting("sfu_write_energy", "0.02");
    section.add_setting("sfu_setup_energy", "0.00");
    section.add_setting("sfu_static_energy", "0.001");
    return section;
}

// Bypass and zero-element boundary: no element, cycle, energy or event (plan table).
void test_bypass_and_zero() {
    sfu_t sfu(reference_section());
    const sfu_invocation_t linear = sfu.elementwise_invocation(SFU_OP_LINEAR, 1000000);
    expect(linear.valid_elements == 0 && linear.busy_cycle == 0.0 &&
           linear.operations == 0 && linear.invocations == 0 &&
           linear.op_energy == 0.0 && linear.read_energy == 0.0 &&
           linear.ingress_bytes == 0 && linear.unpriced.empty(),
           "linear bypass must charge nothing");
    const sfu_invocation_t empty = sfu.elementwise_invocation(SFU_OP_RELU, 0);
    expect(empty.valid_elements == 0 && empty.busy_cycle == 0.0 &&
           empty.operations == 0 && empty.invocations == 0 && empty.op_energy == 0.0,
           "zero valid elements must charge nothing");
}

// Lane-tail boundary: N < lanes, N = lanes, N = lanes + 1 (plan table).
void test_lane_tail() {
    sfu_t sfu(reference_section());
    const sfu_invocation_t under = sfu.elementwise_invocation(SFU_OP_RELU, 15);
    expect(under.chunks == 1, "15 elements on 16 lanes is one chunk");
    expect_near(under.busy_cycle, 1.0, "one-chunk relu is pipeline latency alone");
    expect_near(under.tail_lane_utilization, 15.0/16.0, "15/16 tail utilization");
    const sfu_invocation_t exact = sfu.elementwise_invocation(SFU_OP_RELU, 16);
    expect(exact.chunks == 1, "16 elements on 16 lanes is one chunk");
    expect_near(exact.busy_cycle, 1.0, "exact-chunk relu is pipeline latency alone");
    expect_near(exact.tail_lane_utilization, 1.0, "full tail utilization");
    const sfu_invocation_t over = sfu.elementwise_invocation(SFU_OP_RELU, 17);
    expect(over.chunks == 2, "17 elements on 16 lanes is two chunks");
    expect_near(over.busy_cycle, 2.0, "two-chunk relu is latency + one II");
    expect_near(over.tail_lane_utilization, 1.0/16.0, "1/16 tail utilization");
}

// The plan's invocation formula and per-operation profile independence.
void test_formula_and_profiles() {
    sfu_t sfu(reference_section());
    // 65536 elements / 16 lanes = 4096 chunks: relu 1 + 4095, leaky 2 + 4095.
    const sfu_invocation_t relu = sfu.elementwise_invocation(SFU_OP_RELU, 65536);
    expect_near(relu.busy_cycle, 4096.0, "relu 65536-element window");
    expect(relu.operations == 65536 && relu.invocations == 1 && relu.chunks == 4096,
           "relu event counts");
    expect_near(relu.op_energy, 6553.6, "relu operation energy 65536 x 0.10");
    expect_near(relu.read_energy, 1310.72, "relu ingress energy 65536 x 0.02");
    expect_near(relu.write_energy, 1310.72, "relu egress energy 65536 x 0.02");
    const sfu_invocation_t leaky = sfu.elementwise_invocation(SFU_OP_LEAKY, 65536);
    expect_near(leaky.busy_cycle, 4097.0, "leaky 65536-element window");
    expect_near(leaky.op_energy, 65536*0.18, "leaky operation energy is independent of relu");
    expect(relu.timing_calibrated && leaky.timing_calibrated, "declared profiles are calibrated");
    expect(relu.unpriced.empty() && leaky.unpriced.empty(),
           "fully priced profile has no unpriced events");
}

// Phase-2: commit-event granular invocations. Each final_output_tile event is one SFU
// invocation with its own setup and pipeline fill; the events drain serially.
void test_commit_event_granularity() {
    sfu_t sfu(reference_section());
    // relu (L = II = 1, setup 0): per-event fill adds L - II = 0, so the window is
    // event-count invariant -- 65536 elements in 256 commits must still cost 4096.
    const sfu_invocation_t relu_single = sfu.elementwise_invocation(SFU_OP_RELU, 65536, 1);
    const sfu_invocation_t relu_tiled = sfu.elementwise_invocation(SFU_OP_RELU, 65536, 256);
    expect_near(relu_single.busy_cycle, 4096.0, "single-event relu window");
    expect_near(relu_tiled.busy_cycle, 4096.0,
                "L == II relu window is commit-count invariant");
    expect(relu_tiled.invocations == 256, "one invocation per commit event");
    expect(relu_tiled.chunks == relu_single.chunks,
           "lane-aligned commits conserve chunks");
    expect(relu_tiled.operations == relu_single.operations &&
           relu_tiled.op_energy == relu_single.op_energy,
           "commit granularity never changes operation count or energy");
    // leaky (L = 2, II = 1): every commit pays its own fill, so 64 elements in 4
    // commits of 16 (1 chunk each) cost 4 x 2 = 8 vs the single-event 2 + 3 = 5.
    const sfu_invocation_t leaky_single = sfu.elementwise_invocation(SFU_OP_LEAKY, 64, 1);
    const sfu_invocation_t leaky_tiled = sfu.elementwise_invocation(SFU_OP_LEAKY, 64, 4);
    expect_near(leaky_single.busy_cycle, 5.0, "single-event leaky window");
    expect_near(leaky_tiled.busy_cycle, 8.0, "per-commit pipeline fill accumulates");
    // Uneven split: 17 elements in 2 commits (9 + 8) on 16 lanes -- per-event tail
    // rounding, worst-event tail utilization reported.
    const sfu_invocation_t uneven = sfu.elementwise_invocation(SFU_OP_RELU, 17, 2);
    expect(uneven.chunks == 2 && uneven.ingress_transactions == 2,
           "per-event tails round chunks and transactions up");
    expect_near(uneven.tail_lane_utilization, 8.0/16.0,
                "worst commit's tail utilization is reported");
    // Boundary: more commits than elements clamps to one element per commit.
    const sfu_invocation_t clamped = sfu.elementwise_invocation(SFU_OP_RELU, 2, 10);
    expect(clamped.invocations == 2, "commit count clamps to the element count");
    // Boundary: zero commits behaves as the single aggregate invocation.
    const sfu_invocation_t zero_events = sfu.elementwise_invocation(SFU_OP_RELU, 32, 0);
    expect(zero_events.invocations == 1 && zero_events.busy_cycle == 2.0,
           "zero commit events falls back to one invocation");
}

// Energy perturbation: doubling one operation cost doubles ONLY that operation's energy.
void test_energy_linearity() {
    section_config_t doubled = reference_section();
    // add_setting refuses duplicates; build the x2 fixture from scratch instead.
    section_config_t section("sfu");
    section.add_setting("lanes", "16");
    section.add_setting("relu_pipeline_latency", "1");
    section.add_setting("relu_initiation_interval", "1");
    section.add_setting("sfu_op_energy_relu", "0.20");
    section.add_setting("sfu_read_energy", "0.02");
    section.add_setting("sfu_write_energy", "0.02");
    section.add_setting("sfu_setup_energy", "0.00");
    sfu_t base(reference_section());
    sfu_t perturbed(section);
    const sfu_invocation_t before = base.elementwise_invocation(SFU_OP_RELU, 4096);
    const sfu_invocation_t after = perturbed.elementwise_invocation(SFU_OP_RELU, 4096);
    expect_near(after.op_energy, 2.0*before.op_energy, "2x op energy doubles op energy");
    expect_near(after.read_energy, before.read_energy, "read energy is untouched");
    expect_near(after.write_energy, before.write_energy, "write energy is untouched");
    expect_near(after.busy_cycle, before.busy_cycle, "energy perturbation leaves cycles alone");
    (void)doubled;
}

// Intra-chip unit parallelism: latency follows the busiest unit, work sums.
void test_unit_parallelism() {
    section_config_t section("sfu");
    section.add_setting("num_units_per_chip", "2");
    section.add_setting("lanes", "16");
    section.add_setting("relu_pipeline_latency", "1");
    section.add_setting("relu_initiation_interval", "1");
    section.add_setting("sfu_op_energy_relu", "0.10");
    section.add_setting("sfu_read_energy", "0.02");
    section.add_setting("sfu_write_energy", "0.02");
    section.add_setting("sfu_setup_energy", "0.00");
    sfu_t two_units(section);
    sfu_t one_unit(reference_section());
    const sfu_invocation_t split = two_units.elementwise_invocation(SFU_OP_RELU, 65536);
    const sfu_invocation_t serial = one_unit.elementwise_invocation(SFU_OP_RELU, 65536);
    expect_near(split.busy_cycle, 2048.0, "two units halve the 4096-chunk window");
    expect(split.chunks == serial.chunks, "unit split conserves chunks");
    expect(split.invocations == 2, "one setup per active unit");
    expect_near(split.op_energy, serial.op_energy, "unit split conserves operation energy");
    expect_near(split.read_energy, serial.read_energy, "unit split conserves ingress energy");
    // Tail utilization follows the busiest UNIT's final chunk, not the whole-share
    // remainder: 17 elements on 2 units split 9/8, so the tail is 9/16 -- the pre-fix
    // whole-share formula reported 17 % 16 = 1/16.
    const sfu_invocation_t uneven = two_units.elementwise_invocation(SFU_OP_RELU, 17);
    expect_near(uneven.tail_lane_utilization, 9.0/16.0,
                "multi-unit tail utilization is the busiest unit's final chunk");
}

// Chip-parallel merge: element shares sum EXACTLY to B x K x P x Q, latency is the max.
void test_chip_partition_identity() {
    sfu_t sfu(reference_section());
    const size_t valid_elements = 4*96*55*55;   // AlexNet conv0: B x K x P x Q
    const unsigned chips = 4;
    sfu_invocation_t combined;
    const size_t base = valid_elements/chips;
    const size_t remainder = valid_elements % chips;
    double max_busy = 0.0;
    for(unsigned c = 0; c < chips; ++c) {
        const size_t share = base + (c < remainder ? 1 : 0);
        const sfu_invocation_t inv = sfu.elementwise_invocation(SFU_OP_RELU, share);
        max_busy = std::max(max_busy, inv.busy_cycle);
        combined.merge_parallel(inv);
    }
    expect(combined.valid_elements == valid_elements,
           "chip partition element shares must sum to the layer's valid output count");
    expect_near(combined.busy_cycle, max_busy, "chip-parallel latency is the max");
    expect(combined.operations == valid_elements, "one operation per valid element");
}

// Softmax multi-pass hand calculation (rows=256, N=256, lanes=16, profiles above):
//   max   : 1 + (256*16 - 1)*1 + 256*ceil(log2(16))*1 = 1 + 4095 + 1024 = 5120
//   sub   : 1 + (4096 - 1)*1                           = 4096
//   exp   : 6 + 4095                                   = 4101
//   sum   : 1 + 4095 + 1024                            = 5120
//   recip : 8 + (16 - 1)*2                             = 38
//   mul   : 2 + 4095                                   = 4097
//   total                                              = 22572
void test_softmax_hand_calculation() {
    sfu_t sfu(reference_section());
    const sfu_invocation_t softmax = sfu.softmax_invocation(256, 256);
    expect_near(softmax.busy_cycle, 22572.0, "softmax multi-pass window");
    expect(softmax.valid_elements == 65536, "softmax valid elements = rows x length");
    // ops: max 256*255 + sub 65536 + exp 65536 + sum 256*255 + recip 256 + mul 65536
    expect(softmax.operations == 327424, "softmax scalar operation count");
    // reads: 65536 + (65536+256) + 65536 + 65536 + 256 + (65536+256) = 328448
    expect(softmax.ingress_elements == 328448, "softmax scratchpad reads");
    // writes: 256 + 65536 + 65536 + 256 + 256 + 65536 = 197376
    expect(softmax.egress_elements == 197376, "softmax scratchpad writes");
    expect(softmax.invocations == 6, "six passes, one unit");
    const double expected_op_energy = 65280*0.08 + 65536*0.08 + 65536*0.60 +
                                      65280*0.08 + 256*0.70 + 65536*0.15;
    expect_near(softmax.op_energy, expected_op_energy, "softmax per-pass operation energy");
    expect(softmax.timing_calibrated && softmax.unpriced.empty(),
           "fully profiled softmax is calibrated and priced");

    // Boundary: empty tensors charge nothing; a single-element row skips the reduction
    // passes ENTIRELY -- a length-1 max/sum performs zero operations, so it must charge
    // zero cycles too (every charged cycle needs an event source), and its possibly
    // undeclared vmax profile must not silently produce cycles.
    const sfu_invocation_t empty = sfu.softmax_invocation(0, 256);
    expect(empty.busy_cycle == 0.0 && empty.operations == 0, "empty softmax charges nothing");
    const sfu_invocation_t single = sfu.softmax_invocation(4, 1);
    expect(single.operations == 4*0 + 4 + 4 + 4*0 + 4 + 4,
           "length-1 softmax has no reduction operations");
    // sub 1 + exp 6 + recip 8 + mul 2 = 17; no vmax/vadd reduction cycles.
    expect_near(single.busy_cycle, 17.0, "length-1 softmax runs only the 4 element-wise passes");
    expect(single.invocations == 4, "length-1 softmax sets up 4 passes, not 6");
}

// Reduction-tree fold steps are serially data-dependent, so each step costs the full
// PIPELINE LATENCY, not one initiation interval. rows=1, n=16, lanes=16 with
// vmax latency=4/II=1 and single-cycle other ops:
//   max   : (4 + (1-1)*1) + 1 row * 4 fold steps * 4 latency = 20
//   sub   : 1, exp : 1
//   sum   : (1 + 0) + 1*4*1 = 5
//   recip : 1, mul : 1        -> total 29
void test_softmax_fold_latency() {
    section_config_t section("sfu");
    section.add_setting("lanes", "16");
    section.add_setting("vmax_pipeline_latency", "4");
    section.add_setting("vmax_initiation_interval", "1");
    section.add_setting("vadd_pipeline_latency", "1");
    section.add_setting("vadd_initiation_interval", "1");
    section.add_setting("exp_pipeline_latency", "1");
    section.add_setting("exp_initiation_interval", "1");
    section.add_setting("recip_pipeline_latency", "1");
    section.add_setting("recip_initiation_interval", "1");
    section.add_setting("vmul_pipeline_latency", "1");
    section.add_setting("vmul_initiation_interval", "1");
    section.add_setting("sfu_op_energy_vmax", "0");
    section.add_setting("sfu_op_energy_vadd", "0");
    section.add_setting("sfu_op_energy_exp", "0");
    section.add_setting("sfu_op_energy_recip", "0");
    section.add_setting("sfu_op_energy_vmul", "0");
    section.add_setting("sfu_read_energy", "0");
    section.add_setting("sfu_write_energy", "0");
    section.add_setting("sfu_setup_energy", "0");
    sfu_t sfu(section);
    const sfu_invocation_t softmax = sfu.softmax_invocation(1, 16);
    expect_near(softmax.busy_cycle, 29.0,
                "fold steps cost pipeline latency per serial dependency step");
}

// Serial pass merge: two invocations back to back add their windows.
void test_serial_merge() {
    sfu_t sfu(reference_section());
    sfu_invocation_t first = sfu.elementwise_invocation(SFU_OP_RELU, 65536);
    const sfu_invocation_t second = sfu.elementwise_invocation(SFU_OP_LEAKY, 65536);
    const double sum = first.busy_cycle + second.busy_cycle;
    first.merge_serial(second);
    expect_near(first.busy_cycle, sum, "serial merge adds windows");
    expect(first.valid_elements == 131072, "serial merge sums elements");
}

// Unpriced/uncalibrated detection: an active op with no declared cost or timing must be
// reported, never silently zero (plan principle 7).
void test_unpriced_and_uncalibrated() {
    section_config_t section("sfu");
    section.add_setting("lanes", "16");
    // relu has NO declared latency/II/energy, and no read/write/setup energy either.
    sfu_t sfu(section);
    const sfu_invocation_t inv = sfu.elementwise_invocation(SFU_OP_RELU, 64);
    expect(!inv.timing_calibrated, "defaulted 1/1 profile must be flagged uncalibrated");
    bool op_unpriced = false, read_unpriced = false, write_unpriced = false,
         setup_unpriced = false;
    for(size_t i = 0; i < inv.unpriced.size(); ++i) {
        if(inv.unpriced[i].find("sfu_op_energy_relu") != std::string::npos) op_unpriced = true;
        if(inv.unpriced[i].find("sfu_read_energy") != std::string::npos) read_unpriced = true;
        if(inv.unpriced[i].find("sfu_write_energy") != std::string::npos) write_unpriced = true;
        if(inv.unpriced[i].find("sfu_setup_energy") != std::string::npos) setup_unpriced = true;
    }
    expect(op_unpriced && read_unpriced && write_unpriced && setup_unpriced,
           "every active-but-undeclared cost key must surface as an unpriced event");
    // The defaulted timing still uses the 1/1 formula, so cycles stay well defined.
    expect_near(inv.busy_cycle, 4.0, "64 elements / 16 lanes on the defaulted 1/1 profile");
}

// Nebula activation-name mapping: aliases resolve, unsupported names refuse (no silent
// ReLU fallback).
void test_activation_mapping() {
    sfu_op_t op;
    expect(sfu_t::op_from_name("relu", &op) && op == SFU_OP_RELU, "relu maps");
    expect(sfu_t::op_from_name("logistic", &op) && op == SFU_OP_SIGMOID,
           "Nebula's 'logistic' maps to sigmoid");
    expect(sfu_t::op_from_name("swish", &op) && op == SFU_OP_SILU, "swish maps to silu");
    expect(sfu_t::op_from_name("loggy", &op) && op == SFU_OP_LOGGY, "loggy maps");
    expect(sfu_t::op_from_name("gelu", &op) && op == SFU_OP_GELU, "gelu maps");
    expect(!sfu_t::op_from_name("hardtan", &op), "unsupported activations must refuse");
    expect(!sfu_t::op_from_name("stair", &op), "unsupported activations must refuse");
    expect(!sfu_t::op_from_name("undefined_activation", &op), "undefined must refuse");
}

// Phase-6 item 3: evaluator vs the LINKED Nebula implementation, element for element.
// Covers every operation both sides implement -- including the Phase-6 additions
// (gelu/silu in Nebula's frontend, loggy in the SFU). Nebula's elu is deliberately
// excluded: its implementation returns x*x for x > 0, which is not ELU.
void test_against_nebula_reference() {
    const float inputs[] = {-6.5f, -3.0f, -1.25f, -0.5f, 0.0f, 0.5f, 1.25f, 3.0f, 6.5f};
    const size_t count = sizeof(inputs)/sizeof(inputs[0]);
    struct case_t {
        sfu_op_t op;
        void (*reference)(float*, unsigned);
        const char *label;
    };
    const case_t cases[] = {
        { SFU_OP_RELU,     nebula::relu_activation,     "relu vs Nebula" },
        { SFU_OP_LEAKY,    nebula::leaky_activation,    "leaky vs Nebula" },
        { SFU_OP_SIGMOID,  nebula::logistic_activation, "sigmoid vs Nebula logistic" },
        { SFU_OP_SIGMOID,  nebula::sigmoid_activation,  "sigmoid vs Nebula sigmoid" },
        { SFU_OP_TANH,     nebula::tanh_activation,     "tanh vs Nebula" },
        { SFU_OP_HSIGMOID, nebula::hsigmoid_activation, "hsigmoid vs Nebula" },
        { SFU_OP_HSWISH,   nebula::hswish_activation,   "hswish vs Nebula" },
        { SFU_OP_LOGGY,    nebula::loggy_activation,    "loggy vs Nebula" },
        { SFU_OP_GELU,     nebula::gelu_activation,     "gelu vs Nebula" },
        { SFU_OP_SILU,     nebula::silu_activation,     "silu vs Nebula" },
    };
    for(unsigned c = 0; c < sizeof(cases)/sizeof(cases[0]); ++c) {
        float reference[sizeof(inputs)/sizeof(inputs[0])];
        for(size_t i = 0; i < count; ++i) reference[i] = inputs[i];
        cases[c].reference(reference, static_cast<unsigned>(count));
        for(size_t i = 0; i < count; ++i) {
            expect_close(sfu_t::evaluate(cases[c].op, inputs[i]), reference[i],
                         std::string(cases[c].label) + " at x = " +
                         std::to_string(inputs[i]));
        }
    }
}

// Phase-6/8: architecture allowlist. Declared supported_ops make out-of-list operations
// non-existent (fail-fast, tested via subprocess modes); linear is always implied.
void test_supported_ops_contract() {
    section_config_t section = reference_section();
    section.add_setting("supported_ops", "relu,leaky");
    sfu_t sfu(section);
    expect(sfu.op_supported(SFU_OP_RELU) && sfu.op_supported(SFU_OP_LEAKY),
           "listed operations are supported");
    expect(sfu.op_supported(SFU_OP_LINEAR), "the linear bypass is always implied");
    expect(!sfu.op_supported(SFU_OP_SIGMOID) && !sfu.op_supported(SFU_OP_EXP),
           "unlisted operations are architecturally absent");
    sfu_t open_sfu(reference_section());
    expect(open_sfu.op_supported(SFU_OP_EXP),
           "no declared allowlist means every modeled operation is available");
}

// Phase-6: precision qualification. A profile characterized for a different precision
// than the runtime accumulator (fp32 by default in this test binary) marks every
// invocation UNCALIBRATED; strict mode refuses it outright (subprocess mode below).
void test_precision_mismatch() {
    section_config_t section = reference_section();
    section.add_setting("profile_precision", "int8");
    sfu_t sfu(section);
    const sfu_invocation_t invocation = sfu.elementwise_invocation(SFU_OP_RELU, 64);
    expect(!invocation.timing_calibrated,
           "a precision-mismatched profile must be flagged UNCALIBRATED");
    expect(!sfu.get_precision_note().empty(),
           "the precision mismatch must carry a report note");
    // A matching declaration stays calibrated.
    section_config_t matching = reference_section();
    matching.add_setting("profile_precision", "fp32");
    sfu_t matched(matching);
    expect(matched.elementwise_invocation(SFU_OP_RELU, 64).timing_calibrated,
           "a precision-matched profile stays calibrated");
}

// Functional evaluator against reference math (unit-test-only path; Nebula semantics:
// leaky slope 0.1, logistic sigmoid, relu6-based hard ops).
void test_functional_evaluators() {
    const float inputs[] = {-6.5f, -3.0f, -1.25f, -0.5f, 0.0f, 0.5f, 1.25f, 3.0f, 6.5f};
    for(unsigned i = 0; i < sizeof(inputs)/sizeof(inputs[0]); ++i) {
        const float x = inputs[i];
        expect_close(sfu_t::evaluate(SFU_OP_RELU, x), x > 0.0f ? x : 0.0f, "relu value");
        expect_close(sfu_t::evaluate(SFU_OP_LEAKY, x), x > 0.0f ? x : 0.1f*x, "leaky value");
        expect_close(sfu_t::evaluate(SFU_OP_SIGMOID, x), 1.0/(1.0 + std::exp(-(double)x)),
                     "sigmoid value");
        expect_close(sfu_t::evaluate(SFU_OP_TANH, x), std::tanh((double)x), "tanh value");
        const double relu6 = std::min(std::max((double)x + 3.0, 0.0), 6.0);
        expect_close(sfu_t::evaluate(SFU_OP_HSIGMOID, x), relu6/6.0, "hsigmoid value");
        expect_close(sfu_t::evaluate(SFU_OP_HSWISH, x), x*relu6/6.0, "hswish value");
        expect_close(sfu_t::evaluate(SFU_OP_SILU, x), x/(1.0 + std::exp(-(double)x)),
                     "silu value");
        // GELU tanh approximation vs the exact erf definition: the approximation error
        // bound is ~1e-3 absolute.
        const double exact_gelu = 0.5*x*(1.0 + std::erf(x/std::sqrt(2.0)));
        if(std::fabs(sfu_t::evaluate(SFU_OP_GELU, x) - exact_gelu) > 5e-3) {
            fail("gelu approximation drifted from the erf definition");
        }
    }
    // Softmax: matches the naive reference and sums to 1.
    const float logits[8] = {1.0f, 2.0f, 3.0f, -1.0f, 0.5f, 4.0f, -2.5f, 0.0f};
    float out[8];
    sfu_t::evaluate_softmax(logits, out, 8);
    float max_logit = logits[0];
    for(unsigned i = 1; i < 8; ++i) max_logit = std::max(max_logit, logits[i]);
    double sum = 0.0;
    for(unsigned i = 0; i < 8; ++i) sum += std::exp((double)logits[i] - max_logit);
    double out_sum = 0.0;
    for(unsigned i = 0; i < 8; ++i) {
        const double reference = std::exp((double)logits[i] - max_logit)/sum;
        if(std::fabs(out[i] - reference) > 1e-6) fail("softmax value drifted from reference");
        out_sum += out[i];
    }
    if(std::fabs(out_sum - 1.0) > 1e-5) fail("softmax must sum to 1");
}

// Subprocess modes: each of these must exit non-zero inside sfu_t's fail-fast validation.
int run_invalid_mode(const std::string &mode) {
    section_config_t section("sfu");
    if(mode == "invalid-lanes") {
        section.add_setting("lanes", "0");
    } else if(mode == "invalid-units") {
        section.add_setting("num_units_per_chip", "0");
    } else if(mode == "invalid-queue") {
        // Phase-5: depth >= 2 is the STREAMING contract and must be accepted; only a
        // zero-depth queue (no staging register at all) is rejected.
        section.add_setting("queue_depth", "0");
    } else if(mode == "invalid-residency") {
        section.add_setting("softmax_operand_residency", "sram");
    } else if(mode == "invalid-placement") {
        section.add_setting("placement", "pre_accumulator");
    } else if(mode == "invalid-fusion") {
        section.add_setting("fusion", "materialized");
    } else if(mode == "invalid-zero-ii") {
        section.add_setting("relu_initiation_interval", "0");
    } else if(mode == "invalid-linear-profile") {
        section.add_setting("linear_pipeline_latency", "1");
    } else if(mode == "invalid-supported-ops") {
        section.add_setting("supported_ops", "relu,fancy_op");
    } else if(mode == "unsupported-op") {
        // Architecture allowlist: invoking an out-of-list op must abort.
        section.add_setting("supported_ops", "relu");
        section.add_setting("leaky_pipeline_latency", "2");
        section.add_setting("leaky_initiation_interval", "1");
        sfu_t sfu(section);
        (void)sfu.elementwise_invocation(SFU_OP_LEAKY, 64);   // must abort
        return 0;
    } else if(mode == "unsupported-softmax") {
        // Softmax needs its whole micro-op set; an allowlist without exp must abort.
        section.add_setting("supported_ops", "relu,vmax,vadd,vmul,recip");
        sfu_t sfu(section);
        (void)sfu.softmax_invocation(4, 16);   // must abort
        return 0;
    } else if(mode == "invalid-approximation") {
        // Phase-6: an ALU-class op cannot claim a LUT implementation.
        section.add_setting("relu_approximation", "lut");
    } else if(mode == "invalid-linear-approx") {
        section.add_setting("linear_approximation", "exact");
    } else if(mode == "strict-precision-mismatch") {
        section.add_setting("strict_profiles", "1");
        section.add_setting("profile_precision", "int8");   // runtime accumulator is fp32
    } else if(mode == "strict-missing-approximation") {
        // Phase-6: strict profiles require a transcendental op to state its
        // implementation; a complete latency/II/energy triple alone is half a profile.
        section.add_setting("strict_profiles", "1");
        section.add_setting("sigmoid_pipeline_latency", "6");
        section.add_setting("sigmoid_initiation_interval", "1");
        section.add_setting("sfu_op_energy_sigmoid", "0.6");
        section.add_setting("sfu_read_energy", "0.02");
        section.add_setting("sfu_write_energy", "0.02");
        section.add_setting("sfu_setup_energy", "0");
        sfu_t sfu(section);
        (void)sfu.elementwise_invocation(SFU_OP_SIGMOID, 64);   // must abort
        return 0;
    } else if(mode == "strict-missing-profile") {
        section.add_setting("strict_profiles", "1");
        sfu_t sfu(section);
        (void)sfu.elementwise_invocation(SFU_OP_RELU, 64);   // must abort: no relu profile
        return 0;
    } else {
        std::cerr << "unknown mode " << mode << std::endl;
        return 2;
    }
    sfu_t sfu(section);   // must abort in init()
    (void)sfu;
    return 0;             // reached only if fail-fast did not fire
}

}  // namespace

int main(int argc, char **argv) {
    if(argc > 1) {
        return run_invalid_mode(argv[1]);
    }
    test_bypass_and_zero();
    test_lane_tail();
    test_formula_and_profiles();
    test_commit_event_granularity();
    test_energy_linearity();
    test_unit_parallelism();
    test_chip_partition_identity();
    test_softmax_hand_calculation();
    test_softmax_fold_latency();
    test_serial_merge();
    test_unpriced_and_uncalibrated();
    test_activation_mapping();
    test_functional_evaluators();
    test_against_nebula_reference();
    test_supported_ops_contract();
    test_precision_mismatch();
    std::cout << "SFU validation checks passed" << std::endl;
    return 0;
}
