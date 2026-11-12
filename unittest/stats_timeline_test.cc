// L10: asymmetric multi-chip regression for the stats aggregation and the layer timeline.
//
// A config-only fixture CANNOT make two chips carry different per-datatype cycle profiles:
// global_buffer_t::update_tile_size() assigns every chip the SAME scheduler tile, so the
// live per-chip values are identical by construction and the correct entity/datatype
// reduction order is indistinguishable from the broken one. Making them differ needs
// per-chip mapping support, which is a model feature, not a test fixture.
//
// What IS reachable is to drive the real aggregation path with an asymmetric per-chip state:
// this test builds a stats_t, injects mirrored per-chip GLB profiles, and then runs the
// production scale_serial_repetitions() -> merge_global_buffer_fill() ->
// finalize_layer_timeline() chain. That covers what a pure reduction unit test cannot:
//   * the two GLB sides survive scaling with their DIFFERENT repetition factors (the base
//     scales uniformly, the fill per datatype) without losing the entity or datatype
//     dimension;
//   * the reduction order is the one the live code actually applies at the end of scaling;
//   * a chip that is the bottleneck only after both sides are combined still reaches the
//     stage axis, and from there the critical path.
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "stats.h"

namespace {

unsigned failures = 0;

void expect(bool condition, const std::string &message) {
    if(!condition) {
        std::cerr << "stats timeline test failed: " << message << std::endl;
        ++failures;
    }
}

void expect_close(double measured, double expected, const std::string &message) {
    if(std::fabs(measured - expected) > 1e-6) {
        std::cerr << "stats timeline test failed: " << message
                  << " (measured " << measured << ", expected " << expected << ")" << std::endl;
        ++failures;
    }
}

// Two chips with MIRRORED per-datatype GLB profiles, plus the layer state the timeline
// needs. Repetitions: uniform 2 for the base (GLB read) side, per-datatype (1, 2, 1) for the
// fill (multi-chip -> GLB write) side.
//
//   chip 0 base (10,  0, 0)  fill ( 0, 20, 0)
//   chip 1 base ( 0, 30, 0)  fill (40,  0, 0)
// after scaling:
//   chip 0 base (20,  0, 0)  fill ( 0, 40, 0)   -> per-datatype totals (20, 40, 0)
//   chip 1 base ( 0, 60, 0)  fill (40,  0, 0)   -> per-datatype totals (40, 60, 0)
//
// SEPARATE GLB (per-datatype partitions run concurrently -> max over datatypes):
//   correct: max_chip(max_type(base+fill)) = max(40, 60) = 60
//   broken : max_chip(max_type(base) + max_type(fill)) = max(20+40, 60+40) = 100
// SHARED GLB (one port serializes the streams -> sum over datatypes):
//   both forms agree at max(60, 100) = 100, which is exactly why the separate-buffer error
//   stayed hidden before L1.
void run_case(memory_type_t buffer_type, double expected_axis, const std::string &label) {
    stats_t stats;
    stats.global_buffer_type = buffer_type;
    stats.chip_access_cycle_global_buffer.clear();
    stats.chip_fill_access_cycle_global_buffer.clear();
    const double base[2][3] = {{10.0, 0.0, 0.0}, {0.0, 30.0, 0.0}};
    const double fill[2][3] = {{0.0, 20.0, 0.0}, {40.0, 0.0, 0.0}};
    for(unsigned chip = 0; chip < 2; ++chip) {
        stats.chip_access_cycle_global_buffer.push_back(
            std::vector<double>(base[chip], base[chip] + 3));
        stats.chip_fill_access_cycle_global_buffer.push_back(
            std::vector<double>(fill[chip], fill[chip] + 3));
        for(unsigned type = 0; type < 3; ++type) {
            // The flat per-datatype vectors the report prints are max-across-chips; keep them
            // consistent with the injected matrices so the axis bound (T10) stays meaningful.
            stats.access_cycle_global_buffer[type] =
                std::max(stats.access_cycle_global_buffer[type], base[chip][type]);
            stats.fill_access_cycle_global_buffer[type] =
                std::max(stats.fill_access_cycle_global_buffer[type], fill[chip][type]);
        }
    }
    // A compute schedule and a tile clock, so finalize_layer_timeline() has a real layer to
    // place on the timeline rather than an all-zero one.
    stats.computation_cycle = 100.0;
    stats.timeline_physical_macs = 1.0;
    stats.timeline_boundary_overlap[0] = false;
    stats.timeline_boundary_overlap[1] = false;
    stats.timeline_boundary_overlap[2] = true;
    stats.timeline_boundary_overlap[3] = true;

    std::vector<unsigned> datatype_repetitions(3, 1);
    datatype_repetitions[1] = 2;
    stats.scale_serial_repetitions(2, datatype_repetitions);

    expect_close(stats.stage_axis_access[2], expected_axis,
                 label + ": GLB access axis after scaling");
    // The reduction must be the max ACROSS chips of each chip's own combined total, so it
    // can never exceed the datatype combination of the per-datatype maxima, and it must be
    // at least the busiest single chip-and-datatype value.
    expect(stats.busy_cycle_global_buffer >= expected_axis - 1e-6,
           label + ": GLB busy must cover its access axis");
    // A chip that only becomes the bottleneck once both sides are added still has to reach
    // the critical path.
    expect(stats.layer_latency >= expected_axis - 1e-6,
           label + ": critical path must cover the GLB axis");
    expect(stats.layer_latency >= stats.stage_axis_compute - 1e-6,
           label + ": critical path must cover the compute schedule");
}

// The same injection with the chips SWAPPED must give the same answer: the reduction is a
// max across entities, so the chip order cannot matter.
void run_order_independence() {
    stats_t a, b;
    for(unsigned pass = 0; pass < 2; ++pass) {
        stats_t &stats = (pass == 0) ? a : b;
        stats.global_buffer_type = memory_type_t::SEPARATE;
        const double base[2][3] = {{10.0, 0.0, 0.0}, {0.0, 30.0, 0.0}};
        const double fill[2][3] = {{0.0, 20.0, 0.0}, {40.0, 0.0, 0.0}};
        for(unsigned i = 0; i < 2; ++i) {
            const unsigned chip = (pass == 0) ? i : (1 - i);
            stats.chip_access_cycle_global_buffer.push_back(
                std::vector<double>(base[chip], base[chip] + 3));
            stats.chip_fill_access_cycle_global_buffer.push_back(
                std::vector<double>(fill[chip], fill[chip] + 3));
        }
        stats.computation_cycle = 100.0;
        stats.timeline_physical_macs = 1.0;
        std::vector<unsigned> datatype_repetitions(3, 1);
        datatype_repetitions[1] = 2;
        stats.scale_serial_repetitions(2, datatype_repetitions);
    }
    expect_close(a.stage_axis_access[2], b.stage_axis_access[2],
                 "chip order changed the GLB access axis");
}

// Phase-5 (plan_sfu.md): SFU streaming overlap on the per-tile timeline, against hand
// calculations. Producer: PE-only, per-tile cost p = 10 over T = 4 tiles (pre-scale
// computation_cycle = 10, scaled x4), producer-only latency 40.
//
//   serial (depth 1), S = 20                : latency 40 + 20 = 60, all exposed, stall 20
//   streaming (depth 2), S = 8, s = 2 < p   : latency = T*p + s = 42, exposed 2, hidden 6,
//                                             no back-pressure
//   streaming (depth 2), S = 80, s = 20 > p : latency = p + T*s = 90 (fill + bottleneck),
//                                             exposed 50, hidden 30; the PE start times
//                                             are pushed by finish_sfu(k-2), which delays
//                                             tiles 2 and 3 by 10 each -> stall 20
//   streaming (depth 2), S = 8, but ALL tiles reduce into ONE output (output reps 1):
//                                             the commit releases on the LAST tile, so
//                                             latency = T*p + S = 48 and nothing hides
void run_sfu_streaming_cases() {
    struct sfu_case_t {
        unsigned queue_depth;
        double sfu_busy;
        unsigned output_repetitions;
        double expected_latency;
        double expected_serial;
        double expected_hidden;
        double expected_stall;
        const char *label;
    };
    const sfu_case_t cases[] = {
        {1, 20.0, 4, 60.0, 20.0,  0.0, 20.0, "serial contract appends the SFU window"},
        {2,  8.0, 4, 42.0,  2.0,  6.0,  0.0, "fast streaming SFU hides behind the producer"},
        {2, 80.0, 4, 90.0, 50.0, 30.0, 20.0, "slow streaming SFU back-pressures the producer"},
        {2,  8.0, 1, 48.0,  8.0,  0.0,  0.0, "reduction-bunched outputs release at the end"},
    };
    for(unsigned i = 0; i < sizeof(cases)/sizeof(cases[0]); ++i) {
        const sfu_case_t &sfu_case = cases[i];
        stats_t stats;
        stats.computation_cycle = 10.0;
        stats.timeline_physical_macs = 1.0;

        sfu_invocation_t invocation;
        invocation.operation = "relu";
        invocation.valid_elements = 64;
        invocation.operations = 64;
        invocation.busy_cycle = sfu_case.sfu_busy;
        stats.set_sfu_activation(invocation, 1, 16, 0.0, sfu_case.queue_depth);

        std::vector<unsigned> datatype_repetitions(3, 1);
        datatype_repetitions[data_type_t::OUTPUT] = sfu_case.output_repetitions;
        stats.scale_serial_repetitions(4, datatype_repetitions);

        const std::string label(sfu_case.label);
        expect_close(stats.layer_latency, sfu_case.expected_latency,
                     label + ": critical-path latency");
        expect_close(stats.sfu_busy_cycle, sfu_case.sfu_busy, label + ": SFU busy window");
        expect_close(stats.sfu_serial_cycle, sfu_case.expected_serial,
                     label + ": exposed (critical-path) SFU cycles");
        expect_close(stats.sfu_hidden_cycle, sfu_case.expected_hidden,
                     label + ": hidden SFU cycles");
        expect_close(stats.sfu_stall_cycle, sfu_case.expected_stall,
                     label + ": producer back-pressure stall");
        // Visibility identity: the window is exactly exposed + hidden.
        expect_close(stats.sfu_serial_cycle + stats.sfu_hidden_cycle, sfu_case.sfu_busy,
                     label + ": exposed + hidden must equal the SFU window");
    }
}

// [accelerator] timeline_dram_tiles: the stage totals of a layer whose DRAM-level loop ran `outer` times belong
// on a tile clock of (GLB repetitions x outer) tiles. A DRAM stage of 100 cycles in front of a compute stage of
// 200, both double-buffered, on 2 GLB repetitions:
//   outer 1 -> 2 tiles of (50, 100): the first tile's DRAM fill (50) + 2 x 100             = 250
//   outer 4 -> 8 tiles of (12.5, 25): the same totals, the fill is one eighth (12.5) + 200 = 212.5
// The fill shrinks with the real tile count instead of growing with the DRAM iterations; the totals, and so
// every busy axis, are unchanged.
void run_outer_iterations() {
    double latency[2] = {0.0, 0.0};
    for(unsigned pass = 0; pass < 2; ++pass) {
        stats_t stats;
        stats.computation_cycle = 100.0;                      // x 2 repetitions = 200
        stats.access_cycle_dram[data_type_t::WEIGHT] = 50.0;  // x 2 repetitions = 100
        stats.timeline_physical_macs = 1.0;
        for(unsigned b = 0; b < 4; ++b) stats.timeline_boundary_overlap[b] = true;
        stats.timeline_outer_iterations = (pass == 0) ? 1 : 4;
        std::vector<unsigned> datatype_repetitions(3, 1);
        datatype_repetitions[data_type_t::WEIGHT] = 2;
        stats.scale_serial_repetitions(2, datatype_repetitions);
        expect_close(stats.busy_cycle_dram, 100.0, "outer iterations: DRAM busy total");
        expect_close(stats.stage_axis_compute, 200.0, "outer iterations: compute schedule");
        latency[pass] = stats.layer_latency;
    }
    expect_close(latency[0], 250.0, "outer iterations 1: fill of a 2-tile clock");
    expect_close(latency[1], 212.5, "outer iterations 4: fill of an 8-tile clock");
}

// [accelerator] offchip_window_count / window_switch_cycle: the off-chip stage delivers n windows, each whole
// before its on-chip use, the next fetched meanwhile; w is a compute-side cost per window. With a DRAM stage
// (total S) in front of a compute stage (total C) and nothing in between, the makespan is the two-stage flow
// shop of n identical jobs (d, c) = (S/n, C/n + w):  d + (n - 1) max(d, c) + c.  Expected values below are that
// formula worked by hand; the clock has `tiles` tiles, so a window spans tiles/n of them.
double last_pe_busy = 0.0;      // busy_cycle_pe of the last run_windows() call

double run_windows(unsigned tiles, unsigned windows, double dram_total, double compute_total, double w,
                   double link_total = 0.0, unsigned outstanding = 0, double lb_total = 0.0,
                   bool pe_double_buffer = true) {
    stats_t stats;
    stats.timeline_offchip_outstanding = outstanding;
    stats.entity_combined_access_lb = lb_total/tiles;
    stats.pe_double_buffer = pe_double_buffer;
    stats.computation_cycle = compute_total/tiles;
    stats.access_cycle_dram[data_type_t::WEIGHT] = dram_total/tiles;
    stats.access_cycle_multi_chip[data_type_t::WEIGHT] = link_total/tiles;
    stats.timeline_physical_macs = 1.0;
    for(unsigned b = 0; b < 4; ++b) stats.timeline_boundary_overlap[b] = true;
    stats.timeline_offchip_windows = windows;
    stats.window_switch_cycle = w;
    std::vector<unsigned> datatype_repetitions(3, 1);
    datatype_repetitions[data_type_t::WEIGHT] = tiles;
    stats.scale_serial_repetitions(tiles, datatype_repetitions);
    expect_close(stats.busy_cycle_dram, dram_total, "windows: DRAM busy total is conserved");
    expect_close(stats.busy_cycle_multi_chip, link_total, "windows: multi-chip busy total is conserved");
    expect_close(stats.stage_axis_compute, compute_total + w*windows, "windows: compute schedule = compute + n w");
    expect(stats.timeline_boundary_depth[0] == (outstanding ? outstanding : 2u),
           "windows: DRAM <-> multi-chip depth left as configured");
    last_pe_busy = stats.busy_cycle_pe;
    return stats.layer_latency;
}

void run_window_cases() {
    // Off (0 windows): the off-chip stage serves every tile, 8 jobs of (10, 5): 10 + 7*10 + 5.
    expect_close(run_windows(8, 0, 80.0, 40.0, 0.0), 85.0, "windows off: per-tile service");
    // DMA-bound, 2 windows of (40, 20): 40 + 40 + 20. The first fetch and the last compute are exposed.
    expect_close(run_windows(8, 2, 80.0, 40.0, 0.0), 100.0, "2 windows, DMA-bound");
    // Compute-bound, 2 windows of (20, 80): 20 + 80 + 80.
    expect_close(run_windows(8, 2, 40.0, 160.0, 0.0), 180.0, "2 windows, compute-bound");
    // 3 windows of (10, 40): 10 + 2*40 + 40 -- the third fetch must not wait for the first window's compute.
    expect_close(run_windows(6, 3, 30.0, 120.0, 0.0), 130.0, "3 windows, compute-bound");
    // One window: fetch, then compute, nothing overlaps: 80 + 40.
    expect_close(run_windows(8, 1, 80.0, 40.0, 0.0), 120.0, "1 window: serial");
    // w hidden under DMA slack: 2 windows of (40, 20 + 10): 40 + 40 + 30 -- only the last window's w shows.
    expect_close(run_windows(8, 2, 80.0, 40.0, 10.0), 110.0, "w hidden when DMA-bound");
    // w exposed when compute-bound: 2 windows of (20, 80 + 10): 20 + 90 + 90.
    expect_close(run_windows(8, 2, 40.0, 160.0, 10.0), 200.0, "w exposed when compute-bound");
    // Crossover: (30, 25) is DMA-bound, 30 + 3*30 + 25 = 145; with w = 10 it is (30, 35), compute-bound:
    // 30 + 3*35 + 35 = 170 -- more than the 4*10 of added work minus slack would suggest per window.
    expect_close(run_windows(8, 4, 120.0, 100.0, 0.0), 145.0, "crossover, w = 0");
    expect_close(run_windows(8, 4, 120.0, 100.0, 10.0), 170.0, "crossover, w changes the bottleneck");
    // The count matters, not only the total added cost (Codex's counterexample, S 8, C 4):
    //   4 windows of (2, 1 + 1): 2 + 3*2 + 2 = 10;   2 windows of (4, 2 + 2), the same n w = 4: 4 + 4 + 4 = 12.
    expect_close(run_windows(8, 4, 8.0, 4.0, 1.0), 10.0, "4 windows, w 1");
    expect_close(run_windows(8, 2, 8.0, 4.0, 2.0), 12.0, "2 windows, w 2: same total, different latency");
    // Partial windows: 5 tiles in 2 windows span 3 and 2 tiles (compute 10 per tile), each fetch 40:
    //   window 0: fetch 0-40, compute 40-70; window 1: fetch 40-80, compute 80-100.
    expect_close(run_windows(5, 2, 80.0, 50.0, 0.0), 100.0, "uneven windows (3 + 2 tiles)");
    // A multi-chip link behind the DRAM carries the same window at the same time (cut-through), so a window's
    // fetch is max(DRAM, link), not the sum:  DRAM 40, link 24 per window -> as without the link, 100;
    // link-bound, DRAM 20, link 40, compute 20 per window: 40 + 40 + 20 = 100 (store-and-forward: 120).
    expect_close(run_windows(8, 2, 80.0, 40.0, 0.0, 48.0), 100.0, "2 windows, link hidden behind DRAM");
    expect_close(run_windows(8, 2, 40.0, 40.0, 0.0, 80.0), 100.0, "2 windows, link-bound");
    // Without windows the two stages pipeline per tile as before: 8 jobs through (10, 6, 5): 10 + 6 + 7*10 + 5.
    expect_close(run_windows(8, 0, 80.0, 40.0, 0.0, 48.0), 91.0, "windows off: DRAM and link pipeline per tile");
    // w is one fixed cost per window, not a share per tile (Codex's case): 5 tiles in windows of 3 and 2,
    // compute 10 per tile, fetch 40 each, w 5:  window 0 fetch 0-40, compute 40-75 (30 + 5);
    //   window 1 fetch 40-80, compute 80-105 (20 + 5).  Spread over the tiles it would be 104.
    expect_close(run_windows(5, 2, 80.0, 50.0, 5.0), 105.0, "uneven windows pay the same w each");
    // w joins the compute axis of its tile, inside the PE stage's own combination. 5 tiles (3 + 2), fetch 40
    // per window, LB transfer 10 per tile:
    //   compute 1 per tile, w 5 (Codex's case): every tile is LB-bound (10 > 1 + 5), w fully hidden:
    //     window 0 40-70, window 1 80-100.
    expect_close(run_windows(5, 2, 80.0, 5.0, 5.0, 0.0, 0, 50.0), 100.0, "w hidden under a longer LB axis");
    //   compute 8 per tile, w 5: the two last tiles are max(13, 10), the others 10: 40-73, then 80-103.
    expect_close(run_windows(5, 2, 80.0, 40.0, 5.0, 0.0, 0, 50.0), 103.0, "w partly exposed over the LB axis");
    // ... and the PE stage's reported busy time is what was scheduled, 10 + 10 + 13 + 10 + 13, not max(45, 50).
    expect_close(last_pe_busy, 56.0, "PE busy = sum of the scheduled per-tile service");
    run_windows(5, 2, 80.0, 5.0, 5.0, 0.0, 0, 50.0);
    expect_close(last_pe_busy, 50.0, "PE busy when w is fully hidden under the LB axis");
    run_windows(5, 2, 80.0, 40.0, 0.0, 0.0, 0, 50.0);
    expect_close(last_pe_busy, 50.0, "PE busy without a window cost: max of the axis totals, as before");
    //   single-buffered LB (transfer serializes with compute): tiles 11, 11, 16 | 11, 16: 40-78, then 80-107.
    expect_close(run_windows(5, 2, 80.0, 5.0, 5.0, 0.0, 0, 50.0, false), 107.0, "single-buffered LB: w serial");
    // As many windows as tiles is the per-tile service: 4 jobs of (10, 5): 10 + 3*10 + 5.
    // (More windows than tiles is refused with an error, as is a single-buffered or bypassed chip ingress.)
    expect_close(run_windows(4, 4, 40.0, 20.0, 0.0), 45.0, "one window per tile");
    // An explicit outstanding-request limit on the DRAM <-> multi-chip boundary is kept as configured; the
    // window double buffer sits at the chip ingress, so the DMA-bound 2-window case is still 100.
    expect_close(run_windows(8, 2, 80.0, 40.0, 0.0, 0.0, 1), 100.0, "explicit outstanding limit 1 respected");
    expect_close(run_windows(8, 2, 80.0, 40.0, 0.0, 0.0, 16), 100.0, "explicit outstanding limit 16 respected");
}

// What [accelerator] fused_output_consumer applies: apply_graph_residency(input off chip, output on chip). The
// output's off-chip traffic (multi-chip link and DRAM: requests, transfers, payload transactions, cycles) is
// removed; the input's and the weight's are untouched; everything on chip -- including the output's own GLB
// accesses -- is untouched. No capacity enters: the caller decides, here unconditionally.
void run_fused_output_consumer() {
    stats_t stats;
    for(unsigned t = 0; t < data_type_t::NUM_DATA_TYPES; ++t) {
        stats.num_data_transfer_dram[t] = 100 + t;
        stats.payload_link_transactions_dram[t] = 200 + t;
        stats.access_cycle_dram[t] = 300.0 + t;
        stats.transfer_cycle_dram[t] = 400.0 + t;
        stats.payload_link_transactions_multi_chip[t] = 500 + t;
        stats.access_cycle_multi_chip[t] = 600.0 + t;
        stats.access_cycle_global_buffer[t] = 700.0 + t;
        stats.fill_access_cycle_global_buffer[t] = 800.0 + t;
    }
    stats.apply_graph_residency(false, true, "test");
    const unsigned in = data_type_t::INPUT, wt = data_type_t::WEIGHT, out = data_type_t::OUTPUT;
    for(unsigned t : {in, wt}) {
        expect(stats.num_data_transfer_dram[t] == 100 + t && stats.payload_link_transactions_dram[t] == 200 + t,
               "fused output consumer: input / weight DRAM transfers unchanged");
        expect_close(stats.access_cycle_dram[t], 300.0 + t, "fused output consumer: input / weight DRAM cycles");
        expect_close(stats.transfer_cycle_dram[t], 400.0 + t, "fused output consumer: input / weight DRAM link");
        expect(stats.payload_link_transactions_multi_chip[t] == 500 + t, "fused output consumer: input / weight link");
        expect_close(stats.access_cycle_multi_chip[t], 600.0 + t, "fused output consumer: input / weight link cycles");
        expect_close(stats.fill_access_cycle_global_buffer[t], 800.0 + t, "fused output consumer: input fill kept");
    }
    expect(stats.num_data_transfer_dram[out] == 0 && stats.payload_link_transactions_dram[out] == 0,
           "fused output consumer: no output DRAM transfers");
    expect_close(stats.access_cycle_dram[out] + stats.transfer_cycle_dram[out], 0.0, "fused output consumer: no output DRAM cycles");
    expect(stats.payload_link_transactions_multi_chip[out] == 0, "fused output consumer: no output link transfers");
    expect_close(stats.access_cycle_multi_chip[out], 0.0, "fused output consumer: no output link cycles");
    for(unsigned t = 0; t < data_type_t::NUM_DATA_TYPES; ++t) {
        expect_close(stats.access_cycle_global_buffer[t], 700.0 + t, "fused output consumer: on-chip GLB accesses kept");
    }
}

} // namespace

int main() {
    run_case(memory_type_t::SEPARATE, 60.0, "separate GLB");
    run_case(memory_type_t::SHARED, 100.0, "shared GLB");
    run_order_independence();
    run_sfu_streaming_cases();
    run_outer_iterations();
    run_window_cases();
    run_fused_output_consumer();
    if(failures != 0) {
        std::cerr << failures << " stats timeline check(s) FAILED" << std::endl;
        return 1;
    }
    std::cout << "Asymmetric multi-chip stats/timeline checks passed" << std::endl;
    return 0;
}
