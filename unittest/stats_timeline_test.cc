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

} // namespace

int main() {
    run_case(memory_type_t::SEPARATE, 60.0, "separate GLB");
    run_case(memory_type_t::SHARED, 100.0, "shared GLB");
    run_order_independence();
    if(failures != 0) {
        std::cerr << failures << " stats timeline check(s) FAILED" << std::endl;
        return 1;
    }
    std::cout << "Asymmetric multi-chip stats/timeline checks passed" << std::endl;
    return 0;
}
