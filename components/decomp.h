#ifndef __DECOMP_H__
#define __DECOMP_H__

#include <cstddef>
#include <string>
#include <vector>

#include "config.h"

// Weight-decompression engine (evaluation.md Sec 4, modular-extension case study).
//
// An OPT-IN per-chip component placed on the WEIGHT supply path, between DRAM and the
// on-chip scratchpad (GLB):
//
//     DRAM (compressed weight) -> DMA -> [decompression engine] -> scratchpad -> array
//
// It models two coupled effects that decide whether compression actually helps:
//   1. weight DRAM traffic drops to 1/CR (only the compressed bytes are fetched), which
//      relieves a memory-bound layer;
//   2. the decoder produces dense weight at a bounded throughput, so a slow decoder
//      becomes the new bottleneck -- the DRAM->decoder->scratchpad pipeline (with a
//      bounded input queue and optional overlap) is modeled on the same per-tile
//      timeline the SFU streaming path uses.
//
// Design mirrors the SFU (components/sfu.{h,cc}): a config section, fail-fast validation,
// per-layer invocation events, and -- critically for the modularity claim -- ZERO change
// to the PE/scheduler. Without a [decomp] section the numbers are bit-identical to the
// dense baseline. All cost knobs are opt-in: an active-but-unpriced decoder surfaces as
// UNPRICED/UNCALIBRATED rather than a silent zero.
//
// Every quantity is at WEIGHT precision (the decoder emits dense weight in the runtime
// weight format). Input/output tensors are untouched -- only weight is compressed.

// One decompression event for a layer's weight supply. Values are FINAL (post repetition
// scaling): the layer's whole dense weight footprint arrives in `tiles` decoder tiles.
struct decomp_invocation_t {
    decomp_invocation_t();

    bool active;                    // a [decomp] section priced this layer's weight
    bool bypassed;                  // compressed >= dense for this layer -> stored dense
    size_t tiles;                   // weight tiles the decoder processes (queue units)
    size_t dense_weight_bytes;      // dense weight the decoder emits (runtime format)
    size_t compressed_weight_bytes; // what DRAM actually stores/fetches (dense/CR + meta)
    double effective_ratio;         // dense/compressed as realized (>= 1 or 1 if bypassed)

    double decoder_cycles;          // dense_bytes / decoder throughput (+ startup)
    double dram_weight_cycles_dense;    // weight DRAM transfer cost WITHOUT compression
    double dram_weight_cycles_compressed;   // ... WITH compression (the saving)

    double decoder_energy;          // per dense byte
    bool timing_calibrated;         // decoder throughput/startup were declared
    std::vector<std::string> unpriced;  // active decoder events with no declared cost
};

class decomp_t {

public:
    decomp_t(section_config_t m_section_config);
    ~decomp_t();

    // Parse and fail-fast validate the [decomp] section.
    void init(section_config_t m_section_config);
    // Print the decompression-engine specification.
    void print_specification();
    // Per-layer reset (the engine keeps no cross-layer state; kept for symmetry).
    void reset();

    // Build the decompression event for a layer whose dense weight footprint is
    // m_dense_weight_bytes, delivered over the DRAM link at m_dram_bytes_per_cycle
    // (bytes/cycle), with the compute schedule spanning m_compute_cycles (used only to
    // report the decoder ratio; the timeline integration is done by stats_t). A tile is
    // m_tile_bytes of dense weight (the decoder/queue granularity).
    decomp_invocation_t decompress(size_t m_dense_weight_bytes, double m_dram_bytes_per_cycle,
                                   double m_compute_cycles, size_t m_tile_bytes) const;

    bool enabled() const { return true; }
    double get_compression_ratio() const { return compression_ratio; }
    double get_decoder_bytes_per_cycle() const { return decoder_bytes_per_cycle; }
    double get_decoder_ratio() const { return decoder_ratio; }
    double get_startup_cycles() const { return startup_cycles; }
    unsigned get_queue_depth() const { return input_queue_depth; }
    bool get_overlap() const { return overlap; }
    double get_static_energy_per_cycle() const { return u_static_energy; }
    bool static_energy_declared() const { return static_energy_is_declared; }
    const std::string &get_profile_reference() const { return profile_reference; }

    unsigned index;

private:
    // Bytes the compressed weight occupies in DRAM for a dense footprint, honoring the
    // metadata/scale overhead and the bypass rule (never larger than dense).
    size_t compressed_bytes(size_t m_dense_bytes, bool *m_bypassed) const;

    double compression_ratio;       // CR = dense / (compressed values + metadata + scale)
    double metadata_scale_bytes_per_tile;   // fixed per-tile overhead (metadata + scale)
    double decoder_bytes_per_cycle; // dense-output throughput of the decoder (ABSOLUTE)
    // Decoder throughput RELATIVE to this layer's weight demand (evaluation.md Sec 4:
    // "decoder ratio = BW_decoder,dense-output / BW_weight-required"). When > 0, the
    // absolute throughput is derived per layer as decoder_ratio * (dense_weight_bytes /
    // compute_cycles), so one workload-independent knob sweeps every shape. Mutually
    // exclusive with decoder_bytes_per_cycle. ratio 1 = keeps pace with weight demand,
    // < 1 = decoder-limited, > 1 = decoder slack.
    double decoder_ratio;
    double startup_cycles;          // decoder pipeline fill before the first dense byte
    unsigned input_queue_depth;     // compressed-tile queue between DMA and decoder
    unsigned output_buffer_tiles;   // dense-tile buffer between decoder and scratchpad
    bool overlap;                   // DRAM/decoder/compute pipelining enabled
    std::string granularity;        // "weight_tile" only (initial contract)

    double u_decoder_energy;        // per dense byte
    double u_static_energy;         // pJ/cycle leakage of the physical engine
    bool decoder_energy_is_declared;
    bool static_energy_is_declared;

    std::string profile_reference;  // provenance of the timing/energy numbers
    bool strict_profile;            // fail instead of running an undeclared profile
};

#endif
