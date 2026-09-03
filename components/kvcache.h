#ifndef __KVCACHE_H__
#define __KVCACHE_H__

#include <cstddef>
#include <string>
#include <vector>

#include "config.h"

// KV-cache read-traffic model (evaluation.md Sec 4, decode memory-bound realism).
//
// In autoregressive decode every step reloads the whole KV cache from DRAM to run
// attention over the full context. That read is the dominant DRAM traffic of long-context
// decode -- and it is what actually makes decode memory-bound. NPUsim's modeled scope is
// projection GEMMs (attention score/context are out of scope), so the KV read is not
// produced by any modeled layer. This OPT-IN component injects it as extra DRAM READ
// traffic attached to a decode step, so the projection layer it runs with sees a realistic
// memory-bound operating point WITHOUT artificially slowing the DRAM device.
//
// It is deliberately a SECOND independent component alongside the weight-decompression
// engine (components/decomp.{h,cc}): decompression shrinks WEIGHT traffic, the KV read is a
// separate, incompressible read. Studying them together shows the Amdahl interaction --
// weight compression only speeds up the weight fraction of total DRAM traffic.
//
// Read bytes = kv_bytes_per_token * context_length. kv_bytes_per_token is the per-token KV
// footprint of ONE layer = 2 (K and V) * n_kv_heads * head_dim * precision_bytes; a GQA
// model makes it small, an MHA model large. Both are config knobs, so context_length sweeps
// cleanly. Without a [kvcache] section nothing is added and the numbers are unchanged.

struct kvcache_invocation_t {
    kvcache_invocation_t();

    bool active;                    // a [kvcache] section injected a read this step
    size_t read_bytes;              // KV cache bytes read this decode step (one layer)
    double dram_read_cycles;        // read_bytes * measured DRAM cost/byte (set by stats_t)
    double dram_read_energy;        // read_bytes * u_read_energy
    bool priced;                    // read energy was declared
    std::vector<std::string> unpriced;  // active read with no declared cost
};

class kvcache_t {

public:
    kvcache_t(section_config_t m_section_config);
    ~kvcache_t();

    void init(section_config_t m_section_config);
    void print_specification();
    void reset();

    // Total KV bytes read for one decode step of one layer.
    size_t read_bytes() const;
    // Read energy for m_read_bytes, and whether it was priced.
    double read_energy(size_t m_read_bytes, bool *m_priced) const;

    // KV-cache compression (evaluation.md Sec 4): the cache is stored/fetched compressed
    // (KV quantization / low-rank / eviction), so DRAM fetches dense/CR bytes; a decoder
    // then reconstitutes dense KV. Mirrors the weight-decompression engine but on the KV
    // read. compression_ratio 1.0 = no compression (bypass); a tile that would not shrink is
    // never enlarged.
    size_t dense_read_bytes() const;        // logical KV read (kv_bytes_per_token * context)
    size_t compressed_read_bytes() const;   // what DRAM actually fetches (dense/CR, bypassed)
    bool bypassed() const;                  // compression_ratio <= 1 -> stored dense
    // Decoder cost (cycles) to reconstitute dense KV from the compressed read. Absolute
    // throughput (dense bytes/cycle); 0 means an ideal decoder (no throughput limit),
    // reported UNCALIBRATED unless declared.
    double decoder_cycles(size_t m_dense_bytes) const;
    bool decoder_calibrated() const { return kv_decoder_bytes_per_cycle > 0.0; }

    bool enabled() const { return true; }
    size_t get_bytes_per_token() const { return kv_bytes_per_token; }
    size_t get_context_length() const { return context_length; }
    double get_compression_ratio() const { return kv_compression_ratio; }
    const std::string &get_schedule() const { return kv_schedule; }
    size_t get_tile_bytes() const { return kv_tile_bytes; }
    unsigned get_buffer_tiles() const { return kv_buffer_tiles; }
    const std::string &get_profile_reference() const { return profile_reference; }

    unsigned index;

private:
    size_t kv_bytes_per_token;      // per-token KV footprint of one layer (K+V, precision)
    size_t context_length;          // tokens already in the cache (the read length)
    double kv_compression_ratio;    // CR = dense / compressed (>= 1; 1 = no compression)
    double kv_decoder_bytes_per_cycle;  // dense-output decoder throughput (0 = ideal)
    // KV supply scheduling (evaluation discussion 2026-09-03 Sec 7). "aggregate" (default)
    // folds the read into the layer's DRAM access axis -- the traffic-sensitivity model.
    // blocking / streaming / double_buffered turn the read into a tile-level supply stage
    // against the layer window (stats_t::finalize_layer_timeline), with a bounded KV tile
    // buffer, so the same KV traffic yields different effective performance per schedule.
    std::string kv_schedule;        // aggregate | blocking | streaming | double_buffered
    size_t kv_tile_bytes;           // dense KV bytes per supply tile (streaming modes)
    unsigned kv_buffer_tiles;       // bounded KV tile buffer (default 1; double_buffered 2)

    double u_read_energy;           // per byte (of the COMPRESSED bytes fetched)
    bool read_energy_is_declared;

    std::string profile_reference;  // provenance of the traffic/energy numbers
};

#endif
