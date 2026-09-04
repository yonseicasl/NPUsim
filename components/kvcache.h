#ifndef __KVCACHE_H__
#define __KVCACHE_H__

#include <cstddef>
#include <string>
#include <vector>

#include "config.h"

// KV-cache TRAFFIC / SUPPLY-SCHEDULING PROXY (evaluation.md Sec 4; evaluation discussion
// 2026-09-03 Secs 6-7). Every claim built on this component must call it a proxy, not an
// architecture-level KV-cache model.
//
// In autoregressive decode every step reloads the whole KV cache from DRAM to run
// attention over the full context. That read is the dominant DRAM traffic of long-context
// decode. NPUsim's modeled scope is projection GEMMs (attention score/context are out of
// scope), so the KV read is not produced by any modeled layer. This OPT-IN component
// injects it as extra DRAM READ traffic attached to a decode step (and, in the scheduled
// modes, as a tile-level supply stage), so the projection layer it runs with sees a
// realistic memory-bound operating point WITHOUT artificially slowing the DRAM device.
//
// WHAT THE PROXY MODES DO NOT MODEL (the `attention = 1` consumer below closes most of
// this; without it every claim stays a proxy claim):
//   - no attention consumer: the scheduled modes pipeline KV tiles against the mapped
//     layer's window as a stand-in; QK / softmax / AV never execute or depend on tiles;
//   - no KV read ADDRESS stream: bytes-per-token x context is an aggregate volume, so
//     bank/row behavior of the actual KV layout is invisible;
//   - no KV WRITE and no per-token cache state: the append of the current token's K/V and
//     the cache's occupancy/eviction are absent;
//   - no layer-type awareness: the read attaches to WHATEVER mapped conv/FC layer runs, so
//     pairing [kvcache] with a transformer projection workload (one mapped layer per
//     decode step) is the CONFIG AUTHOR'S contract -- on a CNN it would be meaningless;
//   - DRAM cost per byte is estimated from the layer's own measured weight-DRAM access
//     (stats_t::apply_kv_cache_read), not from a dedicated KV stream on the link/device
//     model, so link bandwidth is only approximately reflected;
//   - no cross-layer state: each layer's read is independent; nothing persists.
//
// ATTENTION CONSUMER MODE (`attention = 1` -- evaluation discussion 2026-09-03 Sec 7):
// runs one decode step's attention as an explicit consumer appended after the mapped
// projection layer (Q depends on the projections, so the step cannot start earlier; K/V
// PREFETCH during the projection window is what double_buffered captures):
//   - QK / softmax / AV execute with declared geometry (n_q_heads, n_kv_heads, head_dim)
//     and declared throughputs; KV tiles carry a real dependency -- tile k's QK/AV cannot
//     start before tile k arrives (pipeline over the bounded kv buffer);
//   - attention_algorithm: "online" (flash-style single interleaved K/V pass, no barrier)
//     or "two_pass" (K pass -> softmax BARRIER over all scores -> V pass);
//   - the KV read/write is costed as a DEDICATED stream on the live DRAM model (device
//     read/write cycles, off-chip link, open-page row activations) -- not the weight
//     estimate. Address model: two contiguous sequential streams (K, then V); irregular
//     layouts are not modeled and the report says so;
//   - KV WRITE of the current token's K/V (the cache append) is charged, and the cache
//     occupancy after the append is reported (fail-fast if a declared kv_capacity_bytes
//     is exceeded);
//   - still one decode step at the declared context: per-layer caches are independent
//     (transformer-correct), but token-by-token cache GROWTH across steps is out of scope.
//
// It is deliberately a SECOND independent component alongside the weight-decompression
// engine (components/decomp.{h,cc}): decompression shrinks WEIGHT traffic, the KV read is a
// separate read stream. Studying them together shows the Amdahl interaction -- weight
// compression only speeds up the weight fraction of total DRAM traffic.
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

    /* Attention consumer mode */
    bool attention_enabled() const { return attention; }
    unsigned get_n_q_heads() const { return n_q_heads; }
    unsigned get_n_kv_heads() const { return n_kv_heads; }
    unsigned get_head_dim() const { return head_dim; }
    const std::string &get_attention_algorithm() const { return attention_algorithm; }
    double get_attention_macs_per_cycle() const { return attention_macs_per_cycle; }
    double get_softmax_cycles_per_element() const { return softmax_cycles_per_element; }
    size_t get_kv_capacity_bytes() const { return kv_capacity_bytes; }
    // The current token's K+V append (the cache write of one decode step).
    size_t kv_write_bytes() const { return kv_bytes_per_token; }
    // QK^T MACs of one decode step (M=1): every query head scores every cached token.
    // AV is the same count, so total attention MACs = 2x this.
    size_t attention_qk_macs() const {
        return static_cast<size_t>(n_q_heads)*head_dim*context_length;
    }
    size_t attention_softmax_elements() const {
        return static_cast<size_t>(n_q_heads)*context_length;
    }
    // Attention compute energy for one step; m_priced = both fired event classes priced.
    double attention_compute_energy(bool *m_mac_priced, bool *m_softmax_priced) const {
        if(m_mac_priced) *m_mac_priced = attention_mac_energy_is_declared;
        if(m_softmax_priced) *m_softmax_priced = softmax_energy_is_declared;
        return 2.0*static_cast<double>(attention_qk_macs())*u_attention_mac_energy +
               static_cast<double>(attention_softmax_elements())*u_softmax_energy;
    }

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

    // Attention consumer mode (attention = 1). Geometry defines both the KV footprint
    // (kv_bytes_per_token = 2 * n_kv_heads * head_dim * kv_precision_bytes -- derived, and
    // cross-checked against an explicit declaration) and the compute volume.
    bool attention;
    unsigned n_q_heads;
    unsigned n_kv_heads;
    unsigned head_dim;
    size_t kv_precision_bytes;          // bytes per K/V element (int8 = 1, fp16 = 2)
    double attention_macs_per_cycle;    // QK/AV datapath rate (0 = ideal -> UNCALIBRATED)
    double softmax_cycles_per_element;  // per attention score (0 = ideal -> UNCALIBRATED)
    std::string attention_algorithm;    // online (no barrier) | two_pass (softmax barrier)
    size_t kv_capacity_bytes;           // optional cache capacity; occupancy fail-fast

    double u_read_energy;           // per byte (of the COMPRESSED bytes fetched)
    bool read_energy_is_declared;
    double u_attention_mac_energy;  // per QK/AV MAC (attention mode)
    bool attention_mac_energy_is_declared;
    double u_softmax_energy;        // per attention score element (attention mode)
    bool softmax_energy_is_declared;

    std::string profile_reference;  // provenance of the traffic/energy numbers
};

#endif
