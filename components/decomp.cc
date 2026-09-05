#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "decomp.h"
#include "energy_units.h"

namespace {
size_t ceil_div(size_t m_value, size_t m_divisor) {
    return m_divisor == 0 ? 0 : (m_value + m_divisor - 1)/m_divisor;
}

// Deterministic unit-interval hash (splitmix64 finalizer): the per-tile compression
// factors must be exactly reproducible across runs and platforms, so no libc rand().
double unit_hash(unsigned long long m_seed, unsigned long long m_index) {
    unsigned long long z = m_seed + 0x9e3779b97f4a7c15ULL*(m_index + 1ULL);
    z = (z ^ (z >> 30))*0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27))*0x94d049bb133111ebULL;
    z ^= z >> 31;
    return static_cast<double>(z >> 11)/9007199254740992.0;   // [0, 1), 53-bit
}
}  // namespace

decomp_invocation_t::decomp_invocation_t() :
    active(false),
    bypassed(false),
    tiles(0),
    bypassed_tiles(0),
    dense_weight_bytes(0),
    compressed_weight_bytes(0),
    effective_ratio(1.0),
    decoder_cycles(0.0),
    dram_weight_cycles_dense(0.0),
    dram_weight_cycles_compressed(0.0),
    decoder_energy(0.0),
    timing_calibrated(true) {
}

decomp_t::decomp_t(section_config_t m_section_config) :
    index(0),
    compression_ratio(2.0),
    metadata_scale_bytes_per_tile(0.0),
    decoder_bytes_per_cycle(0.0),
    decoder_ratio(0.0),
    startup_cycles(0.0),
    tile_ratio_cv(0.0),
    tile_ratio_seed(1),
    input_queue_depth(4),
    output_buffer_tiles(2),
    overlap(true),
    granularity("weight_tile"),
    u_decoder_energy(0.0),
    u_static_energy(0.0),
    decoder_energy_is_declared(false),
    static_energy_is_declared(false),
    strict_profile(false) {

    init(m_section_config);
}

decomp_t::~decomp_t() {
}

void decomp_t::init(section_config_t m_section_config) {
    m_section_config.get_setting("compression_ratio", &compression_ratio);
    m_section_config.get_setting("metadata_scale_bytes_per_tile", &metadata_scale_bytes_per_tile);
    m_section_config.get_setting("decoder_bytes_per_cycle", &decoder_bytes_per_cycle);
    m_section_config.get_setting("decoder_ratio", &decoder_ratio);
    m_section_config.get_setting("startup_cycles", &startup_cycles);
    m_section_config.get_setting("tile_ratio_cv", &tile_ratio_cv);
    m_section_config.get_setting("tile_ratio_seed", &tile_ratio_seed);
    m_section_config.get_setting("input_queue_depth", &input_queue_depth);
    m_section_config.get_setting("output_buffer_tiles", &output_buffer_tiles);
    m_section_config.get_setting("overlap", &overlap);
    m_section_config.get_setting("granularity", &granularity);
    m_section_config.get_setting("strict_profile", &strict_profile);

    // Fail-fast contract checks (evaluation.md: knob validity is stated, not defaulted).
    if(compression_ratio <= 0.0) {
        std::cerr << "Error: [decomp] compression_ratio must be positive" << std::endl;
        exit(1);
    }
    if(metadata_scale_bytes_per_tile < 0.0) {
        std::cerr << "Error: [decomp] metadata_scale_bytes_per_tile must not be negative"
                  << std::endl;
        exit(1);
    }
    if(startup_cycles < 0.0) {
        std::cerr << "Error: [decomp] startup_cycles must not be negative" << std::endl;
        exit(1);
    }
    if(decoder_ratio < 0.0) {
        std::cerr << "Error: [decomp] decoder_ratio must not be negative" << std::endl;
        exit(1);
    }
    // cv >= 1 would allow a non-positive per-tile compressed size (factor 1 - cv <= 0).
    if(tile_ratio_cv < 0.0 || tile_ratio_cv >= 1.0) {
        std::cerr << "Error: [decomp] tile_ratio_cv must be in [0, 1)" << std::endl;
        exit(1);
    }
    // decoder_ratio (relative to weight demand) and decoder_bytes_per_cycle (absolute) are
    // two ways to state the SAME quantity; declaring both is ambiguous.
    if(decoder_ratio > 0.0 && decoder_bytes_per_cycle > 0.0) {
        std::cerr << "Error: [decomp] set either decoder_ratio (relative to weight demand)"
                  << " or decoder_bytes_per_cycle (absolute), not both" << std::endl;
        exit(1);
    }
    if(input_queue_depth == 0) {
        std::cerr << "Error: [decomp] input_queue_depth must be at least 1" << std::endl;
        exit(1);
    }
    if(output_buffer_tiles == 0) {
        std::cerr << "Error: [decomp] output_buffer_tiles must be at least 1" << std::endl;
        exit(1);
    }
    if(granularity != "weight_tile") {
        std::cerr << "Error: [decomp] granularity = '" << granularity
                  << "' is not supported; the initial contract is weight_tile" << std::endl;
        exit(1);
    }
    // decoder_bytes_per_cycle <= 0 means "no throughput limit modeled": the decoder keeps
    // up with any weight demand. That is a legitimate ideal (infinite decoder), but under
    // strict_profile it must be declared explicitly.
    if(strict_profile && decoder_bytes_per_cycle <= 0.0 && decoder_ratio <= 0.0) {
        std::cerr << "Error: [decomp] strict_profile = 1 requires a positive"
                  << " decoder_bytes_per_cycle or decoder_ratio" << std::endl;
        exit(1);
    }

    decoder_energy_is_declared = m_section_config.get_setting("decomp_decoder_energy",
                                                              &u_decoder_energy);
    static_energy_is_declared = m_section_config.get_setting("decomp_static_energy",
                                                             &u_static_energy);
    m_section_config.get_setting("profile_reference", &profile_reference);
}

void decomp_t::reset() {
    // The engine carries no cross-layer state.
}

size_t decomp_t::compressed_bytes(size_t m_dense_bytes, size_t m_tiles,
                                  bool *m_bypassed) const {
    if(m_dense_bytes == 0) {
        if(m_bypassed) *m_bypassed = false;
        return 0;
    }
    // Compressed values plus the fixed metadata/scale overhead of EVERY tile (the knob is
    // per-tile, so a layer of N tiles carries N overheads -- charging it once per layer
    // understated the overhead by the tile count and pushed the break-even point of
    // low-CR configurations the wrong way). An effective compression_ratio that already
    // folds metadata in (evaluation.md: CR includes metadata + scale) simply leaves the
    // knob at its 0 default.
    const double values = static_cast<double>(m_dense_bytes)/compression_ratio;
    const double total = values +
                         metadata_scale_bytes_per_tile*static_cast<double>(std::max<size_t>(1, m_tiles));
    const bool bypass = total >= static_cast<double>(m_dense_bytes);
    if(m_bypassed) *m_bypassed = bypass;
    // Bypass: weight that does not shrink is stored dense (evaluation.md).
    return bypass ? m_dense_bytes : static_cast<size_t>(std::ceil(total));
}

decomp_invocation_t decomp_t::decompress(size_t m_dense_weight_bytes,
                                         double m_dram_bytes_per_cycle,
                                         double m_compute_cycles,
                                         size_t m_tile_bytes) const {
    decomp_invocation_t inv;
    inv.active = true;
    if(m_dense_weight_bytes == 0) {
        return inv;
    }
    inv.dense_weight_bytes = m_dense_weight_bytes;
    const size_t tile_bytes = std::max<size_t>(1, m_tile_bytes);
    inv.tiles = ceil_div(m_dense_weight_bytes, tile_bytes);

    if(tile_ratio_cv > 0.0 && inv.tiles > 1) {
        // Per-tile compression model (evaluation discussion 2026-09-03 Sec 5): each weight
        // tile compresses by its own deterministic factor around the layer mean, the layer
        // TOTAL is pinned to dense/CR pre-bypass (so uniform vs varied runs move identical
        // DRAM bytes), and a tile whose compressed size would reach dense is stored dense
        // (the PER-TILE bypass evaluation.md specifies -- the layer-level path below can
        // only bypass whole layers). The resulting non-uniform arrival is what a bounded
        // decoder queue turns into backpressure.
        const size_t last_tile = m_dense_weight_bytes - (inv.tiles - 1)*tile_bytes;
        const double target_total = static_cast<double>(m_dense_weight_bytes)/compression_ratio;
        double raw_sum = 0.0;
        for(size_t i = 0; i < inv.tiles; ++i) {
            const double dense_i = static_cast<double>(i + 1 == inv.tiles ? last_tile : tile_bytes);
            raw_sum += dense_i/compression_ratio*
                       (1.0 + tile_ratio_cv*(2.0*unit_hash(tile_ratio_seed, i) - 1.0));
        }
        const double scale = raw_sum > 0.0 ? target_total/raw_sum : 1.0;
        const size_t chunks = std::min<size_t>(inv.tiles, 4096);
        inv.tile_supply_fraction.assign(chunks, 0.0);
        double realized = 0.0;
        for(size_t i = 0; i < inv.tiles; ++i) {
            const double dense_i = static_cast<double>(i + 1 == inv.tiles ? last_tile : tile_bytes);
            double comp_i = dense_i/compression_ratio*
                            (1.0 + tile_ratio_cv*(2.0*unit_hash(tile_ratio_seed, i) - 1.0))*scale;
            if(comp_i >= dense_i) {   // per-tile bypass: this tile did not shrink
                comp_i = dense_i;
                ++inv.bypassed_tiles;
            }
            realized += comp_i;
            inv.tile_supply_fraction[i*chunks/inv.tiles] += comp_i;
        }
        const double total = realized +
            metadata_scale_bytes_per_tile*static_cast<double>(inv.tiles);
        inv.bypassed = total >= static_cast<double>(m_dense_weight_bytes);
        inv.compressed_weight_bytes = inv.bypassed
            ? m_dense_weight_bytes : static_cast<size_t>(std::ceil(total));
        if(inv.bypassed || realized <= 0.0) {
            inv.tile_supply_fraction.clear();
        } else {
            for(size_t j = 0; j < chunks; ++j) inv.tile_supply_fraction[j] /= realized;
        }
    } else {
        bool bypassed = false;
        inv.compressed_weight_bytes = compressed_bytes(m_dense_weight_bytes, inv.tiles,
                                                       &bypassed);
        inv.bypassed = bypassed;
    }
    inv.effective_ratio = static_cast<double>(m_dense_weight_bytes)/
                          static_cast<double>(std::max<size_t>(1, inv.compressed_weight_bytes));

    // DRAM weight transfer cost: dense vs compressed, at the DRAM link rate. The saving
    // is the difference -- what the compression buys on the memory side.
    if(m_dram_bytes_per_cycle > 0.0) {
        inv.dram_weight_cycles_dense =
            static_cast<double>(m_dense_weight_bytes)/m_dram_bytes_per_cycle;
        inv.dram_weight_cycles_compressed =
            static_cast<double>(inv.compressed_weight_bytes)/m_dram_bytes_per_cycle;
    }

    // Decoder cost: it emits the full dense weight, plus a one-time pipeline fill. A
    // bypassed layer runs no decoder (the weight is already dense).
    //
    // Two throughput modes:
    //   - ABSOLUTE (decoder_bytes_per_cycle > 0): decoder_cycles computed here.
    //   - RELATIVE (decoder_ratio > 0): the throughput is decoder_ratio x this layer's
    //     weight demand (evaluation.md Sec 4). Weight demand = dense_bytes / compute
    //     schedule, and the FINAL compute schedule is only known after repetition scaling
    //     and timeline assembly -- so this mode's decoder_cycles is filled in by
    //     stats_t::finalize_layer_timeline() (which owns the final compute schedule). The
    //     engine still reports the mode as calibrated because the ratio was declared.
    if(!inv.bypassed) {
        if(decoder_bytes_per_cycle > 0.0) {
            inv.decoder_cycles = startup_cycles +
                static_cast<double>(m_dense_weight_bytes)/decoder_bytes_per_cycle;
        } else if(decoder_ratio > 0.0) {
            inv.decoder_cycles = 0.0;   // deferred to finalize_layer_timeline (see above)
        } else {
            // Ideal (no limit) and undeclared: flag as not calibration-grade.
            inv.timing_calibrated = false;
        }
    }

    (void)m_compute_cycles;   // relative-throughput mode uses the final schedule in stats_t

    // Decoder energy: per dense byte emitted (bypassed layers decode nothing).
    if(!inv.bypassed) {
        inv.decoder_energy = static_cast<double>(m_dense_weight_bytes)*u_decoder_energy;
        if(!decoder_energy_is_declared) {
            inv.unpriced.push_back("decoder decode fired but [decomp] decomp_decoder_energy"
                                   " is not declared");
        }
    }

    return inv;
}

void decomp_t::print_specification() {
    std::cout << "======= Decompression engine ========" << std::endl;
    std::cout << "Compression ratio  :" << std::setw(20) << std::setprecision(2)
              << std::fixed << compression_ratio << "x"
              << (tile_ratio_cv > 0.0 ? "  (per-tile variation below)" : "") << std::endl;
    if(tile_ratio_cv > 0.0) {
        std::cout << "Per-tile variation :" << std::setw(19) << std::setprecision(2)
                  << tile_ratio_cv << " +-relative spread (seed " << tile_ratio_seed
                  << ", layer total pinned)" << std::endl;
    }
    if(decoder_ratio > 0.0) {
        std::cout << "Decoder throughput :" << std::setw(18) << std::setprecision(2)
                  << decoder_ratio << "x weight demand (per-layer)" << std::endl;
    } else {
        std::cout << "Decoder throughput :" << std::setw(18) << std::setprecision(2)
                  << decoder_bytes_per_cycle << " dense B/cycle"
                  << (decoder_bytes_per_cycle <= 0.0 ? "  (ideal: no throughput limit)" : "")
                  << std::endl;
    }
    std::cout << "Startup / queue    :" << std::setw(18) << std::setprecision(0)
              << startup_cycles << " cyc / " << input_queue_depth << " tiles" << std::endl;
    std::cout << "Output buffer      :" << std::setw(18) << output_buffer_tiles
              << " tiles" << std::endl;
    std::cout << "Overlap            :" << std::setw(18)
              << (overlap ? "enabled" : "disabled") << std::endl;
    std::cout << "Metadata/scale     :" << std::setw(18) << std::setprecision(1)
              << metadata_scale_bytes_per_tile << " B/tile" << std::endl;
    std::cout << "Decoder energy     :" << std::setw(18) << std::setprecision(4)
              << u_decoder_energy << "/" << u_static_energy << " (byte/static) "
              << energy_units().label() << std::endl;
    std::cout << "Provenance         : "
              << (profile_reference.empty()
                  ? "NOT DECLARED (timing/energy not calibration-grade)"
                  : profile_reference) << std::endl;
    std::cout << std::endl;
}
