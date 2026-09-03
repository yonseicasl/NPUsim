#include <iomanip>
#include <iostream>

#include "kvcache.h"
#include "energy_units.h"

kvcache_invocation_t::kvcache_invocation_t() :
    active(false),
    read_bytes(0),
    dram_read_cycles(0.0),
    dram_read_energy(0.0),
    priced(true) {
}

kvcache_t::kvcache_t(section_config_t m_section_config) :
    index(0),
    kv_bytes_per_token(0),
    context_length(0),
    kv_compression_ratio(1.0),
    kv_decoder_bytes_per_cycle(0.0),
    u_read_energy(0.0),
    read_energy_is_declared(false) {

    init(m_section_config);
}

kvcache_t::~kvcache_t() {
}

void kvcache_t::init(section_config_t m_section_config) {
    m_section_config.get_setting("kv_bytes_per_token", &kv_bytes_per_token);
    m_section_config.get_setting("context_length", &context_length);
    m_section_config.get_setting("kv_compression_ratio", &kv_compression_ratio);
    m_section_config.get_setting("kv_decoder_bytes_per_cycle", &kv_decoder_bytes_per_cycle);

    // Fail-fast contract (evaluation.md: knob validity is stated, not defaulted). A KV read
    // that is present must describe a real read -- both the per-token footprint and the
    // context length have to be positive, or the section does nothing and silently reads as
    // "no memory pressure", which is exactly the failure this component exists to avoid.
    if(kv_bytes_per_token == 0) {
        std::cerr << "Error: [kvcache] kv_bytes_per_token must be positive" << std::endl;
        exit(1);
    }
    if(context_length == 0) {
        std::cerr << "Error: [kvcache] context_length must be positive" << std::endl;
        exit(1);
    }
    if(kv_compression_ratio < 1.0) {
        std::cerr << "Error: [kvcache] kv_compression_ratio must be >= 1"
                  << " (1 = no compression)" << std::endl;
        exit(1);
    }
    if(kv_decoder_bytes_per_cycle < 0.0) {
        std::cerr << "Error: [kvcache] kv_decoder_bytes_per_cycle must not be negative"
                  << std::endl;
        exit(1);
    }

    read_energy_is_declared = m_section_config.get_setting("kvcache_read_energy", &u_read_energy);
    m_section_config.get_setting("profile_reference", &profile_reference);
}

void kvcache_t::reset() {
    // The component carries no cross-layer state.
}

size_t kvcache_t::read_bytes() const {
    // Compatibility accessor: the bytes DRAM actually fetches (post-compression).
    return compressed_read_bytes();
}

size_t kvcache_t::dense_read_bytes() const {
    return kv_bytes_per_token*context_length;
}

bool kvcache_t::bypassed() const {
    return kv_compression_ratio <= 1.0;
}

size_t kvcache_t::compressed_read_bytes() const {
    const size_t dense = dense_read_bytes();
    if(bypassed()) return dense;
    // Never enlarge: a compressed read that would exceed dense stores/fetches dense.
    const double compressed = static_cast<double>(dense)/kv_compression_ratio;
    return compressed >= static_cast<double>(dense) ? dense
                                                    : static_cast<size_t>(compressed);
}

double kvcache_t::decoder_cycles(size_t m_dense_bytes) const {
    if(bypassed() || kv_decoder_bytes_per_cycle <= 0.0) return 0.0;
    return static_cast<double>(m_dense_bytes)/kv_decoder_bytes_per_cycle;
}

double kvcache_t::read_energy(size_t m_read_bytes, bool *m_priced) const {
    if(m_priced) *m_priced = read_energy_is_declared;
    return static_cast<double>(m_read_bytes)*u_read_energy;
}

void kvcache_t::print_specification() {
    std::cout << "========= KV-cache read =========" << std::endl;
    std::cout << "Per-token KV       :" << std::setw(18) << kv_bytes_per_token
              << " B/token/layer" << std::endl;
    std::cout << "Context length     :" << std::setw(18) << context_length
              << " tokens" << std::endl;
    std::cout << "Compression        :" << std::setw(18) << std::setprecision(2)
              << std::fixed << kv_compression_ratio << "x"
              << (bypassed() ? "  (uncompressed)" : "") << std::endl;
    std::cout << "Read per step      :" << std::setw(18) << std::setprecision(2)
              << static_cast<double>(dense_read_bytes())/1.0e6 << " MB dense / "
              << static_cast<double>(compressed_read_bytes())/1.0e6 << " MB fetched" << std::endl;
    std::cout << "Decoder throughput :" << std::setw(18) << std::setprecision(2)
              << kv_decoder_bytes_per_cycle << " dense B/cycle"
              << (kv_decoder_bytes_per_cycle <= 0.0 ? "  (ideal)" : "") << std::endl;
    std::cout << "Read energy        :" << std::setw(18) << std::setprecision(4)
              << u_read_energy << " /byte "
              << (read_energy_is_declared ? energy_units().label() : "(UNPRICED)")
              << std::endl;
    std::cout << "Provenance         : "
              << (profile_reference.empty()
                  ? "NOT DECLARED (traffic/energy not calibration-grade)"
                  : profile_reference) << std::endl;
    std::cout << std::endl;
}
