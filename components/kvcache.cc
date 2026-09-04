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
    kv_schedule("aggregate"),
    kv_tile_bytes(0),
    kv_buffer_tiles(0),
    attention(false),
    n_q_heads(0),
    n_kv_heads(0),
    head_dim(0),
    kv_precision_bytes(1),
    attention_macs_per_cycle(0.0),
    softmax_cycles_per_element(0.0),
    attention_algorithm("online"),
    kv_capacity_bytes(0),
    u_read_energy(0.0),
    read_energy_is_declared(false),
    u_attention_mac_energy(0.0),
    attention_mac_energy_is_declared(false),
    u_softmax_energy(0.0),
    softmax_energy_is_declared(false) {

    init(m_section_config);
}

kvcache_t::~kvcache_t() {
}

void kvcache_t::init(section_config_t m_section_config) {
    m_section_config.get_setting("kv_bytes_per_token", &kv_bytes_per_token);
    m_section_config.get_setting("context_length", &context_length);
    m_section_config.get_setting("kv_compression_ratio", &kv_compression_ratio);
    m_section_config.get_setting("kv_decoder_bytes_per_cycle", &kv_decoder_bytes_per_cycle);
    m_section_config.get_setting("kv_schedule", &kv_schedule);
    m_section_config.get_setting("kv_tile_bytes", &kv_tile_bytes);
    m_section_config.get_setting("kv_buffer_tiles", &kv_buffer_tiles);
    m_section_config.get_setting("attention", &attention);
    m_section_config.get_setting("n_q_heads", &n_q_heads);
    m_section_config.get_setting("n_kv_heads", &n_kv_heads);
    m_section_config.get_setting("head_dim", &head_dim);
    m_section_config.get_setting("kv_precision_bytes", &kv_precision_bytes);
    m_section_config.get_setting("attention_macs_per_cycle", &attention_macs_per_cycle);
    m_section_config.get_setting("softmax_cycles_per_element", &softmax_cycles_per_element);
    m_section_config.get_setting("attention_algorithm", &attention_algorithm);
    m_section_config.get_setting("kv_capacity_bytes", &kv_capacity_bytes);

    // Attention consumer mode: the geometry defines the per-token KV footprint, so derive
    // it -- and refuse a contradictory explicit declaration instead of silently picking one.
    if(attention) {
        if(n_q_heads == 0 || n_kv_heads == 0 || head_dim == 0) {
            std::cerr << "Error: [kvcache] attention = 1 requires positive n_q_heads,"
                      << " n_kv_heads and head_dim" << std::endl;
            exit(1);
        }
        if(n_q_heads % n_kv_heads != 0) {
            std::cerr << "Error: [kvcache] n_q_heads (" << n_q_heads << ") must be a"
                      << " multiple of n_kv_heads (" << n_kv_heads << ") -- GQA grouping"
                      << std::endl;
            exit(1);
        }
        if(kv_precision_bytes == 0) {
            std::cerr << "Error: [kvcache] kv_precision_bytes must be positive" << std::endl;
            exit(1);
        }
        if(attention_algorithm != "online" && attention_algorithm != "two_pass") {
            std::cerr << "Error: [kvcache] attention_algorithm must be online or two_pass"
                      << " (got '" << attention_algorithm << "')" << std::endl;
            exit(1);
        }
        // The attention consumer IS the schedule's consumer; folding the KV read into the
        // mapped layer's DRAM axis (aggregate) would leave the consumer with nothing.
        if(kv_schedule == "aggregate") {
            std::cerr << "Error: [kvcache] attention = 1 requires kv_schedule = blocking,"
                      << " streaming, or double_buffered (aggregate has no consumer)"
                      << std::endl;
            exit(1);
        }
        const size_t derived = 2ULL*n_kv_heads*head_dim*kv_precision_bytes;
        if(kv_bytes_per_token != 0 && kv_bytes_per_token != derived) {
            std::cerr << "Error: [kvcache] kv_bytes_per_token = " << kv_bytes_per_token
                      << " contradicts the attention geometry (2 * " << n_kv_heads
                      << " * " << head_dim << " * " << kv_precision_bytes << " = "
                      << derived << "); drop one of the two" << std::endl;
            exit(1);
        }
        kv_bytes_per_token = derived;
    }

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
    // Fail-fast schedule contract: a typo'd schedule silently falling back to aggregate
    // would read as "streaming made no difference" -- the failure this section exists to
    // avoid. streaming needs a tile granularity to pipeline over; double_buffered defaults
    // its buffer to 2 tiles (that IS the double buffer), streaming to 1.
    if(kv_schedule != "aggregate" && kv_schedule != "blocking" &&
       kv_schedule != "streaming" && kv_schedule != "double_buffered") {
        std::cerr << "Error: [kvcache] kv_schedule must be aggregate, blocking, streaming,"
                  << " or double_buffered (got '" << kv_schedule << "')" << std::endl;
        exit(1);
    }
    if((kv_schedule == "streaming" || kv_schedule == "double_buffered") && kv_tile_bytes == 0) {
        std::cerr << "Error: [kvcache] kv_schedule = " << kv_schedule
                  << " requires a positive kv_tile_bytes (the KV supply tile)" << std::endl;
        exit(1);
    }
    if(kv_buffer_tiles == 0) {
        // Both overlapped schedules need >= 2 tile slots to actually stream: with a single
        // buffer the consumer must drain a tile before the next can land, which degenerates
        // to a serial handoff (the pipeline model's depth-1 semantics -- physically, one
        // shared buffer cannot be filled and drained at once). An explicit kv_buffer_tiles
        // = 1 is honored as exactly that single-buffer ablation.
        kv_buffer_tiles = (kv_schedule == "streaming" || kv_schedule == "double_buffered")
                          ? 2 : 1;
    }
    // Cache state: this decode step appends the current token's K/V, so the occupancy after
    // the step is (context + 1) tokens. A declared capacity it does not fit is a config
    // error, not a silent overflow.
    if(kv_capacity_bytes > 0) {
        const size_t occupancy = kv_bytes_per_token*(context_length + 1);
        if(occupancy > kv_capacity_bytes) {
            std::cerr << "Error: [kvcache] cache occupancy after this step ("
                      << occupancy << " B = " << kv_bytes_per_token << " B/token x "
                      << (context_length + 1) << " tokens) exceeds kv_capacity_bytes = "
                      << kv_capacity_bytes << std::endl;
            exit(1);
        }
    }

    read_energy_is_declared = m_section_config.get_setting("kvcache_read_energy", &u_read_energy);
    attention_mac_energy_is_declared =
        m_section_config.get_setting("kvcache_attention_mac_energy", &u_attention_mac_energy);
    softmax_energy_is_declared =
        m_section_config.get_setting("kvcache_softmax_energy", &u_softmax_energy);
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
    if(attention) {
        std::cout << "==== KV cache + attention step ====" << std::endl;
        std::cout << "Model scope        : attention consumer (QK/softmax/AV, "
                  << attention_algorithm << "), dedicated KV DRAM stream, cache append;"
                  << " one decode step at the declared context (no multi-step growth)"
                  << std::endl;
        std::cout << "Geometry           :" << std::setw(11) << n_q_heads << " q-heads / "
                  << n_kv_heads << " kv-heads / head_dim " << head_dim << " / "
                  << kv_precision_bytes << " B/elem" << std::endl;
        std::cout << "Attention rate     :" << std::setw(18) << std::setprecision(2)
                  << std::fixed << attention_macs_per_cycle << " MAC/cycle"
                  << (attention_macs_per_cycle <= 0.0 ? "  (ideal: UNCALIBRATED)" : "")
                  << std::endl;
        std::cout << "Softmax rate       :" << std::setw(18) << std::setprecision(2)
                  << softmax_cycles_per_element << " cycles/score"
                  << (softmax_cycles_per_element <= 0.0 ? "  (ideal: UNCALIBRATED)" : "")
                  << std::endl;
    } else {
        std::cout << "==== KV-cache read (PROXY) ====" << std::endl;
        std::cout << "Model scope        : traffic/supply-scheduling proxy -- no attention"
                  << " consumer, no KV address stream, no KV write/cache state" << std::endl;
    }
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
    std::cout << "Schedule           :" << std::setw(18) << kv_schedule;
    if(kv_schedule == "streaming" || kv_schedule == "double_buffered") {
        std::cout << "  (tile " << kv_tile_bytes << " B, buffer " << kv_buffer_tiles
                  << " tiles)";
    }
    std::cout << std::endl;
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
