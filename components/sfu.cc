#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

#include "sfu.h"
#include "datatype.h"
#include "energy_units.h"

namespace {

size_t ceil_div(size_t m_value, size_t m_divisor) {
    return (m_value + m_divisor - 1)/m_divisor;
}

// Depth of a binary reduction tree over m_width live lanes.
size_t ceil_log2(size_t m_width) {
    size_t depth = 0;
    size_t reach = 1;
    while(reach < m_width) { reach *= 2; ++depth; }
    return depth;
}

}  // namespace

sfu_invocation_t::sfu_invocation_t() :
    operation(),
    valid_elements(0),
    operations(0),
    invocations(0),
    chunks(0),
    busy_cycle(0.0),
    tail_lane_utilization(0.0),
    ingress_elements(0),
    egress_elements(0),
    ingress_bytes(0),
    egress_bytes(0),
    ingress_transactions(0),
    egress_transactions(0),
    op_energy(0.0),
    read_energy(0.0),
    write_energy(0.0),
    setup_energy(0.0),
    timing_calibrated(true) {
}

void sfu_invocation_t::merge_parallel(const sfu_invocation_t &m_other) {
    // Chips execute their output partitions concurrently: the SFU window is the slowest
    // chip's window, while events, traffic and energy are physical work and sum.
    if(operation.empty()) operation = m_other.operation;
    valid_elements += m_other.valid_elements;
    operations += m_other.operations;
    invocations += m_other.invocations;
    chunks += m_other.chunks;
    busy_cycle = std::max(busy_cycle, m_other.busy_cycle);
    tail_lane_utilization = std::max(tail_lane_utilization, m_other.tail_lane_utilization);
    ingress_elements += m_other.ingress_elements;
    egress_elements += m_other.egress_elements;
    ingress_bytes += m_other.ingress_bytes;
    egress_bytes += m_other.egress_bytes;
    ingress_transactions += m_other.ingress_transactions;
    egress_transactions += m_other.egress_transactions;
    op_energy += m_other.op_energy;
    read_energy += m_other.read_energy;
    write_energy += m_other.write_energy;
    setup_energy += m_other.setup_energy;
    timing_calibrated = timing_calibrated && m_other.timing_calibrated;
    for(unsigned i = 0; i < m_other.unpriced.size(); ++i) {
        if(std::find(unpriced.begin(), unpriced.end(), m_other.unpriced[i]) == unpriced.end()) {
            unpriced.push_back(m_other.unpriced[i]);
        }
    }
}

void sfu_invocation_t::merge_serial(const sfu_invocation_t &m_other) {
    // Softmax passes execute back to back on the same pipelines: windows add.
    const double parallel_busy = busy_cycle;
    merge_parallel(m_other);
    busy_cycle = parallel_busy + m_other.busy_cycle;
}

sfu_t::sfu_t(section_config_t m_section_config) :
    index(0),
    num_units(1),
    lanes(1),
    queue_depth(1),
    setup_cycle(0.0),
    strict_profiles(false),
    placement("post_accumulator"),
    fusion("fused"),
    u_read_energy(0.0),
    u_write_energy(0.0),
    u_setup_energy(0.0),
    u_static_energy(0.0),
    read_energy_is_declared(false),
    write_energy_is_declared(false),
    setup_energy_is_declared(false),
    static_energy_is_declared(false) {

    init(m_section_config);
}

sfu_t::~sfu_t() {
}

void sfu_t::init(section_config_t m_section_config) {

    // Profile table. Name doubles as the config key stem (<name>_pipeline_latency,
    // <name>_initiation_interval, sfu_op_energy_<name>).
    static const struct { const char *name; bool bypass; bool micro_op; } op_table[NUM_SFU_OPS] = {
        { "linear",   true,  false },
        { "relu",     false, false },
        { "leaky",    false, false },
        { "elu",      false, false },
        { "sigmoid",  false, false },
        { "tanh",     false, false },
        { "hsigmoid", false, false },
        { "hswish",   false, false },
        { "gelu",     false, false },
        { "silu",     false, false },
        { "exp",      false, true  },
        { "recip",    false, true  },
        { "vmax",     false, true  },
        { "vadd",     false, true  },
        { "vmul",     false, true  },
    };
    for(unsigned i = 0; i < NUM_SFU_OPS; ++i) {
        profiles[i].name = op_table[i].name;
        profiles[i].bypass = op_table[i].bypass;
        profiles[i].micro_op = op_table[i].micro_op;
        profiles[i].pipeline_latency = 1.0;
        profiles[i].initiation_interval = 1.0;
        profiles[i].op_energy = 0.0;
        profiles[i].latency_declared = false;
        profiles[i].ii_declared = false;
        profiles[i].energy_declared = false;
    }

    m_section_config.get_setting("num_units_per_chip", &num_units);
    m_section_config.get_setting("lanes", &lanes);
    m_section_config.get_setting("queue_depth", &queue_depth);
    m_section_config.get_setting("setup_cycle", &setup_cycle);
    m_section_config.get_setting("strict_profiles", &strict_profiles);
    m_section_config.get_setting("placement", &placement);
    m_section_config.get_setting("fusion", &fusion);

    // Fail-fast contract checks (plan: 0 lane, 0 II, unsupported placement/fusion are
    // config errors, not silent defaults).
    if(num_units == 0) {
        std::cerr << "Error: [sfu] num_units_per_chip must be non-zero" << std::endl;
        exit(1);
    }
    if(lanes == 0) {
        std::cerr << "Error: [sfu] lanes must be non-zero" << std::endl;
        exit(1);
    }
    if(queue_depth != 1) {
        std::cerr << "Error: [sfu] queue_depth = " << queue_depth
                  << " is not supported yet; the initial SFU contract is the serial"
                  << " queue_depth = 1 model (streaming overlap is a later phase)" << std::endl;
        exit(1);
    }
    if(setup_cycle < 0.0) {
        std::cerr << "Error: [sfu] setup_cycle must not be negative" << std::endl;
        exit(1);
    }
    if(placement != "post_accumulator") {
        std::cerr << "Error: [sfu] placement = '" << placement
                  << "' is not supported; the initial contract is post_accumulator"
                  << " (activation before the output-format cast)" << std::endl;
        exit(1);
    }
    if(fusion != "fused") {
        std::cerr << "Error: [sfu] fusion = '" << fusion
                  << "' is not supported; the initial contract is fused (no intermediate"
                  << " tensor materialization)" << std::endl;
        exit(1);
    }

    // Per-operation profiles.
    for(unsigned i = 0; i < NUM_SFU_OPS; ++i) {
        sfu_op_profile_t &p = profiles[i];
        const std::string latency_key = std::string(p.name) + "_pipeline_latency";
        const std::string ii_key = std::string(p.name) + "_initiation_interval";
        const std::string energy_key = std::string("sfu_op_energy_") + p.name;
        double value = 0.0;
        if(m_section_config.get_setting(latency_key, &value)) {
            if(p.bypass) {
                std::cerr << "Error: [sfu] " << latency_key << " -- linear is a bypass"
                          << " operation and cannot carry a timing profile" << std::endl;
                exit(1);
            }
            if(value < 0.0) {
                std::cerr << "Error: [sfu] " << latency_key << " must not be negative" << std::endl;
                exit(1);
            }
            p.pipeline_latency = value;
            p.latency_declared = true;
        }
        if(m_section_config.get_setting(ii_key, &value)) {
            if(p.bypass) {
                std::cerr << "Error: [sfu] " << ii_key << " -- linear is a bypass"
                          << " operation and cannot carry a timing profile" << std::endl;
                exit(1);
            }
            if(value <= 0.0) {
                std::cerr << "Error: [sfu] " << ii_key
                          << " must be positive (a 0 initiation interval means infinite"
                          << " throughput)" << std::endl;
                exit(1);
            }
            p.initiation_interval = value;
            p.ii_declared = true;
        }
        if(m_section_config.get_setting(energy_key, &value)) {
            if(p.bypass) {
                std::cerr << "Error: [sfu] " << energy_key << " -- linear is a bypass"
                          << " operation and is modeled free by contract" << std::endl;
                exit(1);
            }
            // Negative/non-finite values are already rejected by the energy schema
            // validator; the assignment here only records the declared cost.
            p.op_energy = value;
            p.energy_declared = true;
        }
    }

    read_energy_is_declared = m_section_config.get_setting("sfu_read_energy", &u_read_energy);
    write_energy_is_declared = m_section_config.get_setting("sfu_write_energy", &u_write_energy);
    setup_energy_is_declared = m_section_config.get_setting("sfu_setup_energy", &u_setup_energy);
    static_energy_is_declared = m_section_config.get_setting("sfu_static_energy", &u_static_energy);
}

void sfu_t::reset() {
    // The SFU carries no cross-layer state: every invocation result is returned to the
    // caller and recorded in stats_t.
}

bool sfu_t::op_from_name(const std::string &m_name, sfu_op_t *m_op) {
    static const struct { const char *name; sfu_op_t op; } aliases[] = {
        { "linear",   SFU_OP_LINEAR },
        { "relu",     SFU_OP_RELU },
        { "leaky",    SFU_OP_LEAKY },
        { "elu",      SFU_OP_ELU },
        { "sigmoid",  SFU_OP_SIGMOID },
        { "logistic", SFU_OP_SIGMOID },   // Nebula's name for the same operation
        { "tanh",     SFU_OP_TANH },
        { "hsigmoid", SFU_OP_HSIGMOID },
        { "hswish",   SFU_OP_HSWISH },
        { "gelu",     SFU_OP_GELU },
        { "silu",     SFU_OP_SILU },
        { "swish",    SFU_OP_SILU },
    };
    for(unsigned i = 0; i < sizeof(aliases)/sizeof(aliases[0]); ++i) {
        if(m_name == aliases[i].name) {
            *m_op = aliases[i].op;
            return true;
        }
    }
    return false;
}

const char *sfu_t::op_name(sfu_op_t m_op) {
    static const char *names[NUM_SFU_OPS] = {
        "linear", "relu", "leaky", "elu", "sigmoid", "tanh", "hsigmoid", "hswish",
        "gelu", "silu", "exp", "recip", "vmax", "vadd", "vmul",
    };
    return names[m_op];
}

double sfu_t::invocation_cycles(const sfu_op_profile_t &m_profile, size_t m_chunks) const {
    if(m_chunks == 0) return 0.0;
    return setup_cycle + m_profile.pipeline_latency +
           static_cast<double>(m_chunks - 1)*m_profile.initiation_interval;
}

void sfu_t::require_profile(const sfu_op_profile_t &m_profile) const {
    if(!strict_profiles) return;
    if(!m_profile.latency_declared || !m_profile.ii_declared || !m_profile.energy_declared) {
        std::cerr << "Error: [sfu] strict_profiles = 1 and operation '" << m_profile.name
                  << "' is active without a complete profile; declare "
                  << m_profile.name << "_pipeline_latency, "
                  << m_profile.name << "_initiation_interval and sfu_op_energy_"
                  << m_profile.name << std::endl;
        exit(1);
    }
}

void sfu_t::note_profile_use(const sfu_op_profile_t &m_profile, size_t m_operations,
                             sfu_invocation_t *m_invocation) const {
    if(m_operations == 0) return;
    if(!m_profile.latency_declared || !m_profile.ii_declared) {
        m_invocation->timing_calibrated = false;
    }
    if(!m_profile.energy_declared) {
        m_invocation->unpriced.push_back(std::string("SFU ") + m_profile.name +
            " operation fired but [sfu] sfu_op_energy_" + m_profile.name +
            " is not declared");
    }
    m_invocation->operations += m_operations;
    m_invocation->op_energy += static_cast<double>(m_operations)*m_profile.op_energy;
}

void sfu_t::charge_traffic(size_t m_reads, size_t m_writes, size_t m_setups,
                           sfu_invocation_t *m_invocation) const {
    if(m_reads > 0) {
        m_invocation->ingress_elements += m_reads;
        m_invocation->ingress_bytes += runtime_datatypes().accumulator_storage_bytes(m_reads);
        m_invocation->ingress_transactions += ceil_div(m_reads, lanes);
        m_invocation->read_energy += static_cast<double>(m_reads)*u_read_energy;
        if(!read_energy_is_declared) {
            const std::string message = "SFU ingress read fired but [sfu] sfu_read_energy"
                                        " is not declared";
            if(std::find(m_invocation->unpriced.begin(), m_invocation->unpriced.end(),
                         message) == m_invocation->unpriced.end()) {
                m_invocation->unpriced.push_back(message);
            }
        }
    }
    if(m_writes > 0) {
        m_invocation->egress_elements += m_writes;
        m_invocation->egress_bytes += runtime_datatypes().accumulator_storage_bytes(m_writes);
        m_invocation->egress_transactions += ceil_div(m_writes, lanes);
        m_invocation->write_energy += static_cast<double>(m_writes)*u_write_energy;
        if(!write_energy_is_declared) {
            const std::string message = "SFU egress write fired but [sfu] sfu_write_energy"
                                        " is not declared";
            if(std::find(m_invocation->unpriced.begin(), m_invocation->unpriced.end(),
                         message) == m_invocation->unpriced.end()) {
                m_invocation->unpriced.push_back(message);
            }
        }
    }
    if(m_setups > 0) {
        m_invocation->invocations += m_setups;
        m_invocation->setup_energy += static_cast<double>(m_setups)*u_setup_energy;
        if(!setup_energy_is_declared) {
            const std::string message = "SFU pipeline setup fired but [sfu] sfu_setup_energy"
                                        " is not declared";
            if(std::find(m_invocation->unpriced.begin(), m_invocation->unpriced.end(),
                         message) == m_invocation->unpriced.end()) {
                m_invocation->unpriced.push_back(message);
            }
        }
    }
}

sfu_invocation_t sfu_t::elementwise_invocation(sfu_op_t m_op, size_t m_valid_elements) const {
    sfu_invocation_t invocation;
    const sfu_op_profile_t &profile = profiles[m_op];
    invocation.operation = profile.bypass ? std::string(profile.name) + " (bypass)"
                                          : std::string(profile.name);
    // Plan bypass/boundary contract: LINEAR and an empty output charge no element, no
    // cycle, no energy and no event.
    if(profile.bypass || m_valid_elements == 0) {
        return invocation;
    }
    require_profile(profile);

    invocation.valid_elements = m_valid_elements;
    // Output chunks are distributed evenly across the chip's parallel units; the chip's
    // SFU window is the busiest unit's window (plan cycle model).
    const size_t base = m_valid_elements/num_units;
    const size_t remainder = m_valid_elements % num_units;
    const size_t max_share = base + (remainder != 0 ? 1 : 0);
    size_t total_chunks = 0;
    for(unsigned u = 0; u < num_units; ++u) {
        const size_t share = base + (u < remainder ? 1 : 0);
        total_chunks += ceil_div(share, lanes);
    }
    invocation.chunks = total_chunks;
    invocation.busy_cycle = invocation_cycles(profile, ceil_div(max_share, lanes));
    // Tail utilization of the busiest unit's final chunk -- the elements were just split
    // across units, so the whole-share remainder would misstate multi-unit occupancy.
    const size_t tail = max_share % lanes;
    invocation.tail_lane_utilization = max_share == 0 ? 0.0
        : (tail == 0 ? 1.0 : static_cast<double>(tail)/static_cast<double>(lanes));

    const size_t active_units = base > 0 ? num_units : remainder;
    note_profile_use(profile, m_valid_elements, &invocation);
    // Fused post-accumulator contract: every valid element is read once from the final
    // accumulator value and written once toward the output cast, at ACCUMULATOR
    // precision -- the cast to the output format happens after the SFU.
    charge_traffic(m_valid_elements, m_valid_elements, active_units, &invocation);
    return invocation;
}

sfu_invocation_t sfu_t::softmax_invocation(size_t m_rows, size_t m_row_length) const {
    sfu_invocation_t invocation;
    invocation.operation = "softmax (max-sub-exp-sum-recip-mul multi-pass)";
    if(m_rows == 0 || m_row_length == 0) {
        return invocation;
    }
    const sfu_op_profile_t &vmax = profiles[SFU_OP_VMAX];
    const sfu_op_profile_t &vadd = profiles[SFU_OP_VADD];
    const sfu_op_profile_t &vexp = profiles[SFU_OP_EXP];
    const sfu_op_profile_t &recip = profiles[SFU_OP_RECIP];
    const sfu_op_profile_t &vmul = profiles[SFU_OP_VMUL];
    // The reduction passes only exist for rows longer than one element (a length-1 max
    // or sum performs zero operations, so it must charge zero cycles as well -- every
    // charged cycle needs an event source).
    const bool reduction_passes = m_row_length > 1;
    if(reduction_passes) require_profile(vmax);
    require_profile(vadd);
    require_profile(vexp);
    require_profile(recip);
    require_profile(vmul);

    invocation.valid_elements = m_rows*m_row_length;

    // Rows are distributed across the chip's parallel units and every unit runs the
    // passes serially over its own rows; the chip window follows the busiest unit.
    const size_t base_rows = m_rows/num_units;
    const size_t remainder = m_rows % num_units;
    const size_t max_rows = base_rows + (remainder != 0 ? 1 : 0);
    const size_t active_units = base_rows > 0 ? num_units : remainder;
    const size_t n = m_row_length;
    const size_t row_chunks = ceil_div(n, lanes);
    // Lane-tree fold tail per row reduction: after the streaming pass each row holds
    // min(lanes, n) partial results folding in ceil(log2()) steps. Each fold step is
    // data-dependent on the previous step's result, so a pipelined unit cannot overlap
    // them -- the serial chain costs the full PIPELINE LATENCY per step, not one II.
    const size_t fold_steps = ceil_log2(std::min<size_t>(lanes, n));
    // Tail utilization of the busiest unit's final element-wise chunk.
    const size_t tail = (max_rows*n) % lanes;
    invocation.tail_lane_utilization = max_rows == 0 ? 0.0
        : (tail == 0 ? 1.0 : static_cast<double>(tail)/static_cast<double>(lanes));

    struct pass_cost_t { double busy; size_t chunks; size_t reads; size_t writes; size_t ops; };
    const size_t total_rows = m_rows;
    const size_t elements_total = total_rows*n;
    pass_cost_t max_pass = {0.0, 0, 0, 0, 0};
    pass_cost_t sum_pass = {0.0, 0, 0, 0, 0};

    if(reduction_passes) {
        // Pass 1: max reduction -- every element streams through once per row, then the
        // per-row lane partials fold serially.
        max_pass.busy = invocation_cycles(vmax, max_rows*row_chunks) +
                        static_cast<double>(max_rows)*static_cast<double>(fold_steps)*
                        vmax.pipeline_latency;
        max_pass.chunks = total_rows*row_chunks;
        max_pass.reads = elements_total;
        max_pass.writes = total_rows;
        max_pass.ops = total_rows*(n - 1);

        // Pass 4: sum reduction (same shape, vadd profile).
        sum_pass.busy = invocation_cycles(vadd, max_rows*row_chunks) +
                        static_cast<double>(max_rows)*static_cast<double>(fold_steps)*
                        vadd.pipeline_latency;
        sum_pass.chunks = total_rows*row_chunks;
        sum_pass.reads = elements_total;
        sum_pass.writes = total_rows;
        sum_pass.ops = total_rows*(n - 1);
    }

    // Pass 2: subtract the row max (element-wise, the max value is re-read per row).
    pass_cost_t sub_pass;
    sub_pass.busy = invocation_cycles(vadd, ceil_div(max_rows*n, lanes));
    sub_pass.chunks = ceil_div(elements_total, lanes);
    sub_pass.reads = elements_total + total_rows;
    sub_pass.writes = elements_total;
    sub_pass.ops = elements_total;

    // Pass 3: exponential (element-wise).
    pass_cost_t exp_pass;
    exp_pass.busy = invocation_cycles(vexp, ceil_div(max_rows*n, lanes));
    exp_pass.chunks = ceil_div(elements_total, lanes);
    exp_pass.reads = elements_total;
    exp_pass.writes = elements_total;
    exp_pass.ops = elements_total;

    // Pass 5: reciprocal of the row sums.
    pass_cost_t recip_pass;
    recip_pass.busy = invocation_cycles(recip, ceil_div(max_rows, lanes));
    recip_pass.chunks = ceil_div(total_rows, lanes);
    recip_pass.reads = total_rows;
    recip_pass.writes = total_rows;
    recip_pass.ops = total_rows;

    // Pass 6: normalize (element-wise multiply by the row reciprocal, re-read per row).
    pass_cost_t mul_pass;
    mul_pass.busy = invocation_cycles(vmul, ceil_div(max_rows*n, lanes));
    mul_pass.chunks = ceil_div(elements_total, lanes);
    mul_pass.reads = elements_total + total_rows;
    mul_pass.writes = elements_total;
    mul_pass.ops = elements_total;

    const pass_cost_t passes[6] = {max_pass, sub_pass, exp_pass, sum_pass, recip_pass,
                                   mul_pass};
    size_t reads = 0, writes = 0, executed_passes = 0;
    for(unsigned p = 0; p < 6; ++p) {
        if(passes[p].ops == 0) continue;
        invocation.busy_cycle += passes[p].busy;
        invocation.chunks += passes[p].chunks;
        reads += passes[p].reads;
        writes += passes[p].writes;
        ++executed_passes;
    }

    note_profile_use(vmax, max_pass.ops, &invocation);
    note_profile_use(vadd, sub_pass.ops + sum_pass.ops, &invocation);
    note_profile_use(vexp, exp_pass.ops, &invocation);
    note_profile_use(recip, recip_pass.ops, &invocation);
    note_profile_use(vmul, mul_pass.ops, &invocation);

    charge_traffic(reads, writes, executed_passes*active_units, &invocation);
    return invocation;
}

float sfu_t::evaluate(sfu_op_t m_op, float m_x) {
    switch(m_op) {
        case SFU_OP_LINEAR:   return m_x;
        case SFU_OP_RELU:     return m_x > 0.0f ? m_x : 0.0f;
        case SFU_OP_LEAKY:    return m_x > 0.0f ? m_x : 0.1f*m_x;
        case SFU_OP_ELU:      return m_x >= 0.0f ? m_x : std::exp(m_x) - 1.0f;
        case SFU_OP_SIGMOID:  return 1.0f/(1.0f + std::exp(-m_x));
        case SFU_OP_TANH:     return std::tanh(m_x);
        case SFU_OP_HSIGMOID: return std::min(std::max(m_x + 3.0f, 0.0f), 6.0f)/6.0f;
        case SFU_OP_HSWISH:   return m_x*std::min(std::max(m_x + 3.0f, 0.0f), 6.0f)/6.0f;
        case SFU_OP_GELU:
            // tanh approximation (Hendrycks & Gimpel).
            return 0.5f*m_x*(1.0f + std::tanh(0.7978845608f*(m_x + 0.044715f*m_x*m_x*m_x)));
        case SFU_OP_SILU:     return m_x/(1.0f + std::exp(-m_x));
        case SFU_OP_EXP:      return std::exp(m_x);
        case SFU_OP_RECIP:    return 1.0f/m_x;
        default:              return m_x;
    }
}

void sfu_t::evaluate_softmax(const float *m_input, float *m_output, size_t m_size) {
    if(m_size == 0) return;
    float max_value = m_input[0];
    for(size_t i = 1; i < m_size; ++i) max_value = std::max(max_value, m_input[i]);
    float sum = 0.0f;
    for(size_t i = 0; i < m_size; ++i) {
        m_output[i] = std::exp(m_input[i] - max_value);
        sum += m_output[i];
    }
    const float reciprocal = 1.0f/sum;
    for(size_t i = 0; i < m_size; ++i) m_output[i] *= reciprocal;
}

void sfu_t::print_specification() {
    std::cout << "============ SFU specification =============" << std::endl;
    std::cout << "Units per chip     :" << std::setw(24) << num_units << std::endl;
    std::cout << "Lanes per unit     :" << std::setw(24) << lanes << std::endl;
    std::cout << "Queue depth        :" << std::setw(24) << queue_depth
              << " (serial contract)" << std::endl;
    std::cout << "Placement          :" << std::setw(24) << placement << std::endl;
    std::cout << "Fusion             :" << std::setw(24) << fusion << std::endl;
    std::cout << "Setup cycle        :" << std::setw(24) << std::setprecision(1)
              << std::fixed << setup_cycle << std::endl;
    std::cout << "Clock              :" << std::setw(24) << "chip clock"
              << " (post-accumulator on-chip datapath)" << std::endl;
    std::cout << "Access energy (in/out/setup/static)" << std::endl;
    std::cout << " * per element     :" << std::setw(16) << std::setprecision(3)
              << u_read_energy << "/" << u_write_energy << "/" << u_setup_energy
              << "/" << u_static_energy << " " << energy_units().label() << std::endl;
    std::cout << "Operation profiles (latency/II/energy; '-' = not declared, defaults 1/1/unpriced)"
              << std::endl;
    for(unsigned i = 0; i < NUM_SFU_OPS; ++i) {
        const sfu_op_profile_t &p = profiles[i];
        if(p.bypass) {
            std::cout << " * " << std::left << std::setw(9) << p.name << std::right
                      << ": bypass (free by contract)" << std::endl;
            continue;
        }
        if(!p.latency_declared && !p.ii_declared && !p.energy_declared) continue;
        auto format_value = [](double value) {
            char buffer[32];
            std::snprintf(buffer, sizeof(buffer), "%.6g", value);
            return std::string(buffer);
        };
        std::cout << " * " << std::left << std::setw(9) << p.name << std::right << ": "
                  << (p.latency_declared ? format_value(p.pipeline_latency) : std::string("-"))
                  << "/"
                  << (p.ii_declared ? format_value(p.initiation_interval) : std::string("-"))
                  << "/"
                  << (p.energy_declared ? format_value(p.op_energy) : std::string("-"))
                  << (p.micro_op ? "  [softmax micro-op]" : "") << std::endl;
    }
    std::cout << std::endl;
}
