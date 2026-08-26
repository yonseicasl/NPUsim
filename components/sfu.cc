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
    softmax_operand_residency("dram"),
    precision_mismatch(false),
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
    // <name>_initiation_interval, sfu_op_energy_<name>, <name>_approximation).
    // Phase-6: `transcendental` marks LUT/polynomial hardware (vs piecewise ALU ops)
    // and decides which approximation declarations are legal.
    static const struct {
        const char *name; bool bypass; bool micro_op; bool transcendental;
    } op_table[NUM_SFU_OPS] = {
        { "linear",   true,  false, false },
        { "relu",     false, false, false },
        { "leaky",    false, false, false },
        { "elu",      false, false, true  },
        { "sigmoid",  false, false, true  },
        { "tanh",     false, false, true  },
        { "hsigmoid", false, false, false },
        { "hswish",   false, false, false },
        { "gelu",     false, false, true  },
        { "silu",     false, false, true  },
        { "loggy",    false, false, true  },
        { "exp",      false, true,  true  },
        { "recip",    false, true,  true  },
        { "vmax",     false, true,  false },
        { "vadd",     false, true,  false },
        { "vmul",     false, true,  false },
    };
    for(unsigned i = 0; i < NUM_SFU_OPS; ++i) {
        profiles[i].name = op_table[i].name;
        profiles[i].bypass = op_table[i].bypass;
        profiles[i].micro_op = op_table[i].micro_op;
        profiles[i].transcendental = op_table[i].transcendental;
        profiles[i].pipeline_latency = 1.0;
        profiles[i].initiation_interval = 1.0;
        profiles[i].op_energy = 0.0;
        profiles[i].latency_declared = false;
        profiles[i].ii_declared = false;
        profiles[i].energy_declared = false;
        profiles[i].approximation.clear();
    }

    m_section_config.get_setting("num_units_per_chip", &num_units);
    m_section_config.get_setting("lanes", &lanes);
    m_section_config.get_setting("queue_depth", &queue_depth);
    m_section_config.get_setting("setup_cycle", &setup_cycle);
    m_section_config.get_setting("strict_profiles", &strict_profiles);
    m_section_config.get_setting("placement", &placement);
    m_section_config.get_setting("fusion", &fusion);
    m_section_config.get_setting("softmax_operand_residency", &softmax_operand_residency);
    m_section_config.get_setting("profile_reference", &profile_reference);

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
    // Phase-5 streaming contract: depth 1 keeps the serial model, depth >= 2 lets the
    // producer run ahead by depth-1 output tiles. Depth 0 would mean no staging register
    // at all, which cannot transfer a tile.
    if(queue_depth == 0) {
        std::cerr << "Error: [sfu] queue_depth must be at least 1 (1 = serial contract,"
                  << " >= 2 = streaming overlap through a bounded output-tile queue)" << std::endl;
        exit(1);
    }
    if(softmax_operand_residency != "dram" && softmax_operand_residency != "glb") {
        std::cerr << "Error: [sfu] softmax_operand_residency = '" << softmax_operand_residency
                  << "' is not supported; use 'dram' (materialized round trip) or 'glb'"
                  << " (on-chip retained, must fit the global buffer)" << std::endl;
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
        // Phase-6: approximation mode. ALU-class ops are piecewise/exact by nature;
        // transcendental ops are implemented as a LUT, a polynomial expansion, or a
        // piecewise-linear approximation. Anything outside the op's legal set -- and any
        // declaration on the linear bypass -- fails fast (plan Phase-6 gate).
        const std::string approximation_key = std::string(p.name) + "_approximation";
        std::string approximation;
        if(m_section_config.get_setting(approximation_key, &approximation)) {
            if(p.bypass) {
                std::cerr << "Error: [sfu] " << approximation_key << " -- linear is a"
                          << " bypass and implements no approximation" << std::endl;
                exit(1);
            }
            const bool legal = p.transcendental
                ? (approximation == "lut" || approximation == "polynomial" ||
                   approximation == "piecewise")
                : (approximation == "exact" || approximation == "piecewise");
            if(!legal) {
                std::cerr << "Error: [sfu] " << approximation_key << " = '" << approximation
                          << "' is not supported for '" << p.name << "' ("
                          << (p.transcendental ? "transcendental: lut, polynomial, piecewise"
                                               : "ALU-class: exact, piecewise")
                          << ")" << std::endl;
                exit(1);
            }
            p.approximation = approximation;
        }
    }

    read_energy_is_declared = m_section_config.get_setting("sfu_read_energy", &u_read_energy);
    write_energy_is_declared = m_section_config.get_setting("sfu_write_energy", &u_write_energy);
    setup_energy_is_declared = m_section_config.get_setting("sfu_setup_energy", &u_setup_energy);
    static_energy_is_declared = m_section_config.get_setting("sfu_static_energy", &u_static_energy);

    // Phase-6: precision qualification. The SFU processes the final ACCUMULATOR value
    // (post-accumulator, pre-cast), so a declared profile precision is checked against
    // the runtime accumulator format. parse_data_format() fails fast on an unknown name
    // and canonicalizes aliases (float32 -> fp32) so the comparison is on canonical
    // names. strict_profiles turns a mismatch into a config error; otherwise every
    // invocation is marked UNCALIBRATED and the note reaches the layer report.
    if(m_section_config.get_setting("profile_precision", &profile_precision)) {
        const tensor_format_t declared = parse_data_format(profile_precision);
        profile_precision = declared.name;
        const std::string &runtime = runtime_datatypes().accumulator_format().name;
        if(profile_precision != runtime) {
            if(strict_profiles) {
                std::cerr << "Error: [sfu] profile_precision = " << profile_precision
                          << " but the runtime accumulator format is " << runtime
                          << "; strict_profiles requires a profile characterized for the"
                          << " precision actually running" << std::endl;
                exit(1);
            }
            precision_mismatch = true;
            precision_note = "SFU profile characterized for " + profile_precision +
                             " but the runtime accumulator is " + runtime +
                             " -- timing UNCALIBRATED for this precision";
        }
    }
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
        { "loggy",    SFU_OP_LOGGY },
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
        "gelu", "silu", "loggy", "exp", "recip", "vmax", "vadd", "vmul",
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
    // Phase-6: under strict profiles a transcendental op must also state HOW it is
    // implemented -- a latency without its approximation mode is half a profile.
    if(m_profile.transcendental && m_profile.approximation.empty()) {
        std::cerr << "Error: [sfu] strict_profiles = 1 and transcendental operation '"
                  << m_profile.name << "' is active without a declared "
                  << m_profile.name << "_approximation (lut, polynomial or piecewise)"
                  << std::endl;
        exit(1);
    }
}

void sfu_t::note_profile_use(const sfu_op_profile_t &m_profile, size_t m_operations,
                             sfu_invocation_t *m_invocation) const {
    if(m_operations == 0) return;
    if(!m_profile.latency_declared || !m_profile.ii_declared) {
        m_invocation->timing_calibrated = false;
    }
    // Phase-6: a profile characterized for a different precision than the runtime
    // accumulator cannot ground a calibrated cycle claim.
    if(precision_mismatch) {
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

sfu_invocation_t sfu_t::elementwise_invocation(sfu_op_t m_op, size_t m_valid_elements,
                                               size_t m_commit_events) const {
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
    // Phase-2: the valid elements arrive in m_commit_events final_output_tile events
    // (evenly split -- at most two event sizes), each a serial SFU invocation with its
    // own setup and pipeline fill. Within an event the chip's parallel units split the
    // event's elements; the event's window follows the busiest unit.
    const size_t events = std::max<size_t>(1, std::min(m_commit_events, m_valid_elements));
    const size_t event_base = m_valid_elements/events;
    const size_t event_remainder = m_valid_elements % events;
    double worst_tail_utilization = 1.0;
    size_t setups = 0;
    // At most two distinct event sizes: `count` events of `elements` each.
    const struct { size_t elements; size_t count; } event_shapes[2] = {
        { event_base + 1, event_remainder },
        { event_base,     events - event_remainder },
    };
    for(unsigned s = 0; s < 2; ++s) {
        const size_t elements = event_shapes[s].elements;
        const size_t count = event_shapes[s].count;
        if(count == 0 || elements == 0) continue;
        const size_t unit_base = elements/num_units;
        const size_t unit_remainder = elements % num_units;
        const size_t max_share = unit_base + (unit_remainder != 0 ? 1 : 0);
        size_t event_chunks = 0;
        for(unsigned u = 0; u < num_units; ++u) {
            const size_t share = unit_base + (u < unit_remainder ? 1 : 0);
            event_chunks += ceil_div(share, lanes);
        }
        invocation.chunks += count*event_chunks;
        invocation.busy_cycle += static_cast<double>(count)*
                                 invocation_cycles(profile, ceil_div(max_share, lanes));
        const size_t active_units = unit_base > 0 ? num_units : unit_remainder;
        setups += count*active_units;
        // Tail utilization of the busiest unit's final chunk; report the WORST event
        // shape so per-tile rounding loss stays visible.
        const size_t tail = max_share % lanes;
        const double utilization = tail == 0 ? 1.0
            : static_cast<double>(tail)/static_cast<double>(lanes);
        worst_tail_utilization = std::min(worst_tail_utilization, utilization);
    }
    invocation.tail_lane_utilization = worst_tail_utilization;

    note_profile_use(profile, m_valid_elements, &invocation);
    // Fused post-accumulator contract: every valid element is read once from the final
    // accumulator value and written once toward the output cast, at ACCUMULATOR
    // precision -- the cast to the output format happens after the SFU.
    charge_traffic(m_valid_elements, m_valid_elements, setups, &invocation);
    // Internal transactions are lane-wide transfers, one per issued chunk: per-event
    // tail rounding makes this >= ceil(elements/lanes), so mirror the chunk count.
    invocation.ingress_transactions = invocation.chunks;
    invocation.egress_transactions = invocation.chunks;
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
        case SFU_OP_LOGGY:    return 2.0f/(1.0f + std::exp(-m_x)) - 1.0f;
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
              << (queue_depth == 1 ? " (serial contract)"
                                   : " (streaming: producer may run ahead)") << std::endl;
    std::cout << "Softmax operand    :" << std::setw(24) << softmax_operand_residency
              << (softmax_operand_residency == "dram" ? " (materialized round trip)"
                                                      : " (on-chip retained)") << std::endl;
    std::cout << "Profile provenance : "
              << (profile_reference.empty()
                  ? "NOT DECLARED (timing numbers are not calibration-grade)"
                  : profile_reference) << std::endl;
    std::cout << "Profile precision  : "
              << (profile_precision.empty() ? "not declared" : profile_precision)
              << (precision_mismatch ? "  [MISMATCH vs runtime accumulator -- UNCALIBRATED]"
                                     : "") << std::endl;
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
                  << " [" << (p.approximation.empty() ? "approx?" : p.approximation) << "]"
                  << (p.micro_op ? "  [softmax micro-op]" : "") << std::endl;
    }
    std::cout << std::endl;
}
