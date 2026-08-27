#include "energy_units.h"

#include <algorithm>
#include <iostream>
#include <vector>

#include "datatype.h"

energy_units_t::energy_units_t() : unit(energy_unit_t::UNSPECIFIED), reference() {}

void energy_units_t::configure(section_config_t &m_section_config) {
    std::string unit_str;
    if(m_section_config.get_setting("energy_unit", &unit_str)) {
        for(unsigned i = 0; i < unit_str.size(); ++i) {
            unit_str[i] = static_cast<char>(tolower(unit_str[i]));
        }
        if(unit_str == "pj" || unit_str == "picojoule") {
            unit = energy_unit_t::PICOJOULE;
        } else if(unit_str == "normalized" || unit_str == "relative") {
            unit = energy_unit_t::NORMALIZED;
        } else {
            std::cerr << "Error: energy_unit must be pJ or normalized (got '" << unit_str
                      << "')" << std::endl;
            exit(1);
        }
    }
    m_section_config.get_setting("energy_reference", &reference);
}

const char *energy_units_t::label() const {
    switch(unit) {
        case energy_unit_t::PICOJOULE:  return "pJ";
        case energy_unit_t::NORMALIZED: return "normalized";
        default:                        return "uncalibrated";
    }
}

bool energy_units_t::is_absolute() const {
    // RE2: BOTH conditions. A declared pJ unit with no provenance is an unverifiable claim, and an
    // undeclared unit is not a claim at all -- neither earns an absolute total or a wattage.
    return unit == energy_unit_t::PICOJOULE && !reference.empty();
}

std::string energy_units_t::calibration_note() const {
    switch(unit) {
        case energy_unit_t::PICOJOULE:
            return reference.empty()
                ? "energy_unit = pJ is declared but energy_reference is not, so the absolute scale"
                  " is an unverifiable claim"
                : "declared picojoule scale with provenance; absolute totals and power additionally"
                  " require complete pricing of every active event";
        case energy_unit_t::NORMALIZED:
            return "relative to the declared reference, so absolute totals and watts are not"
                   " meaningful";
        default:
            return "no energy_unit declared, so the numbers are relative event x unit-cost sums of"
                   " unknown scale";
    }
}

std::string energy_units_t::describe() const {
    std::string text(label());
    // E20-1: this description is printed before a layer runs, so it cannot know whether an active
    // setup, accumulator, cast, format, reduction or row event will be unpriced. Calling totals and
    // power meaningful here contradicted the later UNDERCOUNT gate in exactly those runs. State the
    // unit-level qualification only; stats_t owns the run-level decision.
    text += is_absolute() ? " (declared absolute-unit candidate; run-level completeness required)"
                          : " (ESTIMATED scale; absolute totals and power are NOT meaningful)";
    text += " -- ";
    text += calibration_note();
    if(!reference.empty()) {
        text += "; provenance: " + reference;
    }
    return text;
}

energy_units_t &energy_units() {
    static energy_units_t instance;
    return instance;
}

// RE5: the declared energy schema. A key containing "energy" that is not here is a typo, not a
// component with no cost -- the difference used to be invisible.
//
// This table was first built from a survey of the keys the SHIPPED CONFIGS declared, which left out
// every energy key the code reads that no config had ever set -- format_payload_energy,
// format_metadata_energy, lb_static_energy, adder_energy and the nop_energy alias. The omission
// only surfaced when a new fixture tried to declare one and was rejected as a typo. The schema must
// be complete against what the CODE READS, not against what the configs happen to use;
// validation/knobs KN10 now asserts exactly that.
static const energy_key_schema_t g_energy_keys[] = {
    { "energy_unit",               ENERGY_KEY_PROVENANCE },
    { "energy_reference",          ENERGY_KEY_PROVENANCE },
    // per-datatype costs: exactly one value per datatype
    { "read_energy",               ENERGY_KEY_DATATYPE_VECTOR },
    { "write_energy",              ENERGY_KEY_DATATYPE_VECTOR },
    { "static_energy",             ENERGY_KEY_DATATYPE_VECTOR },
    { "mac_read_energy",           ENERGY_KEY_DATATYPE_VECTOR },
    { "mac_write_energy",          ENERGY_KEY_DATATYPE_VECTOR },
    { "lb_read_energy",            ENERGY_KEY_DATATYPE_VECTOR },
    { "lb_write_energy",           ENERGY_KEY_DATATYPE_VECTOR },
    { "pe_array_static_energy",    ENERGY_KEY_DATATYPE_VECTOR },
    { "output_cast_energy",        ENERGY_KEY_DATATYPE_VECTOR },
    { "accumulator_spill_energy",  ENERGY_KEY_DATATYPE_VECTOR },
    { "format_payload_energy",     ENERGY_KEY_DATATYPE_VECTOR },
    { "format_metadata_energy",    ENERGY_KEY_DATATYPE_VECTOR },
    { "lb_static_energy",          ENERGY_KEY_DATATYPE_VECTOR },
    // scalar costs
    { "transfer_energy",           ENERGY_KEY_SCALAR },
    { "transfer_energy_pe",        ENERGY_KEY_SCALAR },
    { "noc_energy",                ENERGY_KEY_SCALAR },
    { "pe_array_read_energy",      ENERGY_KEY_SCALAR },
    { "pe_array_write_energy",     ENERGY_KEY_SCALAR },
    { "computation_energy",        ENERGY_KEY_SCALAR },
    { "weight_fold_fill_energy",   ENERGY_KEY_SCALAR },
    { "layer_setup_energy",        ENERGY_KEY_SCALAR },
    { "stripe_transition_energy",  ENERGY_KEY_SCALAR },
    { "mac_reduction_energy",      ENERGY_KEY_SCALAR },
    { "row_miss_energy",           ENERGY_KEY_SCALAR },
    { "adder_energy",              ENERGY_KEY_SCALAR },
    // multi_chip.cc accepts `nop_energy` as an alias for `noc_energy`, so the schema must know it
    // or a config using the documented alias would be rejected.
    { "nop_energy",                ENERGY_KEY_SCALAR },
    // SFU (plan/plan_sfu.md): internal element access, per-invocation setup and pJ/cycle
    // leakage of one physical SFU unit. Per-operation dynamic costs form a prefix family
    // (sfu_op_energy_relu, sfu_op_energy_exp, ...) so ReLU and softmax micro-ops are never
    // priced by one shared scalar.
    { "sfu_read_energy",           ENERGY_KEY_SCALAR },
    { "sfu_write_energy",          ENERGY_KEY_SCALAR },
    { "sfu_setup_energy",          ENERGY_KEY_SCALAR },
    { "sfu_static_energy",         ENERGY_KEY_SCALAR },
    { "sfu_op_energy_",            ENERGY_KEY_PREFIX_SCALAR },
    // RE4/precision family: mac_energy_<input>_<weight>
    { "mac_energy_",               ENERGY_KEY_PREFIX_SCALAR },
};

const energy_key_schema_t *energy_key_schema(unsigned &count) {
    count = sizeof(g_energy_keys)/sizeof(g_energy_keys[0]);
    return g_energy_keys;
}

// Edit distance, capped -- only used to name the nearest known key in a typo message.
static unsigned key_distance(const std::string &a, const std::string &b) {
    std::vector<unsigned> previous(b.size() + 1), current(b.size() + 1);
    for(size_t j = 0; j <= b.size(); ++j) previous[j] = static_cast<unsigned>(j);
    for(size_t i = 1; i <= a.size(); ++i) {
        current[0] = static_cast<unsigned>(i);
        for(size_t j = 1; j <= b.size(); ++j) {
            const unsigned substitute = previous[j-1] + (a[i-1] == b[j-1] ? 0u : 1u);
            current[j] = std::min(substitute, std::min(previous[j] + 1u, current[j-1] + 1u));
        }
        previous = current;
    }
    return previous[b.size()];
}

// Returns the schema entry for a key, or 0 when the key is not declared.
static const energy_key_schema_t *lookup_energy_key(const std::string &name) {
    unsigned count = 0;
    const energy_key_schema_t *keys = energy_key_schema(count);
    for(unsigned i = 0; i < count; ++i) {
        if(keys[i].kind == ENERGY_KEY_PREFIX_SCALAR) {
            const std::string prefix(keys[i].name);
            if(name.size() > prefix.size() && name.compare(0, prefix.size(), prefix) == 0) {
                return &keys[i];
            }
        } else if(name == keys[i].name) {
            return &keys[i];
        }
    }
    return 0;
}

static std::string nearest_energy_key(const std::string &name) {
    unsigned count = 0;
    const energy_key_schema_t *keys = energy_key_schema(count);
    std::string best;
    unsigned best_distance = 4;   // beyond this a "did you mean" is noise, not help
    for(unsigned i = 0; i < count; ++i) {
        const unsigned distance = key_distance(name, keys[i].name);
        if(distance < best_distance) { best_distance = distance; best = keys[i].name; }
    }
    return best;
}

std::string validate_energy_settings(config_t &m_config) {
    for(unsigned i = 0; i < m_config.sections.size(); ++i) {
        section_config_t &section = m_config.sections[i];
        const std::map<std::string, std::string> &settings = section.all_settings();
        for(std::map<std::string, std::string>::const_iterator it = settings.begin();
            it != settings.end(); ++it) {
            if(it->first.find("energy") == std::string::npos) continue;
            // RE5: an undeclared energy key is a typo. Left unchecked it leaves the component it
            // was meant to price at zero cost, which the report cannot distinguish from a
            // deliberate modeled zero.
            const energy_key_schema_t *schema = lookup_energy_key(it->first);
            if(schema == 0) {
                const std::string suggestion = nearest_energy_key(it->first);
                return "[" + section.name + "] " + it->first + " is not a declared energy key" +
                       (suggestion.empty() ? std::string()
                                           : " (did you mean '" + suggestion + "'?)") +
                       "; an unrecognized energy key silently leaves its component at zero cost";
            }
            // energy_unit / energy_reference are provenance strings, not costs.
            if(schema->kind == ENERGY_KEY_PROVENANCE) continue;
            // A cost is either a scalar or a colon-separated per-datatype vector.
            std::stringstream fields(it->second);
            std::string field;
            bool parsed_any = false;
            unsigned field_count = 0;
            while(std::getline(fields, field, ':')) {
                ++field_count;
                // RE5: an empty field is not "skip me" -- it is a missing per-datatype cost that
                // used to silently default that datatype to zero.
                if(field.empty()) {
                    return "[" + section.name + "] " + it->first + " = '" + it->second +
                           "' has an empty field " + std::to_string(field_count) +
                           "; declare 0 for a modeled zero cost";
                }
                std::stringstream value_stream(field);
                double value = 0.0;
                if(!(value_stream >> value)) {
                    return "[" + section.name + "] " + it->first + " = '" + it->second +
                           "' is not a number";
                }
                value_stream >> std::ws;
                if(!value_stream.eof()) {
                    return "[" + section.name + "] " + it->first + " = '" + it->second +
                           "' has trailing characters";
                }
                if(value < 0.0) {
                    return "[" + section.name + "] " + it->first + " = '" + it->second +
                           "' is negative; an energy unit cost cannot be negative";
                }
                if(value != value || value > 1e300) {
                    return "[" + section.name + "] " + it->first + " = '" + it->second +
                           "' is not finite";
                }
                parsed_any = true;
            }
            if(!parsed_any) {
                return "[" + section.name + "] " + it->first + " is declared with no value; use "
                       "0 for a modeled zero cost or omit the key";
            }
            // RE5: arity. A short per-datatype vector left the trailing datatype(s) at zero, so a
            // two-field `read_energy` priced input and weight and silently made output free.
            const unsigned expected = (schema->kind == ENERGY_KEY_DATATYPE_VECTOR)
                                    ? static_cast<unsigned>(data_type_t::NUM_DATA_TYPES) : 1u;
            if(field_count != expected) {
                return "[" + section.name + "] " + it->first + " = '" + it->second + "' has " +
                       std::to_string(field_count) + " value(s) but the schema requires " +
                       std::to_string(expected) +
                       (expected == 1u ? " (a scalar cost)"
                                       : " (one cost per datatype: input:weight:output)");
            }
        }
        // RE3: an energy needs an event source. `layer_setup_energy > 0` with no setup cycle
        // charges energy for a setup that never executes, which decouples the latency and energy
        // event sources -- the report would show a setup cost for a schedule with no setup in it.
        double setup_energy = 0.0, setup_cycle = 0.0;
        const bool has_setup_energy = section.get_setting("layer_setup_energy", &setup_energy);
        section.get_setting("layer_setup_cycle", &setup_cycle);
        if(has_setup_energy && setup_energy > 0.0 && setup_cycle <= 0.0) {
            return "[" + section.name + "] layer_setup_energy = " + std::to_string(setup_energy) +
                   " but layer_setup_cycle is 0; a setup energy with no setup execution has no "
                   "event source (declare layer_setup_cycle, or drop the energy)";
        }
        double stripe_energy = 0.0, stripe_cycle = 0.0;
        const bool has_stripe_energy =
            section.get_setting("stripe_transition_energy", &stripe_energy);
        section.get_setting("stripe_transition_cycle", &stripe_cycle);
        if(has_stripe_energy && stripe_energy > 0.0 && stripe_cycle == 0.0) {
            return "[" + section.name + "] stripe_transition_energy = " +
                   std::to_string(stripe_energy) +
                   " but stripe_transition_cycle is 0; a transition energy with no modeled "
                   "transition has no event source (declare stripe_transition_cycle, or drop it)";
        }
    }
    return std::string();
}

// RE5: the energy keys each report row needs, and where they live. The PE-array config section
// (`[systolic_array]` / `[spatial_arch]` / `[adder_tree]`) carries three rows' costs, so the
// required-key sets are per ROW rather than per section.
namespace {

struct component_requirement_t {
    energy_component_t component;
    const char *section;       // "pe_array", "global_buffer", "multi_chip" or "dram"
    const char *keys[5];       // required; absence means the component is only PARTIALLY costed
    // Optional-feature costs. These price a feature the architecture may not model at all, so
    // their absence is not a gap -- but a declared non-zero one DOES make the component costed.
    // Without this distinction a config that prices only its layer setup showed a non-zero
    // subtotal annotated "modeled zero".
    const char *optional[5];
};

const component_requirement_t g_requirements[] = {
    { ENERGY_COMPONENT_MAC,           "pe_array",
      { "computation_energy", "mac_read_energy", "mac_write_energy", 0, 0 },
      // accumulator_spill_energy: without edge_accumulation the accumulator sits in the MAC
      // and its (declared, possibly calibrated) energy lands in this row -- so the row must
      // read PARTIAL, not "NOT MODELED", when that is its only declared cost. The key stays
      // in the PE roster too for the edge_accumulation ownership.
      { "mac_energy_", "mac_reduction_energy", "accumulator_spill_energy", 0, 0 } },
    { ENERGY_COMPONENT_PE,            "pe_array",
      { "lb_read_energy", "lb_write_energy", "static_energy", "transfer_energy_pe", 0 },
      { "accumulator_spill_energy", 0, 0, 0, 0 } },
    { ENERGY_COMPONENT_PE_ARRAY,      "pe_array",
      { "pe_array_read_energy", "pe_array_write_energy", "noc_energy", 0, 0 },
      { "pe_array_static_energy", "layer_setup_energy", "weight_fold_fill_energy",
        "stripe_transition_energy", 0 } },
    { ENERGY_COMPONENT_GLOBAL_BUFFER, "global_buffer",
      { "read_energy", "write_energy", "static_energy", "transfer_energy", 0 },
      { 0, 0, 0, 0, 0 } },
    { ENERGY_COMPONENT_MULTI_CHIP,    "multi_chip",
      { "read_energy", "write_energy", "static_energy", "noc_energy", 0 },
      { "output_cast_energy", 0, 0, 0, 0 } },
    { ENERGY_COMPONENT_DRAM,          "dram",
      { "read_energy", "write_energy", "transfer_energy", 0, 0 },
      { "row_miss_energy", 0, 0, 0, 0 } },
    // SFU: every invocation moves elements through its input/output ports, so those two
    // costs are required; the per-operation family, setup and leakage price optional
    // features (their absence on an ACTIVE event still surfaces as UNPRICED_ACTIVE).
    { ENERGY_COMPONENT_SFU,           "sfu",
      { "sfu_read_energy", "sfu_write_energy", 0, 0, 0 },
      { "sfu_op_energy_", "sfu_setup_energy", "sfu_static_energy", 0 } },
};

// Which report-row section a config section belongs to; 0 for sections that carry no energy.
const char *section_role(const std::string &m_name) {
    if(m_name == "systolic_array" || m_name == "spatial_arch" || m_name == "adder_tree") {
        return "pe_array";
    }
    if(m_name == "separate" || m_name == "shared") return "global_buffer";
    if(m_name == "multi_chip") return "multi_chip";
    if(m_name == "dram") return "dram";
    if(m_name == "sfu") return "sfu";
    return 0;
}

// True when any field of a declared cost is non-zero.
bool any_nonzero(const std::string &m_value) {
    std::stringstream fields(m_value);
    std::string field;
    while(std::getline(fields, field, ':')) {
        std::stringstream value_stream(field);
        double value = 0.0;
        if((value_stream >> value) && value != 0.0) return true;
    }
    return false;
}

// Looks a key up in every section playing the given role. Returns whether it was declared and
// increments the non-zero tally when it carries a value. `mac_energy_` matches by prefix, since it
// names a family of per-precision costs rather than one key.
bool declared_nonzero(config_t &m_config, const std::string &m_role, const std::string &m_key,
                      unsigned &m_nonzero) {
    const bool prefix = !m_key.empty() && m_key[m_key.size()-1] == '_';
    bool found = false;
    for(unsigned i = 0; i < m_config.sections.size(); ++i) {
        const char *role = section_role(m_config.sections[i].name);
        if(role == 0 || m_role != role) continue;
        const std::map<std::string, std::string> &settings = m_config.sections[i].all_settings();
        for(std::map<std::string, std::string>::const_iterator it = settings.begin();
            it != settings.end(); ++it) {
            const bool match = prefix ? (it->first.size() > m_key.size() &&
                                         it->first.compare(0, m_key.size(), m_key) == 0)
                                      : (it->first == m_key);
            if(!match) continue;
            found = true;
            if(any_nonzero(it->second)) ++m_nonzero;
        }
    }
    return found;
}

}  // namespace

energy_cost_schema_t::energy_cost_schema_t() : sfu_declared(false) {
    for(unsigned i = 0; i < NUM_ENERGY_COMPONENTS; ++i) {
        state[i] = ENERGY_COST_NOT_MODELED;
    }
}

bool energy_cost_schema_t::is_declared(const char *m_section, const char *m_key) const {
    return declared.count(std::string(m_section) + "/" + m_key) != 0;
}

void energy_cost_schema_t::configure(config_t &m_config) {
    // Record every energy key the config declares, keyed by the role of the section it is in.
    declared.clear();
    sfu_declared = false;
    for(unsigned i = 0; i < m_config.sections.size(); ++i) {
        if(m_config.sections[i].name == "sfu") sfu_declared = true;
    }
    for(unsigned i = 0; i < m_config.sections.size(); ++i) {
        const char *role = section_role(m_config.sections[i].name);
        if(role == 0) continue;
        const std::map<std::string, std::string> &settings = m_config.sections[i].all_settings();
        for(std::map<std::string, std::string>::const_iterator it = settings.begin();
            it != settings.end(); ++it) {
            if(it->first.find("energy") != std::string::npos) {
                declared.insert(std::string(role) + "/" + it->first);
            }
        }
    }

    const unsigned requirements = sizeof(g_requirements)/sizeof(g_requirements[0]);
    for(unsigned r = 0; r < requirements; ++r) {
        const component_requirement_t &need = g_requirements[r];
        unsigned declared = 0, required = 0, nonzero = 0;
        std::string missing_keys;
        for(unsigned k = 0; k < 5 && need.keys[k] != 0; ++k) {
            ++required;
            if(declared_nonzero(m_config, need.section, need.keys[k], nonzero)) { ++declared; }
            else {
                if(!missing_keys.empty()) missing_keys += ", ";
                missing_keys += need.keys[k];
            }
        }
        // 2026-08-26: count optional declarations too (and scan the widened 5-slot list --
        // the bound had stayed at the old array size). The struct's own contract says a
        // declared optional cost MAKES the component costed; the code only honored that when
        // a required key was also declared, so a component whose ONLY declared cost was an
        // optional one (nvdla_small_cacti22: the CACTI-calibrated accumulator on a MAC with
        // no compute cost) printed "NOT MODELED" while carrying calibrated energy. With an
        // optional declared and required keys missing it is a PARTIAL undercount, which is
        // the truthful state.
        unsigned optional_declared = 0;
        for(unsigned k = 0; k < 5 && need.optional[k] != 0; ++k) {
            if(declared_nonzero(m_config, need.section, need.optional[k], nonzero)) {
                ++optional_declared;
            }
        }
        missing[need.component] = missing_keys;
        if(declared == 0 && optional_declared == 0)
                                      state[need.component] = ENERGY_COST_NOT_MODELED;
        else if(declared < required)  state[need.component] = ENERGY_COST_PARTIAL;
        else if(nonzero == 0)         state[need.component] = ENERGY_COST_MODELED_ZERO;
        else                          state[need.component] = ENERGY_COST_CALIBRATED;
    }
}

std::string energy_cost_schema_t::annotate(energy_component_t m_component,
                                           bool m_has_energy) const {
    switch(state[m_component]) {
        case ENERGY_COST_NOT_MODELED:
            return "   <- NOT MODELED: no energy cost declared for this component";
        case ENERGY_COST_PARTIAL:
            return "   <- PARTIAL cost (UNDERCOUNT): missing " + missing[m_component];
        case ENERGY_COST_MODELED_ZERO:
            return "   <- modeled zero: every cost is declared and declared 0";
        default:
            // RE5: a priced component at zero means the layer exercised it zero times. That is a
            // property of the workload, and saying so keeps it from reading as a missing cost.
            return m_has_energy ? std::string()
                                : std::string("   <- costed, but NO ACTIVITY in this layer");
    }
}

unsigned energy_cost_schema_t::undercounted() const {
    unsigned count = 0;
    for(unsigned i = 0; i < NUM_ENERGY_COMPONENTS; ++i) {
        // The SFU row exists only for configs that model an SFU at all: a config without an
        // [sfu] section has activation OUT OF SCOPE (stated in the report), not undercounted.
        if(i == ENERGY_COMPONENT_SFU && !sfu_declared) continue;
        if(state[i] == ENERGY_COST_NOT_MODELED || state[i] == ENERGY_COST_PARTIAL) ++count;
    }
    return count;
}

unsigned energy_cost_schema_t::components_in_scope() const {
    return sfu_declared ? NUM_ENERGY_COMPONENTS : NUM_ENERGY_COMPONENTS - 1;
}

energy_cost_schema_t &energy_cost_schema() {
    static energy_cost_schema_t instance;
    return instance;
}
