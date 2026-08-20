#ifndef __ENERGY_UNITS_H__
#define __ENERGY_UNITS_H__

#include <set>
#include <string>
#include "config.h"

// E7/E8: what the energy numbers in the report actually MEAN, and whether they were
// calibrated.
//
// The problem this solves: eyeriss_energy.cfg deliberately uses the published NORMALIZED
// access costs (MAC = 1, RF = 1x, GIN = 2x, GLB = 6x, DRAM = 200x per 16-bit word). That is a
// valid fixture for relative breakdown analysis -- but every energy value in the report was
// labelled "pJ", so a relative result read as an absolute one. Nothing in the config or the
// output distinguished "1 pJ, measured" from "1, normalized to a MAC".
//
// A config therefore declares its unit and, optionally, where the numbers came from:
//   energy_unit      = pJ | normalized      (default: unspecified/uncalibrated)
//   energy_reference = free text provenance (e.g. a paper table, a tool version)
// The report prints both, so a reader -- and a checker -- can tell whether an absolute total
// is meaningful. `normalized` also means absolute-power output is not meaningful, which is
// where this will matter again once power lands.
// RE2: `UNSPECIFIED` is the default, NOT pJ. Most configs declare no unit at all, and defaulting
// them to pJ made every one of them eligible for absolute totals and average power in mW -- a
// report could say "MAC energy basis: UNCALIBRATED" and "provenance: not declared" and still print
// a milliwatt figure. The arithmetic was right; the QUALIFICATION was wrong.
enum class energy_unit_t { UNSPECIFIED, PICOJOULE, NORMALIZED };

struct energy_units_t {
    energy_unit_t unit;
    std::string reference;

    energy_units_t();
    // Parse `energy_unit` / `energy_reference` from an [accelerator] section.
    void configure(section_config_t &m_section_config);
    // Short unit label for report headers.
    const char *label() const;
    // RE2: true only when the config DECLARES absolute picojoules AND says where they came from.
    // A unit alone is not calibration: `energy_unit = pJ` with no provenance is an assertion
    // nobody can check, so it does not qualify for absolute output.
    bool is_absolute() const;
    // Why absolute energy/power is or is not available, for the report.
    std::string calibration_note() const;
    // One-line provenance for the report; states explicitly when none was declared. A pJ
    // declaration is only a candidate absolute scale here: run-level active-event completeness is
    // checked later by stats_t before a total or wattage is qualified.
    std::string describe() const;
};

energy_units_t &energy_units();

// E8: reject an energy unit cost that cannot mean anything, for EVERY energy setting in the
// config rather than the handful of keys that happened to guard themselves.
//
// A negative unit cost produces negative energy, which no downstream consumer can interpret; a
// non-finite one poisons every total it enters. Both used to pass silently for most of the 22
// energy keys, and a typo in a key name still leaves that component at zero cost -- which is
// why the report now also states how many components ended up with no declared energy at all
// (see stats_t::print_energy_summary()).
//
// Returns an empty string when the config is acceptable, otherwise a message naming the
// section, the key and the offending value.
//
// RE5: value sanity is not enough -- validation is now SCHEMA-level. The old check only looked at
// the values of keys that happened to be spelled with "energy" in them, so it accepted
//   * a misspelled key (`computaton_energy`), which silently leaves that component at zero cost;
//   * a per-datatype vector of the wrong length (`read_energy = 1.0:2.0`), whose missing field
//     defaulted to zero for one datatype only;
//   * an empty middle field (`1.0::2.0`), which parsed as two values;
//   * an energy declared with no event source that can produce it (RE3: `layer_setup_energy > 0`
//     with `layer_setup_cycle = 0` charges energy for a setup that never runs).
// Every energy key is therefore declared below with its arity, and anything outside the schema is
// an error naming the nearest known key rather than a silent zero.
enum energy_key_kind_t {
    ENERGY_KEY_SCALAR,            // one value
    ENERGY_KEY_DATATYPE_VECTOR,   // exactly one value per datatype, colon separated
    ENERGY_KEY_PROVENANCE,        // free text, not a cost
    ENERGY_KEY_PREFIX_SCALAR      // a family of scalars sharing a prefix (mac_energy_<in>_<wt>)
};

struct energy_key_schema_t {
    const char *name;
    energy_key_kind_t kind;
};

// The declared energy schema, and the arity a key must have. Exposed so the unit test can assert
// that every energy key any shipped config uses is in fact declared here.
const energy_key_schema_t *energy_key_schema(unsigned &count);

std::string validate_energy_settings(config_t &m_config);

// RE5: what a zero in the energy breakdown MEANS. `Uncosted components` used to be computed from
// the result -- a subtotal of 0 was reported as "no energy cost declared" -- which conflated four
// completely different situations:
//   NOT_MODELED   the architecture declares no energy key for this component at all
//   PARTIAL       some of its costs are declared and others are missing (a typo looks like this,
//                 and so does a half-finished calibration) -- the subtotal is an UNDERCOUNT
//   MODELED_ZERO  every cost is declared and every one of them is 0 -- a deliberate modeled free
//   CALIBRATED    every cost is declared and at least one is non-zero
// A CALIBRATED component whose subtotal is 0 has simply had NO ACTIVITY in this layer, which is a
// statement about the workload, not about the config. The state is therefore derived from the
// DECLARATION and reported alongside the activity.
enum energy_cost_state_t {
    ENERGY_COST_NOT_MODELED,
    ENERGY_COST_PARTIAL,
    ENERGY_COST_MODELED_ZERO,
    ENERGY_COST_CALIBRATED
};

// The six rows of the energy breakdown.
enum energy_component_t {
    ENERGY_COMPONENT_MAC,
    ENERGY_COMPONENT_PE,
    ENERGY_COMPONENT_PE_ARRAY,
    ENERGY_COMPONENT_GLOBAL_BUFFER,
    ENERGY_COMPONENT_MULTI_CHIP,
    ENERGY_COMPONENT_DRAM,
    NUM_ENERGY_COMPONENTS
};

// E20-1/E20-2: the fifth state, and the one that matters most for an absolute claim.
//
// The four states above are properties of the CONFIG. This one is a property of the RUN: an event
// fired, and the key that prices it is absent from the config. That is not a modeled zero and not a
// half-finished component -- it is a charge the total is missing, of a size nobody declared. It
// cannot be seen in an energy breakdown, because the axis reads 0 either way; it needs the event
// COUNT next to the declaration.
//
// `gemmini_cacti22` is the case that forced this: its SRAM and DRAM costs are externally derived,
// so it satisfied every earlier absolute-energy condition and reported 418.685 mW -- while the same
// report showed a layer setup, 1,048,560 bytes of accumulator reload, 1,048,576 of spill and 4,096
// of final cast, all at 0.00 pJ because those keys were never declared.
struct unpriced_event_t {
    const char *event;      // what fired
    const char *key;        // the cost key that would price it
    const char *section;    // where that key belongs
};

struct energy_cost_schema_t {
    energy_cost_state_t state[NUM_ENERGY_COMPONENTS];
    std::string missing[NUM_ENERGY_COMPONENTS];   // keys the component needs but does not declare
    // Every energy key the config actually declares, as "role/key". Declared-as-zero is present
    // here; absent is not. The distinction is the whole point of UNPRICED_ACTIVE.
    std::set<std::string> declared;

    energy_cost_schema_t();
    // Derive every component's declaration state from the loaded config.
    void configure(config_t &m_config);
    // One-phrase annotation for the breakdown row, given whether the layer had any activity.
    std::string annotate(energy_component_t m_component, bool m_has_energy) const;
    // RE5: counted from the DECLARATION -- components whose subtotal cannot be trusted because a
    // cost is missing, not components that merely came out at zero.
    unsigned undercounted() const;
    // True when the config declares this key at all (any value, including 0).
    bool is_declared(const char *m_section, const char *m_key) const;
};

energy_cost_schema_t &energy_cost_schema();

#endif
