#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <iostream>
#include <algorithm>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <type_traits>
#include "utils.h"

class section_config_t {
public:
    section_config_t(std::string m_name);
    ~section_config_t();

    // Add (key, value) pair to the latest section settings
    bool add_setting(std::string m_key, std::string m_value);
    // Check if a setting exists.
    bool exists(std::string m_key);
    // Get the setting value, Return true if found.
    // A key that IS present but whose value cannot be parsed into the requested type is a CONFIG
    // ERROR, not an absent setting. Returning false there made the declared value vanish and left
    // the caller's default in place, with nothing printed. Two instances of that had shipped:
    //   * `[multi_chip] input_size = 4.5` -- the field was `unsigned`, so the fractional value
    //     failed the eof() check and the buffer stayed at 0. eyerissv2.cfg could not run at all for
    //     as long as it shipped, and nothing noticed because no gate ran that config.
    //   * `exist_temporal_buffer = 2` -- a bool extraction of "2" fails, so the flag silently kept
    //     its default. This is also what made a x2 perturbation of a boolean knob look like a dead
    //     knob in validation/knobs.
    // Both are now hard errors naming the section, the key and the value.
    template <typename T>
    bool get_setting(std::string m_key, T *m_var) {
        std::map<std::string, std::string>::iterator it = settings.find(lowercase(m_key));
        if(it == settings.end() || it->second.empty()) return false;
        if(std::is_unsigned<T>::value && it->second[0] == '-') {
            reject_setting(m_key, it->second, "must not be negative");
        }
        T parsed;
        std::stringstream ss(it->second);
        if(!(ss >> parsed)) {
            reject_setting(m_key, it->second, "is not a valid value for this setting's type"
                                              " (a boolean flag accepts only 0 or 1)");
        }
        ss >> std::ws;
        if(!ss.eof()) {
            reject_setting(m_key, it->second, "has trailing characters -- an integer setting given"
                                              " a fractional value lands here");
        }
        *m_var = parsed;
        return true;
    }

    // A string setting takes the WHOLE value, spaces included: a provenance string is not a token.
    // Stream extraction stopped at the first space and then failed the eof() check, so a
    // multi-word `energy_reference` would have been dropped in silence.
    bool get_setting(std::string m_key, std::string *m_var);

    template <typename T>
    bool get_vector_setting(std::string m_key, std::vector<T> *m_vector) {
        std::map<std::string, std::string>::iterator it = settings.find(lowercase(m_key));
        if(it == settings.end() || m_vector->empty()) return false;
        std::stringstream ss(it->second);
        T first_value;
        // Same contract as the scalar form: a declared vector that does not parse is an error, not
        // an absent setting. The old code could also leave the vector HALF updated before giving
        // up, so a malformed value produced a mix of declared and default entries.
        if(!(ss >> first_value)) {
            reject_setting(m_key, it->second, "is not a valid per-datatype value list");
        }
        std::vector<T> parsed(m_vector->size(), first_value);
        ss >> std::ws;
        if(ss.eof()) {
            *m_vector = parsed;
            return true;
        }
        for(unsigned i = 1; i < m_vector->size(); i++) {
            char delimiter;
            if(!(ss >> delimiter) || delimiter != ':') {
                reject_setting(m_key, it->second, "expects values separated by ':'");
            }
            T temp_value;
            if(!(ss >> temp_value)) {
                reject_setting(m_key, it->second, "has a field that is not a valid value");
            }
            parsed.at(i) = temp_value;
        }
        ss >> std::ws;
        if(!ss.eof()) {
            reject_setting(m_key, it->second, "has more fields than this setting has datatypes");
        }
        *m_vector = parsed;
        return true;
    }

    std::string name;

    // Abort with a message naming this section, the key and the value. Never returns.
    void reject_setting(const std::string &m_key, const std::string &m_value,
                        const char *m_why) const;

    // E8: read-only view of this section's raw (key, value) pairs, so a whole-config validator
    // can check every setting of a KIND -- e.g. every energy unit cost -- without each
    // component having to remember to guard its own keys. Only 3 of the 22 energy keys had a
    // non-negative guard before this.
    const std::map<std::string, std::string> &all_settings() const { return settings; }

private:
    // Section name
    std::map<std::string, std::string> settings;

};

// Configuration.
class config_t {
public:
    config_t();
    ~config_t();

    // Parse configuration file
    void parse(std::string m_config_name);
    std::vector<section_config_t> sections;
};

#endif
