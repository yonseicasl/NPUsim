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
    template <typename T>
    bool get_setting(std::string m_key, T *m_var) {
        std::map<std::string, std::string>::iterator it = settings.find(lowercase(m_key));
        if(it == settings.end() || it->second.empty()) return false;
        if(std::is_unsigned<T>::value && it->second[0] == '-') return false;
        T parsed;
        std::stringstream ss(it->second);
        if(!(ss >> parsed)) return false;
        ss >> std::ws;
        if(!ss.eof()) return false;
        *m_var = parsed;
        return true;
    }

    template <typename T>
    bool get_vector_setting(std::string m_key, std::vector<T> *m_vector) {
        std::map<std::string, std::string>::iterator it = settings.find(lowercase(m_key));
        if(it == settings.end() || m_vector->empty()) return false;
        std::stringstream ss(it->second);
        T first_value;
        if(!(ss >> first_value)) return false;
        m_vector->at(0) = first_value;
        ss >> std::ws;
        if(ss.eof()) {
            std::fill(m_vector->begin(), m_vector->end(), first_value);
            return true;
        }
        for(unsigned i = 1; i < m_vector->size(); i++) {
            char delimiter;
            if(!(ss >> delimiter) || delimiter != ':') return false;
            T temp_value;
            if(!(ss >> temp_value)) return false;
            m_vector->at(i) = temp_value;
        }
        ss >> std::ws;
        return ss.eof();
    }

    std::string name;

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
