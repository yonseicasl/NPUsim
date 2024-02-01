#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "utils.h"

class section_config_t {
public:
    section_config_t(std::string m_name);
    ~section_config_t();

    // Add (key, value) pair to the latest section settings
    void add_setting(std::string m_key, std::string m_value);
    // Check if a setting exists.
    bool exists(std::string m_key);
    // Get the setting value, Return true if found.
    template <typename T>
    bool get_setting(std::string m_key, T *m_var) {
        std::map<std::string, std::string>::iterator it = settings.find(lowercase(m_key));
        if(it == settings.end()) return false;
        std::stringstream ss; ss.str(it->second);
        ss >> *m_var; return true;
    }

    template <typename T>
        bool get_vector_setting(std::string m_key, std::vector<T> *m_vector) {
            std::map<std::string, std::string>::iterator it = settings.find(lowercase(m_key));
            if(it == settings.end()) return false;
            std::stringstream ss(it->second);
            T temp_value;
            for(unsigned i = 0; i < m_vector->size(); i++) {
                ss >> temp_value;
                m_vector->at(i) = temp_value;
                if(ss.peek() == ':') {
                    ss.ignore();
                }
            }
            return true;
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
