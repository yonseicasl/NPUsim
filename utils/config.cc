#include <algorithm>
#include <cctype>
#include <iostream>
#include <locale>
#include <fstream>
#include <cstdlib>
#include "config.h"

section_config_t::section_config_t(std::string m_name) :
    name(m_name) {

}

section_config_t::~section_config_t() {
    settings.clear();
}

bool section_config_t::add_setting(std::string m_key, std::string m_value) {
    lowercase(m_key);
    return settings.insert(std::pair<std::string, std::string>(m_key, m_value)).second;
}

bool section_config_t::exists(std::string m_key) {
    return settings.find(lowercase(m_key)) != settings.end();
}

config_t::config_t() {

}

config_t::~config_t() {
    sections.clear();
}

void config_t::parse(std::string m_config_name) {
    sections.clear();
    std::fstream file_stream;
    file_stream.open(m_config_name.c_str(), std::fstream::in);
    if(!file_stream.is_open()) {
        std::cerr << "Error: failed to open " << m_config_name << std::endl;
        exit(1);
    }

    std::string line;
    unsigned line_number = 0;

    while(getline(file_stream, line)) {
        // Erase all spaces
        line_number++;
        line.erase(remove_if(line.begin(), line.end(),
                             [](unsigned char c) { return std::isspace(c); }), line.end());
        // Skip blank lines or comments
        if(!line.size() || (line[0] == '#')) continue;
        // Beginning of [section]
        if(line[0] == '[') {
            if(line.size() < 3 || line.back() != ']') {
                std::cerr << "Error: malformed section at " << m_config_name << ":" << line_number << std::endl;
                exit(1);
            }
            std::string section_name = line.substr(1, line.size()-2);
            lowercase(section_name);
            if(section_name.empty()) {
                std::cerr << "Error: empty section at " << m_config_name << ":" << line_number << std::endl;
                exit(1);
            }
            sections.push_back(section_config_t(section_name));
        }
        else {
            if(sections.empty()) {
                std::cerr << "Error: setting before section at " << m_config_name << ":" << line_number << std::endl;
                exit(1);
            }
            size_t eq = line.find('=');
            if(eq == std::string::npos || eq == 0 || eq + 1 == line.size()) {
                std::cerr << "Error: invalid setting at " << m_config_name << ":" << line_number
                          << ": " << line << std::endl;
                exit(1);
            }
            // Save (key, value) pair in the latest section setting.
            std::string key   = line.substr(0, eq);
            std::string value = line.substr(eq+1, line.size()-1);
            if(!sections.back().add_setting(key, value)) {
                std::cerr << "Error: duplicate key '" << key << "' at "
                          << m_config_name << ":" << line_number << std::endl;
                exit(1);
            }
        }
    }
}




