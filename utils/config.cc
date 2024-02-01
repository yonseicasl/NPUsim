#include <algorithm>
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

void section_config_t::add_setting(std::string m_key, std::string m_value) {
    settings.insert(std::pair<std::string, std::string>(lowercase(m_key), m_value));
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
    std::fstream file_stream;
    file_stream.open(m_config_name.c_str(), std::fstream::in);
    if(!file_stream.is_open()) {
        std::cerr << "Error: failed to open " << m_config_name << std::endl;
        exit(1);
    }

    std::string line;
    while(getline(file_stream, line)) {
        // Erase all spaces
        line.erase(remove(line.begin(),line.end(),' '),line.end());
        // Skip blank lines or comments
        if(!line.size() || (line[0] == '#')) continue;
        // Beginning of [section]
        if(line[0] == '[') {
            std::string section_name = line.substr(1, line.size()-2);
            sections.push_back(section_config_t(section_name));
        }
        else {
            size_t eq = line.find('=');
            if(eq == std::string::npos) {
                std::cerr << "Error: invalid config" << std::endl << line << std::endl;
                exit(1);
            }
            // Save (key, value) pair in the latest section setting.
            std::string key   = line.substr(0, eq);
            std::string value = line.substr(eq+1, line.size()-1);
            sections[sections.size()-1].add_setting(key, value);
        }
    }
}




