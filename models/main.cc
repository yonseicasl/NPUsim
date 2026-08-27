#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "npu.h"

namespace {

bool file_exists(const std::string &path) {
    std::ifstream input(path.c_str());
    return input.good();
}

std::string file_stem(const std::string &path) {
    const std::string::size_type slash = path.find_last_of("/\\");
    std::string name = slash == std::string::npos ? path : path.substr(slash + 1);
    const std::string::size_type dot = name.find_last_of('.');
    if(dot != std::string::npos) name = name.substr(0, dot);
    if(name.empty()) name = "workload";
    for(char &value : name) {
        const unsigned char byte = static_cast<unsigned char>(value);
        if(!std::isalnum(byte) && value != '_' && value != '-') value = '_';
    }
    return name;
}

void require_file(const std::string &path) {
    if(!file_exists(path)) {
        std::cerr << "Configuration/artifact file not found: " << path << std::endl;
        exit(1);
    }
}

void print_usage(const char *program) {
    std::cerr << "Usage:\n"
              << "  " << program
              << " run [accelerator name] [network name] [mapping name]\n"
              << "  " << program
              << " run-ir [accelerator config path] [executable IR path]"
              << " [mapping path] [optional result name]" << std::endl;
}

} // namespace

int main(int argc, char **argv) {
    if(argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    const std::string run_type = argv[1];
    std::string accelerator_config;
    std::string network_config;
    std::string mapping_config;
    std::string accelerator_label;
    std::string network_label;

    if(run_type == "run") {
        if(argc != 5) {
            print_usage(argv[0]);
            return 1;
        }
        const std::string accelerator = argv[2];
        const std::string network = argv[3];
        const std::string mapping = argv[4];
        const char *config_root_env = std::getenv("NPUSIM_CONFIG_ROOT");
        const std::string config_root = config_root_env ? config_root_env : "../configs";
        accelerator_config = config_root + "/accelerators/" + accelerator + ".cfg";
        network_config = config_root + "/networks/" + network + ".cfg";
        mapping_config = config_root + "/mappings/" + accelerator + "/" +
            network + "/" + mapping + ".map";
        accelerator_label = accelerator;
        network_label = network;
    } else if(run_type == "run-ir") {
        if(argc != 5 && argc != 6) {
            print_usage(argv[0]);
            return 1;
        }
        accelerator_config = argv[2];
        network_config = argv[3];
        mapping_config = argv[4];
        accelerator_label = file_stem(accelerator_config);
        network_label = argc == 6 ? file_stem(argv[5]) : file_stem(network_config);
    } else {
        std::cerr << "Unknown run type " << run_type << std::endl;
        print_usage(argv[0]);
        return 1;
    }

    require_file(accelerator_config);
    require_file(network_config);
    require_file(mapping_config);

    npu_t npu;
    npu.init(accelerator_config, network_config, mapping_config);
    npu.run(accelerator_label, network_label);
    return 0;
}
