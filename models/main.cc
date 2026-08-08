#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "convolutional.h"
#include "npu.h"

int main(int argc, char **argv) {
    if(argc < 5) {
        std::cerr << "Usage : " << argv[0] 
                  << " run "                    
                  << " [accelerator config]"        // Type Accelerator specification.
                  << " [network config]"            // Type Network configuration.
				  << " [mapping table config]"      // Type accelerator cost.
                  << std::endl;
        exit(1);
    }

    std::string run_type = argv[1];                 
    std::string accelerator = argv[2];
    std::string network= argv[3];
	std::string mapping = argv[4];

    const char *config_root_env = std::getenv("NPUSIM_CONFIG_ROOT");
    const std::string config_root = config_root_env ? config_root_env : "../configs";
    const std::string accelerator_config = config_root + "/accelerators/" + accelerator + ".cfg";
    const std::string network_config = config_root + "/networks/" + network + ".cfg";
    const std::string mapping_config = config_root + "/mappings/" + accelerator + "/" + network + "/" + mapping + ".map";

    if(run_type != "run") {
        std::cerr << "Unknown run type " << run_type << std::endl;
        exit(1);
    }

    const std::string configs[] = {accelerator_config, network_config, mapping_config};
    for(const std::string &config : configs) {
        std::ifstream input(config.c_str());
        if(!input.good()) {
            std::cerr << "Configuration file not found: " << config << std::endl;
            return 1;
        }
    }

    npu_t *npu = new npu_t();
    // Initialize the accelerator.
    npu->init(accelerator_config, network_config, mapping_config);
    // Run the accelerator.
    npu->run(accelerator, network);

	delete npu;

    return 0;
}
