#include <cassert>
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
	std::string mapping= argv[4];

    std::string accelerator_config = "accelerators/" + accelerator + ".cfg";
    std::string network_config = "networks/" + network + ".cfg";
    std::string mapping_config = "mappings/" + accelerator + "/" + network +"/" + mapping + ".map";

    if(run_type != "run") {
        std::cerr << "Unknown run type " << run_type << std::endl;
        exit(1);
    }

    npu_t *npu = new npu_t();
    // Initialize the accelerator.
    npu->init(accelerator_config, network_config, mapping_config);
    // Run the accelerator.
    npu->run(accelerator, network);

    std::string command = "./move.sh " + accelerator + " " + network + " " + mapping;
    std::system(command.c_str());
    
	delete npu;

    return 0;
}
