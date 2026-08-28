#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>

#include "workload_graph.h"

int main(int argc, char **argv) {
    if(argc != 2) {
        std::cerr << "usage: workload_graph_test EXECUTABLE_IR" << std::endl;
        return 2;
    }

    workload_graph_t graph;
    graph.load(argv[1]);
    assert(graph.schema_version == "npusim.exec.v1");
    assert(graph.operations.size() == 1);
    const std::vector<std::string> mapping_ids = graph.mapping_operation_ids();
    assert(mapping_ids.size() == 1);
    const std::string legacy = graph.legacy_network_config();

    // 2026-08-28: the test dispatches on the loaded model so BOTH checked-in fixtures --
    // the Linear one and the Conv2d one -- pin their loader/adapter contract here.
    if(graph.model_name == "linear-relu-fixture") {
        assert(graph.operations[0].id == "linear");
        assert(graph.operations[0].kind == WORKLOAD_LINEAR);
        assert(graph.operations[0].activation == "relu");
        assert(graph.operations[0].geometry.batch == 64);
        assert(graph.operations[0].geometry.input_features == 64);
        assert(graph.operations[0].geometry.output_features == 64);
        assert(mapping_ids[0] == "linear");
        assert(legacy.find("[connected]") != std::string::npos);
        assert(legacy.find("output=64") != std::string::npos);
        assert(legacy.find("activation=relu") != std::string::npos);
    } else if(graph.model_name == "conv-relu-fixture") {
        assert(graph.operations[0].id == "conv");
        assert(graph.operations[0].kind == WORKLOAD_CONV2D);
        assert(graph.operations[0].activation == "relu");
        assert(graph.operations[0].geometry.input_channels == 3);
        assert(graph.operations[0].geometry.output_channels == 8);
        assert(graph.operations[0].geometry.filter_height == 3);
        assert(graph.operations[0].geometry.output_height == 6);
        assert(mapping_ids[0] == "conv");
        assert(legacy.find("[convolutional]") != std::string::npos);
        assert(legacy.find("filters=8") != std::string::npos);
        assert(legacy.find("channels=3") != std::string::npos);
        assert(legacy.find("activation=relu") != std::string::npos);
    } else {
        std::cerr << "workload_graph_test: unknown fixture model " << graph.model_name
                  << std::endl;
        return 2;
    }
    std::cout << "workload_graph_test: PASS (" << graph.model_name << ")" << std::endl;
    return 0;
}
