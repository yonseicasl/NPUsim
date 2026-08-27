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
    assert(graph.model_name == "linear-relu-fixture");
    assert(graph.operations.size() == 1);
    assert(graph.operations[0].id == "linear");
    assert(graph.operations[0].kind == WORKLOAD_LINEAR);
    assert(graph.operations[0].activation == "relu");
    assert(graph.operations[0].geometry.batch == 64);
    assert(graph.operations[0].geometry.input_features == 64);
    assert(graph.operations[0].geometry.output_features == 64);

    const std::vector<std::string> mapping_ids = graph.mapping_operation_ids();
    assert(mapping_ids.size() == 1);
    assert(mapping_ids[0] == "linear");

    const std::string legacy = graph.legacy_network_config();
    assert(legacy.find("[connected]") != std::string::npos);
    assert(legacy.find("output=64") != std::string::npos);
    assert(legacy.find("activation=relu") != std::string::npos);
    std::cout << "workload_graph_test: PASS" << std::endl;
    return 0;
}
