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
    const std::vector<std::string> mapping_ids = graph.mapping_operation_ids();
    const std::string legacy = graph.legacy_network_config();

    // 2026-08-28: the test dispatches on the loaded model so BOTH checked-in fixtures --
    // the Linear one and the Conv2d one -- pin their loader/adapter contract here.
    if(graph.model_name == "linear-relu-fixture") {
        assert(graph.operations.size() == 1);
        assert(mapping_ids.size() == 1);
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
        assert(graph.operations.size() == 1);
        assert(mapping_ids.size() == 1);
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
    } else if(graph.model_name == "residual-dag-fixture") {
        assert(graph.operations.size() == 6);
        assert(mapping_ids.size() == 1 && mapping_ids[0] == "conv");
        assert(graph.operations[0].kind == WORKLOAD_CONV2D);
        assert(graph.operations[1].kind == WORKLOAD_BATCH_NORM);
        assert(graph.operations[2].kind == WORKLOAD_ELEMENTWISE);
        assert(graph.operations[2].inputs[1] == "input");
        assert(graph.operations[3].kind == WORKLOAD_POOL2D);
        assert(graph.operations[3].geometry.mode == "max");
        assert(graph.operations[4].kind == WORKLOAD_POOL2D);
        assert(graph.operations[4].geometry.mode == "average");
        assert(graph.operations[5].kind == WORKLOAD_CONCAT);
        assert(graph.operations[5].geometry.axis == 1);
        assert(legacy.find("[convolutional]") != std::string::npos);

        std::map<std::string, size_t> bytes;
        for(const workload_tensor_t &tensor : graph.tensors) bytes[tensor.id] = tensor.elements();
        workload_lifetime_t lifetime(graph, 1024, bytes);
        workload_residency_plan_t plan = lifetime.plan(0);
        assert(plan.inputs[0] == WORKLOAD_RESIDENCY_DRAM && plan.retain_output);
        assert(plan.retain_inputs[0]);
        lifetime.commit(0, &plan);
        assert(lifetime.resident("input") && lifetime.resident("conv_output"));
        plan = lifetime.plan(1);
        assert(plan.inputs[0] == WORKLOAD_RESIDENCY_GLB && plan.retain_output);
        lifetime.commit(1, &plan);
        assert(!lifetime.resident("conv_output") && lifetime.resident("bn_output"));
        plan = lifetime.plan(2);
        assert(plan.inputs[0] == WORKLOAD_RESIDENCY_GLB);
        assert(plan.inputs[1] == WORKLOAD_RESIDENCY_GLB);
        lifetime.commit(2, &plan);
        assert(!lifetime.resident("input") && !lifetime.resident("bn_output"));
        plan = lifetime.plan(3);
        assert(plan.inputs[0] == WORKLOAD_RESIDENCY_GLB);
        lifetime.commit(3, &plan);
        assert(lifetime.resident("residual") && lifetime.resident("max_pool"));
        plan = lifetime.plan(4);
        lifetime.commit(4, &plan);
        assert(!lifetime.resident("residual"));
        plan = lifetime.plan(5);
        assert(plan.inputs[0] == WORKLOAD_RESIDENCY_GLB);
        assert(plan.inputs[1] == WORKLOAD_RESIDENCY_GLB);
        assert(!plan.retain_output);
        lifetime.commit(5, &plan);
        assert(lifetime.occupied_bytes() == 0);

        // With space for only one tensor, the current output has priority over
        // pinning the skip input; the latter remains explicitly DRAM-backed.
        workload_lifetime_t constrained(graph, 256, bytes);
        plan = constrained.plan(0);
        assert(plan.retain_output);
        assert(!plan.retain_inputs[0]);
        constrained.commit(0, &plan);
        assert(constrained.resident("conv_output"));
        assert(!constrained.resident("input"));
        assert(constrained.occupied_bytes() == 256);
    } else if(graph.model_name == "lenet-fixture") {
        // Classifier chain across an elided flatten: the linear reads pool storage
        // through an alias view, with no operation and no separate allocation.
        assert(graph.operations.size() == 7);
        assert(graph.operations[0].kind == WORKLOAD_CONV2D);
        assert(graph.operations[1].kind == WORKLOAD_POOL2D);
        assert(graph.operations[2].kind == WORKLOAD_CONV2D);
        assert(graph.operations[3].kind == WORKLOAD_POOL2D);
        assert(graph.operations[4].kind == WORKLOAD_LINEAR);
        assert(graph.operations[5].kind == WORKLOAD_LINEAR);
        assert(graph.operations[6].kind == WORKLOAD_SOFTMAX);
        assert(mapping_ids.size() == 4);
        assert(mapping_ids[0] == "conv2d" && mapping_ids[3] == "linear_1");
        assert(graph.tensor("flatten").alias_of == "max_pool2d_1");
        assert(graph.storage_tensor("flatten").id == "max_pool2d_1");
        assert(graph.operations[4].inputs[0] == "flatten");
        assert(graph.operations[4].geometry.input_features == 256);

        std::map<std::string, size_t> bytes;
        for(const workload_tensor_t &tensor_value : graph.tensors) bytes[tensor_value.id] = tensor_value.elements();
        workload_lifetime_t lifetime(graph, 1 << 20, bytes);
        for(size_t index = 0; index < 4; ++index) {
            workload_residency_plan_t plan = lifetime.plan(index);
            lifetime.commit(index, &plan);
        }
        // After the second pool, its storage must stay resident: the alias view
        // "flatten" still has the linear as a future consumer.
        assert(lifetime.resident("max_pool2d_1"));
        assert(lifetime.resident("flatten"));
        workload_residency_plan_t plan = lifetime.plan(4);
        assert(plan.inputs[0] == WORKLOAD_RESIDENCY_GLB);
        lifetime.commit(4, &plan);
        // The linear was the storage's last reader; both names now agree it is gone.
        assert(!lifetime.resident("max_pool2d_1"));
        assert(!lifetime.resident("flatten"));
        plan = lifetime.plan(5);
        lifetime.commit(5, &plan);
        plan = lifetime.plan(6);
        assert(!plan.retain_output);
        lifetime.commit(6, &plan);
        assert(lifetime.occupied_bytes() == 0);
    } else if(graph.model_name == "pool-chain-fixture") {
        assert(graph.operations.size() == 3);
        assert(graph.operations[0].kind == WORKLOAD_CONV2D);
        assert(graph.operations[1].kind == WORKLOAD_POOL2D);
        assert(graph.operations[1].geometry.mode == "max");
        assert(graph.operations[2].kind == WORKLOAD_POOL2D);
        assert(graph.operations[2].geometry.mode == "average");
        // pools take no MAC mapping: only the conv appears in the mapping id list.
        assert(mapping_ids.size() == 1);
        assert(mapping_ids[0] == "conv2d");
    } else {
        std::cerr << "workload_graph_test: unknown fixture model " << graph.model_name
                  << std::endl;
        return 2;
    }
    std::cout << "workload_graph_test: PASS (" << graph.model_name << ")" << std::endl;
    return 0;
}
