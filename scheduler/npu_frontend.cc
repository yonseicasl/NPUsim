#include "npu.h"

#include <cstdlib>
#include <iostream>
#include <map>
#include <set>

void npu_t::bind_executable_mappings() {
    if(workload == NULL) return;

    if(sfus.empty()) {
        for(const workload_operation_t &operation : workload->operations) {
            if(operation.kind == WORKLOAD_SOFTMAX || operation.activation != "linear") {
                std::cerr << "Error: executable operation " << operation.id
                          << " requires an SFU, but the accelerator has no [sfu] section"
                          << std::endl;
                exit(1);
            }
        }
    }

    const std::vector<std::string> expected_ids = workload->mapping_operation_ids();
    if(mapping_tables.size() != expected_ids.size()) {
        std::cerr << "Error: executable IR requires " << expected_ids.size()
                  << " MAC mappings, but the mapping file contains "
                  << mapping_tables.size() << std::endl;
        exit(1);
    }

    size_t named_mappings = 0;
    for(mapping_table_t *mapping : mapping_tables) {
        if(!mapping->get_operation_id().empty()) ++named_mappings;
    }
    if(named_mappings == 0) {
        std::cerr << "Warning: executable IR mapping file has no op_id fields; "
                  << "using legacy ordinal binding" << std::endl;
        return;
    }
    if(named_mappings != mapping_tables.size()) {
        std::cerr << "Error: executable IR mapping file mixes named and ordinal sections; "
                  << "declare op_id on every section or none" << std::endl;
        exit(1);
    }

    std::map<std::string, mapping_table_t*> by_id;
    for(mapping_table_t *mapping : mapping_tables) {
        if(!by_id.insert(std::make_pair(mapping->get_operation_id(), mapping)).second) {
            std::cerr << "Error: duplicate mapping op_id " << mapping->get_operation_id()
                      << std::endl;
            exit(1);
        }
    }

    std::vector<mapping_table_t*> ordered;
    ordered.reserve(expected_ids.size());
    std::set<std::string> consumed;
    for(const std::string &operation_id : expected_ids) {
        const std::map<std::string, mapping_table_t*>::const_iterator found = by_id.find(operation_id);
        if(found == by_id.end()) {
            std::cerr << "Error: no mapping section has op_id " << operation_id << std::endl;
            exit(1);
        }

        const workload_operation_t *operation = NULL;
        for(const workload_operation_t &candidate : workload->operations) {
            if(candidate.id == operation_id) {
                operation = &candidate;
                break;
            }
        }
        if(operation == NULL) {
            std::cerr << "Error: internal executable mapping lookup failed for "
                      << operation_id << std::endl;
            exit(1);
        }
        const std::string expected_kind = operation->kind == WORKLOAD_LINEAR
            ? "connected" : "convolutional";
        if(found->second->get_mapping_kind() != expected_kind) {
            std::cerr << "Error: mapping op_id " << operation_id << " has section kind "
                      << found->second->get_mapping_kind() << ", expected "
                      << expected_kind << std::endl;
            exit(1);
        }
        ordered.push_back(found->second);
        consumed.insert(operation_id);
    }
    if(consumed.size() != by_id.size()) {
        std::cerr << "Error: mapping file contains an op_id not used by the executable IR"
                  << std::endl;
        exit(1);
    }
    mapping_tables.swap(ordered);
}

void npu_t::print_workload_provenance(std::ofstream &m_output) const {
    if(workload == NULL) return;
    m_output << "Frontend schema: " << workload->schema_version << std::endl;
    m_output << "Frontend model: " << workload->model_name << std::endl;
    m_output << "Frontend executable SHA256 (declared): "
             << workload->executable_sha256 << std::endl;
    m_output << "Frontend source graph: " << workload->source_schema_version
             << " sha256=" << workload->source_sha256 << std::endl;
    m_output << "Frontend mode: timing-only; accelerator config owns runtime datatypes"
             << std::endl << std::endl;
}
