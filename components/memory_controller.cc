#include <iostream>
#include "memory_controller.h"

#define BASE_ADDR 0x80000000  //2GB
#define CAPACITY 0X200000000  //8GB

memory_controller_t::memory_controller_t(const std::string &m_config_file, const std::string &m_output_dir) :
    callback_counter(0),
    request_counter(0) {

    dramsim = dramsim3::GetMemorySystem(m_config_file, m_output_dir,
                      std::bind(&memory_controller_t::read_callback, this, std::placeholders::_1),
                      std::bind(&memory_controller_t::write_callback, this, std::placeholders::_1));

}

memory_controller_t::~memory_controller_t() {
    delete dramsim;
}

void memory_controller_t::send_request(data_t *m_data, mapping_table_t *m_mapping_table, unsigned m_offset, data_type_t m_data_type, action_type_t m_action_type) {

    std::vector<unsigned> global_buffer_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);
    std::vector<unsigned> dram_param(parameter_type_t::NUM_PARAMETER_TYPES, 1);

    global_buffer_param = m_mapping_table->calculate_parameter_size(component_type_t::GLOBAL_BUFFER);
    dram_param = m_mapping_table->calculate_parameter_size(component_type_t::DRAM);

    callback_counter = 0, request_counter = 0;
    if(m_data_type == data_type_t::INPUT) {
        for(unsigned b = 0; b < global_buffer_param[parameter_type_t::BATCH_SIZE]; b++) {
            for(unsigned c = 0; c < global_buffer_param[parameter_type_t::INPUT_CHANNEL]; c++) {
                for(unsigned h = 0; h < global_buffer_param[parameter_type_t::INPUT_HEIGHT]; h++) {
                    for(unsigned w = 0; w < global_buffer_param[parameter_type_t::INPUT_WIDTH]; w++) {
                        unsigned index = m_offset + b*dram_param[parameter_type_t::INPUT_CHANNEL]
                                         *dram_param[parameter_type_t::INPUT_HEIGHT]*dram_param[parameter_type_t::INPUT_WIDTH]
                                       + c*dram_param[parameter_type_t::INPUT_HEIGHT]*dram_param[parameter_type_t::INPUT_WIDTH]
                                       + h*dram_param[parameter_type_t::INPUT_WIDTH] + w;

                        request_t request;
                        request.address = npu_mmu::v2p((uint64_t)&m_data[index]), request.is_write = false;
                        insert_pending_queue(request);
                    }
                }
            }
        }
    }
    else if(m_data_type == data_type_t::WEIGHT) {
        for(unsigned k = 0; k < global_buffer_param[parameter_type_t::OUTPUT_CHANNEL]; k++) {
            for(unsigned c = 0; c < global_buffer_param[parameter_type_t::INPUT_CHANNEL]; c++) {
                for(unsigned r = 0; r < global_buffer_param[parameter_type_t::FILTER_HEIGHT]; r++) {
                    for(unsigned s = 0; s < global_buffer_param[parameter_type_t::FILTER_WIDTH]; s++) {
                        unsigned index = m_offset + k*dram_param[parameter_type_t::INPUT_CHANNEL]*dram_param[parameter_type_t::FILTER_HEIGHT]*dram_param[parameter_type_t::FILTER_WIDTH]
                                       + c*dram_param[parameter_type_t::FILTER_HEIGHT]*dram_param[parameter_type_t::FILTER_WIDTH]
                                       + r*dram_param[parameter_type_t::FILTER_WIDTH] + s;
                        request_t request;
                        request.address = npu_mmu::v2p((uint64_t)&m_data[index]), request.is_write = false;
                        insert_pending_queue(request);
                    }
                }
            }
        }
    }
    else if(m_data_type == data_type_t::OUTPUT) {
        for(unsigned b = 0; b < global_buffer_param[parameter_type_t::BATCH_SIZE]; b++) {
            for(unsigned k = 0; k < global_buffer_param[parameter_type_t::OUTPUT_CHANNEL]; k++) {
                for(unsigned p = 0; p < global_buffer_param[parameter_type_t::OUTPUT_HEIGHT]; p++) {
                    for(unsigned q = 0; q < global_buffer_param[parameter_type_t::OUTPUT_WIDTH]; q++) {
                        unsigned index = m_offset + b*dram_param[parameter_type_t::OUTPUT_CHANNEL]*dram_param[parameter_type_t::OUTPUT_HEIGHT]*dram_param[parameter_type_t::OUTPUT_WIDTH]
                                       + k*dram_param[parameter_type_t::OUTPUT_HEIGHT]*dram_param[parameter_type_t::OUTPUT_WIDTH]
                                       + p*dram_param[parameter_type_t::OUTPUT_WIDTH] + q;
                        
                        request_t request;
                        request.address = npu_mmu::v2p((uint64_t)&m_data[index]);
                        if(m_action_type == action_type_t::LOAD) {request.is_write = false;}
                        else if(m_action_type == action_type_t::STORE) {request.is_write = true;}

                        insert_pending_queue(request);
                    }
                }
            }
        }
    }
    while(callback_counter != request_counter || !pending_queue.empty()) {
        if(!pending_queue.empty() && will_accept_transaction(pending_queue.front())) {
            add_transaction(pending_queue.front());
        }
        tick();
    }
}

void memory_controller_t::insert_pending_queue(request_t m_request) {
    unsigned mask_bits = 0;
    unsigned bits = get_bus_bits()*get_burst_length();
    while(bits > 8) {
        bits /= 2;
        mask_bits++;
    }
    request_t request;
    request.address = (m_request.address >> mask_bits) << mask_bits, request.is_write = m_request.is_write;
    if(request.address != pending_queue.back().address) { pending_queue.push_back(request); }
}

void memory_controller_t::tick() {
    dramsim->ClockTick();
}

void memory_controller_t::read_callback(uint64_t m_address) {

    for(auto it = request_queue.cbegin(); it != request_queue.cend(); ++it) {
        if(it->address == m_address) {
            callback_counter++;
            request_queue.erase(it);
            return ;
        }
    }
}

void memory_controller_t::write_callback(uint64_t m_address) {
    for(auto it = request_queue.cbegin(); it != request_queue.cend(); ++it) {
        if(it->address == m_address) {
            callback_counter++;
            request_queue.erase(it);
            return ;
        }
    }
}

unsigned int memory_controller_t::get_bus_bits() const {
    return dramsim->GetBusBits();
}

unsigned int memory_controller_t::get_burst_length() const {
    return dramsim->GetBurstLength();
}


bool memory_controller_t::will_accept_transaction(request_t m_request) const {

    return dramsim->WillAcceptTransaction(m_request.address, m_request.is_write);
}

void memory_controller_t::add_transaction(request_t m_request) {
    dramsim->AddTransaction(m_request.address, m_request.is_write);
    request_queue.push_back(m_request);
    request_counter++;

    for(auto it = pending_queue.cbegin(); it != pending_queue.cend(); ++it) {
        if(it->address == m_request.address) {
            pending_queue.erase(it);
            return ;
        }
    }
}

void memory_controller_t::print_stats() const {
    dramsim->PrintStats();
}
