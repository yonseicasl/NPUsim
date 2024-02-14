#ifndef __MEMORY_CONTROLLER_H__
#define __MEMORY_CONTROLLER_H__

#include <list>
#include "dramsim3.h"
#include "mapping_table.h"
#include "def.h"
#include "data-def.h"

class request_t {
public:
    request_t() : address(0), is_write(0) {}
    request_t(uint64_t m_address, bool m_is_write) : address(m_address), is_write(m_is_write) {}
    uint64_t address;
    bool is_write;
};

class memory_controller_t {

public:
    memory_controller_t(const std::string &m_config_file, const std::string &m_output_dir);
                
    ~memory_controller_t();

    void send_request(data_t *m_data, mapping_table_t *m_mapping_table, unsigned m_offset, data_type_t m_data_type, action_type_t m_action_type);

    uint64_t v2p(uint64_t m_ptr);

    void insert_pending_queue(request_t m_request);

    void tick();

    void read_callback(uint64_t m_address);
    void write_callback(uint64_t m_address);

    void register_callback(std::function<void(uint64_t)> m_read_callback,
                           std::function<void(uint64_t)> m_write_callback);
    

    unsigned int get_bus_bits() const;

    unsigned int get_burst_length() const;

    void print_stats() const;

    void reset_stats();

    // Check whether the controller can accept a new packet or not.
    bool will_accept_transaction(request_t m_request) const;
    // when will_accept_transaction() function returns true,
    // Enqueue the packet.
    void add_transaction(request_t m_request);

    bool done();


private: 
    dramsim3::MemorySystem *dramsim;

    std::list<request_t> pending_queue;
    std::list<request_t> request_queue;

    unsigned callback_counter;
    unsigned request_counter;
};

#endif
