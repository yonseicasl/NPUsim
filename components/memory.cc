#include <iostream>

#include "memory.h"

#define BASE_ADDR 0X80000000
#define CAPACITY 0x200000000

npu_mmu *npu_mmu::memory=0;
npu_segment npu_mmu::base_addr;

npu_segment::npu_segment() :
    ptr(0),
    addr(BASE_ADDR),
    valid(false) {

}

npu_segment::~npu_segment() {

}

npu_mmu::npu_mmu() {
}

npu_mmu::~npu_mmu() {
    delete memory;
}

void npu_mmu::init() {
    if(!memory) {
        memory = new npu_mmu();
    }
}

// Setting virtual memory address and physical memory address.
void npu_mmu::npu_malloc(uint64_t m_ptr) {
    // Set the smallest virtual memory address as the baseline virtual memory address.
    if(base_addr.valid) {
        if(m_ptr < base_addr.ptr) {
            base_addr.ptr = m_ptr;
        }
    }
    else {
        base_addr.ptr = m_ptr, base_addr.addr = BASE_ADDR, base_addr.valid = true;
    }
}


void npu_mmu::npu_free(uint64_t m_ptr) {
}

uint64_t npu_mmu::v2p(uint64_t m_ptr) {
    
    uint64_t offset;

    offset = (m_ptr - base_addr.ptr) % CAPACITY;
    

    if(m_ptr >= base_addr.ptr && base_addr.addr + offset < CAPACITY) {
        //std::cout << "Virtual address : 0x" << std::hex << m_ptr << std::endl;
        //std::cout << "Physical address : 0x" << std::hex << base_addr.addr + offset << std::endl;
        return base_addr.addr + offset;
    }
    else if(m_ptr < base_addr.ptr) {
        std::cerr << "Virtual address : 0x" << std::hex << m_ptr << " is smaller than the baseline virtual address : 0x" << std::hex << base_addr.ptr << std::endl;
        exit(1);
    }
    else if(offset > CAPACITY) {
        std::cout << "Physical address : 0x" << std::hex << base_addr.addr+offset << " is bigger than the capacity 0x" << CAPACITY << std::endl;
        std::cerr << "Physical address : 0x" << std::hex << base_addr.addr + offset << " is bigger than the capacity 0x" << CAPACITY << std::endl;
        exit(1);
    }
    return 0;
}
