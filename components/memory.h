#ifndef __NPU_MEMORY_H__
#define __NPU_MEMORY_H__

#include <cstdlib>
#include <cstdint>
#include <set>

class npu_segment {

public:
    //npu_segment(uint64_t m_ptr, uint64_t m_addr, size_t m_size);
    npu_segment();
    ~npu_segment();

    uint64_t ptr;           // Virtual memory address
    uint64_t addr;          // Physical memory address
    bool valid;         
};

class npu_mmu {
public:

    npu_mmu();
    ~npu_mmu();

    static void init();
    static void npu_malloc(uint64_t m_ptr);
    static void npu_free(uint64_t m_ptr);

    // Virtual to physical address translation.
    static uint64_t v2p(uint64_t m_ptr);            

private:

    static npu_mmu *memory;
    static npu_segment base_addr;
};

#endif
