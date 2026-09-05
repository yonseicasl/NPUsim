#ifndef __USER_DEF_H__ 
#define __USER_DEF_H__

#ifdef USER_INTEGER
    #define DATA_BIT 3
#elif defined(USER_FLOAT)
    #define DATA_BIT 16
    #define SIGN_BIT 1
    #define EXPONENT 5

#elif defined(FUNCTIONAL)
    // Functional simulation moves and computes on Nebula's tensors, and Nebula stores
    // float. The former uint8_t default made every `(data_t*)layer->...` cast REINTERPRET
    // the float arrays byte-by-byte -- indices covered a quarter of each tensor and the
    // values were garbage, so functional results were silently meaningless. data_t must
    // match the framework element type unless a quantized mode (USER_INTEGER/USER_FLOAT)
    // explicitly redefines it.
    typedef float data_t;
#else
    // <cstdint> is the header that actually declares uint8_t. Relying on <iostream> to
    // pull it in transitively compiles on some libstdc++ versions and breaks on others.
    #include <cstdint>
    typedef std::uint8_t data_t;
#endif
#endif
