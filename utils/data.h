#ifndef __DATA_H__
#define __DATA_H__

#include "user-def.h"

#ifdef USER_INTEGER
struct data_t {
    int value : DATA_BIT;
};

#elif defined(USER_FLOAT)
struct float_bits {
    unsigned fraction: DATA_BIT - SIGN_BIT - EXPONENT;
    unsigned exponent:EXPONENT;
    unsigned sign:SIGN_BIT;
};

union data_t {
    float value;
    float_bits bits;
};

#endif
inline bool data_is_nonzero(const data_t &data) {
#ifdef USER_INTEGER
    return data.value != 0;
#elif defined(USER_FLOAT)
    return data.value != 0.0f;
#else
    return data != 0;
#endif
}


inline void data_copy(data_t &destination, const data_t &source) {
#if defined(USER_INTEGER) || defined(USER_FLOAT)
    destination.value = source.value;
#else
    destination = source;
#endif
}

inline void data_clear(data_t &data) {
    data = data_t{};
}

inline void data_accumulate_product(data_t &accumulator, const data_t &input, const data_t &weight) {
#if defined(USER_INTEGER) || defined(USER_FLOAT)
    accumulator.value += input.value * weight.value;
#else
    accumulator += input * weight;
#endif
}

inline void data_relu(data_t &data) {
#if defined(USER_INTEGER) || defined(USER_FLOAT)
    data.value = data.value > 0 ? data.value : 0;
#else
    data = data > 0 ? data : 0;
#endif
}

#endif
