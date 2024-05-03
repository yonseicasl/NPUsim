#ifndef __DATA_H__
#define __DATA_H__

#include "user-def.h"

#ifdef USER_INTEGER
struct data_t {
    int value : DATA_BIT;
};

#elif USER_FLOAT
struct float_bits {
    unsigned fraction: DATA_BIT - SIGN_BIT - EXPONENT;
    unsigned exponent:EXPONENT;
    unsigned sign:SIGN_BIT;
};

union data_t {
    float value;
    float_bits bits;
};

#else

#endif

#endif
