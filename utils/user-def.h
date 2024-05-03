#ifndef __USER_DEF_H__ 
#define __USER_DEF_H__

#ifdef USER_INTEGER
    #define DATA_BIT 3
#elif USER_FLOAT
    #define DATA_BIT 16
    #define SIGN_BIT 1
    #define EXPONENT 5

#else
    #include <iostream>
    typedef uint8_t data_t;
#endif
#endif
