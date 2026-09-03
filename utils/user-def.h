#ifndef __USER_DEF_H__ 
#define __USER_DEF_H__

#ifdef USER_INTEGER
    #define DATA_BIT 3
#elif defined(USER_FLOAT)
    #define DATA_BIT 16
    #define SIGN_BIT 1
    #define EXPONENT 5

#else
    // <cstdint> is the header that actually declares uint8_t. Relying on <iostream> to
    // pull it in transitively compiles on some libstdc++ versions and breaks on others.
    #include <cstdint>
    typedef std::uint8_t data_t;
#endif
#endif
