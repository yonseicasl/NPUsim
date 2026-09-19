#ifndef __UTILS_H__
#define __UTILS_H__

#include <cstddef>
#include <list>
#include <string>
#include "user-def.h"

// Convert string to lowercase.
std::string& lowercase(std::string &m_str);

// Ceiling division; returns 0 when the divisor is 0.
inline size_t ceil_div(size_t m_numerator, size_t m_denominator) {
    return m_denominator ? (m_numerator + m_denominator - 1) / m_denominator : 0;
}

template <typename T>
void move_front(std::list<T> *m_queue) {
    m_queue->push_back(m_queue->front());
    m_queue->pop_front();
}

#endif
