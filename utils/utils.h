#ifndef __UTILS_H__
#define __UTILS_H__

//#include <vector>
#include <list>
#include <string>
#include "data-def.h"

// Convert string to lowercase.
std::string& lowercase(std::string &m_str);

// Convert string to uppercase.
std::string& uppercase(std::string &m_str);

void quantization(data_t *m_dest, const float* m_source, unsigned m_size);

template <typename T>
void move_front(std::list<T> *m_queue) {
    m_queue->push_back(m_queue->front());
    m_queue->pop_front();
}

#endif
