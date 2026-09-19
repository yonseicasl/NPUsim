#include <algorithm>
#include <cctype>
#include "utils.h"

std::string& lowercase(std::string &m_str) {
    transform(m_str.begin(), m_str.end(), m_str.begin(), ::tolower);
    return m_str;
}
