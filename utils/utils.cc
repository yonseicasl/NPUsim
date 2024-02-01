#include <algorithm>
#include <cctype>
#include <cmath>
#include <vector>
#include <iostream>
#include "utils.h"

std::string& lowercase(std::string &m_str) {
    transform(m_str.begin(), m_str.end(), m_str.begin(), ::tolower);
    return m_str;
}

std::string& uppercase(std::string &m_str) {
    transform(m_str.begin(), m_str.end(), m_str.begin(), ::toupper);
    return m_str;
}

void quantization(data_t *m_dest, const float* m_source, unsigned m_size) {

    std::vector<float> t_vector_float = {m_source, m_source + m_size};
    float float_min = *std::min_element(std::begin(t_vector_float), std::end(t_vector_float));
    float float_max = *std::max_element(std::begin(t_vector_float), std::end(t_vector_float));

    float ret_min = std::numeric_limits<data_t>::min();
    float ret_max = std::numeric_limits<data_t>::max();

    float scale = (float_max - float_min) / (ret_max - ret_min);
    float zero_point = floor((float_max*ret_max - float_min*ret_min)/(float_max - float_min) + 0.5);

    for(unsigned i = 0; i < t_vector_float.size(); i++) {
        m_dest[i] = (floor((m_source[i] / scale) + 0.5) + zero_point);
    }
}
