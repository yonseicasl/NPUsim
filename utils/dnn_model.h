#ifndef __NETWORK_H__
#define __NETWORK_H__

#include "def.h"
#include "user-def.h"
#include "python_interface.h"

class layer_t {
public:

    layer_t();
    layer_t(layer_name_t m_layer_type);
    ~layer_t();


    layer_name_t layer_type;

    data_t *input_data;
    data_t *weight;
    data_t *output_data;

};

class network_t {
public:
    network_t();
    ~network_t();

    void init(PyObject *pModule, const std::string m_network_config);

    void load_data(PyObject *pModule, const std::string m_network_config, unsigned m_iteration);

    void forward(PyObject *pModule, unsigned m_index);


    std::vector<layer_t*> layers; 
    PyObject *Pynetwork;
    unsigned num_layers;

    data_t *input_data;

};


#endif
