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

    unsigned batch_size;
    unsigned input_height;
    unsigned input_width;
    unsigned input_channel;
    unsigned input_size;

    unsigned weight_height;
    unsigned weight_width;
    unsigned weight_size;

    unsigned output_height;
    unsigned output_width;
    unsigned output_channel;
    unsigned output_size;

    layer_name_t layer_type;

    PyObject *Pylayer;
    PyObject *input_data;
    PyObject *weight;
    PyObject *output_data;
};

class network_t {
public:
    network_t();
    ~network_t();

    void init(PyObject *pModule, const std::string m_network_config);

    void init_layer(PyObject *pModule);

    void init_weight(PyObject *pModule);

    void load_data(PyObject *pModule, const std::string m_network_config, unsigned m_iteration);

    void forward(PyObject *pModule, unsigned m_iteration, unsigned m_index);

    void print_result(PyObject *pModule);


    std::vector<layer_t*> layers; 
    PyObject *Pynetwork;
    PyObject *Pylayers;
    PyObject *Pylayers_name;
    PyObject *Pydataset;
    PyObject *Pyimage;
    PyObject *Pylabel;
    unsigned num_layers;

    PyObject *input_data;

};


#endif
