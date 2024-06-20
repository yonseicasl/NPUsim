#include "dnn_model.h"

layer_t::layer_t(layer_name_t m_layer_type) :
    batch_size(1),
    input_height(1),
    input_width(1),
    input_channel(1),
    input_size(1),
    weight_height(1),
    weight_width(1),
    weight_size(1),
    output_height(1),
    output_width(1),
    output_channel(1),
    output_size(1),
    layer_type(m_layer_type),
    Pylayer(NULL),
    input_data(NULL),
    weight(NULL), 
    output_data(NULL) {

}

layer_t::~layer_t() {
    delete [] input_data;
    delete [] weight;
    delete [] output_data;

}


network_t::network_t() :
    Pynetwork(NULL),
    Pylayers(NULL),
    Pydataset(NULL),
    Pyimage(NULL),
    Pylabel(NULL),
    num_layers(0),
    input_data(NULL) {
}

network_t::~network_t() {
    delete [] input_data;
}

void network_t::init(PyObject* pModule, const std::string m_network_config) {
#ifdef Pytorch
    if(pModule) {
        PyObject *pFunc, *pArgs, *pValue;
        pFunc = PyObject_GetAttrString(pModule, "init");
      
        char *t_network_config = const_cast<char*>(m_network_config.c_str());

        std::vector<std::string> DNN_layers_name;
        // Produce arguments and pass to PyTorch.
        pArgs = PyTuple_Pack(1, PyUnicode_FromString(t_network_config));
        if(pFunc) { 
            pValue = PyObject_CallObject(pFunc, pArgs);
            if(pValue) {
                Pynetwork = PyTuple_GetItem(pValue, 0);
                Pylayers  = PyTuple_GetItem(pValue, 1);
                Pylayers_name = PyTuple_GetItem(pValue, 2);
                DNN_layers_name = python_list_to_vector(Pylayers_name);
            }
        }
        layers.reserve(DNN_layers_name.size());
        for(unsigned i = 0; i < DNN_layers_name.size(); i++) {
            layer_t *layer = NULL;
            if(DNN_layers_name[i] == "Conv2d") {
                layer = new layer_t(layer_name_t::CONVOLUTIONAL_LAYER);
            } else if(DNN_layers_name[i] == "MaxPool2d") {
                layer = new layer_t(layer_name_t::MAXPOOL_LAYER);
            } else if(DNN_layers_name[i] == "Linear") {
                layer = new layer_t(layer_name_t::CONNECTED_LAYER);
            } else if(DNN_layers_name[i] == "AdaptiveAvgPool2d") {
                layer = new layer_t(layer_name_t::AVGPOOL_LAYER);
            } else {
                layer = new layer_t(layer_name_t::UNDEFINED_LAYER);
            } 
            layers.push_back(layer);
            num_layers++; 
        }
    }
#endif
}

void network_t::init_layer(PyObject *pModule) {
#ifdef Pytorch

#endif

}

void network_t::init_weight(PyObject *pModule) {
#ifdef Pytorch
    if(pModule) {
        PyObject *pFunc, *pArgs, *pValue;
        pFunc = PyObject_GetAttrString(pModule, "init_weight");

        for(unsigned i = 0; i < num_layers; i++) {
            if(layers[i]->layer_type == layer_name_t::CONVOLUTIONAL_LAYER || 
               layers[i]->layer_type == layer_name_t::CONNECTED_LAYER) {
                pArgs = PyTuple_Pack(2, Pylayers, PyLong_FromLong(i));
                if(pFunc) {
                    pValue = PyObject_CallObject(pFunc, pArgs);
                    if(pValue) {
                        layers[i]->weight = pValue;
                    }
                }
            }
        }
    }
#endif

}

void network_t::load_data(PyObject *pModule, const std::string m_network_config, unsigned m_iteration) {

#ifdef Pytorch
    if(pModule) {
        PyObject *pFunc, *pArgs, *pValue;
        pFunc = PyObject_GetAttrString(pModule, "load_data");
      
        char *t_network_config = const_cast<char*>(m_network_config.c_str());

        // Produce arguments and pass to PyTorch.
        pArgs = PyTuple_Pack(2, PyUnicode_FromString(t_network_config), PyLong_FromLong(m_iteration));
        if(pFunc) { 
            pValue = PyObject_CallObject(pFunc, pArgs);
            if(pValue) {
                Pyimage = PyTuple_GetItem(pValue, 0);
                Pylabel = PyTuple_GetItem(pValue, 1);
            }
        }
    }
#endif
}

void network_t::forward(PyObject *pModule, unsigned m_iteration, unsigned m_index) {

#ifdef Pytorch
    if(pModule) {
        PyObject *pFunc, *pArgs, *pValue;
        pFunc = PyObject_GetAttrString(pModule, "layerwise_forward");

        PyObject *input_data = m_index > 0 ? layers[m_index-1]->output_data : Pyimage;
        // Produce arguments and pass to PyTorch.
        pArgs = PyTuple_Pack(4, Pylayers, Pylayers_name, input_data, PyLong_FromLong(m_index));

        if(pFunc) { 
            pValue = PyObject_CallObject(pFunc, pArgs);
            if(pValue) {
                layers[m_index]->output_data = pValue;
            }
        }
    }
#endif

}

void network_t::print_result(PyObject *pModule) {
#ifdef Pytorch
    if(pModule) {
        PyObject *pFunc, *pArgs, *pValue;
        pFunc = PyObject_GetAttrString(pModule, "print_result");
        
        pArgs = PyTuple_Pack(2, layers[num_layers-1]->output_data, Pylabel);

        if(pFunc) {
            pValue = PyObject_CallObject(pFunc, pArgs);
        }
    }
#endif
}
