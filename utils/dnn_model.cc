#include "dnn_model.h"

layer_t::layer_t(layer_name_t m_layer_type) :
    layer_type(m_layer_type),
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
    Pydataset(NULL),
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
        //PyObject *pItem1, *pItem2;
        PyObject *pLayer_name;
        pFunc = PyObject_GetAttrString(pModule, "init");
      
        char *t_network_config = const_cast<char*>(m_network_config.c_str());

        std::vector<std::string> DNN_layers_name;
        // Produce arguments and pass to PyTorch.
        pArgs = PyTuple_Pack(1, PyUnicode_FromString(t_network_config));
        if(pFunc) { 
            pValue = PyObject_CallObject(pFunc, pArgs);
            if(pValue) {
                Pynetwork = PyTuple_GetItem(pValue, 0);
                //pItem1 = PyTuple_GetItem(pValue, 1);
                pLayer_name = PyTuple_GetItem(pValue, 2);
                DNN_layers_name = python_list_to_vector(pLayer_name);
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

void network_t::load_data(PyObject *pModule, const std::string m_network_config) {

#ifdef Pytorch
    if(pModule) {
        PyObject *pFunc, *pArgs, *pValue;
        pFunc = PyObject_GetAttrString(pModule, "load_data");
      
        char *t_network_config = const_cast<char*>(m_network_config.c_str());

        // Produce arguments and pass to PyTorch.
        pArgs = PyTuple_Pack(1, PyUnicode_FromString(t_network_config));
        if(pFunc) { 
            pValue = PyObject_CallObject(pFunc, pArgs);
            if(pValue) {
                Pydataset = pValue;
            }
        }
    }
#endif
}

void network_t::forward(PyObject *pModule, unsigned m_iteration, unsigned m_index) {

#ifdef Pytorch
    if(m_index == num_layers-1) {
        if(pModule) {
            PyObject *pFunc, *pArgs, *pValue;
            pFunc = PyObject_GetAttrString(pModule, "forward");

            // Produce arguments and pass to PyTorch.
            pArgs = PyTuple_Pack(4, Pynetwork, Pydataset, PyLong_FromLong(m_iteration), PyLong_FromLong(m_index));

            if(pFunc) { 
                pValue = PyObject_CallObject(pFunc, pArgs);
                if(pValue) {
                    layers[m_index]->output_data = pValue;
                }
            }
        }
    }
#endif

}

