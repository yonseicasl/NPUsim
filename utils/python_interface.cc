#include "python_interface.h"

#ifdef Pytorch
std::vector<std::string> python_list_to_vector(PyObject *pList) {

    std::vector<std::string> layer_list;

    if(PyList_Check(pList)) {
        Py_ssize_t size = PyList_Size(pList);
        layer_list.reserve(size);

        for(Py_ssize_t i = 0; i < size; i++) {
            PyObject* pItem = PyList_GetItem(pList, i);
            if(PyUnicode_Check(pItem)) {
                layer_list.push_back(static_cast<std::string>(PyUnicode_AsUTF8(pItem)));
            }
        }
    }
    return layer_list;
}

std::vector<PyObject*> python_list_to_pyobject(PyObject *pList) {

    std::vector<PyObject*> layer_list;

    if(PyList_Check(pList)) {
        Py_ssize_t size = PyList_Size(pList);
        layer_list.reserve(size);

        for(Py_ssize_t i = 0; i < size; i++) {
            PyObject* pItem = PyList_GetItem(pList, i);
            if(PyUnicode_Check(pItem)) {
                layer_list.push_back(pItem);
            }
        }
    }
    return layer_list;
}


#endif
