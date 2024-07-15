#ifndef __PYTHON_INTERFACE_H__
#define __PYTHON_INTERFACE_H__

#ifdef PyTorch

#include <Python.h>
#include <vector>
#include <string>
#include <iostream>

#include "user-def.h"

std::vector<std::string> python_list_to_vector(PyObject *pList);

data_t* python_tensor_to_data(PyObject *pList);
#endif
#endif
