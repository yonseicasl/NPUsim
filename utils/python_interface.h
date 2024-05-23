#ifndef __PYTHON_INTERFACE_H__
#define __PYTHON_INTERFACE_H__



#ifdef Pytorch

#include <Python.h>
#include <vector>
#include <string>
#include <iostream>


std::vector<std::string> python_list_to_vector(PyObject *pList);


#endif
#endif
