// #ifndef PYTHONFIELDS_H_
// #define PYTHONFIELDS_H_

// #define PY_SSIZE_T_CLEAN
// #include <Python.h>

#include "field.h"
#include <vector>
#include <iostream>
#include "boost/format.hpp"
#include "pythonFields.h"



TPythonField::TPythonField(const std::string ft){
  
  // PyRun_SimpleString("import sys"); 
  // PyRun_SimpleString("sys.path.append(\".\")");

  // PyObject* magpymodule = PyImport_ImportModule("magpylib");
  // pBFieldFunc = PyObject_GetAttrString(magpymodule, "getB");

  std::cout << " pythonField function loading" << std::endl;

  // boost::python::object bmagpymodule = boost::python::import("magpylib");
  // bpBFieldFunc = bmagpymodule.attr("getB");
  boost::python::object bmagpymodule = boost::python::import("pythonField");
  bpBFieldFunc = bmagpymodule.attr("BField");

  std::cout << " python magnetic source building" << std::endl;
  
  // char* ftc = const_cast<char*>(ft.c_str());
  boost::python::str bftc(ft);
  boost::python::object bmagnetmodule = boost::python::import(bftc);
  boost::python::object buildSourceFunc = bmagnetmodule.attr("buildSource");
  bpSourceObject = buildSourceFunc();
  // bpSourceObject = boost::python::call_method<boost::python::object>(, NULL);
  
  // PyObject* magnetmodule = PyImport_ImportModule(ftc);
  // PyObject* pMagnetFunc = PyObject_GetAttrString(magnetmodule, "BuildMagnet");
  // pMagnetObject = PyObject_CallObject(pMagnetFunc, NULL);
  // Py_DECREF(magnetmodule);
  // Py_DECREF(pMagnetFunc);
      
}

void TPythonField::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{

  const double h = 1e-6;
  int dimarg = (dBidxj != nullptr) ? 7 : 1;
  // PyObject *pArgs = PyTuple_New(2);
  // PyObject *pList = PyList_New(dimarg);
  // PyObject *pxyzList;
  // PyArrayObject* npArray;
  double xyz[3];

  boost::python::list bpList;    

  // std::cout << " building list args" << std::endl;

  for(int i=0; i<dimarg; ++i){

    xyz[0] = x;  // swap x and y. Todo : sly
    xyz[1] = y;
    xyz[2] = z;  // swap x and y. Todo : sly
    
    switch(i) {
    case 0:
      break;
    case 1:
      xyz[0] -= h;
      break;
    case 2:
      xyz[0] += h;
      break;
    case 3:
      xyz[1] -= h;
      break;
    case 4:
      xyz[1] += h;
      break;
    case 5:
      xyz[2] -= h;
      break;
    case 6:
      xyz[2] += h;
      break;
    default:
      printf("out of case");
    }

    // pxyzList = PyList_New(3);
    boost::python::list bpxyzList;

    for(int j=0; j<3; ++j){
      // PyList_SetItem(pxyzList, j, PyFloat_FromDouble(1000 * xyz[j]));
      bpxyzList.append(1000 * xyz[j]);
    }
    // PyList_SetItem(pList, i, pxyzList);
    bpList.append(bpxyzList);
  }

  // PyTuple_SetItem(pArgs, 0, pMagnetObject);
  // PyTuple_SetItem(pArgs, 1, pList);

  // std::cout << "calling getBs" << std::endl;

  // PyGILState_STATE gstate = PyGILState_Ensure();
  // PyThreadState* gstate = PyEval_SaveThread();

  
    
  double Bs[dimarg][3];
  boost::python::tuple args;
  boost::python::object bnpArray;
  // npArray = reinterpret_cast<PyArrayObject*>(PyObject_CallObject(pBFieldFunc, pArgs));
  try
    {
      args = boost::python::make_tuple(bpSourceObject, bpList);
      bnpArray = bpBFieldFunc(bpSourceObject, bpList);
    }
  catch (const std::exception& e)
    {
        // Catch and handle C++ exceptions
      std::cout << "Caught C++ exception: " << e.what() << std::endl;
    }
  catch (...)
    {
      // Catch and handle any other C++ exceptions
      std::cout << "Caught unknown C++ exception" << std::endl;
    }


  // Check if a Python exception occurred
  if (PyErr_Occurred())
    {
      // Fetch the Python exception type, value, and traceback
      PyObject* pType, * pValue, * pTraceback;
      PyErr_Fetch(&pType, &pValue, &pTraceback);
      
      // Convert the Python objects to strings
      PyObject* pTypeStr = PyObject_Str(pType);
      PyObject* pValueStr = PyObject_Str(pValue);
      PyObject* pTracebackStr = PyObject_Str(pTraceback);
      
        // Get the C-style string representation of the Python objects
      const char* pTypeCStr = PyUnicode_AsUTF8(pTypeStr);
      const char* pValueCStr = PyUnicode_AsUTF8(pValueStr);
      const char* pTracebackCStr = PyUnicode_AsUTF8(pTracebackStr);
      
      // Print the Python exception information
      std::cout << "Caught Python exception:" << std::endl;
      std::cout << "Type: " << pTypeCStr << std::endl;
      std::cout << "Value: " << pValueCStr << std::endl;
      std::cout << "Traceback: " << pTracebackCStr << std::endl;
      
      // Clean up the Python objects
      Py_XDECREF(pType);
      Py_XDECREF(pValue);
      Py_XDECREF(pTraceback);
      Py_XDECREF(pTypeStr);
      Py_XDECREF(pValueStr);
      Py_XDECREF(pTracebackStr);
    }
  
  
  // PyGILState_Release(gstate);
  // PyEval_RestoreThread(gstate);
  
  // std::cout << "getting Bs" << std::endl;

  if (dBidxj != nullptr){
    for (int i=0; i<dimarg; ++i) {
      for (int j=0; j<3; ++j) { 
	// Bs[i][j] = 0.001 * (*reinterpret_cast<double*>(PyArray_GETPTR2(npArray, i, j)));
	Bs[i][j] = 0.001 * boost::python::extract<double>(bnpArray[i][j]);
      }
    }
  }
  else{
    for (int j=0; j<3; ++j) {
      Bs[0][j] = 0.001 * boost::python::extract<double>(bnpArray[j]);
    }
  }


  for (int i=0; i<3; ++i) {
    B[i] = Bs[0][i];
  }

  
  // std::cout << "Bi = " << B[0] << ", " << B[1] << ", " << B[2] << std::endl;

  if (dBidxj != nullptr){

    double trace_3;
    double dBi_dxj[3][3];

    for(int i=0; i<3; ++i){
      dBi_dxj[i][0] = Bs[2][i]/(2*h) - Bs[1][i]/(2*h);
      dBi_dxj[i][1] = Bs[4][i]/(2*h) - Bs[3][i]/(2*h);
      dBi_dxj[i][2] = Bs[6][i]/(2*h) - Bs[5][i]/(2*h);
      
    }
    
    trace_3 = (dBi_dxj[0][0] + dBi_dxj[1][1] + dBi_dxj[2][2])/3;
    // std::cout << trace_3 << std::endl;

    
    dBidxj[0][0] = dBi_dxj[0][0] - trace_3;
    dBidxj[1][1] = dBi_dxj[1][1] - trace_3;
    dBidxj[2][2] = dBi_dxj[2][2] - trace_3;

    dBidxj[0][1] = dBidxj[1][0] = (dBi_dxj[0][1] + dBi_dxj[1][0])/2;
    dBidxj[0][2] = dBidxj[2][0] = (dBi_dxj[0][2] + dBi_dxj[2][0])/2;
    dBidxj[1][2] = dBidxj[2][1] = (dBi_dxj[1][2] + dBi_dxj[2][1])/2;

    // std::cout << "dBi_dxj = \n " << std::endl;
    // std::cout << dBidxj[0][0] << ", " << dBidxj[0][1] << ", " << dBidxj[0][2] << std::endl;
    // std::cout << dBidxj[1][0] << ", " << dBidxj[1][1] << ", " << dBidxj[1][2] << std::endl;
    // std::cout << dBidxj[2][0] << ", " << dBidxj[2][1] << ", " << dBidxj[2][2] << "\n" << std::endl;

  }

  // Py_DECREF(npArray);
  // Py_DECREF(pList);
  // PyTuple_SET_ITEM(pArgs, 0, pxyzList);
  // Py_DECREF(pxyzList);
  // Py_DECREF(pArgs);

}
