/**
 * \file
 * Static, analytical B fields
 *
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>

#include "boost/format.hpp"

#include "analyticFields.h"


using namespace std;

//TExponentialBFieldX constructor
TExponentialFieldX::TExponentialFieldX(const double _a1, const double _a2, const double _a3, const double _c1, const double _c2){
  a1 = _a1; a2 = _a2; a3 = _a3;
  c1 = _c1; c2 = _c2;
}

void TExponentialFieldX::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
  B[0] = a1 * exp(- a2* x + a3) + c1; //Bx contribution
  B[1] = y * a1 * a2 / 2 * exp(- a2* x + a3) + c2; //By contribution
  B[2] = z* a1 * a2 / 2 * exp(- a2* x + a3) + c2; //Bz contribution
  if (dBidxj != NULL){
    dBidxj[0][0] = -a2 * a1 * exp(- a2* x + a3); // dBxdx
    dBidxj[0][1] = 0; //dBxdy
    dBidxj[0][2] = 0; //dBxdz
    dBidxj[1][0] = -a2 * y * a1 * a2 / 2 * exp(- a2* x + a3); //dBydx
    dBidxj[1][1] = a1 * a2 / 2 * exp(- a2* x + a3); //dBydy
    dBidxj[1][2] = 0; //dBydz
    dBidxj[2][0] = -a2 * z * a1 * a2 / 2 * exp(- a2* x + a3); //dBzdx
    dBidxj[2][1] = 0; //dBzdy
    dBidxj[2][2] = a1 * a2 / 2 * exp(- a2* x + a3); //dBzdz
  }
}


//TLinearFieldZ constructor
TLinearFieldZ::TLinearFieldZ(const double _a1, const double _a2){
  a1 = _a1; a2 = _a2;
}

void TLinearFieldZ::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
  B[0] = 0; //Bx contribution
  B[1] = 0; //By contribution
  B[2] = (a1*x + a2); //Bz contribution
  if (dBidxj != NULL){
    dBidxj[0][0] = 0; // dBxdx
    dBidxj[0][1] = 0; //dBxdy
    dBidxj[0][2] = 0; //dBxdz
    dBidxj[1][0] = 0; //dBydx
    dBidxj[1][1] = 0; //dBydy
    dBidxj[1][2] = 0; //dBydz
    dBidxj[2][0] = a1; //dBzdx
    dBidxj[2][1] = 0; //dBzdy
    dBidxj[2][2] = 0; //dBzdz
  }
}


//TB0GradZ constructor
TB0GradZ::TB0GradZ(const double _a1, const double _a2, const double _z0){
  a1 = _a1; a2 = _a2; z0 = _z0;
}

void TB0GradZ::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
  B[0] = (-(a1*z + a2)/2 * x); //Bx contribution
  B[1] = (-(a1*z + a2)/2 * y); //By contribution
  B[2] = (a1*z*z/2 + a2*z + z0); //Bz contribution
  if (dBidxj != NULL){
    dBidxj[0][0] = (-(a1*z + a2) / 2); // dBxdx
    dBidxj[0][1] = 0; //dBxdy
    dBidxj[0][2] = (- a1*x / 2); //dBxdz
    dBidxj[1][0] = 0; //dBydx
    dBidxj[1][1] = (-(a1*z + a2) / 2); //dBydy
    dBidxj[1][2] = (- a1*y / 2); //dBydz
    dBidxj[2][0] = 0; //dBzdx
    dBidxj[2][1] = 0; //dBzdy
    dBidxj[2][2] = (a1*z + a2); //dBzdz
  }
}


//TB0GradX2 constructor
TB0GradX2::TB0GradX2(const double _a1, const double _a2, const double _a3, const double _z0){
  a1 = _a1; a2 = _a2; a3 = _a3; z0 = _z0;
}

void TB0GradX2::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
  B[0] = ( -a1/6*x*x*x - a2/4*x*x -a3/2*x ); //Bx contribution
  B[1] = ( -( a1*x*x + a2*x + a3 ) / 2 * y ); //By contribution
  B[2] = ((a1*x*x + a2*x + a3) * z + z0); //Bz contribution
  if (dBidxj != NULL){
    dBidxj[0][0] = (-( a1*x*x + a2*x + a3 ) / 2); // dBxdx
    dBidxj[0][1] = 0; //dBxdy
    dBidxj[0][2] = 0; //dBxdz
    dBidxj[1][0] = (-a1*x -a2/2) * y; //dBydx
    dBidxj[1][1] = (-( a1*x*x + a2*x + a3 ) / 2); //dBydy
    dBidxj[1][2] = 0; //dBydz
    dBidxj[2][0] = (2*a1*x + a2)*z; //dBzdx
    dBidxj[2][1] = 0; //dBzdy
    dBidxj[2][2] = ( a1*x*x + a2*x + a3 ); //dBzdz
  }
}


//TB0GradXY constructor
TB0GradXY::TB0GradXY(const double _a1, const double _a2, const double _z0){
  a1 = _a1; a2 = _a2; z0 = _z0;
}

void TB0GradXY::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
  B[0] = (-a1*x*x*y/4 - a2*x/2); //Bx contribution
  B[1] = (-a1*x*y*y/4 - a2*y/2); //By contribution
  B[2] = (a1*x*y*z + a2*z +z0); //Bz contribution
  if (dBidxj != NULL){
    dBidxj[0][0] = (-(a1*x*y +a2) / 2); // dBxdx
    dBidxj[0][1] = (-a1*x*x/4); //dBxdy
    dBidxj[0][2] = 0; //dBxdz
    dBidxj[1][0] = (-a1*y*y/4); //dBydx
    dBidxj[1][1] = (-(a1*x*y +a2) / 2); //dBydy
    dBidxj[1][2] = 0; //dBydz
    dBidxj[2][0] = a1*y*z; //dBzdx
    dBidxj[2][1] = a1*x*z; //dBzdy
    dBidxj[2][2] = (a1*x*y + a2); //dBzdz
  }
}


//TB0_XY constructor
TB0_XY::TB0_XY(const double _a1, const double _z0) {
  a1 = _a1; z0 = _z0;
}

void TB0_XY::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
  B[0] = (a1 * y * z); //Bx contribution
  B[1] = (a1 * x * z); //By contribution
  B[2] = (a1 * x * y + z0); //Bz contribution
  if (dBidxj != NULL){
    dBidxj[0][0] = 0; // dBxdx
    dBidxj[0][1] = (a1 * z); //dBxdy
    dBidxj[0][2] = (a1 * y); //dBxdz
    dBidxj[1][0] = (a1 * z); //dBydx
    dBidxj[1][1] = 0; //dBydy
    dBidxj[1][2] = (a1 * x); //dBydz
    dBidxj[2][0] = (a1 * y); //dBzdx
    dBidxj[2][1] = (a1 * x); //dBzdy
    dBidxj[2][2] = 0; //dBzdz
  }
}


TCustomBField::TCustomBField(const std::string &_Bx, const std::string &_By, const std::string &_Bz){
  tvar = unique_ptr<double>(new double(0.0));
  xvar = unique_ptr<double>(new double(0.0));
  yvar = unique_ptr<double>(new double(0.0));
  zvar = unique_ptr<double>(new double(0.0));
  exprtk::symbol_table<double> symbol_table;
  symbol_table.add_variable("t",*tvar);
  symbol_table.add_variable("x",*xvar);
  symbol_table.add_variable("y",*yvar);
  symbol_table.add_variable("z",*zvar);
  symbol_table.add_constants();
  exprtk::parser<double> parser;

  std::array<std::string, 3> expr{_Bx, _By, _Bz};
  for (int i = 0; i < 3; ++i){
    Bexpr[i].register_symbol_table(symbol_table);
    if (not parser.compile(expr[i], Bexpr[i])){
      throw std::runtime_error(exprtk::parser_error::to_str(parser.get_error(0).mode) + " while parsing CustomBField formula '" + expr[i] + "': " + parser.get_error(0).diagnostic);
    }
  }
}

void TCustomBField::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
  *xvar = x;
  *yvar = y;
  *zvar = z;
  *tvar = t;
  B[0] = Bexpr[0].value();
  B[1] = Bexpr[1].value();
  B[2] = Bexpr[2].value();
  //	std::cout << B[0] << " " << B[1] << " " << B[2] << " ";
	
  if (dBidxj != nullptr){
    dBidxj[0][0] = exprtk::derivative(Bexpr[0], *xvar);
    dBidxj[0][1] = exprtk::derivative(Bexpr[0], *yvar);
    dBidxj[0][2] = exprtk::derivative(Bexpr[0], *zvar);
    dBidxj[1][0] = exprtk::derivative(Bexpr[1], *xvar);
    dBidxj[1][1] = exprtk::derivative(Bexpr[1], *yvar);
    dBidxj[1][2] = exprtk::derivative(Bexpr[1], *zvar);
    dBidxj[2][0] = exprtk::derivative(Bexpr[2], *xvar);
    dBidxj[2][1] = exprtk::derivative(Bexpr[2], *yvar);
    dBidxj[2][2] = exprtk::derivative(Bexpr[2], *zvar);
    //		std::cout << dBidxj[0][0] << " " << dBidxj[0][1] << " " << dBidxj[0][2] << std::endl;
  }
  //	std::cout << std::endl;

}




/**
 * Load python a function from a python module
 *
 * @param v_r Radial component of vector
 * @param v_phi Azimuthal component of vector
 * @param phi Azimuth of vector origin
 * @param v_x Returns x component of vector
 * @param v_y Returns y component of vector
 */
void GetpFunc(const char *moduleName, const char *funcName, PyObject **pFunc){

  PyObject *pModuleName, *pModule;
  pModuleName = PyUnicode_DecodeFSDefault(moduleName);
  /* TODO Error checking of pName left out */

  pModule = PyImport_Import(pModuleName);
  Py_DECREF(pModuleName);

  if (pModule != NULL) {
    *pFunc = PyObject_GetAttrString(pModule, funcName);
    
    if (!*pFunc && !PyCallable_Check(*pFunc)) {
      if (PyErr_Occurred()){
	PyErr_Print();
	fprintf(stderr, "Failed to load \"%s\"\n", funcName); // TODO : duplicate
	Py_DECREF(pModule);
	Py_DECREF(pModuleName);
	throw std::runtime_error((boost::format("Failed to load python function %1%") % funcName).str());	  
      }
    }
    else{
      fprintf(stderr, "Sucess to load \"%s\", is callable \n", funcName); // TODO : duplicate
    }
  }
  else {
    PyErr_Print();
    fprintf(stderr, "Failed to load \"%s\"\n", moduleName);
    Py_DECREF(pModule);
    Py_DECREF(pModuleName);
    throw std::runtime_error((boost::format("Failed to load python module %1%") % moduleName).str());	  
  }
  
  Py_DECREF(pModule);
  Py_DECREF(pModuleName);
  
}



TMagpy::TMagpy(const std::string ft){
  tvar = unique_ptr<double>(new double(0.0));
  xvar = unique_ptr<double>(new double(0.0));
  yvar = unique_ptr<double>(new double(0.0));
  zvar = unique_ptr<double>(new double(0.0));

  char* ftc = const_cast<char*>(ft.c_str());
  PyObject *pMagnetFunc;

  std::cout << "TMapy : initializing python" << std::endl;
  
  Py_Initialize();

  // PyRun_SimpleString("import sys"); 
  // PyRun_SimpleString("sys.path.append(\".\")");

  // std::cout << "TMapy : Get function BuildMagnet" << std::endl;
  GetpFunc(ftc, "BuildMagnet", &pMagnetFunc);

  // std::cout << "TMapy : get function BField" << std::endl;
  GetpFunc("MagpyField", "BField", &pBFieldFunc);
  
  // std::cout << "TMapy : two functions are loaded." << std::endl;

  pMagnetObject = PyObject_CallObject(pMagnetFunc, NULL);
  // std::cout << "TMapy : magnet loaded." << std::endl;
      
  Py_DECREF(pMagnetFunc); // TODO : comment or not for multiple call?
      
}


void TMagpy::GetpBField(double B[3], const double x, const double y, const double z, const double t) const {

  PyObject *pArgs, *pBfieldValue;

  // std::cout << "GetpBField" << std::endl;

  pArgs = PyTuple_New(4);
  PyTuple_SetItem(pArgs, 0, pMagnetObject);
  PyTuple_SetItem(pArgs, 1, PyFloat_FromDouble(x));
  PyTuple_SetItem(pArgs, 2, PyFloat_FromDouble(y));
  PyTuple_SetItem(pArgs, 3, PyFloat_FromDouble(z));
  
  pBfieldValue = PyObject_CallObject(pBFieldFunc, pArgs);

  // std::cout << "Looping for GetpBField" << std::endl;

  if (pBfieldValue != NULL) {
    for (int i=0; i<3; i++) {
      B[i] = PyFloat_AsDouble(PyList_GetItem(pBfieldValue, i));
      // printf("B[%d] = %f, ", i, B[i]);
    }
    // printf("\n");
  }

  // Py_DECREF(pBfieldValue); TODO : why cause issue next time same function is called?

}


void TMagpy::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const{
  *xvar = x;
  *yvar = y;
  *zvar = z;
  *tvar = t;

  // PyObject *pArgs, *pBfieldValue;
  
  // pArgs = PyTuple_New(4);
  // PyTuple_SetItem(pArgs, 0, pMagnetObject);
  // PyTuple_SetItem(pArgs, 1, PyFloat_FromDouble(x));
  // PyTuple_SetItem(pArgs, 2, PyFloat_FromDouble(y));
  // PyTuple_SetItem(pArgs, 3, PyFloat_FromDouble(z));
 
  // pBfieldValue = PyObject_CallObject(pBFieldFunc, pArgs);

  // if (pBfieldValue != NULL) {
  //   for (int i=0; i<3; i++) {
  //     B[i] = PyFloat_AsDouble(PyList_GetItem(pBfieldValue, i));
  //     printf("B[%d] = %f, ", i, B[i]);
  //   }
  //   printf("\n");
  //   // Py_DECREF(pBfieldValue); TODO : why cause issue next time same function is called?
  // }

  // std::cout << "computing B field" << std::endl;

  TMagpy::GetpBField(B, x, y, z, t);

  // std::cout << "Bi = " << B[0] << ", " << B[1] << ", " << B[2] << "\n" << std::endl;

  // std::cout << "computing dBidxj" << std::endl;

  if (dBidxj != nullptr){

    // std::cout << "inside dBidxj" << std::endl;

    const double h = 1e-6;
    double trace_3;

    double B_xph[3];
    double B_yph[3];
    double B_zph[3];
    double B_xmh[3];
    double B_ymh[3];
    double B_zmh[3];

    double dBi_dxj[3][3];
    
    TMagpy::GetpBField(B_xph, x+h, y, z, t);
    TMagpy::GetpBField(B_xmh, x-h, y, z, t);
    TMagpy::GetpBField(B_yph, x, y+h, z, t);
    TMagpy::GetpBField(B_ymh, x, y-h, z, t);
    TMagpy::GetpBField(B_zph, x, y, z+h, t);
    TMagpy::GetpBField(B_zmh, x, y, z-h, t);

    // std::cout << "Computing dBi_dxj" << std::endl;

    for(int i=0; i<3; ++i){
      dBi_dxj[i][0] = (B_xph[i] - B_xmh[i])/(2*h);
      dBi_dxj[i][1] = (B_yph[i] - B_ymh[i])/(2*h);
      dBi_dxj[i][2] = (B_zph[i] - B_zmh[i])/(2*h);
    }
    
    trace_3 = (dBi_dxj[0][0] + dBi_dxj[1][1] + dBi_dxj[2][2])/3;
    
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

 

}
