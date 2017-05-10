// DO NOT EDIT
// Generated automatically using libclang using the command :
// c++2py.py ../c++/pomerol_ed.hpp -I../../../pomerol/installed/include -I/usr/include/eigen3 -I../c++ --only_converters -p -mpytriqs.applications.impurity_solvers.pomerol2triqs -o pomerol2triqs --moduledoc "TRIQS wrapper around Pomerol ED library"


// --- C++ Python converter for g2_iw_l_lp_params_t
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/string.hpp>
#include <algorithm>

namespace triqs { namespace py_tools {

template <> struct py_converter<g2_iw_l_lp_params_t> {
 static PyObject *c2py(g2_iw_l_lp_params_t const & x) {
  PyObject * d = PyDict_New();
  PyDict_SetItemString( d, "gf_struct"  , convert_to_python(x.gf_struct));
  PyDict_SetItemString( d, "beta"       , convert_to_python(x.beta));
  PyDict_SetItemString( d, "channel"    , convert_to_python(x.channel));
  PyDict_SetItemString( d, "block_order", convert_to_python(x.block_order));
  PyDict_SetItemString( d, "blocks"     , convert_to_python(x.blocks));
  PyDict_SetItemString( d, "n_iw"       , convert_to_python(x.n_iw));
  PyDict_SetItemString( d, "n_l"        , convert_to_python(x.n_l));
  PyDict_SetItemString( d, "n_inu_sum"  , convert_to_python(x.n_inu_sum));
  PyDict_SetItemString( d, "inu_sum_tol", convert_to_python(x.inu_sum_tol));
  return d;
 }

 template <typename T, typename U> static void _get_optional(PyObject *dic, const char *name, T &r, U const &init_default) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = init_default;
 }

 template <typename T> static void _get_optional(PyObject *dic, const char *name, T &r) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = T{};
 }

 static g2_iw_l_lp_params_t py2c(PyObject *dic) {
  g2_iw_l_lp_params_t res;
  res.gf_struct = convert_from_python<gf_struct_t>(PyDict_GetItemString(dic, "gf_struct"));
  res.beta = convert_from_python<double>(PyDict_GetItemString(dic, "beta"));
  _get_optional(dic, "channel"    , res.channel       ,PH);
  _get_optional(dic, "block_order", res.block_order   ,AABB);
  _get_optional(dic, "blocks"     , res.blocks        ,g2_blocks_t{});
  _get_optional(dic, "n_iw"       , res.n_iw          ,30);
  _get_optional(dic, "n_l"        , res.n_l           ,20);
  _get_optional(dic, "n_inu_sum"  , res.n_inu_sum     ,500);
  _get_optional(dic, "inu_sum_tol", res.inu_sum_tol   ,1e-6);
  return res;
 }

 template <typename T>
 static void _check(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!convertible_from_python<T>(PyDict_GetItemString(dic, name), false))
   fs << "\n" << ++err << " The parameter " << name << " does not have the right type : expecting " << tname
      << " in C++, but got '" << PyDict_GetItemString(dic, name)->ob_type->tp_name << "' in Python.";
 }

 template <typename T>
 static void _check_mandatory(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!PyDict_Contains(dic, pyref::string(name)))
   fs << "\n" << ++err << " Mandatory parameter " << name << " is missing.";
  else _check<T>(dic,fs,err,name,tname);
 }

 template <typename T>
 static void _check_optional(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (PyDict_Contains(dic, pyref::string(name))) _check<T>(dic, fs, err, name, tname);
 }

 static bool is_convertible(PyObject *dic, bool raise_exception) {
  if (dic == nullptr or !PyDict_Check(dic)) {
   if (raise_exception) { PyErr_SetString(PyExc_TypeError, "The function must be called with named arguments");}
   return false;
  }
  std::stringstream fs, fs2; int err=0;

#ifndef TRIQS_ALLOW_UNUSED_PARAMETERS
  std::vector<std::string> ks, all_keys = {"gf_struct","beta","channel","block_order","blocks","n_iw","n_l","n_inu_sum","inu_sum_tol"};
  pyref keys = PyDict_Keys(dic);
  if (!convertible_from_python<std::vector<std::string>>(keys, true)) {
   fs << "\nThe dict keys are not strings";
   goto _error;
  }
  ks = convert_from_python<std::vector<std::string>>(keys);
  for (auto & k : ks)
   if (std::find(all_keys.begin(), all_keys.end(), k) == all_keys.end())
    fs << "\n"<< ++err << " The parameter '" << k << "' is not recognized.";
#endif

  _check_mandatory<gf_struct_t                 >(dic, fs, err, "gf_struct"  , "gf_struct_t");
  _check_mandatory<double                      >(dic, fs, err, "beta"       , "double");
  _check_optional <pomerol2triqs::channel_t    >(dic, fs, err, "channel"    , "pomerol2triqs::channel_t");
  _check_optional <pomerol2triqs::block_order_t>(dic, fs, err, "block_order", "pomerol2triqs::block_order_t");
  _check_optional <g2_blocks_t                 >(dic, fs, err, "blocks"     , "g2_blocks_t");
  _check_optional <int                         >(dic, fs, err, "n_iw"       , "int");
  _check_optional <int                         >(dic, fs, err, "n_l"        , "int");
  _check_optional <int                         >(dic, fs, err, "n_inu_sum"  , "int");
  _check_optional <double                      >(dic, fs, err, "inu_sum_tol", "double");
  if (err) goto _error;
  return true;

 _error:
   fs2 << "\n---- There " << (err > 1 ? "are " : "is ") << err<< " error"<<(err >1 ?"s" : "")<< " in Python -> C++ transcription for the class g2_iw_l_lp_params_t\n" <<fs.str();
   if (raise_exception) PyErr_SetString(PyExc_TypeError, fs2.str().c_str());
  return false;
 }
};

}}


// --- C++ Python converter for g2_iw_inu_inup_params_t
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/string.hpp>
#include <algorithm>

namespace triqs { namespace py_tools {

template <> struct py_converter<g2_iw_inu_inup_params_t> {
 static PyObject *c2py(g2_iw_inu_inup_params_t const & x) {
  PyObject * d = PyDict_New();
  PyDict_SetItemString( d, "gf_struct"  , convert_to_python(x.gf_struct));
  PyDict_SetItemString( d, "beta"       , convert_to_python(x.beta));
  PyDict_SetItemString( d, "channel"    , convert_to_python(x.channel));
  PyDict_SetItemString( d, "block_order", convert_to_python(x.block_order));
  PyDict_SetItemString( d, "blocks"     , convert_to_python(x.blocks));
  PyDict_SetItemString( d, "n_iw"       , convert_to_python(x.n_iw));
  PyDict_SetItemString( d, "n_inu"      , convert_to_python(x.n_inu));
  return d;
 }

 template <typename T, typename U> static void _get_optional(PyObject *dic, const char *name, T &r, U const &init_default) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = init_default;
 }

 template <typename T> static void _get_optional(PyObject *dic, const char *name, T &r) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = T{};
 }

 static g2_iw_inu_inup_params_t py2c(PyObject *dic) {
  g2_iw_inu_inup_params_t res;
  res.gf_struct = convert_from_python<gf_struct_t>(PyDict_GetItemString(dic, "gf_struct"));
  res.beta = convert_from_python<double>(PyDict_GetItemString(dic, "beta"));
  _get_optional(dic, "channel"    , res.channel       ,PH);
  _get_optional(dic, "block_order", res.block_order   ,AABB);
  _get_optional(dic, "blocks"     , res.blocks        ,g2_blocks_t{});
  _get_optional(dic, "n_iw"       , res.n_iw          ,30);
  _get_optional(dic, "n_inu"      , res.n_inu         ,30);
  return res;
 }

 template <typename T>
 static void _check(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!convertible_from_python<T>(PyDict_GetItemString(dic, name), false))
   fs << "\n" << ++err << " The parameter " << name << " does not have the right type : expecting " << tname
      << " in C++, but got '" << PyDict_GetItemString(dic, name)->ob_type->tp_name << "' in Python.";
 }

 template <typename T>
 static void _check_mandatory(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!PyDict_Contains(dic, pyref::string(name)))
   fs << "\n" << ++err << " Mandatory parameter " << name << " is missing.";
  else _check<T>(dic,fs,err,name,tname);
 }

 template <typename T>
 static void _check_optional(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (PyDict_Contains(dic, pyref::string(name))) _check<T>(dic, fs, err, name, tname);
 }

 static bool is_convertible(PyObject *dic, bool raise_exception) {
  if (dic == nullptr or !PyDict_Check(dic)) {
   if (raise_exception) { PyErr_SetString(PyExc_TypeError, "The function must be called with named arguments");}
   return false;
  }
  std::stringstream fs, fs2; int err=0;

#ifndef TRIQS_ALLOW_UNUSED_PARAMETERS
  std::vector<std::string> ks, all_keys = {"gf_struct","beta","channel","block_order","blocks","n_iw","n_inu"};
  pyref keys = PyDict_Keys(dic);
  if (!convertible_from_python<std::vector<std::string>>(keys, true)) {
   fs << "\nThe dict keys are not strings";
   goto _error;
  }
  ks = convert_from_python<std::vector<std::string>>(keys);
  for (auto & k : ks)
   if (std::find(all_keys.begin(), all_keys.end(), k) == all_keys.end())
    fs << "\n"<< ++err << " The parameter '" << k << "' is not recognized.";
#endif

  _check_mandatory<gf_struct_t                 >(dic, fs, err, "gf_struct"  , "gf_struct_t");
  _check_mandatory<double                      >(dic, fs, err, "beta"       , "double");
  _check_optional <pomerol2triqs::channel_t    >(dic, fs, err, "channel"    , "pomerol2triqs::channel_t");
  _check_optional <pomerol2triqs::block_order_t>(dic, fs, err, "block_order", "pomerol2triqs::block_order_t");
  _check_optional <g2_blocks_t                 >(dic, fs, err, "blocks"     , "g2_blocks_t");
  _check_optional <int                         >(dic, fs, err, "n_iw"       , "int");
  _check_optional <int                         >(dic, fs, err, "n_inu"      , "int");
  if (err) goto _error;
  return true;

 _error:
   fs2 << "\n---- There " << (err > 1 ? "are " : "is ") << err<< " error"<<(err >1 ?"s" : "")<< " in Python -> C++ transcription for the class g2_iw_inu_inup_params_t\n" <<fs.str();
   if (raise_exception) PyErr_SetString(PyExc_TypeError, fs2.str().c_str());
  return false;
 }
};

}}