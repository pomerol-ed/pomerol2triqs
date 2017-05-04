# Generated automatically using the command :
# c++2py.py ../c++/main.hpp -I../../../pomerol/installed/include -I/usr/include/eigen3 -I../c++ -p -mpytriqs.applications.impurity_solvers.pomerol2triqs -o pomerol2triqs --moduledoc "TRIQS wrapper around Pomerol ED library"
from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.applications.impurity_solvers.pomerol2triqs", doc = "TRIQS wrapper around Pomerol ED library")

# All the triqs C++/Python modules

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/main.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/pair.hpp>
#include <triqs/python_tools/converters/map.hpp>
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/variant_int_string.hpp>
#include <triqs/python_tools/converters/gf.hpp>
#include <triqs/python_tools/converters/block_gf.hpp>
""")

# The class pomerol2triqs
c = class_(
        py_type = "Pomerol2TRIQS",  # name of the python class
        c_type = "pomerol2triqs",   # name of the C++ class
        doc = r"",   # doc of the C++ class
)

c.add_method("""block_gf<triqs::gfs::imfreq> G_iw (gf_struct_t gf_struct, int n_iw)""",
             doc = """Green\'s function in Matsubara frequencies """)

c.add_method("""block_gf<triqs::gfs::imtime> G_tau (gf_struct_t gf_struct, int n_tau)""",
             doc = """Green\'s function in imaginary time """)

module.add_class(c)

module.generate_code()
