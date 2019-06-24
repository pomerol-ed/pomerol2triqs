# Generated automatically using the command :
# c++2py ../c++/pomerol_ed.hpp -p -m pomerol2triqs -o pomerol2triqs --moduledoc="TRIQS wrapper around Pomerol ED library" --cxxflags="-std=c++14" -C pytriqs --only="pomerol_ed block_order_t channel_t" -N pomerol2triqs -I../c++ -I/usr/include/eigen3 -I${POMEROL_DIR}/include
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "pomerol2triqs", doc = "TRIQS wrapper around Pomerol ED library", app_name = "pomerol2triqs")

# Imports
import pytriqs.gf
import pytriqs.operators

# Add here all includes
module.add_include("pomerol_ed.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <cpp2py/converters/map.hpp>
#include <cpp2py/converters/pair.hpp>
#include <cpp2py/converters/set.hpp>
#include <cpp2py/converters/string.hpp>
#include <cpp2py/converters/tuple.hpp>
#include <cpp2py/converters/vector.hpp>
#include <cpp2py/converters/variant.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>

using namespace pomerol2triqs;
using namespace Pomerol;
""")

module.add_enum("spin",          ['down', 'up'], "Pomerol", "Spin projection")
module.add_enum("block_order_t", ['AABB', 'ABBA'], "pomerol2triqs", """Order of block indices for Block2Gf objects""")
module.add_enum("channel_t", ['PP', 'PH', 'AllFermionic'], "pomerol2triqs", """Channel in which Matsubara frequency representation is defined""")

# The class pomerol_ed
c = class_(
        py_type = "PomerolED",  # name of the python class
        c_type = "pomerol2triqs::pomerol_ed",   # name of the C++ class
        doc = """Main solver class of pomerol2triqs""",   # doc of the C++ class
)

c.add_constructor("""(index_converter_t index_converter, bool verbose = false)""", doc = """Create a new solver object""")

c.add_property(name = "rho_threshold",
               getter = cfunction("double get_rho_threshold()"),
               setter = cfunction("void set_rho_threshold(double threshold)"),
               doc = """Truncation threshold for density matrix elements""")

c.add_method("""void diagonalize (many_body_op_t hamiltonian, bool ignore_symmetries = false)""",
             doc = """Diagonalize Hamiltonian optionally employing conservation of N and S_z""")

c.add_method("""void diagonalize (many_body_op_t hamiltonian, std::vector<many_body_op_t> integrals_of_motion)""",
             doc = """Diagonalize Hamiltonian using provided integrals of motion""")

c.add_method("""std::complex<double> ensemble_average(indices_t i, indices_t j, double beta)""",
             doc = """Compute the ensemble average of c^+_i c_j""")

c.add_method("""block_gf<triqs::gfs::imfreq> G_iw (gf_struct_t gf_struct, double beta, int n_iw)""",
             doc = """Green\'s function in Matsubara frequencies""")

c.add_method("""block_gf<triqs::gfs::imtime> G_tau (gf_struct_t gf_struct, double beta, int n_tau)""",
             doc = """Green\'s function in imaginary time""")

c.add_method("""block_gf<triqs::gfs::refreq> G_w (gf_struct_t gf_struct, double beta, std::pair<double,double> energy_window, int n_w, double im_shift = 0)""",
             doc = """Retarded Green\'s function on real energy axis""")

c.add_method("""block2_gf<w_nu_nup_t,tensor_valued<4>> G2_iw_inu_inup (**pomerol2triqs::g2_iw_inu_inup_params_t)""",
             doc = """Two-particle Green\'s function, Matsubara frequencies
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| Parameter Name | Type                         | Default            | Documentation                                                              |
+================+==============================+====================+============================================================================+
| gf_struct      | gf_struct_t                  |                    | Structure of G^2 blocks.                                                   |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| beta           | double                       |                    | Inverse temperature                                                        |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| channel        | pomerol2triqs::channel_t     | PH                 | Channel in which Matsubara frequency representation is defined.            |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| block_order    | pomerol2triqs::block_order_t | AABB               | Order of block indices in the definition of G^2.                           |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| blocks         | g2_blocks_t                  | measure all blocks | List of block index pairs of G^2 to measure.                               |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| n_iw           | int                          | 30                 | Number of bosonic Matsubara frequencies.                                   |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| n_inu          | int                          | 30                 | Number of fermionic Matsubara frequencies                                  |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+""")

c.add_method("""block2_gf<w_l_lp_t,tensor_valued<4>> G2_iw_l_lp (**pomerol2triqs::g2_iw_l_lp_params_t)""",
             doc = """Two-particle Green\'s function, bosonic Matsubara frequency + Legendre coefficients
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| Parameter Name | Type                         | Default            | Documentation                                                              |
+================+==============================+====================+============================================================================+
| gf_struct      | gf_struct_t                  |                    | Structure of G^2 blocks.                                                   |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| beta           | double                       |                    | Inverse temperature                                                        |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| channel        | pomerol2triqs::channel_t     | PH                 | Channel in which Matsubara frequency representation is defined.            |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| block_order    | pomerol2triqs::block_order_t | AABB               | Order of block indices in the definition of G^2.                           |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| blocks         | g2_blocks_t                  | measure all blocks | List of block index pairs of G^2 to measure.                               |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| n_iw           | int                          | 30                 | Number of bosonic Matsubara frequencies.                                   |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| n_l            | int                          | 20                 | Number of Legendre coefficients.                                           |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| n_inu_sum      | int                          | 500                | Maximum number of positive Matsubara frequencies in summation.             |
+----------------+------------------------------+--------------------+----------------------------------------------------------------------------+
| inu_sum_tol    | double                       | 1e-6               | Tolerance for Matsubara frequency summation.                               |
+----------------+------------------------------+--------------------+--------------------------------------------------------------------------0-+""")

module.add_class(c)


# Converter for g2_iw_inu_inup_params_t
c = converter_(
        c_type = "pomerol2triqs::g2_iw_inu_inup_params_t",
        doc = """Arguments of G2_iw_inu_inup()""",
)
c.add_member(c_name = "gf_struct",
             c_type = "gf_struct_t",
             initializer = """  """,
             doc = """Structure of G^2 blocks.""")

c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = """Inverse temperature""")

c.add_member(c_name = "channel",
             c_type = "pomerol2triqs::channel_t",
             initializer = """ PH """,
             doc = """Channel in which Matsubara frequency representation is defined.""")

c.add_member(c_name = "block_order",
             c_type = "pomerol2triqs::block_order_t",
             initializer = """ AABB """,
             doc = """Order of block indices in the definition of G^2.""")

c.add_member(c_name = "blocks",
             c_type = "g2_blocks_t",
             initializer = """ g2_blocks_t{} """,
             doc = """List of block index pairs of G^2 to measure.\n     default: measure all blocks""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 30 """,
             doc = """Number of bosonic Matsubara frequencies.""")

c.add_member(c_name = "n_inu",
             c_type = "int",
             initializer = """ 30 """,
             doc = """Number of fermionic Matsubara frequencies.""")

module.add_converter(c)

# Converter for g2_iw_l_lp_params_t
c = converter_(
        c_type = "pomerol2triqs::g2_iw_l_lp_params_t",
        doc = """Arguments of G2_iw_l_lp()""",
)
c.add_member(c_name = "gf_struct",
             c_type = "gf_struct_t",
             initializer = """  """,
             doc = """Structure of G^2 blocks.""")

c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = """Inverse temperature""")

c.add_member(c_name = "channel",
             c_type = "pomerol2triqs::channel_t",
             initializer = """ PH """,
             doc = """Channel in which Matsubara frequency representation is defined.""")

c.add_member(c_name = "block_order",
             c_type = "pomerol2triqs::block_order_t",
             initializer = """ AABB """,
             doc = """Order of block indices in the definition of G^2.""")

c.add_member(c_name = "blocks",
             c_type = "g2_blocks_t",
             initializer = """ g2_blocks_t{} """,
             doc = """List of block index pairs of G^2 to measure.\n     default: measure all blocks""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 30 """,
             doc = """Number of bosonic Matsubara frequencies.""")

c.add_member(c_name = "n_l",
             c_type = "int",
             initializer = """ 20 """,
             doc = """Number of Legendre coefficients.""")

c.add_member(c_name = "n_inu_sum",
             c_type = "int",
             initializer = """ 500 """,
             doc = """Maximum number of positive Matsubara frequencies in summation.""")

c.add_member(c_name = "inu_sum_tol",
             c_type = "double",
             initializer = """ 1e-6 """,
             doc = """Tolerance for Matsubara frequency summation.""")

module.add_converter(c)


module.generate_code()
