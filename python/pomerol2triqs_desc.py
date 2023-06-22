# Generated automatically using the command :
# c++2py ./../c++/pomerol_ed.hpp -p -m pomerol2triqs -o pomerol2triqs --moduledoc="TRIQS wrapper around Pomerol ED library" --cxxflags="-std=c++20 -stdlib=libc++" -C triqs --only="pomerol_ed block_order_t channel_t spin" -N Pomerol::LatticePresets -N pomerol2triqs -I../c++ -I/usr/include/eigen3 -I${POMEROL_DIR}/include
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "pomerol2triqs", doc = r"TRIQS wrapper around Pomerol ED library", app_name = "pomerol2triqs")

# Imports
module.add_imports(*['triqs.gf', 'triqs.gf.meshes', 'triqs.operators'])

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
#include <cpp2py/converters/variant.hpp>
#include <cpp2py/converters/vector.hpp>
#include <triqs/cpp2py_converters/gf.hpp>
#include <triqs/cpp2py_converters/mesh.hpp>
#include <triqs/cpp2py_converters/operators_real_complex.hpp>
#include <triqs/cpp2py_converters/real_or_complex.hpp>

using namespace Pomerol::LatticePresets;
using namespace pomerol2triqs;
""")

module.add_enum("spin", ['undef', 'down', 'up'], "Pomerol::LatticePresets", doc = r"""Possible values of spin-1/2 z-projection.""")
module.add_enum("block_order_t", ['AABB', 'ABBA'], "pomerol2triqs", doc = r"""Order of block indices for Block2Gf objects""")
module.add_enum("channel_t", ['PP', 'PH', 'xPH', 'AllFermionic'], "pomerol2triqs", doc = r"""Channel in which Matsubara frequency representation is defined""")

# The class pomerol_ed
c = class_(
        py_type = "PomerolED",  # name of the python class
        c_type = "pomerol2triqs::pomerol_ed",   # name of the C++ class
        doc = r"""Main solver class of pomerol2triqs""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(pomerol2triqs::index_converter_t index_converter, bool verbose = false)""", doc = r"""Create a new solver object""")

c.add_method("""void diagonalize (pomerol2triqs::many_body_op_t hamiltonian, bool ignore_symmetries = false)""",
             doc = r"""Diagonalize Hamiltonian optionally employing its symmetries""")

c.add_method("""std::complex<double> ensemble_average (pomerol2triqs::indices_t i, pomerol2triqs::indices_t j, double beta)""",
             doc = r"""Compute the ensemble average of c^+_i c_j""")

c.add_method("""block_gf<mesh::imfreq> G_iw (triqs::gfs::gf_struct_t gf_struct, double beta, int n_iw)""",
             doc = r"""Green's function in Matsubara frequencies""")

c.add_method("""block_gf<mesh::imtime> G_tau (triqs::gfs::gf_struct_t gf_struct, double beta, int n_tau)""",
             doc = r"""Green's function in imaginary time""")

c.add_method("""block_gf<mesh::refreq> G_w (triqs::gfs::gf_struct_t gf_struct, double beta, std::pair<double, double> energy_window, int n_w, double im_shift = 0)""",
             doc = r"""Retarded Green's function on real energy axis""")

c.add_method("""block2_gf<pomerol2triqs::w_nu_nup_t, tensor_valued<4>> G2_iw_inu_inup (**pomerol2triqs::g2_iw_inu_inup_params_t)""",
             doc = r"""Two-particle Green's function, Matsubara frequencies



+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| Parameter Name | Type                         | Default            | Documentation                                                    |
+================+==============================+====================+==================================================================+
| gf_struct      | triqs::gfs::gf_struct_t      | --                 | Structure of G^2 blocks.                                         |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| beta           | double                       | --                 | Inverse temperature                                              |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| channel        | pomerol2triqs::channel_t     | PH                 | Channel in which Matsubara frequency representation is defined.  |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| block_order    | pomerol2triqs::block_order_t | AABB               | Order of block indices in the definition of G^2.                 |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| blocks         | pomerol2triqs::g2_blocks_t   | measure all blocks | List of block index pairs of G^2 to measure.                     |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| n_iw           | int                          | 30                 | Number of bosonic Matsubara frequencies.                         |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| n_inu          | int                          | 30                 | Number of fermionic Matsubara frequencies.                       |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
""")

c.add_method("""block2_gf<pomerol2triqs::w_l_lp_t, tensor_valued<4>> G2_iw_l_lp (**pomerol2triqs::g2_iw_l_lp_params_t)""",
             doc = r"""Two-particle Green's function, bosonic Matsubara frequency + Legendre coefficients



+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| Parameter Name | Type                         | Default            | Documentation                                                    |
+================+==============================+====================+==================================================================+
| gf_struct      | triqs::gfs::gf_struct_t      | --                 | Structure of G^2 blocks.                                         |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| beta           | double                       | --                 | Inverse temperature                                              |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| channel        | pomerol2triqs::channel_t     | PH                 | Channel in which Matsubara frequency representation is defined.  |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| block_order    | pomerol2triqs::block_order_t | AABB               | Order of block indices in the definition of G^2.                 |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| blocks         | pomerol2triqs::g2_blocks_t   | measure all blocks | List of block index pairs of G^2 to measure.                     |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| n_iw           | int                          | 30                 | Number of bosonic Matsubara frequencies.                         |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| n_l            | int                          | 20                 | Number of Legendre coefficients.                                 |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| n_inu_sum      | int                          | 500                | Maximum number of positive Matsubara frequencies in summation.   |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| inu_sum_tol    | double                       | 1e-6               | Tolerance for Matsubara frequency summation.                     |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
""")

c.add_method("""gf<mesh::imtime, triqs::gfs::scalar_valued> chi_tau (pomerol2triqs::indices_t i, pomerol2triqs::indices_t j, pomerol2triqs::indices_t k, pomerol2triqs::indices_t l, double beta, int n_tau, bool connected = false, channel_t channel = PH)""",
             doc = r"""Dynamical susceptibility <T c^+_{i_1}(\tau) c_{j_1}(\tau) c^+_{i_2}(0) c_{j_2}(0)> or its connected part""")

c.add_method("""gf<mesh::imfreq, triqs::gfs::scalar_valued> chi_iw (pomerol2triqs::indices_t i, pomerol2triqs::indices_t j, pomerol2triqs::indices_t k, pomerol2triqs::indices_t l, double beta, int n_iw, bool connected = false, channel_t channel = PH)""",
             doc = r"""Dynamical susceptibility <T c^+_{i_1}(\tau) c_{j_1}(\tau) c^+_{i_2}(0) c_{j_2}(0)> or its connected part in Matsubara frequencies""")

c.add_method("""block2_gf<pomerol2triqs::w_nu_t, tensor_valued<4>> chi3_iw_inu (**pomerol2triqs::chi3_iw_inu_params_t)""",
             doc = r"""3-point fermion-boson susceptibility



+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| Parameter Name | Type                         | Default            | Documentation                                                    |
+================+==============================+====================+==================================================================+
| gf_struct      | triqs::gfs::gf_struct_t      | --                 | Structure of \chi^3 blocks.                                      |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| beta           | double                       | --                 | Inverse temperature                                              |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| channel        | pomerol2triqs::channel_t     | PH                 | Channel in which Matsubara frequency representation is defined.  |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| block_order    | pomerol2triqs::block_order_t | AABB               | Order of block indices in the definition of \chi^3.              |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| blocks         | pomerol2triqs::chi3_blocks_t | measure all blocks | List of block index pairs of \chi^3 to measure.                  |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| n_iw           | int                          | 100                | Number of bosonic Matsubara frequencies.                         |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
| n_inu          | int                          | 100                | Number of fermionic Matsubara frequencies.                       |
+----------------+------------------------------+--------------------+------------------------------------------------------------------+
""")

c.add_property(name = "rho_threshold",
               getter = cfunction("double get_rho_threshold ()"),
               setter = cfunction("void set_rho_threshold (double threshold)"),
               doc = r"""Truncation threshold for density matrix elements""")

module.add_class(c)


# Converter for g2_iw_inu_inup_params_t
c = converter_(
        c_type = "pomerol2triqs::g2_iw_inu_inup_params_t",
        doc = r"""Arguments of G2_iw_inu_inup()""",
)
c.add_member(c_name = "gf_struct",
             c_type = "triqs::gfs::gf_struct_t",
             initializer = """  """,
             doc = r"""Structure of G^2 blocks.""")

c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = r"""Inverse temperature""")

c.add_member(c_name = "channel",
             c_type = "pomerol2triqs::channel_t",
             initializer = """ PH """,
             doc = r"""Channel in which Matsubara frequency representation is defined.""")

c.add_member(c_name = "block_order",
             c_type = "pomerol2triqs::block_order_t",
             initializer = """ AABB """,
             doc = r"""Order of block indices in the definition of G^2.""")

c.add_member(c_name = "blocks",
             c_type = "pomerol2triqs::g2_blocks_t",
             initializer = """ g2_blocks_t{} """,
             doc = r"""List of block index pairs of G^2 to measure.
     default: measure all blocks""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 30 """,
             doc = r"""Number of bosonic Matsubara frequencies.""")

c.add_member(c_name = "n_inu",
             c_type = "int",
             initializer = """ 30 """,
             doc = r"""Number of fermionic Matsubara frequencies.""")

module.add_converter(c)

# Converter for g2_iw_l_lp_params_t
c = converter_(
        c_type = "pomerol2triqs::g2_iw_l_lp_params_t",
        doc = r"""Arguments of G2_iw_l_lp()""",
)
c.add_member(c_name = "gf_struct",
             c_type = "triqs::gfs::gf_struct_t",
             initializer = """  """,
             doc = r"""Structure of G^2 blocks.""")

c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = r"""Inverse temperature""")

c.add_member(c_name = "channel",
             c_type = "pomerol2triqs::channel_t",
             initializer = """ PH """,
             doc = r"""Channel in which Matsubara frequency representation is defined.""")

c.add_member(c_name = "block_order",
             c_type = "pomerol2triqs::block_order_t",
             initializer = """ AABB """,
             doc = r"""Order of block indices in the definition of G^2.""")

c.add_member(c_name = "blocks",
             c_type = "pomerol2triqs::g2_blocks_t",
             initializer = """ g2_blocks_t{} """,
             doc = r"""List of block index pairs of G^2 to measure.
     default: measure all blocks""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 30 """,
             doc = r"""Number of bosonic Matsubara frequencies.""")

c.add_member(c_name = "n_l",
             c_type = "int",
             initializer = """ 20 """,
             doc = r"""Number of Legendre coefficients.""")

c.add_member(c_name = "n_inu_sum",
             c_type = "int",
             initializer = """ 500 """,
             doc = r"""Maximum number of positive Matsubara frequencies in summation.""")

c.add_member(c_name = "inu_sum_tol",
             c_type = "double",
             initializer = """ 1e-6 """,
             doc = r"""Tolerance for Matsubara frequency summation.""")

module.add_converter(c)

# Converter for chi3_iw_inu_params_t
c = converter_(
        c_type = "pomerol2triqs::chi3_iw_inu_params_t",
        doc = r"""Arguments of chi3_iw_inu()""",
)
c.add_member(c_name = "gf_struct",
             c_type = "triqs::gfs::gf_struct_t",
             initializer = """  """,
             doc = r"""Structure of \chi^3 blocks.""")

c.add_member(c_name = "beta",
             c_type = "double",
             initializer = """  """,
             doc = r"""Inverse temperature""")

c.add_member(c_name = "channel",
             c_type = "pomerol2triqs::channel_t",
             initializer = """ PH """,
             doc = r"""Channel in which Matsubara frequency representation is defined.""")

c.add_member(c_name = "block_order",
             c_type = "pomerol2triqs::block_order_t",
             initializer = """ AABB """,
             doc = r"""Order of block indices in the definition of \chi^3.""")

c.add_member(c_name = "blocks",
             c_type = "pomerol2triqs::chi3_blocks_t",
             initializer = """ chi3_blocks_t{} """,
             doc = r"""List of block index pairs of \chi^3 to measure.
     default: measure all blocks""")

c.add_member(c_name = "n_iw",
             c_type = "int",
             initializer = """ 100 """,
             doc = r"""Number of bosonic Matsubara frequencies.""")

c.add_member(c_name = "n_inu",
             c_type = "int",
             initializer = """ 100 """,
             doc = r"""Number of fermionic Matsubara frequencies.""")

module.add_converter(c)


module.generate_code()
