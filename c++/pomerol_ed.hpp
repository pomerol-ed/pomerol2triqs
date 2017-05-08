#pragma once

#include <utility>
#include <vector>
#include <memory>
#include <tuple>
#include <functional>

#include <pomerol.h>
#include <triqs/gfs.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/utility/optional_compat.hpp>

#include "g2_parameters.hpp"

namespace pomerol2triqs {

#ifdef POMEROL_COMPLEX_MATRIX_ELEMENTS
using h_scalar_t = std::complex<double>;
#else
using h_scalar_t = double;
#endif

// Operator with real or complex value
using many_body_op_t = triqs::operators::many_body_operator_generic<h_scalar_t>;

using namespace triqs::gfs;
using triqs::hilbert_space::gf_struct_t;

using indices_t = triqs::hilbert_space::fundamental_operator_set::indices_t;
using pomerol_indices_t = std::tuple<std::string, unsigned int, Pomerol::spin>;
using index_converter_t = std::map<indices_t, pomerol_indices_t>;

class pomerol_ed {

 boost::mpi::communicator comm;

 const bool verbose;
 index_converter_t index_converter;
 Pomerol::Lattice bare_lattice;
 Pomerol::IndexClassification index_info;

 many_body_op_t h;

 std::unique_ptr<Pomerol::Lattice> lattice;
 std::unique_ptr<Pomerol::IndexHamiltonian> storage;
 std::unique_ptr<Pomerol::Symmetrizer> symm;
 std::unique_ptr<Pomerol::StatesClassification> states_class;
 std::unique_ptr<Pomerol::Hamiltonian> matrix_h;
 std::unique_ptr<Pomerol::DensityMatrix> rho;
 std::set<Pomerol::ParticleIndex> computed_ops;
 std::unique_ptr<Pomerol::FieldOperatorContainer> ops_container;
 std::unique_ptr<Pomerol::GFContainer> gf_container;
 g2_blocks_t computed_g2_blocks;
 std::unique_ptr<Pomerol::TwoParticleGFContainer> g2_container;

 Pomerol::Lattice init();
 Pomerol::ParticleIndex lookup_pomerol_index(indices_t const& i) const;
 std::set<Pomerol::ParticleIndex> gf_struct_to_pomerol_indices(gf_struct_t const& gf_struct) const;
 void compute_rho(double beta);
 void compute_field_operators(gf_struct_t const& gf_struct);
 void compute_gfs();
 void compute_g2(gf_struct_t const& gf_struct, g2_blocks_t g2_blocks);
 template<typename Mesh, typename Filler>
 block_gf<Mesh> fill_gf(gf_struct_t const& gf_struct, gf_mesh<Mesh> const& mesh, Filler filler) const;

public:

 pomerol_ed(index_converter_t const& index_converter, bool verbose = false);

 void diagonalize(many_body_op_t const& hamiltonian, bool ignore_symmetries = false);

 many_body_op_t const& hamiltonian() const { return h; }

 /// Green's function in Matsubara frequencies
 block_gf<imfreq> G_iw(gf_struct_t const& gf_struct, double beta, int n_iw);

 /// Green's function in imaginary time
 block_gf<imtime> G_tau(gf_struct_t const& gf_struct, double beta, int n_tau);

 /// Retarded Green's function on real energy axis
 block_gf<refreq> G_w(gf_struct_t const& gf_struct,
                      double beta,
                      std::pair<double, double> const& energy_window,
                      int n_w,
                      double im_shift = 0
                     );

 /// Two-particle Green's function, Matsubara frequencies
 TRIQS_WRAP_ARG_AS_DICT
 block2_gf<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> G2_iw(g2_parameters_t const& p);

};

}
