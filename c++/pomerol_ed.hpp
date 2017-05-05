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
using triqs::hilbert_space::fundamental_operator_set;

using indices_t = fundamental_operator_set::indices_t;
using pomerol_indices_t = std::tuple<std::string, unsigned int, Pomerol::spin>;
using index_converter_t = std::map<indices_t, pomerol_indices_t>;

class pomerol_ed {

 boost::mpi::communicator comm;

 const bool verbose;
 index_converter_t index_converter;
 Pomerol::Lattice bare_lattice;
 Pomerol::IndexClassification index_info;

 many_body_op_t h;

 Pomerol::Lattice init();

 std::unique_ptr<Pomerol::Lattice> lattice;
 std::unique_ptr<Pomerol::IndexHamiltonian> storage;
 std::unique_ptr<Pomerol::Symmetrizer> symm;
 std::unique_ptr<Pomerol::StatesClassification> states_class;
 std::unique_ptr<Pomerol::Hamiltonian> matrix_h;
 std::unique_ptr<Pomerol::DensityMatrix> rho;

public:

 pomerol_ed(index_converter_t const& index_converter, bool verbose = false);

 void diagonalize(many_body_op_t const& hamiltonian, bool ignore_symmetries = false);

 many_body_op_t const& hamiltonian() const { return h; }

 /// Green's function in Matsubara frequencies
 block_gf<imfreq> G_iw(gf_struct_t const& gf_struct, int n_iw) const;

 /// Green's function in imaginary time
 block_gf<imtime> G_tau(gf_struct_t const& gf_struct, int n_tau) const;

};

}
