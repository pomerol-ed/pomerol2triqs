#pragma once

#include <pomerol.h>
#include <triqs/gfs.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>

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

class pomerol2triqs {

 many_body_op_t h;
 fundamental_operator_set fops;

public:

 pomerol2triqs(many_body_op_t const& hamiltonian);

 many_body_op_t const& hamiltonian() const { return h; }

 /// Green's function in Matsubara frequencies
 block_gf<imfreq> G_iw(gf_struct_t const& gf_struct, int n_iw) const;

 /// Green's function in imaginary time
 block_gf<imtime> G_tau(gf_struct_t const& gf_struct, int n_tau) const;

};
