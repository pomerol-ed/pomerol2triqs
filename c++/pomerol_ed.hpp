/**
 * pomerol2triqs
 *
 * Copyright (C) 2017-2022 Igor Krivenko <igor.s.krivenko @ gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include <triqs/utility/first_include.hpp>

#include <utility>
#include <vector>
#include <memory>
#include <tuple>
#include <functional>
#include <type_traits>
#include <variant>

#include <pomerol.hpp>
#include <triqs/gfs.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>

#include "parameters.hpp"

namespace pomerol2triqs {

  // Operator with real or complex value (switchable at runtime)
  using many_body_op_t = triqs::operators::many_body_operator;

  using namespace triqs::gfs;
  namespace mesh = triqs::mesh;
  using triqs::hilbert_space::gf_struct_t;

  using indices_t              = triqs::hilbert_space::fundamental_operator_set::indices_t;
  using pomerol_indices_t      = std::tuple<std::string, unsigned short, Pomerol::LatticePresets::spin>;
  using index_classification_t = Pomerol::IndexClassification<std::string, unsigned short, Pomerol::LatticePresets::spin>;
  using hilbert_space_t        = Pomerol::HilbertSpace<std::string, unsigned short, Pomerol::LatticePresets::spin>;

  static_assert(std::is_same<pomerol_indices_t, Pomerol::LatticePresets::RealExpr::index_types>::value);
  static_assert(std::is_same<pomerol_indices_t, Pomerol::LatticePresets::ComplexExpr::index_types>::value);

  // We use 'unsigned int' instead of 'unsigned short' in this type because
  // cpp2py does not know how to convert the latter integer type.
  using index_converter_t = std::map<indices_t, std::tuple<std::string, unsigned int, Pomerol::LatticePresets::spin>>;

  using nu_nup_t = mesh::prod<mesh::imfreq, mesh::imfreq>;
  using w_nu_nup_t = mesh::prod<mesh::imfreq, mesh::imfreq, mesh::imfreq>;
  using w_l_lp_t   = mesh::prod<mesh::imfreq, mesh::legendre, mesh::legendre>;

  /// Main solver class of pomerol2triqs
  class pomerol_ed {

    mpi::communicator comm;

    const bool verbose;
    index_converter_t index_converter;
    index_classification_t index_info;

    double rho_threshold = 0;

    using h_expr_t = std::variant<Pomerol::LatticePresets::RealExpr, Pomerol::LatticePresets::ComplexExpr>;
    std::unique_ptr<h_expr_t> h_expr;
    std::unique_ptr<hilbert_space_t> hs;
    std::unique_ptr<Pomerol::StatesClassification> states_class;
    std::unique_ptr<Pomerol::Hamiltonian> matrix_h;
    std::unique_ptr<Pomerol::DensityMatrix> rho;
    std::set<Pomerol::ParticleIndex> computed_ops;
    std::unique_ptr<Pomerol::FieldOperatorContainer> ops_container;

    Pomerol::ParticleIndex lookup_pomerol_index(indices_t const &i) const;
    std::set<Pomerol::ParticleIndex> gf_struct_to_pomerol_indices(gf_struct_t const &gf_struct) const;
    template <typename HExprType> void diagonalize_prepare_impl(many_body_op_t const &hamiltonian);
    void diagonalize_prepare(many_body_op_t const &hamiltonian);

    void compute_rho(double beta);
    void compute_field_operators(gf_struct_t const &gf_struct);
    template <typename Mesh, typename Filler> block_gf<Mesh> compute_gf(gf_struct_t const &gf_struct, Mesh const &mesh, Filler filler) const;
    template <typename Mesh, typename Filler>
    gf<Mesh, scalar_valued> compute_chi(indices_t const &i1, indices_t const &j1, indices_t const &i2, indices_t const &j2, bool connected,
                                        Mesh const &mesh, Filler filler) const;

    template <typename Mesh, typename Filler>
    block2_gf<Mesh, tensor_valued<4>> compute_g2(gf_struct_t const &gf_struct, Mesh const &mesh, block_order_t block_order,
                                                 g2_blocks_t const &g2_blocks, Filler filler) const;

    template <typename Mesh, typename Filler>
    block2_gf<Mesh, tensor_valued<4>> compute_chi3(gf_struct_t const &gf_struct, Mesh const &mesh, block_order_t block_order,
                                                   channel_t channel, chi3_blocks_t const &chi3_blocks, Filler filler) const;

    public:
    /// Create a new solver object
    pomerol_ed(index_converter_t const &index_converter, bool verbose = false);

    /// Diagonalize Hamiltonian optionally employing its symmetries
    void diagonalize(many_body_op_t const &hamiltonian, bool ignore_symmetries = false);

    /// Compute the ensemble average of c^+_i c_j
    std::complex<double> ensemble_average(indices_t const &i, indices_t const &j, double beta);

    /// Green's function in Matsubara frequencies
    block_gf<mesh::imfreq> G_iw(gf_struct_t const &gf_struct, double beta, int n_iw);

    /// Green's function in imaginary time
    block_gf<mesh::imtime> G_tau(gf_struct_t const &gf_struct, double beta, int n_tau);

    /// Retarded Green's function on real energy axis
    block_gf<mesh::refreq> G_w(gf_struct_t const &gf_struct, double beta, std::pair<double, double> const &energy_window, int n_w,
                               double im_shift = 0);

    /// Two-particle Green's function, Matsubara frequencies
    CPP2PY_ARG_AS_DICT
    block2_gf<w_nu_nup_t, tensor_valued<4>> G2_iw_inu_inup(g2_iw_inu_inup_params_t const &p);

    /// Two-particle Green's function, bosonic Matsubara frequency + Legendre coefficients
    CPP2PY_ARG_AS_DICT
    block2_gf<w_l_lp_t, tensor_valued<4>> G2_iw_l_lp(g2_iw_l_lp_params_t const &p);

    /// Dynamical susceptibility <T c^+_{i_1}(\tau) c_{j_1}(\tau) c^+_{i_2}(0) c_{j_2}(0)> or its connected part
    gf<mesh::imtime, scalar_valued> chi_tau(indices_t const &i1, indices_t const &j1, indices_t const &i2, indices_t const &j2, double beta,
                                            int n_tau, bool connected = false);

    /// Dynamical susceptibility <T c^+_{i_1}(\tau) c_{j_1}(\tau) c^+_{i_2}(0) c_{j_2}(0)> or its connected part in Matsubara frequencies
    gf<mesh::imfreq, scalar_valued> chi_inu(indices_t const &i1, indices_t const &j1, indices_t const &i2, indices_t const &j2, double beta,
                                            int n_inu, bool connected = false);

    /// 3-point fermion-boson susceptibility
    CPP2PY_ARG_AS_DICT
    block2_gf<nu_nup_t, tensor_valued<4>> chi3_inu_inup(chi3_inu_inup_params_t const& p);

    /// Get truncation threshold for density matrix elements
    double get_rho_threshold() const { return rho_threshold; }

    /// Set truncation threshold for density matrix elements
    void set_rho_threshold(double threshold) { rho_threshold = threshold; }
  };
} // namespace pomerol2triqs
