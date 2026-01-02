/**
 * pomerol2triqs
 *
 * Copyright (C) 2017-2026 Igor Krivenko
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

#include <triqs/hilbert_space/fundamental_operator_set.hpp>
using triqs::hilbert_space::gf_struct_t;

namespace pomerol2triqs {

  /// Order of block indices for Block2Gf objects
  enum block_order_t { AABB, ABBA };
  /// Channel in which Matsubara frequency representation is defined
  enum channel_t { PP, PH, xPH, AllFermionic };

  using g2_blocks_t   = std::set<std::pair<std::string, std::string>>;
  using chi3_blocks_t = std::set<std::pair<std::string, std::string>>;

  /// Arguments of G2_iw_inu_inup()
  struct g2_iw_inu_inup_params_t {

    /// Structure of G^2 blocks.
    gf_struct_t gf_struct;

    /// Inverse temperature
    double beta;

    /// Channel in which Matsubara frequency representation is defined.
    channel_t channel = PH;

    /// Order of block indices in the definition of G^2.
    block_order_t block_order = AABB;

    /// List of block index pairs of G^2 to measure.
    /// default: measure all blocks
    g2_blocks_t blocks = g2_blocks_t{};

    /// Number of bosonic Matsubara frequencies.
    int n_iw = 30;

    /// Number of fermionic Matsubara frequencies.
    int n_inu = 30;

    /// Lehmann representation: Maximal distance between energy poles to be consider coinciding.
    double pole_res = 1e-8;

    /// Lehmann representation: Maximal magnitude of a term coefficient to be considered negligible.
    double coeff_tol = 1e-16;

    g2_iw_inu_inup_params_t() {}
    g2_iw_inu_inup_params_t(gf_struct_t const &gf_struct, double beta) : gf_struct(gf_struct), beta(beta) {}
  };

  /// Arguments of G2_iw_l_lp()
  struct g2_iw_l_lp_params_t {

    /// Structure of G^2 blocks.
    gf_struct_t gf_struct;

    /// Inverse temperature
    double beta;

    /// Channel in which Matsubara frequency representation is defined.
    channel_t channel = PH;

    /// Order of block indices in the definition of G^2.
    block_order_t block_order = AABB;

    /// List of block index pairs of G^2 to measure.
    /// default: measure all blocks
    g2_blocks_t blocks = g2_blocks_t{};

    /// Number of bosonic Matsubara frequencies.
    int n_iw = 30;

    /// Number of Legendre coefficients.
    int n_l = 20;

    /// Maximum number of positive Matsubara frequencies in summation.
    int n_inu_sum = 500;

    /// Tolerance for Matsubara frequency summation.
    double inu_sum_tol = 1e-6;

    /// Lehmann representation: Maximal distance between energy poles to be consider coinciding.
    double pole_res = 1e-8;

    /// Lehmann representation: Maximal magnitude of a term coefficient to be considered negligible.
    double coeff_tol = 1e-16;

    g2_iw_l_lp_params_t() {}
    g2_iw_l_lp_params_t(gf_struct_t const &gf_struct, double beta) : gf_struct(gf_struct), beta(beta) {}
  };

  /// Arguments of chi3_iw_inu()
  struct chi3_iw_inu_params_t {

    /// Structure of \chi^3 blocks.
    gf_struct_t gf_struct;

    /// Inverse temperature
    double beta;

    /// Channel in which Matsubara frequency representation is defined.
    channel_t channel = PH;

    /// Order of block indices in the definition of \chi^3.
    block_order_t block_order = AABB;

    /// List of block index pairs of \chi^3 to measure.
    /// default: measure all blocks
    chi3_blocks_t blocks = chi3_blocks_t{};

    /// Number of bosonic Matsubara frequencies.
    int n_iw = 100;

    /// Number of fermionic Matsubara frequencies.
    int n_inu = 100;

    /// Lehmann representation: Maximal distance between energy poles to be consider coinciding.
    double pole_res = 1e-8;

    /// Lehmann representation: Maximal magnitude of a term coefficient to be considered negligible.
    double coeff_tol = 1e-16;

    chi3_iw_inu_params_t() {}
    chi3_iw_inu_params_t(gf_struct_t const &gf_struct, double beta) : gf_struct(gf_struct), beta(beta) {}
  };

} // namespace pomerol2triqs
