/**
 * pomerol2triqs
 *
 * Copyright (C) 2017-2025 Igor Krivenko
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

#include "pomerol_ed.hpp"

///////////////////////////////////////////////////////
// Methods to compute one-particle Green's functions //
///////////////////////////////////////////////////////

namespace pomerol2triqs {

  template <typename Mesh, typename Filler>
  block_gf<Mesh> pomerol_ed::compute_gf(gf_struct_t const &gf_struct, Mesh const &mesh, Filler filler, bool anomalous,
                                        double pole_res, double coeff_tol) const {

    if (!states_class || !matrix_h || !rho || !ops_container) TRIQS_RUNTIME_ERROR << "compute_gf: Internal error!";

    std::vector<std::string> block_names;
    std::vector<gf<Mesh>> g_blocks;

    for (auto const &bl : gf_struct) {
      block_names.push_back(bl.first);
      int n = bl.second;

      g_blocks.push_back(gf<Mesh>{mesh, {n, n}});

      auto &g = g_blocks.back();

      for (int i1 : range(n)) {
        Pomerol::ParticleIndex pom_i1 = lookup_pomerol_index({bl.first, i1});
        for (int i2 : range(n)) {
          Pomerol::ParticleIndex pom_i2 = lookup_pomerol_index({bl.first, i2});

          if (verbose && !comm.rank())
            std::cout << "fill_gf: Filling GF component (" << bl.first << "," << i1 << ")(" << bl.first << "," << i2 << ")" << std::endl;
          auto g_el = slice_target_to_scalar(g, i1, i2);

          auto pom_g = anomalous ? Pomerol::GreensFunction(*states_class,
                                                           *matrix_h,
                                                           ops_container->getAnnihilationOperator(pom_i1),
                                                           ops_container->getAnnihilationOperator(pom_i2),
                                                           *rho) :
                                   Pomerol::GreensFunction(*states_class,
                                                           *matrix_h,
                                                           ops_container->getAnnihilationOperator(pom_i1),
                                                           ops_container->getCreationOperator(pom_i2),
                                                           *rho);
          pom_g.PoleResolution = pole_res;
          pom_g.CoefficientTolerance = coeff_tol;
          pom_g.prepare();
          pom_g.compute();

          filler(g_el, pom_g);
        }
      }
    }
    return make_block_gf(block_names, std::move(g_blocks));
  }

  block_gf<imfreq> pomerol_ed::G_iw(gf_struct_t const &gf_struct, double beta, int n_iw, double pole_res, double coeff_tol) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "G_iw: No Hamiltonian has been diagonalized";
    compute_rho(beta);
    compute_field_operators(gf_struct);

    auto filler = [](gf_view<imfreq, scalar_valued> g_el, Pomerol::GreensFunction const &pom_g) {
      for (auto iw : g_el.mesh()) g_el[iw] = pom_g(std::complex<double>(iw));
    };
    return compute_gf<imfreq>(gf_struct, {beta, Fermion, n_iw}, filler, false, pole_res, coeff_tol);
  }

  block_gf<imtime> pomerol_ed::G_tau(gf_struct_t const &gf_struct, double beta, int n_tau, double pole_res, double coeff_tol) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "G_tau: No Hamiltonian has been diagonalized";
    compute_rho(beta);
    compute_field_operators(gf_struct);

    auto filler = [](gf_view<imtime, scalar_valued> g_el, Pomerol::GreensFunction const &pom_g) {
      for (auto tau : g_el.mesh()) g_el[tau] = pom_g.of_tau(tau);
    };
    return compute_gf<imtime>(gf_struct, {beta, Fermion, n_tau}, filler, false, pole_res, coeff_tol);
  }

  block_gf<refreq> pomerol_ed::G_w(gf_struct_t const &gf_struct, double beta, std::pair<double, double> const &energy_window, int n_w,
                                   double im_shift, double pole_res, double coeff_tol) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "G_w: No Hamiltonian has been diagonalized";
    compute_rho(beta);
    compute_field_operators(gf_struct);

    auto filler = [im_shift](gf_view<refreq, scalar_valued> g_el, Pomerol::GreensFunction const &pom_g) {
      for (auto w : g_el.mesh()) g_el[w] = pom_g(double(w) + 1i * im_shift);
    };
    return compute_gf<refreq>(gf_struct, {energy_window.first, energy_window.second, n_w}, filler, false, pole_res, coeff_tol);
  }

  block_gf<imfreq> pomerol_ed::F_iw(gf_struct_t const &gf_struct, double beta, int n_iw, double pole_res, double coeff_tol) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "F_iw: No Hamiltonian has been diagonalized";
    compute_rho(beta);
    compute_field_operators(gf_struct);

    auto filler = [](gf_view<imfreq, scalar_valued> g_el, Pomerol::GreensFunction const &pom_g) {
      for (auto iw : g_el.mesh()) g_el[iw] = pom_g(std::complex<double>(iw));
    };
    return compute_gf<imfreq>(gf_struct, {beta, Fermion, n_iw}, filler, true, pole_res, coeff_tol);
  }

  block_gf<imtime> pomerol_ed::F_tau(gf_struct_t const &gf_struct, double beta, int n_tau, double pole_res, double coeff_tol) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "F_tau: No Hamiltonian has been diagonalized";
    compute_rho(beta);
    compute_field_operators(gf_struct);

    auto filler = [](gf_view<imtime, scalar_valued> g_el, Pomerol::GreensFunction const &pom_g) {
      for (auto tau : g_el.mesh()) g_el[tau] = pom_g.of_tau(tau);
    };
    return compute_gf<imtime>(gf_struct, {beta, Fermion, n_tau}, filler, true, pole_res, coeff_tol);
  }

  block_gf<refreq> pomerol_ed::F_w(gf_struct_t const &gf_struct, double beta, std::pair<double, double> const &energy_window, int n_w,
                                      double im_shift, double pole_res, double coeff_tol) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "F_w: No Hamiltonian has been diagonalized";
    compute_rho(beta);
    compute_field_operators(gf_struct);

    auto filler = [im_shift](gf_view<refreq, scalar_valued> g_el, Pomerol::GreensFunction const &pom_g) {
      for (auto w : g_el.mesh()) g_el[w] = pom_g(double(w) + 1i * im_shift);
    };
    return compute_gf<refreq>(gf_struct, {energy_window.first, energy_window.second, n_w}, filler, true, pole_res, coeff_tol);
  }

} // namespace pomerol2triqs
