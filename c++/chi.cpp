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

#include <limits>

//////////////////////////////////////////////////////////////////////////////////
// Methods to compute averages and correlation functions of quadratic operators //
//////////////////////////////////////////////////////////////////////////////////

namespace pomerol2triqs {

  std::complex<double> pomerol_ed::ensemble_average(indices_t const &i, indices_t const &j,
                                                    double beta, std::tuple<bool, bool> const& dagger) {
    Pomerol::ParticleIndex pom_i = lookup_pomerol_index(i);
    if (pom_i == -1) TRIQS_RUNTIME_ERROR << "ensemble_average: Unexpected index i = " << i;
    Pomerol::ParticleIndex pom_j = lookup_pomerol_index(j);
    if (pom_j == -1) TRIQS_RUNTIME_ERROR << "ensemble_average: Unexpected index j = " << j;

    if (!matrix_h) TRIQS_RUNTIME_ERROR << "ensemble_average: No Hamiltonian has been diagonalized";
    compute_rho(beta);

    Pomerol::QuadraticOperator op(index_info, *hs, *states_class, *matrix_h, pom_i, pom_j, dagger);
    op.prepare(*hs);
    op.compute(ops_melem_tol);

    Pomerol::EnsembleAverage EA(op, *rho);
    EA.compute();

    return EA();
  }

  std::complex<double> pomerol_ed::ensemble_average(indices_t const &i, indices_t const &j, indices_t const &k, indices_t const &l,
                                                    double beta, std::tuple<bool, bool, bool, bool> const& dagger) {
    auto checked_lookup = [&](indices_t const &i, std::string const & i_name) {
      Pomerol::ParticleIndex pom_i = lookup_pomerol_index(i);
      if (pom_i == -1) TRIQS_RUNTIME_ERROR << "compute_chi: Unexpected index " << i_name << " = " << i;
      return pom_i;
    };

    Pomerol::ParticleIndex pom_i = checked_lookup(i, "i");
    Pomerol::ParticleIndex pom_j = checked_lookup(j, "j");
    Pomerol::ParticleIndex pom_k = checked_lookup(k, "k");
    Pomerol::ParticleIndex pom_l = checked_lookup(l, "l");

    if (!matrix_h) TRIQS_RUNTIME_ERROR << "ensemble_average: No Hamiltonian has been diagonalized";
    compute_rho(beta);

    Pomerol::QuarticOperator op(index_info, *hs, *states_class, *matrix_h, pom_i, pom_j, pom_k, pom_l, dagger);
    op.prepare(*hs);
    op.compute(ops_melem_tol);

    Pomerol::EnsembleAverage EA(op, *rho);
    EA.compute();

    return EA();
  }

  template <typename Mesh, typename Filler>
  gf<Mesh, scalar_valued> pomerol_ed::compute_chi(indices_t const &i, indices_t const &j, indices_t const &k, indices_t const &l, bool connected,
                                                  Mesh const &mesh, Filler filler, channel_t channel,
                                                  double pole_res, double coeff_tol) const {
    if (!states_class || !matrix_h || !rho) TRIQS_RUNTIME_ERROR << "compute_chi: Internal error!";

    auto checked_lookup = [&](indices_t const &i) {
      Pomerol::ParticleIndex pom_i = lookup_pomerol_index(i);
      if (pom_i == -1) TRIQS_RUNTIME_ERROR << "compute_chi: Unexpected index " << i;
      return pom_i;
    };

    Pomerol::ParticleIndex pom_i = checked_lookup(i);
    Pomerol::ParticleIndex pom_j = checked_lookup(j);
    Pomerol::ParticleIndex pom_k = checked_lookup(k);
    Pomerol::ParticleIndex pom_l = checked_lookup(l);

    std::tuple<bool, bool> adag;
    std::tuple<bool, bool> bdag;
    Pomerol::ParticleIndex a1, a2;
    Pomerol::ParticleIndex b1, b2;
    double sign;

    switch(channel) {
      case PP:
        adag = {true, true};
        bdag = {false, false};
        a1 = pom_i;
        a2 = pom_k;
        b1 = pom_j;
        b2 = pom_l;
        sign = -1.0;
        break;
      case PH:
        adag = {true, false};
        bdag = {true, false};
        a1 = pom_i;
        a2 = pom_j;
        b1 = pom_k;
        b2 = pom_l;
        sign = 1.0;
        break;
      case xPH:
        adag = {true, false};
        bdag = {true, false};
        a1 = pom_i;
        a2 = pom_l;
        b1 = pom_k;
        b2 = pom_j;
        sign = -1.0;
        break;
      default:
        TRIQS_RUNTIME_ERROR << "compute_chi: AllFermionic channel is not supported";
    }

    Pomerol::QuadraticOperator A(index_info, *hs, *states_class, *matrix_h, a1, a2, adag);
    Pomerol::QuadraticOperator B(index_info, *hs, *states_class, *matrix_h, b1, b2, bdag);

    A.prepare(*hs);
    A.compute(ops_melem_tol);
    B.prepare(*hs);
    B.compute(ops_melem_tol);

    Pomerol::Susceptibility pom_chi(*states_class, *matrix_h, A, B, *rho);
    pom_chi.PoleResolution = pole_res;
    pom_chi.CoefficientTolerance = coeff_tol;
    pom_chi.prepare();
    pom_chi.compute();
    if (connected) pom_chi.subtractDisconnected();

    gf<Mesh, scalar_valued> chi(mesh);
    filler(chi, pom_chi);
    chi *= sign;

    return chi;
  }

  gf<imtime, scalar_valued> pomerol_ed::chi_tau(indices_t const &i, indices_t const &j, indices_t const &k, indices_t const &l, double beta,
                                                int n_tau, bool connected, channel_t channel,
                                                double pole_res, double coeff_tol) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "chi_tau: No Hamiltonian has been diagonalized";
    compute_rho(beta);

    auto filler = [](gf_view<imtime, scalar_valued> chi, Pomerol::Susceptibility const &pom_chi) {
      for (auto tau : chi.mesh()) chi[tau] = pom_chi.of_tau(double(tau));
    };
    return compute_chi<imtime>(i, j, k, l, connected, {beta, Boson, n_tau}, filler, channel, pole_res, coeff_tol);
  }

  gf<imfreq, scalar_valued> pomerol_ed::chi_iw(indices_t const &i, indices_t const &j, indices_t const &k, indices_t const &l, double beta,
                                               int n_iw, bool connected, channel_t channel,
                                               double pole_res, double coeff_tol) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "chi_iw: No Hamiltonian has been diagonalized";
    compute_rho(beta);

    auto filler = [](gf_view<imfreq, scalar_valued> chi, Pomerol::Susceptibility const &pom_chi) {
      for (auto inu : chi.mesh()) chi[inu] = pom_chi(std::complex<double>(inu));
    };
    return compute_chi<imfreq>(i, j, k, l, connected, {beta, Boson, n_iw}, filler, channel, pole_res, coeff_tol);
  }

  gf<refreq, scalar_valued> pomerol_ed::chi_w(indices_t const &i, indices_t const &j, indices_t const &k, indices_t const &l, double beta,
                                              std::pair<double, double> const &energy_window, int n_w, double im_shift,
                                              bool connected, channel_t channel,
                                              double pole_res, double coeff_tol) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "chi_w: No Hamiltonian has been diagonalized";
    compute_rho(beta);

    auto filler = [im_shift](gf_view<refreq, scalar_valued> chi, Pomerol::Susceptibility const &pom_chi) {
      for (auto w : chi.mesh()) {
        auto w_val = std::complex<double>(w);
        // The imaginary shift is not applied to w = 0 so that the result at this special point
        // coincides with that at the zero Matsubara frequency
        chi[w] = (std::abs(std::real(w_val)) <= std::numeric_limits<double>::epsilon()) ? pom_chi(.0) : pom_chi(w_val + 1i * im_shift);
      }
    };
    return compute_chi<refreq>(i, j, k, l, connected, {energy_window, n_w}, filler, channel, pole_res, coeff_tol);
  }

} // namespace pomerol2triqs
