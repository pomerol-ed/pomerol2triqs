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

#include "pomerol_ed.hpp"

//////////////////////////////////////////////////////////////////////////////////
// Methods to compute averages and correlation functions of quadratic operators //
//////////////////////////////////////////////////////////////////////////////////

namespace pomerol2triqs {

  std::complex<double> pomerol_ed::ensemble_average(indices_t const &i, indices_t const &j, double beta) {

    Pomerol::ParticleIndex pom_i = lookup_pomerol_index(i);
    if (pom_i == -1) TRIQS_RUNTIME_ERROR << "ensemble_average: Unexpected index i = " << i;
    Pomerol::ParticleIndex pom_j = lookup_pomerol_index(j);
    if (pom_j == -1) TRIQS_RUNTIME_ERROR << "ensemble_average: Unexpected index j = " << j;

    if (!matrix_h) TRIQS_RUNTIME_ERROR << "ensemble_average: No Hamiltonian has been diagonalized";
    compute_rho(beta);

    Pomerol::QuadraticOperator op(index_info, *hs, *states_class, *matrix_h, pom_i, pom_j);
    op.prepare(*hs);
    op.compute();

    Pomerol::EnsembleAverage EA(op, *rho);
    EA.compute();

    return EA();
  }

  template <typename Mesh, typename Filler>
  gf<Mesh, scalar_valued> pomerol_ed::compute_chi(indices_t const &i1, indices_t const &j1, indices_t const &i2, indices_t const &j2, bool connected,
                                                  Mesh const &mesh, Filler filler) const {
    if (!states_class || !matrix_h || !rho) TRIQS_RUNTIME_ERROR << "compute_chi: Internal error!";

    auto checked_lookup = [&](indices_t const &i) {
      Pomerol::ParticleIndex pom_i = lookup_pomerol_index(i);
      if (pom_i == -1) TRIQS_RUNTIME_ERROR << "compute_chi: Unexpected index " << i;
      return pom_i;
    };

    Pomerol::ParticleIndex pom_i1 = checked_lookup(i1);
    Pomerol::ParticleIndex pom_j1 = checked_lookup(j1);
    Pomerol::ParticleIndex pom_i2 = checked_lookup(i2);
    Pomerol::ParticleIndex pom_j2 = checked_lookup(j2);

    Pomerol::QuadraticOperator A(index_info, *hs, *states_class, *matrix_h, pom_i1, pom_j1);
    Pomerol::QuadraticOperator B(index_info, *hs, *states_class, *matrix_h, pom_i2, pom_j2);

    A.prepare(*hs);
    A.compute();
    B.prepare(*hs);
    B.compute();

    Pomerol::Susceptibility pom_chi(*states_class, *matrix_h, A, B, *rho);
    pom_chi.prepare();
    pom_chi.compute();
    if (connected) pom_chi.subtractDisconnected();

    gf<Mesh, scalar_valued> chi(mesh);
    filler(chi, pom_chi);

    return chi;
  }

  gf<imtime, scalar_valued> pomerol_ed::chi_tau(indices_t const &i1, indices_t const &j1, indices_t const &i2, indices_t const &j2, double beta,
                                                int n_tau, bool connected) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "chi_tau: No Hamiltonian has been diagonalized";
    compute_rho(beta);

    auto filler = [](gf_view<imtime, scalar_valued> chi, Pomerol::Susceptibility const &pom_chi) {
      for (auto tau : chi.mesh()) chi[tau] = pom_chi.of_tau(double(tau));
    };
    return compute_chi<imtime>(i1, j1, i2, j2, connected, {beta, Boson, n_tau}, filler);
  }

  gf<imfreq, scalar_valued> pomerol_ed::chi_inu(indices_t const &i1, indices_t const &j1, indices_t const &i2, indices_t const &j2, double beta,
                                                int n_inu, bool connected) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "chi_inu: No Hamiltonian has been diagonalized";
    compute_rho(beta);

    auto filler = [](gf_view<imfreq, scalar_valued> chi, Pomerol::Susceptibility const &pom_chi) {
      for (auto inu : chi.mesh()) chi[inu] = pom_chi(std::complex<double>(inu));
    };
    return compute_chi<imfreq>(i1, j1, i2, j2, connected, {beta, Boson, n_inu}, filler);
  }

} // namespace pomerol2triqs
