/**
 * pomerol2triqs
 *
 * Copyright (C) 2017-2019 Igor Krivenko <igor.s.krivenko @ gmail.com>
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

  std::complex<double> pomerol_ed::ensemble_average(indices_t const &i,
                                                    indices_t const &j,
                                                    double beta) {
    Pomerol::ParticleIndex pom_i = lookup_pomerol_index(i);
    if (pom_i == -1)
      TRIQS_RUNTIME_ERROR << "ensemble_average: unexpected first index " << i;
    Pomerol::ParticleIndex pom_j = lookup_pomerol_index(j);
    if (pom_j == -1)
      TRIQS_RUNTIME_ERROR << "ensemble_average: unexpected second index " << j;

    if (!matrix_h)
      TRIQS_RUNTIME_ERROR << "ensemble_average: no Hamiltonian has been diagonalized";
    compute_rho(beta);

    Pomerol::QuadraticOperator op(index_info, *states_class, *matrix_h, pom_i, pom_j);
    op.prepare();
    op.compute();

    Pomerol::EnsembleAverage EA(*states_class, *matrix_h, op, *rho);
    EA.prepare();

    return EA.getResult();
  }

} // namespace pomerol2triqs
