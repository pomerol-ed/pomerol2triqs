/**
 * pomerol2triqs
 *
 * Copyright (C) 2017-2021 Igor Krivenko <igor.s.krivenko @ gmail.com>
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

namespace pomerol2triqs {

  pomerol_ed::pomerol_ed(index_converter_t const &index_converter, bool verbose)
     : verbose(verbose), index_converter(index_converter) {
    for (auto const &ind : index_converter) {
      index_info.addInfo(ind.second);
    }
    if (verbose && !pMPI::rank(comm)) {
      std::cout << "Pomerol: operator indices" << std::endl;
      std::cout << index_info << std::endl;
    }
  }

  template<typename HExprType>
  void pomerol_ed::diagonalize_prepare_impl(many_body_op_t const &hamiltonian) {
    h_expr.reset(new h_expr_t(HExprType()));

    for (auto const &term : hamiltonian) {
      HExprType pom_term(static_cast<typename HExprType::scalar_type>(term.coef));
      for(auto const& o : term.monomial) {
        auto it = index_converter.find(o.indices);
        if (it == index_converter.end())
          TRIQS_RUNTIME_ERROR << "diagonalize: Invalid Hamiltonian, unexpected operator indices " << o.indices;

        pom_term = pom_term * (o.dagger ?
            Pomerol::Operators::c_dag(std::get<0>(it->second), (unsigned short)std::get<1>(it->second), std::get<2>(it->second)) :
            Pomerol::Operators::c(std::get<0>(it->second), (unsigned short)std::get<1>(it->second), std::get<2>(it->second))
          );
      }
      std::get<HExprType>(*h_expr) += pom_term;
    }

    if (verbose && !pMPI::rank(comm)) {
      std::cout << "Pomerol: Hamiltonian" << std::endl;
      std::cout << *h_expr << std::endl;
    }
  }

  void pomerol_ed::diagonalize_prepare(many_body_op_t const &hamiltonian) {
    using namespace Pomerol::LatticePresets;

    if(std::all_of(hamiltonian.cbegin(),
                   hamiltonian.cend(),
                   [](auto const& term) { return term.coef.is_real(); })
    )
      diagonalize_prepare_impl<RealExpr>(hamiltonian);
    else
      diagonalize_prepare_impl<ComplexExpr>(hamiltonian);
  }

  void pomerol_ed::diagonalize(many_body_op_t const &hamiltonian, bool ignore_symmetries) {
    diagonalize_prepare(hamiltonian);

    // Create Hilbert space
    std::visit([&](auto const& h) { hs.reset(new hilbert_space_t(index_info, h)); },
               *h_expr);
    if(!ignore_symmetries)
      hs->compute();

    // Classify many-body states
    states_class.reset(new Pomerol::StatesClassification());
    states_class->compute(*hs);

    // Matrix representation of the Hamiltonian
    matrix_h.reset(new Pomerol::Hamiltonian(*states_class));
    std::visit([&](auto const& h) { matrix_h->prepare(h, *hs, comm); }, *h_expr);
    matrix_h->compute(comm);

    // Get ground state energy
    if (verbose && !pMPI::rank(comm))
      std::cout << "Pomerol: Ground state energy is " << matrix_h->getGroundEnergy() << std::endl;

    // Reset containers, we will compute them later if needed
    rho.release();
    ops_container.release();
  }

  Pomerol::ParticleIndex pomerol_ed::lookup_pomerol_index(indices_t const &i) const {
    auto it = index_converter.find(i);
    if (it == index_converter.end()) return -1;
    return index_info.getIndex(it->second);
  }

  // Translate gf_struct into a set of ParticleIndex
  std::set<Pomerol::ParticleIndex> pomerol_ed::gf_struct_to_pomerol_indices(gf_struct_t const &gf_struct) const {
    std::set<Pomerol::ParticleIndex> indices;
    for (auto const &b : gf_struct) {
      for (auto const &i : b.second) {
        auto pom_ind = lookup_pomerol_index({b.first, i});
        if (pom_ind == -1) TRIQS_RUNTIME_ERROR << "gf_struct_to_pomerol_indices: Unexpected GF index " << b.first << "," << i;
        indices.insert(pom_ind);
      }
    }
    return indices;
  }

  // Create the Density Matrix.
  void pomerol_ed::compute_rho(double beta) {
    if (!states_class || !matrix_h) TRIQS_RUNTIME_ERROR << "compute_rho: Internal error!";

    if (!rho || rho->beta != beta) {
      if (verbose && !pMPI::rank(comm))
        std::cout << "Pomerol: Computing density matrix for \\beta = " << beta << std::endl;
      rho.reset(new Pomerol::DensityMatrix(*states_class, *matrix_h, beta));
      rho->prepare();
      rho->compute();
      if(rho_threshold > 0) rho->truncateBlocks(rho_threshold, verbose);
    }
  }

  void pomerol_ed::compute_field_operators(gf_struct_t const &gf_struct) {
    if (!states_class || !matrix_h) TRIQS_RUNTIME_ERROR << "compute_field_operators: Internal error!";

    auto new_ops = gf_struct_to_pomerol_indices(gf_struct);
    if (!ops_container || computed_ops != new_ops) {
      if (verbose && !pMPI::rank(comm)) {
        std::cout << "Pomerol: Computing field operators with indices ";
        bool comma = false;
        for (auto i : new_ops) {
          std::cout << (comma ? ", " : "") << i;
          comma = true;
        }
        std::cout << std::endl;
      }

      ops_container.reset(new Pomerol::FieldOperatorContainer(index_info, *hs, *states_class, *matrix_h));
      ops_container->prepareAll(*hs);
      ops_container->computeAll();

      computed_ops = new_ops;
    }
  }

} // namespace pomerol2triqs
