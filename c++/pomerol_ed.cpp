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

#include "pomerol_ed.hpp"

#include <algorithm>

namespace pomerol2triqs {

  pomerol_ed::pomerol_ed(index_converter_t const &index_converter, bool verbose) : verbose(verbose), index_converter(index_converter) {

    if (!mpi::is_initialized()) TRIQS_RUNTIME_ERROR << "pomerol2triqs does not support running in the serial (no MPI) mode";

    for (auto const &ind : index_converter) { index_info.addInfo(ind.second); }
    if (verbose && !comm.rank()) {
      std::cout << "Pomerol: operator indices" << std::endl;
      std::cout << index_info << std::endl;
    }
  }

  template <typename HExprType> void pomerol_ed::diagonalize_prepare_impl(many_body_op_t const &hamiltonian) {
    h_expr.reset(new h_expr_t(HExprType()));

    for (auto const &term : hamiltonian) {
      HExprType pom_term(static_cast<typename HExprType::scalar_type>(term.coef));
      for (auto const &o : term.monomial) {
        auto it = index_converter.find(o.indices);
        if (it == index_converter.end()) TRIQS_RUNTIME_ERROR << "diagonalize: Invalid Hamiltonian, unexpected operator indices " << o.indices;

        pom_term = pom_term
           * (o.dagger ? Pomerol::Operators::c_dag(std::get<0>(it->second), (unsigned short)std::get<1>(it->second), std::get<2>(it->second)) :
                         Pomerol::Operators::c(std::get<0>(it->second), (unsigned short)std::get<1>(it->second), std::get<2>(it->second)));
      }
      std::get<HExprType>(*h_expr) += pom_term;
    }

    if (verbose && !comm.rank()) {
      std::cout << "Pomerol: Hamiltonian" << std::endl;
      std::cout << *h_expr << std::endl;
    }
  }

  void pomerol_ed::diagonalize_prepare(many_body_op_t const &hamiltonian) {
    using namespace Pomerol::LatticePresets;

    if (std::all_of(hamiltonian.cbegin(), hamiltonian.cend(), [](auto const &term) { return term.coef.is_real(); }))
      diagonalize_prepare_impl<RealExpr>(hamiltonian);
    else
      diagonalize_prepare_impl<ComplexExpr>(hamiltonian);
  }

  void pomerol_ed::diagonalize(many_body_op_t const &hamiltonian, bool ignore_symmetries) {
    diagonalize_prepare(hamiltonian);

    // Create Hilbert space
    std::visit([&](auto const &h) { hs.reset(new hilbert_space_t(index_info, h)); }, *h_expr);
    if (!ignore_symmetries) hs->compute();

    // Classify many-body states
    states_class.reset(new Pomerol::StatesClassification());
    states_class->compute(*hs);

    // Matrix representation of the Hamiltonian
    matrix_h.reset(new Pomerol::Hamiltonian(*states_class));
    std::visit([&](auto const &h) { matrix_h->prepare(h, *hs, comm.get()); }, *h_expr);
    matrix_h->compute(comm.get());

    // Get ground state energy
    if (verbose && !comm.rank()) std::cout << "Pomerol: Ground state energy is " << matrix_h->getGroundEnergy() << std::endl;

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
      for (auto const &i : range(b.second)) {
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
      if (verbose && !comm.rank()) std::cout << "Pomerol: Computing density matrix for \\beta = " << beta << std::endl;
      rho.reset(new Pomerol::DensityMatrix(*states_class, *matrix_h, beta));
      rho->prepare();
      rho->compute();
      if (rho_threshold > 0) rho->truncateBlocks(rho_threshold, verbose);
    }
  }

  void pomerol_ed::compute_field_operators(gf_struct_t const &gf_struct) {
    if (!states_class || !matrix_h) TRIQS_RUNTIME_ERROR << "compute_field_operators: Internal error!";

    auto new_ops = gf_struct_to_pomerol_indices(gf_struct);
    if (!ops_container || computed_ops != new_ops) {
      if (verbose && !comm.rank()) {
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
      ops_container->computeAll(ops_melem_tol);

      computed_ops = new_ops;
    }
  }

  std::uint64_t pomerol_ed::get_full_hilbert_space_dim() const {
    if (!states_class)
      TRIQS_RUNTIME_ERROR << "get_full_hilbert_space_dim: No Hamiltonian has been diagonalized";
    return states_class->getNumberOfStates();
  }

  std::uint64_t pomerol_ed::get_n_subspaces() const {
    if (!states_class)
      TRIQS_RUNTIME_ERROR << "get_n_subspaces: No Hamiltonian has been diagonalized";
    return states_class->getNumberOfBlocks();
  }

  std::vector<std::uint64_t> pomerol_ed::get_subspace_dims() const {
    if (!states_class)
      TRIQS_RUNTIME_ERROR << "get_subspace_dims: No Hamiltonian has been diagonalized";
    auto n_subspaces = states_class->getNumberOfBlocks();
    std::vector<std::uint64_t> dims(n_subspaces);
    for (auto sp : range(n_subspaces))
      dims[sp] = states_class->getBlockSize(sp);
    return dims;
  }

  std::uint64_t pomerol_ed::get_subspace_dim(std::uint64_t sp) const {
    if (!states_class)
      TRIQS_RUNTIME_ERROR << "get_subspace_dim: No Hamiltonian has been diagonalized";
    return states_class->getBlockSize(sp);
  }

  std::vector<std::vector<std::uint64_t>> pomerol_ed::get_fock_states() const {
    if (!states_class)
      TRIQS_RUNTIME_ERROR << "get_fock_states: No Hamiltonian has been diagonalized";
    auto n_subspaces = states_class->getNumberOfBlocks();
    std::vector<std::vector<std::uint64_t>> fock_states;
    fock_states.reserve(n_subspaces);
    for (auto sp : range(n_subspaces))
      fock_states.emplace_back(std::move(states_class->getFockStates(sp)));
    return fock_states;
  }

  std::vector<std::uint64_t> pomerol_ed::get_subspace_fock_states(std::uint64_t sp) const {
    if (!states_class)
      TRIQS_RUNTIME_ERROR << "get_subspace_fock_states: No Hamiltonian has been diagonalized";
    return states_class->getFockStates(sp);
  }

  std::vector<nda::vector<double>> pomerol_ed::get_energies() const {
    if (!states_class || !matrix_h)
      TRIQS_RUNTIME_ERROR << "get_energies: No Hamiltonian has been diagonalized";
    auto n_subspaces = states_class->getNumberOfBlocks();
    std::vector<nda::vector<double>> res;
    res.reserve(n_subspaces);
    for (auto sp : range(n_subspaces)) {
      res.emplace_back(get_subspace_energies(sp));
    }
    return res;
  }

  nda::vector<double> pomerol_ed::get_subspace_energies(std::uint64_t sp) const {
    if (!states_class || !matrix_h)
      TRIQS_RUNTIME_ERROR << "get_subspace_energies: No Hamiltonian has been diagonalized";
    auto const& eigenvalues = matrix_h->getEigenValues(sp);
    nda::vector<double> res(eigenvalues.size());
    std::copy(eigenvalues.begin(), eigenvalues.end(), res.begin());
    return res;
  }

  std::vector<pomerol_ed::rc_matrix_t> pomerol_ed::get_unitary_matrices() const {
    if (!states_class || !matrix_h)
      TRIQS_RUNTIME_ERROR << "get_unitary_matrices: No Hamiltonian has been diagonalized";
    auto n_subspaces = states_class->getNumberOfBlocks();
    std::vector<rc_matrix_t> res;
    res.reserve(n_subspaces);
    for (auto sp : range(n_subspaces)) {
      res.emplace_back(get_subspace_unitary_matrix(sp));
    }
    return res;
  }

  pomerol_ed::rc_matrix_t pomerol_ed::get_subspace_unitary_matrix(std::uint64_t sp) const {
    if (!states_class || !matrix_h)
      TRIQS_RUNTIME_ERROR << "get_subspace_unitary_matrix: No Hamiltonian has been diagonalized";
    auto const& h_part = matrix_h->getPart(sp);
    auto size = h_part.getSize();
    if(h_part.isComplex()) {
      nda::matrix<dcomplex> res(size, size);
      auto const& mat = h_part.getMatrix<true>();
      for (auto i : range(size)) {
        for (auto j : range(size)) {
          res(i, j) = mat(i, j);
        }
      }
      return res;
    } else {
      nda::matrix<double> res(size, size);
      auto const& mat = h_part.getMatrix<false>();
      for (auto i : range(size)) {
        for (auto j : range(size)) {
          res(i, j) = mat(i, j);
        }
      }
      return res;
    }
  }

} // namespace pomerol2triqs
