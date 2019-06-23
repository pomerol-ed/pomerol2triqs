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

namespace pomerol2triqs {

  Pomerol::Lattice pomerol_ed::init() {

    Pomerol::Lattice l;

    std::map<std::string, int> site_max_orb;
    for (auto const &ind : index_converter) {
      std::string pomerol_site;
      int pomerol_orb;
      std::tie(pomerol_site, pomerol_orb, std::ignore) = ind.second;

      auto it = site_max_orb.find(pomerol_site);
      if (it == site_max_orb.end())
        site_max_orb[pomerol_site] = pomerol_orb;
      else
        it->second = std::max(it->second, pomerol_orb);
    }

    for (auto const &site_orb : site_max_orb) l.addSite(new Pomerol::Lattice::Site(site_orb.first, site_orb.second + 1, 2));

    return l;
  }

  pomerol_ed::pomerol_ed(index_converter_t const &index_converter, bool verbose)
     : verbose(verbose), index_converter(index_converter), bare_lattice(init()), index_info(bare_lattice.getSiteMap()) {
    index_info.prepare();
    if (verbose && !comm.rank()) {
      std::cout << "Pomerol: lattice sites" << std::endl;
      bare_lattice.printSites();
      std::cout << "Pomerol: operator indices" << std::endl;
      index_info.printIndices();
    }
  }

  double pomerol_ed::diagonalize_prepare(many_body_op_t const &hamiltonian) {
    // Workaround for the broken std::vector<bool>
    struct bool_ {
      bool b;
    };

    std::vector<bool_> OperatorSequence;
    std::vector<std::string> SiteLabels;
    std::vector<unsigned short> Orbitals;
    std::vector<unsigned short> Spins;

    lattice.reset(new Pomerol::Lattice(bare_lattice));

    double gs_shift = 0;

    for (auto const &term : hamiltonian) {
      if (term.monomial.empty()) {
        gs_shift = std::real(term.coef);
        continue; // Constant term is unphysical anyway ...
      }
      OperatorSequence.clear();
      SiteLabels.clear();
      Orbitals.clear();
      Spins.clear();

      for (auto o : term.monomial) {
        auto it = index_converter.find(o.indices);
        if (it == index_converter.end()) TRIQS_RUNTIME_ERROR << "diagonalize: invalid Hamiltonian, unexpected operator indices " << o.indices;

        OperatorSequence.push_back({o.dagger});

        std::string site;
        unsigned short orb;
        Pomerol::spin s;
        std::tie(site, orb, s) = it->second;
        SiteLabels.push_back(site);
        Orbitals.push_back(orb);
        Spins.push_back(s);
      }

      lattice->addTerm(new Pomerol::Lattice::Term(term.monomial.size(), reinterpret_cast<bool *>(OperatorSequence.data()), term.coef,
                                                  SiteLabels.data(), Orbitals.data(), Spins.data()));
    }

    storage.reset(new Pomerol::IndexHamiltonian(lattice.get(), index_info));
    storage->prepare();

    if (verbose && !comm.rank()) {
      std::cout << "Pomerol: terms of Hamiltonian" << std::endl;
      std::cout << *storage << std::endl;
    }

    return gs_shift;
  }

  void pomerol_ed::diagonalize_main(double gs_shift) {

    // Classify many-body states
    states_class.reset(new Pomerol::StatesClassification(index_info, *symm));
    states_class->compute();

    // Matrix representation of the Hamiltonian
    matrix_h.reset(new Pomerol::Hamiltonian(index_info, *storage, *states_class));
    matrix_h->prepare(comm);
    matrix_h->compute(comm);

    // Get ground state energy
    if (verbose && !comm.rank()) { std::cout << "Pomerol: ground state energy is " << matrix_h->getGroundEnergy() + gs_shift << std::endl; }

    // Reset containers, we will compute them later if needed
    rho.release();
    ops_container.release();
  }

  void pomerol_ed::diagonalize(many_body_op_t const &hamiltonian, bool ignore_symmetries) {

    double gs_shift = diagonalize_prepare(hamiltonian);

    // Check the Hamiltonian commutes with the total number of particles
    Pomerol::OperatorPresets::N N(index_info.getIndexSize());
    if (!storage->commutes(N)) TRIQS_RUNTIME_ERROR << "diagonalize: Hamiltonian does not conserve the total number of particles";

    // Construct Symmetrizer
    symm.reset(new Pomerol::Symmetrizer(index_info, *storage));
    symm->compute(ignore_symmetries);

    diagonalize_main(gs_shift);
  }

  void pomerol_ed::diagonalize(many_body_op_t const &hamiltonian, std::vector<many_body_op_t> const& integrals_of_motion) {

    double gs_shift = diagonalize_prepare(hamiltonian);

    std::vector<Pomerol::Operator> iom;
    for(auto const& op : integrals_of_motion) {
      Pomerol::Operator pom_op;
      for (auto const &term : op) {
        Pomerol::Operator pom_term;
        pom_term += term.coef;
        for (auto o : term.monomial) {
          Pomerol::ParticleIndex pom_ind = lookup_pomerol_index(o.indices);

          if (pom_ind == -1)
            TRIQS_RUNTIME_ERROR << "diagonalize: invalid integral of motion, unexpected operator indices " << o.indices;

          if(o.dagger)
            pom_term *= Pomerol::OperatorPresets::c_dag(pom_ind);
          else
            pom_term *= Pomerol::OperatorPresets::c(pom_ind);
        }

        pom_op += pom_term;
      }
      iom.push_back(pom_op);
    }

    // Construct Symmetrizer
    symm.reset(new Pomerol::Symmetrizer(index_info, *storage));
    symm->compute(iom);

    diagonalize_main(gs_shift);
  }

  Pomerol::ParticleIndex pomerol_ed::lookup_pomerol_index(indices_t const &i) const {
    auto it = index_converter.find(i);
    if (it == index_converter.end()) return -1;
    std::string site;
    unsigned short orb;
    Pomerol::spin s;
    std::tie(site, orb, s) = it->second;
    return index_info.getIndex(site, orb, s);
  }

  // Translate gf_struct into a set of ParticleIndex
  std::set<Pomerol::ParticleIndex> pomerol_ed::gf_struct_to_pomerol_indices(gf_struct_t const &gf_struct) const {
    std::set<Pomerol::ParticleIndex> indices;
    for (auto const &b : gf_struct) {
      for (auto const &i : b.second) {
        auto pom_ind = lookup_pomerol_index({b.first, i});
        if (pom_ind == -1) TRIQS_RUNTIME_ERROR << "gf_struct_to_pomerol_indices: unexpected GF index " << b.first << "," << i;
        indices.insert(pom_ind);
      }
    }
    return indices;
  }

  // Create the Density Matrix.
  void pomerol_ed::compute_rho(double beta) {
    if (!states_class || !matrix_h) TRIQS_RUNTIME_ERROR << "compute_rho: internal error!";

    if (!rho || rho->beta != beta) {
      if (verbose && !comm.rank()) std::cout << "Pomerol: computing density matrix for \\beta = " << beta << std::endl;
      rho.reset(new Pomerol::DensityMatrix(*states_class, *matrix_h, beta));
      rho->prepare();
      rho->compute();
      if(rho_threshold > 0) rho->truncateBlocks(rho_threshold, verbose);
    }
  }

  void pomerol_ed::compute_field_operators(gf_struct_t const &gf_struct) {
    if (!states_class || !matrix_h) TRIQS_RUNTIME_ERROR << "compute_field_operators: internal error!";

    auto new_ops = gf_struct_to_pomerol_indices(gf_struct);
    if (!ops_container || computed_ops != new_ops) {
      if (verbose && !comm.rank()) {
        std::cout << "Pomerol: computing field operators with indices ";
        bool comma = false;
        for (auto i : new_ops) {
          std::cout << (comma ? ", " : "") << i;
          comma = true;
        }
        std::cout << std::endl;
      }

      ops_container.reset(new Pomerol::FieldOperatorContainer(index_info, *states_class, *matrix_h));
      ops_container->prepareAll(new_ops);
      ops_container->computeAll();

      computed_ops = new_ops;
    }
  }
}
