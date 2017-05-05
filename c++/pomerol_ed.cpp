#include "pomerol_ed.hpp"

namespace pomerol2triqs {

Pomerol::Lattice pomerol_ed::init() {

 Pomerol::Lattice l;

 std::map<std::string, int> site_max_orb;
 for(auto const& ind : index_converter) {
  std::string pomerol_site;
  int pomerol_orb;
  std::tie(pomerol_site, pomerol_orb, std::ignore) = ind.second;

  auto it = site_max_orb.find(pomerol_site);
  if(it == site_max_orb.end())
   site_max_orb[pomerol_site] = pomerol_orb;
  else
   it->second = std::max(it->second, pomerol_orb);
 }

 for(auto const& site_orb : site_max_orb)
  l.addSite(new Pomerol::Lattice::Site(site_orb.first, site_orb.second+1, 2));

 return l;
}

pomerol_ed::pomerol_ed(index_converter_t const& index_converter, bool verbose) :
 verbose(verbose),
 index_converter(index_converter),
 bare_lattice(init()),
 index_info(bare_lattice.getSiteMap()) {
 index_info.prepare();
 if(verbose && !comm.rank()) {
  std::cout << "Pomerol: lattice sites" << std::endl;
  bare_lattice.printSites();
  std::cout << "Pomerol: operator indices" << std::endl;
  index_info.printIndices();
 }
}

void pomerol_ed::diagonalize(many_body_op_t const& hamiltonian, bool ignore_symmetries) {

 // Workaround for the broken std::vector<bool>
 struct bool_ { bool b; };

 std::vector<bool_> OperatorSequence;
 std::vector<std::string> SiteLabels;
 std::vector<unsigned short> Orbitals;
 std::vector<unsigned short> Spins;

 lattice.reset(new Pomerol::Lattice(bare_lattice));

 double gs_shift = 0;

 for(auto const& term : hamiltonian) {
  if(term.monomial.empty()) {
   gs_shift = term.coef.real();
   continue;    // Constant term is unphysical anyway ...
  }
  OperatorSequence.clear();
  SiteLabels.clear();
  Orbitals.clear();
  Spins.clear();

  for(auto o : term.monomial) {
   auto it = index_converter.find(o.indices);
   if(it == index_converter.end())
    TRIQS_RUNTIME_ERROR << "diagonalize: invalid Hamiltonian, unexpected operator indices "
                        << o.indices;

   OperatorSequence.push_back({o.dagger});

   std::string site;
   unsigned short orb;
   Pomerol::spin s;
   std::tie(site, orb, s) = it->second;
   SiteLabels.push_back(site);
   Orbitals.push_back(orb);
   Spins.push_back(s);
  }

  lattice->addTerm(new Pomerol::Lattice::Term(term.monomial.size(),
                                              reinterpret_cast<bool*>(OperatorSequence.data() ),
                                              term.coef,
                                              SiteLabels.data(),
                                              Orbitals.data(),
                                              Spins.data()));
 }

 storage.reset(new Pomerol::IndexHamiltonian(lattice.get(), index_info));
 storage->prepare();

 if(verbose && !comm.rank()) {
  std::cout << "Pomerol: terms of Hamiltonian" << std::endl;
  std::cout << *storage << std::endl;
 }

 // Check the Hamiltonian commutes with the total number of particles
 Pomerol::OperatorPresets::N N(index_info.getIndexSize());
 if(!storage->commutes(N))
  TRIQS_RUNTIME_ERROR << "diagonalize: Hamiltonian does not conserve the total number of particles";

 // Construct Symmetrizer
 symm.reset(new Pomerol::Symmetrizer(index_info, *storage));
 symm->compute(ignore_symmetries);

 // Classify many-body states
 states_class.reset(new Pomerol::StatesClassification(index_info, *symm));
 states_class->compute();

 // Matrix representation of the Hamiltonian
 matrix_h.reset(new Pomerol::Hamiltonian(index_info, *storage, *states_class));
 matrix_h->prepare();
 matrix_h->compute(comm);

 // Get ground state energy
 if(verbose && !comm.rank()) {
  std::cout << "Pomerol: ground state energy is " << matrix_h->getGroundEnergy() + gs_shift << std::endl;
 }

 // Reset density matrix, we will compute it later if needed
 rho.release();
}

block_gf<imfreq> pomerol_ed::G_iw(gf_struct_t const& gf_struct, int n_iw) const {
 // TODO
}

block_gf<imtime> pomerol_ed::G_tau(gf_struct_t const& gf_struct, int n_tau) const {
 // TODO
}

}
