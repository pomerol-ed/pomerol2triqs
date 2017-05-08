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

 // Reset containers, we will compute them later if needed
 rho.release();
 ops_container.release();
 gf_container.release();
 g2_container.release();
}

Pomerol::ParticleIndex pomerol_ed::lookup_pomerol_index(indices_t const& i) const {
 auto it = index_converter.find(i);
 if(it == index_converter.end()) return -1;
 std::string site;
 unsigned short orb;
 Pomerol::spin s;
 std::tie(site, orb, s) = it->second;
 return index_info.getIndex(site, orb, s);
}

// Translate gf_struct into a set of ParticleIndex
std::set<Pomerol::ParticleIndex> pomerol_ed::gf_struct_to_pomerol_indices(gf_struct_t const& gf_struct) const {
 std::set<Pomerol::ParticleIndex> indices;
 for(auto const& b : gf_struct) {
  for(auto const& i : b.second) {
   auto pom_ind = lookup_pomerol_index({b.first, i});
   if(pom_ind == -1)
    TRIQS_RUNTIME_ERROR << "gf_struct_to_pomerol_indices: unexpected GF index " << b.first << "," << i;
   indices.insert(pom_ind);
  }
 }
 return indices;
}

// Create the Density Matrix.
void pomerol_ed::compute_rho(double beta) {
 if(!states_class || !matrix_h)
  TRIQS_RUNTIME_ERROR << "compute_rho: internal error!";

 if(!rho || rho->beta != beta) {
  if(verbose && !comm.rank())
   std::cout << "Pomerol: computing density matrix for \\beta = " << beta << std::endl;
  rho.reset(new Pomerol::DensityMatrix(*states_class, *matrix_h, beta));
  rho->prepare();
  rho->compute();
 }
}

void pomerol_ed::compute_field_operators(gf_struct_t const& gf_struct) {
 if(!states_class || !matrix_h)
  TRIQS_RUNTIME_ERROR << "compute_field_operators: internal error!";

 auto new_ops = gf_struct_to_pomerol_indices(gf_struct);
 if(!ops_container || computed_ops != new_ops) {
  if(verbose && !comm.rank()) {
   std::cout << "Pomerol: computing field operators with indices ";
   bool comma = false;
   for(auto i : new_ops) {
    std::cout << (comma ? ", " : "") << i;
    comma = true;
   }
   std::cout << std::endl;
  }

  ops_container.reset(new Pomerol::FieldOperatorContainer(index_info,
                                                          *states_class,
                                                          *matrix_h));
  ops_container->prepareAll(new_ops);
  ops_container->computeAll();

  computed_ops = new_ops;
  gf_container.release();
  g2_container.release();
 }
}

void pomerol_ed::compute_gfs() {
 if(!states_class || !matrix_h || !rho || !ops_container)
  TRIQS_RUNTIME_ERROR << "compute_gfs: internal error!";

 if(!gf_container) {
  gf_container.reset(new Pomerol::GFContainer(index_info,
                                              *states_class,
                                              *matrix_h,
                                              *rho,
                                              *ops_container));

  std::set<Pomerol::IndexCombination2> in;
  for(auto i1 : computed_ops)
  for(auto i2 : computed_ops)
   in.insert({i1,i2});

  gf_container->prepareAll(in);
  gf_container->computeAll();
 }
}

void pomerol_ed::compute_g2(gf_struct_t const& gf_struct, g2_blocks_t g2_blocks) {
 if(!states_class || !matrix_h || !rho || !ops_container)
  TRIQS_RUNTIME_ERROR << "compute_g2: internal error!";

 if(g2_blocks.empty()) {
  for(auto const& bl1 : gf_struct)
  for(auto const& bl2 : gf_struct)
   g2_blocks.insert({bl1.first, bl2.first});
 }

 if(!g2_container || computed_g2_blocks != g2_blocks) {
  g2_container.reset(new Pomerol::TwoParticleGFContainer(index_info,
                                                         *states_class,
                                                         *matrix_h,
                                                         *rho,
                                                         *ops_container));

  std::set<Pomerol::IndexCombination4> in;

  auto get_inner_indices = [&gf_struct](std::string const& bn) -> gf_struct_t:: mapped_type const& {
   auto it = gf_struct.find(bn);
   if(it == gf_struct.end())
    TRIQS_RUNTIME_ERROR << "compute_g2: block " << bn << " is not in gf_struct";
   return it->second;
  };

  for(auto const& bl : g2_blocks) {
   std::string const& A = bl.first;
   std::string const& B = bl.second;
   auto const& A_inner = get_inner_indices(A);
   auto const& B_inner = get_inner_indices(B);

   for(auto a : A_inner)
   for(auto b : A_inner)
   for(auto c : B_inner)
   for(auto d : B_inner) {
    in.insert({lookup_pomerol_index({A, a}),
               lookup_pomerol_index({A, b}),
               lookup_pomerol_index({B, c}),
               lookup_pomerol_index({B, d})});
   }
  }

  g2_container->prepareAll(in);
  g2_container->computeAll(false, {}, comm);

  computed_g2_blocks = g2_blocks;
 }
}

template<typename Mesh, typename Filler>
block_gf<Mesh> pomerol_ed::fill_gf(gf_struct_t const& gf_struct, gf_mesh<Mesh> const& mesh, Filler filler) const {

 struct index_visitor  {
  std::vector<std::string> indices;
  void operator()(int i) { indices.push_back(std::to_string(i)); }
  void operator()(std::string s) { indices.push_back(s); }
 };

 std::vector<std::string> block_names;
 std::vector<gf<Mesh>> g_blocks;

 for (auto const& bl : gf_struct) {
  block_names.push_back(bl.first);
  int n = bl.second.size();

  index_visitor iv;
  for (auto & ind: bl.second) { apply_visitor(iv, ind); }
  std::vector<std::vector<std::string>> indices{{iv.indices,iv.indices}};

  g_blocks.push_back(gf<Mesh>{mesh, {n, n}, indices});
  auto & g = g_blocks.back();

  for(int i1 : range(n)) {
   Pomerol::ParticleIndex pom_i1 = lookup_pomerol_index({bl.first, bl.second[i1]});
   for(int i2 : range(n)) {
    Pomerol::ParticleIndex pom_i2 = lookup_pomerol_index({bl.first, bl.second[i2]});

    if(verbose && !comm.rank())
     std::cout << "fill_gf: Filling GF component (" << bl.first << ","
                                                    << bl.second[i1] << ")("
                                                    << bl.first << ","
                                                    << bl.second[i2] << ")"
                                                    << std::endl;
     auto g_el = slice_target_to_scalar(g, i1, i2);
     auto const& pom_g = (*gf_container)({pom_i1, pom_i2});
     filler(g_el, pom_g);
     g_el.singularity()(1) = (i1 == i2) ? 1 : 0;
   }
  }
 }
 return make_block_gf(block_names, std::move(g_blocks));
}

block_gf<imfreq> pomerol_ed::G_iw(gf_struct_t const& gf_struct, double beta, int n_iw) {
 if(!matrix_h) TRIQS_RUNTIME_ERROR << "G_iw: no Hamiltonian has been diagonalized";
 compute_rho(beta);
 compute_field_operators(gf_struct);
 compute_gfs();

 auto filler = [](gf_view<imfreq, scalar_valued> g_el, Pomerol::GreensFunction const& pom_g) {
  for(auto iw : g_el.mesh()) g_el[iw] = pom_g(std::complex<double>(iw));
 };
 return fill_gf<imfreq>(gf_struct, {beta, Fermion, n_iw}, filler);
}

block_gf<imtime> pomerol_ed::G_tau(gf_struct_t const& gf_struct, double beta, int n_tau) {
 if(!matrix_h) TRIQS_RUNTIME_ERROR << "G_tau: no Hamiltonian has been diagonalized";
 compute_rho(beta);
 compute_field_operators(gf_struct);
 compute_gfs();

 auto filler = [](gf_view<imtime, scalar_valued> g_el, Pomerol::GreensFunction const& pom_g) {
  for(auto tau : g_el.mesh()) g_el[tau] = pom_g.of_tau(tau);
 };
 return fill_gf<imtime>(gf_struct, {beta, Fermion, n_tau}, filler);
}

block_gf<refreq> pomerol_ed::G_w(gf_struct_t const& gf_struct,
                                 double beta,
                                 std::pair<double, double> const& energy_window,
                                 int n_w,
                                 double im_shift) {
 if(!matrix_h) TRIQS_RUNTIME_ERROR << "G_w: no Hamiltonian has been diagonalized";
 compute_rho(beta);
 compute_field_operators(gf_struct);
 compute_gfs();

 auto filler = [im_shift](gf_view<refreq, scalar_valued> g_el, Pomerol::GreensFunction const& pom_g) {
  for(auto w : g_el.mesh()) g_el[w] = pom_g(double(w) + 1_j*im_shift);
 };
 return fill_gf<refreq>(gf_struct, {energy_window.first, energy_window.second, n_w}, filler);
}

block2_gf<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>> pomerol_ed::G2_iw(g2_parameters_t const& p) {
 if(!matrix_h) TRIQS_RUNTIME_ERROR << "G2_iw: no Hamiltonian has been diagonalized";
 compute_rho(p.beta);
 compute_field_operators(p.gf_struct);
 compute_g2(p.gf_struct, p.blocks);

 std::vector<std::vector<gf<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>>>> gf_vecvec;
 std::vector<std::string> block_names;

 gf_mesh<cartesian_product<imfreq, imfreq, imfreq>> mesh{{p.beta, Fermion, p.n_iw},
                                                         {p.beta, Fermion, p.n_inu},
                                                         {p.beta, Fermion, p.n_inu}};

 for (auto const& bl1 : p.gf_struct) {
  auto & A  = bl1.first;
  int bl1_size = bl1.second.size();
  block_names.push_back(A);

  std::vector<gf<cartesian_product<imfreq, imfreq, imfreq>, tensor_valued<4>>> gf_vec;
  for (auto const& bl2 : p.gf_struct) {
   auto & B  = bl2.first;
   int bl2_size = bl2.second.size();

   gf_vec.emplace_back(mesh, make_shape(bl1_size, bl1_size, bl2_size, bl2_size));

   if(computed_g2_blocks.count({A, B})) {
    auto const& A_inner = bl1.second;
    auto const& B_inner = bl2.second;

    auto & g2_block = gf_vec.back();

    // TODO
   }
  }
  gf_vecvec.emplace_back(std::move(gf_vec));
 }

 return make_block2_gf(block_names, block_names, std::move(gf_vecvec));
}

}
