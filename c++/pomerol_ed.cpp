#include "pomerol_ed.hpp"

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/constants/constants.hpp>
#include <triqs/arrays.hpp>

namespace pomerol2triqs {

  // Generalization of triqs::utility::legendre_T
  // See eq. (4.5) in Lewin Boehnke's thesis
  inline std::complex<double> t_bar(int o, int l) {
    if (o == 0) return l == 0 ? 1 : 0;

    const double pi = boost::math::constants::pi<double>();
    bool neg_o      = false;
    if (o < 0) {
      neg_o = true;
      o     = -o;
    }

    std::complex<double> res = (sqrt(2 * l + 1) / sqrt(o)) * std::pow(1_j, o + l) * boost::math::cyl_bessel_j(l + 0.5, o * pi / 2);
    // \bar T_{-ol} = \bar T_{ol}^*
    return neg_o ? std::conj(res) : res;
  }

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

  //////////////////////
  // Green's function //
  //////////////////////

  template <typename Mesh, typename Filler>
  block_gf<Mesh> pomerol_ed::compute_gf(gf_struct_t const &gf_struct, gf_mesh<Mesh> const &mesh, Filler filler) const {

    if (!states_class || !matrix_h || !rho || !ops_container) TRIQS_RUNTIME_ERROR << "compute_gf: internal error!";

    struct index_visitor {
      std::vector<std::string> indices;
      void operator()(int i) { indices.push_back(std::to_string(i)); }
      void operator()(std::string s) { indices.push_back(s); }
    };

    std::vector<std::string> block_names;
    std::vector<gf<Mesh>> g_blocks;

    for (auto const &bl : gf_struct) {
      block_names.push_back(bl.first);
      int n = bl.second.size();

      index_visitor iv;
      for (auto &ind : bl.second) { apply_visitor(iv, ind); }
      std::vector<std::vector<std::string>> indices{{iv.indices, iv.indices}};

      g_blocks.push_back(gf<Mesh>{mesh, {n, n}, indices});
      auto &g = g_blocks.back();

      for (int i1 : range(n)) {
        Pomerol::ParticleIndex pom_i1 = lookup_pomerol_index({bl.first, bl.second[i1]});
        for (int i2 : range(n)) {
          Pomerol::ParticleIndex pom_i2 = lookup_pomerol_index({bl.first, bl.second[i2]});

          if (verbose && !comm.rank())
            std::cout << "fill_gf: Filling GF component (" << bl.first << "," << bl.second[i1] << ")(" << bl.first << "," << bl.second[i2] << ")"
                      << std::endl;
          auto g_el = slice_target_to_scalar(g, i1, i2);

          Pomerol::GreensFunction pom_g(*states_class, *matrix_h, ops_container->getAnnihilationOperator(pom_i1),
                                        ops_container->getCreationOperator(pom_i2), *rho);
          pom_g.prepare();
          pom_g.compute();

          filler(g_el, pom_g);
          g_el.singularity()(1) = (i1 == i2) ? 1 : 0;
        }
      }
    }
    return make_block_gf(block_names, std::move(g_blocks));
  }

  block_gf<imfreq> pomerol_ed::G_iw(gf_struct_t const &gf_struct, double beta, int n_iw) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "G_iw: no Hamiltonian has been diagonalized";
    compute_rho(beta);
    compute_field_operators(gf_struct);

    auto filler = [](gf_view<imfreq, scalar_valued> g_el, Pomerol::GreensFunction const &pom_g) {
      for (auto iw : g_el.mesh()) g_el[iw] = pom_g(std::complex<double>(iw));
    };
    return compute_gf<imfreq>(gf_struct, {beta, Fermion, n_iw}, filler);
  }

  block_gf<imtime> pomerol_ed::G_tau(gf_struct_t const &gf_struct, double beta, int n_tau) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "G_tau: no Hamiltonian has been diagonalized";
    compute_rho(beta);
    compute_field_operators(gf_struct);

    auto filler = [](gf_view<imtime, scalar_valued> g_el, Pomerol::GreensFunction const &pom_g) {
      for (auto tau : g_el.mesh()) g_el[tau] = pom_g.of_tau(tau);
    };
    return compute_gf<imtime>(gf_struct, {beta, Fermion, n_tau}, filler);
  }

  block_gf<refreq> pomerol_ed::G_w(gf_struct_t const &gf_struct, double beta, std::pair<double, double> const &energy_window, int n_w,
                                   double im_shift) {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "G_w: no Hamiltonian has been diagonalized";
    compute_rho(beta);
    compute_field_operators(gf_struct);

    auto filler = [im_shift](gf_view<refreq, scalar_valued> g_el, Pomerol::GreensFunction const &pom_g) {
      for (auto w : g_el.mesh()) g_el[w] = pom_g(double(w) + 1_j * im_shift);
    };
    return compute_gf<refreq>(gf_struct, {energy_window.first, energy_window.second, n_w}, filler);
  }

  ///////////////////////////////////
  // Two-particle Green's function //
  ///////////////////////////////////

  template <typename Mesh, typename Filler>
  auto pomerol_ed::compute_g2(gf_struct_t const &gf_struct, gf_mesh<Mesh> const &mesh, block_order_t block_order, g2_blocks_t const &g2_blocks,
                              Filler filler) const -> block2_gf<Mesh, tensor_valued<4>> {

    if (!states_class || !matrix_h || !rho || !ops_container) TRIQS_RUNTIME_ERROR << "compute_g2: internal error!";

    bool compute_all_blocks = g2_blocks.empty();

    std::vector<std::vector<gf<Mesh, tensor_valued<4>>>> gf_vecvec;
    std::vector<std::string> block_names;

    for (auto const &bl1 : gf_struct) {
      auto &A    = bl1.first;
      int A_size = bl1.second.size();
      int s1     = A_size;
      block_names.push_back(A);

      std::vector<gf<Mesh, tensor_valued<4>>> gf_vec;
      for (auto const &bl2 : gf_struct) {
        auto &B    = bl2.first;
        int B_size = bl2.second.size();
        int s3     = B_size;

        int s2 = block_order == AABB ? s1 : s3;
        int s4 = block_order == AABB ? s3 : s1;

        gf_vec.emplace_back(mesh, make_shape(s1, s2, s3, s4));

        if (compute_all_blocks || g2_blocks.count({A, B})) {
          auto const &A_inner = bl1.second;
          auto const &B_inner = bl2.second;

          auto &g2_block = gf_vec.back();

          for (int a : range(A_size))
            for (int b : range(A_size))
              for (int c : range(B_size))
                for (int d : range(B_size)) {

                  if (verbose && !comm.rank()) {
                    std::cout << "compute_g2: filling G^2 element ";
                    if (block_order == AABB) {
                      std::cout << "(" << A << "," << a << ")";
                      std::cout << "(" << A << "," << b << ")";
                      std::cout << "(" << B << "," << c << ")";
                      std::cout << "(" << B << "," << d << ")";
                    } else {
                      std::cout << "(" << A << "," << a << ")";
                      std::cout << "(" << B << "," << d << ")";
                      std::cout << "(" << B << "," << c << ")";
                      std::cout << "(" << A << "," << b << ")";
                    }
                    std::cout << std::endl;
                  }

                  auto g2_el = block_order == AABB ? slice_target_to_scalar(g2_block, a, b, c, d) : slice_target_to_scalar(g2_block, a, d, c, b);

                  Pomerol::ParticleIndex pom_i1 = lookup_pomerol_index({A, A_inner[b]});
                  Pomerol::ParticleIndex pom_i2 = lookup_pomerol_index({B, B_inner[d]});
                  Pomerol::ParticleIndex pom_i3 = lookup_pomerol_index({A, A_inner[a]});
                  Pomerol::ParticleIndex pom_i4 = lookup_pomerol_index({B, B_inner[c]});

                  Pomerol::TwoParticleGF pom_g2(*states_class, *matrix_h, ops_container->getAnnihilationOperator(pom_i1),
                                                ops_container->getAnnihilationOperator(pom_i2), ops_container->getCreationOperator(pom_i3),
                                                ops_container->getCreationOperator(pom_i4), *rho);
                  pom_g2.prepare();
                  pom_g2.compute(false, {}, comm);

                  filler(g2_el, pom_g2);
                }
        }
      }
      gf_vecvec.emplace_back(std::move(gf_vec));
    }

    return make_block2_gf(block_names, block_names, std::move(gf_vecvec));
  }

  auto pomerol_ed::G2_iw_inu_inup(g2_iw_inu_inup_params_t const &p) -> block2_gf<w_nu_nup_t, tensor_valued<4>> {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "G2_iw_inu_inup: no Hamiltonian has been diagonalized";
    compute_rho(p.beta);
    compute_field_operators(p.gf_struct);

    if (verbose && !comm.rank()) std::cout << "G2_iw_inu_inup: filling output container" << std::endl;

    auto filler = [&p, this](gf_view<w_nu_nup_t, scalar_valued> g2_el, auto const &pom_g2) {
      long mesh_index = 0;
      for (auto w_nu_nup : g2_el.mesh()) {
        if ((mesh_index++) % comm.size() != comm.rank()) continue;

        if (p.channel == AllFermionic) {

          int n1 = std::get<0>(w_nu_nup).index();
          int n2 = std::get<1>(w_nu_nup).index();
          int n3 = std::get<2>(w_nu_nup).index();

          if (p.block_order == AABB)
            g2_el[w_nu_nup] = -pom_g2(n2, n1 + n3 - n2, n1);
          else
            g2_el[w_nu_nup] = +pom_g2(n1 + n3 - n2, n2, n1);

        } else { // p.channel == PH or PP

          int w_n   = std::get<0>(w_nu_nup).index();
          int nu_n  = std::get<1>(w_nu_nup).index();
          int nup_n = std::get<2>(w_nu_nup).index();

          int W_n = p.channel == PH ? w_n + nu_n : w_n - nup_n - 1;

          if (p.block_order == AABB) {
            g2_el[w_nu_nup] = -pom_g2(W_n, nup_n, nu_n);
          } else {
            g2_el[w_nu_nup] = +pom_g2(nup_n, W_n, nu_n);
          }
        }
      }
    };

    gf_mesh<imfreq> mesh_b{p.beta, Boson, p.n_iw};
    gf_mesh<imfreq> mesh_f{p.beta, Fermion, p.n_inu};

    gf_mesh<w_nu_nup_t> mesh_bff{mesh_b, mesh_f, mesh_f};
    gf_mesh<w_nu_nup_t> mesh_fff{mesh_f, mesh_f, mesh_f};

    block2_gf<w_nu_nup_t, tensor_valued<4>> g2;

    if (p.channel == AllFermionic)
      g2 = compute_g2<w_nu_nup_t>(p.gf_struct, mesh_fff, p.block_order, p.blocks, filler);
    else
      g2 = compute_g2<w_nu_nup_t>(p.gf_struct, mesh_bff, p.block_order, p.blocks, filler);

    g2() = mpi_all_reduce(g2(), comm);

    return g2;
  }

  auto pomerol_ed::G2_iw_l_lp(g2_iw_l_lp_params_t const &p) -> block2_gf<w_l_lp_t, tensor_valued<4>> {
    if (!matrix_h) TRIQS_RUNTIME_ERROR << "G2_iw_l_lp: no Hamiltonian has been diagonalized";
    compute_rho(p.beta);
    compute_field_operators(p.gf_struct);

    gf_mesh<w_l_lp_t> mesh{{p.beta, Boson, p.n_iw}, {p.beta, Fermion, static_cast<size_t>(p.n_l)}, {p.beta, Fermion, static_cast<size_t>(p.n_l)}};

    if (verbose && !comm.rank()) std::cout << "G2_iw_l_lp: filling output container" << std::endl;

    auto filler = [&p, this](gf_view<w_l_lp_t, scalar_valued> g2_el, auto const &pom_g2) {

      auto get_g2_iw_inu_inup_val = [&p, &pom_g2](long w_m, long nu_n, long nup_n) {
        int W_n = p.channel == PH ? w_m + nu_n : w_m - nup_n - 1;
        if (p.block_order == AABB)
          return -pom_g2(W_n, nup_n, nu_n);
        else
          return +pom_g2(nup_n, W_n, nu_n);
      };

      array<std::complex<double>, 2> border_contrib(p.n_l, p.n_l);
      array<bool, 2> llp_element_converged(p.n_l, p.n_l);
      int n_llp_elements_converged;

      long mesh_index = 0;
      for (auto iw : std::get<0>(g2_el.mesh())) {
        if((mesh_index++) % comm.size() != comm.rank()) continue;

        int w_m = iw.index();

        llp_element_converged()  = false;
        n_llp_elements_converged = 0;

        // Summation over n and n' is done by adding new border points to
        // a square summation domain.
        for (int r = 0; r < p.n_inu_sum; ++r) { // r is current size of the domain
          if (n_llp_elements_converged == p.n_l * p.n_l) break;

          border_contrib() = 0;
          for (int n = -r; n <= r; ++n) {
            auto iw_val_1 = get_g2_iw_inu_inup_val(w_m, r, n);
            auto iw_val_2 = get_g2_iw_inu_inup_val(w_m, -r - 1, n - 1);
            auto iw_val_3 = get_g2_iw_inu_inup_val(w_m, n - 1, r);
            auto iw_val_4 = get_g2_iw_inu_inup_val(w_m, n, -r - 1);

            for (auto l : std::get<1>(g2_el.mesh()))
              for (auto lp : std::get<2>(g2_el.mesh())) {
                if (llp_element_converged(l, lp)) continue;

                using std::conj;
                std::complex<double> val = 0;
                val += t_bar(2 * r + w_m + 1, l) * iw_val_1 * conj(t_bar(2 * n + w_m + 1, lp));
                val += t_bar(2 * (-r - 1) + w_m + 1, l) * iw_val_2 * conj(t_bar(2 * (n - 1) + w_m + 1, lp));
                val += t_bar(2 * (n - 1) + w_m + 1, l) * iw_val_3 * conj(t_bar(2 * r + w_m + 1, lp));
                val += t_bar(2 * n + w_m + 1, l) * iw_val_4 * conj(t_bar(2 * (-r - 1) + w_m + 1, lp));

                g2_el[{iw, l, lp}] += val;
                border_contrib(l, lp) += val;
                if (std::abs(border_contrib(l, lp)) < p.inu_sum_tol) {
                  llp_element_converged(l, lp) = true;
                  ++n_llp_elements_converged;
                }
              }
          }
        }
      }
    };

    auto g2 = compute_g2<w_l_lp_t>(p.gf_struct, mesh, p.block_order, p.blocks, filler);
    g2() = mpi_all_reduce(g2(), comm);

    return g2;
  }
}
