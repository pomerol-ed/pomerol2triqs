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

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/constants/constants.hpp>

#include <nda/nda.hpp>

///////////////////////////////////////////////////////
// Methods to compute two-particle Green's functions //
///////////////////////////////////////////////////////

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

    std::complex<double> res = (sqrt(2 * l + 1) / sqrt(o)) * std::pow(dcomplex{1i}, o + l) * boost::math::cyl_bessel_j(l + 0.5, o * pi / 2);
    // \bar T_{-ol} = \bar T_{ol}^*
    return neg_o ? std::conj(res) : res;
  }

  template <typename Mesh, typename Filler>
  auto pomerol_ed::compute_g2(gf_struct_t const &gf_struct, Mesh const &mesh, block_order_t block_order, g2_blocks_t const &g2_blocks,
                              Filler filler) const -> block2_gf<Mesh, tensor_valued<4>> {

    if (!states_class || !matrix_h || !rho || !ops_container) TRIQS_RUNTIME_ERROR << "compute_g2: Internal error!";

    bool compute_all_blocks = g2_blocks.empty();

    std::vector<std::vector<gf<Mesh, tensor_valued<4>>>> gf_vecvec;
    std::vector<std::string> block_names;

    for (auto const &bl1 : gf_struct) {
      auto &A    = bl1.first;
      int A_size = bl1.second;
      int s1     = A_size;
      block_names.push_back(A);

      std::vector<gf<Mesh, tensor_valued<4>>> gf_vec;
      for (auto const &bl2 : gf_struct) {
        auto &B    = bl2.first;
        int B_size = bl2.second;
        int s3     = B_size;

        int s2 = block_order == AABB ? s1 : s3;
        int s4 = block_order == AABB ? s3 : s1;

        gf_vec.emplace_back(mesh, make_shape(s1, s2, s3, s4));

        if (compute_all_blocks || g2_blocks.count({A, B})) {
          auto &g2_block = gf_vec.back();

          for (int a : range(A_size))
            for (int b : range(A_size))
              for (int c : range(B_size))
                for (int d : range(B_size)) {

                  if (verbose && !comm.rank()) {
                    std::cout << "compute_g2: Filling G^2 element ";
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

                  Pomerol::ParticleIndex pom_i1 = lookup_pomerol_index({A, b});
                  Pomerol::ParticleIndex pom_i2 = lookup_pomerol_index({B, d});
                  Pomerol::ParticleIndex pom_i3 = lookup_pomerol_index({A, a});
                  Pomerol::ParticleIndex pom_i4 = lookup_pomerol_index({B, c});

                  Pomerol::TwoParticleGF pom_g2(*states_class, *matrix_h, ops_container->getAnnihilationOperator(pom_i1),
                                                ops_container->getAnnihilationOperator(pom_i2), ops_container->getCreationOperator(pom_i3),
                                                ops_container->getCreationOperator(pom_i4), *rho);
                  pom_g2.prepare();
                  pom_g2.compute(false, {}, comm.get());

                  filler(g2_el, pom_g2);
                }
        }
      }
      gf_vecvec.emplace_back(std::move(gf_vec));
    }

    return make_block2_gf(block_names, block_names, std::move(gf_vecvec));
  }

  auto pomerol_ed::G2_iw_inu_inup(g2_iw_inu_inup_params_t const &p) -> block2_gf<w_nu_nup_t, tensor_valued<4>> {
    if(p.channel == xPH) TRIQS_RUNTIME_ERROR << "G2_iw_inu_inup: Crossed particle-hole channel is not supported";

    if (!matrix_h) TRIQS_RUNTIME_ERROR << "G2_iw_inu_inup: No Hamiltonian has been diagonalized";
    compute_rho(p.beta);
    compute_field_operators(p.gf_struct);

    if (verbose && !comm.rank()) std::cout << "G2_iw_inu_inup: Filling output container" << std::endl;

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

    mesh::imfreq mesh_b{p.beta, Boson, p.n_iw};
    mesh::imfreq mesh_f{p.beta, Fermion, p.n_inu};

    w_nu_nup_t mesh_bff{mesh_b, mesh_f, mesh_f};
    w_nu_nup_t mesh_fff{mesh_f, mesh_f, mesh_f};

    block2_gf<w_nu_nup_t, tensor_valued<4>> g2;

    if (p.channel == AllFermionic)
      g2 = compute_g2<w_nu_nup_t>(p.gf_struct, mesh_fff, p.block_order, p.blocks, filler);
    else
      g2 = compute_g2<w_nu_nup_t>(p.gf_struct, mesh_bff, p.block_order, p.blocks, filler);

    g2() = mpi::all_reduce(g2(), comm);

    return g2;
  }

  auto pomerol_ed::G2_iw_l_lp(g2_iw_l_lp_params_t const &p) -> block2_gf<w_l_lp_t, tensor_valued<4>> {
    if(p.channel == xPH) TRIQS_RUNTIME_ERROR << "G2_iw_l_lp: Crossed particle-hole channel is not supported";

    if (!matrix_h) TRIQS_RUNTIME_ERROR << "G2_iw_l_lp: No Hamiltonian has been diagonalized";
    compute_rho(p.beta);
    compute_field_operators(p.gf_struct);

    w_l_lp_t mesh{{p.beta, Boson, p.n_iw}, {p.beta, Fermion, p.n_l}, {p.beta, Fermion, p.n_l}};

    if (verbose && !comm.rank()) std::cout << "G2_iw_l_lp: Filling output container" << std::endl;

    auto filler = [&p, this](gf_view<w_l_lp_t, scalar_valued> g2_el, auto const &pom_g2) {
      auto get_g2_iw_inu_inup_val = [&p, &pom_g2](long w_m, long nu_n, long nup_n) {
        int W_n = p.channel == PH ? w_m + nu_n : w_m - nup_n - 1;
        if (p.block_order == AABB)
          return -pom_g2(W_n, nup_n, nu_n);
        else
          return +pom_g2(nup_n, W_n, nu_n);
      };

      nda::array<std::complex<double>, 2> border_contrib(p.n_l, p.n_l);
      nda::array<bool, 2> llp_element_converged(p.n_l, p.n_l);
      int n_llp_elements_converged;

      long mesh_index = 0;
      for (auto iw : std::get<0>(g2_el.mesh())) {
        if ((mesh_index++) % comm.size() != comm.rank()) continue;

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

            for (auto l_mp : std::get<1>(g2_el.mesh()))
              for (auto lp_mp : std::get<2>(g2_el.mesh())) {
                auto l  = l_mp.index();
                auto lp = lp_mp.index();

                if (llp_element_converged(l, lp)) continue;

                using std::conj;
                std::complex<double> val = 0;
                val += t_bar(2 * r + w_m + 1, l) * iw_val_1 * conj(t_bar(2 * n + w_m + 1, lp));
                val += t_bar(2 * (-r - 1) + w_m + 1, l) * iw_val_2 * conj(t_bar(2 * (n - 1) + w_m + 1, lp));
                val += t_bar(2 * (n - 1) + w_m + 1, l) * iw_val_3 * conj(t_bar(2 * r + w_m + 1, lp));
                val += t_bar(2 * n + w_m + 1, l) * iw_val_4 * conj(t_bar(2 * (-r - 1) + w_m + 1, lp));

                g2_el[iw, l_mp, lp_mp] += val;
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
    g2()    = mpi::all_reduce(g2(), comm);

    return g2;
  }
} // namespace pomerol2triqs
