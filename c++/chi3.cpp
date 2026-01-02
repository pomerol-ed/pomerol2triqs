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

#include <array>
#include <utility>

/////////////////////////////////////////////////////////////
// Methods to compute 3-point fermion-boson susceptibility //
/////////////////////////////////////////////////////////////

namespace pomerol2triqs {

template <typename Mesh, typename Filler>
auto pomerol_ed::compute_chi3(gf_struct_t const &gf_struct,
                              Mesh const &mesh,
                              block_order_t block_order,
                              channel_t channel,
                              chi3_blocks_t const &chi3_blocks,
                              Filler filler,
                              double pole_res, double coeff_tol) const -> block2_gf<Mesh, tensor_valued<4>> {

    if (!states_class || !matrix_h || !rho || !ops_container) TRIQS_RUNTIME_ERROR << "compute_chi3: Internal error!";

    bool compute_all_blocks = chi3_blocks.empty();

    std::vector<std::vector<gf<Mesh, tensor_valued<4>>>> gf_vecvec;
    std::vector<std::string> block_names;

    static const std::array<Pomerol::Channel, 3> pom_channels = {Pomerol::PP, Pomerol::PH, Pomerol::xPH};

    for (auto const &[A, A_size] : gf_struct) {
      block_names.push_back(A);

      std::vector<gf<Mesh, tensor_valued<4>>> gf_vec;
      for (auto const &[B, B_size] : gf_struct) {

        int s2 = block_order == AABB ? A_size : B_size;
        int s4 = block_order == AABB ? B_size : A_size;

        gf_vec.emplace_back(mesh, make_shape(A_size, s2, B_size, s4));

        if (compute_all_blocks || chi3_blocks.count({A, B})) {
          auto &chi3_block = gf_vec.back();

          for (int a : range(A_size))
            for (int b : range(A_size))
              for (int c : range(B_size))
                for (int d : range(B_size)) {

                  if (verbose && !comm.rank()) {
                    std::cout << "compute_chi3: Filling \\chi^3 element ";
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

                  auto chi3_el = block_order == AABB ? slice_target_to_scalar(chi3_block, a, b, c, d) :
                                                       slice_target_to_scalar(chi3_block, a, d, c, b);

                  Pomerol::ParticleIndex pom_CX1_i = lookup_pomerol_index({A, a});
                  Pomerol::ParticleIndex pom_C2_i = block_order == AABB ? lookup_pomerol_index({A, b}) :
                                                                          lookup_pomerol_index({B, d});
                  Pomerol::ParticleIndex pom_CX3_i = lookup_pomerol_index({B, c});
                  Pomerol::ParticleIndex pom_C4_i = block_order == AABB ? lookup_pomerol_index({B, d}) :
                                                                          lookup_pomerol_index({A, b});

                  Pomerol::Channel pom_channel = pom_channels[int(channel)];

                  Pomerol::ThreePointSusceptibility pom_chi3(pom_channel,
                                                             *states_class,
                                                             *matrix_h,
                                                             ops_container->getCreationOperator(pom_CX1_i),
                                                             ops_container->getAnnihilationOperator(pom_C2_i),
                                                             ops_container->getCreationOperator(pom_CX3_i),
                                                             ops_container->getAnnihilationOperator(pom_C4_i),
                                                             *rho);

                  pom_chi3.PoleResolution = pole_res;
                  pom_chi3.CoefficientTolerance = coeff_tol;
                  pom_chi3.prepare();
                  pom_chi3.compute(false, {}, comm.get());

                  filler(chi3_el, pom_chi3);
                }
        }
      }

      gf_vecvec.emplace_back(std::move(gf_vec));
    }

    return make_block2_gf(block_names, block_names, std::move(gf_vecvec));
}

auto pomerol_ed::chi3_iw_inu(chi3_iw_inu_params_t const& p) -> block2_gf<w_nu_t, tensor_valued<4>> {
  if(p.channel == AllFermionic) TRIQS_RUNTIME_ERROR << "chi3_iw_inu: AllFermionic channel is not supported";

  if (!matrix_h) TRIQS_RUNTIME_ERROR << "chi3_iw_inu: No Hamiltonian has been diagonalized";
  compute_rho(p.beta);
  compute_field_operators(p.gf_struct);

  if (verbose && !comm.rank()) std::cout << "chi3_iw_inu: Filling output container" << std::endl;

  auto filler = [&p, this](gf_view<w_nu_t, scalar_valued> chi3_el, auto const &pom_chi3) {
    long mesh_index = 0;
    for (auto [w, nu] : chi3_el.mesh()) {
      if ((mesh_index++) % comm.size() != comm.rank()) continue;

      if(p.channel == PP)
        chi3_el[w, nu] = pom_chi3(nu, w - nu);
      else // PH and xPH
        chi3_el[w, nu] = pom_chi3(nu, w + nu);
    }
  };

  mesh::imfreq mesh_b{p.beta, Boson, p.n_iw};
  mesh::imfreq mesh_f{p.beta, Fermion, p.n_inu};
  w_nu_t mesh_bf{mesh_b, mesh_f};

  auto chi3 = compute_chi3<w_nu_t>(p.gf_struct, mesh_bf, p.block_order, p.channel, p.blocks, filler,
                                   p.pole_res, p.coeff_tol);

  chi3() = mpi::all_reduce(chi3(), comm);

  return chi3;
}

} // namespace pomerol2triqs
