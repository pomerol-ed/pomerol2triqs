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
                              Filler filler) const -> block2_gf<Mesh, tensor_valued<4>> {

    if (!states_class || !matrix_h || !rho || !ops_container) TRIQS_RUNTIME_ERROR << "compute_chi3: Internal error!";

    bool compute_all_blocks = chi3_blocks.empty();

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

                  switch(channel) {
                    case PP: {
                      Pomerol::ParticleIndex pom_CX1_i = lookup_pomerol_index({A, a});
                      Pomerol::ParticleIndex pom_CX3_i = lookup_pomerol_index({B, c});
                      Pomerol::ParticleIndex pom_Delta_i1 = block_order == AABB ? lookup_pomerol_index({A, b}) :
                                                                                  lookup_pomerol_index({B, d});
                      Pomerol::ParticleIndex pom_Delta_i2 = block_order == AABB ? lookup_pomerol_index({B, d}) :
                                                                                  lookup_pomerol_index({A, b});

                      Pomerol::QuadraticOperator Delta(index_info, *hs, *states_class, *matrix_h, pom_Delta_i1, pom_Delta_i2,
                                                       {false, false});
                      Delta.prepare(*hs);
                      Delta.compute();

                      Pomerol::ThreePointSusceptibility pom_chi3(*states_class, *matrix_h,
                                                                 ops_container->getCreationOperator(pom_CX1_i),
                                                                 ops_container->getCreationOperator(pom_CX3_i),
                                                                 Delta, *rho);
                      pom_chi3.prepare();
                      pom_chi3.compute(false, {}, comm.get());

                      filler(chi3_el, pom_chi3);
                      break;
                    }
                    case PH: {
                      Pomerol::ParticleIndex pom_CX_i = lookup_pomerol_index({A, a});
                      Pomerol::ParticleIndex pom_C_i = block_order == AABB ? lookup_pomerol_index({A, b}) :
                                                                             lookup_pomerol_index({B, d});
                      Pomerol::ParticleIndex pom_N_i1 = lookup_pomerol_index({B, c});
                      Pomerol::ParticleIndex pom_N_i2 = block_order == AABB ? lookup_pomerol_index({B, d}) :
                                                                              lookup_pomerol_index({A, b});

                      Pomerol::QuadraticOperator N(index_info, *hs, *states_class, *matrix_h, pom_N_i1, pom_N_i2);
                      N.prepare(*hs);
                      N.compute();

                      Pomerol::ThreePointSusceptibility pom_chi3(*states_class, *matrix_h,
                                                                 ops_container->getCreationOperator(pom_CX_i),
                                                                 ops_container->getAnnihilationOperator(pom_C_i),
                                                                 N, *rho);
                      pom_chi3.prepare();
                      pom_chi3.compute(false, {}, comm.get());

                      filler(chi3_el, pom_chi3);
                      break;
                    }
                    case xPH: {
                      Pomerol::ParticleIndex pom_CX_i = lookup_pomerol_index({A, a});
                      Pomerol::ParticleIndex pom_C_i = block_order == AABB ? lookup_pomerol_index({B, d}) :
                                                                             lookup_pomerol_index({A, b});
                      Pomerol::ParticleIndex pom_N_i1 = lookup_pomerol_index({B, c});
                      Pomerol::ParticleIndex pom_N_i2 = block_order == AABB ? lookup_pomerol_index({A, b}) :
                                                                              lookup_pomerol_index({B, d});

                      Pomerol::QuadraticOperator N(index_info, *hs, *states_class, *matrix_h, pom_N_i1, pom_N_i2);
                      N.prepare(*hs);
                      N.compute();

                      Pomerol::ThreePointSusceptibility pom_chi3(*states_class, *matrix_h,
                                                                 ops_container->getCreationOperator(pom_CX_i),
                                                                 ops_container->getAnnihilationOperator(pom_C_i),
                                                                 N, *rho);
                      pom_chi3.prepare();
                      pom_chi3.compute(false, {}, comm.get());

                      filler(chi3_el, pom_chi3);
                      break;
                    }
                    default:
                      TRIQS_RUNTIME_ERROR << "compute_chi3: Internal error!";
                  }
                }
        }
      }

      gf_vecvec.emplace_back(std::move(gf_vec));
    }

    return make_block2_gf(block_names, block_names, std::move(gf_vecvec));
}

auto pomerol_ed::chi3_inu_inup(chi3_inu_inup_params_t const& p) -> block2_gf<nu_nup_t, tensor_valued<4>> {
  if(p.channel == AllFermionic) TRIQS_RUNTIME_ERROR << "chi3_inu_inup: AllFermionic channel is not supported";

  if (!matrix_h) TRIQS_RUNTIME_ERROR << "chi3_inu_inup: No Hamiltonian has been diagonalized";
  compute_rho(p.beta);
  compute_field_operators(p.gf_struct);

  if (verbose && !comm.rank()) std::cout << "chi3_inu_inup: Filling output container" << std::endl;

  auto filler = [&p, this](gf_view<nu_nup_t, scalar_valued> chi3_el, auto const &pom_chi3) {
    long mesh_index = 0;
    for (auto nu_nup : chi3_el.mesh()) {
      if ((mesh_index++) % comm.size() != comm.rank()) continue;

      auto w1 = std::complex<double>(std::get<0>(nu_nup));
      auto w2 = std::complex<double>(std::get<1>(nu_nup));

      if(p.channel == xPH)
        chi3_el[nu_nup] = -pom_chi3(w1, w2); // Extra minus sign from crossing-symmetry relation
      else
        chi3_el[nu_nup] = +pom_chi3(w1, w2);
    }
  };

  mesh::imfreq mesh_f{p.beta, Fermion, p.n_inu};
  nu_nup_t mesh_ff{mesh_f, mesh_f};

  auto chi3 = compute_chi3<nu_nup_t>(p.gf_struct, mesh_ff, p.block_order, p.channel, p.blocks, filler);

  chi3() = mpi::all_reduce(chi3(), comm);

  return chi3;
}

} // namespace pomerol2triqs
