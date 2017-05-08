#pragma once

#include <triqs/hilbert_space/fundamental_operator_set.hpp>
using triqs::hilbert_space::gf_struct_t;

namespace pomerol2triqs {

  enum block_order_t { AABB, ABBA };
  enum channel_t { PP, PH };

   using g2_blocks_t = std::set<std::pair<std::string, std::string>>;

  struct g2_parameters_t {

    /// Structure of G^2 blocks.
    gf_struct_t gf_struct;

    /// Inverse temperature
    double beta;

    /// Channel in which Matsubara frequency representation is defined.
    channel_t channel = PH;

    /// Order of block indices in the definition of G^2.
    block_order_t block_order = AABB;

    /// List of block index pairs of G^2 to measure.
    /// default: measure all blocks
    g2_blocks_t blocks = g2_blocks_t{};

    /// Number of bosonic Matsubara frequencies.
    int n_iw = 30;

    /// Number of fermionic Matsubara frequencies.
    int n_inu = 30;

    g2_parameters_t() {}
    g2_parameters_t(gf_struct_t const& gf_struct, double beta) : gf_struct(gf_struct), beta(beta) {}
  };
}
