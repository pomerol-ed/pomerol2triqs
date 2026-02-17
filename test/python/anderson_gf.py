################################################################################
#
# pomerol2triqs
#
# Copyright (C) 2017-2026 Igor Krivenko
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

from h5 import HDFArchive
from triqs.gf import *
from triqs.operators import Operator, c, c_dag, n
from triqs.utility import mpi
from triqs.utility.comparison_tests import *
from pomerol2triqs import PomerolED
import numpy as np
from itertools import product

# Single orbital Anderson model

####################
# Input parameters #
####################

beta = 10.0             # Inverse temperature
U = 2.0                 # Coulomb repulsion
mu = 1.0                # Chemical potential

# Levels of the bath sites
epsilon = [-1.0, 1.0]
# Hopping amplitudes
V = [0.5, 0.5]

spin_names = ("up", "dn")

# Number of Matsubara frequencies for GF calculation
n_iw = 200

# Number of imaginary time slices for GF calculation
n_tau = 1001

# Energy window for real frequency GF calculation
energy_window = (-5, 5)
# Number of frequency points for real frequency GF calculation
n_w = 1000

# GF structure
gf_struct = [['up', 1], ['dn', 1]]

# Conversion from TRIQS to Pomerol notation for operator indices
index_converter = {}
index_converter.update({(sn, 0) : ("loc", 0, "down" if sn == "dn" else "up") for sn in spin_names})
index_converter.update({("B%i_%s" % (k, sn), 0) : ("bath" + str(k), 0, "down" if sn == "dn" else "up")
                        for k, sn in product(list(range(len(epsilon))), spin_names)})

# Make PomerolED solver object
ed = PomerolED(index_converter, verbose = True)

# Number of particles on the impurity
H_loc = -mu*(n('up', 0) + n('dn', 0)) + U * n('up', 0) * n('dn', 0)

# Bath Hamiltonian
H_bath = sum(eps*n("B%i_%s" % (k, sn), 0)
             for sn, (k, eps) in product(spin_names, enumerate(epsilon)))

# Hybridization Hamiltonian
H_hyb = Operator()
for k, v in enumerate(V):
    H_hyb += sum(        v   * c_dag("B%i_%s" % (k, sn), 0) * c(sn, 0) +
                 np.conj(v)  * c_dag(sn, 0) * c("B%i_%s" % (k, sn), 0)
                 for sn in spin_names)

# Complete Hamiltonian
H = H_loc + H_hyb + H_bath

# Diagonalize H
ed.diagonalize(H)

# Structure of the Hilbert space
assert ed.full_hilbert_space_dim == 64
assert ed.n_subspaces == 16
fock_states_ref = frozenset([
    frozenset([0]),                                     # N_dn = 0, N_up = 0
    frozenset([1, 4, 16]),                              # N_dn = 1, N_up = 0
    frozenset([2, 8, 32]),                              # N_dn = 0, N_up = 1
    frozenset([5, 17, 20]),                             # N_dn = 2, N_up = 0
    frozenset([10, 34, 40]),                            # N_dn = 0, N_up = 2
    frozenset([3, 6, 9, 12, 18, 24, 33, 36, 48]),       # N_dn = 1, N_up = 1
    frozenset([21]),                                    # N_dn = 3, N_up = 0
    frozenset([42]),                                    # N_dn = 0, N_up = 3
    frozenset([7, 13, 19, 22, 25, 28, 37, 49, 52]),     # N_dn = 2, N_up = 1
    frozenset([11, 14, 26, 35, 38, 41, 44, 50, 56]),    # N_dn = 1, N_up = 2
    frozenset([23, 29, 53]),                            # N_dn = 3, N_up = 1
    frozenset([43, 46, 58]),                            # N_dn = 1, N_up = 3
    frozenset([15, 27, 30, 39, 45, 51, 54, 57, 60]),    # N_dn = 2, N_up = 2
    frozenset([47, 59, 62]),                            # N_dn = 2, N_up = 3
    frozenset([31, 55, 61]),                            # N_dn = 3, N_up = 2
    frozenset([63])                                     # N_dn = 3, N_up = 3
])

assert sorted(ed.subspace_dims) == sorted(map(len, fock_states_ref))
assert sorted(ed.subspace_dim(sp) for sp in range(16)) == sorted(map(len, fock_states_ref))
assert frozenset(frozenset(s) for s in ed.fock_states) == fock_states_ref
assert frozenset(frozenset(ed.subspace_fock_states(sp)) for sp in range(16)) == fock_states_ref

# Compute G(i\omega)
G_iw = ed.G_iw(gf_struct, beta, n_iw)

# Compute G(\tau)
G_tau = ed.G_tau(gf_struct, beta, n_tau)

# Compute G(\omega)
G_w = ed.G_w(gf_struct, beta, energy_window, n_w, 0.01)

if mpi.is_master_node():
    with HDFArchive('anderson_gf.out.h5', 'w') as ar:
        ar['H'] = H
        ar['G_iw'] = G_iw
        ar['G_tau'] = G_tau
        ar['G_w'] = G_w

    with HDFArchive("anderson_gf.ref.h5", 'r') as ar:
        assert (ar['H'] - H).is_zero()
        assert_block_gfs_are_close(ar['G_iw'], G_iw)
        assert_block_gfs_are_close(ar['G_tau'], G_tau)
        assert_block_gfs_are_close(ar['G_w'], G_w)
