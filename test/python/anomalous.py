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

# Single orbital Anderson model with superconducting bath

####################
# Input parameters #
####################

beta = 10.0             # Inverse temperature
U = 2.0                 # Coulomb repulsion
mu = 1.0                # Chemical potential
Delta = 0.1             # Bath pairing amplitude

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
gf_struct = [('up_dn', 2)]

# Conversion from TRIQS to Pomerol notation for operator indices
index_converter = {('up_dn', 0): ("loc", 0, "up"),
                   ('up_dn', 1): ("loc", 0, "down")}
index_converter.update({("B%i_%s" % (k, sn), 0):
                        ("bath" + str(k), 0, "down" if sn == "dn" else "up")
                        for k, sn in product(range(len(epsilon)), spin_names)})

# Make PomerolED solver object
ed = PomerolED(index_converter, verbose = True)

# Number of particles on the impurity
H_loc = -mu * (n('up_dn', 0) + n('up_dn', 1)) + U * n('up_dn', 0) * n('up_dn', 1)

# Bath Hamiltonian
H_bath = sum(eps*n("B%i_%s" % (k, sn), 0)
             for sn, (k, eps) in product(spin_names, enumerate(epsilon)))

# Hybridization Hamiltonian
H_hyb = Operator()
for k, v in enumerate(V):
    H_hyb += sum(        v   * c_dag("B%i_%s" % (k, sn), 0) * c("up_dn", s) +
                 np.conj(v)  * c_dag("up_dn", s) * c("B%i_%s" % (k, sn), 0)
                 for s, sn in enumerate(spin_names))

# Pairing terms
H_sc = Delta * sum(c_dag("B%i_up" % k, 0) * c_dag("B%i_dn" % k, 0) +
                   c("B%i_dn" % k, 0) * c("B%i_up" % k, 0)
                   for k in range(len(epsilon)))

# Complete Hamiltonian
H = H_loc + H_hyb + H_bath + H_sc

# Diagonalize H
ed.diagonalize(H)

# Compute superconducting order parameter \phi = <c_\up c_\dn>
phi = ed.ensemble_average(("up_dn", 0), ("up_dn", 1), beta, (False, False))

#
# Compute normal GFs
#

# Compute G(i\omega)
G_iw = ed.G_iw(gf_struct, beta, n_iw)

# Compute G(\tau)
G_tau = ed.G_tau(gf_struct, beta, n_tau)

# Compute G(\omega)
G_w = ed.G_w(gf_struct, beta, energy_window, n_w, 0.01)

#
# Compute anomalous GFs
#

# Compute F(i\omega)
F_iw = ed.F_iw(gf_struct, beta, n_iw)

# Compute F(\tau)
F_tau = ed.F_tau(gf_struct, beta, n_tau)

# Compute F(\omega)
F_w = ed.F_w(gf_struct, beta, energy_window, n_w, 0.01)

# Check that F_{\up,\dn}(\tau=0) = -<c_\up c_\dn>
assert np.isclose(F_tau["up_dn"][0, 1].data[0], -phi, atol=1e-10)

if mpi.is_master_node():
    with HDFArchive('anomalous.out.h5', 'w') as ar:
        ar['H'] = H
        ar['phi'] = phi
        ar['G_iw'] = G_iw
        ar['G_tau'] = G_tau
        ar['G_w'] = G_w
        ar['F_iw'] = F_iw
        ar['F_tau'] = F_tau
        ar['F_w'] = F_w

    with HDFArchive("anomalous.ref.h5", 'r') as ar:
        assert (ar['H'] - H).is_zero()
        assert_block_gfs_are_close(ar['G_iw'], G_iw)
        assert_block_gfs_are_close(ar['G_tau'], G_tau)
        assert_block_gfs_are_close(ar['G_w'], G_w)
        assert_block_gfs_are_close(ar['F_iw'], F_iw)
        assert_block_gfs_are_close(ar['F_tau'], F_tau)
        assert_block_gfs_are_close(ar['F_w'], F_w)
