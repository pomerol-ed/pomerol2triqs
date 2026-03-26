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
from pomerol2triqs import PomerolED, BosonParams
import numpy as np
from itertools import product

# Single orbital Anderson model + bosons

####################
# Input parameters #
####################

beta = 10.0             # Inverse temperature
U = 2.0                 # Coulomb repulsion
mu = 1.0                # Chemical potential
w0_ch = 0.5             # Frequency of the phonon coupled to the charge
w0_sp = 0.4             # Frequency of the phonon coupled to the spin
g_ch = 0.2              # Strength of the charge-phonon coupling
g_sp = 0.1              # Strength of the spin-phonon coupling
n_bits = 2              # Consider occupations of bosonic modes 0 ... 2^2 - 1

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

# Bosonic contributions to the Hamiltonian
#
# Each bosonic mode contributes the terms w0 a^\dagger a + \hat O a^\dagger + \hat O^\dagger a,
# where \hat O is a fermionic coupling operator
bosons = [
    # Boson coupled to charge
    BosonParams(frequency=w0_ch, coupling=g_ch * (n('up', 0) + n('dn', 0)), n_bits=n_bits),
    # Boson coupled to spin
    BosonParams(frequency=w0_sp, coupling=g_sp * (n('up', 0) - n('dn', 0)), n_bits=n_bits)
]

# Make PomerolED solver object
ed = PomerolED(index_converter, verbose = True)

# pomerol_index()
for pom_i, (k, sn) in enumerate(product(range(len(epsilon)), ['dn', 'up'])):
    assert ed.pomerol_index(("B%i_%s" % (k, sn), 0)) == pom_i
assert ed.pomerol_index(("dn", 0)) == 2 * len(epsilon)
assert ed.pomerol_index(("up", 0)) == 2 * len(epsilon) + 1

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
ed.diagonalize(H, bosons=bosons)

# Compute G(i\omega)
G_iw = ed.G_iw(gf_struct, beta, n_iw)

# Compute G(\tau)
G_tau = ed.G_tau(gf_struct, beta, n_tau)

# Compute G(\omega)
G_w = ed.G_w(gf_struct, beta, energy_window, n_w, 0.01)

if mpi.is_master_node():
    with HDFArchive('bosons.out.h5', 'w') as ar:
        ar['H'] = H
        ar['frequencies'] = [b.frequency for b in bosons]
        ar['couplings'] = [b.coupling for b in bosons]
        ar['n_bits'] = [b.n_bits for b in bosons]
        ar['G_iw'] = G_iw
        ar['G_tau'] = G_tau
        ar['G_w'] = G_w

    with HDFArchive("bosons.ref.h5", 'r') as ar:
        assert (ar['H'] - H).is_zero()
        assert ar['frequencies'] == [b.frequency for b in bosons]
        assert ar['couplings'] == [b.coupling for b in bosons]
        assert ar['n_bits'] == [b.n_bits for b in bosons]
        assert_block_gfs_are_close(ar['G_iw'], G_iw)
        assert_block_gfs_are_close(ar['G_tau'], G_tau)
        assert_block_gfs_are_close(ar['G_w'], G_w)
