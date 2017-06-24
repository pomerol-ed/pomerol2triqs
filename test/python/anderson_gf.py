from pytriqs.archive import HDFArchive
from pytriqs.gf import *
from pytriqs.operators import Operator, c, c_dag, n
from pytriqs.utility import mpi
from pytriqs.applications.impurity_solvers.pomerol2triqs import PomerolED
from pytriqs.utility.comparison_tests import *
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
gf_struct = {'up' : [0], 'dn' : [0]}

# Conversion from TRIQS to Pomerol notation for operator indices
index_converter = {}
index_converter.update({(sn, 0) : ("loc", 0, "down" if sn == "dn" else "up") for sn in spin_names})
index_converter.update({("B%i_%s" % (k, sn), 0) : ("bath" + str(k), 0, "down" if sn == "dn" else "up")
                        for k, sn in product(range(len(epsilon)), spin_names)})

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
