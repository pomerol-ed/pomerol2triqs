from pytriqs.archive import HDFArchive
from pytriqs.gf import *
from pytriqs.operators import Operator, c, c_dag, n
from pytriqs.utility import mpi
from pytriqs.applications.impurity_solvers.pomerol2triqs import PomerolED
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
epsilon = [-1.0, 0, 1.0]
# Hopping amplitudes
V = [0.5, 0.5, 0.5]

spin_names = ("up", "dn")

# Number of Matsubara frequencies for GF calculation
n_iw = 1024

# Number of imaginary time slices for GF calculation
n_tau = 10001

# Energy window for real frequency GF calculation
energy_window = (-5, 5)
# Number of frequency points for real frequency GF calculation
n_w = 1000

# Number of bosonic Matsubara frequencies for G^2 calculations
g2_n_iw = 5
# Number of fermionic Matsubara frequencies for G^2 calculations
g2_n_inu = 10
# Number of Legendre coefficients for G^2 calculations
g2_n_l = 10
# Block index combinations for G^2 calculations
g2_blocks = set([("up", "up"), ("up", "dn"), ("dn", "up")])

gf_struct = {'up' : [0], 'dn' : [0]}

# Conversion from TRIQS to Pomerol notation for operator indices
# TRIQS: block_name, inner_index
# Pomerol: site_label, orbital_index, spin_name
index_converter = {}

# Local degrees of freedom
index_converter.update({(sn, 0) : ("loc", 0, "down" if sn == "dn" else "up") for sn in spin_names})
# Bath degrees of freedom
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

###########
# G^{(2)} #
###########

common_g2_params = {'gf_struct' : gf_struct,
                    'beta' : beta,
                    'blocks' : g2_blocks,
                    'n_iw' : g2_n_iw}

###############################
# G^{(2)}(i\omega;i\nu,i\nu') #
###############################

# Compute G^{(2),ph}(i\omega;i\nu,i\nu'), AABB block order
G2_iw_inu_inup_ph_AABB = ed.G2_iw_inu_inup(channel = "PH",
                                           block_order = "AABB",
                                           n_inu = g2_n_inu,
                                           **common_g2_params)

# Compute G^{(2),ph}(i\omega;i\nu,i\nu'), ABBA block order
G2_iw_inu_inup_ph_ABBA = ed.G2_iw_inu_inup(channel = "PH",
                                           block_order = "ABBA",
                                           n_inu = g2_n_inu,
                                           **common_g2_params)

# Compute G^{(2),pp}(i\omega;i\nu,i\nu'), AABB block order
G2_iw_inu_inup_pp_AABB = ed.G2_iw_inu_inup(channel = "PP",
                                           block_order = "AABB",
                                           n_inu = g2_n_inu,
                                           **common_g2_params)

# Compute G^{(2),pp}(i\omega;i\nu,i\nu'), ABBA block order
G2_iw_inu_inup_pp_ABBA = ed.G2_iw_inu_inup(channel = "PP",
                                           block_order = "ABBA",
                                           n_inu = g2_n_inu,
                                           **common_g2_params)

################
# Save results #
################

if mpi.is_master_node():
    with HDFArchive('anderson.h5', 'w') as ar:
        ar['G_iw'] = G_iw
        ar['G_tau'] = G_tau
        ar['G_w'] = G_w
        ar['G2_iw_inu_inup_ph_AABB'] = G2_iw_inu_inup_ph_AABB
        ar['G2_iw_inu_inup_ph_ABBA'] = G2_iw_inu_inup_ph_ABBA
        ar['G2_iw_inu_inup_pp_AABB'] = G2_iw_inu_inup_pp_AABB
        ar['G2_iw_inu_inup_pp_ABBA'] = G2_iw_inu_inup_pp_ABBA
