from pytriqs.archive import HDFArchive
from pytriqs.gf import *
from pytriqs.operators import Operator, c, c_dag, n
from pytriqs.operators.util.hamiltonians import h_int_kanamori
from pytriqs.utility import mpi
from pytriqs.applications.impurity_solvers.pomerol2triqs import PomerolED
import numpy as np
from itertools import product

# 2-orbital impurity Anderson model (bath: 1 site * 2 orbitals)

####################
# Input parameters #
####################

beta = 10.0             # Inverse temperature
num_orb = 2             # Number of orbitals
U = 2.0                 # Coulomb repulsion
mu = 1.5                # Chemical potential
J = 0.2                 # Hund coupling

# Levels of the bath sites
epsilon = np.array([-0.2, 0.2])
# Hopping matrix
V = 0.7*np.eye(num_orb) + 0.1*(np.ones((num_orb, num_orb)) - np.eye(num_orb))

spin_names = ("up", "dn")
orb_names = range(num_orb)

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

gf_struct = {"up" : orb_names, "dn" : orb_names}
print "Block structure of single-particle Green's functions:", gf_struct

# Conversion from TRIQS to Pomerol notation for operator indices
# TRIQS: block_name, inner_index
# Pomerol: site_label, orbital_index, spin_name
index_converter = {}

# Local degrees of freedom
index_converter.update({(sn, o) : ("loc", o, "down" if sn == "dn" else "up")
                        for sn, o in product(spin_names, orb_names)})
# Bath degrees of freedom
index_converter.update({("B_" + sn, o) : ("bath", o, "down" if sn == "dn" else "up")
                        for sn, o in product(spin_names, orb_names)})

# Make PomerolED solver object
ed = PomerolED(index_converter, verbose = True)

# Number of particles on the impurity
N = sum(n(sn, o) for sn, o in product(spin_names, orb_names))

# Local Hamiltonian
H_loc = h_int_kanamori(spin_names, orb_names,
                       np.array([[0, U-3*J], [U-3*J, 0]]),
                       np.array([[U, U-2*J], [U-2*J, U]]),
                       J, True)
H_loc -= mu*N

# Bath Hamiltonian
H_bath = sum(epsilon[o] * n("B_" + sn, o) for sn, o in product(spin_names, orb_names))

# Hybridization Hamiltonian
H_hyb = sum(        V[o1,o2]  * c_dag("B_" + sn, o1) * c(sn, o2) +
            np.conj(V[o2,o1]) * c_dag(sn, o1) * c("B_" + sn, o2)
            for sn, o1, o2 in product(spin_names, orb_names, orb_names))

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

#########################
# G^{(2)}(i\omega;l,l') #
#########################

# Compute G^{(2),ph}(i\omega;l,l'), AABB block order
G2_iw_l_lp_ph_AABB = ed.G2_iw_l_lp(channel = "PH",
                                   block_order = "AABB",
                                   n_l = g2_n_l,
                                   **common_g2_params)

# Compute G^{(2),ph}(i\omega;l,l'), ABBA block order
G2_iw_l_lp_ph_ABBA = ed.G2_iw_l_lp(channel = "PH",
                                   block_order = "ABBA",
                                   n_l = g2_n_l,
                                   **common_g2_params)

# Compute G^{(2),pp}(i\omega;l,l'), AABB block order
G2_iw_l_lp_pp_AABB = ed.G2_iw_l_lp(channel = "PP",
                                   block_order = "AABB",
                                   n_l = g2_n_l,
                                   **common_g2_params)

# Compute G^{(2),pp}(i\omega;l,l'), ABBA block order
G2_iw_l_lp_pp_ABBA = ed.G2_iw_l_lp(channel = "PP",
                                   block_order = "ABBA",
                                   n_l = g2_n_l,
                                   **common_g2_params)

################
# Save results #
################

if mpi.is_master_node():
    with HDFArchive('2band.h5', 'w') as ar:
        ar['G_iw'] = G_iw
        ar['G_tau'] = G_tau
        ar['G_w'] = G_w
        ar['G2_iw_inu_inup_ph_AABB'] = G2_iw_inu_inup_ph_AABB
        ar['G2_iw_inu_inup_ph_ABBA'] = G2_iw_inu_inup_ph_ABBA
        ar['G2_iw_inu_inup_pp_AABB'] = G2_iw_inu_inup_pp_AABB
        ar['G2_iw_inu_inup_pp_ABBA'] = G2_iw_inu_inup_pp_ABBA
        ar['G2_iw_l_lp_ph_AABB'] = G2_iw_l_lp_ph_AABB
        ar['G2_iw_l_lp_ph_ABBA'] = G2_iw_l_lp_ph_ABBA
        ar['G2_iw_l_lp_pp_AABB'] = G2_iw_l_lp_pp_AABB
        ar['G2_iw_l_lp_pp_ABBA'] = G2_iw_l_lp_pp_ABBA
