from h5 import HDFArchive
from triqs.gf import *
from triqs.operators import Operator, c, c_dag, n
from triqs.utility import mpi
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

# Number of Matsubara frequencies for \chi^3 calculations
chi3_n_inu = 10
# Block index combinations for \chi^3 calculations
chi3_blocks = set([("up", "up"), ("up", "dn"), ("dn", "up")])

gf_struct = [('up', 1), ('dn', 1)]

# Conversion from TRIQS to Pomerol notation for operator indices
# TRIQS: block_name, inner_index
# Pomerol: site_label, orbital_index, spin_name
index_converter = {}

# Local degrees of freedom
index_converter.update({(sn, 0) : ("loc", 0, "down" if sn == "dn" else "up") for sn in spin_names})
# Bath degrees of freedom
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

# Compute occupations
occ = [ed.ensemble_average((s, 0), (s, 0), beta).real for s in spin_names]

# Compute G(i\omega)
G_iw = ed.G_iw(gf_struct, beta, n_iw)

# Compute G(\tau)
G_tau = ed.G_tau(gf_struct, beta, n_tau)

# Compute G(\omega)
G_w = ed.G_w(gf_struct, beta, energy_window, n_w, 0.01)

# Compute \chi(\tau) = <n_{up}(\tau) n_{dn}(0)>
chi_tau = ed.chi_tau(('up',0), ('up',0), ('dn',0), ('dn',0), beta, n_tau)

# Compute \chi(i\nu)
chi_inu = ed.chi_inu(('up',0), ('up',0), ('dn',0), ('dn',0), beta, n_iw)

# Compute \chi_c(\tau) = <n_{up}(\tau) n_{dn}(0)> - <n_{up}><n_{dn}>
chi_tau_c = ed.chi_tau(('up',0), ('up',0), ('dn',0), ('dn',0), beta, n_tau, True)

# Compute \chi_c(i\nu)
chi_inu_c = ed.chi_inu(('up',0), ('up',0), ('dn',0), ('dn',0), beta, n_iw, True)

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

##########################
# \chi^{(3)}(i\nu,i\nu') #
##########################

common_chi3_params = {'gf_struct' : gf_struct,
                      'beta' : beta,
                      'blocks' : chi3_blocks,
                      'n_inu' : chi3_n_inu}

# Compute \chi^{(3),pp}(i\nu,i\nu'), AABB block order
chi3_iw_inu_pp_AABB = ed.chi3_iw_inu(channel = "PP",
                                         block_order = "AABB",
                                         **common_chi3_params)

# Compute \chi^{(3),pp}(i\nu,i\nu'), ABBA block order
chi3_iw_inu_pp_ABBA = ed.chi3_iw_inu(channel = "PP",
                                         block_order = "ABBA",
                                         **common_chi3_params)

# Compute \chi^{(3),ph}(i\nu,i\nu'), AABB block order
chi3_iw_inu_ph_AABB = ed.chi3_iw_inu(channel = "PH",
                                         block_order = "AABB",
                                         **common_chi3_params)

# Compute \chi^{(3),ph}(i\nu,i\nu'), ABBA block order
chi3_iw_inu_ph_ABBA = ed.chi3_iw_inu(channel = "PH",
                                         block_order = "ABBA",
                                         **common_chi3_params)

# Compute \chi^{(3),xph}(i\nu,i\nu'), AABB block order
chi3_iw_inu_xph_AABB = ed.chi3_iw_inu(channel = "xPH",
                                          block_order = "AABB",
                                          **common_chi3_params)

# Compute \chi^{(3),xph}(i\nu,i\nu'), ABBA block order
chi3_iw_inu_xph_ABBA = ed.chi3_iw_inu(channel = "xPH",
                                          block_order = "ABBA",
                                          **common_chi3_params)

################
# Save results #
################

if mpi.is_master_node():
    with HDFArchive('anderson.h5', 'w') as ar:
        ar['occ'] = occ
        ar['G_iw'] = G_iw
        ar['G_tau'] = G_tau
        ar['G_w'] = G_w
        ar['chi_tau'] = chi_tau
        ar['chi_inu'] = chi_inu
        ar['chi_tau_c'] = chi_tau_c
        ar['chi_inu_c'] = chi_inu_c
        ar['G2_iw_inu_inup_ph_AABB'] = G2_iw_inu_inup_ph_AABB
        ar['G2_iw_inu_inup_ph_ABBA'] = G2_iw_inu_inup_ph_ABBA
        ar['G2_iw_inu_inup_pp_AABB'] = G2_iw_inu_inup_pp_AABB
        ar['G2_iw_inu_inup_pp_ABBA'] = G2_iw_inu_inup_pp_ABBA
        ar['chi3_iw_inu_pp_AABB'] = chi3_iw_inu_pp_AABB
        ar['chi3_iw_inu_pp_ABBA'] = chi3_iw_inu_pp_ABBA
        ar['chi3_iw_inu_ph_AABB'] = chi3_iw_inu_ph_AABB
        ar['chi3_iw_inu_ph_ABBA'] = chi3_iw_inu_ph_ABBA
        ar['chi3_iw_inu_xph_AABB'] = chi3_iw_inu_xph_AABB
        ar['chi3_iw_inu_xph_ABBA'] = chi3_iw_inu_xph_ABBA
