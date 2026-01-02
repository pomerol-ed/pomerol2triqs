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
from triqs.operators.util.hamiltonians import h_int_kanamori
from triqs.utility import mpi
from pomerol2triqs import PomerolED
import numpy as np
from itertools import product

# 2-orbital Hubbard-Kanamori atom

####################
# Input parameters #
####################

beta = 10.0         # Inverse temperature
num_orb = 2         # Number of orbitals
U = 2.0             # Coulomb repulsion
mu = 1.5            # Chemical potential
J = 0.2             # Hund coupling

spin_names = ("up", "dn")
orb_names = list(range(num_orb))

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

# Number of bosonic Matsubara frequencies for \chi^3 calculations
chi3_n_iw = 10
# Number of fermionic Matsubara frequencies for \chi^3 calculations
chi3_n_inu = 10
# Block index combinations for \chi^3 calculations
chi3_blocks = set([("up", "up"), ("up", "dn"), ("dn", "up")])

gf_struct = [("up", num_orb), ("dn", num_orb)]
print("Block structure of single-particle Green's functions:", gf_struct)

# Conversion from TRIQS to Pomerol notation for operator indices
# TRIQS: block_name, inner_index
# Pomerol: site_label, orbital_index, spin_name
index_converter = {(sn, o) : ("loc", o, "down" if sn == "dn" else "up")
                   for sn, o in product(spin_names, orb_names)}

# Make PomerolED solver object
ed = PomerolED(index_converter, verbose = True)

# Number of particles on the impurity
N = sum(n(sn, o) for sn, o in product(spin_names, orb_names))

# Hamiltonian
H = h_int_kanamori(spin_names, num_orb,
                   np.array([[0, U-3*J], [U-3*J, 0]]),
                   np.array([[U, U-2*J], [U-2*J, U]]),
                   J, True)
H -= mu*N

# Diagonalize H
ed.diagonalize(H)

# Compute occupations
occ = [ed.ensemble_average(i, i, beta).real for i in product(spin_names, orb_names)]

# Compute G(i\omega)
G_iw = ed.G_iw(gf_struct, beta, n_iw)

# Compute G(\tau)
G_tau = ed.G_tau(gf_struct, beta, n_tau)

# Compute G(\omega)
G_w = ed.G_w(gf_struct, beta, energy_window, n_w, 0.01)

# Compute \chi(\tau) = <n_{up,0}(\tau) n_{dn,0}(0)>
chi_tau = ed.chi_tau(('up',0), ('up',0), ('dn',0), ('dn',0), beta, n_tau)

# Compute \chi(i\nu)
chi_iw = ed.chi_iw(('up',0), ('up',0), ('dn',0), ('dn',0), beta, n_iw)

# Compute \chi(\omega)
chi_w = ed.chi_w(('up',0), ('up',0), ('dn',0), ('dn',0), beta, energy_window, n_w, 0.01)

# Compute \chi_c(\tau) = <n_{up,0}(\tau) n_{dn,0}(0)> - <n_{up,0}><n_{dn,0}>
chi_tau_c = ed.chi_tau(('up',0), ('up',0), ('dn',0), ('dn',0), beta, n_tau, True)

# Compute \chi_c(i\nu)
chi_iw_c = ed.chi_iw(('up',0), ('up',0), ('dn',0), ('dn',0), beta, n_iw, True)

# Compute \chi_c(\omega)
chi_w_c = ed.chi_w(('up',0), ('up',0), ('dn',0), ('dn',0), beta, energy_window, n_w, 0.01, True)

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

############################
# \chi^{(3)}(i\omega,i\nu) #
############################

common_chi3_params = {'gf_struct' : gf_struct,
                      'beta' : beta,
                      'blocks' : chi3_blocks,
                      'n_iw' : chi3_n_iw,
                      'n_inu' : chi3_n_inu}

# Compute \chi^{(3),pp}(i\omega,i\nu), AABB block order
chi3_iw_inu_pp_AABB = ed.chi3_iw_inu(channel = "PP",
                                     block_order = "AABB",
                                     **common_chi3_params)

# Compute \chi^{(3),pp}(i\omega,i\nu), ABBA block order
chi3_iw_inu_pp_ABBA = ed.chi3_iw_inu(channel = "PP",
                                     block_order = "ABBA",
                                     **common_chi3_params)

# Compute \chi^{(3),ph}(i\omega,i\nu), AABB block order
chi3_iw_inu_ph_AABB = ed.chi3_iw_inu(channel = "PH",
                                     block_order = "AABB",
                                     **common_chi3_params)

# Compute \chi^{(3),ph}(i\omega,i\nu), ABBA block order
chi3_iw_inu_ph_ABBA = ed.chi3_iw_inu(channel = "PH",
                                     block_order = "ABBA",
                                     **common_chi3_params)

# Compute \chi^{(3),xph}(i\omega,i\nu), AABB block order
chi3_iw_inu_xph_AABB = ed.chi3_iw_inu(channel = "xPH",
                                      block_order = "AABB",
                                      **common_chi3_params)

# Compute \chi^{(3),xph}(i\omega,i\nu), ABBA block order
chi3_iw_inu_xph_ABBA = ed.chi3_iw_inu(channel = "xPH",
                                      block_order = "ABBA",
                                      **common_chi3_params)

################
# Save results #
################

if mpi.is_master_node():
    with HDFArchive('2band.atom.h5', 'w') as ar:
        ar['occ'] = occ
        ar['G_iw'] = G_iw
        ar['G_tau'] = G_tau
        ar['G_w'] = G_w
        ar['chi_tau'] = chi_tau
        ar['chi_iw'] = chi_iw
        ar['chi_w'] = chi_w
        ar['chi_tau_c'] = chi_tau_c
        ar['chi_iw_c'] = chi_iw_c
        ar['chi_w_c'] = chi_w_c
        ar['G2_iw_inu_inup_ph_AABB'] = G2_iw_inu_inup_ph_AABB
        ar['G2_iw_inu_inup_ph_ABBA'] = G2_iw_inu_inup_ph_ABBA
        ar['G2_iw_inu_inup_pp_AABB'] = G2_iw_inu_inup_pp_AABB
        ar['G2_iw_inu_inup_pp_ABBA'] = G2_iw_inu_inup_pp_ABBA
        ar['G2_iw_l_lp_ph_AABB'] = G2_iw_l_lp_ph_AABB
        ar['G2_iw_l_lp_ph_ABBA'] = G2_iw_l_lp_ph_ABBA
        ar['G2_iw_l_lp_pp_AABB'] = G2_iw_l_lp_pp_AABB
        ar['G2_iw_l_lp_pp_ABBA'] = G2_iw_l_lp_pp_ABBA
        ar['chi3_iw_inu_pp_AABB'] = chi3_iw_inu_pp_AABB
        ar['chi3_iw_inu_pp_ABBA'] = chi3_iw_inu_pp_ABBA
        ar['chi3_iw_inu_ph_AABB'] = chi3_iw_inu_ph_AABB
        ar['chi3_iw_inu_ph_ABBA'] = chi3_iw_inu_ph_ABBA
        ar['chi3_iw_inu_xph_AABB'] = chi3_iw_inu_xph_AABB
        ar['chi3_iw_inu_xph_ABBA'] = chi3_iw_inu_xph_ABBA
