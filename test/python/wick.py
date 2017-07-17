from pytriqs.archive import HDFArchive
from pytriqs.gf import *
from pytriqs.operators import Operator, c, c_dag, n
from pytriqs.utility import mpi
from pytriqs.applications.impurity_solvers.pomerol2triqs import PomerolED
from pytriqs.utility.comparison_tests import *
import numpy as np
from itertools import product

# Non-interacting impurity in a bath

####################
# Input parameters #
####################

beta = 10.0             # Inverse temperature
e_d = 0.1               # Local level
h = 0.15                # Magnetic field

# Levels of the bath sites
epsilon = [-1.0, 1.0]
# Hopping amplitudes
V = [0.5, 0.5]

spin_names = ("up", "dn")

# Number of Matsubara frequencies for GF calculation
n_iw = 200

# Number of bosonic Matsubara frequencies for G^2 calculations

g2_n_iw = 5

# Number of fermionic Matsubara frequencies for G^2 calculations
g2_n_inu = 5

# GF structure
gf_struct = {'up' : [0], 'dn' : [0]}

# Conversion from TRIQS to Pomerol notation for operator indices
index_converter = {}
index_converter.update({(sn, 0) : ("loc", 0, "down" if sn == "dn" else "up") for sn in spin_names})
index_converter.update({("B%i_%s" % (k, sn), 0) : ("bath" + str(k), 0, "down" if sn == "dn" else "up")
                        for k, sn in product(range(len(epsilon)), spin_names)})

# Make PomerolED solver object
ed = PomerolED(index_converter, verbose = True)

# Local Hamiltonian
H_loc = e_d*(n('up', 0) + n('dn', 0)) + h*(n('up', 0) - n('dn', 0))

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

###########
# G^{(2)} #
###########

common_g2_params = {'gf_struct' : gf_struct,
                    'beta' : beta,
                    'n_iw' : g2_n_iw}

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

b_mesh = MeshImFreq(beta = beta, S = 'Boson', n_max = g2_n_iw)
f_mesh = MeshImFreq(beta = beta, S = 'Fermion', n_max = g2_n_inu)
g2_mesh = MeshProduct(b_mesh, f_mesh, f_mesh)

# Check meshes of G2 containers
for bn1, bn2 in product(spin_names, spin_names):
    assert G2_iw_inu_inup_ph_AABB[bn1, bn2].mesh == g2_mesh
    assert G2_iw_inu_inup_ph_ABBA[bn1, bn2].mesh == g2_mesh
    assert G2_iw_inu_inup_pp_AABB[bn1, bn2].mesh == g2_mesh
    assert G2_iw_inu_inup_pp_ABBA[bn1, bn2].mesh == g2_mesh

G2_iw_inu_inup_ph_AABB_wick = G2_iw_inu_inup_ph_AABB.copy()
G2_iw_inu_inup_ph_ABBA_wick = G2_iw_inu_inup_ph_AABB.copy()
G2_iw_inu_inup_pp_AABB_wick = G2_iw_inu_inup_ph_AABB.copy()
G2_iw_inu_inup_pp_ABBA_wick = G2_iw_inu_inup_ph_AABB.copy()

G = lambda s, i: G_iw[s].data[i + n_iw, 0, 0]

g2_mesh_enum = product(enumerate(b_mesh), enumerate(f_mesh), enumerate(f_mesh))
for (m, w), (i, nu), (ip, nup) in g2_mesh_enum:
    m -= g2_n_iw - 1
    i -= g2_n_inu
    ip -= g2_n_inu
    d_w_0 = int(m == 0)
    d_w_nu_nup = int(m == i + ip + 1)
    d_nu_nup = int(i == ip)


    for s1, s2 in product(spin_names, spin_names):
        d_s1_s2 = int(s1 == s2)

        G2_iw_inu_inup_ph_AABB_wick[s1, s2].data[m + (g2_n_iw - 1), i + g2_n_inu, ip + g2_n_inu] = \
            beta * d_w_0 * G(s1, i) * G(s2, ip) - beta * d_nu_nup * d_s1_s2 * G(s1, m + i) * G(s2, i)
        G2_iw_inu_inup_ph_ABBA_wick[s1, s2].data[m + (g2_n_iw - 1), i + g2_n_inu, ip + g2_n_inu] = \
            beta * d_w_0 * d_s1_s2 * G(s1, i) * G(s2, ip) - beta * d_nu_nup * G(s2, m + i) * G(s1, i)
        G2_iw_inu_inup_pp_AABB_wick[s1, s2].data[m + (g2_n_iw - 1), i + g2_n_inu, ip + g2_n_inu] = \
            beta * d_w_nu_nup * G(s1, i) * G(s2, ip) - beta * d_nu_nup * d_s1_s2 * G(s1, m - i - 1) * G(s2, i)
        G2_iw_inu_inup_pp_ABBA_wick[s1, s2].data[m + (g2_n_iw - 1), i + g2_n_inu, ip + g2_n_inu] = \
            beta * d_w_nu_nup * d_s1_s2 * G(s1, i) * G(s2, ip) - beta * d_nu_nup * G(s2, m - i - 1) * G(s1, i)

for s1, s2 in product(spin_names, spin_names):
    assert_gfs_are_close(G2_iw_inu_inup_ph_AABB_wick[s1, s2], G2_iw_inu_inup_ph_AABB[s1, s2])
    assert_gfs_are_close(G2_iw_inu_inup_ph_ABBA_wick[s1, s2], G2_iw_inu_inup_ph_ABBA[s1, s2])
    assert_gfs_are_close(G2_iw_inu_inup_pp_AABB_wick[s1, s2], G2_iw_inu_inup_pp_AABB[s1, s2])
    assert_gfs_are_close(G2_iw_inu_inup_pp_ABBA_wick[s1, s2], G2_iw_inu_inup_pp_ABBA[s1, s2])
