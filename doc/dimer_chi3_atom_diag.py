# Call atom_diag to generate reference data for chi3.py unit test

from h5 import HDFArchive
from triqs.operators import c, c_dag, n
from triqs.gf import *
from triqs.atom_diag import *
from itertools import product
from copy import deepcopy
import numpy as np

####################
# Input parameters #
####################

beta = 2.0              # Inverse temperature
eps = [-1.9, -2.1]      # Energy levels of the atoms
t = 0.5                 # Hopping matrix element
U = 4.0                 # Coulomb repulsion
h_field = 0.05          # Magnetic field

spin_names = ("up", "dn")
atoms = (0, 1)

fops = [(sn, a) for sn, a in product(spin_names, atoms)]

# Number of bosonic Matsubara frequencies for susceptibility calculation
n_iw = 2
# Number of fermionic Matsubara frequencies for susceptibility calculation
n_inu = 2

# Perform computation for both atoms in the 'up' block and only for the 0-th
# atom in the 'down' block.
bs = {"up": 2, "dn": 1}
def make_block2_gf(AABB):
    mesh_b = MeshImFreq(beta=beta, S="Boson", n_iw=n_iw)
    mesh_f = MeshImFreq(beta=beta, S="Fermion", n_iw=n_inu)
    mesh = MeshProduct(mesh_b, mesh_f)
    print(mesh)
    blocks = []
    for sp1 in spin_names:
        blocks.append([])
        for sp2 in spin_names:
            if AABB:
                shape = (bs[sp1], bs[sp1], bs[sp2], bs[sp2])
            else:
                shape = (bs[sp1], bs[sp2], bs[sp2], bs[sp1])
            blocks[-1].append(Gf(mesh=mesh, target_shape=shape))
    return Block2Gf(spin_names, spin_names, blocks)

# Hamiltonian
H = sum(e * n(s, a) for (a, e), s in product(zip(atoms, eps), spin_names))
H += sum(-h_field * (n('up', a) - n('dn', a)) for a in atoms)
H += U * sum(n('up', a) * n('dn', a) for a in atoms)
H += t * sum((c_dag(sp, 0) * c(sp, 1) + c_dag(sp, 1) * c(sp, 0)) for sp in spin_names)

# Exact diagonalization without partitioning into sectors
ad = AtomDiag(H, fops, [])

# Extract energy levels and statistical weights
E = ad.energies[0]
w = np.diag(atomic_density_matrix(ad, beta)[0])

# Function f_{234} from Dominik Kiese's notes
def f(i, j, k, w1, w2):
    res = (w[j] + w[i]) / (E[i] - E[j] - 1j*w1)
    if np.isclose(E[i], E[k], atol=1e-13):
        res += beta * w[i] * np.isclose(w1, -w2, atol=1e-13)
    else:
        res += (w[k] - w[i]) / (E[i] - E[k] - 1j*w1 - 1j*w2)
    res /= -(E[j] - E[k] - 1j*w2)
    return res

# \chi^{(3)}_{pp}
def chi3_pp(x1p, x1, x2p, x2, w1, w2):
    Delta = ad.c_matrix(fops.index(x1), 0) @ ad.c_matrix(fops.index(x2), 0)
    c_dag_1 = ad.cdag_matrix(fops.index(x1p), 0)
    c_dag_2 = ad.cdag_matrix(fops.index(x2p), 0)
    res = 0
    for i, j, k in product(range(ad.get_subspace_dim(0)), repeat=3):
        res += f(i, j, k, w1, w2) * c_dag_1[i, j] * c_dag_2[j, k] * Delta[k, i]
        res += -f(i, j, k, w2, w1) * c_dag_2[i, j] * c_dag_1[j, k] * Delta[k, i]
    return res

# \chi^{(3)}_{ph}
def chi3_ph(x1p, x1, x2p, x2, w1, w2):
    N = ad.cdag_matrix(fops.index(x2p), 0) @ ad.c_matrix(fops.index(x2), 0)
    Cdag = ad.cdag_matrix(fops.index(x1p), 0)
    C = ad.c_matrix(fops.index(x1), 0)
    res = 0
    for i, j, k in product(range(ad.get_subspace_dim(0)), repeat=3):
        res += -f(i, j, k, w1, -w2) * Cdag[i, j] * C[j, k] * N[k, i]
        res += f(i, j, k, -w2, w1) * C[i, j] * Cdag[j, k] * N[k, i]
    return res

# \chi^{(3)}_{xph}
def chi3_xph(x1p, x1, x2p, x2, w1, w2):
    barN = -ad.cdag_matrix(fops.index(x2p), 0) @ ad.c_matrix(fops.index(x1), 0)
    Cdag = ad.cdag_matrix(fops.index(x1p), 0)
    C = ad.c_matrix(fops.index(x2), 0)
    res = 0
    for i, j, k in product(range(ad.get_subspace_dim(0)), repeat=3):
        res += -f(i, j, k, w1, -w2) * Cdag[i, j] * C[j, k] * barN[k, i]
        res += f(i, j, k, -w2, w1) * C[i, j] * Cdag[j, k] * barN[k, i]
    return res

def make_chi3(f, PP, AABB):
    print("PP =", PP, "AABB =", AABB)
    chi3 = make_block2_gf(AABB)
    for sp1, sp2 in product(spin_names, spin_names):
        block = chi3[sp1, sp2]
        for a1, a2, a3, a4 in product(*map(range, block.target_shape)):
            if AABB:
                indices = [(sp1, a1), (sp1, a2), (sp2, a3), (sp2, a4)]
            else:
                indices = [(sp1, a1), (sp2, a2), (sp2, a3), (sp1, a4)]
            print(indices)
            for iw, inu in block.mesh:
                print(iw, inu)
                nu1 = inu.imag
                nu2 = (iw.imag - nu1) if PP else (iw.imag + nu1)
                block[a1, a2, a3, a4][iw, inu] = f(*indices, nu1, nu2)
    return chi3

chi3_pp_AABB = make_chi3(chi3_pp, True, True)
chi3_pp_ABBA = make_chi3(chi3_pp, True, False)
chi3_ph_AABB = make_chi3(chi3_ph, False, True)
chi3_ph_ABBA = make_chi3(chi3_ph, False, False)
chi3_xph_AABB = make_chi3(chi3_xph, False, True)
chi3_xph_ABBA = make_chi3(chi3_xph, False, False)

with HDFArchive('dimer_chi3.ref.h5', 'w') as ar:
    ar['H'] = H
    ar['chi3_pp_AABB'] = chi3_pp_AABB
    ar['chi3_pp_ABBA'] = chi3_pp_ABBA
    ar['chi3_ph_AABB'] = chi3_ph_AABB
    ar['chi3_ph_ABBA'] = chi3_ph_ABBA
    ar['chi3_xph_AABB'] = chi3_xph_AABB
    ar['chi3_xph_ABBA'] = chi3_xph_ABBA
