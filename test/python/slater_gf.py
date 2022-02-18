from h5 import HDFArchive
from triqs.gf import *
from triqs.operators import *
from triqs.operators.util.op_struct import set_operator_structure, get_mkind
from triqs.operators.util.U_matrix import U_matrix
from triqs.operators.util.hamiltonians import h_int_slater
from triqs.operators.util.observables import N_op, S_op, L_op
from triqs.utility import mpi
from triqs.utility.comparison_tests import *
from pomerol2triqs import PomerolED
from itertools import product

# 5-orbital atom with Slater interaction term

####################
# Input parameters #
####################

L = 2                   # Angular momentum
beta = 10.0             # Inverse temperature
mu = 32.5               # Chemical potential (3 electrons in 5 bands)
# Slater parameters
U = 5.0
J = 0.1
F0 = U
F2 = J*(14.0/(1.0 + 0.625))
F4 = F2*0.625

spin_names = ("up", "dn")
num_orb = 2*L+1
orb_names = list(range(num_orb))
U_mat = U_matrix(L, radial_integrals = [F0,F2,F4], basis='spherical')

# Do not split H into blocks
ignore_symmetries = False

# Number of Matsubara frequencies for GF calculation
n_iw = 200

# Number of imaginary time slices for GF calculation
n_tau = 1001

# Energy window for real frequency GF calculation
energy_window = (-5, 5)
# Number of frequency points for real frequency GF calculation
n_w = 1000

# GF structure
gf_struct = set_operator_structure(spin_names, num_orb, False)

mkind = get_mkind(False, None)
# Conversion from TRIQS to Pomerol notation for operator indices
index_converter = {mkind(sn, on) : ("atom", on, "down" if sn == "dn" else "up")
                   for sn, on in product(spin_names, orb_names)}

# Make PomerolED solver object
ed = PomerolED(index_converter, verbose = True)

# Hamiltonian
H = h_int_slater(spin_names, orb_names, U_mat, False)

# Number of particles
N = N_op(spin_names, orb_names, False)

# z-component of spin
Sz = S_op('z', spin_names, orb_names, False)

# z-component of angular momentum
Lz = L_op('z', spin_names, orb_names, off_diag = False, basis = 'spherical')

# Double check that we are actually using integrals of motion
h_comm = lambda op: H*op - op*H
assert h_comm(N).is_zero()
assert h_comm(Sz).is_zero()
assert h_comm(Lz).is_zero()

# Diagonalize H
ed.diagonalize(H, ignore_symmetries)

# Compute G(i\omega)
G_iw = ed.G_iw(gf_struct, beta, n_iw)

# Compute G(\tau)
G_tau = ed.G_tau(gf_struct, beta, n_tau)

if mpi.is_master_node():
    with HDFArchive('slater_gf.out.h5', 'w') as ar:
        ar['H'] = H
        ar['G_iw'] = G_iw
        ar['G_tau'] = G_tau

    with HDFArchive("slater_gf.ref.h5", 'r') as ar:
        assert (ar['H'] - H).is_zero()
        assert_block_gfs_are_close(ar['G_iw'], G_iw)
        assert_block_gfs_are_close(ar['G_tau'], G_tau)

