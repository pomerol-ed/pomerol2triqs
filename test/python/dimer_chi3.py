from h5 import HDFArchive
from triqs.operators import c, c_dag, n
from triqs.utility import mpi
from triqs.utility.comparison_tests import *
from pomerol2triqs import PomerolED
from itertools import product

# Hubbard dimer

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

# Number of Matsubara frequencies for susceptibility calculation
n_iw = 10

# Block structure of \chi^3
gf_struct = [['up', 2], ['dn', 2]]

# Conversion from TRIQS to Pomerol notation for operator indices
index_converter = {}
index_converter.update({(sn, 0) : ("A", 0, "down" if sn == "dn" else "up") for sn in spin_names})
index_converter.update({(sn, 1) : ("B", 0, "down" if sn == "dn" else "up") for sn in spin_names})

# Make PomerolED solver object
ed = PomerolED(index_converter, verbose = True)

# Hamiltonian
H = sum(e * n(s, a) for (a, e), s in product(zip(atoms, eps), spin_names))
H += sum(-h_field * (n('up', a) - n('dn', a)) for a in atoms)
H += U * sum(n('up', a) * n('dn', a) for a in atoms)
H += t * sum((c_dag(sp, 0) * c(sp, 1) + c_dag(sp, 1) * c(sp, 0)) for sp in spin_names)

# Diagonalize H
ed.diagonalize(H)

##################
# Compute \chi^3 #
##################

params = {'gf_struct': gf_struct, 'beta': beta, 'n_inu': n_inu}

# Particle-particle channel
chi3_pp_AABB = ed.chi3_iw_inu(**params, channel='PP', block_order='AABB')
chi3_pp_ABBA = ed.chi3_iw_inu(**params, channel='PP', block_order='ABBA')

# Particle-hole channel
chi3_ph_AABB = ed.chi3_iw_inu(**params, channel='PH', block_order='AABB')
chi3_ph_ABBA = ed.chi3_iw_inu(**params, channel='PH', block_order='ABBA')

# Crossed particle-hole channel
chi3_xph_AABB = ed.chi3_iw_inu(**params, channel='xPH', block_order='AABB')
chi3_xph_ABBA = ed.chi3_iw_inu(**params, channel='xPH', block_order='ABBA')

if mpi.is_master_node():
    with HDFArchive('dimer_chi3_np%i.out.h5' % mpi.size, 'w') as ar:
        ar['H'] = H
        ar['chi3_pp_AABB'] = chi3_pp_AABB
        ar['chi3_pp_ABBA'] = chi3_pp_ABBA
        ar['chi3_ph_AABB'] = chi3_ph_AABB
        ar['chi3_ph_ABBA'] = chi3_ph_ABBA
        ar['chi3_xph_AABB'] = chi3_xph_AABB
        ar['chi3_xph_ABBA'] = chi3_xph_ABBA

    with HDFArchive("dimer_chi3.ref.h5", 'r') as ar:
        assert (ar['H'] - H).is_zero()
        for bn1, bn2 in product(spin_names, spin_names):
            assert_gfs_are_close(ar['chi3_pp_AABB'][bn1, bn2], chi3_pp_AABB[bn1, bn2])
            assert_gfs_are_close(ar['chi3_pp_ABBA'][bn1, bn2], chi3_pp_ABBA[bn1, bn2])
            assert_gfs_are_close(ar['chi3_ph_AABB'][bn1, bn2], chi3_ph_AABB[bn1, bn2])
            assert_gfs_are_close(ar['chi3_ph_ABBA'][bn1, bn2], chi3_ph_ABBA[bn1, bn2])
            assert_gfs_are_close(ar['chi3_xph_AABB'][bn1, bn2], chi3_xph_AABB[bn1, bn2])
            assert_gfs_are_close(ar['chi3_xph_ABBA'][bn1, bn2], chi3_xph_ABBA[bn1, bn2])
