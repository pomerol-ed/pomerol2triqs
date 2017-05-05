from pytriqs.archive import HDFArchive
from pytriqs.gf import *
from pytriqs.operators import Operator, c, c_dag, n
from pytriqs.operators.util.op_struct import set_operator_structure, get_mkind
from pytriqs.operators.util.hamiltonians import h_int_kanamori
from pytriqs.utility import mpi
from pytriqs.applications.impurity_solvers.pomerol2triqs import PomerolED
import numpy as np
from itertools import product

# Input parameters
beta = 10.0
num_orb = 2
mu = 1.0
U = 2.0
J = 0.2
epsilon = [-2.3, 2.3]
V = [1.0*np.eye(num_orb) + 0.1*(np.ones((num_orb,num_orb)) - np.eye(num_orb))]*2

spin_names = ("up", "dn")
orb_names = range(num_orb)
n_iw = 1024

g2_n_iw = 5
g2_n_inu = 10
g2_n_l = 4
g2_blocks = set([("up","up"),("up","dn"),("dn","up")])

gf_struct = set_operator_structure(spin_names, orb_names, True)

# Conversion from TRIQS to Pomerol notation for operator indices
index_converter = {}
# Local degrees of freedom
index_converter.update({(sn, o) : ("loc", o, "down" if sn == "dn" else "up")
                        for sn, o in product(spin_names, orb_names)})
# Bath degrees of freedom
index_converter.update({("B%i_%s" % (k,sn), o) : ("bath" + str(k), o, "down" if sn == "dn" else "up")
                        for k, sn, o in product(range(len(epsilon)), spin_names, orb_names)})

ed = PomerolED(index_converter, verbose = True)

# Number of particles on the impurity
N = sum(n(sn, o) for sn, o in product(spin_names, orb_names))

# Local Hamiltonian
H_loc = h_int_kanamori(spin_names,orb_names,
                       np.array([[0,U-3*J],[U-3*J,0]]),
                       np.array([[U,U-2*J],[U-2*J,U]]),
                       J, True)
H_loc -= mu*N

# Bath Hamiltonian
H_bath = sum(eps*n("B%i_%s" % (k, sn), o)
             for sn, (k, eps), o in product(spin_names, enumerate(epsilon), orb_names))

# Hybridization Hamiltonian
H_hyb = Operator()
for k, v in enumerate(V):
    H_hyb = sum(        v[o1,o2]   * c_dag("B%i_%s" % (k, sn), o1) * c(sn, o2) +
                np.conj(v[o2,o1])  * c_dag(sn, o1) * c("B%i_%s" % (k, sn), o2)
                for sn, o1, o2 in product(spin_names, orb_names, orb_names))

# Complete Hamiltonian
H = H_loc + H_hyb + H_bath

# Diagonalize H
ed.diagonalize(H)
