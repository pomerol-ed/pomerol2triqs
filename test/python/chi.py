from h5 import HDFArchive
from triqs.gf import *
from triqs.operators import Operator, c, c_dag, n
from triqs.utility import mpi
from triqs.utility.comparison_tests import *
from pomerol2triqs import PomerolED
import numpy as np
from itertools import product

# Single orbital Anderson atom

####################
# Input parameters #
####################

beta = 10.0             # Inverse temperature
U = 1.0                 # Coulomb repulsion
mu = 0.4                # Chemical potential
h_field = 0.01

spin_names = ("up", "dn")

# Conversion from TRIQS to Pomerol notation for operator indices
index_converter = {}
index_converter.update({(sn, 0) : ("A", 0, "down" if sn == "dn" else "up") for sn in spin_names})

# Make PomerolED solver object
ed = PomerolED(index_converter, verbose = True)

# Hamiltonian
H = -mu*(n('up', 0) + n('dn', 0)) - h_field*(n('up', 0) - n('dn', 0))
H += U * n('up', 0) * n('dn', 0)

# Diagonalize H
ed.diagonalize(H)

# Compute reference values of <n('up', 0)> and <n('dn', 0)>
w = np.array([1.0,
              np.exp(beta*(mu+h_field)),
              np.exp(beta*(mu-h_field)),
              np.exp(-beta*(-2*mu+U))])
Z = np.sum(w)
w /= Z
n_up_ref, n_dn_ref = complex(w[1] + w[3]), complex(w[2] + w[3])

n_up = ed.ensemble_average(('up', 0), ('up', 0), beta)
n_dn = ed.ensemble_average(('dn', 0), ('dn', 0), beta)
S_p = ed.ensemble_average(('up', 0), ('dn', 0), beta)
S_m = ed.ensemble_average(('dn', 0), ('up', 0), beta)

assert abs(n_up - n_up_ref) < 1e-10
assert abs(n_dn - n_dn_ref) < 1e-10
assert abs(S_p) < 1e-10
assert abs(S_m) < 1e-10

# Compute 3 susceptibilities
# < n_up ; n_up >
# < n_up ; n_dn >
# < S_+ ; S_- >

# Number of Matsubara frequencies for susceptibility calculation
n_iw = 200

chi_up_up = ed.chi_iw(('up',0), ('up',0), ('up',0), ('up',0), beta, n_iw, True)
chi_up_dn = ed.chi_iw(('up',0), ('up',0), ('dn',0), ('dn',0), beta, n_iw, True)
chi_Sp_Sm = ed.chi_iw(('up',0), ('dn',0), ('dn',0), ('up',0), beta, n_iw, True)

for inu in chi_up_up.mesh:
  zero_freq = abs(complex(inu)) < 1e-10

  chi_up_up_ref = (w[1] + w[3]) * (1 - w[1] - w[3]) * beta if zero_freq else 0
  assert abs(chi_up_up[inu] - chi_up_up_ref) < 1e-10

  chi_up_dn_ref = (w[3] - (w[1] + w[3]) * (w[2] + w[3])) * beta if zero_freq else 0
  assert abs(chi_up_dn[inu] - chi_up_dn_ref) < 1e-10

  if h_field == 0:
    chi_Sp_Sm_ref = w[1] * beta if zero_freq else 0
  else:
    chi_Sp_Sm_ref = -(w[1] - w[2]) / (inu - 2*h_field)
  assert abs(chi_Sp_Sm[inu] - chi_Sp_Sm_ref) < 1e-10

# Number of time slices for susceptibility calculation
n_tau = 200

chi_up_up = ed.chi_tau(('up',0), ('up',0), ('up',0), ('up',0), beta, n_tau, True)
chi_up_dn = ed.chi_tau(('up',0), ('up',0), ('dn',0), ('dn',0), beta, n_tau, True)

for tau in chi_up_up.mesh:
  chi_up_up_ref = (w[1] + w[3]) * (1 - w[1] - w[3])
  assert abs(chi_up_up[tau] - chi_up_up_ref) < 1e-10

  chi_up_dn_ref = (w[3] - (w[1] + w[3]) * (w[2] + w[3]))
  assert abs(chi_up_dn[tau] - chi_up_dn_ref) < 1e-10
