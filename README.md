pomerol2triqs
=============

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5735413.svg)](https://doi.org/10.5281/zenodo.5735413)
[![CI](https://github.com/pomerol-ed/pomerol2triqs/actions/workflows/CI.yml/badge.svg)](https://github.com/pomerol-ed/pomerol2triqs/actions/workflows/CI.yml)


Quick and dirty TRIQS wrapper around the Pomerol exact diagonalization library.

To learn how to use this wrapper, see `example` subdir in the source directory.

Features
--------

* Diagonalization of finite fermionic models with Hamiltonians written in terms of second quantization operators.
* Calculation of single-particle Green's functions: $G(\tau)$, $G(i\omega_n)$, $G(\omega)$.
* Calculation of two-particle Green's functions: $G(\omega;\nu,\nu')$ and $G(\omega;\ell,\ell')$.
* Calculation of ensemble averages of quadratic operators, $\langle c^\dagger_i c_j \rangle$.
* Calculation of dynamic susceptibilities
  $\langle \mathbb{T} c^\dagger_i(\tau) c_j(0^+) c^\dagger_k(\tau) c_l(0) \rangle$,
  $\langle \mathbb{T} c^\dagger_i(\tau) c_j(\tau) c^\dagger_k(0^+) c_l(0) \rangle$,
  $\langle \mathbb{T} c^\dagger_i(\tau) c_j(0) c^\dagger_k(0^+) c_l(\tau) \rangle$
  in the Matsubara frequency representation.
* Calculation of 3-point fermion-boson susceptibilities
  $\langle \mathbb{T} c^\dagger_i(\tau) c_j(0^+) c^\dagger_k(\tau') c_l(0) \rangle$,
  $\langle \mathbb{T} c^\dagger_i(\tau) c_j(\tau') c^\dagger_k(0^+) c_l(0) \rangle$,
  $\langle \mathbb{T} c^\dagger_i(\tau) c_j(0) c^\dagger_k(0^+) c_l(\tau') \rangle$
  in the Matsubara frequency representation.

Notation for the two-particle Green's functions is adopted from the
[PhD thesis of Lewin Boehnke](http://ediss.sub.uni-hamburg.de/volltexte/2015/7325/pdf/Dissertation.pdf).

Installation
------------

- Install the [Pomerol](http://pomerol-ed.github.io/pomerol/) exact diagonalization library.
- Install the [TRIQS](http://triqs.github.io/triqs/latest/install.html) library.
- `source <path_to_triqs_install_dir>/share/triqs/triqsvars.sh`
  (skip this step if TRIQS has been installed in a system location, e.g. `/usr`
  or `/usr/local` on Linux).
- `git clone https://github.com/pomerol-ed/pomerol2triqs.git pomerol2triqs.git`
- `mkdir pomerol2triqs.build && cd pomerol2triqs.build`
- `cmake ../pomerol2triqs.git -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<path_to_install_dir> -DPOMEROL_PATH=<path_to_pomerol_install_dir>`
- `make`
- `make test`
- `make install`

Pomerol and TRIQS version compatibility
---------------------------------------

| pomerol2triqs release | Pomerol version | TRIQS version |
|-----------------------|-----------------|---------------|
| v0.9                  | 2.1             | 3.2.x, 3.3.x  |
| v0.8                  | 2.1             | 3.1.x         |
| v0.7                  | 2.0             | 3.1.x         |
| v0.6                  | 2.0             | 3.0.x         |
| v0.5                  | 1.3             | 3.0.x         |
| v0.4                  | 1.3             | 2.2.x         |
| v0.3                  | 1.3             | 2.1.x         |

License
-------

Copyright (C) 2017-2024 Igor Krivenko <igor.s.krivenko @ gmail.com>

With contributions from

* Hugo U.R. Strand
* Nils Wentzell
* Dominik Kiese

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
