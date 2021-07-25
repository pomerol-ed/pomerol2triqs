pomerol2triqs
=============

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3908394.svg)](https://doi.org/10.5281/zenodo.3908394)
[![CI](https://github.com/krivenko/pomerol2triqs/actions/workflows/CI.yml/badge.svg)](https://github.com/krivenko/pomerol2triqs/actions/workflows/CI.yml)


Quick and dirty TRIQS wrapper around the Pomerol exact diagonalization library.

To learn how to use this wrapper, see `example` subdir in the source directory.

Features
--------

* Diagonalization of finite fermionic models with Hamiltonians written in terms of second quantization operators.
* Calculation of single-particle Green's functions: `G(\tau)`, `G(i\omega_n)`, `G(\omega)`.
* Calculation of two-particle Green's functions: `G(\omega;\nu,\nu')` and `G(\omega;\ell,\ell')`.
* Calculation of ensemble averages of quadratic operators, `\langle c^\dagger_i c_j \rangle`.
* Calculation of dynamic susceptibilities, `\langle T c^\dagger_{i_1}(\tau) c_{j_1}(\tau) c^\dagger_{i_2}(0) c_{j_2}(0) \rangle`.

Notation for the two-particle Green's functions is adopted from the
[PhD thesis of Lewin Boehnke](http://ediss.sub.uni-hamburg.de/volltexte/2015/7325/pdf/Dissertation.pdf).

Installation
------------

- Install the latest version of [Pomerol](http://aeantipov.github.io/pomerol/) exact diagonalization library (`master` branch).
- Install the [TRIQS](http://triqs.github.io/triqs/3.0.x/install.html) library version 3.0.x.
- `source <path_to_triqs_install_dir>/share/triqsvars.sh`
- `git clone https://github.com/krivenko/pomerol2triqs.git pomerol2triqs.git`
- `mkdir pomerol2triqs.build && cd pomerol2triqs.build`
- `cmake ../pomerol2triqs.git -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<path_to_install_dir> -DPOMEROL_PATH=<path_to_pomerol_install_dir>`
- `make`
- `make test`
- `make install`

License
-------

Copyright (C) 2017-2020 Igor Krivenko <igor.s.krivenko @ gmail.com>

With contributions from Hugo U.R. Strand and Nils Wentzell.

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
