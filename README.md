pomerol2triqs
=============

[![Build Status](https://travis-ci.org/krivenko/pomerol2triqs.svg?branch=master)](https://travis-ci.org/krivenko/pomerol2triqs)

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
- Install the [TRIQS](http://triqs.github.io/triqs/2.1.x/install.html) library version 2.1.0.
- `source <path_to_triqs_install_dir>/share/cpp2pyvars.sh`
- `source <path_to_triqs_install_dir>/share/triqsvars.sh`
- `git clone https://github.com/krivenko/pomerol2triqs.git pomerol2triqs.git`
- `mkdir pomerol2triqs.build && cd pomerol2triqs.build`
- `cmake ../pomerol2triqs.git -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<path_to_install_dir> -DPOMEROL_PATH=<path_to_pomerol_install_dir>`
- `make`
- `make test`
- `make install`

License
-------

Copyright (C) 2017-2019 Igor Krivenko <igor.s.krivenko @ gmail.com>

With contributions from Hugo U.R. Strand

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
