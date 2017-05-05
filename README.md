Quick and dirty TRIQS wrapper around the Pomerol exact diagonalization library.

Copyright (C) 2017 I. Krivenko

Installation
------------

- Install the latest version of [Pomerol](http://aeantipov.github.io/pomerol/) exact diagonalization library (`master` branch).
- Install the latest version of [TRIQS](https://triqs.ipht.cnrs.fr/1.x/install.html) library (**`unstable`** branch).
- `git clone git@bitbucket.org:igork/pomerol2triqs.git pomerol2triqs.git`.
- `mkdir pomerol2triqs.build && cd pomerol2triqs.build`.
- `cmake ../pomerol2triqs.git -DCMAKE_BUILD_TYPE=Release -DTRIQS_PATH=<path_to_triqs_install_dir> -DPOMEROL_PATH=<path_to_pomerol_install_dir>`.
- `make`.
- `make install`.
