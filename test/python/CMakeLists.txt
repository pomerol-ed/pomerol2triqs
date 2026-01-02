################################################################################
#
# pomerol2triqs
#
# Copyright (C) 2017-2026 Igor Krivenko
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

# Copy h5 files to binary dir
file(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

macro(add_mpi_test basename filename numproc)
  add_test(NAME "${basename}_np${numproc}"
           COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numproc} ${CPP2PY_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${filename})
  set_property(TEST "${basename}_np${numproc}" PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH})
endmacro()

# Single MPI rank tests
foreach(t anderson_gf slater_gf wick chi anomalous)
  add_mpi_test(${t} ${t}.py 1)
endforeach()

# Multiple MPI ranks tests
foreach(t anderson_g2_matsubara dimer_chi3)
  foreach(numproc 1 2 4)
    add_mpi_test(${t} ${t}.py ${numproc})
  endforeach()
endforeach()
