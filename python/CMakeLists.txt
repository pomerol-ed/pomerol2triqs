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

# Build the python module
add_cpp2py_module(pomerol2triqs)

target_link_libraries(pomerol2triqs triqs pomerol2triqs_c)
include_directories(${pomerol_INCLUDE_DIRS})

# Configure the version
configure_file(version.py.in version.py)

# Copy __init__.py
file(COPY __init__.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Install python module to proper location
set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/pomerol2triqs)
install(TARGETS pomerol2triqs DESTINATION ${PYTHON_LIB_DEST})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/__init__.py ${CMAKE_CURRENT_BINARY_DIR}/version.py
        DESTINATION ${PYTHON_LIB_DEST})
