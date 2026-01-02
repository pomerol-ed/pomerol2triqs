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

add_library(pomerol2triqs_c pomerol_ed.cpp g.cpp g2.cpp chi.cpp chi3.cpp)

target_link_libraries(pomerol2triqs_c PUBLIC triqs pomerol::pomerol)
target_compile_options(pomerol2triqs_c PRIVATE
  -DPROJECT_GIT_HASH=${PROJECT_GIT_HASH}
  -DTRIQS_GIT_HASH=${TRIQS_GIT_HASH}
)

install(TARGETS pomerol2triqs_c DESTINATION lib EXPORT pomerol2triqs-targets)
install(EXPORT pomerol2triqs-targets DESTINATION share/cmake)

install(FILES pomerol_ed.hpp parameters.hpp DESTINATION include/pomerol2triqs)
