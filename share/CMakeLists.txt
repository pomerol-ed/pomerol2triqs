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

if(NOT IS_SUBPROJECT
  AND NOT CMAKE_INSTALL_PREFIX STREQUAL TRIQS_ROOT
  AND NOT CMAKE_INSTALL_PREFIX STREQUAL "/usr"
  AND NOT CMAKE_INSTALL_PREFIX STREQUAL "/usr/local"
)

  if(PythonSupport)
    set(EXPORT_PYTHON_PATH "export PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${CPP2PY_PYTHON_LIB_DEST_ROOT}:$PYTHONPATH")
    set(MODFILE_PYTHON_PATH "prepend-path    PYTHONPATH          $root/${CPP2PY_PYTHON_LIB_DEST_ROOT}")
  endif()

  configure_file(${PROJECT_NAME}.modulefile.in ${PROJECT_NAME}.modulefile @ONLY)
  configure_file(${PROJECT_NAME}vars.sh.in ${PROJECT_NAME}vars.sh @ONLY)

  install(
    FILES
      ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.modulefile
      ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}vars.sh
    DESTINATION share/${PROJECT_NAME}
  )

  message(STATUS "***************************************************************")
  message(STATUS "* Custom install Location. Use:                                ")
  message(STATUS "*                                                              ")
  message(STATUS "*   source ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/${PROJECT_NAME}vars.sh      ")
  message(STATUS "*                                                              ")
  message(STATUS "* to set up the environment variables                          ")
  if(DEFINED ENV{MODULEPATH})
    message(STATUS "*                                                                          ")
    message(STATUS "* Consider copying ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.modulefile  ")
    message(STATUS "*   into your environment module directories                               ")
  endif()
  message(STATUS "***************************************************************")

endif()
