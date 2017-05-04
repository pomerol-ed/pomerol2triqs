# where will the python end up in triqs?
set(python_destination pytriqs/applications/impurity_solvers/pomerol2triqs)

# make the build_xxx, install python files...
triqs_prepare_local_pytriqs(${python_destination})

# Build the python module for the solver
triqs_python_extension(pomerol2triqs ${python_destination})
target_link_libraries(pomerol2triqs pomerol2triqs_c)
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR} ${pomerol_INCLUDE_DIRS})
triqs_set_rpath_for_target(pomerol2triqs)
install(TARGETS pomerol2triqs DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})
