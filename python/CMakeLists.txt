# where will the python end up in triqs?
set(python_destination pytriqs/applications/impurity_solvers/pomerol2triqs)

# make the build_xxx, install python files...
triqs_prepare_local_pytriqs(${python_destination})

# sitecustomize.py for build
set(package_name "pytriqs.applications.impurity_solvers")
configure_file(${CMAKE_SOURCE_DIR}/cmake/sitecustomize.py ${CMAKE_CURRENT_BINARY_DIR}/sitecustomize.py @ONLY)

# Build the python module for the solver
triqs_python_extension(pomerol2triqs ${python_destination})
target_link_libraries(pomerol2triqs pomerol2triqs_c)
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR} ${pomerol_INCLUDE_DIRS})
triqs_set_rpath_for_target(pomerol2triqs)
install(TARGETS pomerol2triqs DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})
