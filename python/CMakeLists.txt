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
