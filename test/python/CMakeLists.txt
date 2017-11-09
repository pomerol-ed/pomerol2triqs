# Copy h5 files to binary dir
file(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

foreach(t anderson_gf slater_gf wick)
  add_test(NAME ${t} COMMAND ${CPP2PY_PYTHON_INTERPRETER} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.py)
  set_property(TEST ${t} PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} )
endforeach()

# MPI tests
macro(add_mpi_test basename filename numproc)
  add_test(NAME "${basename}_${numproc}"
           COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numproc} ${CPP2PY_PYTHON_INTERPRETER} ${CMAKE_CURRENT_SOURCE_DIR}/${filename})
  set_property(TEST "${basename}_${numproc}" PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH})
endmacro()

foreach(numproc 1 2 4)
  add_mpi_test(anderson_g2_matsubara anderson_g2_matsubara.py ${numproc})
endforeach()
