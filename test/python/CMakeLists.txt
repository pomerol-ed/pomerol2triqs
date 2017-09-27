find_package(TriqsTest)

# Copy h5 files to binary dir
FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

triqs_add_python_test(anderson_gf)
triqs_add_python_test(slater_gf)
triqs_add_python_test(wick)

foreach(TEST_MPI_NUMPROC 1 2 4)
    triqs_add_python_test(anderson_g2_matsubara)
endforeach()
