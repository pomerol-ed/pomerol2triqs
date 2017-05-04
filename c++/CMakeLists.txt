add_library(pomerol2triqs_c main.cpp)

target_link_libraries(pomerol2triqs_c ${TRIQS_LIBRARY_ALL} ${pomerol_LIBRARIES})
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR} ${pomerol_INCLUDE_DIRS})
triqs_set_rpath_for_target(pomerol2triqs_c)

install(TARGETS pomerol2triqs_c DESTINATION lib)
install(FILES main.hpp DESTINATION include/pomerol2triqs)
