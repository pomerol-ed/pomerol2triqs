add_library(pomerol2triqs_c pomerol_ed.cpp)

target_link_libraries(pomerol2triqs_c ${TRIQS_LIBRARY_ALL} ${pomerol_LIBRARIES})
include_directories(${TRIQS_INCLUDE_ALL} ${pomerol_INCLUDE_DIRS})
triqs_set_rpath_for_target(pomerol2triqs_c)

install(TARGETS pomerol2triqs_c DESTINATION lib)
install(FILES pomerol_ed.hpp DESTINATION include/pomerol2triqs)
