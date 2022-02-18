add_library(pomerol2triqs_c pomerol_ed.cpp g.cpp g2.cpp chi.cpp)

target_link_libraries(pomerol2triqs_c PUBLIC triqs pomerol::pomerol)
target_compile_options(pomerol2triqs_c PRIVATE
  -DPROJECT_GIT_HASH=${PROJECT_GIT_HASH}
  -DTRIQS_GIT_HASH=${TRIQS_GIT_HASH}
)

install(TARGETS pomerol2triqs_c DESTINATION lib EXPORT pomerol2triqs-targets)
install(EXPORT pomerol2triqs-targets DESTINATION share/cmake)

install(FILES pomerol_ed.hpp g2_parameters.hpp DESTINATION include/pomerol2triqs)
