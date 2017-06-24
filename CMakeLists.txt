# Start configuration
cmake_minimum_required(VERSION 2.8.7)

# Version number of the application
project(pomerol2triqs CXX)
set(POMEROL2TRIQS_VERSION "0.1")

# Build the optimized version by default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# We use shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Append TRIQS installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)
# Append Pomerol installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${POMEROL_PATH}/share/pomerol/cmake)

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

if(NOT TRIQS_VERSION EQUAL 1.5)
    message(FATAL_ERROR "The application requires the TRIQS library version 1.5 (got ${TRIQS_VERSION})")
endif()

# Load Pomerol
find_package(pomerol REQUIRED)

# Find MPI
find_package(MPI)

# Get hash
triqs_get_git_hash(${CMAKE_SOURCE_DIR} "POMEROL2TRIQS")
if(${GIT_RESULT} EQUAL 0)
    message(STATUS "Hash: ${POMEROL2TRIQS_GIT_HASH}")
endif(${GIT_RESULT} EQUAL 0)

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

message(STATUS "TRIQS : Adding compilation flags detected by the library")
add_definitions(${TRIQS_CXX_DEFINITIONS})

option(Tests "Enable Tests" ON)

include_directories(c++)

# Compile C++ code
add_subdirectory(c++)

# Python interface
if (${TRIQS_WITH_PYTHON_SUPPORT})
    add_subdirectory(python)
    if (${Tests})
        enable_testing()
        add_subdirectory(test)
    endif()
endif()
