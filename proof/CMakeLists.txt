cmake_minimum_required(VERSION 3.12)

project(proof)
# set( CMAKE_VERBOSE_MAKEFILE on )

SET(LIBFQFFT ${CMAKE_SOURCE_DIR}/libfqfft/install/lib/libff.a)
SET(LIBFQFFTH ${CMAKE_SOURCE_DIR}/libfqfft/install/include)

SET(GCC_COVERAGE_COMPILE_FLAGS "-std=c++17 -fopenmp -pthread -Ofast -march=native -frename-registers -flto")
SET(GCC_COVERAGE_LINK_FLAGS    "")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")

enable_testing()
add_subdirectory(client)
add_subdirectory(main)

