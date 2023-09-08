
cmake_minimum_required(VERSION 3.0)

project(Alya2pos VERSION 0.1 LANGUAGES C CXX Fortran)
set (CMAKE_CXX_STANDARD 11)

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release")

# Include Directories (for all targets)
include_directories(./ ${CMAKE_BINARY_DIR})

# wildcard all the sources in src
file(GLOB SOURCES ./*.f90)

add_executable(alya2pos ${SOURCES})

target_compile_options(alya2pos PUBLIC -Mfree)
target_compile_options(alya2pos PUBLIC -cpp)
target_compile_options(alya2pos PUBLIC -Mbackslash)
