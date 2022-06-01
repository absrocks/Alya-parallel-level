
cmake_minimum_required(VERSION 3.0)

project(PLEPP VERSION 0.1 LANGUAGES C CXX Fortran)
set (CMAKE_CXX_STANDARD 11)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release")

# Include Directories (for all targets)
include_directories(../LIBPLE/include include ${CMAKE_BINARY_DIR})

# wildcard all the sources in src
file(GLOB SOURCES src/*.cpp src/*.c)

add_library(commdom ${SOURCES})

