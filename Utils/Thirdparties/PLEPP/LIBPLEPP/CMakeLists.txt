
cmake_minimum_required(VERSION 3.0)

project(PLEPP VERSION 0.1 LANGUAGES C CXX Fortran)
set (CMAKE_CXX_STANDARD 11)

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release")

# Include Directories (for all targets)
include_directories(../LIBPLE ../LIBPLE/src Wrappers Include ${CMAKE_BINARY_DIR})

# wildcard all the sources in src
file(GLOB SOURCES 
	Wrappers/*.cpp
	Include/commdom.cpp
	Include/build_octree.c
       	)

add_library(commdom ${SOURCES})

