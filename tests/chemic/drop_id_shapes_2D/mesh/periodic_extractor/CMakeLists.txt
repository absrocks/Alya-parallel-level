cmake_minimum_required (VERSION 2.6)
project (periodic_extractor)

add_library(D3class D3class.cpp)

add_executable(per.x getPerNodes.cpp)
target_link_libraries (per.x D3class)
