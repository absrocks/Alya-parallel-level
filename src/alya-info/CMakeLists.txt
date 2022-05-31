cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(alya-info LANGUAGES Fortran C)

set(GITINFO_FILE "gitinfo.f90")

add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${GITINFO_FILE}" "${CMAKE_CURRENT_BINARY_DIR}/_${GITINFO_FILE}"
  COMMAND ${CMAKE_COMMAND}
    -Dlocal_dir="${CMAKE_SOURCE_DIR}/cmake"
    -Doutput_dir="${CMAKE_CURRENT_BINARY_DIR}"
    -Droot_dir="${CMAKE_SOURCE_DIR}"
    -P "${CMAKE_SOURCE_DIR}/cmake/gitinfo.cmake"
)

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} "${CMAKE_CURRENT_BINARY_DIR}/${GITINFO_FILE}")

include(compiler)
set_flags()
set_defs()
set_openmpi()
include(metis)
metis()
include(gmsh)
gmsh()
include(commdom)
commdom()
include(cantera)
cantera()
include(talp)
talp()

set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(${PROJECT_NAME}_MODULE_DIR ${CMAKE_Fortran_MODULE_DIRECTORY} CACHE INTERNAL "${PROJECT_NAME}: module directory" FORCE)

install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION modules)
