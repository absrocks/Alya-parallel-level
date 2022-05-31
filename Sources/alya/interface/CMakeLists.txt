cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(alya-module-interface LANGUAGES Fortran)

file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.f90 ${CMAKE_CURRENT_LIST_DIR}/*.c)
add_library(${PROJECT_NAME} ${LIBRARY_TYPE} ${SRC_FILES})

include(compiler)
set_flags()
set_defs()
set_openmpi()
include(modules)

include_directories(${alya-core_MODULE_DIR} ${MODULES_MODULE_DIR})

target_link_libraries(${PROJECT_NAME} alya-core ${MODULES_LIB})

set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(${PROJECT_NAME}_MODULE_DIR ${CMAKE_Fortran_MODULE_DIRECTORY} CACHE INTERNAL "${PROJECT_NAME}: module directory" FORCE)

install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION modules)
