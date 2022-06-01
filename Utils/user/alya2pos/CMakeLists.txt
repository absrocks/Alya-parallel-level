cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(alya2pos LANGUAGES Fortran)

file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.f90 ${CMAKE_CURRENT_LIST_DIR}/*.c)
add_executable(${PROJECT_NAME} ${SRC_FILES})

set(WITH_STD2008 OFF)

include(compiler)
set_flags()
set_defs()

set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
