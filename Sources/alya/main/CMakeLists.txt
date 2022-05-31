cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(alya LANGUAGES Fortran)

file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.f90 ${CMAKE_CURRENT_LIST_DIR}/*.c)
add_executable(${PROJECT_NAME} ${SRC_FILES})

include(compiler)
set_flags()
set_defs()
set_openmpi()
include(commdom)
commdom()
include(cantera)
cantera()

include_directories(${alya-core_MODULE_DIR} ${alya-module-interface_MODULE_DIR})

target_link_libraries(${PROJECT_NAME} alya-core alya-module-interface c)

set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
