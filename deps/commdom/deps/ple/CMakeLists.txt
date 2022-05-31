cmake_minimum_required(VERSION 3.1)

project(ple LANGUAGES C)

file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*.c)
set(HEADER_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
file(GLOB HEADER_FILES ${HEADER_DIR}/*.h)

include_directories(SYSTEM ${HEADER_DIR})

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} ${SRC_FILES} ${HEADER_FILES})

include(compiler)
set_flags()
set_defs()
set_mpi_c()

set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(${PROJECT_NAME}_INCLUDE_DIR ${HEADER_DIR} CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(FILES ${HEADER_FILES} DESTINATION include)
