cmake_minimum_required(VERSION 3.1)

project(cantera_fwrapper LANGUAGES CXX)

file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/src/*)

include_directories(SYSTEM)

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} ${SRC_FILES})

set(CUSTOM_CXX_FLAGS "${CUSTOM_CXX_FLAGS} -std=c++0x")

include(compiler)
set_flags()
set_defs()
include(cantera)
cantera()

set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(${PROJECT_NAME}_INCLUDE_DIR ${HEADER_DIR} CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(FILES ${HEADER_FILES} DESTINATION include)

