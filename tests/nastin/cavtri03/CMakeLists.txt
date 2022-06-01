cmake_minimum_required(VERSION 3.1)

project(cavtri03)

add_test(NAME "${PROJECT_NAME}_sequential" COMMAND ${CMAKE_SOURCE_DIR}/test/alya/alya ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_BINARY_DIR}/sequential" "${CMAKE_BINARY_DIR}/src/alya/alya" ${PROJECT_NAME})
