cmake_minimum_required(VERSION 3.16)
add_library(unit_tests "")

target_sources(unit_tests
  PRIVATE
    Cell_Bucket_test.cpp
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/unit_tests.hpp
  )

target_include_directories(unit_tests 
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )




