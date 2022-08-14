cmake_minimum_required(VERSION 3.20)
add_library(unit_tests STATIC)

target_sources(unit_tests
  PRIVATE
    Cell_Bucket_test.cpp
  PUBLIC
    unit_tests.hpp
  )

target_include_directories(unit_tests PUBLIC .)




