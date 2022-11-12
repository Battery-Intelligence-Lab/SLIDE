cmake_minimum_required(VERSION 3.17)
add_library(unit STATIC)

target_sources(unit
  PRIVATE
    Cell_Bucket_test.cpp
  PUBLIC
    unit_tests.hpp
  )

target_include_directories(unit PUBLIC .)




