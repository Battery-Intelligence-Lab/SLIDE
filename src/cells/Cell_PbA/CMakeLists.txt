cmake_minimum_required(VERSION 3.20)
add_library(Cell_PbA STATIC)

target_sources(Cell_PbA
  PRIVATE
  Cell_PbA.cpp
  PUBLIC
    Cell_PbA.hpp
    State_PbA.hpp
  )

target_include_directories(Cell_PbA  PUBLIC  .)