cmake_minimum_required(VERSION 3.16)
add_library(Cell_PbA "")

target_sources(Cell_PbA
  PRIVATE
  Cell_PbA.cpp
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/Cell_PbA.hpp
    ${CMAKE_CURRENT_LIST_DIR}/State_PbA.hpp
  )

target_include_directories(Cell_Bucket
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )