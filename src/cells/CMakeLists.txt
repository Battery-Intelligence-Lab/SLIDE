cmake_minimum_required(VERSION 3.17)
add_library(cells "")

add_subdirectory(Cell_SPM)

target_sources(cells
  PRIVATE
  cell_util.cpp
  PUBLIC
    cell_util.hpp
    cell_limits.hpp
  )

target_include_directories(cells  PUBLIC .)

target_link_libraries(cells
  PRIVATE
    Cell_SPM
  )