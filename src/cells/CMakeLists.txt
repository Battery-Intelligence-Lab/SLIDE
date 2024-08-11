cmake_minimum_required(VERSION 3.21)
add_library(cells OBJECT)
set_target_properties(cells PROPERTIES LINKER_LANGUAGE CXX)

add_subdirectory(Cell_SPM)

target_sources(cells
  PUBLIC
    cell_limits.hpp
    cells.hpp
  )

target_include_directories(cells PUBLIC .)

target_link_libraries(cells
  PUBLIC
    Cell_SPM
  )