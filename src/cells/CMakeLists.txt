cmake_minimum_required(VERSION 3.16)
add_library(cells "")

add_subdirectory(Cell_SPM)

target_sources(cells
  PRIVATE
  cell_util.cpp
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/cell_util.hpp
    ${CMAKE_CURRENT_LIST_DIR}/cell_limits.hpp
  )

target_include_directories(cells
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )

target_link_libraries(cells
  PRIVATE
    Cell_SPM
  )