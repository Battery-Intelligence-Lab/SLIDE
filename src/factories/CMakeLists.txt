cmake_minimum_required(VERSION 3.16)
add_library(factories OBJECT)

target_sources(factories
  PRIVATE
    makeBattery.cpp
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/makeBattery.hpp
  )

target_include_directories(factories 
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )

  target_link_libraries(factories
  PRIVATE
    cells
    modules
  )
