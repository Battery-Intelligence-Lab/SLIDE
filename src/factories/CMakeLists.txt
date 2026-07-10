cmake_minimum_required(VERSION 3.31)
add_library(factories OBJECT)

target_sources(factories
  PRIVATE
    makeBattery.cpp
  PUBLIC
    makeBattery.hpp
  )

target_include_directories(factories PUBLIC .)

target_link_libraries(factories
PUBLIC
  cells
  modules
  procedures
  settings
)
