cmake_minimum_required(VERSION 3.21)
add_library(factories OBJECT)

target_sources(factories
  PUBLIC
    makeBattery.cpp
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
