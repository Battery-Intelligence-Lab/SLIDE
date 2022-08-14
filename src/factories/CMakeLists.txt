cmake_minimum_required(VERSION 3.20)
add_library(factories STATIC)

target_sources(factories
  PRIVATE
    makeBattery.cpp
  PUBLIC
    makeBattery.hpp
  )

target_include_directories(factories  PUBLIC .)

target_link_libraries(factories
PRIVATE
  cells
  modules
)
