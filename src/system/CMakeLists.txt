cmake_minimum_required(VERSION 3.17)
add_library(system STATIC)

target_sources(system
  PRIVATE
    Battery.cpp
  PUBLIC
    Battery.hpp
  )

target_include_directories(system PUBLIC .)
target_link_libraries(system PRIVATE power_conversion)