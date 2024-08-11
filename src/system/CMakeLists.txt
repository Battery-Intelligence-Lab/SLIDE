cmake_minimum_required(VERSION 3.21)
add_library(system OBJECT)

target_sources(system
  PUBLIC
    Battery.cpp
    Battery.hpp
  )

target_include_directories(system PUBLIC .)
target_link_libraries(system PUBLIC power_conversion settings)