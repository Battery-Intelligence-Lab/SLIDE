cmake_minimum_required(VERSION 3.17)
add_library(power_conversion STATIC)

target_sources(power_conversion
  PRIVATE
    Converter.cpp
  PUBLIC
    Converter.hpp
  )

target_include_directories(power_conversion PUBLIC .)