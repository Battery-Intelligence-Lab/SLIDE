cmake_minimum_required(VERSION 3.16)
add_library(power_conversion OBJECT)

target_sources(power_conversion
  PRIVATE
    Converter.cpp
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/Converter.hpp
  )

target_include_directories(power_conversion 
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )




