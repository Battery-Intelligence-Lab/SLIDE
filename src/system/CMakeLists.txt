cmake_minimum_required(VERSION 3.16)
add_library(system OBJECT)

target_sources(system
  PRIVATE
    Battery.cpp
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/Battery.hpp
  )

target_include_directories(system 
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )




