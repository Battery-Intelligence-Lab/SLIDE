cmake_minimum_required(VERSION 3.16)
add_library(modules OBJECT)

target_sources(modules
  PRIVATE
    Module_p.cpp
    Module_s.cpp
    Module.cpp
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/Module_p.hpp
    ${CMAKE_CURRENT_LIST_DIR}/Module_s.hpp
    ${CMAKE_CURRENT_LIST_DIR}/Module.hpp
  )

target_include_directories(modules 
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  )




