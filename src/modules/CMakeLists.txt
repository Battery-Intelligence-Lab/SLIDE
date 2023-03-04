cmake_minimum_required(VERSION 3.17)
add_library(modules STATIC)

target_sources(modules
  PRIVATE
    Module_p.cpp
    Module_s.cpp
    Module.cpp
  PUBLIC
    Module_p.hpp
    Module_s.hpp
    Module.hpp
  )

target_include_directories(modules PUBLIC .)
target_link_libraries(modules PRIVATE cooling)




