cmake_minimum_required(VERSION 3.21)
add_library(modules OBJECT)

target_sources(modules
  PUBLIC
    Module_p_impl.cpp
    Module_p.cpp
    Module_s.cpp
    Module.cpp
    Module_p.hpp
    Module_s.hpp
    Module.hpp
  )

target_include_directories(modules PUBLIC .)
target_link_libraries(modules PUBLIC cooling cells settings Eigen3::Eigen Boost::boost)