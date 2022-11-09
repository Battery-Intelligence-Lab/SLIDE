cmake_minimum_required(VERSION 3.17)
add_library(Cell_SPM STATIC)

target_sources(Cell_SPM
  PRIVATE
  Cell_SPM.cpp
  Cell_SPM_degradation.cpp
  Cell_SPM_dstate.cpp
  Cell_SPM_thermal.cpp
  Cell_SPM_fitting.cpp
  Cell_SPM_diffusion.cpp
  PUBLIC
  Cell_SPM.hpp
)
target_include_directories(Cell_SPM PUBLIC .)
 