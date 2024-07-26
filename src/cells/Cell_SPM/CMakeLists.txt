cmake_minimum_required(VERSION 3.21)
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
  Model_SPM.hpp
)
target_include_directories(Cell_SPM PUBLIC .)
target_link_libraries(Cell_SPM PUBLIC Eigen3::Eigen)
 