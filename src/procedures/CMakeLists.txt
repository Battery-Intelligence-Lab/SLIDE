cmake_minimum_required(VERSION 3.16)
add_library(procedures OBJECT)

target_sources(procedures
  PUBLIC
  Procedure.cpp
  Cycler.cpp

  # determine_OCV.cpp
  PUBLIC
  ${CMAKE_CURRENT_LIST_DIR}/Procedure.hpp
  ${CMAKE_CURRENT_LIST_DIR}/Cycler.hpp

  # ${CMAKE_CURRENT_LIST_DIR}/determine_OCV.hpp
)

target_include_directories(procedures
  PUBLIC
  ${CMAKE_CURRENT_LIST_DIR}
)
