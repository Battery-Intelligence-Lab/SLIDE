add_library(procedures STATIC)

target_sources(procedures
  PUBLIC
  Procedure.cpp
  Cycler.cpp
  determine_OCV.cpp
  PUBLIC
  Procedure.hpp
  Cycler.hpp
  determine_OCV.hpp
)

target_include_directories(procedures PUBLIC .)
