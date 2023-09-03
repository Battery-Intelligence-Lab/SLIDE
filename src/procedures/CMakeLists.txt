add_library(procedures STATIC)

target_sources(procedures
  PUBLIC
  Procedure.cpp
  Cycler.cpp
  determine_OCV.cpp
  determine_characterisation.cpp
  #-------------------------
  Procedure.hpp
  Cycler.hpp
  determine_OCV.hpp
  determine_characterisation.hpp
)

target_include_directories(procedures PUBLIC .)
target_link_libraries(procedures PUBLIC types)