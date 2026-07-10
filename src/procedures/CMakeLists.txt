add_library(procedures OBJECT)

target_sources(procedures
  PRIVATE
  Procedure.cpp
  Cycler.cpp
  determine_OCV.cpp
#  determine_characterisation.cpp
  #-------------------------
  PUBLIC
  Procedure.hpp
  Cycler.hpp
  determine_OCV.hpp
  determine_characterisation.hpp
)

target_include_directories(procedures PUBLIC .)
target_link_libraries(procedures PUBLIC types cells modules system)
