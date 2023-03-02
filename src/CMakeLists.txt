cmake_minimum_required(VERSION 3.17)

add_library(src STATIC
            paperCode.cpp)

#add_definitions(-DSLIDE_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../")

add_subdirectory(cells)
add_subdirectory(modules)
add_subdirectory(cooling)
add_subdirectory(power_conversion)
add_subdirectory(system)
add_subdirectory(procedures)
add_subdirectory(factories)
add_subdirectory(settings)

target_link_libraries(src
  PRIVATE
  cells
  modules
  cooling
  power_conversion
  system
  procedures
  factories
  settings
 # pthread #-> For WSL.
)