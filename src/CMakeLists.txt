cmake_minimum_required(VERSION 3.20)

add_library(src STATIC)

#add_definitions(-DSLIDE_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../")

# Set global compiler warnings
# if(MSVC)
#   message(STATUS "Building for MSVC")
#   set(common_options /W2 /permissive-) # /WX
#   set(debug_options ${common_options} /D_ITERATOR_DEBUG_LEVEL=0 /D_SCL_SECURE_NO_WARNINGS)
#   set(not_debug_options ${common_options} /D_ITERATOR_DEBUG_LEVEL=0 /D_SCL_SECURE_NO_WARNINGS)

#   add_compile_options( # REMOVE THIS #CHECK
#     "$<$<CONFIG:Debug>:${debug_options}>$<$<NOT:$<CONFIG:Debug>>:${not_debug_options}>"
#   )
# else()
#   message(STATUS "Building for non-MSVC")
#   add_compile_options("-Wall"
#     "$<$<CONFIG:Release>:-Ofast;-march=native>") # "-Weffc++" -Ofast -march=native -g -fno-omit-frame-pointer -gdwarf-2 (flto not good) -Wextra -pedantic
# endif()


add_subdirectory(cells)
add_subdirectory(modules)
add_subdirectory(cooling)
add_subdirectory(power_conversion)
add_subdirectory(system)
add_subdirectory(procedures)
add_subdirectory(factories)

# add_subdirectory(unit_tests)
# add_subdirectory(settings)
# target_sources(src
#   PRIVATE
#   slide.hpp
# )

target_link_libraries(src
  PRIVATE
  cells
  modules
  cooling
  power_conversion
  system
  procedures
  factories

  # unit_tests
  # pthread -> For WSL.
)