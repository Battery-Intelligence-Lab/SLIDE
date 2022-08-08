cmake_minimum_required(VERSION 3.16)

set(PROJECT_NAME "SLIDE")

set(SLIDE_MAJOR_VERSION 3)
set(SLIDE_MINOR_VERSION 0)
set(SLIDE_PATCH_VERSION 0)
set(SLIDE_VERSION
  ${SLIDE_MAJOR_VERSION}.${SLIDE_MINOR_VERSION}.${SLIDE_PATCH_VERSION})

project(${PROJECT_NAME} VERSION ${SLIDE_VERSION} LANGUAGES CXX)
add_library(${PROJECT_NAME} "")

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
  EXPORT_COMPILE_COMMANDS ON)

# ensure C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # export compiler flags for code completion engines

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

# default build type.
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build." FORCE)

  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

add_definitions(-DSLIDE_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${SLIDE_ROOT_DIR}/exe/Debug)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${SLIDE_ROOT_DIR}/exe/Release)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${SLIDE_ROOT_DIR}/exe/RelWithDebInfo)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL ${SLIDE_ROOT_DIR}/exe/MinSizeRel)

# Set global compiler warnings
if(MSVC)
  message(STATUS "Building for MSVC")
  set(common_options /W2 /permissive-) # /WX
  set(debug_options ${common_options} /D_ITERATOR_DEBUG_LEVEL=0 /D_SCL_SECURE_NO_WARNINGS)
  set(not_debug_options ${common_options} /D_ITERATOR_DEBUG_LEVEL=0 /D_SCL_SECURE_NO_WARNINGS)

  add_compile_options(
    "$<$<CONFIG:Debug>:${debug_options}>$<$<NOT:$<CONFIG:Debug>>:${not_debug_options}>"
  )
else()
  message(STATUS "Building for non-MSVC")
  add_compile_options("-Wall"
    "$<$<CONFIG:Release>:-Ofast;-march=native>") # "-Weffc++" -Ofast -march=native -g -fno-omit-frame-pointer -gdwarf-2 (flto not good) -Wextra -pedantic
endif()

# guard against in-source builds -> from Eigen
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

message(STATUS "Host system: ${CMAKE_HOST_SYSTEM}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# set(slide_headers
# # determine_OCV.h
# # BasicCycler.hpp
# # cycling.h
# # degradation.h
# # cell_fit.hpp
# # cell_user.hpp
# # determine_characterisation.h
# settings/settings.hpp
# # currentWork.hpp
# # converter.hpp
# StorageUnit.hpp
# slide.hpp
# )

# set (slide_source
# # basic_cycler.cpp
# # degradation.cpp
# # cell_user.cpp
# # determine_characterisation.cpp
# # cell_fit.cpp
# # determine_OCV.cpp
# # cycling.cpp
# # currentWork.cpp
# # converter_test.cpp
# )
add_subdirectory(cells)
add_subdirectory(modules)
add_subdirectory(cooling)
add_subdirectory(power_conversion)
add_subdirectory(system)
add_subdirectory(procedures)
add_subdirectory(factories)

# add_subdirectory(unit_tests)
# add_subdirectory(settings)
target_sources(${PROJECT_NAME}
  PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/slide.hpp
)

target_link_libraries(${PROJECT_NAME}
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