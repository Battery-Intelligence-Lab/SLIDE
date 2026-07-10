cmake_minimum_required(VERSION 3.31)
project(SLIDEExternalConsumer LANGUAGES CXX)

if(NOT DEFINED SLIDE_SOURCE_DIR)
  message(FATAL_ERROR "Configure with -DSLIDE_SOURCE_DIR=<SLIDE checkout>")
endif()

set(SLIDE_CORE_ONLY ON CACHE BOOL "" FORCE)
set(SLIDE_BUILD_CORE_TESTS OFF CACHE BOOL "" FORCE)
set(SLIDE_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
set(SLIDE_WITH_CUDA OFF CACHE BOOL "" FORCE)
set(SLIDE_WITH_MATLAB OFF CACHE BOOL "" FORCE)
set(SLIDE_WITH_ZSTD OFF CACHE BOOL "" FORCE)
set(SLIDE_WITH_ARROW OFF CACHE BOOL "" FORCE)

add_subdirectory("${SLIDE_SOURCE_DIR}" slide)

add_executable(slide_external_consumer ../core_only_smoke.cpp)
target_link_libraries(slide_external_consumer PRIVATE slide_core)

enable_testing()
add_test(NAME slide_external_consumer COMMAND slide_external_consumer)
