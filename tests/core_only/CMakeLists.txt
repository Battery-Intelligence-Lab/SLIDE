if(SLIDE_WITH_CUDA OR SLIDE_WITH_MATLAB OR SLIDE_WITH_ZSTD OR SLIDE_WITH_ARROW)
  message(FATAL_ERROR "The core-only smoke gate requires every optional toolchain disabled")
endif()

add_executable(slide_core_only_smoke core_only_smoke.cpp)
target_link_libraries(slide_core_only_smoke
  PRIVATE
  slide_core
  project_options
  project_warnings)
add_test(NAME slide_core_only_smoke COMMAND slide_core_only_smoke)
