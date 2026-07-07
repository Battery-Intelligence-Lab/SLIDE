# v4 core (slide::core) — grows alongside legacy src/ (strangler pattern, PLAN.md §5).
# Header-only for now; legacy targets do NOT link this.
add_library(slide_core INTERFACE)
target_include_directories(slide_core INTERFACE ${CMAKE_SOURCE_DIR}/src)
target_compile_features(slide_core INTERFACE cxx_std_20)
