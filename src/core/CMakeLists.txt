# v4 core (slide::core) — grows alongside legacy src/ (strangler pattern, PLAN.md §5).
# The registry is compiled once so downstream callers do not instantiate every common composition.
add_library(slide_core STATIC PackSolver.cpp PackTopology.cpp SpmFactory.cpp Simulation.cpp)
target_include_directories(slide_core PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_compile_features(slide_core PUBLIC cxx_std_20)
target_link_libraries(slide_core PUBLIC Eigen3::Eigen)
