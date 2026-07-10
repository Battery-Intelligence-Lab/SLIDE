include_guard(GLOBAL)

# Keep the dependency-light core source set in one place. Fuzzing builds a
# separately instrumented copy so sanitizer and Windows runtime/iterator ABI
# settings can never leak into the production slide_core target.
function(slide_add_core_library target)
  add_library(${target} STATIC
    "${PROJECT_SOURCE_DIR}/src/core/AsyncRecorder.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/CudaSpmBatch.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/Experiment.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/ForwardSensitivity.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/NetlistCsv.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/PackSolver.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/PackStepper.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/PackTopology.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/ParameterSet.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/Recorder.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/SpmFactory.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/Simulation.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/ThreadPool.cpp")
  set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_include_directories(${target}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>)
  target_compile_features(${target} PUBLIC cxx_std_20)
  target_link_libraries(${target} PUBLIC Eigen3::Eigen)
endfunction()
