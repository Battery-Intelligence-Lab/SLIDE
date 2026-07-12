include_guard(GLOBAL)

# Keep translation units whose correctness contract requires IEEE NaN/Inf
# generation or propagation outside the Release finite-math contract. Adding
# the strict option last preserves ordinary optimization at those boundaries.
function(slide_use_strict_fp)
  if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
    set(strict_fp_option "/fp:strict")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(strict_fp_option "-fno-fast-math")
  else()
    message(FATAL_ERROR
      "slide_use_strict_fp needs a strict floating-point option for ${CMAKE_CXX_COMPILER_ID}")
  endif()
  set_source_files_properties(${ARGN}
    PROPERTIES COMPILE_OPTIONS "${strict_fp_option}")
endfunction()

# Keep the dependency-light core source set in one place. Fuzzing builds a
# separately instrumented copy so sanitizer and Windows runtime/iterator ABI
# settings can never leak into the production slide_core target.
function(slide_add_core_library target)
  # Imported targets are directory-scoped. Resolve Threads where each core
  # target is created so sibling consumers such as tests/fuzz see the target.
  find_package(Threads REQUIRED)
  add_library(${target} STATIC
    "${PROJECT_SOURCE_DIR}/src/core/AsyncRecorder.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/CudaSpmBatch.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/Experiment.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/ForwardSensitivity.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/NetlistCsv.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/PackSolver.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/PackSolverValidation.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/PackStepper.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/PackTopology.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/ParameterSet.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/Recorder.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/SpmFactory.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/Simulation.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/ThreadPool.cpp")
  slide_use_strict_fp(
    "${PROJECT_SOURCE_DIR}/src/core/ForwardSensitivity.cpp"
    "${PROJECT_SOURCE_DIR}/src/core/PackSolverValidation.cpp")
  set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_include_directories(${target}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>)
  target_compile_features(${target} PUBLIC cxx_std_20)
  target_link_libraries(${target} PUBLIC Eigen3::Eigen Threads::Threads)
endfunction()
