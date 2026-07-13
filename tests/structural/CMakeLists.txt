add_test(
  NAME structural_test_core_PC10SingleSource
  COMMAND "${CMAKE_COMMAND}"
    "-DSLIDE_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/pc10_single_source.cmake")
set_tests_properties(structural_test_core_PC10SingleSource
  PROPERTIES
    LABELS "structural;pc10"
    TIMEOUT 10)

add_test(
  NAME structural_test_core_9C2AgeingKernel
  COMMAND "${CMAKE_COMMAND}"
    "-DSLIDE_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/p9c_architecture.cmake")
set_tests_properties(structural_test_core_9C2AgeingKernel
  PROPERTIES
    LABELS "structural;9c2;9c3"
    TIMEOUT 10)

# Compiler-verified, so it needs the real compiler rather than a token scan.
add_test(
  NAME structural_test_core_9C5HeaderSelfContained
  COMMAND "${CMAKE_COMMAND}"
    "-DSLIDE_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
    "-DSLIDE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/p9c5_header_selfcontained.cmake")
set_tests_properties(structural_test_core_9C5HeaderSelfContained
  PROPERTIES
    LABELS "structural;9c5"
    TIMEOUT 300)
