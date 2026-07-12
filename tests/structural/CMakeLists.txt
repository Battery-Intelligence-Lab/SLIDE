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
