add_test(
  NAME structural_test_core_PC10SingleSource
  COMMAND "${CMAKE_COMMAND}"
    "-DSLIDE_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/pc10_single_source.cmake")
set_tests_properties(structural_test_core_PC10SingleSource
  PROPERTIES
    LABELS "structural;pc10"
    TIMEOUT 10)
