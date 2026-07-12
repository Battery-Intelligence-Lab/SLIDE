include_guard(GLOBAL)

option(SLIDE_ENABLE_COVERAGE
  "Enable the dedicated Clang source-based coverage lane" OFF)

function(_slide_require_matching_llvm_tool executable label expected_major)
  execute_process(
    COMMAND "${executable}" --version
    RESULT_VARIABLE tool_result
    OUTPUT_VARIABLE tool_stdout
    ERROR_VARIABLE tool_stderr
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE)
  if(NOT tool_result EQUAL 0)
    message(FATAL_ERROR
      "${label} failed its version probe (${tool_result}): ${tool_stderr}")
  endif()

  string(TOLOWER "${tool_stdout}\n${tool_stderr}" tool_version)
  string(REGEX MATCH "version[ \t]+([0-9]+)" _ "${tool_version}")
  if(NOT CMAKE_MATCH_1)
    message(FATAL_ERROR
      "Could not determine ${label}'s major version from: ${tool_stdout} ${tool_stderr}")
  endif()
  if(NOT CMAKE_MATCH_1 STREQUAL expected_major)
    message(FATAL_ERROR
      "Coverage profile formats require matching LLVM tools: Clang is major "
      "${expected_major}, but ${label} is major ${CMAKE_MATCH_1} (${executable})")
  endif()
endfunction()

# Create the reusable coverage policy before any production or test targets are
# created. The policy is attached to slide_core and every registered test TU;
# legacy physics is intentionally outside P9-G3's src/core failure-site scope.
function(slide_enable_coverage)
  if(NOT SLIDE_ENABLE_COVERAGE)
    return()
  endif()
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(FATAL_ERROR
      "SLIDE_ENABLE_COVERAGE requires Clang source-based coverage; found "
      "${CMAKE_CXX_COMPILER_ID}")
  endif()
  if(CMAKE_CONFIGURATION_TYPES)
    message(FATAL_ERROR
      "SLIDE_ENABLE_COVERAGE requires a single-config Debug build (use Ninja)")
  endif()
  if(NOT CMAKE_GENERATOR MATCHES "^Ninja")
    message(FATAL_ERROR
      "SLIDE_ENABLE_COVERAGE requires Ninja for reproducible freshness checks")
  endif()
  if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(FATAL_ERROR
      "SLIDE_ENABLE_COVERAGE requires -DCMAKE_BUILD_TYPE=Debug")
  endif()
  if(ENABLE_IPO)
    message(FATAL_ERROR
      "Coverage builds require -DENABLE_IPO=OFF so mappings stay inspectable")
  endif()
  if(SLIDE_CORE_ONLY)
    message(FATAL_ERROR
      "SLIDE_ENABLE_COVERAGE requires the full test graph; SLIDE_CORE_ONLY is unsupported")
  endif()
  if(SLIDE_BUILD_CORE_TESTS OR SLIDE_BUILD_FUZZERS)
    message(FATAL_ERROR
      "SLIDE_ENABLE_COVERAGE cannot include unregistered core-smoke or fuzz targets")
  endif()
  if(SLIDE_WITH_CUDA OR SLIDE_WITH_ARROW OR SLIDE_WITH_ZSTD)
    message(FATAL_ERROR
      "P9-G3 coverage uses the optional-off core configuration")
  endif()

  set(enabled_sanitizers)
  foreach(sanitizer IN ITEMS
      ADDRESS LEAK UNDEFINED_BEHAVIOR THREAD MEMORY)
    if(ENABLE_SANITIZER_${sanitizer})
      list(APPEND enabled_sanitizers "${sanitizer}")
    endif()
  endforeach()
  if(enabled_sanitizers)
    list(JOIN enabled_sanitizers ", " sanitizer_list)
    message(FATAL_ERROR
      "Coverage must use its own build tree; enabled sanitizers: ${sanitizer_list}")
  endif()

  string(REGEX MATCH "^[0-9]+" clang_major "${CMAKE_CXX_COMPILER_VERSION}")
  if(NOT clang_major)
    message(FATAL_ERROR
      "Could not determine Clang's major version from ${CMAKE_CXX_COMPILER_VERSION}")
  endif()
  get_filename_component(clang_binary_directory "${CMAKE_CXX_COMPILER}" DIRECTORY)
  find_program(SLIDE_LLVM_PROFDATA_EXECUTABLE
    NAMES "llvm-profdata-${clang_major}" llvm-profdata
    HINTS "${clang_binary_directory}"
    REQUIRED)
  find_program(SLIDE_LLVM_COV_EXECUTABLE
    NAMES "llvm-cov-${clang_major}" llvm-cov
    HINTS "${clang_binary_directory}"
    REQUIRED)
  _slide_require_matching_llvm_tool(
    "${SLIDE_LLVM_PROFDATA_EXECUTABLE}" llvm-profdata "${clang_major}")
  _slide_require_matching_llvm_tool(
    "${SLIDE_LLVM_COV_EXECUTABLE}" llvm-cov "${clang_major}")
  mark_as_advanced(
    SLIDE_LLVM_PROFDATA_EXECUTABLE
    SLIDE_LLVM_COV_EXECUTABLE)

  set(SLIDE_COVERAGE_OUTPUT_DIRECTORY
      "${CMAKE_BINARY_DIR}/coverage" CACHE INTERNAL
      "Build-local SLIDE coverage artifact directory" FORCE)
  set(SLIDE_COVERAGE_RAW_DIRECTORY
      "${CMAKE_BINARY_DIR}/coverage/raw" CACHE INTERNAL
      "Build-local SLIDE raw-profile directory" FORCE)
  set(SLIDE_COVERAGE_TARGET_MANIFEST
      "${CMAKE_BINARY_DIR}/coverage/targets.tsv" CACHE INTERNAL
      "Generated target-to-binary coverage manifest" FORCE)

  set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL
      "Export exact coverage preprocessing commands" FORCE)
  add_library(slide_coverage_options INTERFACE)
  target_compile_options(slide_coverage_options INTERFACE
    -fprofile-instr-generate
    -fcoverage-mapping
    -fprofile-update=atomic
    -O0
    -g)
  target_link_options(slide_coverage_options INTERFACE
    -fprofile-instr-generate)
  if(UNIX AND NOT APPLE)
    find_program(SLIDE_COVERAGE_GOLD_LINKER NAMES ld.gold)
    if(SLIDE_COVERAGE_GOLD_LINKER)
      # Gold handles the many mapping-heavy standalone test links materially
      # faster than BFD on DrvFS and accepts the same ELF coverage sections.
      target_link_options(slide_coverage_options INTERFACE -fuse-ld=gold)
      mark_as_advanced(SLIDE_COVERAGE_GOLD_LINKER)
    endif()
  endif()

  message(STATUS
    "SLIDE Clang source coverage: ${SLIDE_LLVM_PROFDATA_EXECUTABLE}; "
    "${SLIDE_LLVM_COV_EXECUTABLE}; linker ${SLIDE_COVERAGE_GOLD_LINKER}")
endfunction()

function(slide_instrument_coverage_target target_name)
  if(NOT SLIDE_ENABLE_COVERAGE)
    return()
  endif()
  if(NOT TARGET ${target_name})
    message(FATAL_ERROR "Coverage target does not exist: ${target_name}")
  endif()
  target_link_libraries(${target_name} PRIVATE slide_coverage_options)
endfunction()

# All existing test declarations pass through this helper. In a coverage lane,
# every process gets a target-, PID-, and binary-signature-specific raw profile
# under the current build tree, so parallel CTest cannot overwrite evidence.
function(add_executable_with_coverage_and_test target_name source_name)
  add_executable(${target_name} ${source_name})
  target_link_libraries(${target_name} PRIVATE src Catch2::Catch2WithMain)
  add_test(NAME ${target_name} COMMAND ${target_name} WORKING_DIRECTORY bin)

  if(SLIDE_ENABLE_COVERAGE)
    target_link_libraries(${target_name} PRIVATE slide_coverage_options)
    set_property(GLOBAL APPEND PROPERTY
      SLIDE_COVERAGE_TEST_TARGETS "${target_name}")
    set_tests_properties(${target_name} PROPERTIES
      ENVIRONMENT
        "LLVM_PROFILE_FILE=${SLIDE_COVERAGE_RAW_DIRECTORY}/${target_name}-%p-%m.profraw")
  endif()
endfunction()

# llvm-cov cannot read an ar archive directly. This build-only executable forces
# every slide_core object into one valid coverage-mapping container. Test
# executables remain separate objects in the manifest so their header/template
# instantiations and execution counts are retained by the reporter.
function(slide_finalize_coverage production_target)
  if(NOT SLIDE_ENABLE_COVERAGE)
    return()
  endif()
  if(NOT TARGET ${production_target})
    message(FATAL_ERROR
      "Coverage production target does not exist: ${production_target}")
  endif()

  get_property(test_targets GLOBAL PROPERTY SLIDE_COVERAGE_TEST_TARGETS)
  if(NOT test_targets)
    message(FATAL_ERROR "Coverage enabled, but no test targets were registered")
  endif()
  list(REMOVE_DUPLICATES test_targets)

  add_executable(slide_coverage_anchor
    "${PROJECT_SOURCE_DIR}/tests/coverage/coverage_anchor.cpp")
  target_link_libraries(slide_coverage_anchor PRIVATE
    slide_coverage_options
    "$<LINK_LIBRARY:WHOLE_ARCHIVE,${production_target}>")
  # The anchor's only job is retaining production mappings. Suppress coverage
  # for its trivial main so the reporter need not execute it merely to avoid a
  # mismatched/missing function record.
  target_compile_options(slide_coverage_anchor PRIVATE
    -fno-profile-instr-generate
    -fno-coverage-mapping)

  # Keep benchmarks and user-facing executables outside the evidence build.
  # They are neither CTest producers nor objects consumed by the reporter.
  add_custom_target(slide_coverage_binaries
    DEPENDS slide_coverage_anchor ${test_targets})

  file(MAKE_DIRECTORY "${SLIDE_COVERAGE_OUTPUT_DIRECTORY}")
  set(manifest
    "anchor\tslide_coverage_anchor\t$<TARGET_FILE:slide_coverage_anchor>\n")
  foreach(test_target IN LISTS test_targets)
    string(APPEND manifest
      "test\t${test_target}\t$<TARGET_FILE:${test_target}>\n")
  endforeach()
  file(GENERATE
    OUTPUT "${SLIDE_COVERAGE_TARGET_MANIFEST}"
    CONTENT "${manifest}")

  list(LENGTH test_targets test_target_count)
  message(STATUS
    "SLIDE coverage manifest: 1 production anchor + ${test_target_count} test binaries")
endfunction()
