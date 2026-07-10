if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  message(FATAL_ERROR "SLIDE_BUILD_FUZZERS requires Clang/libFuzzer")
endif()

set(SLIDE_FUZZ_RUNTIME_SANITIZERS "address,undefined")

include(${PROJECT_SOURCE_DIR}/cmake/SlideCoreTarget.cmake)
slide_add_core_library(slide_core_fuzz)

target_compile_options(slide_core_fuzz PRIVATE
  -fsanitize=fuzzer-no-link,${SLIDE_FUZZ_RUNTIME_SANITIZERS}
  -fno-omit-frame-pointer)
if(WIN32)
  # LLVM's prebuilt Windows libFuzzer runtime uses the release STL ABI even
  # when the instrumented target keeps Debug symbols and optimisation level.
  target_compile_definitions(slide_core_fuzz PRIVATE
    _ITERATOR_DEBUG_LEVEL=0 _HAS_ITERATOR_DEBUGGING=0)
  set_property(TARGET slide_core_fuzz PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded")
endif()
set_target_properties(slide_core_fuzz PROPERTIES
  CXX_SCAN_FOR_MODULES OFF
  INTERPROCEDURAL_OPTIMIZATION OFF)

function(add_slide_parser_fuzzer target source)
  add_executable(${target} ${source})
  target_link_libraries(${target} PRIVATE slide_core_fuzz)
  target_compile_options(${target} PRIVATE
    -fsanitize=fuzzer-no-link,${SLIDE_FUZZ_RUNTIME_SANITIZERS}
    -fno-omit-frame-pointer)
  if(WIN32)
    target_compile_definitions(${target} PRIVATE
      _ITERATOR_DEBUG_LEVEL=0 _HAS_ITERATOR_DEBUGGING=0)
    set_property(TARGET ${target} PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded")
  endif()
  target_link_options(${target} PRIVATE
    -fsanitize=fuzzer,${SLIDE_FUZZ_RUNTIME_SANITIZERS}
    -fno-omit-frame-pointer)
  set_target_properties(${target} PROPERTIES
    CXX_SCAN_FOR_MODULES OFF
    INTERPROCEDURAL_OPTIMIZATION OFF
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/fuzz"
    RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/fuzz"
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/fuzz")
endfunction()

add_slide_parser_fuzzer(fuzz_experiment fuzz_Experiment.cpp)
add_slide_parser_fuzzer(fuzz_bpx fuzz_Bpx.cpp)
add_slide_parser_fuzzer(fuzz_netlist_csv fuzz_NetlistCsv.cpp)
