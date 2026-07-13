# M0.9 / 9C-5: every api and support header compiles on its own.
#
# Adversarial review found `SeiParams.hpp` and `SurfaceCrackParams.hpp` using `slide::Status`
# without including `../types/Status.hpp`: they compiled only because `SpmFactory.hpp` happened
# to include a sibling that pulled Status in first. An include reorder, an IWYU pass, or a user
# writing `#include <core/SeiParams.hpp>` as their first include would have broken the build.
# A token gate cannot see that. The compiler can, so this gate asks it.

if(NOT DEFINED SLIDE_SOURCE_DIR OR NOT DEFINED SLIDE_CXX_COMPILER)
  message(FATAL_ERROR "SLIDE_SOURCE_DIR and SLIDE_CXX_COMPILER are required")
endif()

set(P9C5_PUBLIC_HEADERS
  AsyncRecorder CellDesign CudaSpmBatch EulerLegacy Experiment ExponentialModal
  ForwardSensitivity NetlistCsv PackSolver PackStepper PackTopology ParameterSet
  Recorder Simulation SpmFactory ThreadPool
  AgeingModelMask BatchBuilder BatchView CompiledCurve LamParams LithiumPlatingParams
  Numeric SeiParams SpmBatchLayout SpmScalarKernels StateArena SurfaceCrackParams)

set(scratch "${CMAKE_CURRENT_BINARY_DIR}/p9c5_selfcontained")
file(MAKE_DIRECTORY "${scratch}")

set(broken "")
foreach(header IN LISTS P9C5_PUBLIC_HEADERS)
  set(source "${scratch}/${header}_alone.cpp")
  file(WRITE "${source}" "#include <core/${header}.hpp>\nint main() { return 0; }\n")
  execute_process(
    COMMAND "${SLIDE_CXX_COMPILER}" -std=c++20 "-I${SLIDE_SOURCE_DIR}/src"
            -fsyntax-only "${source}"
    RESULT_VARIABLE status
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err)
  if(NOT status EQUAL 0)
    string(REGEX MATCH "error:[^\n]*" first_error "${err}")
    list(APPEND broken "${header}.hpp (${first_error})")
  endif()
endforeach()

if(NOT broken STREQUAL "")
  string(REPLACE ";" "\n  " report "${broken}")
  message(FATAL_ERROR
    "9C-5 self-containment: these public headers do not compile on their own:\n  ${report}\n"
    "A header on the public surface must not depend on the include order of whoever reaches it first.")
endif()

message(STATUS "9C-5 header self-containment gate passed (28 headers)")
