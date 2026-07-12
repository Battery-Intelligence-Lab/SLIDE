# M0.7 / 9C-3: one concept per cold translation unit, with one shared
# implementation for reusable parameter curves and Experiment semantics.

if(NOT DEFINED SLIDE_SOURCE_DIR)
  message(FATAL_ERROR "SLIDE_SOURCE_DIR is required")
endif()

function(p9c3_load_compact relative_path output)
  set(path "${SLIDE_SOURCE_DIR}/${relative_path}")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "9C-3 source is missing: ${relative_path}")
  endif()
  file(READ "${path}" content)
  string(REGEX REPLACE "//[^\r\n]*" "" content "${content}")
  string(REGEX REPLACE "/\\*([^*]|\\*+[^*/])*\\*+/" "" content "${content}")
  string(REGEX REPLACE "[ \t\r\n]" "" content "${content}")
  set(${output} "${content}" PARENT_SCOPE)
endfunction()

function(p9c3_require_tokens label content_variable)
  foreach(token IN ITEMS ${ARGN})
    string(FIND "${${content_variable}}" "${token}" position)
    if(position EQUAL -1)
      message(FATAL_ERROR
        "9C-3 ${label}: required token is absent: ${token}")
    endif()
  endforeach()
endfunction()

function(p9c3_forbid_tokens label content_variable)
  foreach(token IN ITEMS ${ARGN})
    string(FIND "${${content_variable}}" "${token}" position)
    if(NOT position EQUAL -1)
      message(FATAL_ERROR
        "9C-3 ${label}: forbidden ownership/dependency token remains: ${token}")
    endif()
  endforeach()
endfunction()

function(p9c3_require_token_count label content_variable token expected)
  string(LENGTH "${${content_variable}}" original_length)
  string(REPLACE "${token}" "" without_token "${${content_variable}}")
  string(LENGTH "${without_token}" stripped_length)
  string(LENGTH "${token}" token_length)
  if(token_length EQUAL 0)
    message(FATAL_ERROR "9C-3 internal error: empty count token for ${label}")
  endif()
  math(EXPR count
    "(${original_length} - ${stripped_length}) / ${token_length}")
  if(NOT count EQUAL expected)
    message(FATAL_ERROR
      "9C-3 ${label}: expected ${expected} occurrences of ${token}, found ${count}")
  endif()
endfunction()

function(p9c3_require_max_lines relative_path maximum)
  set(path "${SLIDE_SOURCE_DIR}/${relative_path}")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "9C-3 source is missing: ${relative_path}")
  endif()
  file(READ "${path}" content)
  string(REGEX REPLACE "[^\n]" "" newline_characters "${content}")
  string(LENGTH "${newline_characters}" line_count)
  string(LENGTH "${content}" content_length)
  if(content_length GREATER 0)
    math(EXPR last_index "${content_length} - 1")
    string(SUBSTRING "${content}" ${last_index} 1 last_character)
    if(NOT last_character STREQUAL "\n")
      math(EXPR line_count "${line_count} + 1")
    endif()
  endif()
  if(line_count GREATER maximum)
    message(FATAL_ERROR
      "9C-3 ${relative_path}: ${line_count} physical lines exceeds ceiling ${maximum}")
  endif()
endfunction()

# Production concepts and their deliberately non-public sharing seams.
p9c3_load_compact("src/core/ParameterSet.cpp" parameter_set)
p9c3_load_compact("src/core/BpxParameterReader.cpp" bpx_reader)
p9c3_load_compact("src/core/BpxExpression.cpp" bpx_expression)
p9c3_load_compact("src/core/StrictJson.cpp" strict_json)
p9c3_load_compact("src/core/detail/ParameterCurve.hpp" parameter_curve)
p9c3_load_compact("src/core/detail/BpxExpression.hpp" bpx_expression_interface)
p9c3_load_compact("src/core/detail/StrictJson.hpp" strict_json_interface)
p9c3_load_compact("src/core/Experiment.cpp" experiment_parser)
p9c3_load_compact("src/core/CyclerV2.cpp" experiment_runner)
p9c3_load_compact("src/core/detail/ExperimentSemantics.hpp" experiment_semantics)
p9c3_load_compact("src/core/ParameterSet.hpp" parameter_set_public)
p9c3_load_compact("src/core/Experiment.hpp" experiment_public)
p9c3_load_compact("cmake/SlideCoreTarget.cmake" core_target)
p9c3_load_compact("tests/unit/CMakeLists.txt" unit_cmake)
p9c3_load_compact("tests/unit/core_Experiment_test.cpp" experiment_test)

# The centralized list feeds production, fuzz, coverage, core-only, and nested
# consumers. Each moved implementation must occur exactly once there.
foreach(source IN ITEMS
    BpxExpression.cpp
    BpxParameterReader.cpp
    CyclerV2.cpp
    Experiment.cpp
    ParameterSet.cpp
    StrictJson.cpp)
  set(source_token "\"\${PROJECT_SOURCE_DIR}/src/core/${source}\"")
  p9c3_require_token_count(
    "central core source list (${source})" core_target "${source_token}" 1)
endforeach()

# The value-semantic ParameterSet owns absorption and SPM mapping, not wire
# formats, scanners, or file I/O.
p9c3_require_tokens("ParameterSet absorption" parameter_set
  "#include\"ParameterSet.hpp\""
  "#include\"detail/ParameterCurve.hpp\""
  "graphiteOcp("
  "nmcOcp("
  "requiredScalar("
  "ParameterSet::chen2020("
  "ParameterSet::toSpmInput(")
foreach(method IN ITEMS
    canonicalName set update contains find findScalar findCurve describe
    chen2020 toSpmInput)
  p9c3_require_token_count(
    "ParameterSet method ownership (${method})"
    parameter_set "ParameterSet::${method}(" 1)
endforeach()
p9c3_require_token_count(
  "ParameterSet excludes BPX JSON" parameter_set
  "ParameterSet::fromBpxJson(" 0)
p9c3_require_token_count(
  "ParameterSet excludes BPX files" parameter_set
  "ParameterSet::fromBpxFile(" 0)
p9c3_require_tokens("ParameterSet allocation Status boundary" parameter_set
  "[[nodiscard]]slide::StatusallocationFailureStatus()noexcept{returnslide::Status::Numerical_failure;}")
p9c3_require_token_count(
  "ParameterSet allocation mapper ownership"
  parameter_set "allocationFailureStatus()" 8)
p9c3_require_token_count(
  "ParameterSet bad_alloc boundary"
  parameter_set "catch(conststd::bad_alloc&)" 4)
p9c3_require_token_count(
  "ParameterSet length_error boundary"
  parameter_set "catch(conststd::length_error&)" 3)
p9c3_require_token_count(
  "ParameterSet allocation delegation"
  parameter_set "returnallocationFailureStatus();" 7)
p9c3_forbid_tokens("ParameterSet dependency boundary" parameter_set
  "BoundedFileReader.hpp"
  "BpxExpression.hpp"
  "StrictJson.hpp"
  "JsonParser"
  "max_bpx_json_bytes")

# The BPX adapter is the sole owner of wire semantics and atomic publication.
p9c3_require_tokens("BPX reader" bpx_reader
  "#include\"ParameterSet.hpp\""
  "#include\"BoundedFileReader.hpp\""
  "#include\"detail/BpxExpression.hpp\""
  "#include\"detail/ParameterCurve.hpp\""
  "#include\"detail/StrictJson.hpp\""
  "max_bpx_json_bytes"
  "jsonPath("
  "jsonConstant("
  "jsonCurve("
  "parseStrictJson("
  "ParameterSet::fromBpxJson("
  "ParameterSet::fromBpxFile(")
p9c3_require_token_count(
  "BPX JSON API ownership" bpx_reader "ParameterSet::fromBpxJson(" 1)
p9c3_require_token_count(
  "BPX file API ownership" bpx_reader "ParameterSet::fromBpxFile(" 1)
foreach(seam IN ITEMS
    evaluateBpxExpressionSamples
    sampleBpxExpressionCurve
    parseStrictJson)
  p9c3_require_token_count(
    "BPX reader consumes internal seam (${seam})"
    bpx_reader "detail::${seam}(" 1)
endforeach()
foreach(method IN ITEMS
    canonicalName set update contains find findScalar findCurve describe
    chen2020 toSpmInput)
  p9c3_require_token_count(
    "BPX reader excludes ParameterSet core method (${method})"
    bpx_reader "ParameterSet::${method}(" 0)
endforeach()
p9c3_forbid_tokens("BPX reader excludes private parsers" bpx_reader
  "classBpxExpression"
  "classJsonParser")

# Expression parsing/evaluation and strict JSON parsing stay independent of
# ParameterSet and of one another. Their headers are narrow internal seams.
p9c3_require_tokens("BPX expression implementation" bpx_expression
  "#include\"detail/BpxExpression.hpp\""
  "evaluateBpxExpressionSamples("
  "sampleBpxExpressionCurve(")
p9c3_require_tokens("BPX expression interface" bpx_expression_interface
  "evaluateBpxExpressionSamples("
  "sampleBpxExpressionCurve(")
foreach(seam IN ITEMS
    evaluateBpxExpressionSamples
    sampleBpxExpressionCurve)
  p9c3_require_token_count(
    "BPX expression definition (${seam})" bpx_expression "${seam}(" 1)
  p9c3_require_token_count(
    "BPX expression declaration (${seam})"
    bpx_expression_interface "${seam}(" 1)
endforeach()
p9c3_forbid_tokens("BPX expression DAG" bpx_expression
  "ParameterSet"
  "StrictJson"
  "JsonParser"
  "BoundedFileReader")
p9c3_forbid_tokens("BPX expression interface DAG" bpx_expression_interface
  "ParameterSet"
  "StrictJson"
  "BoundedFileReader")

p9c3_require_tokens("strict JSON implementation" strict_json
  "#include\"detail/StrictJson.hpp\""
  "JsonParser"
  "max_json_values"
  "parseStrictJson("
  "StrictJsonValuecandidate;"
  "static_assert(std::is_nothrow_move_assignable_v<StrictJsonValue>);"
  "output=std::move(candidate);")
p9c3_require_token_count(
  "strict JSON parser boundary" strict_json "parseStrictJson(" 1)
p9c3_forbid_tokens("strict JSON atomic publication" strict_json
  "parser.parse(output)")
p9c3_forbid_tokens("strict JSON DAG" strict_json
  "ParameterSet"
  "BpxExpression"
  "BoundedFileReader"
  "ParameterCurve")
p9c3_require_tokens("strict JSON interface" strict_json_interface
  "StrictJsonValue"
  "parseStrictJson(")
p9c3_require_token_count(
  "strict JSON declaration" strict_json_interface "parseStrictJson(" 1)
p9c3_forbid_tokens("strict JSON interface DAG" strict_json_interface
  "ParameterSet"
  "BpxExpression"
  "BoundedFileReader")

p9c3_require_tokens("shared parameter curve" parameter_curve
  "validParameterCurve("
  "sampleParameterCurve("
  "template<")
p9c3_forbid_tokens("shared parameter curve DAG" parameter_curve
  "ParameterSet"
  "BpxExpression"
  "StrictJson"
  "BoundedFileReader")

# Parser and runner share normalization and semantic admission from one inline
# definition, while retaining distinct public-method ownership.
p9c3_require_tokens("shared Experiment semantics" experiment_semantics
  "inlinestd::stringnormalized("
  "inlineboolnormalizedEqual("
  "inlineboolvalidDirection("
  "inlineboolvalidSegment(")
p9c3_forbid_tokens("shared Experiment semantics boundary" experiment_semantics
  "CyclerV2::"
  "Experiment::parse("
  "ParameterSet"
  "Bpx")

p9c3_require_tokens("Experiment parser" experiment_parser
  "#include\"Experiment.hpp\""
  "#include\"detail/ExperimentSemantics.hpp\""
  "detail::normalized("
  "detail::validSegment("
  "Experiment::parse(")
p9c3_require_token_count(
  "Experiment parser ownership" experiment_parser "Experiment::parse(" 1)
p9c3_require_token_count(
  "Experiment parser excludes runner" experiment_parser "CyclerV2::" 0)
p9c3_forbid_tokens("Experiment parser has no private semantic copy"
  experiment_parser
  "std::stringnormalized("
  "boolvalidDirection("
  "boolvalidSegment(")

p9c3_require_tokens("Cycler runner" experiment_runner
  "#include\"Experiment.hpp\""
  "#include\"detail/ExperimentSemantics.hpp\""
  "detail::normalizedEqual("
  "detail::validSegment(")
p9c3_require_token_count(
  "Cycler runner excludes parser" experiment_runner "Experiment::parse(" 0)
foreach(method IN ITEMS
    configure
    registerDriveCycle
    voltageAt
    currentForVoltage
    currentForPower
    evaluateFunction
    currentForCustom
    advance
    findDriveCycle
    driveCurrent
    run)
  p9c3_require_token_count(
    "CyclerV2 method ownership (${method})"
    experiment_runner "CyclerV2::${method}(" 1)
endforeach()
p9c3_require_token_count(
  "CyclerV2 exact method count" experiment_runner "CyclerV2::" 11)
p9c3_require_tokens("Cycler allocation Status boundary" experiment_runner
  "[[nodiscard]]slide::StatusallocationFailureStatus()noexcept{returnslide::Status::Numerical_failure;}")
p9c3_require_token_count(
  "Cycler allocation mapper ownership"
  experiment_runner "allocationFailureStatus()" 7)
p9c3_require_token_count(
  "Cycler bad_alloc boundaries"
  experiment_runner "catch(conststd::bad_alloc&)" 4)
p9c3_require_token_count(
  "Cycler length_error boundaries"
  experiment_runner "catch(conststd::length_error&)" 4)
p9c3_require_token_count(
  "Cycler allocation delegation"
  experiment_runner "returnallocationFailureStatus();" 6)
p9c3_require_token_count(
  "Cycler callback allocation rethrows"
  experiment_runner "throw;" 2)
p9c3_require_token_count(
  "Cycler snapshot restore ownership"
  experiment_runner "restoreBatchSnapshot(" 3)
p9c3_require_token_count(
  "Cycler run rollback guards"
  experiment_runner "run_snapshot_ready_&&batch_!=nullptr" 2)
p9c3_forbid_tokens("Cycler runner has no private semantic copy"
  experiment_runner
  "std::stringnormalized("
  "boolvalidDirection("
  "boolvalidSegment(")

p9c3_require_tokens("Cycler stale-scratch coverage witness" experiment_test
  "[core][experiment][configuration][coverage]"
  "core::SpmModelOptions{.nch=8}"
  "replacement.state().raw().size()!=batch.state().raw().size()"
  "cycler.run(experiment,1.0,output)==Status::Numerical_failure"
  "std::ranges::equal(batch.state().raw(),state_before)"
  "std::ranges::equal(batch.derivative().raw(),derivative_before)")

# Internal sharing must not leak into the two public description headers or
# any other top-level core header.
p9c3_forbid_tokens("ParameterSet public header" parameter_set_public "detail/")
p9c3_forbid_tokens("Experiment public header" experiment_public "detail/")
file(GLOB p9c3_public_headers "${SLIDE_SOURCE_DIR}/src/core/*.hpp")
foreach(public_header_path IN LISTS p9c3_public_headers)
  file(RELATIVE_PATH public_header_relative
    "${SLIDE_SOURCE_DIR}" "${public_header_path}")
  p9c3_load_compact("${public_header_relative}" public_header)
  p9c3_forbid_tokens(
    "public header isolation (${public_header_relative})" public_header
    "#include\"detail/"
    "#include<detail/")
endforeach()

# Exact fixtures remain in existing test executables, and configuration rather
# than __FAST_MATH__ identifies Release because the Experiment test is strict-FP.
p9c3_require_token_count(
  "ParameterSet recorded fixture wiring" unit_cmake
  "slide_configure_recorded_scalar_fixture(unit_test_core_ParameterSet)" 1)
p9c3_require_token_count(
  "Experiment recorded fixture wiring" unit_cmake
  "slide_configure_recorded_scalar_fixture(unit_test_core_Experiment)" 1)
p9c3_require_tokens("recorded fixture configuration" unit_cmake
  "\"$<$<CONFIG:Release>:SLIDE_TEST_RELEASE=1>\""
  "SLIDE_TEST_IPO=1")

# Physical-line ceilings implement MC-1 without allowing whitespace removal to
# disguise an oversized concept.
p9c3_require_max_lines("src/core/ParameterSet.cpp" 450)
p9c3_require_max_lines("src/core/BpxParameterReader.cpp" 500)
p9c3_require_max_lines("src/core/BpxExpression.cpp" 400)
p9c3_require_max_lines("src/core/StrictJson.cpp" 400)
p9c3_require_max_lines("src/core/Experiment.cpp" 450)
# Reviewed MC-1 exception: the 771-line runner is one transaction boundary.
# A further TU split would change same-TU non-IPO inlining and could perturb the
# frozen Release trace; revisit only with a stable runner transaction seam.
p9c3_require_max_lines("src/core/CyclerV2.cpp" 800)
p9c3_require_max_lines("src/core/detail/ParameterCurve.hpp" 200)
p9c3_require_max_lines("src/core/detail/BpxExpression.hpp" 200)
p9c3_require_max_lines("src/core/detail/StrictJson.hpp" 200)
p9c3_require_max_lines("src/core/detail/ExperimentSemantics.hpp" 200)

message(STATUS "9C-3 cold-file split structural gate passed")
