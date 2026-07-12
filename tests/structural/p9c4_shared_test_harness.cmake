# M0.8 / 9C-4: one narrow, caller-owned SPM test harness with explicit units.

if(NOT DEFINED SLIDE_SOURCE_DIR)
  message(FATAL_ERROR "SLIDE_SOURCE_DIR is required")
endif()

function(p9c4_load_compact relative_path output)
  set(path "${SLIDE_SOURCE_DIR}/${relative_path}")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "9C-4 source is missing: ${relative_path}")
  endif()
  file(READ "${path}" content)
  string(REGEX REPLACE "//[^\r\n]*" "" content "${content}")
  string(REGEX REPLACE "/\\*([^*]|\\*+[^*/])*\\*+/" "" content "${content}")
  string(REGEX REPLACE "[ \t\r\n]" "" content "${content}")
  set(${output} "${content}" PARENT_SCOPE)
endfunction()

function(p9c4_require_tokens label content_variable)
  foreach(token IN ITEMS ${ARGN})
    string(FIND "${${content_variable}}" "${token}" position)
    if(position EQUAL -1)
      message(FATAL_ERROR
        "9C-4 ${label}: required discriminator is absent: ${token}")
    endif()
  endforeach()
endfunction()

function(p9c4_forbid_tokens label content_variable)
  foreach(token IN ITEMS ${ARGN})
    string(FIND "${${content_variable}}" "${token}" position)
    if(NOT position EQUAL -1)
      message(FATAL_ERROR
        "9C-4 ${label}: forbidden ownership/default token remains: ${token}")
    endif()
  endforeach()
endfunction()

function(p9c4_require_token_count label content_variable token expected)
  string(LENGTH "${${content_variable}}" original_length)
  string(REPLACE "${token}" "" without_token "${${content_variable}}")
  string(LENGTH "${without_token}" stripped_length)
  string(LENGTH "${token}" token_length)
  if(token_length EQUAL 0)
    message(FATAL_ERROR "9C-4 internal error: empty count token for ${label}")
  endif()
  math(EXPR count
    "(${original_length} - ${stripped_length}) / ${token_length}")
  if(NOT count EQUAL expected)
    message(FATAL_ERROR
      "9C-4 ${label}: expected ${expected} occurrences of ${token}, found ${count}")
  endif()
endfunction()

function(p9c4_require_ordered_tokens label content_variable)
  set(search_offset 0)
  foreach(token IN LISTS ARGN)
    string(SUBSTRING "${${content_variable}}" ${search_offset} -1 remaining)
    string(FIND "${remaining}" "${token}" relative_position)
    if(relative_position EQUAL -1)
      message(FATAL_ERROR
        "9C-4 ${label}: ordered discriminator is absent after offset ${search_offset}: ${token}")
    endif()
    string(LENGTH "${token}" token_length)
    math(EXPR search_offset
      "${search_offset} + ${relative_position} + ${token_length}")
  endforeach()
endfunction()

function(p9c4_require_max_lines relative_path maximum)
  set(path "${SLIDE_SOURCE_DIR}/${relative_path}")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "9C-4 source is missing: ${relative_path}")
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
      "9C-4 ${relative_path}: ${line_count} physical lines exceeds ceiling ${maximum}")
  endif()
endfunction()

p9c4_load_compact("tests/support/CoreSpmTestHarness.hpp" harness)
p9c4_load_compact("tests/unit/core_CoreSpmTestHarness_test.cpp" self_test)
p9c4_load_compact("tests/unit/CMakeLists.txt" unit_cmake)
string(REPLACE "static_cast<" "" harness_without_casts "${harness}")

# The build seam returns a move-only batch to its caller. Observation and trace
# seams accept all remaining storage from the caller, with units in the types.
p9c4_require_tokens("public seams" harness
  "#include\"../../src/core/Numeric.hpp\""
  "structCurrentA{explicitCurrentA(std::span<constcore::real_t>lane_values)noexcept"
  "structCurrentDensityApm2{explicitCurrentDensityApm2(std::span<constcore::real_t>lane_values)noexcept"
  "[[nodiscard]]inlinecore::SpmBatchrequireSpmBatch(constcore::SpmFactoryInput&input,constcore::SpmModelOptions&options,intn_lanes)"
  "core::buildSpmBatch(input,options,n_lanes,batch)==Status::Success"
  "returnbatch"
  "inlinevoidrequireTerminalVoltage(core::SpmBatch&batch,CurrentDensityApm2current_density,constcore::real_t&time_s,std::span<core::real_t>terminal_voltage_V)"
  "inlinevoidrequireTerminalVoltage(core::SpmBatch&batch,CurrentAcurrent,constcore::real_t&time_s,std::span<core::real_t>current_density_scratch_Apm2,std::span<core::real_t>terminal_voltage_V)"
  "inlinevoidrequireConstantCurrentTrace(core::SpmBatch&batch,core::ExponentialModal&stepper,CurrentDensityApm2current_density,std::span<constcore::real_t>sample_time_s,std::span<core::real_t>sample_major_terminal_voltage_V)"
  "[[nodiscard]]inlineboolcomputeVoltageError(std::span<constcore::real_t>actual_voltage_V,std::span<constcore::real_t>expected_voltage_V,VoltageError&output)noexcept")
p9c4_require_token_count(
  "single successful build assertion" harness
  "core::buildSpmBatch(input,options,n_lanes,batch)==Status::Success" 1)
p9c4_require_token_count(
  "single build seam" harness "requireSpmBatch(" 1)
p9c4_require_token_count(
  "single trace seam" harness "requireConstantCurrentTrace(" 1)
p9c4_require_token_count(
  "single metric seam" harness "computeVoltageError(" 1)
p9c4_require_token_count(
  "only returned batch is built internally" harness "core::SpmBatchbatch" 1)

# Amperes are admitted for every lane before any quotient is evaluated. The
# quotient expression occurs once, is bit-checked, then is published once.
p9c4_require_ordered_tokens("current conversion" harness
  "inlinevoidrequireTerminalVoltage(core::SpmBatch&batch,CurrentAcurrent,"
  "constcore::real_tarea=batch.electrode_area()"
  "REQUIRE(core::is_finite(area))"
  "REQUIRE(area>0.0)"
  "for(constcore::real_t&current_A:current.lane_values){"
  "REQUIRE(core::is_finite(current_A))"
  "REQUIRE(detail::canDivideByPositiveNormal(current_A,area))"
  "for(std::size_tlane=0;lane<lane_count;++lane){"
  "constcore::real_tdensity_Apm2=current.lane_values[lane]/area"
  "REQUIRE(core::is_finite(density_Apm2))"
  "current_density_scratch_Apm2[lane]=density_Apm2"
  "requireTerminalVoltage(batch,CurrentDensityApm2{current_density_scratch_Apm2},time_s,terminal_voltage_V)")
p9c4_require_tokens("division overflow preflight" harness
  "inlineboolcanDivideByPositiveNormal("
  "denominator_bits&exponent_mask"
  "returnnumerator_magnitude<std::numeric_limits<core::real_t>::max()*positive_denominator")
p9c4_require_token_count(
  "one evaluated current quotient" harness
  "current.lane_values[lane]/area" 1)
p9c4_require_token_count(
  "one quotient publication" harness
  "current_density_scratch_Apm2[lane]=density_Apm2" 1)
p9c4_require_token_count(
  "opaque-reference density validation" harness
  "for(constcore::real_t&density_Apm2:current_density.lane_values)" 2)
p9c4_require_token_count(
  "opaque-reference voltage validation" harness
  "for(constcore::real_t&voltage_V:terminal_voltage_V)" 1)
p9c4_forbid_tokens("no evaluated division in preflight" harness
  "numerator/positive_denominator"
  "canDivideByPositiveFinite")

# Every writable span is disjoint from spans which remain live. Pointer range
# comparison uses the standard-library total order, including unrelated arrays.
p9c4_require_tokens("portable overlap contract" harness
  "inlineboolspansOverlap(std::span<constcore::real_t>left,std::span<constcore::real_t>right)noexcept"
  "constautoless=std::less<constcore::real_t*>{}"
  "returnless(left.data(),right_end)&&less(right.data(),left_end)"
  "detail::spansOverlap(current_density.lane_values,terminal_voltage_V)"
  "detail::spansOverlap(current.lane_values,current_density_scratch_Apm2)"
  "detail::spansOverlap(current.lane_values,terminal_voltage_V)"
  "detail::spansOverlap(current_density_scratch_Apm2,terminal_voltage_V)"
  "detail::spansOverlap(current_density.lane_values,sample_major_terminal_voltage_V)"
  "detail::spansOverlap(sample_time_s,sample_major_terminal_voltage_V)")
p9c4_require_token_count(
  "all mutable-span overlap guards" harness "detail::spansOverlap(" 6)

# Full-grid admission precedes t0 observation; the stepper is configured only
# after that observation and receives the left endpoint plus adjacent delta.
p9c4_require_ordered_tokens("trace framing" harness
  "inlinevoidrequireConstantCurrentTrace("
  "for(std::size_tsample=1;sample<sample_time_s.size();++sample)REQUIRE((core::is_finite(sample_time_s[sample])&&sample_time_s[sample]>sample_time_s[sample-1]))"
  "requireTerminalVoltage(batch,current_density,sample_time_s.front(),sample_major_terminal_voltage_V.first(lane_count))"
  "stepper.configure(batch)==Status::Success"
  "for(std::size_tsample=1;sample<sample_time_s.size();++sample){"
  "constcore::real_tstart_time_s=sample_time_s[sample-1]"
  "constcore::real_tdt_s=sample_time_s[sample]-start_time_s"
  "stepper.step(batch,current_density.lane_values,start_time_s,dt_s)==Status::Success"
  "constautovoltage_V=stepper.terminalVoltage()"
  "REQUIRE(voltage_V.size()==lane_count)"
  "std::ranges::copy(voltage_V,sample_major_terminal_voltage_V.subspan(sample*lane_count,lane_count).begin())"
  "std::ranges::all_of(sample_major_terminal_voltage_V,[](constcore::real_t&voltage_V){returncore::is_finite(voltage_V);})")
p9c4_require_tokens("observation time forwarding" harness
  "constcore::StepCtxobservation{.time=time_s,.dt=0.0,.i_app=current_density.lane_values,}")
p9c4_forbid_tokens("no hidden trace origin" harness
  "sample_time_s.front()==0.0")
p9c4_require_token_count(
  "one grid interval admission" harness
  "sample_time_s[sample]>sample_time_s[sample-1]" 1)
p9c4_require_token_count(
  "one accepted step per interval" harness
  "stepper.step(batch,current_density.lane_values,start_time_s,dt_s)==Status::Success" 1)

# Metric arithmetic mirrors P7 exactly: sample zero initializes both reductions,
# then samples 1..N-1 accumulate in order. Each risky operation is preflighted.
p9c4_require_ordered_tokens("literal metric" harness
  "detail::canEvaluateAbsoluteDifference(actual_voltage_V[0],expected_voltage_V[0])"
  "maximum_absolute_V=std::abs(actual_voltage_V[0]-expected_voltage_V[0])"
  "core::is_finite(maximum_absolute_V)&&detail::canSquare(maximum_absolute_V)"
  "square_sum_V2=maximum_absolute_V*maximum_absolute_V"
  "core::is_finite(square_sum_V2)"
  "for(std::size_tsample=1;sample<actual_voltage_V.size();++sample){"
  "constcore::real_t&actual_V=actual_voltage_V[sample]"
  "constcore::real_t&expected_V=expected_voltage_V[sample]"
  "detail::canEvaluateAbsoluteDifference(actual_V,expected_V)"
  "constcore::real_terror_V=std::abs(actual_V-expected_V)"
  "core::is_finite(error_V)"
  "maximum_absolute_V=std::max(maximum_absolute_V,error_V)"
  "detail::canSquare(error_V)&&detail::canAddSquare(square_sum_V2,error_V)"
  "square_sum_V2+=error_V*error_V"
  "core::is_finite(square_sum_V2)"
  "std::sqrt(square_sum_V2/static_cast<core::real_t>(actual_voltage_V.size()))"
  "core::is_finite(rms_V)"
  "output={maximum_absolute_V,rms_V}")
p9c4_require_tokens("metric overflow preflights" harness
  "inlineboolcanEvaluateAbsoluteDifference("
  "(left_bits^right_bits)&sign_mask"
  "returnlarger<std::numeric_limits<core::real_t>::max()-smaller"
  "inlineboolcanSquare("
  "returnvalue<std::numeric_limits<core::real_t>::max()/value"
  "inlineboolcanAddSquare("
  "std::numeric_limits<core::real_t>::max()-sum"
  "if(value<=1.0)returnroom>value"
  "returnvalue<room/value")
p9c4_require_token_count(
  "metric output is atomically published" harness
  "output={maximum_absolute_V,rms_V}" 1)
p9c4_forbid_tokens("no sample-zero fold" harness
  "for(std::size_tsample=0")
p9c4_forbid_tokens("transparent metric" harness
  "longdouble"
  "scale_V"
  "scaled_square_sum"
  "std::hypot"
  "Catch::Approx"
  "Approx("
  "tolerance"
  "epsilon"
  "margin")

# The harness must not own fixtures, references, storage, steppers, policy, or
# hidden process state. These tokens are negative mutation discriminators.
p9c4_forbid_tokens("narrow ownership and policy" harness
  "Kokam"
  "chen2020"
  "ParameterSet"
  "capacity_Ah"
  "c_rate"
  "C_rate"
  "SpmModelOptions{"
  "=default"
  ".nch="
  ".thermal="
  "std::vector"
  "std::array"
  "std::unique_ptr"
  "std::shared_ptr"
  "make_unique"
  "make_shared"
  "filesystem"
  "fstream"
  "ifstream"
  "ofstream"
  "SLIDE_SOURCE_DIR"
  "thread_local"
  "mutable"
  "core::ExponentialModalstepper")
p9c4_forbid_tokens("no static process state" harness_without_casts "static")

# Independent coverage pins configuration, both current units, caller storage,
# nonuniform/partial framing, literal metric order, final max, and IEEE rejection.
p9c4_require_tokens("independent self-test" self_test
  "constexprcore::SpmModelOptionsharness_options{.nch=8,.thermal=true}"
  "constexprintharness_lanes=3"
  "constexprcore::SpmModelOptionsalternate_harness_options{.nch=12,.thermal=false}"
  "constexprintalternate_harness_lanes=2"
  "input.design.electrode_area=0.25"
  "test_support::requireSpmBatch(input,harness_options,harness_lanes)"
  "input.design.electrode_area=0.4"
  "test_support::requireSpmBatch(input,alternate_harness_options,alternate_harness_lanes)"
  "test_support::CurrentA{std::span<constcore::real_t>{current_A}}"
  "test_support::CurrentDensityApm2{std::span<constcore::real_t>{expected_density_Apm2}}"
  "expected_density_Apm2{1.0,-2.0,0.5}"
  "expected_density_Apm2{1.0,-0.5}"
  "density_scratch_Apm2==expected_density_Apm2"
  "constexprstd::array<core::real_t,5>sample_time_s{2.0,2.125,2.5,2.75,2.875}"
  ".time=sample_time_s.front()"
  ".dt=0.0"
  "manual_stepper.step(manual_batch,current_density_Apm2,start_time_s,dt_s)"
  "batch.state().raw().data()==state_storage"
  "actual_voltage_V.data()==output_storage"
  "actual_voltage_V==expected_voltage_V"
  "literal_error_V{1.0,2.0,-1.0}"
  "final_maximum_error_V{1.0,-1.0,2.0}"
  "error.maximum_absolute_V==2.0"
  "error.rms_V==std::sqrt(2.0)"
  "UINT64_C(0x7ff8000000000000)"
  "UINT64_C(0x7ff0000000000000)"
  "std::span<constcore::real_t>{literal_error_V}.first(2)"
  "std::span<constcore::real_t>{},std::span<constcore::real_t>{},error"
  "volatilestd::uint64_truntime_bits=bits"
  "test_support::detail::spansOverlap(first,overlapping)"
  "test_support::detail::spansOverlap(first,separate_storage)"
  "canDivideByPositiveNormal(maximum,0.5)"
  "UINT64_C(0x0000000000000001)"
  "canAddSquare(maximum,minimum_subnormal)"
  "UINT64_C(0x7fefffffffffffff)"
  "UINT64_C(0xffefffffffffffff)"
  "UINT64_C(0x5ff0000000000000)"
  "UINT64_C(0x5fe0000000000000)"
  "core::try_multiply_nonnegative(two_to_511,two_to_511,witness_square)")
p9c4_require_token_count(
  "three distinct fast-math overflow tests" self_test "[fast-math]" 3)
p9c4_require_token_count(
  "both explicit build tuples are exercised" self_test
  "test_support::requireSpmBatch(" 3)
p9c4_forbid_tokens("self-test independence" self_test
  "KokamSpmFixture"
  "make_legacy_kokam"
  "Catch::Approx"
  "Approx(")

p9c4_require_token_count(
  "self-test target registration" unit_cmake
  "add_executable_with_coverage_and_test(unit_test_core_CoreSpmTestHarnesscore_CoreSpmTestHarness_test.cpp)" 1)
p9c4_require_token_count(
  "self-test core-only link" unit_cmake
  "target_link_libraries(unit_test_core_CoreSpmTestHarnessPRIVATEslide_core)" 1)
p9c4_require_token_count(
  "self-test exercises fast-math validation" unit_cmake
  "slide_use_strict_fp(core_CoreSpmTestHarness_test.cpp)" 0)

p9c4_require_max_lines("tests/support/CoreSpmTestHarness.hpp" 360)
p9c4_require_max_lines("tests/unit/core_CoreSpmTestHarness_test.cpp" 480)
p9c4_require_max_lines("tests/structural/p9c4_shared_test_harness.cmake" 380)

message(STATUS "9C-4 shared test harness structural gate passed")
