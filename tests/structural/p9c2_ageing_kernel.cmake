# M0.6 / 9C-2: four named ageing equation bodies, one production scaffold.

if(NOT DEFINED SLIDE_SOURCE_DIR)
  message(FATAL_ERROR "SLIDE_SOURCE_DIR is required")
endif()

function(load_compact relative_path output)
  set(path "${SLIDE_SOURCE_DIR}/${relative_path}")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "9C-2 source is missing: ${relative_path}")
  endif()
  file(READ "${path}" content)
  string(REGEX REPLACE "//[^\r\n]*" "" content "${content}")
  string(REGEX REPLACE "/\\*([^*]|\\*+[^*/])*\\*+/" "" content "${content}")
  string(REGEX REPLACE "[ \t\r\n]" "" content "${content}")
  set(${output} "${content}" PARENT_SCOPE)
endfunction()

function(require_tokens label content_variable)
  foreach(token IN ITEMS ${ARGN})
    string(FIND "${${content_variable}}" "${token}" position)
    if(position EQUAL -1)
      message(FATAL_ERROR "9C-2 ${label}: required token is absent: ${token}")
    endif()
  endforeach()
endfunction()

function(forbid_tokens label content_variable)
  foreach(token IN ITEMS ${ARGN})
    string(FIND "${${content_variable}}" "${token}" position)
    if(NOT position EQUAL -1)
      message(FATAL_ERROR "9C-2 ${label}: private scaffold remains: ${token}")
    endif()
  endforeach()
endfunction()

function(require_token_count label content_variable token expected)
  string(LENGTH "${${content_variable}}" original_length)
  string(REPLACE "${token}" "" without_token "${${content_variable}}")
  string(LENGTH "${without_token}" stripped_length)
  string(LENGTH "${token}" token_length)
  math(EXPR count "(${original_length} - ${stripped_length}) / ${token_length}")
  if(NOT count EQUAL expected)
    message(
      FATAL_ERROR
      "9C-2 ${label}: expected ${expected} occurrences of ${token}, found ${count}")
  endif()
endfunction()

load_compact("src/core/Sei.hpp" sei)
load_compact("src/core/AgeingKernel.hpp" scaffold)
load_compact("src/core/SurfaceCrack.hpp" crack)
load_compact("src/core/Lam.hpp" lam)
load_compact("src/core/LithiumPlating.hpp" plating)
load_compact("src/core/SpmStress.hpp" stress)
load_compact("src/core/SpmPipeline.hpp" pipeline)
load_compact("src/core/SpmFactory.cpp" factory)
load_compact("src/cells/Cell_SPM/Cell_SPM_degradation.cpp" legacy)

foreach(consumer IN ITEMS sei crack lam plating stress pipeline)
  require_tokens("${consumer}" ${consumer} "#include\"AgeingKernel.hpp\"")
endforeach()

require_tokens("shared scaffold" scaffold
  "for(unsignedmodel=1;model<=ModelCount;++model)"
  "SLIDE_AGEING_FORCE_INLINEslide::Statusfor_each_enabled_ageing_model_lane("
  "SLIDE_AGEING_FORCE_INLINEslide::Statusevaluate_ageing_stage(")

require_tokens("SEI" sei
  "ageing_model_bit<4>(model)"
  "valid_ageing_model_mask<4>(p.model_mask)"
  "detail::AgeingScratchStorage<Real,2>"
  "detail::clear_ageing_fields<Real,2>("
  "detail::for_each_enabled_ageing_model_lane<4>("
  "detail::for_each_ageing_lane_while_success("
  "detail::for_each_ageing_lane("
  "return{storage_.field(0),storage_.field(1)}")

require_tokens("surface crack" crack
  "ageing_model_bit<5>(model)"
  "valid_ageing_model_mask<5>(p.model_mask)"
  "detail::AgeingScratchStorage<Real,3>"
  "detail::clear_ageing_fields<Real,3>("
  "detail::for_each_enabled_ageing_model_lane<5>("
  "detail::for_each_ageing_lane_while_success("
  "detail::for_each_ageing_lane("
  "return{storage_.field(0),storage_.field(1),storage_.field(2)}")

require_tokens("LAM" lam
  "ageing_model_bit<4>(model)"
  "valid_ageing_model_mask<4>(p.model_mask)"
  "detail::AgeingScratchStorage<Real,6>"
  "detail::clear_ageing_fields<Real,6>("
  "detail::for_each_enabled_ageing_model_lane<4>("
  "detail::for_each_ageing_lane_while_success("
  "detail::for_each_ageing_lane("
  "(*field)[domain_index(domain)]=storage_.field(cursor)"
  "++cursor")

require_tokens("lithium plating" plating
  "BasicLithiumPlatingOutput"
  "detail::AgeingScratchStorage<Real,1>"
  "return{storage_.field(0)}"
  "detail::for_each_ageing_lane_while_success("
  "detail::for_each_ageing_lane(")

require_tokens("stress scratch" stress "detail::AgeingScratchStorage<Real,3>")
require_token_count("pipeline stages" pipeline "detail::evaluate_ageing_stage(" 4)
require_tokens("pipeline disabled-SEI dependency" pipeline
  "detail::clear_ageing_fields<real_t,2>(")
require_tokens("factory masks" factory
  "valid_optional_ageing_model_mask<4>(options.sei_model_mask)"
  "valid_optional_ageing_model_mask<5>(options.surface_crack_model_mask)"
  "valid_optional_ageing_model_mask<4>(options.lam_model_mask)")

foreach(mechanism IN ITEMS sei crack lam plating)
  forbid_tokens("${mechanism}" ${mechanism}
    "std::vector<Real>storage_"
    "std::fill("
    "for(unsignedmodel=1;model<="
    "for(intlane=0;lane<lanes;++lane)")
endforeach()

# The legacy Cell_SPM equations are the independent parity oracle. Sharing the
# production helper here would turn mechanism parity into a common-mode test.
forbid_tokens("legacy oracle" legacy
  "AgeingKernel.hpp"
  "AgeingScratchStorage"
  "for_each_enabled_ageing_model_lane"
  "evaluate_ageing_stage")

message(STATUS "9C-2 shared ageing-kernel structural gate passed")
