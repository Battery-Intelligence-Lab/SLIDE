# PC-10 is an architectural invariant: production adapters may own indexing,
# validation, and status publication, but not private copies of SPM algebra.

if(NOT DEFINED SLIDE_SOURCE_DIR)
  message(FATAL_ERROR "SLIDE_SOURCE_DIR is required")
endif()

function(load_compact relative_path output)
  set(path "${SLIDE_SOURCE_DIR}/${relative_path}")
  if(NOT EXISTS "${path}")
    message(FATAL_ERROR "PC-10 source is missing: ${relative_path}")
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
      message(FATAL_ERROR "PC-10 ${label}: required token is absent: ${token}")
    endif()
  endforeach()
endfunction()

function(forbid_tokens label content_variable)
  foreach(token IN ITEMS ${ARGN})
    string(FIND "${${content_variable}}" "${token}" position)
    if(NOT position EQUAL -1)
      message(FATAL_ERROR "PC-10 ${label}: private physics copy remains: ${token}")
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
      "PC-10 ${label}: expected ${expected} occurrences of ${token}, found ${count}")
  endif()
endfunction()

function(slice_between label content_variable begin_marker end_marker output)
  string(FIND "${${content_variable}}" "${begin_marker}" begin)
  if(begin EQUAL -1)
    message(FATAL_ERROR "PC-10 ${label}: begin marker is absent: ${begin_marker}")
  endif()
  string(SUBSTRING "${${content_variable}}" ${begin} -1 tail)
  string(FIND "${tail}" "${end_marker}" length)
  if(length EQUAL -1)
    message(FATAL_ERROR "PC-10 ${label}: end marker is absent: ${end_marker}")
  endif()
  string(SUBSTRING "${tail}" 0 ${length} result)
  set(${output} "${result}" PARENT_SCOPE)
endfunction()

# M0.10 / 9C-6 moved the diffusion RHS mapping into its own header, beside every other
# mechanism's. The PC-10 invariant is untouched -- every CPU diffusion leaf still calls the one
# shared scalar kernel -- so the count is taken over the pipeline and the mapping together, and
# the expected total (3) is unchanged.
load_compact("src/core/SpmPipeline.hpp" pipeline)
load_compact("src/core/SpmDiffusionRhs.hpp" diffusion_rhs)
string(APPEND pipeline "${diffusion_rhs}")
load_compact("src/core/SpmObservables.hpp" observables)
load_compact("src/core/ForwardSensitivity.cpp" dual)
load_compact("src/core/CudaSpmRuntime.cu" cuda)
load_compact("src/core/CompiledCurve.hpp" curves)
load_compact("src/core/SpectralDiffusion.hpp" spectral)
load_compact("src/core/SpectralDiffusionLegacy.hpp" spectral_legacy)
load_compact("src/core/SpmScalarKernels.hpp" scalar_kernels)
load_compact("src/core/Sei.hpp" sei)

foreach(consumer IN ITEMS pipeline observables dual cuda curves spectral)
  require_tokens("${consumer}" ${consumer} "#include\"SpmScalarKernels.hpp\"")
endforeach()

slice_between(
  "CPU modal slice"
  pipeline
  "[[nodiscard]]slide::StatusadvanceExponential("
  "[[nodiscard]]slide::StatusobserveTerminalVoltage("
  pipeline_modal)
require_tokens("CPU modal" pipeline_modal "SLIDE_SPM_ADVANCE_MODAL_STD(")
forbid_tokens("CPU modal" pipeline_modal
  "std::expm1(" "std::exp(x)" "1.0/6.0" "/24.0")

require_token_count(
  "CPU diffusion leaves" pipeline "spm_scalar::diffusionRate(" 3)

require_tokens("CPU observables" observables
  "spm_scalar::arrheniusFactor("
  "spm_scalar::activatedValue("
  "spm_scalar::fluxDenominator("
  "spm_scalar::molarFlux("
  "spm_scalar::concentrationOutput("
  "spm_scalar::surfaceStoichiometry("
  "spm_scalar::exchangeCurrent("
  "spm_scalar::activationArgument("
  "spm_scalar::activationOverpotential("
  "spm_scalar::activeArea("
  "spm_scalar::seriesResistance("
  "spm_scalar::cellOpenCircuitVoltage("
  "spm_scalar::terminalVoltage(")
forbid_tokens("CPU observables" observables
  "std::exp(" "std::sqrt(" "std::asinh(")

require_tokens("Dual" dual
  "SLIDE_SPM_ADVANCE_MODAL_ADL("
  "spm_scalar::arrheniusFactor("
  "spm_scalar::activatedValue("
  "spm_scalar::fluxDenominator("
  "spm_scalar::molarFlux("
  "spm_scalar::concentrationOutput("
  "spm_scalar::surfaceStoichiometry("
  "spm_scalar::exchangeCurrent("
  "spm_scalar::activationArgument("
  "spm_scalar::activationOverpotential("
  "spm_scalar::activeArea("
  "spm_scalar::seriesResistance("
  "spm_scalar::cellOpenCircuitVoltage("
  "spm_scalar::terminalVoltage(")
forbid_tokens("Dual" dual
  "std::exp(" "std::expm1(" "std::sqrt(" "std::asinh("
  "1.0/6.0" "/24.0")

require_tokens("CUDA" cuda
  "SLIDE_SPM_ADVANCE_MODAL_CUDA("
  "spm_scalar::linearInterpolate("
  "spm_scalar::arrheniusFactor("
  "spm_scalar::activatedValue("
  "spm_scalar::fluxDenominator("
  "spm_scalar::molarFlux("
  "spm_scalar::concentrationOutput("
  "spm_scalar::surfaceStoichiometry("
  "spm_scalar::exchangeCurrent("
  "spm_scalar::activationArgument("
  "spm_scalar::activationOverpotential("
  "spm_scalar::activeArea("
  "spm_scalar::seriesResistance("
  "spm_scalar::cellOpenCircuitVoltage("
  "spm_scalar::terminalVoltage(")
forbid_tokens("CUDA" cuda
  "expm1(" "sqrt(" "asinh(" "1.0/6.0" "/24.0")

slice_between(
  "indexed curve slice"
  curves
  "template<classReal>Realeval(constReal&x)const"
  "real_tderivative(constreal_t&query)const"
  indexed_curve)
require_tokens("indexed curve" indexed_curve "spm_scalar::linearInterpolate(")
forbid_tokens("indexed curve" indexed_curve
  "y_[i+1]-y_[i]" "x_[i+1]-x_[i]")

require_tokens("spectral diffusion" spectral
  "spm_scalar::arrheniusFactor("
  "spm_scalar::activatedValue("
  "spm_scalar::fluxDenominator("
  "spm_scalar::molarFlux("
  "SLIDE_SPM_DIFFUSION_RATE(")
forbid_tokens("spectral diffusion" spectral
  "std::exp("
  "(1.0/p_.T_ref-1.0/T[c])/p_.Rg"
  "p_.a[d]*p_.n*p_.F*p_.thick[d]"
  "sgnd*i_app[c]/flux_den"
  "De[c]*Ak*z[c]+Bk*fl[c]")

forbid_tokens("spectral legacy oracle" spectral_legacy
  "SpmScalarKernels.hpp" "spm_scalar::")

slice_between(
  "scalar ageing owners"
  scalar_kernels
  "template<classScalar>SLIDE_SPM_HOST_DEVICEinlineScalararrheniusFactor("
  "template<classArea,classElectronCount,classFaraday,classThickness>"
  scalar_ageing)
slice_between(
  "SEI consumer"
  sei
  "template<classReal>[[nodiscard]]slide::StatuscomputeSei("
  "template<intNCH>structSeiRhsParams"
  sei_consumer)

require_token_count(
  "scalar ageing owner" scalar_ageing "SLIDE_SPM_SEI_KINETIC_CURRENT(" 1)
require_token_count(
  "SEI activation consumers" sei_consumer "spm_scalar::activatedValue(" 5)
require_token_count(
  "SEI kinetic consumers" sei_consumer "SLIDE_SPM_SEI_KINETIC_CURRENT(" 3)
set(sei_owner_and_consumer "${scalar_ageing}${sei_consumer}")
require_token_count(
  "SEI kinetic owner and consumers" sei_owner_and_consumer
  "SLIDE_SPM_SEI_KINETIC_CURRENT(" 4)
require_tokens("SEI direct scalar owner" sei
  "#include\"SpmScalarKernels.hpp\"")
forbid_tokens("SEI private activation copies" sei_consumer
  "p.model1_k*exp(p.model1_k_activation*arrhenius)"
  "p.model2_k*exp(p.model2_k_activation*arrhenius)"
  "p.model2_D*exp(p.model2_D_activation*arrhenius)"
  "k_ref*exp(k_activation*arrhenius)"
  "D_ref*exp(D_activation*arrhenius)")
forbid_tokens("SEI private kinetic copies" sei_consumer
  "p.n_sei*p.F*kt*exp(")

message(STATUS "PC-10 single-source structural gate passed")
