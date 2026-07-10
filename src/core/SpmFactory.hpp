/**
 * @file SpmFactory.hpp
 * @brief Cold SPM parameter compiler and D-02 batch-composition registry.
 */

#pragma once

#include "SpmPipeline.hpp"

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace slide::core {

inline constexpr std::array registered_spm_nch{ 5, 8, 12 };

enum class SpmComposition : unsigned char {
  isothermal,
  thermal,
  isothermal_ageing,
  thermal_ageing
};

struct SpmModelOptions
{
  int nch{ 5 };
  bool thermal{};
  std::uint8_t sei_model_mask{};
  bool sei_porosity{};
  std::uint8_t surface_crack_model_mask{};
  bool surface_crack_diffusivity{};
  std::uint8_t lam_model_mask{};
  bool lithium_plating{};

  bool has_ageing() const
  {
    return sei_model_mask != 0 || surface_crack_model_mask != 0
           || lam_model_mask != 0 || lithium_plating;
  }
};

/** Cold, value-semantic inputs not already represented by CellDesign. */
struct SpmFactoryInput
{
  CellDesign design{};
  OCVCurve total_entropic_coefficient{};
  OCVCurve negative_entropic_coefficient{};
  OCVCurve negative_laresgoiti_stress{};

  real_t initial_soc{ 0.5 };
  real_t initial_temperature{}; //!< zero selects design.thermal.reference_temperature
  real_t initial_sei_thickness{ 1e-9 };
  real_t initial_lost_lithium{};
  real_t initial_crack_surface_fraction{ 0.01 };
  real_t initial_plated_lithium_thickness{};
  PerDomain<real_t> initial_specific_resistance{};
  real_t initial_current_collector_resistance{};
  real_t initial_stress_interval{ 1.0 };
  real_t sei_resistivity_area{};

  SeiParams sei{};
  SurfaceCrackParams surface_crack{};
  LamParams lam{};
  LithiumPlatingParams lithium_plating{};
};

/**
 * Move-only runtime batch. `rhs()` performs exactly one indirect call, after which the concrete
 * compile-time pipeline executes with no virtual/function-pointer work inside its lane sweeps.
 */
class SpmBatch
{
public:
  SpmBatch() = default;
  ~SpmBatch();
  SpmBatch(const SpmBatch &) = delete;
  SpmBatch &operator=(const SpmBatch &) = delete;
  SpmBatch(SpmBatch &&other) noexcept;
  SpmBatch &operator=(SpmBatch &&other) noexcept;

  bool valid() const { return implementation_ != nullptr; }
  int nch() const { return nch_; }
  int n_lanes() const { return state_.n_lanes(); }
  real_t capacity_Ah() const { return capacity_Ah_; }
  real_t electrode_area() const { return electrode_area_; }
  SpmComposition composition() const { return composition_; }
  const SpmPipelineLayout &layout() const { return layout_; }
  std::span<const StateRole> roles() const { return roles_; }

  StateArena &state() { return state_; }
  const StateArena &state() const { return state_; }
  StateArena &derivative() { return derivative_; }
  const StateArena &derivative() const { return derivative_; }

  [[nodiscard]] slide::Status rhs(std::span<const real_t> trial_state,
                                  std::span<real_t> trial_derivative,
                                  const StepCtx &ctx);
  [[nodiscard]] slide::Status evaluate(const StepCtx &ctx);
  bool hasFusedEuler() const { return fused_euler_ != nullptr; }
  [[nodiscard]] slide::Status fusedEuler(const StepCtx &ctx, real_t dt,
                                         std::span<real_t> terminal_voltage);
  [[nodiscard]] slide::Status exponentialStep(const StepCtx &ctx, real_t dt,
                                              std::span<real_t> terminal_voltage);
  [[nodiscard]] slide::Status terminalVoltage(const StepCtx &ctx,
                                              std::span<real_t> output);
  /** Recompute voltage lazily from a full recorded arena snapshot. */
  [[nodiscard]] slide::Status terminalVoltageAt(
    std::span<const real_t> snapshot,
    std::span<const real_t> current_density,
    std::span<real_t> output);
  [[nodiscard]] slide::Status linearizeThevenin(
    std::span<const real_t> current,
    std::span<real_t> intercept_ocv,
    std::span<real_t> resistance);
  /** Cold pack compiler hook: enable an already-validated repeated-lane brick. */
  [[nodiscard]] slide::Status setTrustedLanePeriod(int maximum_period);
  int trustedLanePeriod() const;
  [[nodiscard]] slide::Status storeStressHistory(real_t interval);

private:
  using EvaluateFn = slide::Status (*)(void *, RhsViews &, const StepCtx &);
  using ObserveVoltageFn = slide::Status (*)(void *, const ConstBatchView &,
                                             const StepCtx &, std::span<real_t>);
  using LinearizeTheveninFn = slide::Status (*)(void *, const ConstBatchView &,
                                                std::span<const real_t>,
                                                std::span<real_t>,
                                                std::span<real_t>);
  using SetLanePeriodFn = slide::Status (*)(void *, const ConstBatchView &, int);
  using GetLanePeriodFn = int (*)(const void *);
  using StoreStressFn = void (*)(void *, BatchView, real_t);
  using FusedEulerFn = slide::Status (*)(void *, BatchView, const StepCtx &, real_t,
                                         std::span<real_t>);
  using ExponentialFn = slide::Status (*)(void *, BatchView, BatchView,
                                          const StepCtx &, real_t,
                                          std::span<real_t>);
  using DestroyFn = void (*)(void *);

  SpmBatch(void *implementation,
           EvaluateFn evaluate,
           ObserveVoltageFn observe_voltage,
           LinearizeTheveninFn linearize_thevenin,
           SetLanePeriodFn set_lane_period,
           GetLanePeriodFn get_lane_period,
           StoreStressFn store_stress,
           FusedEulerFn fused_euler,
           ExponentialFn exponential,
           DestroyFn destroy,
           int nch,
           real_t capacity_Ah,
           real_t electrode_area,
           SpmComposition composition,
           SpmPipelineLayout layout,
           StateArena state,
           StateArena derivative,
           std::vector<StateRole>
             roles);
  void reset() noexcept;

  void *implementation_{};
  EvaluateFn evaluate_{};
  ObserveVoltageFn observe_voltage_{};
  LinearizeTheveninFn linearize_thevenin_{};
  SetLanePeriodFn set_lane_period_{};
  GetLanePeriodFn get_lane_period_{};
  StoreStressFn store_stress_{};
  FusedEulerFn fused_euler_{};
  ExponentialFn exponential_{};
  DestroyFn destroy_{};
  int nch_{};
  real_t capacity_Ah_{};
  real_t electrode_area_{};
  SpmComposition composition_{ SpmComposition::isothermal };
  SpmPipelineLayout layout_{};
  StateArena state_{};
  StateArena derivative_{};
  std::vector<StateRole> roles_{};

  friend struct SpmBatchFactoryAccess;
};

/** Build a homogeneous batch; `output` is unchanged if any cold validation gate fails. */
[[nodiscard]] slide::Status buildSpmBatch(const SpmFactoryInput &input,
                                          const SpmModelOptions &options,
                                          int n_lanes,
                                          SpmBatch &output);

} // namespace slide::core
