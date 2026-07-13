/**
 * @file SpmBatch.cpp
 * @brief The runtime batch object: move semantics and the one indirect call per operation.
 *
 * Owns: every `SpmBatch` member. Implements PLAN.md §3.2. Hot: `rhs`/`evaluate`/`exponentialStep`
 * dispatch once through a function pointer into the compiled pipeline. This file deliberately
 * does NOT include `SpmPipeline.hpp` -- the batch reaches its implementation only through the
 * type-erased pointers the factory installs, which is what makes the erasure real.
 */

#include "SpmFactory.hpp"

#include <utility>

namespace slide::core {

SpmBatch::SpmBatch(void *implementation,
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
                   SpmBatchLayout layout,
                   StateArena state,
                   StateArena derivative,
                   std::vector<StateRole>
                     roles)
  : implementation_{ implementation }, evaluate_{ evaluate },
    observe_voltage_{ observe_voltage },
    linearize_thevenin_{ linearize_thevenin },
    set_lane_period_{ set_lane_period },
    get_lane_period_{ get_lane_period },
    store_stress_{ store_stress }, fused_euler_{ fused_euler },
    exponential_{ exponential }, destroy_{ destroy }, nch_{ nch },
    capacity_Ah_{ capacity_Ah }, electrode_area_{ electrode_area },
    composition_{ composition }, layout_{ layout }, state_{ std::move(state) },
    derivative_{ std::move(derivative) }, roles_{ std::move(roles) }
{}

SpmBatch::~SpmBatch()
{
  reset();
}

SpmBatch::SpmBatch(SpmBatch &&other) noexcept
  : implementation_{ std::exchange(other.implementation_, nullptr) },
    evaluate_{ std::exchange(other.evaluate_, nullptr) },
    observe_voltage_{ std::exchange(other.observe_voltage_, nullptr) },
    linearize_thevenin_{ std::exchange(other.linearize_thevenin_, nullptr) },
    set_lane_period_{ std::exchange(other.set_lane_period_, nullptr) },
    get_lane_period_{ std::exchange(other.get_lane_period_, nullptr) },
    store_stress_{ std::exchange(other.store_stress_, nullptr) },
    fused_euler_{ std::exchange(other.fused_euler_, nullptr) },
    exponential_{ std::exchange(other.exponential_, nullptr) },
    destroy_{ std::exchange(other.destroy_, nullptr) },
    nch_{ std::exchange(other.nch_, 0) },
    capacity_Ah_{ std::exchange(other.capacity_Ah_, 0.0) },
    electrode_area_{ std::exchange(other.electrode_area_, 0.0) },
    composition_{ other.composition_ },
    layout_{ other.layout_ }, state_{ std::move(other.state_) },
    derivative_{ std::move(other.derivative_) }, roles_{ std::move(other.roles_) }
{}

SpmBatch &SpmBatch::operator=(SpmBatch &&other) noexcept
{
  if (this != &other) {
    reset();
    implementation_ = std::exchange(other.implementation_, nullptr);
    evaluate_ = std::exchange(other.evaluate_, nullptr);
    observe_voltage_ = std::exchange(other.observe_voltage_, nullptr);
    linearize_thevenin_ = std::exchange(other.linearize_thevenin_, nullptr);
    set_lane_period_ = std::exchange(other.set_lane_period_, nullptr);
    get_lane_period_ = std::exchange(other.get_lane_period_, nullptr);
    store_stress_ = std::exchange(other.store_stress_, nullptr);
    fused_euler_ = std::exchange(other.fused_euler_, nullptr);
    exponential_ = std::exchange(other.exponential_, nullptr);
    destroy_ = std::exchange(other.destroy_, nullptr);
    nch_ = std::exchange(other.nch_, 0);
    capacity_Ah_ = std::exchange(other.capacity_Ah_, 0.0);
    electrode_area_ = std::exchange(other.electrode_area_, 0.0);
    composition_ = other.composition_;
    layout_ = other.layout_;
    state_ = std::move(other.state_);
    derivative_ = std::move(other.derivative_);
    roles_ = std::move(other.roles_);
  }
  return *this;
}

void SpmBatch::reset() noexcept
{
  if (implementation_ != nullptr)
    destroy_(implementation_);
  implementation_ = nullptr;
  evaluate_ = nullptr;
  observe_voltage_ = nullptr;
  linearize_thevenin_ = nullptr;
  set_lane_period_ = nullptr;
  get_lane_period_ = nullptr;
  store_stress_ = nullptr;
  fused_euler_ = nullptr;
  exponential_ = nullptr;
  destroy_ = nullptr;
  nch_ = 0;
  capacity_Ah_ = 0.0;
  electrode_area_ = 0.0;
  state_ = {};
  derivative_ = {};
  roles_.clear();
}

slide::Status SpmBatch::rhs(std::span<const real_t> trial_state,
                            std::span<real_t>
                              trial_derivative,
                            const StepCtx &ctx)
{
  if (!valid() || trial_state.size() != state_.size()
      || trial_derivative.size() != derivative_.size()
      || static_cast<int>(ctx.i_app.size()) != n_lanes()
      || !is_finite(ctx.time) || !is_finite(ctx.dt))
    return slide::Status::Invalid_parameters;
  RhsViews views{ BatchShape::from(state_) };
  views.rebind(trial_state, trial_derivative);
  return evaluate_(implementation_, views, ctx);
}

slide::Status SpmBatch::evaluate(const StepCtx &ctx)
{
  return rhs(std::span<const real_t>{ state_.raw() }, derivative_.raw(), ctx);
}

slide::Status SpmBatch::fusedEuler(const StepCtx &ctx, real_t dt,
                                   std::span<real_t> terminal_voltage)
{
  if (!valid() || fused_euler_ == nullptr || static_cast<int>(ctx.i_app.size()) != n_lanes()
      || static_cast<int>(terminal_voltage.size()) != n_lanes()
      || !is_finite(ctx.time) || !is_finite(ctx.dt) || !is_finite(dt) || !(dt > 0.0))
    return slide::Status::Invalid_parameters;
  BatchView state_view{ BatchShape::from(state_), state_.raw() };
  return fused_euler_(implementation_, state_view, ctx, dt, terminal_voltage);
}

slide::Status SpmBatch::exponentialStep(const StepCtx &ctx, real_t dt,
                                        std::span<real_t> terminal_voltage)
{
  if (!valid() || exponential_ == nullptr
      || static_cast<int>(ctx.i_app.size()) != n_lanes()
      || static_cast<int>(terminal_voltage.size()) != n_lanes()
      || !is_finite(ctx.time) || !is_finite(ctx.dt)
      || !is_finite(dt) || !(dt > 0.0))
    return slide::Status::Invalid_parameters;
  BatchView state_view{ BatchShape::from(state_), state_.raw() };
  BatchView derivative_view{ BatchShape::from(derivative_), derivative_.raw() };
  return exponential_(implementation_, state_view, derivative_view, ctx, dt, terminal_voltage);
}

slide::Status SpmBatch::terminalVoltage(const StepCtx &ctx,
                                        std::span<real_t>
                                          output)
{
  if (!valid() || static_cast<int>(ctx.i_app.size()) != n_lanes()
      || static_cast<int>(output.size()) != n_lanes()
      || !is_finite(ctx.time) || !is_finite(ctx.dt))
    return slide::Status::Invalid_parameters;
  const ConstBatchView state_view{ BatchShape::from(state_),
                                   std::span<const real_t>{ state_.raw() } };
  return observe_voltage_(implementation_, state_view, ctx, output);
}

slide::Status SpmBatch::terminalVoltageAt(
  std::span<const real_t> snapshot,
  std::span<const real_t>
    current_density,
  std::span<real_t>
    output)
{
  if (!valid() || observe_voltage_ == nullptr || snapshot.size() != state_.size()
      || static_cast<int>(current_density.size()) != n_lanes()
      || output.size() != current_density.size())
    return slide::Status::Invalid_parameters;
  const ConstBatchView state_view{ BatchShape::from(state_), snapshot };
  const StepCtx ctx{ .time = 0.0, .dt = 0.0, .i_app = current_density };
  return observe_voltage_(implementation_, state_view, ctx, output);
}

slide::Status SpmBatch::linearizeThevenin(
  std::span<const real_t> current,
  std::span<real_t>
    intercept_ocv,
  std::span<real_t>
    resistance)
{
  if (!valid() || linearize_thevenin_ == nullptr
      || static_cast<int>(current.size()) != n_lanes()
      || intercept_ocv.size() != current.size()
      || resistance.size() != current.size())
    return slide::Status::Invalid_parameters;
  const ConstBatchView state_view{ BatchShape::from(state_),
                                   std::span<const real_t>{ state_.raw() } };
  return linearize_thevenin_(implementation_, state_view, current, intercept_ocv, resistance);
}

slide::Status SpmBatch::setTrustedLanePeriod(int maximum_period)
{
  if (!valid() || set_lane_period_ == nullptr || maximum_period <= 0
      || n_lanes() % maximum_period != 0)
    return slide::Status::Invalid_parameters;
  const ConstBatchView state_view{ BatchShape::from(state_),
                                   std::span<const real_t>{ state_.raw() } };
  return set_lane_period_(implementation_, state_view, maximum_period);
}

int SpmBatch::trustedLanePeriod() const
{
  return valid() && get_lane_period_ != nullptr
           ? get_lane_period_(implementation_)
           : n_lanes();
}

slide::Status SpmBatch::storeStressHistory(real_t interval)
{
  if (!valid() || !is_finite(interval) || !(interval > 0.0))
    return slide::Status::Invalid_parameters;
  BatchView view{ BatchShape::from(state_), state_.raw() };
  store_stress_(implementation_, view, interval);
  return slide::Status::Success;
}

} // namespace slide::core
