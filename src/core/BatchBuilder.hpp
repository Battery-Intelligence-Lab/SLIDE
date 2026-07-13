/**
 * @file BatchBuilder.hpp
 * @brief Build-time state layout for one batch archetype (PLAN.md §3.1).
 *
 * Components (electrode diffusion, thermal, each ageing model) declare StateSpecs at
 * build time; the builder performs layout and hands each component its StateSlice.
 * Adding a model appends rows. Layout is frozen by build(); never recomputed during
 * simulation. The builder itself is cold-path (strings, vectors allowed); nothing of
 * it survives into the hot loop except integer StateSlices.
 *
 * @date 2026-07-07
 * @surface support
 */

#pragma once

#include "StateArena.hpp"

#include <optional>
#include <vector>

namespace slide::core {

class BatchBuilder
{
public:
  //!< Append `spec.rows` rows; returns the component's integer handle.
  //!< Names must be unique within a batch (asserted) — recording looks state up by name.
  StateSlice declare(StateSpec spec)
  {
    assert(!frozen_ && "layout is frozen after build()");
    assert(spec.rows > 0);
    assert(!find(spec.name).has_value() && "duplicate state name in batch");

    const StateSlice s{ next_row_, spec.rows };
    specs_.push_back(spec);
    slices_.push_back(s);
    roles_.insert(roles_.end(), static_cast<std::size_t>(spec.rows), spec.role);
    next_row_ += spec.rows;
    return s;
  }

  /**
   * Q9 seam (PLAN.md §7, DECIDED 2026-07-07): cross-batch thermal coupling.
   *
   * Reserves ONE per-lane row "q_ext" [W] — external heat flow into each cell. The
   * in-batch thermal kernel (ThermalLumped) CONSUMES it as a source term; the pack-level
   * thermal solver (D-21, to be designed before Phase 2 opens) PRODUCES it between batch
   * steps. This is the entire declared interface — no coupling topology, no solver here.
   * Idempotent: repeated calls return the same slice.
   */
  StateSlice reserve_thermal_flux()
  {
    if (!q_ext_.has_value())
      q_ext_ = declare({ "q_ext", 1, Unit::W, StateRole::input });
    return *q_ext_;
  }

  //!< Cold-path lookup (recording, diagnostics, tests).
  std::optional<StateSlice> find(std::string_view name) const
  {
    for (std::size_t i{}; i < specs_.size(); i++)
      if (specs_[i].name == name)
        return slices_[i];
    return std::nullopt;
  }

  //!< Total rows declared so far.
  int rows() const { return next_row_; }

  std::span<const StateSpec> specs() const { return specs_; }

  //!< One entry per arena row. Generic steppers advance only StateRole::ode rows.
  std::span<const StateRole> roles() const { return roles_; }

  bool is_ode_row(int row) const
  {
    assert(0 <= row && row < next_row_);
    return roles_[static_cast<std::size_t>(row)] == StateRole::ode;
  }

  /**
   * Freeze the layout and allocate the arena — the ONLY allocation of state memory
   * (PC-1). May be called more than once (same archetype layout, different lane counts:
   * one arena per batch instance), but the layout cannot grow afterwards.
   */
  StateArena build(int n_lanes)
  {
    assert(next_row_ > 0 && "no states declared");
    frozen_ = true;
    return StateArena{ next_row_, n_lanes };
  }

private:
  std::vector<StateSpec> specs_{};
  std::vector<StateSlice> slices_{};
  std::vector<StateRole> roles_{};
  std::optional<StateSlice> q_ext_{};
  int next_row_{ 0 };
  bool frozen_{ false };
};

} // namespace slide::core
