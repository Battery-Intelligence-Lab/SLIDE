/**
 * @file SpmBatchLayout.hpp
 * @brief Which arena rows an SPM batch owns — the layout `SpmBatch::layout()` returns.
 *
 * Owns: the `StateSlice` descriptors for the SPM, thermal, and stress-history row groups,
 * plus the aggregate an `SpmBatch` carries. Implements PLAN.md §3.1 (SoA arena) and MC-5:
 * the layout is public vocabulary, so it must not live inside the kernel headers that fill it.
 * Cold: constructed once per batch, read by recorders and callers.
 * @surface support
 */

#pragma once

#include "CellDesign.hpp"
#include "StateArena.hpp"

namespace slide::core {

/** Rows of the SPM composition itself (concentrations, temperature, ageing states). */
struct SpmStateLayout
{
  PerDomain<StateSlice> z{};
  StateSlice temperature{};
  StateSlice sei_thickness{};
  StateSlice lost_lithium{};
  StateSlice crack_surface{};
  StateSlice plated_lithium_thickness{};
  PerDomain<StateSlice> active_fraction{};
  PerDomain<StateSlice> diffusion_coefficient{};
  PerDomain<StateSlice> electrode_thickness{};
  PerDomain<StateSlice> specific_surface_area{};
  PerDomain<StateSlice> specific_resistance{};
  StateSlice current_collector_resistance{};
};

/** Rows the lumped thermal component reads and writes. */
struct ThermalLumpedLayout
{
  StateSlice temperature{};
  StateSlice external_heat_flow{};    //!< q_ext, positive into the cell [W]
  StateSlice generated_heat_energy{}; //!< integral of internal heat generation [J]
  StateSlice thermal_elapsed_time{};  //!< time integrated by the thermal component [s]
};

/** Previous-step stress values consumed by crack/LAM rate laws; all checkpointed arena rows. */
struct StressHistoryLayout
{
  PerDomain<StateSlice> previous_dai{};
  StateSlice previous_laresgoiti_negative{};
  StateSlice interval{};
};

/** Every row group an assembled SPM batch owns. */
struct SpmBatchLayout
{
  SpmStateLayout spm{};
  ThermalLumpedLayout thermal{};
  StressHistoryLayout stress_history{};
  StateSlice elapsed_time{};
  StateSlice charge_throughput{};
  StateSlice energy_throughput{};
};

} // namespace slide::core
