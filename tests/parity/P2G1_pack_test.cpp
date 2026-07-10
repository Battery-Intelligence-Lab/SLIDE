/**
 * @file P2G1_pack_test.cpp
 * @brief P2-G1 3s2p compiled-pack parity against legacy Module_s/Module_p.
 */

#include "../../src/core/PackStepper.hpp"
#include "../support/KokamSpmFixture.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>

using namespace slide;

namespace {

struct Drift
{
  double voltage{};
  double current{};
  double state_scaled{};
  double worst_legacy{};
  double worst_core{};
  bool state_ok{ true };

  void state(double legacy, double value)
  {
    const double error = std::abs(value - legacy);
    const double tolerance = 1e-15 + 1e-12 * std::max(std::abs(value), std::abs(legacy));
    if (error / tolerance > state_scaled) {
      state_scaled = error / tolerance;
      worst_legacy = legacy;
      worst_core = value;
    }
    state_ok = state_ok && error <= tolerance;
  }
};

void compareCell(const Cell_SPM &legacy, const core::SpmBatch &batch,
                 int lane, Drift &drift)
{
  auto &legacy_state = const_cast<Cell_SPM &>(legacy).getStateObj();
  const auto &state = batch.state();
  const auto &layout = batch.layout();
  for (const auto legacy_domain : { neg, pos }) {
    const auto domain = legacy_domain == neg ? core::Domain::neg : core::Domain::pos;
    const auto d = core::domain_index(domain);
    for (int mode = 0; mode < static_cast<int>(settings::nch); ++mode)
      drift.state(legacy_state.z(static_cast<std::size_t>(mode), legacy_domain),
                  state.at(layout.spm.z[d], mode, lane));
    drift.state(legacy_state.e(legacy_domain), state.at(layout.spm.active_fraction[d], 0, lane));
    drift.state(legacy_state.D(legacy_domain), state.at(layout.spm.diffusion_coefficient[d], 0, lane));
    drift.state(legacy_state.thick(legacy_domain), state.at(layout.spm.electrode_thickness[d], 0, lane));
    drift.state(legacy_state.a(legacy_domain), state.at(layout.spm.specific_surface_area[d], 0, lane));
    drift.state(legacy_state.rDC(legacy_domain), state.at(layout.spm.specific_resistance[d], 0, lane));
  }
  drift.state(legacy_state.T(), state.at(layout.spm.temperature, 0, lane));
  drift.state(legacy_state.delta(), state.at(layout.spm.sei_thickness, 0, lane));
  drift.state(legacy_state.LLI(), state.at(layout.spm.lost_lithium, 0, lane));
  drift.state(legacy_state.CS(), state.at(layout.spm.crack_surface, 0, lane));
  drift.state(legacy_state.rDCcc(), state.at(layout.spm.current_collector_resistance, 0, lane));
  drift.state(legacy_state.time(), state.at(layout.elapsed_time, 0, lane));
  drift.state(legacy_state.Ah(), state.at(layout.charge_throughput, 0, lane));
  drift.state(legacy_state.Wh(), state.at(layout.energy_throughput, 0, lane));
}

} // namespace

TEST_CASE("P2-G1 compiled 3s2p Kokam trajectory matches legacy modules",
          "[parity][core][pack][P2-G1]")
{
  constexpr int series_count = 3;
  constexpr int parallel_count = 2;
  constexpr int lanes = series_count * parallel_count;
  constexpr double pack_current = 32.0;
  constexpr double dt = 1.0;
  constexpr int steps = 300;

  std::array<Deep_ptr<StorageUnit>, series_count> parallel_modules;
  double initial_temperature{};
  double reference_temperature{};
  for (int s = 0; s < series_count; ++s) {
    std::array<Deep_ptr<StorageUnit>, parallel_count> cells;
    for (auto &storage : cells) {
      auto cell = make<Cell_SPM>();
      REQUIRE(test_support::initialize_legacy_kokam(*cell, 0.55, 0.0)
              == Status::Success);
      if (s == 0 && &storage == &cells.front()) {
        double environment_temperature{};
        cell->getTemperatures(&environment_temperature, &reference_temperature);
        initial_temperature = cell->T();
      }
      storage = std::move(cell);
    }
    auto module = make<Module_p>("p", 298.0, false, false, parallel_count, 1, 1);
    module->setSUs(cells, false, false);
    parallel_modules[static_cast<std::size_t>(s)] = std::move(module);
  }
  auto legacy = make<Module_s>("s", 298.0, false, false, lanes, 1, 1);
  legacy->setSUs(parallel_modules, false, false);
  legacy->setBlockDegAndTherm(true);

  core::SpmBatch batch;
  REQUIRE(core::buildSpmBatch(
            test_support::make_legacy_kokam_input(
              0.55, initial_temperature, reference_temperature),
            {},
            lanes,
            batch)
          == Status::Success);
  core::CompiledPackTopology topology;
  REQUIRE(core::compilePackDescription(
            { .root = core::series(series_count,
                                   core::parallel(parallel_count,
                                                  core::cell({ .archetype = "spm" }))) },
            topology)
          == Status::Success);
  std::array<core::SpmBatch *, 1> batches{ &batch };
  core::PackStepper core_pack;
  REQUIRE(core_pack.configure(topology, batches) == Status::Success);

  Drift drift;
  for (int step = 0; step < steps; ++step) {
    REQUIRE(legacy->setCurrent(pack_current, false, false) == Status::Success);
    REQUIRE(core_pack.step(pack_current, step * dt, dt) == Status::Success);
    drift.voltage = std::max(drift.voltage,
                             std::abs(core_pack.solution().terminal_voltage - legacy->V()));

    for (int s = 0; s < series_count; ++s) {
      auto *module = dynamic_cast<Module_p *>(legacy->getSUs()[static_cast<std::size_t>(s)].get());
      REQUIRE(module != nullptr);
      for (int p = 0; p < parallel_count; ++p) {
        auto *cell = dynamic_cast<Cell_SPM *>(module->getSUs()[static_cast<std::size_t>(p)].get());
        REQUIRE(cell != nullptr);
        const int lane = s * parallel_count + p;
        drift.current = std::max(drift.current,
                                 std::abs(core_pack.solution().cell_current[static_cast<std::size_t>(lane)]
                                          - cell->I()));
      }
    }

    legacy->timeStep_CC(dt);
    for (int s = 0; s < series_count; ++s) {
      auto *module = dynamic_cast<Module_p *>(legacy->getSUs()[static_cast<std::size_t>(s)].get());
      for (int p = 0; p < parallel_count; ++p) {
        auto *cell = dynamic_cast<Cell_SPM *>(module->getSUs()[static_cast<std::size_t>(p)].get());
        compareCell(*cell, batch, s * parallel_count + p, drift);
      }
    }
  }

  std::printf("P2-G1 3s2p: dV=%.17g dI=%.17g state_scaled=%.17g worst=(%.17g,%.17g) factors=%d\n",
              drift.voltage,
              drift.current,
              drift.state_scaled,
              drift.worst_legacy,
              drift.worst_core,
              core_pack.solver().workspace().numericFactorizations());
  CAPTURE(drift.voltage, drift.current, drift.state_scaled);
  // PLAN §5.2 registers voltage digit parity at 1e-12; Phase 9A explicitly
  // adopts the same scale for branch current. Do not add silent headroom.
  REQUIRE(drift.voltage <= 1e-12);
  REQUIRE(drift.current <= 1e-12);
  REQUIRE(drift.state_ok);
}
