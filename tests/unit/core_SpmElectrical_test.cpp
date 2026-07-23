/**
 * @file core_SpmElectrical_test.cpp
 * @brief Phase-1 parity gate for kinetics, OCV, resistance, voltage and heat observables.
 *
 * REGISTERED BEFORE FIRST RUN: for {0 A at 298.15 K, +20 A at 298.15 K, -10 A at 310 K},
 * core OCV, total resistance and terminal voltage match legacy Cell_SPM within 1e-12 absolute.
 * Nonzero-current cases must produce nonzero electrode overpotentials and finite heat terms.
 */

#include "../../src/slide.hpp"
#include "../../src/core/SpmObservables.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <bit>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <span>

using namespace slide;

namespace {

constexpr std::size_t core_index(slide::Domain domain)
{
  return core::domain_index(domain == pos ? core::Domain::pos : core::Domain::neg);
}

double runtime_real_from_bits(std::uint64_t bits)
{
  volatile std::uint64_t opaque_bits = bits;
  const std::uint64_t copied_bits = opaque_bits;
  return std::bit_cast<double>(copied_bits);
}

template <int NCH>
core::SpmElectricalParams<NCH> make_params(Model_SPM<NCH> &model, OCVcurves &curves,
                                           double reference_temperature)
{
  core::SpmElectricalParams<NCH> p;
  p.reference_temperature = reference_temperature;
  p.electrode_area = 0.1 * 0.2 * 31;
  p.sei_resistivity_area = 2037.4 * 50.0;

  auto &concentration = p.concentration;
  concentration.T_ref = reference_temperature;
  concentration.D_T[core_index(pos)] = 29000.0;
  concentration.D_T[core_index(neg)] = 35000.0 / 5.0;
  for (auto legacy_domain : { pos, neg }) {
    const auto d = core_index(legacy_domain);
    concentration.R[d] = legacy_domain == pos ? model.Rp : model.Rn;
    for (int node = 0; node < NCH + 1; ++node) {
      concentration.Dout[d][node] = model.D[legacy_domain](node);
      for (int mode = 0; mode < NCH; ++mode)
        concentration.C[d][node][mode] = model.C[legacy_domain](node, mode);
    }
  }
  for (int node = 0; node < NCH + 1; ++node)
    concentration.Cc[node] = model.Cc(node);
  concentration.cc_coeff = model.cc_coeff;

  auto &negative = core::domain_value(p.electrode, core::Domain::neg);
  negative.cs_max = 30555.0;
  negative.reaction_rate_ref = 1.764e-11;
  negative.reaction_activation = 20000.0;
  auto &positive = core::domain_value(p.electrode, core::Domain::pos);
  positive.cs_max = 51385.0;
  positive.reaction_rate_ref = 5e-11;
  positive.reaction_activation = 58000.0;

  REQUIRE(core::domain_value(p.electrode_ocv, core::Domain::neg)
            .build(curves.OCV_neg.x, curves.OCV_neg.y)
          == Status::Success);
  REQUIRE(core::domain_value(p.electrode_ocv, core::Domain::pos)
            .build(curves.OCV_pos.x, curves.OCV_pos.y)
          == Status::Success);
  REQUIRE(p.total_entropic_coefficient.build(curves.dOCV_tot.x, curves.dOCV_tot.y)
          == Status::Success);
  REQUIRE(p.negative_entropic_coefficient.build(curves.dOCV_neg.x, curves.dOCV_neg.y)
          == Status::Success);
  return p;
}

template <int NCH>
void copy_state(State_SPM &legacy, const core::SpmStateLayout &layout,
                core::StateArena &arena)
{
  for (auto legacy_domain : { pos, neg }) {
    const auto d = core_index(legacy_domain);
    for (int mode = 0; mode < NCH; ++mode)
      arena.at(layout.z[d], mode, 0) = legacy.z(mode, legacy_domain);
    arena.at(layout.diffusion_coefficient[d], 0, 0) = legacy.D(legacy_domain);
    arena.at(layout.electrode_thickness[d], 0, 0) = legacy.thick(legacy_domain);
    arena.at(layout.specific_surface_area[d], 0, 0) = legacy.a(legacy_domain);
    arena.at(layout.specific_resistance[d], 0, 0) = legacy.rDC(legacy_domain);
  }
  arena.at(layout.temperature, 0, 0) = legacy.T();
  arena.at(layout.sei_thickness, 0, 0) = legacy.delta();
  arena.at(layout.current_collector_resistance, 0, 0) = legacy.rDCcc();
}

} // namespace

TEST_CASE("SPM observable scratch spans its complete named field extent",
          "[core][observables][scratch]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  constexpr int lanes = 3;
  constexpr std::size_t expected_values =
    static_cast<std::size_t>(lanes)
    * (2 * (static_cast<std::size_t>(NCH) + 2) + 2 * 6 + 9);

  core::SpmObservableScratch<NCH> scratch{ lanes };
  const auto output = scratch.view();
  const auto *first =
    output.concentration[core::domain_index(core::Domain::neg)].data();
  const auto *last = output.total_heat.data() + output.total_heat.size();

  REQUIRE(output.concentration[core::domain_index(core::Domain::neg)].size()
          == static_cast<std::size_t>(NCH + 2) * lanes);
  REQUIRE(output.total_heat.size() == static_cast<std::size_t>(lanes));
  REQUIRE(static_cast<std::size_t>(last - first) == expected_values);
}

TEST_CASE("SPM electrical observable stage matches legacy Cell_SPM", "[core][observables]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  struct Scenario
  {
    double current;
    double temperature;
  };
  constexpr std::array scenarios{
    Scenario{ 0.0, 298.15 }, Scenario{ 20.0, 298.15 }, Scenario{ -10.0, 310.0 }
  };

  for (const auto scenario : scenarios) {
    Cell_SPM cell;
    cell.setBlockDegAndTherm(true);
    cell.setT(scenario.temperature);
    REQUIRE(cell.setCurrent(scenario.current, false, false) == Status::Success);
    double environment{}, reference_temperature{};
    cell.getTemperatures(&environment, &reference_temperature);
    auto *model = Model_SPM<>::makeModel();
    auto curves = OCVcurves::makeOCVcurves(cellType::KokamNMC);
    const auto params = make_params<NCH>(*model, curves, reference_temperature);

    core::BatchBuilder builder;
    const auto layout = core::declareSpmState<NCH>(builder);
    auto arena = builder.build(1);
    copy_state<NCH>(cell.getStateObj(), layout, arena);

    const std::array<double, 1> iapp{ scenario.current / params.electrode_area };
    const core::StepCtx ctx{ .time = 0.0, .dt = 0.0, .i_app = iapp };
    const core::ConstBatchView state{ core::BatchShape::from(arena),
                                      std::span<const core::real_t>{ arena.raw() } };
    core::SpmObservableScratch<NCH> scratch{ 1 };
    auto output = scratch.view();
    REQUIRE(core::computeSpmObservables(params, state, layout, ctx, output) == Status::Success);

    const double legacy_ocv = cell.getOCV();
    const double legacy_resistance = cell.getRtot();
    const double legacy_voltage = cell.V();
    const double ocv_error = std::abs(output.open_circuit_voltage[0] - legacy_ocv);
    const double resistance_error = std::abs(output.resistance[0] - legacy_resistance);
    const double voltage_error = std::abs(output.terminal_voltage[0] - legacy_voltage);
    std::printf("SpmElectrical I=%g T=%g: dOCV=%.3e dR=%.3e dV=%.3e Q=%.9g\n",
                scenario.current,
                scenario.temperature,
                ocv_error,
                resistance_error,
                voltage_error,
                output.total_heat[0]);

    CAPTURE(scenario.current, scenario.temperature, legacy_ocv, legacy_resistance, legacy_voltage);
    REQUIRE(ocv_error <= 1e-12);
    REQUIRE(resistance_error <= 1e-12);
    REQUIRE(voltage_error <= 1e-12);
    REQUIRE(std::isfinite(output.total_heat[0]));
    if (scenario.current != 0.0) {
      REQUIRE(core::domain_value(output.overpotential, core::Domain::neg)[0] != 0.0);
      REQUIRE(core::domain_value(output.overpotential, core::Domain::pos)[0] != 0.0);
    }
  }
}

TEST_CASE("SPM electrical observable stage returns Status for invalid concentrations",
          "[core][observables]")
{
  constexpr int NCH = static_cast<int>(settings::nch);
  Cell_SPM cell;
  double environment{}, reference_temperature{};
  cell.getTemperatures(&environment, &reference_temperature);
  auto *model = Model_SPM<>::makeModel();
  auto curves = OCVcurves::makeOCVcurves(cellType::KokamNMC);
  const auto params = make_params<NCH>(*model, curves, reference_temperature);

  core::BatchBuilder builder;
  const auto layout = core::declareSpmState<NCH>(builder);
  auto arena = builder.build(1); // zero z => invalid zero surface concentration
  arena.at(layout.temperature, 0, 0) = 298.15;
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    arena.at(layout.electrode_thickness[d], 0, 0) = 1.0;
    arena.at(layout.specific_surface_area[d], 0, 0) = 1.0;
  }
  const std::array<double, 1> iapp{};
  const core::StepCtx ctx{ .i_app = iapp };
  const core::ConstBatchView state{ core::BatchShape::from(arena),
                                    std::span<const core::real_t>{ arena.raw() } };
  core::SpmObservableScratch<NCH> scratch{ 1 };
  REQUIRE(core::computeSpmObservables(params, state, layout, ctx, scratch.view())
          == Status::Invalid_states);

  copy_state<NCH>(cell.getStateObj(), layout, arena);
  for (const auto domain : core::domains) {
    const auto d = core::domain_index(domain);
    REQUIRE(arena.at(layout.diffusion_coefficient[d], 0, 0) > 0.0);
  }
  REQUIRE(core::computeSpmObservables(params, state, layout, ctx, scratch.view())
          == Status::Success);

  auto &modal_state =
    arena.at(layout.z[core::domain_index(core::Domain::neg)], 0, 0);
  const auto valid_modal_state = modal_state;
  constexpr std::array invalid_bits{
    UINT64_C(0x7ff8000000000000),
    UINT64_C(0x7ff0000000000000),
    UINT64_C(0xfff0000000000000),
  };
  for (const auto bits : invalid_bits) {
    modal_state = runtime_real_from_bits(bits);
    REQUIRE(core::computeSpmObservables(params, state, layout, ctx, scratch.view())
            == Status::Invalid_states);
  }
  modal_state = valid_modal_state;
}
