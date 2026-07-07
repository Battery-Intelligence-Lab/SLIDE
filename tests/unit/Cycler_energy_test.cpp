/**
 * @file Cycler_energy_test.cpp
 * @brief Registered test for bug B3 (PLAN.md 2.4): Cycler::CC energy throughput.
 *
 * Cycler::CC accumulated Wh throughput with a single end-of-step voltage sample
 * (`th.Wh() += |I|*dt/3600 * vi`). The in-source TODO says it should use the
 * trapezoid (v_before + v_after)/2. This test drives a cell whose terminal voltage
 * ramps LINEARLY from 4.0 V to 3.0 V over 3600 s at a constant |I| = 1 A, dt = 1 s.
 *
 * REGISTERED PREDICTION (written before running):
 *   - Energy = |I| * integral(V dt) / 3600 = 1 * (avg V 3.5) * 3600 / 3600 = 3.5 Wh.
 *     The composite trapezoid rule is EXACT for a linear integrand, so the fixed
 *     code must give 3.5 Wh. Registered pass band: |Wh - 3.5| <= 1e-6.
 *   - Charge throughput Ah = |I| * 3600 / 3600 = 1.0 Ah (unaffected by the fix).
 *   - For reference, a correct end-of-step Riemann sum would give
 *     sum_{k=1..3600}(1/3600)(4 - k/3600) = 3.49986... Wh (off by ~1.4e-4, i.e. it
 *     would FAIL the 1e-6 band). The unfixed code is in fact worse still: its `vi`
 *     is never assigned by Cycler::setCurrent, so pre-fix Wh == 0.
 */

#include "../../src/slide.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <span>
#include <vector>
#include <string>

using Catch::Matchers::WithinAbs;

namespace {

//!< Minimal StorageUnit whose terminal voltage decreases linearly with elapsed
//!< time, so integral(V dt) has a known closed form. Only the members Cycler::CC
//!< actually touches carry behaviour; the rest are inert stubs.
class RampCell : public slide::StorageUnit
{
  double Icell{ 0.0 };
  double t{ 0.0 }; //!< elapsed time [s]

  static constexpr double Vstart = 4.0, Vend = 3.0, Tramp = 3600.0;

public:
  RampCell() : slide::StorageUnit("rampcell") {}

  double Vof(double tt) const { return Vstart + (Vend - Vstart) * (tt / Tramp); }

  //!< --- behaviour used by Cycler::CC / Cycler::CV ---
  double V() override { return Vof(t); }
  double getOCV() override { return Vof(t); }
  double I() const override { return Icell; }
  slide::Status setCurrent(double Inew, bool = true, bool = true) override
  {
    Icell = Inew;
    return slide::Status::Success;
  }
  //!< Inert voltage regulation: Cycler::CV only needs setVoltage to succeed; the
  //!< terminal voltage keeps following the time ramp so integral(V dt) stays known.
  slide::Status setVoltage(double, bool = true, bool = true) override
  {
    return slide::Status::Success;
  }
  void timeStep_CC(double dt, int steps = 1) override { t += dt * steps; }

  //!< --- inert stubs ---
  double Cap() const override { return 1.0; }
  double getRtot() override { return 0.0; }
  size_t getNcells() override { return 1; }
  void getStates(std::vector<double> &) override {}
  slide::Status setStates(std::span<const double>, int &, bool = true, bool = true) override
  {
    return slide::Status::Success;
  }
  slide::Status checkVoltage(double &v, bool) noexcept override
  {
    v = Vof(t);
    return slide::Status::Success;
  }
  double getVhigh() override { return V(); }
  double getVlow() override { return V(); }
  double Vmin() const override { return 2.5; }
  double VMIN() const override { return 2.0; }
  double Vmax() const override { return 4.2; }
  double VMAX() const override { return 4.3; }
  double T() override { return 298.15; }
  double getThotSpot() override { return 298.15; }
  double getThermalSurface() override { return 1.0; }
  double thermalModel(int, double[], double[], double[], double) override { return 298.15; }
  void setT(double) override {}
  bool validStates(bool = true) override { return true; }
  slide::StorageUnit *copy() override { return new RampCell(*this); }
  void storeData() override {}
  void writeData(const std::string &) override {}
};

} // namespace

TEST_CASE("Cycler::CC energy throughput uses the trapezoid rule", "[Cycler][B3]")
{
  RampCell cell;
  slide::Cycler cyc(&cell, "B3_energy");

  slide::ThroughputData th{};
  const double I = 1.0;      //!< discharge, |I| = 1 A
  const double vlim = 0.0;   //!< unreachable low limit -> terminate on time only
  const double tlim = 3600.0;
  const double dt = 1.0;
  const int ndt_data = 0;    //!< no data storage

  const auto succ = cyc.CC(I, vlim, tlim, dt, ndt_data, th);

  REQUIRE(succ == slide::Status::ReachedTimeLimit);
  REQUIRE_THAT(th.time(), WithinAbs(3600.0, 1e-9));
  REQUIRE_THAT(th.Ah(), WithinAbs(1.0, 1e-9));
  //!< Registered band (linear ramp -> trapezoid is exact):
  REQUIRE_THAT(th.Wh(), WithinAbs(3.5, 1e-6));
}

/**
 * Registered test for P0-C6 (PLAN.md 2.5): Cycler::CV throughput.
 *
 * Cycler::CV accumulated Wh with the END-of-step voltage only (`th.Wh() += dAh * vi`
 * with vi sampled after timeStep_CC) and never accumulated th.time() at all.
 * Same 4.0 -> 3.0 V linear ramp over 3600 s at |I| = 1 A, dt = 1 s (the stub's
 * setVoltage is inert, so CV reduces to time integration of a known V(t)).
 *
 * REGISTERED PREDICTION (written before running):
 *   - time = 3600 s (pre-fix: 0 s, since CV never touched th.time() -> FAILS).
 *   - Ah   = 1.0 (I is held constant during each step, |I|*dt is exact pre- and post-fix).
 *   - Wh   = 3.5 +- 1e-6 (trapezoid exact on a linear ramp). The pre-fix end-of-step
 *     Riemann sum gives sum_{k=1..3600}(1/3600)(4 - k/3600) = 3.49986... Wh,
 *     which MISSES the band by ~1.4e-4 -> FAILS.
 */
TEST_CASE("Cycler::CV throughput uses the trapezoid rule and counts time", "[Cycler][P0-C6]")
{
  RampCell cell;
  cell.setCurrent(1.0); //!< CV loop reads su->I(); keep |I| = 1 A throughout

  slide::Cycler cyc(&cell, "P0C6_energy");

  slide::ThroughputData th{};
  const double Vset = 3.6;  //!< inert for the stub; any value inside the ramp
  const double Ilim = 0.5;  //!< |I| stays 1 A -> current limit never reached
  const double tlim = 3600.0;
  const double dt = 1.0;
  const int ndt_data = 0;

  const auto succ = cyc.CV(Vset, Ilim, tlim, dt, ndt_data, th);

  REQUIRE(succ == slide::Status::ReachedTimeLimit);
  REQUIRE_THAT(th.time(), WithinAbs(3600.0, 1e-9));
  REQUIRE_THAT(th.Ah(), WithinAbs(1.0, 1e-9));
  //!< Registered band (linear ramp -> trapezoid is exact):
  REQUIRE_THAT(th.Wh(), WithinAbs(3.5, 1e-6));
}
