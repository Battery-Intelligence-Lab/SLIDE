/**
 * @file Experiment.hpp
 * @brief PyBaMM-style experiment grammar and event-aligned core cycler.
 */

#pragma once

#include "EulerLegacy.hpp"
#include "ExponentialModal.hpp"

#include <span>
#include <string>
#include <vector>

namespace slide::core {

enum class ControlMode : unsigned char { current,
                                         voltage,
                                         power,
                                         rest,
                                         drive_cycle };
enum class Direction : signed char { charge = -1,
                                     none = 0,
                                     discharge = 1 };

struct ExperimentSegment
{
  ControlMode mode{ ControlMode::rest };
  Direction direction{ Direction::none };
  real_t value{}; //!< A, C-rate, V, or W according to mode
  bool value_is_c_rate{};
  real_t duration{};       //!< seconds; zero means event-terminated
  real_t voltage_limit{};  //!< V; zero means absent
  real_t current_cutoff{}; //!< A or C-rate magnitude; zero means absent
  bool cutoff_is_c_rate{};
  std::string drive_cycle{};
  std::string source{};
};

struct ParseDiagnostic
{
  std::size_t step{};
  std::size_t offset{};
  std::string message{};
};

struct Experiment
{
  std::vector<ExperimentSegment> segments{};

  [[nodiscard]] static slide::Status parse(
    std::span<const std::string> steps,
    Experiment &output,
    ParseDiagnostic &diagnostic);
};

enum class TerminationReason : unsigned char { event,
                                               limit,
                                               error,
                                               final_time };

struct ExperimentSolution
{
  std::vector<real_t> time{};
  std::vector<real_t> voltage{};
  std::vector<real_t> current{};
  TerminationReason reason{ TerminationReason::final_time };
  slide::Status status{ slide::Status::Success };
  std::size_t segment{};
};

struct DriveCycle
{
  std::string name{};
  std::vector<real_t> time{};
  std::vector<real_t> current{};
};

enum class CyclerIntegrator : unsigned char { euler_legacy,
                                              exponential };

class CyclerV2
{
public:
  [[nodiscard]] slide::Status configure(
    SpmBatch &batch,
    CyclerIntegrator integrator = CyclerIntegrator::exponential);
  [[nodiscard]] slide::Status registerDriveCycle(DriveCycle cycle);
  [[nodiscard]] slide::Status run(const Experiment &experiment,
                                  real_t sample_step,
                                  ExperimentSolution &output);

private:
  [[nodiscard]] slide::Status currentForVoltage(real_t target, real_t &current);
  [[nodiscard]] slide::Status currentForPower(real_t target_power,
                                              Direction direction,
                                              real_t &current);
  [[nodiscard]] slide::Status advance(real_t current, real_t time, real_t dt);
  [[nodiscard]] slide::Status voltageAt(real_t current, real_t &voltage);
  const DriveCycle *findDriveCycle(const std::string &name) const;
  real_t driveCurrent(const DriveCycle &cycle, real_t local_time) const;

  SpmBatch *batch_{};
  CyclerIntegrator integrator_{ CyclerIntegrator::exponential };
  EulerLegacy euler_{};
  ExponentialModal exponential_{};
  std::vector<real_t> density_{};
  std::vector<real_t> event_backup_{};
  std::vector<DriveCycle> drive_cycles_{};
};

} // namespace slide::core
