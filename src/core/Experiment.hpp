/**
 * @file Experiment.hpp
 * @brief PyBaMM-style experiment grammar and event-aligned core cycler.
 */

#pragma once

#include "EulerLegacy.hpp"
#include "ExponentialModal.hpp"

#include <functional>
#include <limits>
#include <span>
#include <string>
#include <vector>

namespace slide::core {

enum class ControlMode : unsigned char { current,
                                         voltage,
                                         power,
                                         rest,
                                         drive_cycle,
                                         custom_explicit,
                                         custom_implicit,
                                         custom_differential };
enum class Direction : signed char { charge = -1,
                                     none = 0,
                                     discharge = 1 };

/** Observable subset available to dependency-free custom controls/events. */
struct ExperimentVariables
{
  real_t time{};
  real_t local_time{};
  real_t voltage{};
  real_t current{};
  real_t power{};
};

using ExperimentFunction = std::function<real_t(const ExperimentVariables &)>;

struct CustomTermination
{
  std::string name{};
  ExperimentFunction indicator{}; //!< positive before event, zero/negative at event
};

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
  ExperimentFunction custom_control{};
  std::vector<CustomTermination> custom_terminations{};
  /** Seconds from the first scheduled step; a negative value means start immediately. */
  real_t scheduled_start{ -1.0 };
  /** Per-step recording/control period; a negative value selects run's default. */
  real_t sample_period{ -1.0 };
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
  std::vector<std::size_t> sample_segment{};
  TerminationReason reason{ TerminationReason::final_time };
  slide::Status status{ slide::Status::Success };
  std::size_t segment{};
  std::string termination_name{};
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
  [[nodiscard]] slide::Status registerDriveCycle(const DriveCycle &cycle);
  [[nodiscard]] slide::Status run(const Experiment &experiment,
                                  real_t sample_step,
                                  ExperimentSolution &output);

private:
  [[nodiscard]] slide::Status currentForVoltage(real_t target, real_t &current);
  [[nodiscard]] slide::Status currentForPower(real_t target_power,
                                              Direction direction,
                                              real_t &current);
  [[nodiscard]] slide::Status currentForCustom(
    const ExperimentSegment &segment,
    real_t time,
    real_t local_time,
    real_t &current);
  [[nodiscard]] slide::Status evaluateFunction(
    const ExperimentFunction &function,
    real_t time,
    real_t local_time,
    real_t voltage,
    real_t current,
    real_t &value) const;
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
  std::vector<real_t> event_derivative_backup_{};
  std::vector<real_t> run_state_backup_{};
  std::vector<real_t> run_derivative_backup_{};
  std::vector<DriveCycle> drive_cycles_{};
  bool in_run_{};
  bool run_snapshot_ready_{};
};

} // namespace slide::core
