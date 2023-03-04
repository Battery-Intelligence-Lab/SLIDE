/*
 * Status.hpp
 *
 * A small class for error codes.
 *  Created on: 21 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <cstdint>

namespace slide {
enum class Status : int_fast8_t //!< -128 to 127 = 1 byte.
{
  ReachedCurrentLimit = -6,
  ReachedVoltageLimit,
  ReachedTimeLimit,
  ReachedSmallCurrent, //!< Not sure if this is an error?
  Invalid_Vset,
  SafeVoltage,
  //!< Upper part is from cycler!
  Success, //!< = 0
  Vmin_violation,
  Vmax_violation,
  VMIN_violation,
  VMAX_violation,
  VMINsafety_violation, //!< 0.99 VMIN
  VMAXsafety_violation, //!< 1.01 VMAX
  V_not_calculated,
  SOC_limits_violation,
  Invalid_states,
  Invalid_SUs,
  ParallelUnit_failed, //!< throw 14.
  RedistributeCurrent_failed,
  timeStep_CC_failed,
  setVoltage_not_defined,
  Unknown_problem = 127,

  //!< Auxillary definitions:
  Critical = VMIN_violation, //!< Non-critical status definition.
  NotSafe = VMINsafety_violation

};

const char *getStatusMessage(Status status);
bool inline isStatusSuccessful(Status status) { return status == Status::Success; }
bool inline isStatusFailed(Status status) { return !isStatusSuccessful(status); }
bool inline isStatusOK(Status status) { return (status < Status::Critical); }
bool inline isStatusBad(Status status) { return !isStatusOK(status); }
bool inline isStatusWarning(Status status) { return ((status != Status::Success) && isStatusOK(status)); } //!< Not successful but voltage is in limits.
bool inline isStatusSafe(Status status) { return status < Status::NotSafe; }
bool inline isLimitsReached(Status status) { return status < Status::ReachedTimeLimit; }

bool inline isStatusVoltageLimitsViolated(Status status) { return (Status::Success < status && status <= Status::VMAXsafety_violation); }

bool inline isCCLimitReached(Status status) { return status == Status::ReachedVoltageLimit || status == Status::ReachedTimeLimit; }
bool inline isCVLimitReached(Status status) { return status == Status::ReachedCurrentLimit || status == Status::ReachedTimeLimit; }
bool inline isCurrentLimitReached(Status status) { return status == Status::ReachedCurrentLimit; }
bool inline isVoltageLimitReached(Status status) { return status == Status::ReachedVoltageLimit; }


inline const char *getStatusMessage(Status status)
{
  switch (status) {
  case Status::VMIN_violation:
    return "VMIN is violated!";
  case Status::Vmin_violation:
    return "Vmin is violated!";
  case Status::Success:
    return "Success! Yay!";
  case Status::Vmax_violation:
    return "Vmax is violated!";
  case Status::VMAX_violation:
    return "VMAX is violated!";
  case Status::V_not_calculated:
    return "V could not be calculated at all!";
  case Status::SOC_limits_violation:
    return "SOC_limits_violation!";
  case Status::Invalid_states:
    return "Invalid_states!";
  case Status::Invalid_SUs:
    return "Invalid_SUs!";
  case Status::ReachedCurrentLimit:
    return "ReachedCurrentLimit!";
  case Status::ReachedVoltageLimit:
    return "ReachedVoltageLimit!";
  case Status::ReachedTimeLimit:
    return "ReachedTimeLimit!";
  case Status::ReachedSmallCurrent:
    return "ReachedSmallCurrent!";
  case Status::Invalid_Vset:
    return "Invalid_Vset!";
  case Status::SafeVoltage:
    return "SafeVoltage!";
  case Status::ParallelUnit_failed:
    return "ParallelUnit_failed!";
  case Status::RedistributeCurrent_failed:
    return "RedistributeCurrent_failed!";
  case Status::timeStep_CC_failed:
    return "timeStep_CC_failed!";
  case Status::setVoltage_not_defined:
    return "setVoltage is not defined for this class!";
  case Status::Unknown_problem:
    return "Unknown problem!";
  default:
    return "Unkown status?!";
  }
}
} // namespace slide