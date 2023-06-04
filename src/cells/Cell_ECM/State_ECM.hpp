/**
 * @file State_ECM.hpp
 * @brief State_ECM class definition
 * @author Jorn Reniers, Volkan Kumtepeli
 * @date 14 Apr 2022
 */

#pragma once

#include "../../types/State.hpp"

namespace slide {

/**
 * @brief State_ECM template class for holding the states of the ECM model.
 * @tparam N_RC Number of parallel resistance-capacitance elements (default: 1).
 */
template <size_t N_RC = 1>
struct State_ECM : public State<3 + N_RC>
{

  /**
   * @brief Enumeration of indices for state variables.
   */
  enum Index : size_t {
    i_T,      //!< Cell temperature [K].
    i_SOC,    //!< State of charge [0-1].
    i_I,      //!< Current, [A], + for discharge, - for charge.
    i_Ir,     //!< Current through the parallel resistance, [I].
    N_states, //!< Total number of states (not recommended for usage, use .size() instead).
  };

  /**
   * @brief Get the current value (const method).
   * @return The current value [A], + for discharge, - for charge.
   */
  inline auto I() const { return (*this)[i_I]; } //!< Current, [A], + for discharge, - for charge

  /**
   * @brief Get the current value (non-const method).
   * @return The current value [A], + for discharge, - for charge.
   */
  inline auto &I() { return (*this)[i_I]; }

  /**
   * @brief Get the current through the parallel resistance (non-const method).
   * @param i Index of the parallel resistance (default: 0).
   * @return The current through the parallel resistance [I].
   */
  inline auto &Ir(size_t i = 0) { return (*this)[i_I + i + 1]; }

  /**
   * @brief Get the state of charge value (non-const method).
   * @return The state of charge value [0-1].
   */
  inline auto &SOC() { return (*this)[i_SOC]; }

  /**
   * @brief Get the temperature value (non-const method).
   * @return The temperature value [K].
   */
  inline auto &T() { return (*this)[i_T]; }
};

using State_Bucket = State_ECM<0>; //!< Bucket Cell with no RC pairs.

} // namespace slide