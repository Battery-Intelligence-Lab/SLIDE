/*
 * DataStorage.hpp
 *
 * This class is created to generic interface to store data.
 *  Created on: 28 Aug 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */


#pragma once

struct DataStorage
{
  enum Index : size_t //!< Index variables for:
  {
    i_I,
    i_SOC,
    i_T, //!< cell temperature [K]
    N_states,
  };
};
