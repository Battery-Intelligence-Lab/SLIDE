/**
 * @file cell_limits.hpp
 * @brief  Limits for cell
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 07 Feb 2022
 */

#pragma once

#include "../settings/settings.hpp"

namespace slide {
struct CellLimits
{
  using data_type = double;
  data_type Vmin{ 2.7 }, Vmax{ 4.2 };                                     //!< minimum / maximum voltage  [V]
  data_type VMIN{ 2.0 }, VMAX{ 4.3 };                                     //!< lower   / upper safety cut-off limit [V]
  data_type Tmin{ settings::Tmin_Cell_K }, Tmax{ settings::Tmax_Cell_K }; //!< minimum / maximum temperature;
};

inline constexpr CellLimits defaultCellLimits;
} // namespace slide