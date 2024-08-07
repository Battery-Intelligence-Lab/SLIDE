/**
 * @file interpolation.hpp
 * @brief Groups functions to do linear interpolation.
 * @author Jorn Reniers
 * @author Volkan Kumtepeli
 * @date 09 Aug 2021
 */

#pragma once

#include "../settings/settings.hpp"

#include <vector>
#include <array>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <utility>
#include <cmath>
#include <cstdlib>

namespace slide {

/**
 * @brief Linear interpolation function with exception handling for out-of-bound values.
 *
 * @tparam Tx Type of x data points.
 * @tparam Ty Type of y data points.
 * @param bound Boolean deciding what to do if the value of x is out of range of xdat.
 *              If true, an error is thrown. If false, the value closest to x is returned.
 * @param xdat x data points in strictly increasing order.
 * @param ydat y data points.
 * @param nin Number of data points.
 * @param x x point at which value is needed.
 * @param is_fixed If the difference between values is fixed.
 * @return std::pair<double, int> A pair containing the interpolated y value and the status (0 if successful, 1 if x > first, and -1 if x < last).
 */
template <typename Tx, typename Ty>
auto linInt_noexcept(bool bound, Tx &xdat, Ty &ydat, int nin, double x, bool is_fixed = false)
{
  double yy{ 0.0 }; //!< Some programs depend on 0.0 initial condition when status != 0, do not change.
  int status = 0;   //!< Set the status as inverse of bound. So that first two branches of if are invalid if bound is true.

  //!< Check that x is within the limits of the data points.
  if (bound) {
    if (x < xdat[0]) {
      status = -1;
      return std::make_pair(yy, status);
    } else if (x > xdat[nin - 1]) {
      status = 1;
      return std::make_pair(yy, status);
    }
  }

  if (x <= xdat[0])
    yy = ydat[0]; //!< If x is below the minimum value, return the y value of the first data point.
  else if (x >= xdat[nin - 1])
    yy = ydat[nin - 1]; //!< If x is above the maximum value, return the y value of the last data point.
  else {
    //!< Scan the data points.
    int i_low{ 0 };

    //!< For fixed step, no need to iterate:
    if (is_fixed) {
      double dt = xdat[1] - xdat[0];
      i_low = static_cast<int>((x - xdat[0]) / dt) + 1;
    } else {
      //!< Binary search algorithm:
      const auto it = std::lower_bound(xdat.begin(), xdat.begin() + nin, x);
      i_low = static_cast<int>(it - xdat.begin());
    }

    const double xr = xdat[i_low]; //!< Then that point is the point 'to the right' of x.
    const double yr = ydat[i_low];
    const double xl = xdat[i_low - 1]; //!< While the previous point is the point 'to the left' of x.
    const double yl = ydat[i_low - 1];
    yy = yl + (yr - yl) * (x - xl) / (xr - xl);
  }

  return std::make_pair(yy, status);
}

/**
 * @brief Linear interpolation function with verbosity and exception handling for out-of-bound values.
 *
 * @tparam Tx Type of x data points.
 * @tparam Ty Type of y data points.
 * @param verbose If false, no error message is written (but the error is still thrown). If true, an error message is written.
 * @param bound Boolean deciding what to do if the value of x is out of range of xdat.
 *              If true, an error is thrown. If false, the value closest to x is returned.
 * @param xdat x data points in strictly increasing order.
 * @param ydat y data points.
 * @param nin Number of data points.
 * @param x x point at which value is needed.
 * @param is_fixed If the difference between values is fixed.
 * @return double The interpolated y value.
 * @throws int If bound is true and x is out of bounds, i.e., smaller than the smallest value of xdat or larger than the largest value of xdat.
 */
template <typename Tx, typename Ty>
double linInt(bool verbose, bool bound, Tx &xdat, Ty &ydat, int nin, double x, bool is_fixed = false)
{
  auto [yy, status] = linInt_noexcept(bound, xdat, ydat, nin, x, is_fixed);

  if (status) {
    if (verbose)
      std::cerr << "ERROR in Interpolation::linInt: x is out of bounds. x = " << x << " while xmin = " << xdat[0] << " and xmax is " << xdat[nin - 1] << ".\n";
    throw 1;
  }

  return yy;
}

} // namespace slide
