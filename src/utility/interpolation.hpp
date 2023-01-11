/*
 * interpolation.hpp
 *
 * Groups functions to do linear interpolation
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
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
template <typename Tx, typename Ty>
auto linInt_noexcept(bool bound, Tx &xdat, Ty &ydat, int nin, double x, bool is_fixed = false)
{
  /*
   * function for linear interpolation with the data points provided as two arrays
   *
   * IN
   * bound	boolean deciding what to do if the value of x is out of range of xdat
   * 			if true, an error is thrown
   * 			if false, the value closest to x is returned
   * xdat 	x data points in strictly increasing order
   * ydat 	y data points
   * nin 		number of data points
   * x 		x point at which value is needed
   * is_fixed If the difference between values are fixed.
   *
   * OUT
   * y 		y value corresponding to x
   * status     0 if successful, 1 if x>first and -1 if x < last
   */

  double yy{ 0.0 }; //!< Some programs depend on 0.0 initial condition when status != 0, do not change.
  int status = 0;   //!< Set the status as inverse of bound. So that first two branches of if are invalid if bound is true.
  //!< check that x is within the limits of the data points

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
    yy = ydat[0]; //!< if x is below the minimum value, return the y value of the first data point
  else if (x >= xdat[nin - 1])
    yy = ydat[nin - 1]; //!< if x is above the maximum value, return the y value of the last data point
  else {
    //!< scan the data points
    int i_low{ 0 };
    //!< For fixed step no need to iterate:

    if (is_fixed) {
      double dt = xdat[1] - xdat[0];
      i_low = static_cast<int>((x - xdat[0]) / dt) + 1;
    } else {
      //!< binary search algorithm:
      //!< i_low will be the first index which compares greater than x;
      //!< const auto it = std::find_if(std::begin(xdat), std::begin(xdat) + nin, [x](double element) { return (x < element); }); -> Linear search if needed.
      const auto it = std::lower_bound(xdat.begin(), xdat.begin() + nin, x);
      i_low = static_cast<int>(it - xdat.begin());
    }

    const double xr = xdat[i_low]; //!< then that point is the point 'to the right' of x
    const double yr = ydat[i_low];
    const double xl = xdat[i_low - 1]; //!< while the previous point is the point 'to the left' of x
    const double yl = ydat[i_low - 1];
    yy = yl + (yr - yl) * (x - xl) / (xr - xl);
  }

  return std::make_pair(yy, status);
}

template <typename Tx, typename Ty>
double linInt(bool verbose, bool bound, Tx &xdat, Ty &ydat, int nin, double x, bool is_fixed = false)
{
  /*
   * function for linear interpolation with the data points provided as two arrays
   *
   * IN
   * verbose 	if false, no error message is written (but the error is still thrown)
   * 			if true, an error message is written
   * bound	boolean deciding what to do if the value of x is out of range of xdat
   * 			if true, an error is thrown
   * 			if false, the value closest to x is returned
   * xdat 	x data points in strictly increasing order
   * ydat 	y data points
   * nin 		number of data points
   * x 		x point at which value is needed
   * is_fixed If the difference between values are fixed.
   *
   * OUT
   * y 		y value corresponding to x
   *
   * THROWS
   * 1 		if bound = true && if x is out of bounds, i.e. smaller than the smallest value of xdat or larger than the largest value of xdat
   * 				i.e. bound AND (x < xdat [0] OR x > xdat[end])
   */

  auto [yy, status] = linInt_noexcept(bound, xdat, ydat, nin, x, is_fixed);

  if (status) {
    if (verbose)
      std::cerr << "ERROR in Interpolation::linInt: x is out of bounds. x = " << x << " while xmin = " << xdat[0] << " and xmax is " << xdat[nin - 1] << ".\n";
    throw 1;
  }

  return yy;
}

} // namespace slide
