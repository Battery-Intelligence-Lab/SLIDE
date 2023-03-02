/*
 * XYdata.hpp
 *
 *  Created on: 07 Feb 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "FixedData.hpp"
#include "../utility/interpolation.hpp"
#include "../utility/io/read_CSVfiles.hpp"

#include <stdexcept>
#include <vector>
#include <array>
#include <functional>
#include <cmath>
#include <filesystem>
#include <type_traits>
#include <span>

namespace slide {

template <typename Tx>
bool check_is_fixed(Tx &xdat) //!< Checks if steps are fixed length. #TODO should not be here.
{
  const double dt = xdat[1] - xdat[0];
  const double tol = 0.01 * dt; //!< Tolerance.

  bool is_fixed = true;

  for (size_t i{ 1 }; i < (xdat.size() - 1); i++) {
    const double residual = std::abs(xdat[i + 1] - xdat[i] - dt);
    if (residual > tol) {
      is_fixed = false;
      break;
    }
  }

  return is_fixed;
}

template <typename Tx, typename Ty>
class XYdata
{
  bool is_fixed{ false };

public:
  Tx x;
  Ty y;

  XYdata() = default;
  explicit XYdata(size_t N) : x(N), y(N) {}
  XYdata(Tx &x_, Ty &y_) : x(x_), y(y_) { check_is_fixed(); }

  //!< XYdata(FixedData x, Ty y) : is_fixed(true), x(x), y(y) {} #TODO this should be on but error in GCC

  void reserve(int n) { x.reserve(n), y.reserve(n); }
  void clear() { x.clear(), y.clear(); }

  void resize(size_t n) { x.resize(n), y.resize(n); }

  double interp(double x_i, bool print = false, bool bound = true) // #TODO cannot put const here.
  {
    return linInt(print, bound, x, y, x.size(), x_i, is_fixed);
  }

  auto size() const { return y.size(); }

  void check_is_fixed()
  {
    is_fixed = slide::check_is_fixed(x);
  }

  template <typename Tpath>
  void setCurve(Tpath &&path)
  {
    loadCSV_2col(path, x, y);
    check_is_fixed();
  }
};

using XYdata_ff = XYdata<FixedData<double>, FixedData<double>>;
using XYdata_fv = XYdata<FixedData<double>, std::vector<double>>;
using XYdata_vv = XYdata<std::vector<double>, std::vector<double>>;
using XYdata_ss = XYdata<std::span<double>, std::span<double>>;

template <typename Tpath>
void loadCSV_2col(Tpath &&name, slide::XYdata_vv &data, int n = 0)
{
  //!< slide::XYdata_vv overload.
  loadCSV_2col(name, data.x, data.y, n);
}

} // namespace slide
