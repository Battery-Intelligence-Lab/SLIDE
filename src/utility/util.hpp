/*
 * util.hpp
 *
 * Some utility functions. So the code is less verbose.
 *
 * A cycler implements check-up procedures and degradation procedures.
 * The data from the check-up procedures is written in csv files in the same subfolder as where the cycling data of a cell is written (see BasicCycler.cpp).
 * There is one file per 'type' of check-up (capacity measurement, OCV measurement, CCCV cycles and a pulse discharge).
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

//!< Include other util files.
#pragma once

#include "../settings/settings.hpp"
#include "slide_aux.hpp"
#include "util_debug.hpp"
#include "units.hpp"


#include <string>
#include <iostream>
#include <array>
#include <vector>
#include <cmath>

namespace slide {

constexpr inline auto abs_sqrt(auto x) { return std::sqrt(std::abs(x)); }

constexpr inline auto sqr(auto x) { return x * x; }
constexpr inline auto cube(auto x) { return x * x * x; }


template <typename T>
void output_printer(const std::vector<T> &vec, const auto &save_path)
{
  std::ofstream out_file(save_path, std::ios_base::out);

  for (const auto &state : vec) {
    for (size_t i{ 0 }; i < state.size(); i++) {
      if (i != 0) out_file << ',';
      out_file << state[i];
    }
  }

  out_file.close();
}
} // namespace slide


namespace slide::util {

} // namespace slide::util

namespace slide {
std::vector<double> linstep(double x_min, double x_step, int Nstep); // #TODO not defined.
std::vector<double> logstep(double x_min, double x_step, int Nstep); // #TODO not defined.

//!< FixedData<int> range(int stop); #TODO -> FixedData is not good since it has function
//!< FixedData<int> range(int start, int stop, int step = 1);

FixedData<double> linspace_fix(double x1, double x2, int N);

template <size_t N>
constexpr std::array<double, N> linspace(double x1, double x2)
{
  std::array<double, N> out;

  double dx{ 1 };
  if (N < 1)
    return out;
  else if (N > 1)
    dx = (x2 - x1) / static_cast<double>(N - 1);

  for (size_t n{ 0 }; n < N; n++) {
    out[N - 1 - n] = x2 - n * dx;
  }

  return out;
}

} // namespace slide

namespace slide {

inline FixedData<double> range_fix(double x_min, double x_max, double x_step)
{
  int Nstep = static_cast<int>((x_max - x_min) / x_step) + 1;
  return FixedData<double>(x_min, x_step, Nstep);
}

inline FixedData<double> linstep_fix(double x_min, double x_step, int Nstep)
{
  return FixedData<double>(x_min, x_step, Nstep);
}

inline FixedData<double> logstep_fix(double x_min, double x_step, int Nstep)
{
  auto fun = [](double x_min, double x_step, int i) { return x_min * std::pow(x_step, i); };

  return FixedData<double>(x_min, x_step, Nstep, fun);
}

inline FixedData<double> linspace_fix(double x1, double x2, int N)
{
  if (N < 1)
    return FixedData(x2, 0.0, 0);
  else if (N == 1)
    return FixedData(x2, 0.0, 1);
  else
    return FixedData(x1, (x2 - x1) / static_cast<double>(N - 1), N);
}

inline std::vector<double> linspace(double x1, double x2, int N)
{
  N = std::max(0, N);
  std::vector<double> out(N);

  double dx{ 1 };
  if (N < 1)
    return out;
  else if (N > 1)
    dx = (x2 - x1) / static_cast<double>(N - 1);

  for (int n{ 0 }; n < N; n++) {
    out[N - 1 - n] = x2 - n * dx;
  }

  return out;
}
} // namespace slide