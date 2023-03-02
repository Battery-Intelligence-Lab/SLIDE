/*
 * slide_aux.hpp
 *
 * Auxillary classes and functions. So the code is less verbose.
 *
 * A cycler implements check-up procedures and degradation procedures.
 * The data from the check-up procedures is written in csv files in the same subfolder as where the cycling data of a cell is written (see BasicCycler.cpp).
 * There is one file per 'type' of check-up (capacity measurement, OCV measurement, CCCV cycles and a pulse discharge).
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include "array_util.hpp"
#include "interpolation.hpp"
#include "../types/matrix.hpp"
#include "../types/XYdata.hpp"
#include "../types/FixedData.hpp"

#include <stdexcept>
#include <vector>
#include <array>
#include <functional>
#include <cmath>
#include <filesystem>

namespace slide {
template <typename StoredClass, int N>
struct DataStore
{
  //!< Will be filled in future.
};

template <int N, typename T> //!< Takes vector or array:
double norm_sq(const T &x)
{
  if constexpr (N < -1)
    throw 12345;
  else if (N == -1)
    throw 12344; //!< Infinite norm
  else if (N == 0)
    throw 333456; //!< 0 norm
  else {
    double sum = 0;
    for (const auto &x_i : x)
      sum += std::pow(x_i, N);

    return sum;
  }
  throw 345789; //!< This should not happen.
}

template <int N, typename T> //!< Takes vector or array:
double norm(const T &x)
{
  const double norm_square = norm_sq(x);

  if (N < 2)
    return norm_square;
  else {
    const double N_double = N;
    return std::pow(norm_square, 1.0 / N_double);
  }
}

class card_prod
{
};

struct GammaDensityFunctor
{
  double a, inv_b, scale;
  GammaDensityFunctor() = delete; //!< Even it is delete, it constructs with default values.
  GammaDensityFunctor(double a_, double b_) : a(a_), inv_b(1 / b_), scale(1 / std::tgamma(a_) / std::pow(b_, a_)) {}
  double operator()(double x) { return scale * std::exp(-inv_b * x) * std::pow(x, a - 1); }
};

} // namespace slide
