/**
 * @file slide_aux.hpp
 * @brief Auxillary classes and functions. So the code is less verbose.
 * @author Volkan Kumtepeli
 * @author Jorn Reniers
 * @date 8 Aug 2022
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
#include <algorithm>

namespace slide {

template <int N, typename T> //!< Takes vector or array:
double norm_sq(const T &x)
{
  if constexpr (N < -1)
    throw "N cannot be less than -1";
  else if (N == -1) //!< Infinite norm
    return *std::max_element(x.begin(), x.end(), [](auto a, auto b) {
      return std::abs(a) < std::abs(b);
    });
  else if (N == 0) //!< 0 norm
    return std::count_if(x.begin(), x.end(), [](auto a) {
      return std::abs(a) > 1e-15;
    });
  else {
    double sum = 0;
    for (const auto &x_i : x)
      sum += std::pow(x_i, N);

    return sum;
  }
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
