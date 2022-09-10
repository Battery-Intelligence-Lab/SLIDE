/*
 * matrix.hpp
 *
 * matrix data type
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include "DynamicMatrix.hpp"

#include <cstdlib>
#include <array>
#include <cmath>

namespace slide {
template <typename T, size_t ROW, size_t COL>
using Matrix = std::array<std::array<T, COL>, ROW>; //!< See source: http://cpptruths.blogspot.com/2011/10/multi-dimensional-arrays-in-c11.html

template <size_t N, size_t M = N>
auto eye(double k = 1.0)
{
  auto A = Matrix<double, M, N>{};
  size_t m = std::min(M, N);
  for (size_t i = 0; i < m; i++)
    A[i][i] = k;

  return A;
}

template <size_t N, size_t M = N>
auto zeros()
{
  auto A = Matrix<double, M, N>{};
  return A;
}

template <int N, int M = N>
auto ones(double k = 1.0)
{
  auto A = Matrix<double, M, N>{};

  for (int i = 0; i < M; i++)
    for (int j = 0; j < N; j++)
      A[i][j] = k;

  return A;
}

template <int N>
void cholUpdate(slide::Matrix<double, N, N> &L, std::array<double, N> x, bool isDowndate = false)
{
  //!< In-place rank-one update for Cholesky
  for (int k = 0; k < N; k++) {
    const double Lk = L[k][k];
    const double s = x[k] / Lk;
    double r;
    if (isDowndate)
      r = std::sqrt(Lk * Lk - x[k] * x[k]);
    else
      r = std::hypot(Lk, x[k]);

    const double c = r / Lk;
    L[k][k] = r;

    for (int j = k + 1; (k != (N - 1) && (j < N)); j++) {
      if (isDowndate)
        L[j][k] -= s * x[j] / c;
      else
        L[j][k] += s * x[j] / c;
      x[j] = c * x[j] - s * L[j][k];
    }
  }
}

} // namespace slide