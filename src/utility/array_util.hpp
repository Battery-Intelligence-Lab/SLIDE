/*
 * array_util.cpp
 *
 * Array operations to make summing/subtracting arrays easier.
 *
 * Created on: 04 Jan 2023
 *  Author(s): Jorn Reniers, Volkan Kumtepeli
 */


#pragma once

#include <array>
#include <cmath>

namespace slide {

template <typename T, size_t N>
auto arrSum(const std::array<T, N> &a1, const std::array<T, N> &a2, double b1, double b2)
{
  std::array<T, N> c;

  for (size_t i = 0; i < N; i++)
    c[i] = b1 * a1[i] + b2 * a2[i];

  return c;
}

template <typename T, size_t N>
auto operator+=(std::array<T, N> &c, const std::array<T, N> &b)
{
  for (size_t i = 0; i < N; i++)
    c[i] += b[i];

  return c;
}

template <typename T, size_t N>
auto operator+(const std::array<T, N> &a, const std::array<T, N> &b)
{
  auto c = a;
  c += b;
  return c;
}

template <typename T, size_t N>
auto operator-=(std::array<T, N> &c, const std::array<T, N> &b)
{
  for (size_t i = 0; i < N; i++)
    c[i] -= b[i];

  return c;
}

template <typename T, size_t N>
auto operator-(const std::array<T, N> &a, const std::array<T, N> &b)
{
  auto c = a;
  c -= b;
  return c;
}
} // namespace slide