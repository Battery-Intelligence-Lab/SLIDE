/*
 * FixedData.hpp
 *
 *  Created on: 07 Feb 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <stdexcept>
#include <vector>
#include <array>
#include <functional>
#include <cmath>
#include <filesystem>

namespace slide {
template <typename T, bool extrapolation = true>
class FixedData;

template <typename T, bool extrapolation = true>
class FixedDataIter
{ //!< For more information, see: https://internalpointers.com/post/writing-custom-iterators-modern-cpp
public:
  using iterator_category = std::input_iterator_tag;
  using difference_type = int;
  using value_type = T;
  using pointer = value_type *;
  using reference = value_type &;

  FixedDataIter(FixedData<T, extrapolation> *const f_data_ptr_, int n) : f_data_ptr{ f_data_ptr_ }, n{ n } {}

  FixedDataIter &operator++() //!< Prefix increment
  {
    n++;
    return *this;
  }

  FixedDataIter operator++(int) //!< Postfix increment
  {
    FixedDataIter temp = *this;
    ++(*this);
    return temp;
  }

  FixedDataIter &operator--() //!< Prefix decrement
  {
    n--;
    return *this;
  }

  FixedDataIter operator--(int) //!< Postfix decrement
  {
    FixedDataIter temp = *this;
    --(*this);
    return temp;
  }

  value_type operator[](int i)
  {
    return f_data_ptr->operator[](i);
  }

  const value_type operator*() const { return f_data_ptr->operator[](n); }

  friend bool operator==(const FixedDataIter &a, const FixedDataIter &b) { return a.n == b.n; }
  friend bool operator!=(const FixedDataIter &a, const FixedDataIter &b) { return a.n != b.n; }

  friend FixedDataIter operator+(const FixedDataIter &a, const int b)
  {
    return FixedDataIter{ a.f_data_ptr, a.n + b };
  }

  friend difference_type operator-(const FixedDataIter &a, const FixedDataIter &b)
  {
    return a.n - b.n;
  }

private:
  FixedData<T, extrapolation> *f_data_ptr;
  int n;
};

template <typename T, bool extrapolation = true>
int distance(const FixedDataIter<T, extrapolation> &a, const FixedDataIter<T, extrapolation> &b)
{
  return a - b;
}

template <typename T, bool extrapolation>
class FixedData
{
  using size_type = int;
  T x_min{};
  T dx{};
  size_type n{};
  std::function<T(T, T, int)> F = [](T x_min, T dx, size_type i) { return x_min + i * dx; };

public:
  FixedData() = default;
  FixedData(T x_min, T dx, int n) : x_min(x_min), dx(dx), n(n) {}
  FixedData(T x_min, T dx, int n, std::function<T(T, T, int)> F) : x_min(x_min), dx(dx), n(n), F(F) {}

  T operator[](int i) const
  {

    if constexpr (!extrapolation) //!< Extrapolation will be deprecated.
    {
      if (i >= n || i < 0) {
        std::cout << i << ' ' << n << '\n';
        throw std::out_of_range("fixed data is out of range");
      }
    }

    auto x = F(x_min, dx, i);
    return x;
  }

  T operator()(int i) const //!< Can extrapolate (go out of bounds).
  {
    return F(x_min, dx, i);
  }

  [[nodiscard]] constexpr FixedDataIter<T, extrapolation> begin() noexcept { return FixedDataIter(this, 0); }
  [[nodiscard]] constexpr FixedDataIter<T, extrapolation> end() noexcept { return FixedDataIter(this, n); }

  [[nodiscard]] constexpr const FixedDataIter<T, extrapolation> cbegin() const noexcept { return FixedDataIter(this, 0); }
  [[nodiscard]] constexpr const FixedDataIter<T, extrapolation> cend() const noexcept { return FixedDataIter(this, n); }

  T back() const { return operator[](n - 1); }
  T front() const { return operator[](0); }
  T dstep() const noexcept { return dx; }

  T prev(T x_current) noexcept
  {
    //!< Gets previous point compared to the current point x_current
    auto temp_data = *this;
    temp_data.x_min = x_current;
    return temp_data(-1);
  }

  T next(T x_current) noexcept
  {
    //!< Gets next point compared to the current point x_current
    auto temp_data = *this;
    temp_data.x_min = x_current;
    return temp_data(+1);
  }

  constexpr void reserve(int n) noexcept {}
  constexpr void clear() noexcept
  {
    x_min = 0;
    dx = 0;
    n = 0;
  }
  constexpr auto size() const noexcept { return static_cast<size_t>(n); }
};
} // namespace slide