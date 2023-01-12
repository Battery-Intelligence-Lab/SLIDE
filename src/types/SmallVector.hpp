

#pragma once

#include <array>

namespace slide {

template <typename T, T Nmax>
class SmallVector
{
  std::array<T, Nmax> data{};

public:
  T N{ 0 }; //!< number of degradation models to use (length of SEI_ID) #TODO if default should be 1 and array initialised with 0.

  SmallVector() = default;
  SmallVector(T N) : N(N)
  {
    if (N > Nmax) {
      std::cerr << "Array is very small! Please increase its size from settings!\n";
      throw 1234;
    }

    for (T i = 0; i < N; i++)
      data[i] = 0;
  }

  inline void push_back(const T &&elem)
  {
    if (N > Nmax) {
      std::cerr << "Array is very small! Please increase its size from settings!\n";
      throw 1234;
    }

    data[N] = elem;
    N++;
  }

  inline auto &operator[](T idx)
  {
    if (idx < N)
      return data[idx];
    else {
      std::cerr << "Array is very small! Please increase its size from settings!\n";
      throw 1234;
    }
  }

  inline auto empty() { return N == 0; }

  [[nodiscard]] inline const auto &operator[](T idx) const noexcept { return data[idx]; }

  [[nodiscard]] constexpr auto begin() noexcept { return data.begin(); }
  [[nodiscard]] constexpr auto cbegin() const noexcept { return data.cbegin(); }

  [[nodiscard]] constexpr auto end() { return data.begin() + N; }
  [[nodiscard]] constexpr auto cend() const noexcept { return data.cbegin() + N; }
};

} // namespace slide