/*
 * DynamicMatrix.hpp
 *
 *  Created on: 14 Jul 2022
 *   Author(s): Volkan Kumtepeli
 */

#pragma once

#include <iostream>
#include <vector>

namespace slide {
template <typename Tdata>
class DynamicMatrix
{
  using VarType = int;
  using data_type = Tdata;
  VarType m{}, n{};

public:
  std::vector<Tdata> data;
  DynamicMatrix() = default;
  DynamicMatrix(VarType m_) : m(m_), n(m_), data(m_ * m_) {} //!< Sequare matrix
  DynamicMatrix(VarType m_, VarType n_) : m(m_), n(n_), data(m_ * n_) {}
  DynamicMatrix(VarType m_, VarType n_, Tdata x) : m(m_), n(n_), data(m_ * n_, x) {}

  inline void resize(VarType m_, VarType n_, Tdata x = 0)
  {
    m = m_;
    n = n_;
    data.resize(m_ * n_, x);
  }

  inline void reshape(VarType m_, VarType n_)
  {
    m = m_;
    n = n_;

    //!< #TODO add only one to be changed.
  }

  inline auto rows() const { return m; }
  inline auto cols() const { return n; }
  inline auto size() const { return m * n; }

  inline auto &operator()(VarType i, VarType j)
  {
    return data[i + j * m];
  }

  inline auto operator()(VarType i, VarType j) const
  {
    return data[i + j * m];
  }

  void print(std::ostream &os = std::cout) const
  {
    for (VarType i = 0; i < m; i++) {
      for (VarType j = 0; j < n; j++)
        os << this->operator()(i, j) << ", "; //!< \t\t

      os << '\n';
    }
  }

  template <typename TTdata>
  friend std::ostream &operator<<(std::ostream &os, const DynamicMatrix<TTdata> &M);
};

template <typename Tdata>
std::ostream &operator<<(std::ostream &os, const DynamicMatrix<Tdata> &M) //!< #TODO in future this should be const.
{
  M.print(os);
  return os;
}
} // namespace slide
