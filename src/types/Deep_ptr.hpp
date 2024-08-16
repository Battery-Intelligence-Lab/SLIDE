/**
 * @file Deep_ptr.hpp
 * @brief A pointer class, that deep copies itself when it is copied.
 * @author Volkan Kumtepeli
 * @date 05 Nov 2022
 * @details A pointer class, that deep copies itself when it is copied.
 * Otherwise, adding a unique_ptr member inside a class deletes the
 * default copy constructor.
 *
 * Code is adapted from: https://www.reddit.com/r/cpp_questions/comments/avr9op/specialize_a_copy_constructor_for_unique_ptr/
 * Accessed 2022-11-05
 */

#pragma once

#include <memory>
#include <type_traits> /// For std::is_base_of

namespace slide {

template <typename T>
class Deep_ptr : public std::unique_ptr<T>
{
public:
  Deep_ptr() = default;                                                           /// Default constructor
  explicit Deep_ptr(T *p) noexcept : std::unique_ptr<T>{ p } {}                   /// Constructor taking a raw pointer
  Deep_ptr(Deep_ptr &&other) noexcept : std::unique_ptr<T>{ std::move(other) } {} /// Move constructor

  /// Copy constructor
  Deep_ptr(const Deep_ptr &other) : std::unique_ptr<T>{ bool(other) ? other->copy() : nullptr } {}

  /// Templated copy constructor for implicit conversion from derived class.
  template <typename U, typename = std::enable_if_t<std::is_base_of<T, U>::value>>
  Deep_ptr(const Deep_ptr<U> &other) : std::unique_ptr<T>{ bool(other) ? other->copy() : nullptr } {}


  Deep_ptr(std::unique_ptr<T> &&up) : std::unique_ptr<T>{ std::move(up) } {}
  Deep_ptr &operator=(const Deep_ptr &rhs)
  {
    std::unique_ptr<T>::reset(Deep_ptr{ rhs }.std::unique_ptr<T>::release());
    return *this;
  }

  /// Copy assignment operator
  Deep_ptr &operator=(Deep_ptr &&rhs) noexcept
  {
    std::unique_ptr<T>::reset(rhs.release());
    return *this;
  }

  ~Deep_ptr() = default;
};

/// @brief  Makes storage units.
/// @tparam T
/// @tparam ...Args
/// @param ...args
/// @return
template <typename T, typename... Args>
Deep_ptr<T> make(Args &&...args) { return Deep_ptr<T>(new T(std::forward<Args>(args)...)); } // #TODO implicit conversion loses integer precision: 'const unsigned long long' to 'int' [-Wshorten-64-to-32]

} // namespace slide