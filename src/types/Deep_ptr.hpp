/*
 * Deep_ptr.hpp
 *
 * A pointer class, that deep copies itself when it is copied.
 * Otherwise, adding a unique_ptr member inside a class deletes the
 * default copy constructor.
 *
 * Code is adapted from: https://www.reddit.com/r/cpp_questions/comments/avr9op/specialize_a_copy_constructor_for_unique_ptr/
 * Accessed 2022-11-05
 *
 * Created on: 05 Nov 2022
 * Author(s): Volkan Kumtepeli
 */

#pragma once

#include <memory>

namespace slide {

template <typename T>
class Deep_ptr
{
  std::unique_ptr<T> ptr{};

public:
  Deep_ptr() = default;               //!< Default constructor
  explicit Deep_ptr(T *p) : ptr(p) {} //!< Constructor taking a raw pointer

  Deep_ptr(std::unique_ptr<T> &&uptr) : ptr(std::move(uptr)) {} // Constructor taking a unique_ptr

  /// Copy constructor
  Deep_ptr(const Deep_ptr &other) : ptr(other ? other->copy() : nullptr) {}

  // Templated copy constructor for implicit conversion from derived class.
  template <typename U, typename = std::enable_if_t<std::is_base_of<T, U>::value>>
  Deep_ptr(const Deep_ptr<U> &other) : ptr(other ? other->copy() : nullptr) {}

  // Templated constructor for implicit conversion from derived class for unique_ptr.
  template <typename U, typename = std::enable_if_t<std::is_base_of<T, U>::value>>
  Deep_ptr(const std::unique_ptr<U> &other) : ptr(other ? other->copy() : nullptr) {}

  /// Move constructor
  Deep_ptr(Deep_ptr &&other) noexcept : ptr(std::move(other.ptr)) {}

  // Templated move constructor for implicit conversion
  template <typename U, typename = std::enable_if_t<std::is_base_of<T, U>::value>>
  Deep_ptr(Deep_ptr<U> &&other) noexcept : ptr(std::move(other.ptr)) {}

  // // Templated constructor for implicit conversion from raw pointers
  // template <typename U, typename = std::enable_if_t<std::is_base_of<T, U>::value>>
  // explicit Deep_ptr(U *raw_ptr) : ptr(raw_ptr) {}

  // Declare the templated Deep_ptr class as a friend
  template <typename U>
  friend class Deep_ptr;

  /// Copy assignment operator
  Deep_ptr &operator=(const Deep_ptr &other)
  {
    if (this != &other) ptr = other ? other->copy() : nullptr;

    return *this;
  }

  /// Move assignment operator
  Deep_ptr &operator=(Deep_ptr &&other) noexcept
  {
    if (this != &other) ptr = std::move(other.ptr);

    return *this;
  }

  T *operator->() const { return ptr.operator->(); } //!< Accessor operator
  T &operator*() const { return *ptr; }              //!< Dereference operator
  T *get() const { return ptr.get(); }               //!< Get the raw pointer

  /// Check if the Deep_ptr is not empty
  explicit operator bool() const { return static_cast<bool>(ptr); }
  void reset(T *p = nullptr) { ptr.reset(p); } //!< Reset the Deep_ptr

  // Swap the Deep_ptr with another Deep_ptr
  void swap(Deep_ptr &other) { ptr.swap(other.ptr); }

  // Comparison operators
  [[nodiscard]] friend bool operator==(const Deep_ptr<T> &lhs, std::nullptr_t) noexcept
  {
    return lhs.ptr == nullptr;
  }

  [[nodiscard]] friend bool operator==(std::nullptr_t, const Deep_ptr<T> &rhs) noexcept
  {
    return rhs == nullptr;
  }

  [[nodiscard]] friend bool operator!=(const Deep_ptr<T> &lhs, std::nullptr_t) noexcept
  {
    return !(lhs == nullptr);
  }

  [[nodiscard]] friend bool operator!=(std::nullptr_t, const Deep_ptr<T> &rhs) noexcept
  {
    return rhs != nullptr;
  }
};

// Free functions:
template <typename T, typename... Args>
Deep_ptr<T> make_Deep_ptr(Args &&...args)
{
  return Deep_ptr<T>(new T(std::forward<Args>(args)...));
}

/// @brief  Makes storage units.
/// @tparam T
/// @tparam ...Args
/// @param ...args
/// @return
template <typename T, typename... Args>
Deep_ptr<T> make(Args &&...args) { return Deep_ptr<T>(new T(std::forward<Args>(args)...)); }

} // namespace slide