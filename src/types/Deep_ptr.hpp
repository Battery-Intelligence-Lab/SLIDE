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
  std::unique_ptr<T> ptr;

public:
  Deep_ptr() = default;               //!< Default constructor
  explicit Deep_ptr(T *p) : ptr(p) {} //!< Constructor taking a raw pointer

  /// Copy constructor
  Deep_ptr(const Deep_ptr &other)
  {
    if (other.ptr)
      ptr = std::unique_ptr<T>(other.ptr->copy());
  }

  /// Move constructor
  Deep_ptr(Deep_ptr &&other) noexcept : ptr(std::move(other.ptr)) {}

  /// Copy assignment operator
  Deep_ptr &operator=(const Deep_ptr &other)
  {
    if (this != &other)
      if (other.ptr)
        ptr = std::unique_ptr<T>(other.ptr->copy());
      else
        ptr.reset();

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

  // Check if the Deep_ptr is not empty
  explicit operator bool() const { return static_cast<bool>(ptr); }
  void reset(T *p = nullptr) { ptr.reset(p); } //!< Reset the Deep_ptr

  // Swap the Deep_ptr with another Deep_ptr
  void swap(Deep_ptr &other) { ptr.swap(other.ptr); }
};

} // namespace slide