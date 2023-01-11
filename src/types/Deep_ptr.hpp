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
  explicit Deep_ptr(T *t) : ptr(t) {}

  Deep_ptr(Deep_ptr const &other) : ptr(new T(*other.ptr)) {}
  Deep_ptr(Deep_ptr &&other) : ptr(std::move(other.ptr)) {}
  Deep_ptr &operator=(Deep_ptr const &other) { ptr.reset(new T(*other.ptr)); }
  Deep_ptr &operator=(Deep_ptr &&other) { ptr = std::move(other.ptr); }

  T *release() noexcept { return ptr.release(); }
  void reset(T *t = nullptr) noexcept { ptr.reset(t); }
  void swap(Deep_ptr &other) noexcept { this->ptr.swap(other.ptr); }

  explicit operator bool() const noexcept { return bool(ptr); }

  T *get() const noexcept { return ptr.get(); }
  T &operator*() const { return *ptr.get(); }
  T *operator->() const noexcept { return ptr.get(); }
};
} // namespace slide