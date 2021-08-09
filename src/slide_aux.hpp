/*
 * slide_aux.hpp
 *
 * Auxillary classes and functions. So the code is less verbose. 
 *
 * A cycler implements check-up procedures and degradation procedures.
 * The data from the check-up procedures is written in csv files in the same subfolder as where the cycling data of a cell is written (see BasicCycler.cpp).
 * There is one file per 'type' of check-up (capacity measurement, OCV measurement, CCCV cycles and a pulse discharge).
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include <stdexcept>
#include <vector>
#include <array>
#include <functional>

//#include "interpolation.h"

namespace slide
{

    template <typename Tx, typename Ty>
    class XYdata
    {
    public:
        Tx x;
        Ty y;

        XYdata() = default;
        XYdata(size_t N) : x(N), y(N) {}
        XYdata(Tx x, Ty y) : x(x), y(y) {}

        void reserve(int n)
        {
            x.reserve(n), y.reserve(n);
        }
        void clear()
        {
            x.clear(), y.clear();
        }
        //double interp(double x_i) { return ::linInt(false, true, x, y, x.size(), x_i, false); }
        auto size() const { return y.size(); }
    };

    template <typename T, bool extrapolation = true>
    class fixed_data;

    template <typename T, bool extrapolation = true>
    class FixedDataIter
    { // For more information, see: https://internalpointers.com/post/writing-custom-iterators-modern-cpp
    public:
        using iterator_category = std::input_iterator_tag;
        using difference_type = int;
        using value_type = T;
        using pointer = value_type *;
        using reference = value_type &;

        FixedDataIter(fixed_data<T, extrapolation> *const f_data_ptr, int n) : f_data_ptr{f_data_ptr}, n{n} {}

        FixedDataIter &operator++() // Prefix increment
        {
            n++;
            return *this;
        }

        FixedDataIter operator++(int) // Postfix increment
        {
            FixedDataIter temp = *this;
            ++(*this);
            return temp;
        }

        FixedDataIter &operator--() // Prefix decrement
        {
            n--;
            return *this;
        }

        FixedDataIter operator--(int) // Postfix decrement
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
            return FixedDataIter{a.f_data_ptr, a.n + b};
        }

        friend difference_type operator-(const FixedDataIter &a, const FixedDataIter &b)
        {
            return a.n - b.n;
        }

    private:
        fixed_data<T, extrapolation> *f_data_ptr;
        int n;
    };

    template <typename T, bool extrapolation = true>
    int distance(const FixedDataIter<T, extrapolation> &a, const FixedDataIter<T, extrapolation> &b)
    {
        return a - b;
    }

    template <typename T, bool extrapolation>
    class fixed_data
    {
        using size_type = int;
        T x_min{};
        T dx{};
        size_type n{};
        std::function<T(T, T, int)> F = [](T x_min, T dx, size_type i)
        { return x_min + i * dx; };

    public:
        fixed_data() = default;
        fixed_data(T x_min, T dx, int n) : x_min(x_min), dx(dx), n(n) {}
        fixed_data(T x_min, T dx, int n, std::function<T(T, T, int)> F) : x_min(x_min), dx(dx), n(n), F(F) {}

        T operator[](int i) const
        {

            if constexpr (!extrapolation) // Extrapolation will be deprecated.
            {
                if (i >= n || i < 0)
                {
                    std::cout << i << ' ' << n << '\n';
                    throw std::out_of_range("fixed data is out of range");
                }
            }

            auto x = F(x_min, dx, i);
            return x;
        }

        T operator()(int i) const // Can extrapolate (go out of bounds).
        {
            return F(x_min, dx, i);
        }

        FixedDataIter<T, extrapolation> begin() { return FixedDataIter(this, 0); }
        FixedDataIter<T, extrapolation> end() { return FixedDataIter(this, n); }

        const FixedDataIter<T, extrapolation> cbegin() const { return FixedDataIter(this, 0); }
        const FixedDataIter<T, extrapolation> cend() const { return FixedDataIter(this, n); }

        T back() const { return operator[](n - 1); }
        T front() const { return operator[](0); }
        T dstep() const { return dx; }

        T prev(T x_current)
        {
            // Gets previous point compared to the current point x_current
            auto temp_data = *this;
            temp_data.x_min = x_current;
            return temp_data(-1);
        }

        T next(T x_current)
        {
            // Gets next point compared to the current point x_current
            auto temp_data = *this;
            temp_data.x_min = x_current;
            return temp_data(+1);
        }

        void reserve(int n) {}
        void clear()
        {
            x_min = 0;
            dx = 0;
            n = 0;
        }
        auto size() const { return n; }
    };

    using fixed_XYdata = XYdata<fixed_data<double>, std::vector<double>>;
    using vec_XYdata = XYdata<std::vector<double>, std::vector<double>>;

    template <typename T, size_t ROW, size_t COL>
    using Matrix = std::array<std::array<T, COL>, ROW>; // See source: http://cpptruths.blogspot.com/2011/10/multi-dimensional-arrays-in-c11.html

    class card_prod
    {
    };

}

// Array operations to make summing/subtracting arrays easier.

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
