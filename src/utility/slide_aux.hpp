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
#include <cmath>
#include <filesystem>

#include "interpolation.hpp"
#include "../types/matrix.hpp"
#include "../types/XYdata.hpp"
#include "../types/FixedData.hpp"

namespace slide
{
    template <typename StoredClass, int N>
    struct DataStore
    {
        // Will be filled in future.
    };

    template <int N, typename T> // Takes vector or array:
    double norm_sq(const T &x)
    {
        if constexpr (N < -1)
            throw 12345;
        else if (N == -1)
            throw 12344; // Infinite norm
        else if (N == 0)
            throw 333456; // 0 norm
        else
        {
            double sum = 0;
            for (const auto &x_i : x)
                sum += std::pow(x_i, N);

            return sum;
        }
        throw 345789; // This should not happen.
    }

    template <int N, typename T> // Takes vector or array:
    double norm(const T &x)
    {
        const double norm_square = norm_sq(x);

        if (N < 2)
            return norm_square;
        else
        {
            const double N_double = N;
            return std::pow(norm_square, 1.0 / N_double);
        }
    }

    class card_prod
    {
    };

    struct GammaDensityFunctor
    {
        double a, inv_b, scale;
        GammaDensityFunctor() = delete; // Even it is delete, it constructs with default values.
        GammaDensityFunctor(double a, double b) : a(a), inv_b(1 / b), scale(1 / std::tgamma(a) / std::pow(b, a)) {}
        double operator()(double x) { return scale * std::exp(-inv_b * x) * std::pow(x, a - 1); }
    };

}

// Array operations to make summing/subtracting arrays easier.

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
