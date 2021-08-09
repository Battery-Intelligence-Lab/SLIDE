/*
 * Util.cpp
 *
 * All utility functions. So the code is less verbose. 
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#include <vector>
#include <cmath>

#include "util.hpp"
#include "slide_aux.hpp"

namespace slide
{

    slide::fixed_data<double> range_fix(double x_min, double x_max, double x_step)
    {
        int Nstep = static_cast<int>((x_max - x_min) / x_step) + 1;
        return slide::fixed_data<double>(x_min, x_step, Nstep);
    }

    slide::fixed_data<double> linstep_fix(double x_min, double x_step, int Nstep)
    {
        return slide::fixed_data<double>(x_min, x_step, Nstep);
    }
    slide::fixed_data<double> logstep_fix(double x_min, double x_step, int Nstep)
    {
        auto fun = [](double x_min, double x_step, int i)
        { return x_min * std::pow(x_step, i); };

        return slide::fixed_data<double>(x_min, x_step, Nstep, fun);
    }

    slide::fixed_data<double> linspace_fix(double x1, double x2, int N)
    {
        if (N < 1)
            return fixed_data(x2, 0.0, 0);
        else if (N == 1)
            return fixed_data(x2, 0.0, 1);
        else
            return fixed_data(x1, (x2 - x1) / static_cast<double>(N - 1), N);
    }

    std::vector<double> linspace(double x1, double x2, int N)
    {
        N = std::max(0, N);
        std::vector<double> out(N);

        double dx{1};
        if (N < 1)
            return out;
        else if (N > 1)
            dx = (x2 - x1) / static_cast<double>(N - 1);

        for (int n{0}; n < N; n++)
        {
            out[N - 1 - n] = x2 - n * dx;
        }

        return out;
    }
}
// namespace slide::util