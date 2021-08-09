/*
 * util.hpp
 *
 * All utility functions. So the code is less verbose. 
 *
 * A cycler implements check-up procedures and degradation procedures.
 * The data from the check-up procedures is written in csv files in the same subfolder as where the cycling data of a cell is written (see BasicCycler.cpp).
 * There is one file per 'type' of check-up (capacity measurement, OCV measurement, CCCV cycles and a pulse discharge).
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

// Include other util files.
#pragma once

#include <string>
#include <iostream>
#include <array>
#include <thread>

#include "util_error.hpp"
#include "constants.hpp"
#include "slide_aux.hpp"

namespace slide
{
    template <typename Tfun>
    void run(Tfun task_indv, int i_end, unsigned int numMaxParallelWorkers = settings::numMaxParallelWorkers)
    {

        auto task_par = [&](int i_begin, int i_end, int Nth)
        {
            while (i_begin < i_end)
            {
                task_indv(i_begin);
                i_begin += Nth;
            }
        };

        if constexpr (settings::isParallel)
        {
            if (numMaxParallelWorkers <= 1)
            {
                task_par(0, i_end, 1);
            }
            else
            {
                const unsigned int N_th_max = std::min(numMaxParallelWorkers, std::thread::hardware_concurrency());
                std::vector<std::thread> threads;
                threads.reserve(N_th_max);

                for (unsigned int i_begin = 0; i_begin < N_th_max; i_begin++) // indices for the threads
                {
                    // Multi threaded simul:

                    threads.emplace_back(task_par, i_begin, i_end, N_th_max);
                }

                for (auto &th : threads)
                {
                    if (th.joinable())
                        th.join();
                }
            }
        }
        else
        {
            task_par(0, i_end, 1);
        }
    }
}

namespace slide::util
{
    //template <int pLevel>
    //inline void print(const std::string& message) // Print messages.
    // {
    // if constexpr (settings::verbose >= pLevel)
    //  {
    //       std::cout << message << "\n";
    //    }
    // }

    struct Counter // Counts how many times it is called.
    {

        // To use counter create a static Counter object in the function you would like to measure
        // slide::util::Counter myCounterName(printEveryXiter);  Then just call it
        // myCounterName();
        std::atomic_ullong count{0};
        long long printEveryXiter{500000};

        Counter() = default;
        Counter(unsigned long long printEveryXiter) : printEveryXiter(printEveryXiter){};

        void operator()()
        {
            count++;

            if (count % printEveryXiter == 0)
                std::cout << "Counter: " << count << '\n';
        }
    };

} // namespace slide::util

namespace slide
{

    std::vector<double> linstep(double x_min, double x_step, int Nstep);
    std::vector<double> logstep(double x_min, double x_step, int Nstep);

    slide::fixed_data<double> linstep_fix(double x_min, double x_step, int Nstep);
    slide::fixed_data<double> logstep_fix(double x_min, double x_step, int Nstep);

    slide::fixed_data<double> range_fix(double x_min, double x_max, double dx);

    std::vector<double> linspace(double x1, double x2, int N);
    slide::fixed_data<double> linspace_fix(double x1, double x2, int N);

    template <size_t N>
    std::array<double, N> linspace(double x1, double x2)
    {
        std::array<double, N> out;

        double dx{1};
        if (N < 1)
            return out;
        else if (N > 1)
            dx = (x2 - x1) / static_cast<double>(N - 1);

        for (size_t n{0}; n < N; n++)
        {
            out[N - 1 - n] = x2 - n * dx;
        }

        return out;
    }
}