#pragma once

// This file is for default parameter values.

#include "cell_param.hpp"

namespace settings
{
    namespace def_param
    {
        // fitting parameters of the models
        constexpr SEIparam SEIparam_Kokam{
            0.075e-14,  // sei1k
            130e3,      // sei1k_T
            2.75e-11,   // sei2k
            130e3,      // sei2k_T
            1.125e-14,  // sei2D
            20e3,       // sei2D_T
            1.1458e-15, // sei3k
            65e3,       // sei3k_T
            0.25e-15,   // sei3D
            200e3,      // sei3D_T
            7.5e-7      // sei_porosity
        };

        constexpr SEIparam SEIparam_LGCChemNMC{
            0.075e-14, // sei1k
            130e3,     // sei1k_T
            2.75e-11,  // sei2k
            130e3,     // sei2k_T
            2.5e-15,   // sei2D
            200e3,     // sei2D_T
            1e-11,     // sei3k
            0.0,       // sei3k_T
            1.05e-16,  // sei3D
            20e3,      // sei3D_T
            7.5e-7     // sei_porosity
        };

        constexpr SEIparam SEIparam_User{
            0.075e-14,  // sei1k
            130e3,      // sei1k_T
            2.75e-11,   // sei2k
            130e3,      // sei2k_T
            1.125e-14,  // sei2D
            20e3,       // sei2D_T
            1.1458e-15, // sei3k
            65e3,       // sei3k_T
            0.25e-15,   // sei3D
            200e3,      // sei3D_T
            7.5e-7      // sei_porosity
        };

        constexpr LAMparam LAMparam_Kokam{
            3.4985e-9, // lam1p
            5.88e-13,  // lam1n
            -1.675e-5, // lam2ap
            0.0,       // lam2bp
            -1.675e-5, // lam2an
            0.0,       // am2bn
            54611.0,   // lam2t
            12.5e-6,   // lam3k
            27305.0,   // lam3k_T
            8.3333e-9, // lam4p
            8.3333e-9  // lam4n
        };

        constexpr LAMparam LAMparam_LGCChemNMC{
            2.6031e-9,   // lam1p
            1.0417e-12,  // lam1n
            3.015e-11,   // lam2ap
            -1.72125e-6, // lam2bp
            3.015e-11,   // lam2an
            -1.72125e-6, // lam2bn
            54611.0,     // lam2t
            1.21e-6,     // lam3k
            27305.0,     // lam3k_T
            7.5e-9,      // lam4p
            7.5e-9       // lam4n
        };

        constexpr LAMparam LAMparam_User{
            3.4985e-9, // lam1p
            5.88e-13,  // lam1n
            -1.675e-5, // lam2ap
            0.0,       // lam2bp
            -1.675e-5, // lam2an
            0.0,       // am2bn
            54611.0,   // lam2t
            12.5e-6,   // lam3k
            27305.0,   // lam3k_T
            8.3333e-9, // lam4p
            8.3333e-9  // lam4n
        };

        constexpr StressParam StressParam_Kokam{
            2.1e-6,  // omegap, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
            3.17e-6, // omegan, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
            138.73,  // Ep, from Salco de Vasconcelos, Xu, Li, Zhao, Grid indentation analysis of mechanical properties of composite electrodes in Li-ion batteries, Extreme Mechanics Letters 9 (3), 2016
            10,      // En, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
            0.3,     // nup, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
            0.3,     // nun, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
            false,   // s_dai
            false,   // s_lares
            false,   // s_dai_update
            false,   // s_lares_update
            0,       // s_dai_p
            0,       // s_dai_n
            0,       // s_lares_n
            0,       // s_dai_p_prev
            0,       // s_dai_n_prev
            0        // s_lares_n_prev
        };

        constexpr StressParam StressParam_LGCChemNMC{
            2.1e-6,  // omegap, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
            3.17e-6, // omegan, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
            138.73,  // Ep, from Salco de Vasconcelos, Xu, Li, Zhao, Grid indentation analysis of mechanical properties of composite electrodes in Li-ion batteries, Extreme Mechanics Letters 9 (3), 2016
            10,      // En, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
            0.3,     // nup, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
            0.3,     // nun, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
            false,   // s_dai
            false,   // s_lares
            false,   // s_dai_update
            false,   // s_lares_update
            0,       // s_dai_p
            0,       // s_dai_n
            0,       // s_lares_n
            0,       // s_dai_p_prev
            0,       // s_dai_n_prev
            0        // s_lares_n_prev
        };

        constexpr StressParam StressParam_User{
            2.1e-6,  // omegap, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
            3.17e-6, // omegan, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
            138.73,  // Ep, from Salco de Vasconcelos, Xu, Li, Zhao, Grid indentation analysis of mechanical properties of composite electrodes in Li-ion batteries, Extreme Mechanics Letters 9 (3), 2016
            10,      // En, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
            0.3,     // nup, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
            0.3,     // nun, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
            false,   // s_dai
            false,   // s_lares
            false,   // s_dai_update
            false,   // s_lares_update
            0,       // s_dai_p
            0,       // s_dai_n
            0,       // s_lares_n
            0,       // s_dai_p_prev
            0,       // s_dai_n_prev
            0        // s_lares_n_prev
        };

        constexpr StressParam StressParam_Fit{
            2.1e-6,  // omegap, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
            3.17e-6, // omegan, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
            138.73,  // Ep, from Salco de Vasconcelos, Xu, Li, Zhao, Grid indentation analysis of mechanical properties of composite electrodes in Li-ion batteries, Extreme Mechanics Letters 9 (3), 2016
            10,      // En, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
            0.3,     // nup, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
            0.3,     // nun, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
            false,   // s_dai
            false,   // s_lares
            false,   // s_dai_update
            false,   // s_lares_update
            0,       // s_dai_p
            0,       // s_dai_n
            0,       // s_lares_n
            0,       // s_dai_p_prev
            0,       // s_dai_n_prev
            0        // s_lares_n_prev
        };

        // constexpr slide::State iniState_Kokam(){

        // };

    }

}