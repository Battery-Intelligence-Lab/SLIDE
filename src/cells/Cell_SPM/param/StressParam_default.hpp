/*
 * StressParam_default.hpp
 *
 *
 *
 *  Created on: 28 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 *
 */

#pragma once

#include "StressParam.hpp"

namespace slide::param::def {
constexpr StressParam StressParam_Kokam{
  2.1e-6,  //!< < omegap, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
  3.17e-6, //!< < omegan, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
  138.73,  //!< < Ep, from Salco de Vasconcelos, Xu, Li, Zhao, Grid indentation analysis of mechanical properties of composite electrodes in Li-ion batteries, Extreme Mechanics Letters 9 (3), 2016
  10,      //!< < En, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
  0.3,     //!< < nup, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
  0.3,     //!< < nun, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
  false,   //!< < s_dai
  false,   //!< < s_lares
  false,   //!< < s_dai_update
  false,   //!< < s_lares_update
  0,       //!< < s_dai_p
  0,       //!< < s_dai_n
  0,       //!< < s_lares_n
  0,       //!< < s_dai_p_prev
  0,       //!< < s_dai_n_prev
  0        //!< < s_lares_n_prev
};

constexpr StressParam StressParam_LGCChemNMC{
  2.1e-6,  //!< < omegap, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
  3.17e-6, //!< < omegan, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
  138.73,  //!< < Ep, from Salco de Vasconcelos, Xu, Li, Zhao, Grid indentation analysis of mechanical properties of composite electrodes in Li-ion batteries, Extreme Mechanics Letters 9 (3), 2016
  10,      //!< < En, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
  0.3,     //!< < nup, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
  0.3,     //!< < nun, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
  false,   //!< < s_dai
  false,   //!< < s_lares
  false,   //!< < s_dai_update
  false,   //!< < s_lares_update
  0,       //!< < s_dai_p
  0,       //!< < s_dai_n
  0,       //!< < s_lares_n
  0,       //!< < s_dai_p_prev
  0,       //!< < s_dai_n_prev
  0        //!< < s_lares_n_prev
};

constexpr StressParam StressParam_User{
  2.1e-6,  //!< < omegap, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
  3.17e-6, //!< < omegan, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
  138.73,  //!< < Ep, from Salco de Vasconcelos, Xu, Li, Zhao, Grid indentation analysis of mechanical properties of composite electrodes in Li-ion batteries, Extreme Mechanics Letters 9 (3), 2016
  10,      //!< < En, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
  0.3,     //!< < nup, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
  0.3,     //!< < nun, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
  false,   //!< < s_dai
  false,   //!< < s_lares
  false,   //!< < s_dai_update
  false,   //!< < s_lares_update
  0,       //!< < s_dai_p
  0,       //!< < s_dai_n
  0,       //!< < s_lares_n
  0,       //!< < s_dai_p_prev
  0,       //!< < s_dai_n_prev
  0        //!< < s_lares_n_prev
};

constexpr StressParam StressParam_Fit{
  2.1e-6,  //!< < omegap, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
  3.17e-6, //!< < omegan, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
  138.73,  //!< < Ep, from Salco de Vasconcelos, Xu, Li, Zhao, Grid indentation analysis of mechanical properties of composite electrodes in Li-ion batteries, Extreme Mechanics Letters 9 (3), 2016
  10,      //!< < En, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
  0.3,     //!< < nup, from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
  0.3,     //!< < nun, from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
  false,   //!< < s_dai
  false,   //!< < s_lares
  false,   //!< < s_dai_update
  false,   //!< < s_lares_update
  0,       //!< < s_dai_p
  0,       //!< < s_dai_n
  0,       //!< < s_lares_n
  0,       //!< < s_dai_p_prev
  0,       //!< < s_dai_n_prev
  0        //!< < s_lares_n_prev
};
} // namespace slide::param::def