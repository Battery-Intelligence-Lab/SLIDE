/*
 * tests_source.cpp
 *
 * This file implements test functions
 *
 * Date: 		2022.08.11
 * Author(s): 	Jorn Reniers, Volkan Kumtepeli
 *
 */

//!< Include header files
#include <ctime>
#include <thread>
#include <array>
#include <tuple>
#include <random>
#include <cmath>
#include <iomanip>

#include "slide.hpp"


int main()
{

  //!< print that you start simulations
  //!< slide::unit_tests::test_all();

  std::cout << "Start tests" << std::endl;

  //!< Slide-pack tests:

  //!< *********************************************** END ********************************************************
  //!< Now all the simulations have finished. Print this message, as well as how long it took to do the simulations
  std::cout << "finished all tests in " << clk << ".\n";
}
