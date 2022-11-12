/*
 * unit_tests.hpp
 *
 *  Created on: 26 Nov 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#ifndef SRC_UNITTESTS_H_
#define SRC_UNITTESTS_H_

#include <iostream>

namespace slide::unit_tests {

bool TEST(auto &&fun, auto &&fun_name)
{
  try {
    fun();
  } catch (...) {
    std::cout << fun_name << " is successful!\n";
    return false;
  }
  std::cout << fun_name << " is failed!\n";
  return true;
}

bool test_Cell_Bucket();              //!< Bucket Cell
bool test_Cell_ECM();                 //!< ECM cell
bool test_Cell_SPM();                 //!< SPM cell
bool test_Module_s_testAll();         //!< Module_s (series connected module)
bool test_Module_p_testAll();         //!< Module_p (parallel connected base module)
bool test_Cycler(int coolcontrol);    //!< Cycler
bool test_Procedure(int coolcontrol); //!< procedure
bool test_Converter();                //!< power electronic converter
bool test_Battery();                  //!< Battery


bool inline test_all()
{
  std::cout << "Unit tests are starting!\n";
  return test_Cell_Bucket();
}


} // namespace slide::unit_tests
#endif /* SRC_UNITTESTS_H_ */
