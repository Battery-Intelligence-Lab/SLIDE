/*
 * unit_tests.hpp
 *
 *  Created on: 26 Nov 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#ifndef SRC_UNITTESTS_H_
#define SRC_UNITTESTS_H_

#include <iostream>

namespace slide::unit_tests
{
    void test_Cell_Bucket(bool testErrors);             // Bucket Cell
    void test_Cell_ECM(bool testErrors);                // ECM cell
    void test_Cell_SPM(bool testErrors);                // SPM cell
    void test_Module_s_testAll(bool testErrors);        // Module_s (series connected module)
    void test_Module_p_testAll(bool testErrors);        // Module_p (parallel connected base module)
    void test_Cycler(bool testErrors, int coolcontrol); // Cycler
    void test_Procedure(int coolcontrol);               // procedure
    void test_Converter();                              // power electronic converter
    void test_Battery();                                // Battery

    void inline test_all()
    {
        std::cout << "Unit tests are starting!\n";
        test_Cell_Bucket(true);
    }

}
#endif /* SRC_UNITTESTS_H_ */
