/*
 * Module_base_s_test.cpp
 *
 *  Created on: 9 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "../tests_util.hpp"
#include "../../src/slide.hpp"

#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <memory>
#include <typeinfo>

namespace slide::tests::unit {

//!< ********************************************************** test functions from Module_base *******************************************************************
bool test_BasicGetters()
{
  //!< int Module_base::getNcells()
  //!< int Module_base::getNstates()
  //!< double Module_base::getTmod()

  std::string n = "na";

  std::unique_ptr<StorageUnit> cs[] = { std::make_unique<Cell_Bucket>(), std::make_unique<Cell_Bucket>() };

  auto cp1{ cs[0].get() }, cp2{ cs[1].get() };

  double T = settings::T_ENV;
  bool checkCells = false;

  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);

  mp->setSUs(cs, checkCells, true);
  assert(mp->getNSUs() == std::size(cs));
  assert(mp->getNstates() == std::size(cs) * cp1->getNstates() + 1);
  assert(mp->T() == T);

  return true;
}
bool test_getCellV()
{
  //!< Module_base::getCellVotages

  std::unique_ptr<StorageUnit> cs[];

  auto cs[0] = std::make_unique<Cell_Bucket>();
  auto cs[1] = std::make_unique<Cell_Bucket>();
  std::string n = "na";
  cs[1]->setSOC(0.6); //!< change cell 2 to verify the order is correct

  double v[] = { cs[0]->V(), cs[1]->V() };
  constexpr double T = settings::T_ENV;
  constexpr bool checkCells = false;

  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  assert(mp->Vi(0) == v[0]);
  assert(mp->Vi(1) == v[1]);
}
bool test_getStates()
{
  //!< bool Module_base::getStates(double s[], int nin, int& nout)

  constexpr int nin = settings::STORAGEUNIT_NSTATES_MAX;
  constexpr int nin1 = settings::CELL_NSTATE_MAX;
  double s[nin], s1[nin1], s2[nin1];
  int nout, nout1, nout2;
  auto cp1 = std::make_unique<Cell_Bucket>();
  auto cp2 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> cs[] = { cp1, cp2 };
  std::string n = "na";
  cp2->setSOC(0.6); //!< change cell 2 to verify the order is correct
  cp1->getStates(s1, nin1, nout1);
  cp2->getStates(s2, nin1, nout2);
  double T = settings::T_ENV;
  bool checkCells = false;

  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);
  mp->getStates(s, nin, nout);

  assert(nout == 1 + nout1 + nout2);
  assert(s[nout - 1] == T);
  for (int i = 0; i < nout1; i++) {
    assert(s1[i] == s[i]);
    assert(s2[i] == s[nout1 + i]);
  }
}
bool test_getCells()
{


  int nin1 = settings::CELL_NSTATE_MAX;
  int nn;
  double so1[nin1], so2[nin1]; //!< original states
  double s1[nin1], s2[nin1];   //!< states of returned cells
  double v1, v2;

  auto cp1 = std::make_unique<Cell_Bucket>();
  auto cp2 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> cs[] = { cp1, cp2 };
  std::string n = "na";
  cp2->setSOC(0.6); //!< change cell 2 to verify the order is correct
  cp1->getStates(so1, nin1, nn);
  cp2->getStates(so2, nin1, nn);
  v1 = cp1->V();
  v2 = cp2->V();
  double T = settings::T_ENV;
  bool checkCells = false;

  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);
  std::unique_ptr<StorageUnit> cs2[];
  mp->getSUs(cs2);                 //!< returns storage units, not cells. So getStates is the implementation from SU not Cell_Bucket
  cs2[0]->getStates(s1, nin1, nn); //!< this is SU::getStates, not Cell_Bucket::getStates
  cs2[1]->getStates(s2, nin1, nn);

  //!< validate cell states
  for (int i = 0; i < cp1->getNstates(); i++) {
    assert(s1[i] == so1[i]);
    assert(s2[i] == so2[i]);
  }
  assert(cs2[0]->V() == cp1->V());
  assert(cs2[1]->V() == cp2->V());

  /*
   * Side step about type checking and dynamic casting.
   * See https://stackoverflow.com/questions/1347691/static-vs-dynamic-type-checking-in-c
   */
  /*cout<<"Dynamic type checking: "<<endl;
  std::cout << typeid( cs2[0] ).name() << '\n';   //prints P11StorageUnit
  std::cout << typeid( *cs2[0] ).name() << '\n';   //prints St4Cell

  std::cout<<"remove the first three characters from the string to get the name of the class"<<endl;
  std::string a = typeid( *cs2[0] ).name();
  a.erase(a.begin() + 0, a.begin()+3); 	//!< remove the first 3 random characters, to end up with the string "Class" or "StorageUnit"
  std::cout<<a<<endl;

  std::cout<<"Check typeID with respect to class"<<endl;
  if(typeid(*cs2[0]) == typeid(Cell_Bucket))
          std::cout<<"It is a cell after all"<<endl;
  if(typeid(*cs2[0]) != typeid(Module_base))
          std::cout<<"It is not a base module"<<endl;

  std::cout<<"Dynamic type casting to the correct type: "<<endl;
  Cell_Bucket* cell1 = dynamic_cast<Cell_Bucket*>(cs2[0]);
  std::cout<<cell1->SOC()<<endl;				//!< prints correctly to 0.5
  std::cout<<"Dynamic type casting to the wrong type: "<<endl;
  Module_base_s* mb1 = dynamic_cast<Module_base_s*>(cs2[0]);
  if(mb1 == NULL)
          std::cout<<"wrong type, returned null"<<endl;

  //!< casting smart pointers
  std::unique_ptr<StorageUnit> c2 = c1.copy(); 		//!< c2 is a pointer to a SU
  Cell_Bucket* c22 = dynamic_cast<Cell_Bucket*>(c2.get()); 	//!< Dynamic cast from StorageUnit to Cell_Bucket
   */

  //!< Set the current, ensure it has changed in the cells
  double Inew = 1;
  mp->setCurrent(Inew);
  mp->getSUs(cs2);
  assert(cs2[0]->I() == Inew);
  assert(cs2[1]->I() == Inew);
  assert(cs2[0]->V() < v1);
  assert(cs2[1]->V() < v2);

  //!< charge
  Inew = -1;
  mp->setCurrent(Inew);
  mp->getSUs(cs2);
  assert(cs2[0]->I() == Inew);
  assert(cs2[1]->I() == Inew);
  assert(cs2[0]->V() > v1);
  assert(cs2[1]->V() > v2);

  return true;
}
bool test_setT()
{
  //!< bool Module_base::setT(double Tnew)


  auto cp1 = std::make_unique<Cell_Bucket>();
  auto cp2 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> cs[] = { cp1, cp2 };
  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  double Tnew = 273;
  mp->setT(Tnew);
  assert(mp->T() == Tnew);

  return true;
}
bool test_setStates(bool testError)
{
  //!< double Module_base::setStates(double s[], int nin, bool checkStates, bool print)


  auto cp1 = std::make_unique<Cell_Bucket>();
  auto cp2 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> cs[] = { cp1, cp2 };
  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  //!< Make a new State vector from cells 1 and 3, where 3 has a different SOC
  int nin = settings::STORAGEUNIT_NSTATES_MAX;
  int nin1 = settings::CELL_NSTATE_MAX;
  double s[nin], sout[nin], s1[nin1], s2[nin1];
  int nout, nout1;
  Cell_Bucket c3;
  c3.setSOC(0.4);        //!< change SOC of a cell
  double Tnew = 273 + 1; //!< change T
  cp1->getStates(s1, nin1, nout1);
  c3.getStates(s2, nin1, nout1);
  for (int i = 0; i < nout1; i++) {
    s[i] = s1[i];
    s[nout1 + i] = s2[i];
  }
  s[nout1 * 2] = Tnew;
  mp->setStates(s, 1 + 2 * nout1); //!< this changes the states of m, and should therefore change the states of cell 2 too

  //!< Check cell 2 and the states have changed
  assert(cp2->SOC() == 0.4); //!< mp->setStates invoked cp2->setStates, so also c2 should have changed
  mp->getStates(sout, nin, nout);
  assert(nout == 1 + 2 * nout1);
  assert(sout[nout - 1] == Tnew);
  for (int i = 0; i < nout1; i++) {
    assert(s1[i] == sout[i]);
    assert(s2[i] == sout[nout1 + i]);
  }

  if (testError) {
    //!< try setting an array with the wrong length
    try {
      mp->setStates(s1, nout1);
      assert(false);
    } catch (...) {
    };

    //!< try setting states with different currents
    c3.setCurrent(1);
    cp1->getStates(s1, nin1, nout1);
    c3.getStates(s2, nin1, nout1);
    s[2 * nout1] = Tnew;
    for (int i = 0; i < nout1; i++) {
      s[i] = s1[i];
      s[nout1 + i] = s2[i];
    }
    try {
      mp->setStates(s, 1 + 2 * nout1);
      assert(false);
    } catch (...) {
    };
  }
}
bool test_setCells()
{
  //!< bool Module_base::setCells(Cell_Bucket c[], int nin, bool checkCells, bool print)


  auto cp1 = std::make_unique<Cell_Bucket>();
  auto cp2 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> cs[] = { cp1, cp2 };
  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells);

  //!< set cell with new SOC
  double socnew = 0.4;
  auto cp22 = std::make_unique<Cell_Bucket>();
  cp22->setSOC(socnew);
  cs[1] = cp22;
  mp->setSUs(cs);
  //!< check using getCells
  std::unique_ptr<StorageUnit> cs2[];
  int nn;
  mp->getSUs(cs2);
  Cell_Bucket *cell1 = dynamic_cast<Cell_Bucket *>(cs2[1].get()); //!< Dynamic cast from smart ptr of StorageUnit to regular pointer of Cell_Bucket
  assert(cell1->SOC() == socnew);                                 //!< the new soc should be 0.1
  assert(cp2->SOC() == 0.5);                                      //!< we replaced the pointer in mp->cells[1] which now points to c22 (with new SOC). c2 should not have changed

  //!< try updating the number of cells
  auto cp3 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> cs3[] = { cp1, cp22, cp3 };
  mp->setSUs(cs3);

  //!< check using getCells
  std::unique_ptr<StorageUnit> cs4[];
  mp->getSUs(cs4);
  cell1 = dynamic_cast<Cell_Bucket *>(cs4[1].get()); //!< Dynamic cast from StorageUnit to Cell_Bucket
  assert(cell1->SOC() == socnew);                    //!< cells[1] is cell22 with new soc
  cell1 = dynamic_cast<Cell_Bucket *>(cs4[2].get()); //!< Dynamic cast from StorageUnit to Cell_Bucket
  assert(cell1->SOC() == 0.5);                       //!< cells[2] is cell3 with standard soc
  assert(mp->getNSUs() == std::size(cs3));

  if () {
    //!< try setting cells with different currents
    double Inew = 1;
    cp2->setCurrent(Inew);
    std::unique_ptr<StorageUnit> cs5[] = { cp1, cp2, cp3 };
    try {
      mp->setSUs(cs5); //!< this causes some erorror. 10 is thrown, but we don't arrive in the catch-statement
      assert(false);   //!< so probably one of the smart pointer shits is fucking up
    } catch (...) {
    }

    /*	//!< cells which belong to a different module
* note: that gives errors with the pointers
    auto mp2  = std::make_unique<Module_s>();
    try{
            mp2->setSUs(cs,std::size(cs),false); //this has to throw 10
            assert(false);
    }catch(int e){};*/
  }

  return true;
}

//!< ********************************************************* test functions from Module_base_s *******************************************************************
bool test_Constructor()
{
  //!< Module_base_s::Module_base_s()
  //!< Module_base_s::Module_base_s(int ncellsi, Cell_Bucket ci[], double Ti, bool checkCells, bool print)

  auto mp = std::make_unique<Module_s>();
  assert(mp->getNSUs() == 0);
  assert(mp->T() == settings::T_ENV);


  auto cp1 = std::make_unique<Cell_Bucket>();
  auto cp2 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> cs[] = { cp1, cp2 };
  assert(cp1->getID() == "cell");
  assert(cp1->getFullID() == "cell"); //!< has no parent yet
  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp2 = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp2->setSUs(cs, checkCells, true);
  //!< note: you have to make a new module, else the getParent() fails
  //!< I made an m before, so &m points to the first module_s made

  assert(mp2->getNSUs() == std::size(cs));
  assert(mp2->T() == T);
  assert(mp2->getID() == n);
  assert(cp1->getParent() == mp2);
  assert(cp2->getParent() == mp2);
  assert(cp1->getID() == "cell");
  assert(cp1->getFullID() == "na_cell");
}
bool test_BasicGetters_s()
{
  auto cp1 = std::make_unique<Cell_Bucket>();
  auto cp2 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> cs[] = { cp1, cp2 };
  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  assert(mp->Cap() == cp1->Cap());
  assert(mp->Vmin() == 2 * cp1->Vmin());
  assert(mp->Vmax() == 2 * cp1->Vmax());
  assert(mp->I() == 0);
  assert(mp->V() == 2 * cp1->V());
}
bool test_setI()
{
  //!< double Module_base_s::setCurrent(double Inew, bool checkV, bool print)


  auto cp1 = std::make_unique<Cell_Bucket>();
  auto cp2 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> cs[] = { cp1, cp2 };
  std::string n = "na";
  double v1 = cp1->V();
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);
  assert(mp->I() == 0);
  assert(mp->V() == 2 * v1);

  double Inew = 1;
  double V;
  mp->setCurrent(Inew, true);
  V = mp->V();
  assert(mp->I() == Inew);
  assert(V < 2 * v1); //!< voltage must decrease
  //!< do not check individual cells, that is done in getCells

  Inew = -1;
  mp->setCurrent(Inew, true);
  V = mp->V();
  assert(mp->I() == Inew);
  assert(mp->V() > 2 * v1); //!< voltage must increase

  //!< test things which should break
  Inew = 10000;                         //!< very large current, should give too low voltage
  Status status = mp->setCurrent(Inew); //!< should fail because the current equation cannot be solved
  if (isStatusSuccessful(status)) return false;


  Inew = -10000;                        //!< very large current, should give too low voltage
  Status status = mp->setCurrent(Inew); //!< should fail because the current equation cannot be solved
  if (isStatusSuccessful(status)) return false;

  return true;
}
bool test_validStates()
{
  //!< bool Module_base_s::validStates(double s[], int nin, bool print)


  auto cp1 = std::make_unique<Cell_Bucket>();
  auto cp2 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> cs[] = { cp1, cp2 };
  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  //!< valid states (current states)
  int nin = settings::MODULE_NSUs_MAX;
  double s[nin];
  int nout;
  mp->getStates(s, nin, nout);
  assert(mp->validStates(s, nout));

  //!< valid states (new T)
  s[nout - 1] = 273 + 5;
  assert(mp->validStates(s, nout));

  if () {
    //!< wrong length
    int nc = settings::CELL_NSTATE_MAX;
    double sc[nc];
    int noutc;
    cp1->getStates(sc, nc, noutc);
    assert(!mp->validStates(sc, noutc, false));

    //!< an SOC which is too large
    s[0] = 2; //!< this is the SOC of cell 1
    s[nout - 1] = 273 + 1;
    assert(!mp->validStates(s, nout, false));

    //!< different current values
    s[0] = 0.5;
    s[2] = 5; //!< current of cell 1
    assert(!mp->validStates(s, nout, false));
  }
}
bool test_validCells()
{
  //!< bool Module_base_s::validCells(Cell_Bucket c[], int nin, bool print)

  auto cp1 = std::make_unique<Cell_Bucket>();
  auto cp2 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> cs[] = { cp1, cp2 };
  std::string n = "na";
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  //!< valid cells with the current cells
  std::unique_ptr<StorageUnit> cs2[];
  int nout;
  mp->getSUs(cs2);
  assert(mp->validSUs(cs2));

  //!< valid cells (new cell)
  auto cp3 = std::make_unique<Cell_Bucket>();
  cs2[1] = cp3;
  assert(mp->validSUs(cs2));

  //!< add an additional cell
  std::unique_ptr<StorageUnit> cs3[] = { cp1, cp2, cp3 };
  assert(mp->validSUs(cs3, 3));

  if () {
    //!< different currents
    cp3->setCurrent(1);
    assert(cs3[2]->I() == 1); //!< array of pointers, so it should have been changed automatically in the array
    assert(!mp->validSUs(cs3, 3, false));

    //!< validSU does not check whether cells are outside their voltage range
  }

  return true;
}
bool test_timeStep_CC()
{
  //!< bool Module_base_s::timeStep_CC(double dt)


  auto cp1 = std::make_unique<Cell_Bucket>();
  auto cp2 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> cs[] = { cp1, cp2 };
  std::string n = "na";
  double v1 = cp1->V();
  double soc1 = cp1->SOC();
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);

  //!< time step with 0 current
  double dt = 5;
  mp->timeStep_CC(dt);
  assert(mp->V() == 2 * v1);

  //!< discharge
  double Inew = 1;
  double V, err;
  double tol = 0.002;
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  V = mp->V();
  assert(V < 2 * v1);
  assert(mp->I() == Inew);

  //!< check individual cells
  std::unique_ptr<StorageUnit> cs2[];
  int nout;
  Cell_Bucket *cell1;
  mp->getSUs(cs2, nout);
  for (int i = 0; i < std::size(cs2); i++) {
    assert(cs2[i]->I() == Inew);
    assert(cs2[i]->V() < v1);
    cell1 = dynamic_cast<Cell_Bucket *>(cs2[i].get()); //!< Dynamic cast from StorageUnit to Cell_Bucket
    assert(cell1->SOC() < soc1);
    err = cell1->SOC() - (0.5 - 1.0 * 5.0 / 3600.0 / cell1->Cap());
    assert(err < tol && err > -tol);
  }

  //!< charge
  Inew = -1;
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  mp->getSUs(cs2);
  assert(mp->V() > V);
  assert(mp->I() == Inew);
  for (int i = 0; i < std::size(cs2); i++) {
    assert(cs2[i]->I() == Inew);
    assert(cs2[i]->V() > V / 2.0);
    cell1 = dynamic_cast<Cell_Bucket *>(cs2[i].get()); //!< Dynamic cast from StorageUnit to Cell_Bucket
    err = cell1->SOC() - (0.5);
    assert(err < tol && err > -tol);
  }

  return true;
}

bool test_Modules_s_ECM()
{
  //!< test series modules make out of ECM cells
  //!< this just repeats the other tests but with a different Cell_Bucket type

  int nin = settings::STORAGEUNIT_NSTATES_MAX;
  int nin1 = settings::CELL_NSTATE_MAX;
  double s[nin], sout[nin], s1[nin1], s2[nin1];
  int nout, nout1, nout2;

  auto cp1 = std::make_unique<Cell_ECM>();
  auto cp2 = std::make_unique<Cell_ECM>();
  auto cp3 = std::make_unique<Cell_ECM>();
  auto cp22 = std::make_unique<Cell_ECM>();
  std::unique_ptr<StorageUnit> cs[] = { cp1, cp2 };
  Cell_ECM *cell1;
  std::string n = "na";

  //!< getStates - setStates
  cp2->setSOC(0.6); //!< change cell 2 to verify the order is correct
  cp1->getStates(s1, nin1, nout1);
  cp2->getStates(s2, nin1, nout2);
  double T = settings::T_ENV;
  bool checkCells = true;
  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);
  mp->getStates(s, nin, nout);
  assert(nout == 1 + nout1 + nout2);
  assert(s[nout - 1] == T);
  for (int i = 0; i < nout1; i++) {
    assert(s1[i] == s[i]);
    assert(s2[i] == s[nout1 + i]);
  }
  cp3->setSOC(0.4);      //!< change SOC of a cell
  double Tnew = 273 + 1; //!< change T
  cp1->getStates(s1, nin1, nout1);
  cp3->getStates(s2, nin1, nout1);
  for (int i = 0; i < nout1; i++) {
    s[i] = s1[i];
    s[nout1 + i] = s2[i];
  }
  s[nout1 * 2] = Tnew;
  mp->setStates(s, 1 + 2 * nout1);
  assert(cp2->SOC() == 0.4); //!< mp->setStates invoked cp2->setStates, so also c2 should have changed
  mp->getStates(sout, nin, nout);
  assert(nout == 1 + 2 * nout1);
  assert(sout[nout - 1] == Tnew);
  for (int i = 0; i < nout1; i++) {
    assert(s1[i] == sout[i]);
    assert(s2[i] == sout[nout1 + i]);
  }
  mp->setT(settings::T_ENV); //!< reset the tenperature, else the thermal model in timeStep will freak out

  //!< getCells - setCells
  std::unique_ptr<StorageUnit> cs2[];
  double socnew = 0.4;
  cp22->setSOC(socnew);
  cp1->setSOC(0.5);
  cp2->setSOC(0.5);
  cp3->setSOC(0.5);
  cs[1] = cp22;
  mp->setSUs(cs);
  //!< check using getCells
  int nn;
  mp->getSUs(cs2);
  cell1 = dynamic_cast<Cell_ECM *>(cs2[1].get()); //!< Dynamic cast from StorageUnit to Cell_Bucket
  assert(cell1->SOC() == socnew);                 //!< the new soc should be 0.1
  assert(cp2->SOC() == 0.5);                      //!< we replaced the pointer in mp->cells[1] which now points to c22 (with new SOC). c2 should not have changed
  //!< try updating the number of cells
  std::unique_ptr<StorageUnit> cs3[] = { cp1, cp22, cp3 };
  mp->setSUs(cs3);
  //!< check using getCells
  std::unique_ptr<StorageUnit> cs4[];
  mp->getSUs(cs4);
  cell1 = dynamic_cast<Cell_ECM *>(cs4[1].get()); //!< Dynamic cast from StorageUnit to Cell_Bucket
  assert(cell1->SOC() == socnew);                 //!< cells[1] is cell22 with new soc
  cell1 = dynamic_cast<Cell_ECM *>(cs4[2].get()); //!< Dynamic cast from StorageUnit to Cell_Bucket
  assert(cell1->SOC() == 0.5);                    //!< cells[2] is cell3 with standard soc
  assert(mp->getNSUs() == std::size(cs3));

  //!< timeStep_CC
  auto cp111 = std::make_unique<Cell_ECM>();
  auto cp222 = std::make_unique<Cell_ECM>();
  cp1 = cp111;
  cp2 = cp222;
  cs[0] = cp1;
  cs[1] = cp2;
  double v1 = cp1->V();
  double soc1 = cp1->SOC();
  mp->setSUs(cs, checkCells, true);
  //!< time step with 0 current
  double dt = 5;
  mp->timeStep_CC(dt);
  assert(mp->V() == 2 * v1);
  //!< discharge
  double Inew = 1;
  double tol = 0.002;
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  double V = mp->V();
  double err;
  assert(V < 2 * v1);
  assert(mp->I() == Inew);
  //!< check individual cells
  mp->getSUs(cs2);
  for (int i = 0; i < std::size(cs2); i++) {
    assert(cs2[i]->I() == Inew);
    assert(cs2[i]->V() < v1);
    cell1 = dynamic_cast<Cell_ECM *>(cs2[i].get()); //!< Dynamic cast from StorageUnit to Cell_Bucket
    assert(cell1->SOC() < soc1);
    err = cell1->SOC() - (0.5 - 1.0 * 5.0 / 3600.0 / cell1->Cap());
    assert(err < tol && err > -tol);
  }
  //!< charge
  Inew = -1;
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  mp->getSUs(cs2);
  assert(mp->V() > V);
  assert(mp->I() == Inew);
  for (int i = 0; i < std::size(cs2); i++) {
    assert(cs2[i]->I() == Inew);
    assert(cs2[i]->V() > V / 2.0);
    cell1 = dynamic_cast<Cell_ECM *>(cs2[i].get()); //!< Dynamic cast from StorageUnit to Cell_Bucket
    err = cell1->SOC() - (0.5);
    assert(err < tol && err > -tol);
  }

  return true;
}

bool test_Modules_s_SPM()
{
  //!< test series modules make out of ECM cells
  //!< this just repeats the other tests but with a different Cell_Bucket type

  int nin = settings::STORAGEUNIT_NSTATES_MAX;
  int nin1 = settings::CELL_NSTATE_MAX;
  double s[nin], sout[nin], s1[nin1], s2[nin1];
  int nout, nout1, nout2;

  auto cp1 = std::make_unique<Cell_SPM>();
  auto cp2 = std::make_unique<Cell_SPM>();
  auto cp3 = std::make_unique<Cell_SPM>();
  auto cp22 = std::make_unique<Cell_SPM>();
  std::unique_ptr<StorageUnit> cs[] = { cp1, cp2 };
  Cell_SPM *cell1;
  std::string n = "na";

  //!< getStates - setStates
  cp2->setSOC(0.6); //!< change cell 2 to verify the order is correct
  cp1->getStates(s1, nin1, nout1);
  cp2->getStates(s2, nin1, nout2);
  double T = settings::T_ENV;
  bool checkCells = true;
  auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  mp->setSUs(cs, checkCells, true);
  mp->getStates(s, nin, nout);
  assert(nout == 1 + nout1 + nout2);
  assert(s[nout - 1] == T);
  for (int i = 0; i < nout1; i++) {
    assert(s1[i] == s[i]);
    assert(s2[i] == s[nout1 + i]);
  }
  cp3->setSOC(0.1);      //!< change SOC of a cell
  double Tnew = 273 + 1; //!< change T
  cp1->getStates(s1, nin1, nout1);
  cp3->getStates(s2, nin1, nout1);
  for (int i = 0; i < nout1; i++) {
    s[i] = s1[i];
    s[nout1 + i] = s2[i];
  }
  s[nout1 * 2] = Tnew;
  mp->setStates(s, 1 + 2 * nout1);
  assert(cp2->SOC() == 0.1); //!< mp->setStates invoked cp2->setStates, so also c2 should have changed
  mp->getStates(sout, nin, nout);
  assert(nout == 1 + 2 * nout1);
  assert(sout[nout - 1] == Tnew);
  for (int i = 0; i < nout1; i++) {
    assert(s1[i] == sout[i]);
    assert(s2[i] == sout[nout1 + i]);
  }

  //!< getCells - setCells
  std::unique_ptr<StorageUnit> cs2[];
  /*	double socnew = 0.1;
  cp22->setSOC(socnew);
  cp1->setSOC(0.5);
  cp2->setSOC(0.5);
  cp3->setSOC(0.5);
  cs[1] = cp22;
  mp->setSUs(cs);
  //!< check using getCells
  int nn;
  mp->getSUs(cs2, nn);
  cell1 = dynamic_cast<Cell_SPM*>(cs2[1].get()); 	//!< Dynamic cast from StorageUnit to Cell_Bucket
  assert(cell1->SOC()==socnew);			//!< the new soc should be 0.1
  assert(cp2->SOC()==0.5);					//!< we replaced the pointer in mp->cells[1] which now points to c22 (with new SOC). c2 should not have changed
  //!< try updating the number of cells
  std::unique_ptr<StorageUnit> cs3[] = {cp1, cp22, cp3};
  mp->setSUs(cs3);
  //!< check using getCells
  std::unique_ptr<StorageUnit> cs4[];
  mp->getSUs(cs4, nn);
  cell1 = dynamic_cast<Cell_SPM*>(cs4[1].get()); 	//!< Dynamic cast from StorageUnit to Cell_Bucket
  assert(cell1->SOC()==socnew);		//!< cells[1] is cell22 with new soc
  cell1 = dynamic_cast<Cell_SPM*>(cs4[2].get()); 	//!< Dynamic cast from StorageUnit to Cell_Bucket
  assert(cell1->SOC()==0.5);			//!< cells[2] is cell3 with standard soc
  assert(mp->getNSUs()==std::size(cs4));
*/
  //!< timeStep_CC
  auto cp111 = std::make_unique<Cell_SPM>();
  auto cp222 = std::make_unique<Cell_SPM>();
  cp1 = cp111;
  cp2 = cp222;
  cs[0] = cp1;
  cs[1] = cp2;
  double v1 = cp1->V();
  double soc1 = cp1->SOC();
  mp->setSUs(cs, checkCells, true);
  //!< time step with 0 current
  double dt = 5;
  mp->timeStep_CC(dt);
  assert(mp->V() == 2 * v1);
  //!< discharge
  double Inew = 1;
  double tol = 0.002;
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  double V = mp->V();
  double err;
  assert(V < 2 * v1);
  assert(mp->I() == Inew);
  //!< check individual cells
  mp->getSUs(cs2, nout);
  for (int i = 0; i < std::size(cs2); i++) {
    assert(cs2[i]->I() == Inew);
    assert(cs2[i]->V() < v1);
    cell1 = dynamic_cast<Cell_SPM *>(cs2[i].get()); //!< Dynamic cast from StorageUnit to Cell_Bucket
    assert(cell1->SOC() < soc1);
    err = cell1->SOC() - (0.5 - 1.0 * 5.0 / 3600.0 / cell1->Cap());
    assert(err < tol && err > -tol);
  }
  //!< charge
  Inew = -1;
  mp->setCurrent(Inew);
  mp->timeStep_CC(dt);
  mp->getSUs(cs2, nout);
  assert(mp->V() > V);
  assert(mp->I() == Inew);
  for (int i = 0; i < std::size(cs2); i++) {
    assert(cs2[i]->I() == Inew);
    assert(cs2[i]->V() > V / 2.0);
    cell1 = dynamic_cast<Cell_SPM *>(cs2[i].get()); //!< Dynamic cast from StorageUnit to Cell_Bucket
    err = cell1->SOC() - (0.5);
    assert(err < tol && err > -tol);
  }
}

bool test_Hierarchichal()
{
  //!< test series modules made out of other series modules

  double tol = 0.001;
  std::string n1 = "H1";
  std::string n2 = "H2";
  std::string n3 = "H3";
  auto cp1 = std::make_unique<Cell_Bucket>();
  auto cp2 = std::make_unique<Cell_Bucket>();
  auto cp3 = std::make_unique<Cell_Bucket>();
  auto cp4 = std::make_unique<Cell_Bucket>();
  auto cp5 = std::make_unique<Cell_Bucket>();
  auto cp6 = std::make_unique<Cell_Bucket>();
  auto cp7 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> SU1[] = { cp1, cp2 };
  std::unique_ptr<StorageUnit> SU2[] = { cp3, cp4 };
  std::unique_ptr<StorageUnit> SU3[] = { cp5, cp6, cp7 };
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp1 = std::make_unique<Module_s>(ids[0], T, true, false, std::size(SU1), 1, 2); //!< middle level modules are pass through
  auto mp2 = std::make_unique<Module_s>(ids[1], T, true, false, std::size(SU2), 1, 2);
  auto mp3 = std::make_unique<Module_s>(ids[2], T, true, false, std::size(SU3), 1, 2);
  mp1->setSUs(SU1, checkCells);
  mp2->setSUs(SU2, checkCells);
  mp3->setSUs(SU3, checkCells);

  //!< make the hierarichical module
  int nm = 3;
  std::string n4 = "H4";
  std::unique_ptr<StorageUnit> MU[] = { mp1, mp2, mp3 };
  checkCells = true;
  auto mp = std::make_unique<Module_s>(n4, T, true, false, 7, 1, 1);
  mp->setSUs(MU, checkCells, true);
  double Vini = mp->V();
  assert(Vini == mp1->V() + mp2->V() + mp3->V());
  assert(Vini == 7 * cp1->V());
  assert(mp->getFullID() == "H4");
  assert(mp1->getFullID() == "H4_H1");
  assert(cp1->getFullID() == "H4_H1_cell");
  assert(cp4->getFullID() == "H4_H2_cell");
  assert(cp5->getFullID() == "H4_H3_cell");

  //!< set a CC current
  double Inew = -5;
  mp->setCurrent(Inew);
  assert(mp1->I() == Inew);
  assert(cp3->I() == Inew);
  assert(cp7->I() == Inew);

  //!< time a CC time step
  Vini = mp->V();
  double dt = 5;
  mp->timeStep_CC(dt);
  assert(std::abs(cp1->SOC() - (0.5 - Inew * dt / 3600.0 / cp1->Cap())) < tol); //!< the SOC must have increased (check just 1 cell out of all 7)
  assert(mp->V() > Vini);
}

bool test_Hierarchical_Cross()
{
  //!< test series-modules made out of parallel-modules

  double tol = settings::MODULE_P_I_ABSTOL;
  std::string ids[] = { "1", "2", "3" };
  auto cp1 = std::make_unique<Cell_Bucket>();
  auto cp2 = std::make_unique<Cell_Bucket>();
  auto cp3 = std::make_unique<Cell_Bucket>();
  auto cp4 = std::make_unique<Cell_Bucket>();
  auto cp5 = std::make_unique<Cell_Bucket>();
  auto cp6 = std::make_unique<Cell_Bucket>();
  auto cp7 = std::make_unique<Cell_Bucket>();
  std::unique_ptr<StorageUnit> SU1[] = { cp1, cp2 };
  std::unique_ptr<StorageUnit> SU2[] = { cp3, cp4 };
  std::unique_ptr<StorageUnit> SU3[] = { cp5, cp6, cp7 };
  double T = settings::T_ENV;
  bool checkCells = false;
  auto mp1 = std::make_unique<Module_p>(ids[0], T, true, false, std::size(SU1), 1, 2); //!< middle level modules are pass through
  auto mp2 = std::make_unique<Module_p>(ids[1], T, true, false, std::size(SU2), 1, 2);
  auto mp3 = std::make_unique<Module_p>(ids[2], T, true, false, std::size(SU3), 1, 2);
  mp1->setSUs(SU1, checkCells);
  mp2->setSUs(SU2, checkCells);
  mp3->setSUs(SU3, checkCells);

  //!< make the hierarichical series-module
  int nm = 3;
  std::string n4 = "4";
  std::unique_ptr<StorageUnit> MU[] = { mp1, mp2, mp3 };
  checkCells = true;
  auto mp = std::make_unique<Module_s>(n4, T, true, false, 7, 1, 1);
  mp->setSUs(MU, checkCells, true);
  double Vini = mp->V();
  assert(Vini == mp1->V() + mp2->V() + mp3->V());
  assert(Vini == 3 * cp1->V());

  //!< set a CC current
  double Inew = -5;
  mp->setCurrent(Inew);
  assert(NEAR(mp->I(), Inew, tol));                   //!< no longer exact since it involves parallel modules
  assert(NEAR(cp3->I(), Inew / mp2->getNSUs(), tol)); //!< 2nd module has 2 cells, so current should split in 2
  assert(NEAR(cp7->I(), Inew / mp3->getNSUs(), tol)); //!< 3rd module has 3 cells, so current should split in 3

  //!< time a CC time step
  Vini = mp->V();
  double dt = 5;
  mp->timeStep_CC(dt);
  assert(std::abs(cp1->SOC() - (0.5 - Inew * dt / 3600.0 / cp1->Cap())) < tol); //!< the SOC must have increased (check just 1 cell out of all 7)
  assert(mp->V() > Vini);

  return true;
}

bool test_copy_s()
{
  //!< 	/*
  //!< 	 * test the copy-function
  //!< 	 */

  //!< 	//!< make module
  //!<
  //!< 	auto cp1 = std::make_unique<Cell_Bucket>();
  //!< 	auto cp2 = std::make_unique<Cell_Bucket>();
  //!< 	std::unique_ptr<StorageUnit> cs[] = {cp1, cp2};
  //!< 	std::string n = "na";
  //!< 	double v1 = cp1->V();
  //!< 	double T = settings::T_ENV;
  //!< 	bool checkCells = false;
  //!< 	auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), 1, 1);
  //!< 	mp->setSUs(cs, checkCells, true);

  //!< 	//!< copy this one and check they are identical
  //!< 	std::unique_ptr<StorageUnit> cn = mp->copy();
  //!< 	Module_s *c22 = dynamic_cast<Module_s *>(cn.get()); //!< Dynamic cast from StorageUnit to Cell_Bucket
  //!< 	assert(mp->V() == c22->V());
  //!< 	std::unique_ptr<StorageUnit> corig[], cnew[];
  //!< 	int nout;
  //!< 	mp->getSUs(corig);
  //!< 	c22->getSUs(cnew);
  //!< 	for (int i = 0; i < mp->getNSUs(); i++)
  //!< 		assert(corig[i]->V() == cnew[i]->V());

  //!< 	//!< change the copied version, and ensure the old one is still the same
  //!< 	c22->setCurrent(1 * std::size(cs), false, false); //!< discharge
  //!< 	for (int t = 0; t < 10; t++)
  //!< 		c22->timeStep_CC(2);
  //!< 	mp->getSUs(corig);
  //!< 	c22->getSUs(cnew);
  //!< 	for (int i = 0; i < mp->getNSUs(); i++)
  //!< 	{
  //!< 		assert(corig[i]->V() == v1);
  //!< 		assert(cnew[i]->V() < v1);
  //!< 	}
  return true;
}

bool test_CoolSystem_s()
{
  /*
   * test the cool system design
   *
   * The difference between coolsystems is almost entirely due to different heat generation in cells.
   * Difference in heat generation between cool1 and cool 5 is about
   * 	1% for single cell
   * 	2.9% for simple module (4 cells)
   * 	5.5% for complex module (7 cells)
   *
   * Using print statements in Cell_SPM::dState_thermal for cool 5 vs cool 1
   * 		Ohmic heat generation is exactly the same (0.326687)
   * 		Reaction heat generation is about 1% higher (1.18595 vs 1.17032)
   * 			eta_n 0.7% larger
   * 			eta_p 3.8% larger [in absolute value]
   * 			cell temperature 0.2 % lower
   * 		Reversable heat generation is about 0.1% higher (0.0880738 vs 0.0879128)
   * so due to different cell T (here cool 5 is colder) different overpotentials (here cool 5 larger overpotential) and therefore different heat generation (here cool 5 has a bit more heat generation)
   *
   * for TENV = 30 compared to 20 (cool1)
   * 	T = 30 has 16.5 % less heat generation (in line with the previous ones, lower overpotentials)
   * 	Qcool and Qheat both decrease by same fraction
   */

  //!< Parameters from Cell_SPM which are needed to calculate the heat balance
  double rho = 1626;
  double Cp = 750;
  double L = 1.6850 * pow(10, -4);    //!< thickness of one layer
  double width = 0.1;                 //!< width of the pouch
  double height = 0.2;                //!< height of the pouch
  int nlayers = 31;                   //!< number of layers in the pouch
  double Acell = width * height;      //!< geometric surface area of the pouch
  double elec_surf = Acell * nlayers; //!< total 'unrolled' surface area of the electrodes

  //!< General settings
  double T = settings::T_ENV;
  bool checkCells = true;
  double Icha, Idis;
  double Irest = 0;
  double dt = 2;
  int N = 100; //!< number of steps per charge / discharge / rest

  //!< Loop for each setting of the cool controller
  for (int coolControl = 1; coolControl < 6; coolControl++) {

    //!< ****************************************************************************************************************************************************
    //!< Make a simple module with one SPM cell
    //!< cout<<"module_s_test start coolsystem test with a single cell for cool control setting "<<coolControl<<endl;
    auto cp0 = std::make_unique<Cell_SPM>();
    std::unique_ptr<StorageUnit> cs[] = { cp0 };
    std::string n = "testCoolSystem";
    auto mp = std::make_unique<Module_s>(n, T, true, false, std::size(cs), coolControl, true);
    mp->setSUs(cs, checkCells, true);
    double Tini[1] = { cp0->T() };

    //!< do a few 1C cycles (note just some time steps since we don't have the Cycler
    Icha = -cp0->Cap();
    Idis = cp0->Cap();
    for (int i = 0; i < 5; i++) {
      //!< charge
      mp->setCurrent(Icha);
      for (int t = 0; t < N; t++)
        mp->timeStep_CC(dt);
      //!< rest
      mp->setCurrent(Irest);
      for (int t = 0; t < N; t++)
        mp->timeStep_CC(dt);
      //!< discharge
      mp->setCurrent(Idis);
      for (int t = 0; t < N; t++)
        mp->timeStep_CC(dt);
      //!< rest
      mp->setCurrent(Irest);
      for (int t = 0; t < N; t++)
        mp->timeStep_CC(dt);
    }

    //!< check the energy balance of the outer module
    double Qgen = cp0->thermal_getTotalHeat();         //!< total heat generated by cells
    double Qcool = mp->getCoolSystem()->getHeatEvac(); //!< total heat extracted by the coolsystem from the cells
    double Tnew[1] = { cp0->T() };
    double Qheat = 0; //!< total energy in heating up the cells
    for (int i = 0; i < 1; i++)
      Qheat += (Tnew[i] - Tini[i]) * (rho * Cp * L * elec_surf);
    //!< cout<<"\t Total heat balance of coolsystem single cell "<<coolControl<<" is "<<Qgen<<", "<<Qheat<<", "<<Qcool<<" and error "<<abs(Qgen - Qcool - Qheat)<<endl<<flush;
    //!<  note: Qheat < 0 means the cell has cooled down compared to its initial value
    //!<  	that is possible because cells are made at T_env (which is 20), but the HVAC coolsystem can cool the SU as cold as controlAC_onoff_Toff, which is 15
    //!<  	so cells can cool down 5 degrees
    assert(std::abs(Qgen - Qcool - Qheat) / std::abs(Qgen) < 1e-10);

    //!< **********************************************************************************************************************************************************
    //!< Make a simple module with SPM cells
    //!< cout<<"module_s_test start coolsystem test with a simple module for cool control setting "<<coolControl<<endl;
    auto cp1 = std::make_unique<Cell_SPM>();
    auto cp2 = std::make_unique<Cell_SPM>();
    auto cp3 = std::make_unique<Cell_SPM>();
    auto cp4 = std::make_unique<Cell_SPM>();
    std::unique_ptr<StorageUnit> cs2[] = { cp1, cp2, cp3, cp4 };
    std::string n2 = "testCoolSystem";
    auto mp2 = std::make_unique<Module_s>(n2, T, true, false, std::size(cs2), coolControl, true);
    mp2->setSUs(cs2, checkCells, true);
    double Tini2[4] = { cp1->T(), cp2->T(), cp3->T(), cp4->T() };

    //!< do a few 1C cycles (note just some time steps since we don't have the Cycler
    Icha = -cp1->Cap();
    Idis = cp1->Cap();
    for (int i = 0; i < 5; i++) {
      //!< charge
      mp2->setCurrent(Icha);
      for (int t = 0; t < N; t++)
        mp2->timeStep_CC(dt);
      //!< rest
      mp2->setCurrent(Irest);
      for (int t = 0; t < N; t++)
        mp2->timeStep_CC(dt);
      //!< discharge
      mp2->setCurrent(Idis);
      for (int t = 0; t < N; t++)
        mp2->timeStep_CC(dt);
      //!< rest
      mp2->setCurrent(Irest);
      for (int t = 0; t < N; t++)
        mp2->timeStep_CC(dt);
    }

    //!< check the energy balance of the outer module
    double Qgen2 = cp1->thermal_getTotalHeat() + cp2->thermal_getTotalHeat() + cp3->thermal_getTotalHeat() + cp4->thermal_getTotalHeat(); //!< total heat generated by cells
    double Qcool2 = mp2->getCoolSystem()->getHeatEvac();                                                                                  //!< total heat extracted by the coolsystem from the cells
    double Tnew2[4] = { cp1->T(), cp2->T(), cp3->T(), cp4->T() };
    double Qheat2 = 0; //!< total energy in heating up the cells
    for (int i = 0; i < 4; i++)
      Qheat2 += (Tnew2[i] - Tini2[i]) * (rho * Cp * L * elec_surf);
    //!< cout<<"\t Total heat balance of coolsystem simple module "<<coolControl<<" is "<<Qgen2<<", "<<Qheat2<<", "<<Qcool2<<" and error "<<abs(Qgen2 - Qcool2 - Qheat2)<<endl<<flush;
    //!<  note: Qheat < 0 means the cell has cooled down compared to its initial value
    //!<  	that is possible because cells are made at T_env (which is 20), but the HVAC coolsystem can cool the SU as cold as controlAC_onoff_Toff, which is 15
    //!<  	so cells can cool down 5 degrees
    assert(std::abs(Qgen2 - Qcool2 - Qheat2) / std::abs(Qgen2) < 1e-10);

    //!< ******************************************************************************************************************************************************
    //!< make the hierarchical module
    //!< cout<<"Module_s_test start coolsystem test with a complex module for cool control setting "<<coolControl<<endl;
    std::string n11 = "H1";
    std::string n22 = "H2";
    std::string n33 = "H3";
    auto cp11 = std::make_unique<Cell_SPM>();
    auto cp22 = std::make_unique<Cell_SPM>();
    auto cp33 = std::make_unique<Cell_SPM>();
    auto cp44 = std::make_unique<Cell_SPM>();
    auto cp55 = std::make_unique<Cell_SPM>();
    auto cp66 = std::make_unique<Cell_SPM>();
    auto cp77 = std::make_unique<Cell_SPM>();
    std::unique_ptr<StorageUnit> SU1[] = { cp11, cp22 };
    std::unique_ptr<StorageUnit> SU2[] = { cp33, cp44 };
    std::unique_ptr<StorageUnit> SU3[] = { cp55, cp66, cp77 };
    auto mp11 = std::make_unique<Module_s>(n11, T, true, false, std::size(SU1), coolControl, 2); //!< pass through cool system
    auto mp22 = std::make_unique<Module_s>(n22, T, true, false, std::size(SU2), coolControl, 2);
    auto mp33 = std::make_unique<Module_s>(n33, T, true, false, std::size(SU3), coolControl, 2);
    mp11->setSUs(SU1, checkCells);
    mp22->setSUs(SU2, checkCells);
    mp33->setSUs(SU3, checkCells);
    int nm = 3;
    std::string n44 = "H4";
    std::unique_ptr<StorageUnit> MU[] = { mp11, mp22, mp33 };
    auto mp44 = std::make_unique<Module_s>(n44, T, true, true, 7, coolControl, 1);
    mp44->setSUs(MU, checkCells, true);
    double Tini22[7] = { cp11->T(), cp22->T(), cp33->T(), cp44->T(), cp55->T(), cp66->T(), cp77->T() };

    //!< do a few 1C cycles (note just some time steps since we don't have the Cycler
    Icha = -cp11->Cap();
    Idis = cp11->Cap();
    for (int i = 0; i < 5; i++) {
      //!< charge
      mp44->setCurrent(Icha);
      for (int t = 0; t < N; t++)
        mp44->timeStep_CC(dt);
      //!< rest
      mp44->setCurrent(Irest);
      for (int t = 0; t < N; t++)
        mp44->timeStep_CC(dt);
      //!< discharge
      mp44->setCurrent(Idis);
      for (int t = 0; t < N; t++)
        mp44->timeStep_CC(dt);
      //!< rest
      mp44->setCurrent(Irest);
      for (int t = 0; t < N; t++)
        mp44->timeStep_CC(dt);
    }

    double Qgen3, Qcool3, Qheat3;
    //!< check balance of module mp11
    Qgen3 = cp11->thermal_getTotalHeat() + cp22->thermal_getTotalHeat();                                                     //!< total heat generated by cells
    Qcool3 = mp11->getCoolSystem()->getHeatEvac();                                                                           //!< total heat extracted by the coolsystem from the cells
    Qheat3 = -((Tini22[0] - cp11->T()) * (rho * Cp * L * elec_surf) + (Tini22[1] - cp22->T()) * (rho * Cp * L * elec_surf)); //!< total energy in heating up the cells
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) / std::abs(Qgen3) < 1e-10);
    //!< check balance of module mp22
    Qgen3 = cp33->thermal_getTotalHeat() + cp44->thermal_getTotalHeat();                                                     //!< total heat generated by cells
    Qcool3 = mp22->getCoolSystem()->getHeatEvac();                                                                           //!< total heat extracted by the coolsystem from the cells
    Qheat3 = -((Tini22[2] - cp33->T()) * (rho * Cp * L * elec_surf) + (Tini22[3] - cp44->T()) * (rho * Cp * L * elec_surf)); //!< total energy in heating up the cells
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) / std::abs(Qgen3) < 1e-10);
    //!< check balance of module mp33
    Qgen3 = cp55->thermal_getTotalHeat() + cp66->thermal_getTotalHeat() + cp77->thermal_getTotalHeat();                                                                             //!< total heat generated by cells
    Qcool3 = mp33->getCoolSystem()->getHeatEvac();                                                                                                                                  //!< total heat extracted by the coolsystem from the cells
    Qheat3 = -((Tini22[4] - cp55->T()) * (rho * Cp * L * elec_surf) + (Tini22[5] - cp66->T()) * (rho * Cp * L * elec_surf) + (Tini22[6] - cp77->T()) * (rho * Cp * L * elec_surf)); //!< total energy in heating up the cells
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) / std::abs(Qgen3) < 1e-10);

    //!< check balance of the top level module
    Qgen3 = mp11->getCoolSystem()->getHeatEvac() + mp22->getCoolSystem()->getHeatEvac() + mp33->getCoolSystem()->getHeatEvac();
    Qcool3 = mp44->getCoolSystem()->getHeatEvac();
    Qheat3 = mp11->getCoolSystem()->getHeatabsorbed() + mp22->getCoolSystem()->getHeatabsorbed() + mp33->getCoolSystem()->getHeatabsorbed();
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) / std::abs(Qgen3) < 1e-10);

    //!< check balance of total system
    Qgen3 = cp11->thermal_getTotalHeat() + cp22->thermal_getTotalHeat();
    Qgen3 += cp33->thermal_getTotalHeat() + cp44->thermal_getTotalHeat();
    Qgen3 += cp55->thermal_getTotalHeat() + cp66->thermal_getTotalHeat() + cp77->thermal_getTotalHeat();
    Qcool3 = mp44->getCoolSystem()->getHeatEvac();
    Qheat3 = -((Tini22[0] - cp11->T()) * (rho * Cp * L * elec_surf) + (Tini22[1] - cp22->T()) * (rho * Cp * L * elec_surf));
    Qheat3 += -((Tini22[2] - cp33->T()) * (rho * Cp * L * elec_surf) + (Tini22[3] - cp44->T()) * (rho * Cp * L * elec_surf));
    Qheat3 += -((Tini22[4] - cp55->T()) * (rho * Cp * L * elec_surf) + (Tini22[5] - cp66->T()) * (rho * Cp * L * elec_surf) + (Tini22[6] - cp77->T()) * (rho * Cp * L * elec_surf));
    Qheat3 += mp11->getCoolSystem()->getHeatabsorbed() + mp22->getCoolSystem()->getHeatabsorbed() + mp33->getCoolSystem()->getHeatabsorbed();
    //!< cout<<"\t Total heat balance of coolsystem complex module "<<coolControl<<" is "<<Qgen3<<", "<<Qheat3<<", "<<Qcool3<<" and error "<<abs(Qgen3 - Qcool3 - Qheat3)<<endl;
    assert(std::abs(Qgen3 - Qcool3 - Qheat3) / std::abs(Qgen3) < 1e-10);

    //!< Comparison of cool system performance in the different control strategies: print out the following statement
    //!< cout<<"Total heat balance of coolsystem complex module "<<coolControl<<" is "<<Qgen3<<", "<<Qheat3<<", "<<Qcool3<<" and error "<<abs(Qgen3 - Qcool3 - Qheat3)<<endl;
    //!< note: Qheat < 0 means the cell has cooled down compared to its initial value
    //!< 	that is possible because cells are made at T_env (which is 20), but the HVAC coolsystem can cool the SU as cold as controlAC_onoff_Toff, which is 15
    //!< 	so cells can cool down 5 degrees
  }

  return true;
}

//!< ***************************************************************** test all functions *************************************************************************
int test_all_Module_s()
{
  //!< if we test the errors, suppress error messages
  //!< "Pure" unit tests (series modules with Cells)
  if (!TEST(test_Constructor, "Cell_test")) return 1;
  if (!TEST(test_BasicGetters, "Cell_test")) return 2;
  if (!TEST(test_BasicGetters_s, "Cell_test")) return 3;
  if (!TEST(test_setI, "Cell_test")) return 4;

  if (!TEST(test_getCellV, "Cell_test")) return 5;
  if (!TEST(test_getStates, "Cell_test")) return 6;
  if (!TEST(test_getCells, "Cell_test")) return 7;
  if (!TEST(test_setT, "Cell_test")) return 8;

  if (!TEST(test_setStates, "Cell_test")) return 9;
  if (!TEST(test_validCells, "Cell_test")) return 10; //!< includes setState
  if (!TEST(test_validStates, "Cell_test")) return 11;
  if (!TEST(test_setCells, "Cell_test")) return 12; //!< (includes validCells)

  if (!TEST(test_timeStep_CC, "Cell_test")) return 13;
  if (!TEST(test_copy_s, "Cell_test")) return 14;

  //!< Combinations
  if (!TEST(test_Modules_s_ECM, "Cell_test")) return 15;
  if (!TEST(test_Modules_s_SPM, "Cell_test")) return 16;
  if (!TEST(test_Hierarchichal, "Cell_test")) return 17;      //!< series of series
  if (!TEST(test_Hierarchical_Cross, "Cell_test")) return 18; //!< series of parallel

  //!< coolsystem (includes hierarchical modules and uses SPM cells)
  if (!TEST(test_CoolSystem_s, "Cell_test")) return 19;

  return 0;
}
} // namespace slide::tests::unit

int main() { return slide::tests::unit::test_all_Module_s(); }