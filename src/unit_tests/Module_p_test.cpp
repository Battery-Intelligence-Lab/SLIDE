/*
 * Module_base_p_test.cpp
 *
 *  Created on: 18 Dec 2019
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#include "Module_p.hpp"
#include "Cycler.hpp"

#include "Cell.hpp"
#include "Cell_ECM.hpp"
#include "Cell_SPM.hpp"
#include "unit_tests.hpp"
#include "constants.hpp"
#include "Module_s.hpp"
#include "Interpolation.h"

#include <cassert>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>

namespace slide::unit_tests
{

	void test_Constructor_p()
	{
		// Module_base_s::Module_base_s()
		// Module_base_s::Module_base_s(int ncellsi, Cell ci[], double Ti, bool checkCells, bool print)

		std::unique_ptr<Module_p> mp(new Module_p());
		assert(mp->getNSUs() == 0);
		assert(mp->T() == settings::T_ENV);

		constexpr int ncel = 2;
		std::unique_ptr<Cell> cp1(new Cell);
		std::unique_ptr<Cell> cp2(new Cell);
		std::unique_ptr<StorageUnit> cs[ncel] = {cp1, cp2};
		assert(cp1->getID() == "cell");
		assert(cp1->getFullID() == "cell"); // has no parent yet
		std::string n = "na";
		double T = settings::T_ENV;
		bool checkCells = false;
		std::unique_ptr<Module_p> mp2(new Module_p(n, T, true, false, ncel, 1, 1));
		mp2->setSUs(cs, ncel, checkCells, true);

		assert(mp2->getNSUs() == ncel);
		assert(mp2->T() == T);
	}
	void test_BasicGetters_p()
	{
		// double Module_base_s::Cap()
		// double Module_base_s::Vmin(){
		// double Module_base_s::Vmax()
		// double Module_base_s::I()
		// double double Module_base_s::V()

		constexpr int ncel = 2;
		std::unique_ptr<Cell> cp1(new Cell);
		std::unique_ptr<Cell> cp2(new Cell);
		std::unique_ptr<StorageUnit> cs[ncel] = {cp1, cp2};
		assert(cp1->getID() == "cell");
		assert(cp1->getFullID() == "cell"); // has no parent yet
		std::string n = "na";
		double T = settings::T_ENV;
		bool checkCells = false;
		std::unique_ptr<Module_p> mp(new Module_p(n, T, true, false, ncel, 1, 1));
		mp->setSUs(cs, ncel, checkCells, true);

		assert(mp->Cap() == ncel * cp1->Cap());
		assert(mp->Vmin() == cp1->Vmin());
		assert(mp->Vmax() == cp1->Vmax());
		assert(mp->getVMIN() == cp1->getVMIN());
		assert(mp->getVMAX() == cp1->getVMAX());
		assert(mp->I() == 0);
		assert(mp->V() == cp1->V());
	}
	void test_setI_p(bool testErrors)
	{
		// double Module_base_s::setCurrent(double Inew, bool checkV, bool print)
		double tol = 0.005;
		double Inew, V;

		constexpr int ncel = 2;
		std::unique_ptr<Cell> cp1(new Cell);
		std::unique_ptr<Cell> cp2(new Cell);
		std::unique_ptr<StorageUnit> cs[ncel] = {cp1, cp2};
		assert(cp1->getID() == "cell");
		assert(cp1->getFullID() == "cell"); // has no parent yet
		double v1 = cp1->V();
		std::string n = "na";
		double T = settings::T_ENV;
		bool checkCells = false;
		std::unique_ptr<Module_p> mp(new Module_p(n, T, true, false, ncel, 1, 1));
		mp->setSUs(cs, ncel, checkCells, true);
		assert(mp->I() == 0);
		assert(mp->V() == cp1->V());

		// discharge
		Inew = 1.0 * ncel;
		V = mp->setCurrent(Inew, true);
		assert(std::abs(mp->I() - Inew) < tol);
		assert(V < v1); // voltage must decrease
		// do not check individual cells, that is done in getCells

		// charge
		Inew = -1.0 * ncel;
		V = mp->setCurrent(Inew, true);
		assert(std::abs(mp->I() - Inew) < tol);
		assert(mp->V() > v1); // voltage must increase

		// rest with different SOC values
		Inew = 0;
		cp2->setSOC(0.4);													 // c2 has lower OCV -> should charge
		V = mp->setCurrent(Inew, true);										 // the large change in OCV causes a large voltage change, which cannot be fixed by setCurrent
		assert(std::abs(cp1->V() - cp2->V()) < settings::MODULE_P_V_ABSTOL); // cell voltages are equal
		assert(std::abs(cp1->I() - cp1->I()) < settings::MODULE_P_V_ABSTOL); // cell currents are opposite
		assert(cp1->I() > 0);
		assert(cp2->I() < 0);

		// test things which should break
		if (testErrors)
		{
			Inew = 10000; // very large current, should give too low voltage
			try
			{
				V = mp->setCurrent(Inew); // should fail because the current equation cannot be solved
				assert(false);
			}
			catch (int e)
			{
			};

			Inew = -10000; // very large current, should give too low voltage
			try
			{
				V = mp->setCurrent(Inew); // should fail
				assert(false);
			}
			catch (int e)
			{
			};
		}
	}

	void test_validStates_p(bool testErrors)
	{
		// bool Module_base_s::validStates(double s[], int nin, bool print)

		constexpr int ncel = 2;
		std::unique_ptr<Cell> cp1(new Cell);
		std::unique_ptr<Cell> cp2(new Cell);
		std::unique_ptr<StorageUnit> cs[ncel] = {cp1, cp2};
		assert(cp1->getID() == "cell");
		assert(cp1->getFullID() == "cell"); // has no parent yet
		std::string n = "na";
		double T = settings::T_ENV;
		bool checkCells = false;
		std::unique_ptr<Module_p> mp(new Module_p(n, T, true, false, ncel, 1, 1));
		mp->setSUs(cs, ncel, checkCells, true);

		// valid states (current states)
		int nin = settings::STORAGEUNIT_NSTATES_MAX;
		double s[nin];
		int nout;
		mp->getStates(s, nin, nout);
		assert(mp->validStates(s, nout));

		// valid states (new T)
		s[nout - 1] = 273 + 5;
		assert(mp->validStates(s, nout));

		if (testErrors)
		{

			// wrong length
			int nc = settings::CELL_NSTATE_MAX;
			double sc[nc];
			int noutc;
			cp1->getStates(sc, nc, noutc);
			assert(!mp->validStates(sc, noutc, false));

			// an SOC which is too large
			s[0] = 2; // this is the SOC of cell 1
			s[nout - 1] = 273 + 1;
			assert(!mp->validStates(s, nout, false));

			// different voltage values (by changing the SOCs)
			s[0] = 0.4; // this is the SOC of cell 1 (soc of cell 2 is 0.5)
			assert(!mp->validStates(s, nout, false));
		}
	}
	void test_validCells_p(bool testErrors)
	{
		// bool Module_base_s::validCells(Cell c[], int nin, bool print)

		constexpr int ncel = 2;
		std::unique_ptr<Cell> cp1(new Cell);
		std::unique_ptr<Cell> cp2(new Cell);
		std::unique_ptr<Cell> cp3(new Cell);
		std::unique_ptr<StorageUnit> cs[ncel] = {cp1, cp2};
		assert(cp1->getID() == "cell");
		assert(cp1->getFullID() == "cell"); // has no parent yet
		std::string n = "na";
		double T = settings::T_ENV;
		bool checkCells = false;
		std::unique_ptr<Module_p> mp(new Module_p(n, T, true, false, ncel, 1, 1));
		mp->setSUs(cs, ncel, checkCells, true);

		// valid cells with the current cells
		std::unique_ptr<StorageUnit> cs2[ncel];
		int nout;
		mp->getSUs(cs2, ncel, nout);
		assert(mp->validSUs(cs2, ncel));

		// valid cells (new cell)
		cs2[1] = cp3;
		assert(mp->validSUs(cs2, ncel));

		// add an additional cell
		std::unique_ptr<StorageUnit> cs3[3] = {cp1, cp2, cp3};
		assert(mp->validSUs(cs3, 3));

		if (testErrors)
		{
			// different SOC values -> different voltages
			cp3->setSOC(0.4);
			assert(!mp->validSUs(cs3, 3, false));
		}
	}
	void test_timeStep_CC_p()
	{
		// void Module_base_s::timeStep_CC(double dt)

		double T = settings::T_ENV;
		bool checkCells = false;
		constexpr int ncel = 2;
		std::unique_ptr<Cell> cp1(new Cell);
		std::unique_ptr<Cell> cp2(new Cell);
		std::unique_ptr<StorageUnit> cs[ncel] = {cp1, cp2};
		assert(cp1->getID() == "cell");
		assert(cp1->getFullID() == "cell"); // has no parent yet
		std::string n = "na";
		std::unique_ptr<Module_p> mp(new Module_p(n, T, true, false, ncel, 1, 1));
		mp->setSUs(cs, ncel, checkCells, true);
		double v1 = cp1->V();
		double soc1 = cp1->SOC();

		// time step with 0 current
		double dt = 5;
		mp->timeStep_CC(dt);
		assert(mp->V() == cp1->V());

		// discharge
		double Inew = 1 * ncel;
		double V, err;
		double tol = settings::MODULE_P_I_ABSTOL;
		mp->setCurrent(Inew);
		mp->timeStep_CC(dt);
		V = mp->V();
		assert(V < v1);
		err = mp->I() - Inew;
		assert(err < tol && err > -tol); // note: this is not going to be exact because we solve a nonlinear equation to obain I

		// check individual cells
		std::unique_ptr<StorageUnit> cs2[ncel];
		int nout;
		mp->getSUs(cs2, ncel, nout);
		Cell *cell1;
		for (int i = 0; i < ncel; i++)
		{
			err = cs2[i]->I() - Inew / ncel; // we know the current has to split equally between both cells
			assert(err < tol && err > -tol);
			assert(cs2[i]->V() < v1);
			cell1 = dynamic_cast<Cell *>(cs2[i].get()); // Dynamic cast from StorageUnit to Cell
			assert(cell1->SOC() < soc1);
			err = cell1->SOC() - (0.5 - 1.0 * 5.0 / 3600.0 / cell1->Cap());
			assert(err < tol && err > -tol);
		}

		// charge
		Inew = -1;
		mp->setCurrent(Inew);
		mp->timeStep_CC(dt);
		mp->getSUs(cs2, ncel, nout);
		assert(mp->V() > V);
		err = mp->I() - Inew;
		assert(err < tol && err > -tol);
		for (int i = 0; i < ncel; i++)
		{
			err = cs2[i]->I() - Inew / ncel; // we know the current has to split equally between both cells
			assert(err < tol && err > -tol);
			assert(cs2[i]->V() > V);
			cell1 = dynamic_cast<Cell *>(cs2[i].get()); // Dynamic cast from StorageUnit to Cell
			err = cell1->SOC() - (0.5);
			assert(err < tol && err > -tol);
		}

		// loop to charge and check cell voltages
		cp2->setSOC(0.49); // start off with a slightly different soC
		for (int t = 0; t < 100; t++)
		{
			mp->timeStep_CC(dt);
			assert(mp->validSUs());
		}
	}

	void test_Modules_p_ECM(bool testErrors)
	{
		// test parallel modules with ECM cells
		double tol = settings::MODULE_P_I_ABSTOL;

		// setCurrent
		constexpr int ncel = 2;
		std::unique_ptr<Cell_ECM> cp1(new Cell_ECM);
		std::unique_ptr<Cell_ECM> cp2(new Cell_ECM);
		std::unique_ptr<StorageUnit> cs[ncel] = {cp1, cp2};
		assert(cp1->getID() == "cell");
		assert(cp1->getFullID() == "cell"); // has no parent yet
		std::string n = "na";
		double T = settings::T_ENV;
		bool checkCells = false;
		std::unique_ptr<Module_p> mp(new Module_p(n, T, true, false, ncel, 1, 1));
		mp->setSUs(cs, ncel, checkCells, true);
		double v1 = cp1->V();
		assert(mp->I() == 0);
		assert(mp->V() == cp1->V());
		// discharge
		double Inew = 1.0 * ncel;
		double V;
		V = mp->setCurrent(Inew, true);
		assert(std::abs(mp->I() - Inew) < tol);
		assert(V < v1); // voltage must decrease
		// do not check individual cells, that is done in getCells
		// charge
		Inew = -1.0 * ncel;
		V = mp->setCurrent(Inew, true);
		assert(std::abs(mp->I() - Inew) < tol);
		assert(mp->V() > v1); // voltage must increase
		// rest with different SOC values
		Inew = 0;
		cp2->setSOC(0.4); // c2 has lower OCV -> should charge
		V = mp->setCurrent(Inew, true);
		assert(std::abs(cp1->V() - cp2->V()) < settings::MODULE_P_V_ABSTOL); // cell voltages are equal
		assert(std::abs(cp1->I() + cp2->I()) < tol);						 // cell currents are opposite
		assert(cp1->I() > 0);
		assert(cp2->I() < 0);

		// validCells
		cp1 = std::unique_ptr<Cell_ECM>(new Cell_ECM());
		cp2 = std::unique_ptr<Cell_ECM>(new Cell_ECM());
		std::unique_ptr<Cell_ECM> cp3(new Cell_ECM);
		// valid cells with the current cells
		std::unique_ptr<StorageUnit> cs2[ncel] = {cp1, cp2};
		mp->setSUs(cs2, ncel);
		int nout;
		mp->getSUs(cs2, ncel, nout);
		assert(mp->validSUs(cs2, ncel));
		// valid cells (new cell)
		cs2[1] = cp3;
		assert(mp->validSUs(cs2, ncel));
		// add an additional cell
		std::unique_ptr<StorageUnit> cs3[3] = {cp1, cp2, cp3};
		assert(mp->validSUs(cs3, 3));

		// CC timestep
		cp1 = std::unique_ptr<Cell_ECM>(new Cell_ECM());
		cp2 = std::unique_ptr<Cell_ECM>(new Cell_ECM());
		cs[0] = cp1;
		cs[1] = cp2;
		double soc1 = cp1->SOC();
		v1 = cp1->V();
		mp = std::unique_ptr<Module_p>(new Module_p(n, T, true, false, ncel, 1, 1));
		mp->setSUs(cs, ncel, checkCells);
		// time step with 0 current
		double dt = 5;
		mp->timeStep_CC(dt);
		assert(mp->V() == cp1->V());
		// discharge
		Inew = 1 * ncel;
		double err;
		mp->setCurrent(Inew);
		mp->timeStep_CC(dt);
		V = mp->V();
		assert(V < v1);
		err = mp->I() - Inew;
		assert(err < tol && err > -tol); // note: this is not going to be exact because we solve a nonlinear equation to obain I
		// check individual cells
		mp->getSUs(cs2, ncel, nout);
		Cell_ECM *cell1;
		for (int i = 0; i < ncel; i++)
		{
			err = cs2[i]->I() - Inew / ncel; // we know the current has to split equally between both cells
			assert(err < tol && err > -tol);
			assert(cs2[i]->V() < v1);
			cell1 = dynamic_cast<Cell_ECM *>(cs2[i].get()); // Dynamic cast from StorageUnit to Cell
			assert(cell1->SOC() < soc1);
			err = cell1->SOC() - (0.5 - 1.0 * 5.0 / 3600.0 / cell1->Cap());
			assert(err < tol && err > -tol);
		}
		// charge
		Inew = -1;
		mp->setCurrent(Inew);
		mp->timeStep_CC(dt);
		mp->getSUs(cs2, ncel, nout);
		assert(mp->V() > V);
		err = mp->I() - Inew;
		assert(err < tol && err > -tol);
		for (int i = 0; i < ncel; i++)
		{
			err = cs2[i]->I() - Inew / ncel; // we know the current has to split equally between both cells
			assert(err < tol && err > -tol);
			assert(cs2[i]->V() > V);
			cell1 = dynamic_cast<Cell_ECM *>(cs2[i].get()); // Dynamic cast from StorageUnit to Cell
			err = cell1->SOC() - (0.5);
			assert(err < tol && err > -tol);
		}
	}

	void test_Modules_p_SPM(bool testErrors)
	{
		// test parallel modules with ECM cells
		double tol = settings::MODULE_P_I_ABSTOL;
		Cell_SPM *cell1;

		// setCurrent
		constexpr int ncel = 2;
		std::unique_ptr<Cell> cp1(new Cell);
		std::unique_ptr<Cell> cp2(new Cell);
		std::unique_ptr<StorageUnit> cs[ncel] = {cp1, cp2};
		assert(cp1->getID() == "cell");
		assert(cp1->getFullID() == "cell"); // has no parent yet
		std::string n = "na";
		double T = settings::T_ENV;
		bool checkCells = false;
		std::unique_ptr<Module_p> mp(new Module_p(n, T, true, false, ncel, 1, 1));
		mp->setSUs(cs, ncel, checkCells, true);
		double v1 = cp1->V();
		assert(mp->I() == 0);
		assert(mp->V() == cp1->V());
		// discharge
		double Inew = 1.0 * ncel;
		double V;
		V = mp->setCurrent(Inew, true);
		assert(std::abs(mp->I() - Inew) < tol);
		assert(V < v1); // voltage must decrease
		// do not check individual cells, that is done in getCells
		// charge
		Inew = -1.0 * ncel;
		V = mp->setCurrent(Inew, true);
		assert(std::abs(mp->I() - Inew) < tol);
		assert(mp->V() > v1); // voltage must increase
		// rest with different SOC values
		// -> must be different lithium fractions since SOC does not affect the concentration
		// too complicated so skip this test here

		// validCells
		cp1 = std::unique_ptr<Cell_SPM>(new Cell_SPM());
		cp2 = std::unique_ptr<Cell_SPM>(new Cell_SPM());
		std::unique_ptr<Cell_SPM> cp3(new Cell_SPM);
		std::unique_ptr<StorageUnit> cs2[ncel] = {cp1, cp2};
		mp->setSUs(cs2, ncel);
		int nout;
		// valid cells with the current cells
		std::unique_ptr<StorageUnit> cs22[ncel]; // if this is cs2 and you reuse it in CC time step, for some reason, the code fails
		mp->getSUs(cs22, ncel, nout);			 // similarly, if you comment out the creation of cs22, the code also fails
		assert(mp->validSUs(cs22, ncel));		 // that seems to suggest there is a data-overflow somewhere
		// valid cells (new cell)
		cs22[1] = cp3;
		assert(mp->validSUs(cs22, ncel));
		// add an additional cell
		std::unique_ptr<StorageUnit> cs33[3] = {cp1, cp2, cp3};
		assert(mp->validSUs(cs33, 3));

		// CC timestep
		cp1 = std::unique_ptr<Cell_SPM>(new Cell_SPM());
		cp2 = std::unique_ptr<Cell_SPM>(new Cell_SPM());
		cs[0] = cp1;
		cs[1] = cp2;
		double soc1 = cp1->SOC();
		v1 = cp1->V();
		mp->setSUs(cs, ncel, checkCells, true);
		// time step with 0 current
		double dt = 5;
		mp->timeStep_CC(dt);
		assert(mp->V() == cp1->V());
		// discharge
		Inew = 1 * ncel;
		double err;
		mp->setCurrent(Inew);
		mp->timeStep_CC(dt);
		V = mp->V();
		assert(V < v1);
		err = mp->I() - Inew;
		assert(err < tol && err > -tol); // note: this is not going to be exact because we solve a nonlinear equation to obain I
		// check individual cells
		mp->getSUs(cs2, ncel, nout);
		for (int i = 0; i < ncel; i++)
		{
			err = cs2[i]->I() - Inew / ncel; // we know the current has to split equally between both cells
			assert(err < tol && err > -tol);
			assert(cs2[i]->V() < v1);
			cell1 = dynamic_cast<Cell_SPM *>(cs2[i].get()); // Dynamic cast from StorageUnit to Cell
			err = cell1->SOC() - (0.5 - 1.0 * 5.0 / 3600.0 / cell1->Cap());
			assert(cell1->SOC() < soc1);
			assert(err < tol && err > -tol);
		}
		// charge
		Inew = -1;
		mp->setCurrent(Inew);
		mp->timeStep_CC(dt);
		mp->getSUs(cs2, ncel, nout);
		assert(mp->V() > V);
		err = mp->I() - Inew;
		assert(err < tol && err > -tol);
		for (int i = 0; i < ncel; i++)
		{
			err = cs2[i]->I() - Inew / ncel; // we know the current has to split equally between both cells
			assert(err < tol && err > -tol);
			assert(cs2[i]->V() > V);
			cell1 = dynamic_cast<Cell_SPM *>(cs2[i].get()); // Dynamic cast from StorageUnit to Cell
			err = cell1->SOC() - (0.5);
			assert(err < tol && err > -tol);
		}
	}

	void test_contactR()
	{
		/*
		 * Make a module with 3 cells and a contact resistance
		 */

		double Rc = 0.01;
		double tol = 0.0001;

		int ncel = 3;
		std::unique_ptr<Cell> cp1(new Cell);
		std::unique_ptr<Cell> cp2(new Cell);
		std::unique_ptr<Cell> cp3(new Cell);
		std::unique_ptr<StorageUnit> cs[ncel] = {cp1, cp2, cp3};
		double Rcs[ncel] = {Rc, Rc, Rc};
		std::string n = "na";
		double T = settings::T_ENV;
		bool checkCells = false;
		std::unique_ptr<Module_p> mp(new Module_p(n, T, true, false, ncel, 1, 1));
		mp->setSUs(cs, ncel, checkCells, true);
		mp->setRcontact(Rcs, ncel);

		// total resistance:
		// 		-Rp+-Rp-+-Rp-|
		// 		   |    |    |
		// 		  Rs    Rs   Rs
		// where Rp = contact resistance (value Rc) and Rs = cell resistance = 0.01;
		double Rs = cp1->getRtot();
		double Rq = Rc + (Rs * (Rs + Rc)) / (Rc + 2 * Rs); // resistance of last two branches
		double Rtot = Rc + Rs * Rq / (Rs + Rq);
		assert(std::abs(mp->getRtot() - Rtot) < tol);

		// setCurrent
		// check voltages at each node from the branch going 'down' and the branch going 'right'
		// 	R1*I1 = Rp2*(I2 + I3) + R2*I2
		// 	R2*I2 = (Rp3 + R3)*I3
		// 	where Ri = resistance of cell i
		// 		  Rpi = contact resistance in parallel at cell i
		// 	since all cells have the same OCV
		double I = 20;
		mp->setCurrent(I, true, true);
		double I1 = cp1->I();
		double I2 = cp2->I();
		double I3 = cp3->I();

		// assert the currents of cells further from the terminal are smaller
		assert(std::abs(I1) > std::abs(I2));
		assert(std::abs(I2) > std::abs(I3));

		// Check the voltage at every node
		double Rcell = cp1->getRtot();				  // all cell resistances are the same
		double V11 = Rcell * I1;					  // voltage at the node connecting the first cell, going down [ignoring OCV]
		double V12 = Rcs[1] * (I2 + I3) + Rcell * I2; // voltage at the node connecting the first cell, going right
		// double V13 = Rcs[1] * (I2 + I3) + Rcs[2]*I3 + Rcell*I3;
		double V22 = Rcell * I2;			// voltage at node of 2nd cell going down
		double V23 = (Rcs[2] + Rcell) * I3; // voltage at node of 2nc cell going right
		assert(std::abs(V11 - V12) < settings::MODULE_P_V_ABSTOL);
		assert(std::abs(V22 - V23) < settings::MODULE_P_V_ABSTOL);
		assert(mp->validSUs()); // ensure the SUs have valid voltages

		// check the total voltage
		double V1 = cp1->V() - Rcs[0] * (I1 + I2 + I3);
		double V2 = cp2->V() - Rcs[0] * (I1 + I2 + I3) - Rcs[1] * (I2 + I3);
		double V3 = cp3->V() - Rcs[0] * (I1 + I2 + I3) - Rcs[1] * (I2 + I3) - Rcs[2] * I3;
		assert(std::abs(V1 - V2) < settings::MODULE_P_V_ABSTOL);
		assert(std::abs(V1 - V3) < settings::MODULE_P_V_ABSTOL);
		assert(std::abs(V2 - V3) < settings::MODULE_P_V_ABSTOL);
		assert(std::abs(mp->V() - V2) < settings::MODULE_P_V_ABSTOL);

		assert(std::abs(V1 - mp->Vi(0)) < tol); // these numbers should be exactly the same
		assert(std::abs(V2 - mp->Vi(1)) < tol); // these numbers should be exactly the same
		assert(std::abs(V3 - mp->Vi(2)) < tol); // these numbers should be exactly the same

		// set charging current
		I = -20;
		mp->setCurrent(I, true, true);
		I1 = cp1->I();
		I2 = cp2->I();
		I3 = cp3->I();
		assert(std::abs(I1) > std::abs(I2));
		assert(std::abs(I2) > std::abs(I3));
		Rcell = cp1->getRtot();				   // all cell resistances are the same
		V11 = Rcell * I1;					   // voltage at the node connecting the first cell, going down [ignoring OCV]
		V12 = Rcs[1] * (I2 + I3) + Rcell * I2; // voltage at the node connecting the first cell, going right
		V22 = Rcell * I2;					   // voltage at node of 2nd cell going down
		V23 = (Rcs[2] + Rcell) * I3;		   // voltage at node of 2nc cell going right
		assert(std::abs(V11 - V12) < settings::MODULE_P_V_ABSTOL);
		assert(std::abs(V22 - V23) < settings::MODULE_P_V_ABSTOL);
		assert(mp->validSUs()); // ensure the SUs have valid voltages
		V1 = cp1->V() - Rcs[0] * (I1 + I2 + I3);
		V2 = cp2->V() - Rcs[0] * (I1 + I2 + I3) - Rcs[1] * (I2 + I3);
		V3 = cp3->V() - Rcs[0] * (I1 + I2 + I3) - Rcs[1] * (I2 + I3) - Rcs[2] * I3;
		assert(std::abs(V1 - V2) < settings::MODULE_P_V_ABSTOL);
		assert(std::abs(V1 - V3) < settings::MODULE_P_V_ABSTOL);
		assert(std::abs(V2 - V3) < settings::MODULE_P_V_ABSTOL);
		assert(std::abs(mp->V() - V2) < settings::MODULE_P_V_ABSTOL);

		assert(std::abs(V1 - mp->Vi(0)) < tol); // these numbers should be exactly the same
		assert(std::abs(V2 - mp->Vi(1)) < tol); // these numbers should be exactly the same
		assert(std::abs(V3 - mp->Vi(2)) < tol); // these numbers should be exactly the same
	}

	void test_Hierarchichal_p()
	{
		// test parallel modules made out of other parallel modules

		double tol = settings::MODULE_P_I_ABSTOL;
		constexpr int ncel1 = 2;
		constexpr int ncel2 = 2;
		int ncel3 = 3;
		std::string n1 = "H1";
		std::string n2 = "H2";
		std::string n3 = "H3";
		std::unique_ptr<Cell> cp1(new Cell);
		std::unique_ptr<Cell> cp2(new Cell);
		std::unique_ptr<Cell> cp3(new Cell);
		std::unique_ptr<Cell> cp4(new Cell);
		std::unique_ptr<Cell> cp5(new Cell);
		std::unique_ptr<Cell> cp6(new Cell);
		std::unique_ptr<Cell> cp7(new Cell);
		std::unique_ptr<StorageUnit> SU1[ncel1] = {cp1, cp2};
		std::unique_ptr<StorageUnit> SU2[ncel2] = {cp3, cp4};
		std::unique_ptr<StorageUnit> SU3[ncel3] = {cp5, cp6, cp7};
		double T = settings::T_ENV;
		bool checkCells = false;
		std::unique_ptr<Module_p> mp1(new Module_p(n1, T, true, false, ncel1, 1, 2)); // pass through cool systems
		std::unique_ptr<Module_p> mp2(new Module_p(n2, T, true, false, ncel2, 1, 2));
		std::unique_ptr<Module_p> mp3(new Module_p(n3, T, true, false, ncel3, 1, 2));
		mp1->setSUs(SU1, ncel1, checkCells);
		mp2->setSUs(SU2, ncel2, checkCells);
		mp3->setSUs(SU3, ncel3, checkCells);

		// make the hierarichical module
		int nm = 3;
		std::string n4 = "4";
		std::unique_ptr<StorageUnit> MU[nm] = {mp1, mp2, mp3};
		checkCells = true;
		std::unique_ptr<Module_p> mp(new Module_p(n4, T, true, true, 7, 1, 1));
		mp->setSUs(MU, nm, checkCells, true);
		double Vini = mp->V();
		assert(std::abs(Vini - mp1->V()) < tol);
		assert(std::abs(Vini - cp5->V()) < tol);
		assert(mp->getFullID() == "4");
		assert(mp1->getFullID() == "4_H1");
		assert(cp1->getFullID() == "4_H1_cell");
		assert(cp4->getFullID() == "4_H2_cell");
		assert(cp5->getFullID() == "4_H3_cell");

		// set a CC current
		double Inew = -14; // should give about 2A per cell
		mp->setCurrent(Inew);
		assert(std::abs(mp->I() - Inew) < tol);
		assert(std::abs(mp1->V() - mp2->V()) < settings::MODULE_P_V_ABSTOL);
		assert(std::abs(mp3->V() - mp2->V()) < settings::MODULE_P_V_ABSTOL);
		/*assert( std::abs(cp3->I() + 2)<tol); // this is not necessarily true, the functions only put constraints on the voltage
		assert( std::abs(cp7->I() + 2)<tol);		// so we should not expect certain tolerances on the current. E.g. if R very small, we can get larger current errors
		assert( std::abs(mp1->I() + 2*2)<tol); // m1 has two cells
		assert( std::abs(mp3->I() + 2*3)<tol); // m3 has three cells*/

		// time a CC time step
		Vini = mp->V();
		double dt = 5;
		mp->timeStep_CC(dt);
		assert(std::abs(cp1->SOC() - (0.5 - 2 * dt / 3600.0 / cp1->Cap())) < tol); // the SOC must have increased (check just 1 cell out of all 7)
		assert(mp->V() > Vini);
		assert(std::abs(mp2->V() - mp3->V()) < tol); // submodules must have same voltage
	}

	void test_Hierarchical_cross_p()
	{
		// test parallel module made out of series modules
		// note: series modules must have same number of cells to get the same voltage
		double tol = settings::MODULE_P_I_ABSTOL;
		constexpr int ncel1 = 2;
		constexpr int ncel2 = 2;
		constexpr int ncel3 = 2;
		std::string n1 = "H1";
		std::string n2 = "H2";
		std::string n3 = "H3";
		std::unique_ptr<Cell> cp1(new Cell);
		std::unique_ptr<Cell> cp2(new Cell);
		std::unique_ptr<Cell> cp3(new Cell);
		std::unique_ptr<Cell> cp4(new Cell);
		std::unique_ptr<Cell> cp5(new Cell);
		std::unique_ptr<Cell> cp6(new Cell);
		std::unique_ptr<StorageUnit> SU1[ncel1] = {cp1, cp2};
		std::unique_ptr<StorageUnit> SU2[ncel2] = {cp3, cp4};
		std::unique_ptr<StorageUnit> SU3[ncel3] = {cp5, cp6};
		double T = settings::T_ENV;
		bool checkCells = false;
		std::unique_ptr<Module_s> mp1(new Module_s(n1, T, true, false, ncel1, 1, 2)); // pass through coolsystems
		std::unique_ptr<Module_s> mp2(new Module_s(n2, T, true, false, ncel2, 1, 2));
		std::unique_ptr<Module_s> mp3(new Module_s(n3, T, true, false, ncel3, 1, 2));
		mp1->setSUs(SU1, ncel1, checkCells);
		mp2->setSUs(SU2, ncel2, checkCells);
		mp3->setSUs(SU3, ncel3, checkCells);

		// make the hierarichical module
		int nm = 3;
		std::string n4 = "4";
		std::unique_ptr<StorageUnit> MU[nm] = {mp1, mp2, mp3};
		checkCells = true;
		std::unique_ptr<Module_p> mp(new Module_p(n4, T, true, true, 7, 1, 1));
		mp->setSUs(MU, nm, checkCells, true);
		double Vini = mp->V();
		assert(std::abs(Vini - mp1->V()) < tol);
		assert(std::abs(Vini - cp5->V() * 2) < tol); // one module has 2 cells so voltage should split in 2

		// set a CC current
		double Inew = -6; // should give about 2A per cell
		mp->setCurrent(Inew);
		assert(std::abs(mp->I() - Inew) < tol);
		assert(std::abs(cp3->I() + 2) < tol);
		assert(std::abs(cp6->I() + 2) < tol);
		assert(std::abs(mp1->I() + 2) < tol);		 // m1 has two cells
		assert(std::abs(mp3->I() + 2) < tol);		 // m3 has two cells
		assert(std::abs(mp1->V() - mp3->V()) < tol); // check voltage is equal

		/* the iterative function is now private
			Inew = 5;
			mp->setI_iterative(Inew);
			assert(mp->I() == Inew);
			Inew = -5;
			mp->setI_iterative(Inew);
			assert(mp->I() == Inew);
			Inew = -6; // should give about 2A per cell
			mp->setCurrent(Inew);
		*/

		// time a CC time step
		Vini = mp->V();
		double dt = 5;
		mp->timeStep_CC(dt);
		assert(std::abs(cp1->SOC() - (0.5 - 2 * dt / 3600.0 / cp1->Cap())) < tol); // the SOC must have increased (check just 1 cell out of all 7)
		assert(mp->V() > Vini);
		assert(std::abs(mp2->V() - mp3->V()) < tol); // submodules must have same voltage
													 // note: there is no check on sub-modules with different SOC but I assume that works since it works with sub-cells of different SOC
	}

	// void test_copy_p()
	// {
	// 	/*
	// 	 * test the copy-function
	// 	 */

	// 	// make module
	// 	constexpr int ncel = 2;
	// 	std::unique_ptr<Cell> cp1(new Cell);
	// 	std::unique_ptr<Cell> cp2(new Cell);
	// 	std::unique_ptr<StorageUnit> cs[ncel] = {cp1, cp2};
	// 	std::string n = "na";
	// 	double v1 = cp1->V();
	// 	double T = settings::T_ENV;
	// 	bool checkCells = false;
	// 	std::unique_ptr<Module_p> mp(new Module_p(n, T, true, false, ncel, 1, 1));
	// 	mp->setSUs(cs, ncel, checkCells, true);

	// 	// copy this one and check they are identical
	// 	std::unique_ptr<StorageUnit> cn = mp->copy();
	// 	Module_p *c22 = dynamic_cast<Module_p *>(cn.get()); // Dynamic cast from StorageUnit to Cell
	// 	assert(mp->V() == c22->V());
	// 	std::unique_ptr<StorageUnit> corig[ncel], cnew[ncel];
	// 	int nout;
	// 	mp->getSUs(corig, ncel, nout);
	// 	c22->getSUs(cnew, ncel, nout);
	// 	for (int i = 0; i < mp->getNSUs(); i++)
	// 		assert(corig[i]->V() == cnew[i]->V());

	// 	// change the copied version, and ensure the old one is still the same
	// 	c22->setCurrent(1 * ncel, false, false); // discharge
	// 	for (int t = 0; t < 10; t++)
	// 		c22->timeStep_CC(2, false);
	// 	mp->getSUs(corig, ncel, nout);
	// 	c22->getSUs(cnew, ncel, nout);
	// 	for (int i = 0; i < mp->getNSUs(); i++)
	// 	{
	// 		assert(corig[i]->V() == v1);
	// 		assert(cnew[i]->V() < v1);
	// 	}
	// }

	void test_equaliseV_timing(std::unique_ptr<Module_p> mp, std::unique_ptr<StorageUnit> c[], int nin)
	{
		// test timing
		// IN
		// mp		parallel module
		// SUs 		array with smart pointers to the children of mp

		const auto nout = mp->getNSUs();
		mp->setBlockDegAndTherm(true); // ignore thermal and degradation during this function (we mess with individual cells time time keeping for thermal gives errors)

		// set a 1C current to the individual cells, then redistribute
		double I = mp->Cap();
		double Ii = I / mp->getNcells(); // current per cell
		for (int i = 0; i < nin; i++)
			c[i]->setCurrent(Ii * c[i]->getNcells());

		std::cout << "after setCurrent, the cell voltages are ";
		for (int i = 0; i < nout; i++)
			std::cout << mp->Vi(i) << " ";
		int n1 = mp->redistributeI(true);
		std::cout << "while after " << n1 << " steps in redistributeI, the cell voltages are ";
		for (int i = 0; i < nout; i++)
			std::cout << mp->Vi(i) << " ";
		std::cout << '\n';

		// take 10 time steps, then redistribute again
		int N = 10;
		double dt = 2;
		for (int i = 0; i < nin; i++)
			mp->timeStep_CC_i(i, dt, false, N);

		std::cout << "after 10 time steps, the cell voltages are ";
		for (int i = 0; i < nout; i++)
			std::cout << mp->Vi(i) << " ";

		int n2 = mp->redistributeI(true);

		std::cout << "while after " << n2 << " steps in redistributeI, the cell voltages are ";
		for (int i = 0; i < nout; i++)
			std::cout << mp->Vi(i) << " ";
		std::cout << '\n';

		// go to low SOC (lowest cell voltage is 3V

		auto getSU_Vmin = [&mp, N = nout]()
		{
			double Vmin{1e19};
			for (size_t i = 0; i < N; i++)
				Vmin = std::min(Vmin, mp->Vi(i));

			return Vmin;
		};

		while (getSU_Vmin() > 3.0)
			mp->timeStep_CC(dt, false, 1);

		// take 10 time steps, then redistribute again
		N = 10;
		for (int i = 0; i < nin; i++)
			mp->timeStep_CC_i(i, dt, false, N);

		std::cout << "after discharge, the cell voltages are ";
		for (int i = 0; i < nout; i++)
			std::cout << mp->Vi(i) << " ";

		int n3 = mp->redistributeI(true);
		std::cout << "while after " << n3 << " steps in redistributeI, the cell voltages are ";
		for (int i = 0; i < nout; i++)
			std::cout << mp->Vi(i) << " ";
		std::cout << '\n';

		std::cout << "Number of time steps for mp " << mp->getFullID()
				  << " is for setCurrent " << n1 << ", for timestep " << n2
				  << ", for timestep at low SOC " << n3 << ", now starting a CC cycle\n";

		// Use a Cycler to do a full CC cycle with redistributeI every time step
		Cycler cyc;
		double lim = 0.0;
		int ndata = 0;
		double Ah, Wh, vlim, tlim;
		tlim = 99999999;
		cyc.initialise(mp, mp->getFullID());
		vlim = mp->Vmax() - lim;
		cyc.CC(-I, vlim, tlim, dt, ndata, Ah, Wh); // CC charge
		vlim = mp->Vmin() + lim;
		cyc.CC(I, vlim, tlim, dt, ndata, Ah, Wh); // CC discharge

		std::cout << "Finished CC cycle\n";
	}

	void test_equaliseV()
	{
		// test timing with
		// 		5 identical cells
		// 		5 cells with minor differences
		// 		5 widely different cells
		// 		4 similar and one aged cell

		DEG_ID deg;
		deg.SEI_n = 1;		  // there is 1 SEI model
		deg.SEI_id[0] = 4;	  // chirstensen SEI growth
		deg.SEI_porosity = 0; // don't decrease the porosity (set to 1 if you do want to decrease the porosity)
		deg.CS_n = 1;		  // there is 1 model (that there are no cracks)
		deg.CS_id[0] = 0;	  // no surface cracks
		deg.CS_diffusion = 0; // don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)
		deg.LAM_n = 1;		  // there are 1 LAM model
		deg.LAM_id[0] = 0;	  // no LAM
		deg.pl_id = 0;		  // no litihium plating
		double T2 = settings::T_ENV;
		bool checkCells2 = false;

		// 5 identical cells
		int ncel1 = 5;
		std::string n1 = "mp_identical";
		std::unique_ptr<Cell_SPM> cp7(new Cell_SPM("cell1", deg, 1, 1, 1, 1));
		std::unique_ptr<Cell_SPM> cp8(new Cell_SPM("cell2", deg, 1, 1, 1, 1));
		std::unique_ptr<Cell_SPM> cp9(new Cell_SPM("cell3", deg, 1, 1, 1, 1));
		std::unique_ptr<Cell_SPM> cp10(new Cell_SPM("cell4", deg, 1, 1, 1, 1));
		std::unique_ptr<Cell_SPM> cp11(new Cell_SPM("cell5", deg, 1, 1, 1, 1));
		std::unique_ptr<StorageUnit> cs1[ncel1] = {cp7, cp8, cp9, cp10, cp11};
		std::unique_ptr<Module_p> mpp1(new Module_p(n1, T2, true, false, ncel1, 1, 1)); // no multithreading, nt_Vcheck time steps between checking SU voltage
		mpp1->setSUs(cs1, ncel1, checkCells2, true);
		test_equaliseV_timing(mpp1, cs1, ncel1);

		// 5 cells with small distribution
		default_random_engine gen;
		double std1 = 0;
		double std2 = 0;
		std1 = 0.004;
		std2 = 0.025;
		normal_distribution<double> distr_c(1.0, std1); // normal distribution with mean 1 and std 0.4% for cell capacity
		normal_distribution<double> distr_r(1.0, std2); // normal distribution with mean 1 and std 2.5% for cell resistance
		normal_distribution<double> distr_d(1.0, std2); // normal distribution with mean 1 and std 2.5% for cell degradation rate
		std::string n2 = "mp_variation";
		std::unique_ptr<Cell_SPM> cp16(new Cell_SPM("cell5", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
		std::unique_ptr<Cell_SPM> cp12(new Cell_SPM("cell6", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
		std::unique_ptr<Cell_SPM> cp13(new Cell_SPM("cell7", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
		std::unique_ptr<Cell_SPM> cp14(new Cell_SPM("cell8", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
		std::unique_ptr<Cell_SPM> cp15(new Cell_SPM("cell9", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
		std::unique_ptr<StorageUnit> cs2[ncel1] = {cp16, cp12, cp13, cp14, cp15};
		std::unique_ptr<Module_p> mpp2(new Module_p(n2, T2, true, false, ncel1, 1, 1)); // no multithreading, nt_Vcheck time steps between checking SU voltage
		mpp2->setSUs(cs2, ncel1, checkCells2, true);
		test_equaliseV_timing(mpp2, cs2, ncel1);

		// 5 cells with large distribution
		std1 = 0.1;
		std2 = 0.15;
		normal_distribution<double> distr_c2(1.0, std1); // normal distribution with mean 1 and std 0.4% for cell capacity
		normal_distribution<double> distr_r2(1.0, std2); // normal distribution with mean 1 and std 2.5% for cell resistance
		normal_distribution<double> distr_d2(1.0, std2); // normal distribution with mean 1 and std 2.5% for cell degradation rate
		std::string n3 = "mp_largeVariation";
		std::unique_ptr<Cell_SPM> cp116(new Cell_SPM("cell5", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen)));
		std::unique_ptr<Cell_SPM> cp112(new Cell_SPM("cell6", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen)));
		std::unique_ptr<Cell_SPM> cp113(new Cell_SPM("cell7", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen)));
		std::unique_ptr<Cell_SPM> cp114(new Cell_SPM("cell8", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen)));
		std::unique_ptr<Cell_SPM> cp115(new Cell_SPM("cell9", deg, distr_c2(gen), distr_r2(gen), distr_d2(gen), distr_d2(gen)));
		std::unique_ptr<StorageUnit> cs3[ncel1] = {cp116, cp112, cp113, cp114, cp115};
		std::unique_ptr<Module_p> mpp3(new Module_p(n3, T2, true, false, ncel1, 1, 1)); // no multithreading, nt_Vcheck time steps between checking SU voltage
		mpp3->setSUs(cs3, ncel1, checkCells2, true);
		test_equaliseV_timing(mpp3, cs3, ncel1);

		// 4 similar and one very different
		std::string n4 = "mp_4and1";
		std::unique_ptr<Cell_SPM> cp1116(new Cell_SPM("cell5", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
		std::unique_ptr<Cell_SPM> cp1112(new Cell_SPM("cell6", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
		std::unique_ptr<Cell_SPM> cp1113(new Cell_SPM("cell7", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
		std::unique_ptr<Cell_SPM> cp1114(new Cell_SPM("cell8", deg, distr_c(gen), distr_r(gen), distr_d(gen), distr_d(gen)));
		std::unique_ptr<Cell_SPM> cp1115(new Cell_SPM("cell9", deg, 0.5, 2.0, 1.1, 1.1)); // one with half the capacity and double the resistance
		std::unique_ptr<StorageUnit> cs4[ncel1] = {cp1116, cp1112, cp1113, cp1114, cp1115};
		std::unique_ptr<Module_p> mpp4(new Module_p(n4, T2, true, false, ncel1, 1, 1)); // no multithreading, nt_Vcheck time steps between checking SU voltage
		mpp4->setSUs(cs4, ncel1, checkCells2, true);
		test_equaliseV_timing(mpp4, cs4, ncel1);
	}

	void test_Module_p_testAll(bool testErrors)
	{
		/*
		 * Test the functions from the parallel module
		 * note that we already test the function from the base module in the unit test for the series-connected module so there is no point to repeat them
		 */

		// if we test the errors, suppress error messages

		// "pure" unit tests
		test_Constructor_p();
		test_BasicGetters_p();

		test_setI_p(testErrors);
		test_validCells_p(testErrors);
		test_validStates_p(testErrors);

		test_timeStep_CC_p();
		test_contactR();

		// Combinations
		test_Modules_p_ECM(testErrors); // parallel from ECM cells
		test_Modules_p_SPM(testErrors); // parallel from ECM cells
		test_Hierarchichal_p();			// parallel from parallel
		test_Hierarchical_cross_p();	// parallel from series

		test_copy_p();

		// functions to equalise the voltage
		// test_equaliseV(); 			// this prints how many iterations are needed in various situations but does not check wheather things work correctly
	}
}