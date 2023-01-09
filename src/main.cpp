/*
 * main.cpp
 *
 * This file implements the main-function.
 * Here, the user has to choose what to simulate by commenting things in (by placing \\ at the start of the line) and out (by removing the \\)
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

//!< Include header files
#include "slide.hpp"
#include "../benchmark/benchmarks.hpp"
#include "../examples/examples.hpp"

#include <ctime>
#include <thread>
#include <array>
#include <tuple>
#include <random>
#include <cmath>
#include <iomanip>

int main()
{
  /*
   * This function decides what will be simulated.
   * In the code below, you have to uncomment the line which you want to execute (to uncomment, remove the //!<in front of the line)
   * and comment all the other lines in the blocks with FUNCTION CALLS (to comment in, add //!<in front of the lines)
   */

  //!< print that you start simulations
  //!< slide::tests::unit::test_all();
  std::cout << "Start simulations" << std::endl;

  //!< Make a clock to measure how long the simulation takes
  slide::Clock clk;
  //!< Read the values for the matrices for the spatial discretisation of the solid diffusion PDE.
  //!< these values are calculated by the MATLAB-script 'modelSetup.m', which writes them to csv files.
  //!< Here, we only need to read these csv files

  //!< slide::Model_SPM M; //!<Initialised at constructor.
  //!<  Look for extra settings at Constants.hpp (such as verbosity)

  //!< Choose the prefix, a string with which the names of all subfolders with the outputfiles will start (the folders will be called prefix_xxxx).
  //!< This is typically a number, and the reason why we use it is that results of previous simulations are not overwritten.
  //!< i.e. for one simulation, use pref "0", then the next simulation could have pref "1".
  //!< this ensures the results of the second simulation are written in different files and folders, rather than overwriting the files from the first simulation
  const std::string pref = "0"; //!< todo
                                //!<  will be used to name a folder with the results, so must comply with the naming convention on your operating system
                                //!<  letters, numbers, underscores and dots are fine
                                //!<  don't use spaces
                                //!<  don't use forbidden characters, e.g. for Windows the following are not allowed < > : " \ | / ? *
                                //!<  https://docs.microsoft.com/en-us/windows/desktop/fileio/naming-a-file

  //!< Choose which cell should be used for the simulations
  //!< const int cellType = slide::cellType::KokamNMC;
  //!< which cell to use for the simulation.
  //!< 0 	high power Kokam NMC cell (18650)
  //!< 1 	high energy LG Chem NMC cell (18650)
  //!< 2 	user cell (template class provided for where the user can define his own parameters)

  //!< Choose which degradation models to use for the simulations
  //!< The models are identified by a number, see DEG_ID in DEG_ID.hpp and below
  /* SEI_id 		identification for which model(s) to use for SEI growth
   * 					0	no SEI growth
   * 					1 	kinetic model only (Tafel kinetics)
   * 							ref: Ning & Popov, Journal of the Electrochemical Society 151 (10), 2004
   * 					2 	Pinson&Bazant model: linear diffusion + Tafel kinetics
   * 							ref: Pinson & Bazant, Journal of the Electrochemical society 160 (2), 2013
   * 					3	Christensen and Newman model: diffusion and kinetics model
   * 							ref: Christensen & Newmann, Journal of the Electrochemical Society 152 (4), 2005
   * SEI_porosity	integer deciding whether we reduce the volume fraction of active electrode material due to SEI growth
   *  				0	don't subtract anything
   * 					1	use correlation from Ashwin
   * 							ref: Ashwin, Chung, Wang, Journal of Power Sources 328, 2016
   * CS_id	identification for which model(s) to use for surface crack growth
   *  				0 	no surface cracking
   * 					1 	Laresgoiti's stress + crack growth model
   * 							ref: Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015
   * 					2 	Dai stress model + Laresgoiti crack growth
   * 							ref: Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015
   * 								 Dai, Cai, White, Journal of Power sources 247, 2014
   * 					3 	model based on Deshpande and Bernardi
   * 							ref: Deshpande & Bernardi,Journal of the Electrochemical Society 164 (2), 2017
   * 					4 	model from Barai
   * 							ref: Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
   * 					5 	model from Ekstrom
   * 							ref: Ekstrom and Lindbergh, Journal of the Electrochemical Society 162 (6), 2015
   * CS_diffusion	integer deciding whether to reduce the negative diffusion constant due to surface cracks
   *  				0 	don't decrease diffusion
   * 					1	decrease according to Barai
   * 							ref: Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
   * LAM_id  identification for which model(s) to use for loss of active material
   *  				0 	no LAM
   * 					1	Dai's stress model and Laresgoiti's correlation to relate the stress to LAM
   * 							ref: Dai, Cai, White, Journal of Power sources 247, 2014
   * 								 Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015
   * 					2	delacourt's correlation between abs(current density) and porosity
   * 							ref: Delacourt & Safari, Journal of the Electrochemical Society 159 (8), 2012
   * 					3 	Kindermann's model for cathode dissolution: tafel kinetics for increased porosity
   * 							ref: Kindermann, Keil, Frank, Jossen, Journal of the Electrochemical Society 164 (12), 2017
   * 					4 	Narayanrao's correlation between effective surface area and current density
   * 							ref: Narayanrao, Joglekar, Inguva, Journal of the Electrochemical Society 160 (1), 2012
   * pl_id 	identification for which model to use for li-plating
   *  				0 	no plating
   * 					1	Yang et al's model for thermodynamic plating (Tafel kinetics)
   * 							ref: Yang, Leng, Zhang, Ge, Wang, Journal of Power Sources 360, 2017
   */
  //!< For now, assume we want to include 'Pinson&Bazant'-SEI growth, 'Delacourt'-LAM, 'Kindermann'-LAM and 'Yang'-lithium plating
  // slide::DEG_ID deg;
  // deg.SEI_id.add_model(2); //!< Pinson & Bazant SEI growth
  // deg.SEI_porosity = 0;    //!< don't decrease the porosity (set to 1 if you do want to decrease the porosity)

  // deg.CS_id.add_model(0); //!< no surface cracks
  // deg.CS_diffusion = 0;   //!< don't decrease the diffusion coefficient (set to 1 if you do want to decrease the diffusion)

  // //!< there are 2 LAM models (Delacourt and Kindermann)
  // deg.LAM_id.add_model(2); //!< Delacourt LAM
  // deg.LAM_id.add_model(3); //!< Kindermann LAM

  // deg.pl_id = 1; //!< Yang lithium plating

  // std::cout << "Used ageing model: " << deg.print() << '\n';

  //!< Then the user has to choose what is simulated.
  //!< In the code below, uncomment the line which calls the function you want to execute (uncommenting means removing the //!<in front of the line)
  //!< and comment all the other lines (commenting means putting //!<in front of the line)

  auto numThreads = std::thread::hardware_concurrency(); //!< asd
  std::cout << "Available number of threads : " << numThreads << '\n';


  // Hopefully module is working:

  using namespace slide;

  slide::benchmarks::run_Cell_SPM_1(1);

  // paperCode::paper2022::Vequalisation();
  //   paperCode::paper2022::thermalModel();
  //  paperCode::paper2022::degradation_thermal();

  //!< Benchmarks:

  // slide::benchmarks::run_Cell_Bucket();
  // slide::benchmarks::run_Cell_ECM();
  // slide::benchmarks::run_Cell_SPM();

  // Timings:

  // slide::examples::estimatingOCVparameters();
  // slide::examples::drive_cycle_artemis();


  // std::unique_ptr<slide::StorageUnit> cs[2];
  // cs[0] = std::make_unique<slide::Cell_Bucket>();
  // cs[1] = std::make_unique<slide::Cell_Bucket>();

  // std::string n = "na";
  // double T = 300;
  // bool checkCells = false;

  // std::unique_ptr<slide::Module_s> mp(new slide::Module_s(n, T, true, false, 2, 1, 1));

  // mp->setSUs(cs, checkCells, true);

  // auto tst = mp->Vmin();
  // std::cout << tst << '\n';
  // std::cout << mp->Vmax() << '\n';
  // std::cout << mp->VMIN() << '\n';
  // std::cout << mp->VMAX() << '\n';

  //!< Slide-pack tests:

  //!< std::cout << "Cell_Bucket:\n";
  //!<{
  //!<	std::unique_ptr<slide::StorageUnit> cb_ptr = std::make_unique<slide::Cell_Bucket>();

  //!<	std::cout << "Voltage: " << cb_ptr->V() << '\n';

  //!<	cb_ptr->setCurrent(cb_ptr->Cap() / 2.0);
  //!<	std::cout << "Current: " << cb_ptr->I() << '\n';

  //!<	cb_ptr->timeStep_CC(2, 900);
  //!<	std::cout << "Voltage: " << cb_ptr->V() << '\n';

  //!<	auto u = cb_ptr->viewStates();

  //!<	for (auto x : u)
  //!<		std::cout << x << ' ';

  //!<	std::cout << "\n\n"
  //!<			  << std::endl;
  //!<}

  //!< std::cout << "Cell_Bucket with procedure:\n";
  //!<{
  //!<	bool balance = true;
  //!<	bool capspread = true;
  //!<	bool Rcellspread = true;
  //!<	bool degratespread = true;
  //!<	bool contactR = true;
  //!<	bool CV = false;
  //!<	double Vbal = 3.5; //!<rebalance to 3.1V
  //!<	int ndata = 60;	   //!<how often to store data (if 0, nothing is stored, not even usage stats)
  //!<	slide::Procedure p(balance, Vbal, ndata, false);
  //!<	std::cout << "Size of a Procedure: " << sizeof(p) << '\n';

  //!<	using namespace slide;
  //!<	constexpr int ncp = 2;
  //!<	constexpr int coolControl = 5;
  //!<	double Rc4[] = {0.0075 * 1e-3, 0.0073 * 1e-3};

  //!<	std::unique_ptr<StorageUnit> CinpM[2];

  //!<	CinpM[0] = std::make_unique<Cell_Bucket>();
  //!<	CinpM[1] = std::make_unique<Cell_Bucket>();

  //!<	auto mip = std::make_unique<Module_p>("p1", settings::T_ENV, true, false, ncp, coolControl, 1); //!<print warning messages, single-threaded, pass through coolsystem

  //!<	mip->setSUs(CinpM, false, true);
  //!<	mip->setRcontact(Rc4);

  //!<	p.useCaseAge(mip.get(), CV);
  //!<}

  //!< std::cout << "Cell_ECM:\n";
  //!<{
  //!<	std::unique_ptr<slide::StorageUnit> cb_ptr = std::make_unique<slide::Cell_ECM>();

  //!<	std::cout << "Voltage: " << cb_ptr->V() << '\n';

  //!<	cb_ptr->setCurrent(5);
  //!<	std::cout << "Current: " << cb_ptr->I() << '\n';

  //!<	cb_ptr->timeStep_CC(0.1, 3600);
  //!<	std::cout << "Voltage: " << cb_ptr->V() << '\n';

  //!<	auto u = cb_ptr->viewStates();

  //!<	for (auto x : u)
  //!<		std::cout << x << ' ';

  //!<	std::cout << "\n\n";
  //!<}

  //!< std::cout << "Cell_SPM:\n";
  //!<{
  //!<	auto cb_ptr = std::make_unique<slide::Cell_SPM>(); //!<std::unique_ptr<slide::StorageUnit>

  //!<	cb_ptr->deg_id = deg;

  //!<	std::cout << "Voltage: " << cb_ptr->V() << '\n';

  //!<	cb_ptr->setCurrent(cb_ptr->Cap() / 2.0);
  //!<	std::cout << "Current: " << cb_ptr->I() << '\n';

  //!<	cb_ptr->timeStep_CC(2, 900);
  //!<	std::cout << "Voltage: " << cb_ptr->V() << '\n';

  //!<	auto u = cb_ptr->viewStates();

  //!<	//!<for (auto x : u)
  //!<	//!<	std::cout << x << ' ';

  //!<	std::cout << "\n\n";
  //!<}

  //!< std::cout << "Cell_SPM and Cycler tests:\n";
  //!<{
  //!<	slide::Cycler cyc;

  //!<	auto cb_ptr = std::make_unique<slide::Cell_SPM>(); //!<std::unique_ptr<slide::StorageUnit>

  //!<	cb_ptr->deg_id = deg;

  //!<	std::string pref = "TEST";
  //!<	const bool diagnostic = true;

  //!<	cyc.initialise(cb_ptr.get(), pref);
  //!<	cyc.setDiagnostic(diagnostic);

  //!<	double Ah, Wh, v;
  //!<	const double dt = 1;
  //!<	constexpr unsigned Ncycle = 20; //!<365 * 10; //!<10 year
  //!<	constexpr unsigned Ncheck = 10; //!<do a checkup ever month
  //!<	constexpr unsigned Nbal = 7;	//!<balance every week

  //!<	std::vector<std::array<double, 3>> res(29);

  //!<	int ndt_data = 10;

  //!<	auto u = cb_ptr->viewStates();

  //!<	for (size_t i = 0; i < 29; i++)
  //!<		res[i][0] = u[i];

  //!<	cyc.CC(-cb_ptr->Cap(), cb_ptr->Vmin(), 100, dt, ndt_data, Ah, Wh); //!<4) 1C discharge

  //!<	for (size_t i = 0; i < 29; i++)
  //!<		res[i][1] = u[i];

  //!<	const auto succ = cyc.rest(4 * 3600, dt, ndt_data, Ah, Wh); //!<1) 4h rest
  //!<	if (succ != slide::Status::ReachedTimeLimit)
  //!<	{
  //!<		std::cout << "Error in useAge when resting 1 in cycle.\n";
  //!<	}

  //!<	for (size_t i = 0; i < 29; i++)
  //!<		res[i][2] = u[i];

  //!<	std::cout.setf(std::ios::fixed, std::ios::floatfield);

  //!<	for (auto item : res)
  //!<		std::cout << std::setprecision(10) << item[0] << "\t" << item[1] << "\t" << item[2] << '\n';

  //!<	std::cout << "\n\n";

  //!<	std::cout << "Ah: " << Ah << "   Wh: " << Wh << '\n';

  //!<	cyc.writeData();
  //!<}

  //!<{
  //!<	//!<Simulations with a large battery
  //!<	bool balance = true;
  //!<	bool capspread = true;
  //!<	bool Rcellspread = true;
  //!<	bool degratespread = true;
  //!<	bool contactR = true;
  //!<	bool CV = false;
  //!<	double Vbal = 3.5; //!<rebalance to 3.1V
  //!<	int ndata = 60;	   //!<how often to store data (if 0, nothing is stored, not even usage stats)
  //!<	slide::Procedure p(balance, Vbal, ndata, false);

  //!<	int coolControl = 5;
  //!<	//!<shared_ptr<StorageUnit> su = makeBattery(balance, capspread, Rcellspread, degratespread, contactR,coolControl); 	//!<20s15s9p, i.e. 20 cells in s module, 15 modules in s string, 9 strings in p battery
  //!<	//!<shared_ptr<StorageUnit> su = makeBattery2(balance, capspread, Rcellspread, degratespread, contactR,coolControl);	//!<9ps1520s, i.e. 9 cells in p module, 15 modules in s string, 20 strings in s battery
  //!<	//!<shared_ptr<StorageUnit> su = makeBattery_EPFL(capspread, Rcellspread, degratespread, contactR,coolControl,"",1);	//!<9ps1520s, i.e. 9 cells in p module, 15 modules in s string, 20 strings in s battery
  //!<	//!<p->cycleAge(su, CV);

  //!<	//!<loop for control settings, 1-2-3-4-5
  //!<	//	for(int coolControl = 1; coolControl < 6; coolControl++){
  //!<	//		try{
  //!<	//!<shared_ptr<StorageUnit> su = makeBattery(balance, capspread, Rcellspread, degratespread, contactR, coolControl); 	//!<20s15s9p, i.e. 20 cells in s module, 15 modules in s string, 9 strings in p battery
  //!<	//!<shared_ptr<StorageUnit> su = makeBattery2(balance, capspread, Rcellspread, degratespread, contactR,coolControl);	//!<9ps1520s, i.e. 9 cells in p module, 15 modules in s string, 20 strings in s battery
  //!<	//	auto su = slide::makeBattery_TestParallel(capspread, Rcellspread, degratespread, contactR, coolControl, "useCase", 1); //!<9ps1520s, i.e. 9 cells in p module, 15 modules in s string, 20 strings in s battery
  //!<	//!<p->cycleAge(su, CV);

  //!<	auto su = std::make_unique<slide::Cell_Bucket>(); //!<std::unique_ptr<slide::StorageUnit>

  //!<	p.useCaseAge(su.get(), CV);
  //!<	//		}catch(...){};
  //!<	//	}
  //!<}

  //!< std::cout << "Make battery test:\n";
  //!<{

  //!<	//!<Simulations with a large battery
  //!<	bool balance = true;
  //!<	bool capspread = true;
  //!<	bool Rcellspread = true;
  //!<	bool degratespread = true;
  //!<	bool contactR = true;
  //!<	bool CV = false;
  //!<	double Vbal = 3.5; //!<rebalance to 3.1V
  //!<	int ndata = 60;	   //!<how often to store data (if 0, nothing is stored, not even usage stats)

  //!<	//!<slide::Procedure p(balance, Vbal, ndata, false);

  //!<	int coolControl = 5;
  //!<	//

  //!<	//!<//!<Get the battery from makeBattery
  //!<	//!<int coolControl = 1; //!<just for naming since there is no thermal model
  //!<	//!<bool capspread = false;
  //!<	//!<bool Rcellspread = false;
  //!<	//!<bool degratespread = false;
  //!<	//!<bool contactR = true;
  //!<	//!<todo	std::shared_ptr<StorageUnit> su = makeBattery_EPFL(capspread, Rcellspread, degratespread, contactR,coolControl,IDaddition,1);

  //!<	//!<degrade the battery
  //!<	//!<todo	degrade(su);

  //!<	//!<Sensitivty with 10 times higher Rcontact
  //!<	std::string IDaddition = "_noThermalModel_HighRContact";
  //!<	auto cb_ptr = slide::makeBattery_EPFL(capspread, Rcellspread, degratespread, contactR, coolControl, IDaddition, 1);

  //!<	cb_ptr->setBlockDegAndTherm(true) ;

  //!<	//!<p.useCaseAge(su2.get(), CV);

  //!<	//!<auto cb_ptr = makeBattery(); //!<std::unique_ptr<slide::StorageUnit>

  //!<	//!<cb_ptr->deg_id = deg;

  //!<	std::cout << "Voltage: " << cb_ptr->V() << '\n';

  //!<	cb_ptr->setCurrent(5);
  //!<	std::cout << "Current: " << cb_ptr->I() << '\n';

  //!<	cb_ptr->timeStep_CC(1, 360);
  //!<	std::cout << "Voltage: " << cb_ptr->V() << '\n';

  //!<	//!<auto u = cb_ptr->viewStates();

  //!<	//!<//!<for (auto x : u)
  //!<	//!<//!<	std::cout << x << ' ';

  //!<	//!<std::cout << "\n\n";
  //!<}

  //!< std::cout << "Make battery test with Cyclers:\n";
  //!<{
  //!<	//!<Simulations with a large battery
  //!<	bool balance = true;
  //!<	bool capspread = true;
  //!<	bool Rcellspread = true;
  //!<	bool degratespread = true;
  //!<	bool contactR = true;
  //!<	bool CV = false;
  //!<	double Vbal = 3.5; //!<rebalance to 3.1V
  //!<	int ndata = 60;	   //!<how often to store data (if 0, nothing is stored, not even usage stats)
  //!<	//!<slide::Procedure p(balance, Vbal, ndata, false);
  //!<	int coolControl = 5;
  //!<	//
  //!<	//!<//!<Get the battery from makeBattery
  //!<	//!<int coolControl = 1; //!<just for naming since there is no thermal model
  //!<	//!<bool capspread = false;
  //!<	//!<bool Rcellspread = false;
  //!<	//!<bool degratespread = false;
  //!<	//!<bool contactR = true;
  //!<	//!<todo	std::shared_ptr<StorageUnit> su = makeBattery_EPFL(capspread, Rcellspread, degratespread, contactR,coolControl,IDaddition,1);
  //!<	//!<degrade the battery
  //!<	//!<todo	degrade(su);
  //!<	//!<Sensitivty with 10 times higher Rcontact
  //!<	std::string IDaddition = "_noThermalModel_HighRContact";
  //!<	auto cb_ptr = slide::makeBattery_EPFL(capspread, Rcellspread, degratespread, contactR, coolControl, IDaddition, 1);
  //!<	cb_ptr->setBlockDegAndTherm(true);
  //!<	//!<p.useCaseAge(su2.get(), CV);
  //!<	//!<auto cb_ptr = makeBattery(); //!<std::unique_ptr<slide::StorageUnit>
  //!<	//!<cb_ptr->deg_id = deg;
  //!<	slide::Cycler cyc;
  //!<	cyc.initialise(cb_ptr.get(), "VkCycler");
  //!<	std::cout << "Voltage: " << cb_ptr->V() << '\n';
  //!<	double Ah, Wh;
  //!<	cyc.CC(5, 4, 360, 1, 10, Ah, Wh);
  //!<	//!<cb_ptr->setCurrent(5);
  //!<	//!<std::cout << "Current: " << cb_ptr->I() << '\n';
  //!<	//!<cb_ptr->timeStep_CC(0.1 3600);
  //!<	std::cout << "Voltage: " << cb_ptr->V() << '\n';
  //!<	//!<auto u = cb_ptr->viewStates();
  //!<	//!<//!<for (auto x : u)
  //!<	//!<//!<	std::cout << x << ' ';
  //!<	//!<std::cout << "\n\n";
  //!< }

  //!<{
  //!<	//!<Simulations with a large battery
  //!<	bool balance = true;
  //!<	bool capspread = true;
  //!<	bool Rcellspread = true;
  //!<	bool degratespread = true;
  //!<	bool contactR = true;
  //!<	bool CV = false;
  //!<	double Vbal = 3.5; //!<rebalance to 3.1V
  //!<	int ndata = 60;	   //!<how often to store data (if 0, nothing is stored, not even usage stats)
  //!<	slide::Procedure p(balance, Vbal, ndata, false);

  //!<	int coolControl = 5;
  //!<	//!<shared_ptr<StorageUnit> su = makeBattery(balance, capspread, Rcellspread, degratespread, contactR,coolControl); 	//!<20s15s9p, i.e. 20 cells in s module, 15 modules in s string, 9 strings in p battery
  //!<	//!<shared_ptr<StorageUnit> su = makeBattery2(balance, capspread, Rcellspread, degratespread, contactR,coolControl);	//!<9ps1520s, i.e. 9 cells in p module, 15 modules in s string, 20 strings in s battery
  //!<	//!<shared_ptr<StorageUnit> su = makeBattery_EPFL(capspread, Rcellspread, degratespread, contactR,coolControl,"",1);	//!<9ps1520s, i.e. 9 cells in p module, 15 modules in s string, 20 strings in s battery
  //!<	//!<p->cycleAge(su, CV);

  //!<	//!<loop for control settings, 1-2-3-4-5
  //!<	//	for(int coolControl = 1; coolControl < 6; coolControl++){
  //!<	//		try{
  //!<	//!<shared_ptr<StorageUnit> su = makeBattery(balance, capspread, Rcellspread, degratespread, contactR, coolControl); 	//!<20s15s9p, i.e. 20 cells in s module, 15 modules in s string, 9 strings in p battery
  //!<	//!<shared_ptr<StorageUnit> su = makeBattery2(balance, capspread, Rcellspread, degratespread, contactR,coolControl);	//!<9ps1520s, i.e. 9 cells in p module, 15 modules in s string, 20 strings in s battery
  //!<	auto su = slide::makeBattery_TestParallel(capspread, Rcellspread, degratespread, contactR, coolControl, "useCase", 1); //!<9ps1520s, i.e. 9 cells in p module, 15 modules in s string, 20 strings in s battery
  //!<	//!<p->cycleAge(su, CV);
  //!<	p.useCaseAge(su.get(), CV);
  //!<	//		}catch(...){};
  //!<	//	}
  //!<}

  //!< Current Work section:

  //!<*********************************************** PARAMETRISATION FUNCTION CALLS *********************************************************************
  //!< slide::estimateOCVparameters(); //!<OCV parametrisation
  //!< slide::estimateCharacterisation(); //!<parametrisation of diffusion constant, rate constant and DC resistance

  //!<*********************************************** CYCLING FUNCTION CALLS ********************************************************
  //!< CCCV(M, pref, deg, cellType, settings::verbose); //!<a cell does a few CCCV cycles
  //!< FollowCurrent(M, pref, deg, cellType, settings::verbose); //!<a cell follows the current profile specified in a csv file

  //!<*********************************************** DEGRADATION FUNCTION CALLS ********************************************************
  //!< CalendarAgeing(M, pref, deg, cellType, settings::verbose); //!<simulates a bunch of calendar degradation experiments
  //!< CycleAgeing(M, pref, deg, cellType, settings::verbose); //!<simulates a bunch of cycle degradation experiments
  //!< ProfileAgeing(M, pref, deg, cellType, settings::verbose); //!<simulates a bunch of drive cycle degradation experiments

  //!<*********************************************** END ********************************************************
  //!< Now all the simulations have finished. Print this message, as well as how long it took to do the simulations
  std::cout << "finished all simulations in " << clk << ".\n";
}
