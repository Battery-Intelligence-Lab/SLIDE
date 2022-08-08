/*
 * CellSPM.cpp
 *
 *  Created on: 8 Feb 2020
 *   Author(s): Jorn Reniers
 */

#include "Cell_SPM.hpp"
#include "Interpolation.h"
#include "ReadCSVfiles.h"

#include "cassert"
#include "cmath"
#include "fstream"
#include "iostream"
#include <ctime>

Cell_SPM::Cell_SPM() : Cell()
{

	// ID string
	ID = "cell_SPM";

	// discretisation model
	Model_initialise(M);

	// parameters
	// maximum concentrations
	Cmaxpos = 51385;		  // value for NMC
	Cmaxneg = 30555;		  // value for C
	C_elec_sqrt = sqrt(1000); // standard concentration of 1 molar

	// constants
	n = 1;
	F = 96487;
	Rg = 8.314;

	// Cell parameters
	cap = 16;	// real capacity is 16.9248
	Vmax = 4.2; // value for an NMC/C cell
	VMAX = 4.3;
	Vmin = 2.7; // value for an NMC/C cell
	VMIN = 2.0;

	// thermal parameters
	T_ref = 273 + 25;
	Qch = 45; // 90;						// 90 gives very good cooling, as if there is a fan pointed at the cell. values of 30-60 are representative for a cell on a shelf without forced cooling
	rho = 1626;
	Cp = 750;
	Tmin = 273 + 0;
	Tmax = 273 + 60;

	// geometry
	L = 1.6850 * pow(10, -4);	 // thickness of one layer
	double width = 0.1;			 // width of the pouch
	double height = 0.2;		 // height of the pouch
	int nlayers = 31;			 // number of layers in the pouch
	Acell = width * height;		 // geometric surface area of the pouch
	elec_surf = Acell * nlayers; // total 'unrolled' surface area of the electrodes
	Rp = 8.5 * pow(10, -6);		 // do NOT CHANGE this value, if you do change it, you have to recalculate the spatial discretisation with the supplied Matlab scripts. See the word document '2 overview of the code', section 'Matlab setup before running the C++ code'
	Rn = 1.25 * pow(10, -5);	 // do NOT CHANGE this value, if you do change it, you have to recalculate the spatial discretisation with the supplied Matlab scripts. See the word document '2 overview of the code', section 'Matlab setup before running the C++ code'

	// main Li reaction
	kp = 5 * pow(10, -11); // fitting parameter
	kp_T = 58000;
	kn = 1.7640 * pow(10, -11); // fitting parameter (*2.5 for ageing fit)
	kn_T = 20000;
	// The diffusion coefficients at reference temperature are part of 'State'.
	// The values are set in the block of code below ('Initialise state variables')
	Dp_T = 29000;
	Dn_T = 35000.0 / 5.0; // ageing fit

	// State
	nstates = 2 * CELL_NCH + 18;
	double zp[CELL_NCH], zn[CELL_NCH];
	double fp, fn, delta, LLI, thickp, thickn, ep, en, ap, an, CS, Dp, Dn, rp, rn, rcc, delta_pl;
	fp = 0.689332;						 // lithium fraction in the cathode at 50% soc (3.68136 V) [-]
	fn = 0.479283;						 // lithium fraction in the anode at 50% soc (3.68136 V) [-]
	T = T_ENV;							 // cell temperature
	delta = pow(10, -9);				 // SEI thickness. Start with a fresh cell, which has undergone some formation cycles so it has an initial SEI layer.
										 // never start with a value of 0, because some equations have a term 1/delta, which would give nan or inf
										 // so this will give errors in the code
	LLI = 0;							 // lost lithium. Start with 0 so we can keep track of how much li we lose while cycling the cell
	thickp = 70 * pow(10, -6);			 // thickness of the positive electrode
	thickn = 73.5 * pow(10, -6);		 // thickness of the negative electrode
	ep = 0.5;							 // volume fraction of active material in the cathode
	en = 0.5;							 // volume fraction of active material in the anode
	ap = 3 * ep / Rp;					 // effective surface area of the cathode, the 'real' surface area is the product of the effective surface area (a) with the electrode volume (elec_surf * thick)
	an = 3 * en / Rn;					 // effective surface area of the anode
	CS = 0.01 * an * elec_surf * thickn; // initial crack surface. Start with 1% of the real surface area
	Dp = 8 * pow(10, -14);				 // diffusion constant of the cathode at reference temperature
	Dn = 7 * pow(10, -14);				 // diffusion constant of the anode at reference temperature
	rp = 0.0028;
	rn = 0.0028;
	rcc = 0.0002325; // note: this is divided by geometric surface (elec_surf) instead of effective surf (a*thick*elec_surf)
	delta_pl = 0;	 // thickness of the plated lithium layer. You can start with 0 here
	SoC = 0.5;
	I = 0;
	var_cap = 1;
	var_R = 1;
	var_degSEI = 1;
	var_degLAM = 1;

	// state-like parameters which remember values or cumulative effects
	Therm_Qgen = 0;
	Therm_Qgentot = 0;
	Therm_time = 0;
	Vcell = 0;
	Vcell_valid = false;
	etapcell = 0;
	etancell = 0;
	etacell_valid = 0;
	for (int i = 0; i < CELL_NSTATE_MAX; i++)
		degState[i] = 0;

	// Rdc is set at the end of this function since it requires a calculation for which you need all states
	// zp and zn: set to 0 and after setStates, use setC function to calculate the correct values
	for (int i = 0; i < CELL_NCH; i++)
	{
		zp[i] = 0;
		zn[i] = 0;
	}
	double sini[2 * CELL_NCH + 18];
	for (int i = 0; i > CELL_NCH; i++)
	{
		sini[i] = zp[i];
		sini[CELL_NCH + i] = zn[i];
	}
	sini[2 * CELL_NCH + 0] = delta;
	sini[2 * CELL_NCH + 1] = LLI;
	sini[2 * CELL_NCH + 2] = thickp;
	sini[2 * CELL_NCH + 3] = thickn;
	sini[2 * CELL_NCH + 4] = ep;
	sini[2 * CELL_NCH + 5] = en;
	sini[2 * CELL_NCH + 6] = ap;
	sini[2 * CELL_NCH + 7] = an;
	sini[2 * CELL_NCH + 8] = CS;
	sini[2 * CELL_NCH + 9] = Dp;
	sini[2 * CELL_NCH + 10] = Dn;
	sini[2 * CELL_NCH + 11] = rp;
	sini[2 * CELL_NCH + 12] = rn;
	sini[2 * CELL_NCH + 13] = rcc;
	sini[2 * CELL_NCH + 14] = delta_pl;
	sini[2 * CELL_NCH + 15] = SoC;
	sini[2 * CELL_NCH + 16] = T;
	sini[2 * CELL_NCH + 17] = I;

	setStates(sini, 2 * CELL_NCH + 18, false); // do not check that the states are valid, since we haven't set any parameters yet
	setC(fp, fn);							   // set the lithium concentration

	// OCV curves
	string namepos = "Kokam_OCV_NMC.csv";				 // name of the csv file with the cathode OCV curve
	string nameneg = "Kokam_OCV_C.csv";					 // name of the csv file with the anode OCV curve
	string nameentropicC = "Kokam_entropic_C.csv";		 // name of the csv file with the entropic coefficient curve for the anode OCV curve
	string nameentropicCell = "Kokam_entropic_cell.csv"; // name of the csv file with the entropic coefficient curve for the total cell OCV curve
	OCV_pos_n = 50;
	OCV_neg_n = 63;
	dOCV_neg_n = 11;
	dOCV_tot_n = 11;

	// shorter curves to save time
	OCV_pos_n = 22;						   // 49;
	OCV_neg_n = 24;						   // 63;
	namepos = "Kokam_OCV_NMC_reduced.csv"; // name of the csv file with the cathode OCV curve
	nameneg = "Kokam_OCV_C_reduced.csv";   // name of the csv file with the anode OCV curve

	// check that the OCV curves are shorter than the allowed range
	bool leng = false;
	if (OCV_pos_n >= CELL_NOCV_MAX)
	{
		std::cerr << "ERROR in the constructor of Cell_SPM. The cathode OCV curve has " << OCV_pos_n << " points while the maximum number is " << CELL_NOCV_MAX << ". Throwing an error" << '\n';
		leng = true;
	}
	if (OCV_neg_n >= CELL_NOCV_MAX)
	{
		std::cerr << "ERROR in the constructor of Cell_SPM. The anode OCV curve has " << OCV_neg_n << " points while the maximum number is " << CELL_NOCV_MAX << ". Throwing an error" << '\n';
		leng = true;
	}
	if (dOCV_neg_n >= CELL_NOCV_MAX)
	{
		std::cerr << "ERROR in the constructor of Cell_SPM. The anode entropic coefficient curve has " << dOCV_neg_n << " points while the maximum number is " << CELL_NOCV_MAX << ". Throwing an error" << '\n';
		leng = true;
	}
	if (dOCV_tot_n >= CELL_NOCV_MAX)
	{
		std::cerr << "ERROR in the constructor of Cell_SPM. The cell entropic coefficient curve has " << dOCV_tot_n << " points while the maximum number is " << CELL_NOCV_MAX << ". Throwing an error" << '\n';
		leng = true;
	}
	if (leng)
		throw 111;
	// if you get this error and you do want to use OCV curves which are longer than 1000 points, increase the size of the OCV arrays in Cell.hpp and the value of OCV_maxLength
	// Read the csv files
	try
	{
		loadCSV_2col(nameneg, OCV_neg_n, OCV_neg_x, OCV_neg_y);				// the OCV curve of the anode, the first column gives the lithium fractions (increasing), the 2nd column gives the OCV vs li/li+
		loadCSV_2col(namepos, OCV_pos_n, OCV_pos_x, OCV_pos_y);				// the OCV curve of the cathode, the first column gives the lithium fractions (increasing), the 2nd column gives the OCV vs li/li+
		loadCSV_2col(nameentropicC, dOCV_neg_n, dOCV_neg_x, dOCV_neg_y);	// the entropic coefficient of the anode, the first column gives the lithium fractions (increasing), the 2nd column gives the entropic coefficient [V K-1]
		loadCSV_2col(nameentropicCell, dOCV_tot_n, dOCV_tot_x, dOCV_tot_y); // the entropic coefficient of the entire cell, the first column gives the lithium fractions (increasing), the 2nd column gives the entropic coefficient [V K-1]
	}
	catch (int e)
	{
		cout << "ERROR in Cell_SPM::Cell_SPM when loading the OCV curves from the CSV files: " << e << ". Throwing it on" << endl
			 << flush;
		throw e;
	}

	// check if the inputs to the Matlab code are the same as the ones here in the C++ code
	// input:
	// 		M.Input[0] has to be the same as nch (defined in State.hpp)
	// 		M.Input[1] has to be the same as Rp (defined earlier in this constructor)
	// 		M.Input[2] has to be the same as Rn (defined earlier in this constructor)
	// 		M.input[3] has to give the location of the 0 eigenvalue
	bool Mnch = (M.Input[0] - CELL_NCH) / M.Input[0] > pow(10, -10); // allow a relative difference of e-10 due to numerical errors
	if (Mnch)
		std::cerr << "ERROR in Cell_SPM::Cell_SPM: the value of nch used in the Matlab script " << M.Input[0] << " is not the same as the value of nch used in the c++ code " << CELL_NCH << '\n';
	bool Mrp = (M.Input[1] - Rp) / M.Input[1] > pow(10, -10); // allow a relative difference of e-10 due to numerical errors
	if (Mrp)
		std::cerr << "ERROR in Cell_SPM::Cell_SPM: the value of Rp used in the Matlab script " << M.Input[1] << " is not the same as the value of Rp used in the c++ code " << Rp << '\n';
	bool Mrn = (M.Input[2] - Rn) / M.Input[2] > pow(10, -10); // allow a relative difference of e-10 due to numerical errors
	if (Mrn)
		std::cerr << "ERROR in Cell_SPM::Cell_SPM: the value of Rn used in the Matlab script " << M.Input[2] << " is not the same as the value of Rn used in the c++ code " << Rn << '\n';
	int a = M.Input[3];
	bool Meig = abs(M.An[a]) > pow(10, -10) || abs(M.Ap[a]) > pow(10, -10); // allow a relative difference of e-10 due to numerical errors
	if (Meig)
		std::cerr << "ERROR in Cell_SPM::Cell_SPM: the row of the 0-eigenvalue is " << M.Input[3] << " but that row has a positive eigenvalue of " << M.Ap[a] << " and negative eigenvalue of " << M.An[a] << ". They are not 0." << '\n';
	if (Mnch || Mrp || Mrn || Meig)
	{
		cout << "The Matlab script modelSetup.m produces matrices used by the C++ code for the discretisation of the solid diffusion equation."
				" Matlab needs some input parameters, such as the number of nodes and the radius of the particles. These parameters are specified on top of the Matlab scripts."
				" These values are also defined in the C++ code. Of course, both values have to be the same."
				" It turned out this was not the case, so you either have to change the values in the Matlab script or the ones in the C++ code."
				" We are throwing an error."
			 << endl
			 << flush;
		throw 10;
	}

	// Stress parameters
	omegap = 2.1 * pow(10, -6);	 // value from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
	omegan = 3.17 * pow(10, -6); // value from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
	Ep = 138.73;				 // value from Salco de Vasconcelos, Xu, Li, Zhao, Grid indentation analysis of mechanical properties of composite electrodes in Li-ion batteries, Extreme Mechanics Letters 9 (3), 2016
	En = 10;					 // value from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
	nup = 0.3;					 // value from Wu, Xiao, Wen, Zhang, Three-dimensional finite element study on stress generation in synchrotron X-ray tomography reconstructed nickel-manganese-cobalt based half cell, Journal of Power Sources, 336, 2016
	nun = 0.3;					 // value from Takahashi, Srinivasan, Examination of Graphite Particle Cracking as a Failure Mode in Lithium-Ion Batteries: A Model-Experimental Study, Journal of the Electrochemical Society, 162(4), 2015
	s_dai = false;
	s_lares = false;
	s_dai_update = false;
	s_lares_update = false;
	s_dai_p = 0;
	s_dai_n = 0;
	s_lares_n = 0;
	s_dai_p_prev = 0;
	s_dai_n_prev = 0;
	s_lares_n_prev = 0;
	s_dt = 1;

	// SEI parameters
	nsei = 1;
	alphasei = 1;
	OCVsei = 0.4;
	rhosei = 100 * pow(10, 3);
	rsei = 2037.4 * 50; // ageing fit
	Vmain = 13.0;
	Vsei = 64.39;
	c_elec0 = 4.541 / 1000;
	// fitting parameters of the models
	seiparam.sei1k = 0.075 * pow(10, -14);
	seiparam.sei1k_T = 130000;
	seiparam.sei2k = 2.75 * pow(10, -11);
	seiparam.sei2k_T = 130000;
	seiparam.sei2D = 1.125 * pow(10, -14);
	seiparam.sei2D_T = 20000;
	seiparam.sei3k = 1.1458 * pow(10, -15);
	seiparam.sei3k_T = 65000;
	seiparam.sei3D = 0.25 * pow(10, -15);
	seiparam.sei3D_T = 200000;
	seiparam.sei4k = (0.5 * 1 * pow(10, -14) / 2.0) * 1.5; // ageing fit
	seiparam.sei4k_T = 130000.0 / 1.5;					   // ageing fit
	seiparam.sei4D = 0.5 * pow(10, -16) / 15.0;			   // ageing fit
	seiparam.sei4D_T = 200000.0 / 2.5;					   // ageing fit
	seiparam.sei_porosity = 7.5 * pow(10, -7) * 3;		   // ageing fit

	// surface cracking
	csparam.CS1alpha = 4.25 * pow(10, -5);
	csparam.CS2alpha = 6.3 * pow(10, -7);
	csparam.CS3alpha = 2.31 * pow(10, -16);
	csparam.CS4alpha = 4.3306 * pow(10, -8);
	csparam.CS4Amax = 5 * getAn() * getThickn() * elec_surf; // assume the maximum crack surface is 5 times the initial anode surface
	csparam.CS5k = pow(10, -18);
	csparam.CS5k_T = -127040;
	csparam.CS_diffusion = 2;

	// LAM
	OCVnmc = 4.1;
	lamparam.lam1p = 3.4985 * pow(10, -9) / 10.0;		   // ageing fit
	lamparam.lam1n = 5.88 * pow(10, -13) * 2 * pow(10, 4); // ageing fit
	lamparam.lam2ap = -1.675 * pow(10, -5);
	lamparam.lam2bp = 0;
	lamparam.lam2an = -1.675 * pow(10, -5);
	lamparam.lam2bn = 0;
	lamparam.lam2t = 54611;
	lamparam.lam3k = 12.5 * pow(10, -6);
	lamparam.lam3k_T = 27305;
	lamparam.lam4p = 8.3333 * pow(10, -9);
	lamparam.lam4n = 8.3333 * pow(10, -9);

	// li-plating parameters
	npl = 1;
	alphapl = 1;
	OCVpl = 0;
	rhopl = 10000 * pow(10, 3);
	plparam.pl1k = 4.5 * pow(10, -10);
	plparam.pl1k_T = -2.014008 * pow(10, 5);

	// degradation identifiers: no degradation
	deg_id.SEI_id[0] = 0;	 // no SEI growth
	deg_id.SEI_n = 1;		 // there is 1 SEI model (namely '0')
	deg_id.SEI_porosity = 0; // don't decrease the volume fraction
	deg_id.CS_id[0] = 0;	 // no surface cracks
	deg_id.CS_n = 1;		 // there is 1 model (that there are no cracks)
	deg_id.CS_diffusion = 0; // don't decrease the diffusion coefficient
	deg_id.LAM_id[0] = 0;	 // no LAM
	deg_id.LAM_n = 1;		 // there is 1 LAM model
	deg_id.pl_id = 0;		 // no lithium plating

	// Check if we will have to calculate the stress according to Dai's stress model
	for (int i = 0; i < deg_id.CS_n; i++)
		s_dai = s_dai || deg_id.CS_id[i] == 2;
	for (int i = 0; i < deg_id.LAM_n; i++)
		s_dai = s_dai || deg_id.LAM_id[i] == 1;

	// check if we need to calculate the stress according to Laresgoiti's stress model
	for (int i = 0; i < deg_id.CS_n; i++)
		s_lares = s_lares || deg_id.CS_id[i] == 1;

	Rdc = getRdc(); // get the DC resistance

#if TIMING
	T_dstate = 0;
	T_getV = 0;
	T_getOCV = 0;
	T_validStates = 0;
	T_setStates = 0;
#endif
}
Cell_SPM::Cell_SPM(string IDi, const DEG_ID &degid, double capf, double resf, double degfsei, double degflam) : Cell_SPM()
{
	ID = IDi;

	// store such that the cell can later report its precise variation
	var_cap = capf;
	var_R = resf;
	var_degSEI = degfsei;
	var_degLAM = degflam;

	// Set the resistance
	setrdccc(getrdccc() * resf); // current collector
	setrdcp(getrdcp() * resf);	 // cathode
	setrdcn(getrdcn() * resf);	 // anode
	Rdc = getRdc();				 // set the DC resistance

	// set the capacity
	cap = cap * capf;			  // nominal capacity
	elec_surf = elec_surf * capf; // surface area of the electrodes (current to current density)

	// set the degradation factor
	seiparam.sei1k = seiparam.sei1k * var_degSEI;
	seiparam.sei1k_T = seiparam.sei1k_T * var_degSEI;
	seiparam.sei2k = seiparam.sei2k * var_degSEI;
	seiparam.sei2k_T = seiparam.sei2k_T * var_degSEI;
	seiparam.sei2D = seiparam.sei2D * var_degSEI;
	seiparam.sei2D_T = seiparam.sei2D_T * var_degSEI;
	seiparam.sei3k = seiparam.sei3k * var_degSEI;
	seiparam.sei3k_T = seiparam.sei3k_T * var_degSEI;
	seiparam.sei3D = seiparam.sei3D * var_degSEI;
	seiparam.sei3D_T = seiparam.sei3D_T * var_degSEI;
	seiparam.sei4k = seiparam.sei4k * var_degSEI;
	seiparam.sei4k_T = seiparam.sei4k_T * var_degSEI;
	seiparam.sei4D = seiparam.sei4D * var_degSEI;
	seiparam.sei4D_T = seiparam.sei4D_T * var_degSEI;
	seiparam.sei_porosity = seiparam.sei_porosity * var_degSEI;
	csparam.CS1alpha = csparam.CS1alpha * var_degSEI;
	csparam.CS2alpha = csparam.CS2alpha * var_degSEI;
	csparam.CS3alpha = csparam.CS3alpha * var_degSEI;
	csparam.CS4alpha = csparam.CS4alpha * var_degSEI;
	csparam.CS4Amax = csparam.CS4Amax * var_degSEI;
	csparam.CS5k = csparam.CS5k * var_degSEI;
	csparam.CS5k_T = csparam.CS5k_T * var_degSEI;
	csparam.CS_diffusion = csparam.CS_diffusion * var_degSEI;
	lamparam.lam1p = lamparam.lam1p * var_degLAM;
	lamparam.lam1n = lamparam.lam1n * var_degLAM;
	lamparam.lam2ap = lamparam.lam2ap * var_degLAM;
	lamparam.lam2bp = lamparam.lam2bp * var_degLAM;
	lamparam.lam2an = lamparam.lam2an * var_degLAM;
	lamparam.lam2bn = lamparam.lam2bn * var_degLAM;
	lamparam.lam2t = lamparam.lam2t * var_degLAM;
	lamparam.lam3k = lamparam.lam3k * var_degLAM;
	lamparam.lam3k_T = lamparam.lam3k_T * var_degLAM;
	lamparam.lam4p = lamparam.lam4p * var_degLAM;
	lamparam.lam4n = lamparam.lam4n * var_degLAM;
	plparam.pl1k = plparam.pl1k * var_degSEI;
	plparam.pl1k_T = plparam.pl1k_T * var_degSEI;

	// set the degradation ID and related settings
	deg_id = degid;
	bockDegAndTherm = false;
	// Check if we will have to calculate the stress according to Dai's stress model
	for (int i = 0; i < deg_id.CS_n; i++)
		s_dai = s_dai || deg_id.CS_id[i] == 2;
	for (int i = 0; i < deg_id.LAM_n; i++)
		s_dai = s_dai || deg_id.LAM_id[i] == 1;
	// check if we need to calculate the stress according to Laresgoiti's stress model
	for (int i = 0; i < deg_id.CS_n; i++)
		s_lares = s_lares || deg_id.CS_id[i] == 1;
}

Cell_SPM::~Cell_SPM() {}

double *Cell_SPM::getZp() { return &state[0]; }
double *Cell_SPM::getZn() { return &state[CELL_NCH]; }
double Cell_SPM::getCp_surface(double *zpi)
{
	/*
	 * Get the surface concentration in the cathodic particle
	 *
	 * zpi 	optional argument
	 * 		if none provided, zp is used and the surface concentration of the current state is calculated
	 * 		if a different array is provided, the surface concentration for this different array is calculated
	 * 			i.e. what the surface concentration would be if the (transformed) concentration is set to zpi
	 */

	// Calculate the diffusion constant at the battery temperature using an Arrhenius relation
	double Dpt = getDp() * exp(Dp_T / Rg * (1 / T_ref - 1 / getT())); // diffusion constant of the positive particle [m s-1]

	// Calculate the molar flux on the surfaces
	double jp = -getI() / (getAp() * elec_surf * getThickp()) / (n * F); // molar flux on the positive particle [mol m-2 s-1]

	// Calculate the surface concentration at the positive particle
	// 	cp_surf = M.Cp[0][:] * zp[:] + M.Dp*jp/Dpt
	double cp_surf = 0;
	for (int j = 0; j < CELL_NCH; j++)
		cp_surf += M.Cp[0][j] * zpi[j];
	cp_surf += M.Dp[0] * jp / Dpt;

	// Make the output parameters
	return cp_surf;
}
double Cell_SPM::getCn_surface(double *zni)
{
	/*
	 * Get the surface concentration in the anodic particle
	 */

	// Calculate the diffusion constant at the battery temperature using an Arrhenius relation
	double Dnt = getDn() * exp(Dn_T / Rg * (1 / T_ref - 1 / getT())); // diffusion constant of the negative particle [m s-1]

	// Calculate the molar flux on the surfaces
	double jn = getI() / (getAn() * elec_surf * getThickn()) / (n * F); // molar flux on the negative particle [mol m-2 s-1]

	// Calculate the surface concentration at the negative particle
	// 	cn_surf = M.Cn[0][:] * zn[:] + M.Dn*jn/Dnt
	double cn_surf = 0;
	for (int j = 0; j < CELL_NCH; j++)
		cn_surf += M.Cn[0][j] * zni[j];
	cn_surf += M.Dn[0] * jn / Dnt;

	// Make the output parameters
	return cn_surf;
}
double Cell_SPM::getDelta() { return state[2 * CELL_NCH + 0]; }
double Cell_SPM::getLLI() { return state[2 * CELL_NCH + 1]; }
double Cell_SPM::getThickp() { return state[2 * CELL_NCH + 2]; }
double Cell_SPM::getThickn() { return state[2 * CELL_NCH + 3]; }
double Cell_SPM::getEp() { return state[2 * CELL_NCH + 4]; }
double Cell_SPM::getEn() { return state[2 * CELL_NCH + 5]; }
double Cell_SPM::getAp() { return state[2 * CELL_NCH + 6]; }
double Cell_SPM::getAn() { return state[2 * CELL_NCH + 7]; }
double Cell_SPM::getCS() { return state[2 * CELL_NCH + 8]; }
double Cell_SPM::getDp() { return state[2 * CELL_NCH + 9]; }
double Cell_SPM::getDn() { return state[2 * CELL_NCH + 10]; }
double Cell_SPM::getrdcp() { return state[2 * CELL_NCH + 11]; }
double Cell_SPM::getrdcn() { return state[2 * CELL_NCH + 12]; }
double Cell_SPM::getrdccc() { return state[2 * CELL_NCH + 13]; }
double Cell_SPM::getDelta_pl() { return state[2 * CELL_NCH + 14]; }
double Cell_SPM::getRdc()
{

	// calculate the resistance from every component
	double Rdc_sei = getDelta() * rsei / (getAn() * elec_surf * getThickn());
	double Rdc_p = getrdcp() / (getAp() * elec_surf * getThickp());
	double Rdc_n = getrdcn() / (getAn() * elec_surf * getThickn());
	double Rdc_cc = getrdccc() / (elec_surf);

	return Rdc_sei + Rdc_p + Rdc_n + Rdc_cc;
}
double Cell_SPM::getRtot()
{
	return getRdc();
}
void Cell_SPM::getConcentration(int nin, double cp[], double cn[])
{
	/*
	 * Calculates the lithium concentration at each positive chebyshev node (0 <= x <= 1), including the centre and surface nodes.
	 * Uses the matrices from the state space matrices.
	 * See the equations in the explanatory documents.
	 *
	 * IN
	 * nin 	length of the arrays provided
	 *
	 * OUT
	 * cp 	li-concentration at each chebyshev node in the positive electrode [mol m-3], length of the array should be CELL_NCH+2
	 * 			cp[0]			concentration at the surface of the sphere
	 * 			cp[1 to nch]	concentration at the inner nodes
	 * 			cp[nch + 1]		concentration at the centre of the sphere
	 * cn 	li-concentration at each chebyshev node in the negative electrode [mol m-3], length of the array should be CELL_NCH+2
	 * 			cn[0]			concentration at the surface of the sphere
	 * 			cn[1 to nch]	concentration at the inner nodes
	 * 			cn[nch + 1]		concentration at the centre of the sphere
	 *
	 * THROWS
	 * 100 	the arrays provided have the wrong length
	 */

	if (nin != CELL_NCH + 2)
	{ // nch is the number of inner nodes, +2 for the node at the centre and surface
		std::cerr << "ERROR in Cell::getC. The length of the array provided is " << nin << " instead of " << CELL_NCH + 2 << '\n';
		throw 100;
	}

	// Calculate the diffusion constant at the battery temperature using an Arrhenius relation
	double Dpt = getDp() * exp(Dp_T / Rg * (1 / T_ref - 1 / getT())); // diffusion constant of the positive particle [m s-1]
	double Dnt = getDn() * exp(Dn_T / Rg * (1 / T_ref - 1 / getT())); // diffusion constant of the negative particle [m s-1]

	// Calculate the molar flux on the surfaces
	double jp = -getI() / (getAp() * elec_surf * getThickp()) / (n * F); // molar flux on the positive particle [mol m-2 s-1]
	double jn = getI() / (getAn() * elec_surf * getThickn()) / (n * F);	 // molar flux on the negative particle [mol m-2 s-1]

	// Calculate concentration at the surface and inner nodes using the matrices from the spatial discretisation of the solid diffusion PDE
	// 	cp = M.Cp[:][:] * zp[:] + M.Dp*jp/Dpt
	// 	cn = M.Cn[:][:] * zn[:] + M.Dn*jn/Dnt
	double cpt, cnt;
	for (int i = 0; i < CELL_NCH + 1; i++)
	{ // loop to calculate at each surface + inner node
		cpt = 0.0;
		cnt = 0.0;
		for (int j = 0; j < CELL_NCH; j++)
		{
			cpt += M.Cp[i][j] * getZp()[j];
			cnt += M.Cn[i][j] * getZn()[j];
		}
		cp[i] = cpt + M.Dp[i] * jp / Dpt;
		cn[i] = cnt + M.Dn[i] * jn / Dnt;
	}

	// Calculate the concentration at centre node using the boundary condition (the concentration gradient at the centre has to be 0 due to symmetry)
	// cp_centre = -1/2 (M.Cc[:]*cp +jp*Rp/Dpt)
	// cn_centre = -1/2 (M.Cc[:]*cn +jn*Rn/Dnt)
	double DM = 2.0; // we need a constant of 2 in the equations
	cpt = 0.0;
	cnt = 0.0;
	for (int i = 0; i < CELL_NCH + 1; i++)
	{
		cpt += M.Cc[i] * cp[i];
		cnt += M.Cc[i] * cn[i];
	}
	cp[CELL_NCH + 1] = (-1.0 / DM) * (cpt + jp * Rp / Dpt);
	cn[CELL_NCH + 1] = (-1.0 / DM) * (cnt + jn * Rn / Dnt);
}

void Cell_SPM::getDaiStress(int nin, double *sigma_p, double *sigma_n,
							double sigma_r_p[], double sigma_r_n[], double sigma_t_p[],
							double sigma_t_n[], double sigma_h_p[], double sigma_h_n[])
{
	/*
	 * Calculates the radial and tangential stress for each positive Chebyshev node according to the formula by
	 * Dai, Cai, White, Journal of Power sources 247, 2014
	 *
	 * It takes quite long to calculate the stress, so only call this function when needed.
	 *
	 * IN
	 * nin 			length of the arrays provided
	 *
	 * OUT
	 * sigma_p		maximum hydrostatic stress in the positive particle, can be both positive and negative [Pa]
	 * sigma_n		maximum hydrostatic stress in the negative particle, can be both positive and negative [Pa]
	 * sigma_r_p	array with the radial stress at each positive chebyshev node in the positive electrode, length CELL_NCH+2, [Pa]
	 * sigma_r_n	array with the radial stress at each positive chebyshev node in the negative electrode, length CELL_NCH+2, [Pa]
	 * sigma_t_p	array with the tangential stress at each positive chebyshev node in the positive electrode, length CELL_NCH+2, [Pa]
	 * sigma_t_n	array with the tangential stress at each positive chebyshev node in the negative electrode, length CELL_NCH+2, [Pa]
	 * sigma_h_p	array with the hydrostatic stress at each positive chebyshev node in the positive electrode, length CELL_NCH+2, [Pa]
	 * sigma_h_n	array with the hydrostatic stress at each positive chebyshev node in the negative electrode, length CELL_NCH+2, [Pa]
	 * 				[0]			stress at the surface of the sphere
	 * 				[1 to nch]	stress at the inner nodes
	 * 				[nch + 1]	stress at the centre of the sphere
	 *
	 * THROWS
	 * 100 	the arrays provided have the wrong length
	 */

	// check the arrays have the correct length
	if (nin != CELL_NCH + 2)
	{
		std::cerr << "ERROR in Cell::getDaiStress. The length of the array provided is " << nin << " instead of " << CELL_NCH + 2 << '\n';
		throw 100;
	}

	// Get the locations of the Chebyshev nodes
	double xp[CELL_NCH + 2]; // location (x-value) of the positive Chebyshev nodes
	xp[0] = 1;
	for (int i = 0; i < CELL_NCH; i++)
		xp[i + 1] = M.xch[i];
	xp[CELL_NCH + 1] = 0;
	double xtot[2 * CELL_NCH + 3];		   // location (x-value) of the positive and negative Chebyshev nodes [-surface .. centre .. +surface]
	xtot[CELL_NCH + 1] = xp[CELL_NCH + 1]; // centre node
	for (int i = 0; i < CELL_NCH + 1; i++)
	{
		xtot[i] = -xp[i];						   // negative nodes
		xtot[CELL_NCH + 2 + i] = xp[CELL_NCH - i]; // positive nodes
	}

	// get concentrations at each Chebyshev node
	// Due to symmetry, the concentration at the negative point is the same as the concentration of the positive point: c(-x) = c(x)
	double cp[CELL_NCH + 2], cn[CELL_NCH + 2]; // positive and negative nodes, [+surface .. inner .. centre]
	getConcentration(CELL_NCH + 2, cp, cn);
	double CP[2 * CELL_NCH + 3], CN[2 * CELL_NCH + 3]; // concentrations at all nodes, [-surface .. inner .. centre .. inner .. +surface]
	CP[CELL_NCH + 1] = cp[CELL_NCH + 1];			   // cathode centre node
	CN[CELL_NCH + 1] = cn[CELL_NCH + 1];			   // anode centre node
	for (int i = 0; i < CELL_NCH + 1; i++)
	{
		CP[i] = cp[i];							 // cathode negative points
		CN[i] = cn[i];							 // anode negative points
		CP[CELL_NCH + 2 + i] = cp[CELL_NCH - i]; // cathode positive points
		CN[CELL_NCH + 2 + i] = cn[CELL_NCH - i]; // anode positive points
	}

	// The formula's to calculate the stress have integrals.
	// Integrals of Chebyshev points can be calculated using the Q-matrix in the state space struct (M)
	// The integral from the negative surface to node i is given by row i of the product Q*f
	// 		with Q the Chebyshev integration matrix
	// 			 f the value of the function you want to integrate, evaluated at every node
	// All integrals have to start from the negative surface (because you must cover the entire Chebyshev domain)
	// 	so if you need the integral of a function f from the centre (x = 0) to a positive point in the sphere (x = i)
	// 		int(f, x = 0 .. i) = int(f, x=-1 .. i) - int(f, x=-1 .. 0)
	// E.g. to get the integral of the positive li-concentration from the centre until the 4th positive Chebyshev node:
	// 		F = Q * CP 				[-surface .. centre .. +surface]
	// 		int(c(x), x = 0 .. i(4)) = int(c(x), x=-1 .. i(4)) - int(c(x), x=-1 .. 0)
	// 							  	 = F[nch + 1 + 4] 		   - F[CELL_NCH+1]
	// 		(remember that the centre node is at [CELL_NCH+1])

	// Calculate the matrix-vector product of the required functions as given in the paper by Dai et al. (concentration * radius^2)
	// we need to remember the transformation of variables from x (-1<x<1) to r (-R<r<R)
	// 		int( c * r^2 dr) = int(c * (x*R)^2 * d(R*x)) = int(c x^2 R^3 dx)
	// so the matrix-vector product we need is F = Q * (concentration * x^2 * R^3)
	double Fp[2 * CELL_NCH + 3], Fn[2 * CELL_NCH + 3]; // arrays with the product for the positive and negative electrode
	for (int i = 0; i < 2 * CELL_NCH + 3; i++)
	{			   // loop for each row (one row = one node)
		Fp[i] = 0; // calculate the matrix-vector product for row i as you would do it by hand:
		Fn[i] = 0; // 	F(i) = sum(Q(i,j)*C(j)*x(j)^2*R^3, j=0..2*CELL_NCH+2)
		for (int j = 0; j < 2 * CELL_NCH + 3; j++)
		{													  // loop through the columns to calculate the sum
			Fp[i] += M.Q[i][j] * (CP[j] * xtot[j] * xtot[j]); // 		Q(i,j)*C(j)*x(j)^2
			Fn[i] += M.Q[i][j] * (CN[j] * xtot[j] * xtot[j]);
		}
		Fp[i] = Fp[i] * Rp * Rp * Rp; // *R^3 (R is constant so it can be out of the sum)
		Fn[i] = Fn[i] * Rn * Rn * Rn;
	}

	// Calculate the integral from the centre to the positive surface, which is a constant present in all equations
	double ap = Fp[2 * CELL_NCH + 2] - Fp[CELL_NCH + 1]; // int( cp*r^2, r=0..Rp )
	double an = Fn[2 * CELL_NCH + 2] - Fn[CELL_NCH + 1]; // int( cn*r^2, r=0..Rn )

	// Calculate the equations for all nodes
	double rp;				  // radius of positive node i in the positive particle
	double rn;				  // radius of positive node i in the negative particle
	double bp;				  // integral from the centre to positive node i, int(cp*zp^2, zp=0..rp(i) )
	double bn;				  // integral from the centre to positive node i, int(cn*zn^2, zn=0..rn(i) )
	double srp[CELL_NCH + 2]; // radial stress in the positive particle at the positive nodes [centre .. +surface]
	double srn[CELL_NCH + 2]; // radial stress in the negative particle at the positive nodes [centre .. +surface]
	double stp[CELL_NCH + 2]; // tangential stress in the positive particle at the positive nodes [centre .. +surface]
	double stn[CELL_NCH + 2]; // tangential stress in the negative particle at the positive nodes [centre .. +surface]
	for (int i = 0; i < CELL_NCH + 2; i++)
	{									  // loop for the positive nodes
		rp = Rp * xtot[CELL_NCH + 1 + i]; // r(i) = R * x(i)
		rn = Rn * xtot[CELL_NCH + 1 + i];
		bp = Fp[CELL_NCH + 1 + i] - Fp[CELL_NCH + 1]; // int( cp*z^2, z=0..r(i) ) = F[CELL_NCH+1+i] - F[CELL_NCH+1]
		bn = Fn[CELL_NCH + 1 + i] - Fn[CELL_NCH + 1];

		// Implement the equations from Dai et al.
		if (i == 0)
		{ // centre node -> special formula (31 & 33) in Dai, Cai, White
			srp[i] = 2 * omegap * Ep / (9 * (1 - nup)) * (3 / pow(Rp, 3) * ap - CP[CELL_NCH + 1]);
			srn[i] = 2 * omegan * En / (9 * (1 - nun)) * (3 / pow(Rn, 3) * an - CN[CELL_NCH + 1]);

			stp[i] = 2 * omegap * Ep / (9 * (1 - nup)) * (3 / pow(Rp, 3) * ap - CP[CELL_NCH + 1]);
			stn[i] = 2 * omegan * En / (9 * (1 - nun)) * (3 / pow(Rn, 3) * an - CN[CELL_NCH + 1]);
		}
		else
		{																							  // other nodes -> equation 13 in Dai, Cai, White
			srp[i] = 2 * omegap * Ep / (3 * (1 - nup)) * (1 / pow(Rp, 3) * ap - 1 / pow(rp, 3) * bp); // ap = int (c x^2, x=0..R), bp = int (c x^2 , x=0..r)
			srn[i] = 2 * omegan * En / (3 * (1 - nun)) * (1 / pow(Rn, 3) * an - 1 / pow(rn, 3) * bn);

			stp[i] = omegap * Ep / (3 * (1 - nup)) * (2 / pow(Rp, 3) * ap + 1 / pow(rp, 3) * bp - cp[i]);
			stn[i] = omegan * En / (3 * (1 - nun)) * (2 / pow(Rn, 3) * an + 1 / pow(rn, 3) * bn - cn[i]);
		}
	}

	// Flip all arrays to get the opposite order (now it is [centre .. +surface] and we want [+surface .. centre]
	// and store in the output arrays
	for (int i = 0; i < CELL_NCH + 2; i++)
	{ // loop for the positive nodes
		sigma_r_p[i] = srp[CELL_NCH + 2 - 1 - i];
		sigma_r_n[i] = srn[CELL_NCH + 2 - 1 - i];
		sigma_t_p[i] = stp[CELL_NCH + 2 - 1 - i];
		sigma_t_n[i] = stn[CELL_NCH + 2 - 1 - i];
	}

	// Make the hydrostatic stress sh = (sr + 2sp)/3
	int sp = 0; // node with the maximum hydrostatic stress in the positive particle
	int sn = 0; // node with the maximum hydrostatic stress in the negative
	for (int i = 0; i < CELL_NCH + 2; i++)
	{														  // loop for all nodes
		sigma_h_p[i] = (sigma_r_p[i] + 2 * sigma_t_p[i]) / 3; // calculate hydrostatic stress
		sigma_h_n[i] = (sigma_r_n[i] + 2 * sigma_t_n[i]) / 3;

		// find the maximum (in absolute value) of the stresses
		if (abs(sigma_h_p[i]) > abs(sigma_h_p[sp]))
			sp = i;
		if (abs(sigma_h_n[i]) > abs(sigma_h_n[sp]))
			sn = i;
	}
	*sigma_p = sigma_h_p[sp];
	*sigma_n = sigma_h_n[sn];
}

void Cell_SPM::updateDaiStress()
{
	/*
	 * Function which will update the values stored in the stress variables relating with Dai's stress model
	 */

	// Make variables to store the stress
	double sigma_r_p[CELL_NCH + 2], sigma_r_n[CELL_NCH + 2], sigma_t_p[CELL_NCH + 2], sigma_t_n[CELL_NCH + 2], sigma_h_p[CELL_NCH + 2], sigma_h_n[CELL_NCH + 2];
	double sigma_p, sigma_n;

	// Get the stress
	try
	{
		getDaiStress(CELL_NCH + 2, &sigma_p, &sigma_n, sigma_r_p, sigma_r_n, sigma_t_p, sigma_t_n, sigma_h_p, sigma_h_n);
	}
	catch (int e)
	{
		cout << "Error in Cell::getDaiStress when calling updateDaiStress: " << e << ". throwing it on" << endl
			 << flush;
		throw e;
	}

	// Update the stored values
	s_dai_p = sigma_p;
	s_dai_n = sigma_n;

	// indicate that the values in the class variables are updated
	s_dai_update = true;
}

void Cell_SPM::getLaresgoitiStress(bool print, double *sigma_n)
{
	/*
	 * Calculate the stress according to Laresgoiti's stress model
	 *
	 * IN
	 * print 	boolean indicating if we want to print error messages or not
	 * 				if true, error messages are printed
	 * 				if false no error messages are printed (but the error will still be thrown)
	 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
	 *
	 * OUT
	 * sigma_n	stress in the negative particle [MPa]
	 *
	 * THROWS
	 * 101 		the surface concentration is out of bounds
	 */

	// Arrays with the stress from Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015 (figure 5)
	double xx[11] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};		// li-fraction in the graphite
	double yy[11] = {0.0, 5.5, 8.5, 9.5, 10.0, 10.0, 10.5, 13.0, 16.5, 21.0, 23.5}; // stress [MPa]

	// Get the surface concentration
	double cps = getCp_surface(getZp());
	double cns = getCn_surface(getZn());
	double zn_surf = (cns / Cmaxneg); // lithium fraction on negative surface [0 1]

	// check if the surface concentration is within the allowed range
	// 	0 < cp < Cmaxpos
	// 	0 < cn < Cmaxneg
	if (cps < 0 || cns < 0 || cps > Cmaxpos || cns > Cmaxneg)
	{
		if (print)
		{
			std::cerr << "ERROR in Cell::getLaresgoitiStress: concentration out of bounds. the positive lithium fraction is " << cps / Cmaxpos << " and the negative lithium fraction is " << cns / Cmaxneg;
			std::cerr << "they should both be between 0 and 1" << '\n';
		}
		throw 101;
	}

	// Interpolate linearly to get the stress
	double s;
	try
	{
		s = linInt(print, true, xx, yy, 11, zn_surf); // throw an error if you are out of the allowed range
	}
	catch (int e)
	{
		if (print)
			cout << "Error in Cell::getLaresgoitiStress when interpolating in the arrays " << e << ". Throwing it on" << endl
				 << flush;
		throw e;
	}

	// Make the output variable
	*sigma_n = s;
}

void Cell_SPM::updateLaresgoitiStress(bool print)
{
	/*
	 * Function which will update the values stored in the stress variables relating with Laresgoiti's stress model
	 *
	 * IN
	 * print 	boolean indicating if we want to print error messages or not
	 * 				if true, error messages are printed
	 * 				if false no error messages are printed (but the error will still be thrown)
	 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
	 */

	double s;
	getLaresgoitiStress(print, &s);

	// Update the stored value
	s_lares_n = s;
	s_lares_update = true; // indicate that the values in the class variables are updated
}

void Cell_SPM::setZp(double in[])
{
	// set the (transformed) concentration to a new value
	for (int i = 0; i < CELL_NCH; i++)
		state[i] = in[i];
	Vcell_valid = false;
};
void Cell_SPM::setZn(double in[])
{
	for (int i = 0; i < CELL_NCH; i++)
		state[CELL_NCH + i] = in[i];
	Vcell_valid = false;
};
void Cell_SPM::setdZp(double din[])
{
	// add a value to the (transformed) concentration
	for (int i = 0; i < CELL_NCH; i++)
		state[i] += din[i];
	Vcell_valid = false;
};
void Cell_SPM::setdZn(double din[])
{
	for (int i = 0; i < CELL_NCH; i++)
		state[CELL_NCH + i] += din[i];
	Vcell_valid = false;
};
void Cell_SPM::setDelta(double in)
{
	state[2 * CELL_NCH + 0] = in;
	Vcell_valid = false;
}
void Cell_SPM::setLLI(double in)
{
	state[2 * CELL_NCH + 1] = in;
	Vcell_valid = false;
}
void Cell_SPM::setThickp(double in)
{
	state[2 * CELL_NCH + 2] = in;
	Vcell_valid = false;
}
void Cell_SPM::setThickn(double in)
{
	state[2 * CELL_NCH + 3] = in;
	Vcell_valid = false;
}
void Cell_SPM::setEp(double in)
{
	state[2 * CELL_NCH + 4] = in;
	Vcell_valid = false;
}
void Cell_SPM::setEn(double in)
{
	state[2 * CELL_NCH + 5] = in;
	Vcell_valid = false;
}
void Cell_SPM::setAp(double in)
{
	state[2 * CELL_NCH + 6] = in;
	Vcell_valid = false;
}
void Cell_SPM::setAn(double in)
{
	state[2 * CELL_NCH + 7] = in;
	Vcell_valid = false;
}
void Cell_SPM::setCS(double in) { state[2 * CELL_NCH + 8] = in; }
void Cell_SPM::setDp(double in)
{
	state[2 * CELL_NCH + 9] = in;
	Vcell_valid = false;
}
void Cell_SPM::setDn(double in)
{
	state[2 * CELL_NCH + 10] = in;
	Vcell_valid = false;
}
void Cell_SPM::setrdcp(double in)
{
	state[2 * CELL_NCH + 11] = in;
	Vcell_valid = false;
}
void Cell_SPM::setrdcn(double in)
{
	state[2 * CELL_NCH + 12] = in;
	Vcell_valid = false;
}
void Cell_SPM::setrdccc(double in)
{
	state[2 * CELL_NCH + 13] = in;
	Vcell_valid = false;
}
void Cell_SPM::setDelta_pl(double in)
{
	state[2 * CELL_NCH + 14] = in;
	Vcell_valid = false;
}
double Cell_SPM::setSoC(double SoCnew, bool checkV, bool print)
{
	/*
	 * Overwrite the implementation from Cell since we don't want to throw errors if SoC > 1 or < 0
	 * in Cell, SoC determines OCV so must stay within limits
	 * here, that is done by surface concentration.
	 * SoC is integration of current, and if the nominal capacity or initial SoC is wrong, SoC might become larger than 1 or smaller than 0.
	 * However, it causes no problems in the code.
	 *
	 * In the future, SoC should be linked to the mean concentration, but for now it is completely separate
	 */

	bool verb = print && (verbose_gl >= v_crit_gl); // print if the (global) verbose-setting is above the threshold

	double SoCold = getSoC();
	SoC = SoCnew;
	Vcell_valid = false;

	// check the voltage if desired
	if (checkV)
	{
		double V;
		int err;

		// get the voltage
		try
		{
			err = getVcheck(V, print); // throws error if SoC is outside of range
		}
		catch (int e)
		{
			if (verb)
				std::cerr << "Error in Cell::setSoC when getting the voltage, which is " << V << " restoring the old states and throwing on " << e << endl
						  << flush;
			SoC = SoCold;
			throw e;
		}

		// check if the voltage is valid
		if (err == -2 || err == 2)
		{ // safety limits VMIN and VMAX
			SoC = SoCold;
			if (verb)
				std::cerr << "Error in Cell::setSoC the voltage is " << V << " which is outside the safety limits. restoring the old states and throwing 3 " << endl
						  << flush;
			throw 3;
		}
		if (err == -1 || err == 1)
		{ // valid limits Vmin and Vmax
			if (verb)
				std::cerr << "Error in Cell::setSoC the voltage is " << V << " which is outside the allowed limits. throwing 2 " << endl
						  << flush;
			throw 2; // don't restore old SoC
		}

		return V;
	}
	else
	{
		return 0;
	}
}

double Cell_SPM::setI(double Inew, bool checkV, bool print)
{
	/*
	 * sets the current
	 *
	 * checkV	true, the voltage is checked after setting the current
	 * 				if it is outside the safety limits of the cell, error 3 is thrown and the old current is restored
	 * 				if it is outside the valid limits of the cell, error 2 is thrown but the new current is kept
	 * 				if inside allowed Vrange, it returns the voltage
	 * 			false, the voltage is not checked (function returns 0, no errors are thrown)
	 *  		if no value of checkV is given, it is set to true
	 * print 	controls the printing of error messages
	 * 			if true, error messages are printed (if the global printing variable is high enough)
	 * 			if false, no messages are printed, but the errors are still thrown
	 * 			if no value, the default is true
	 *
	 * returns the voltage if checkV = true, else it returns 0
	 *
	 * THROWS
	 * 2 	checkV is true && the voltage is outside the allowed range but still in the safety range
	 * 			and current is in the correct direction. I.e. if charging and V > Vmax or discharging and V < Vmin
	 * 3 	checkV is true && the voltage is outside the safety limits, old current is restored
	 * 			and current is in the correct direction. I.e. if charging and V > VMAX or discharging and V < VMAX
	 * 		if currents are in the 'wrong' direction (e.g. charging but V < Vmin or V < VMIN) then don't throw errors
	 * 			since this current is helping to rectify the situation
	 */

	double Iold = getI();
	I = Inew;
	Vcell_valid = false;
	etacell_valid = false;

	if (checkV)
	{
		bool verb = print && (verbose_gl >= v_crit_gl); // print if the (global) verbose-setting is above the threshold
		double V;
		int err;

		// get the voltage
		try
		{
			err = getVcheck(V, print); // throws error if SoC is outside of range
		}
		catch (int e)
		{
			if (verb)
				std::cerr << "Error in Cell_SPM::setI when getting the voltage, which is " << V << " restoring the old states and throwing on " << e << endl
						  << flush;
			I = Iold;
			throw e;
		}

		// check if the voltage is valid
		if ((err == -2 && I > 0) || (err == 2 && I < 0))
		{ // outside safety limits (< VMIN and discharge || > VMAX and charge)
			I = Iold;
			if (verb)
				std::cerr << "Error in Cell_SPM::setI the voltage is " << V << " which is outside the safety limits when attempting to set I of " << I << ". restoring the old states and throwing 3 " << endl
						  << flush;
			throw 3;
		}
		if ((err == -1 && I > 0) || (err == 1 && I < 0))
		{ // outside valid limits for the cell
			if (verb)
				std::cerr << "Error in Cell_SPM::setI the voltage is " << V << " which is outside the allowed limits when attempting to set I of " << I << ". throwing 2 " << endl
						  << flush;
			throw 2; // do not restore the old current
		}

		return V;
	}
	else
	{
		return 0;
	}
}

void Cell_SPM::setC(double fp0, double fn0)
{
	/*
	 * Function to set the value of the transformed concentrations to the values
	 * corresponding with a uniform (user-specified) lithium concentration.
	 * The cell's current is set to 0 because the concentration is fully uniform (which means the gradient at the surface is 0, so the current must be 0)
	 *
	 * IN
	 * fp0	lithium fraction in the positive electrode 0 <= cp0 <= 1
	 * fn0	lithium fraction in the negative electrode 0 <= cn0 <= 1
	 *
	 * THROWS
	 * 10 	illegal input lithium fractions
	 */

	bool pp = (fp0 < 0) || (fp0 > 1);
	if (pp)
		std::cerr << "ERROR in Cell_SPM::setC, illegal input li-fraction for the positive electrode : " << fp0 << ". The value has to be between 0 and 1" << '\n';
	bool nn = (fn0 < 0) || (fn0 > 1);
	if (nn)
		std::cerr << "ERROR in Cell_SPM::setC, illegal input li-fraction for the negative electrode : " << fn0 << ". The value has to be between 0 and 1" << '\n';
	if (pp || nn)
		throw 10;

	// Calculate the corresponding li-concentrations in [mol m-3]
	double cp = fp0 * Cmaxpos;
	double cn = fn0 * Cmaxneg;

	// Do the first transformation, to u(i) = radius(i) * concentration = x(i) * R * concentration(i)
	double uneg[CELL_NCH], upos[CELL_NCH];
	for (int i = 0; i < CELL_NCH; i++)
	{
		uneg[i] = cn * M.xch[i] * Rn;
		upos[i] = cp * M.xch[i] * Rp;
	}

	// The second transformation is to the eigenspace: z = V * u with V the inverse of the matrix with the eigenvectors.
	// As explained, we know that there is one eigenvector corresponding to a uniform concentration
	// So we need to calculate only this one nonzero value for the (twice) transformed concentration
	// The location of the uniform eigenvector (which has a 0 eigenvalue) is written in M.Input[3]
	int ind = M.Input[3];
	double znu = 0; // (twice) transformed negative uniform concentration
	double zpu = 0; // (twice) transformed positive uniform concentration
	for (int i = 0; i < CELL_NCH; i++)
	{ // loop to calculate the row of V * u corresponding to the uniform concentration
		znu += M.Vn[ind][i] * uneg[i];
		zpu += M.Vp[ind][i] * upos[i];
	}

	// Make the full arrays for the (twice) transformed concentration
	double zp[CELL_NCH], zn[CELL_NCH];
	for (int i = 0; i < CELL_NCH; i++)
	{ // set all values to 0
		zp[i] = 0;
		zn[i] = 0;
	}
	zp[ind] = zpu; // set the non-zero value
	zn[ind] = znu;
	setZp(zp);
	setZn(zn);
	Vcell_valid = false;
}

void Cell_SPM::getStates(double s[], int nin, int &nout)
{
	/*
	 *	returns the states of the cell
	 *
	 * IN
	 * nin 		length of the array provided, must be correct length
	 *
	 * OUT
	 * states	array in which the battery states will be put
	 * 				zp 			the transformed li concentration at the positive inner nodes of the positive particle (CELL_NCH values)
	 * 				zn			the transformed li concentration at the positive inner nodes of the negative particle (CELL_NCH values)
	 * 				delta 		the thickness of the SEI layer [m]
	 * 				LLI 		the lost lithium [As]
	 * 				thickp 		the thickness of the cathode [m]
	 * 				thickn		the thickness of the anode [m]
	 * 				ep 			the volume fraction of active material in the cathode [-]
	 * 				en 			the volume fraction of active material in the anode [-]
	 * 				ap 			the effective surface area of the cathode [m2 m-3]
	 * 				an			the effective surface area of the anode [m2 m-3]
	 * 				CS			the surface area of the cracks at the surface of the negative particle [m2]
	 * 				Dp			the diffusion constant at reference temperature of the cathode [m s-1]
	 * 				Dn			the diffusion constant at reference temperature of the anode [m s-1]
	 * 				rp 			the specific resistance of the cathode [Ohm m2]
	 * 				rn 			the specific resistance of the anode
	 * 				rcc 		the specific resistance of the separator
	 * 				delta_pl 	the thickness of the plated lithium layer [m]
	 * 				SoC			the state of charge of the cell
	 * 				T 			the cell temperature [K]
	 * 				I 			the current of the cell
	 *
	 * throws
	 * 10 invalid array (wrong length)
	 */

	if (nin < getNstates())
	{
		if (verbose_gl >= v_crit_gl)
			std::cerr << "ERROR in Cell_SPM::getStates, array has the wrong length, should be " << 2 * CELL_NCH + 17 + 1 << "  but is " << nin << '\n';
		throw 10;
	}

	// Put the states in the array
	for (int j = 0; j < 2 * CELL_NCH + 15; j++)
		s[j] = state[j];
	s[2 * CELL_NCH + 15] = getSoC();
	s[2 * CELL_NCH + 16] = getT();
	s[2 * CELL_NCH + 17] = getI();

	nout = 2 * CELL_NCH + 17 + 1;
}

void Cell_SPM::getVariations(double var[], int nin, int &nout)
{
	/*
	 * Return the parameters of this cell's variation
	 *
	 * SPM cells have variation in
	 * 	capacity
	 * 	resistance
	 * 	degradation rate
	 */

	if (nin < 4)
	{
		if (verbose_gl >= v_crit_gl)
			std::cerr << "ERROR in Cell_SPM::getVariations, array has the wrong length, should be " << 4 << "  but is " << nin << '\n';
		throw 10;
	}

	var[0] = var_cap;
	var[1] = var_R;
	var[2] = var_degSEI;
	var[3] = var_degLAM;

	nout = 4;
}

double Cell_SPM::getOCV(bool print)
{
	/*
	 * Calculate the OCV of the cell
	 *
	 * THROWS
	 * 6 	the surface concentration is out of bounds
	 */

#if TIMING
	std::clock_t tstart = std::clock();
#endif
	bool verb = print && (verbose_gl >= v_crit_gl); // print if the (global) verbose-setting is above the threshold

	// Calculcate the lithium fractions at the surface of the particles
	double cps = getCp_surface(getZp()) / Cmaxpos;
	double cns = getCn_surface(getZn()) / Cmaxneg;

	// check if the surface concentration is within the allowed range
	// 	0 < cp < 1
	// 	0 < cn < 1
	// don't allow 0 or 1 because in that case, i0 will be 0, and the overpotentials will have 1/0 = inf or nan
	if (cps <= 0 || cps >= 1)
	{
		if (print)
		{ // print error message unless you want to suppress the output
			std::cerr << "ERROR in Cell_SPM::getOCV: concentration out of bounds. the positive lithium fraction is " << cps << " and the negative lithium fraction is " << cns << flush;
			std::cerr << " they should both be between 0 and 1, OCV_n = " << linInt(verb, true, OCV_neg_x, OCV_neg_y, OCV_neg_n, cns) << ", T = " << getT() - 273 << " centigrade" << '\n';
		}
		throw 6;
	}
	if (cns <= 0 || cns >= 1)
	{
		if (print)
		{ // print error message unless you want to suppress the output
			std::cerr << "ERROR in Cell_SPM::getOCV: concentration out of bounds. the positive lithium fraction is " << cps << " and the negative lithium fraction is " << cns << flush;
			std::cerr << " they should both be between 0 and 1. OCV_p = " << linInt(verb, true, OCV_pos_x, OCV_pos_y, OCV_pos_n, cps) << ", T = " << getT() - 273 << " centigrade" << '\n';
		}
		throw 6;
	}

	// Calculate the electrode potentials
	double OCV_p;	   // cathode potential [V]
	double OCV_n;	   // anode potential [V]
	double dOCV;	   // entropic coefficient of the total cell voltage [V/K]
	bool bound = true; // in linear interpolation, throw an error if you are out of the allowed range
	try
	{
		dOCV = linInt(verb, bound, dOCV_tot_x, dOCV_tot_y, dOCV_tot_n, cps);
		OCV_n = linInt(verb, bound, OCV_neg_x, OCV_neg_y, OCV_neg_n, cns);
		OCV_p = linInt(verb, bound, OCV_pos_x, OCV_pos_y, OCV_pos_n, cps);
	}
	catch (int e)
	{
		if (print)
			cout << "error in Cell_SPM::getOCV when getting the electrode potentials " << e << ". Throwing it up" << endl
				 << flush;
		throw e;
	}

#if TIMING
	T_getOCV += (std::clock() - tstart) / (double)CLOCKS_PER_SEC; // time in seconds
#endif
	return OCV_p - OCV_n + (getT() - T_ref) * dOCV;
}
double Cell_SPM::getOCV(double fps, double fns, bool print)
{
/*
 * Calculate the OCV of the cell
 *
 * IN
 * fps 	lithium FRACTION at the surface of the cathode
 * fns 	lihtium FRACTION at the surface of the anode
 * 		must be between 0 and 1
 *
 * THROWS
 * 6 	the surface concentration is out of bounds
 */
#if TIMING
	std::clock_t tstart = std::clock();
#endif
	bool verb = print && (verbose_gl >= v_crit_gl); // print if the (global) verbose-setting is above the threshold

	// Calculate the electrode potentials
	double OCV_p; // cathode potential [V]
	double OCV_n; // anode potential [V]
	// double dOCV;				// entropic coefficient of the total cell voltage [V/K]
	bool bound = true; // in linear interpolation, throw an error if you are out of the allowed range
	try
	{
		// dOCV = linInt(verb, bound, dOCV_tot_x, dOCV_tot_y, dOCV_tot_n, fps);
		OCV_n = linInt(verb, bound, OCV_neg_x, OCV_neg_y, OCV_neg_n, fns);
		OCV_p = linInt(verb, bound, OCV_pos_x, OCV_pos_y, OCV_pos_n, fps);
	}
	catch (int e)
	{
		if (print)
			cout << "error in Cell_SPM::getOCV when getting the electrode potentials " << e << ". Throwing it up" << endl
				 << flush;
		throw e;
	}

#if TIMING
	T_getOCV += (std::clock() - tstart) / (double)CLOCKS_PER_SEC; // time in seconds
#endif
	return OCV_p - OCV_n; //+ (getT()-T_ref)*dOCV;
}
double Cell_SPM::getV(bool print)
{
/*
 * Function to calculate the cell voltage
 *
 * IN
 * print 	boolean indicating if we want to print error messages or not
 * 				if true, error messages are printed
 * 				if false no error messages are printed (but the error will still be thrown)
 * 			we need this input from higher level functions because at this point we cannot know if an error will be critical or not
 *
 * OUT
 * V 		battery voltage [V]
 * bool 	indicates if the voltage is in the allowed range Vmin <= V <= Vmax
 *
 * THROWS
 */
#if TIMING
	std::clock_t tstart = std::clock();
#endif
	bool verb = print && (verbose_gl >= v_crit_gl); // print if the (global) verbose-setting is above the threshold

	// If the stored value is the most up to date one, then simply return this value
	if (Vcell_valid)
		return Vcell;

	// Calculcate the lithium fractions at the surface of the particles
	double cps = getCp_surface(getZp());
	double cns = getCn_surface(getZn());
	double fps = cps / Cmaxpos;
	double fns = cns / Cmaxneg;

	// check if the surface concentration is within the allowed range
	// 	0 < cp < 1
	// 	0 < cn < 1
	// don't allow 0 or 1 because in that case, i0 will be 0, and the overpotentials will have 1/0 = inf or nan
	if (fps <= 0 || fps >= 1)
	{
		if (print)
		{ // print error message unless you want to suppress the output
			std::cerr << "ERROR in Cell_SPM::getV for cell " << getFullID() << " concentration out of bounds. the positive lithium fraction is " << fps << " and the negative lithium fraction is " << fns << flush;
			std::cerr << " they should both be between 0 and 1, OCV_n = " << linInt(verb, true, OCV_neg_x, OCV_neg_y, OCV_neg_n, fns) << ", T = " << getT() - 273 << " centigrade, I = " << getI() << '\n';
		}
		throw 6;
	}
	if (fns <= 0 || fns >= 1)
	{
		if (print)
		{ // print error message unless you want to suppress the output
			std::cerr << "ERROR in Cell_SPM::getV for cell " << getFullID() << " concentration out of bounds. the positive lithium fraction is " << fps << " and the negative lithium fraction is " << fns << flush;
			std::cerr << " they should both be between 0 and 1. OCV_p = " << linInt(verb, true, OCV_pos_x, OCV_pos_y, OCV_pos_n, fps) << ", T = " << getT() - 273 << " centigrade, I = " << getI() << '\n';
		}
		throw 6;
	}

	// get the OCV
	double OCV;
	try
	{
		OCV = getOCV(fps, fns, print);
		// OCV = getOCV(print);
	}
	catch (int e)
	{
		if (verb)
			cout << "error in Cell_SPM when getting the open circuit voltage, passing it on" << endl
				 << flush;
		throw e;
	}

	// Calculate the overpotential if needed
	if (!etacell_valid)
	{
		// Calculate the rate constants at the cell's temperature using an Arrhenius relation
		double kpt = kp * exp(kp_T / Rg * (1 / T_ref - 1 / getT()));
		double knt = kn * exp(kn_T / Rg * (1 / T_ref - 1 / getT()));

		// Calculate the overpotential using the Bulter-Volmer equation
		// 		if alpha is 0.5, the Bulter-Volmer relation can be inverted to eta = 2RT / (nF) asinh(x)
		double i_app = getI() / elec_surf;															  // current density on the electrodes [I m-2]
		double i0p = kpt * n * F * C_elec_sqrt * sqrt(cps) * sqrt(Cmaxpos - cps);					  // exchange current density of the positive electrode
		double i0n = knt * n * F * C_elec_sqrt * sqrt(cns) * sqrt(Cmaxneg - cns);					  // exchange current density of the negative electrode
		etapcell = (2 * Rg * getT()) / (n * F) * asinh(-0.5 * i_app / (getAp() * getThickp()) / i0p); // cathode overpotential [V], < 0 on discharge
		etancell = (2 * Rg * getT()) / (n * F) * asinh(0.5 * i_app / (getAn() * getThickn()) / i0n);  // anode overpotential [V],  > 0 on discharge
		etacell_valid = true;
	}

	// Calculate the cell voltage
	Vcell = OCV + (etapcell - etancell) - getRdc() * getI();
	Vcell_valid = true; // we now have the most up to date value stored

#if TIMING
	T_getV += (std::clock() - tstart) / (double)CLOCKS_PER_SEC; // time in seconds
#endif

	return Vcell;
}

double Cell_SPM::setStates(double s[], int nin, bool checkV, bool print)
{
/*
 *
 * THROWS
 * 2 	checkV is true && the voltage is outside the allowed range
 * 10 	illegal state proposed, or not all fields are present
 * 		old state is restored
 */
#if TIMING
	std::clock_t tstart = std::clock();
#endif
	bool verb = print && (verbose_gl >= v_crit_gl); // print if the (global) verbose-setting is above the threshold

	if (nin < getNstates())
	{
		if (verb)
			std::cerr << "ERROR in Cell_SPM::setStates, array has the wrong length, should be " << 2 * CELL_NCH + 18 << "  but is " << nin << '\n';
		throw 10;
	}

	// if we check the voltage, first check the states are valid
	// else you might get weird errors in getV() if the states are illega;
	if (checkV)
	{
		bool val = validStates(s, nin, print);
		if (!val)
		{
			if (verb)
				cout << "Error in Cell_SPM::setStates, the state vector is not allowed as determined by validStates" << endl
					 << flush;
			throw 10;
		}
	}

	// Use the setters of every property
	setZp(&s[0]);
	setZn(&s[CELL_NCH]);
	setDelta(s[2 * CELL_NCH + 0]);
	setLLI(s[2 * CELL_NCH + 1]);
	setThickp(s[2 * CELL_NCH + 2]);
	setThickn(s[2 * CELL_NCH + 3]);
	setEp(s[2 * CELL_NCH + 4]);
	setEn(s[2 * CELL_NCH + 5]);
	setAp(s[2 * CELL_NCH + 6]);
	setAn(s[2 * CELL_NCH + 7]);
	setCS(s[2 * CELL_NCH + 8]);
	setDp(s[2 * CELL_NCH + 9]);
	setDn(s[2 * CELL_NCH + 10]);
	setrdcp(s[2 * CELL_NCH + 11]);
	setrdcn(s[2 * CELL_NCH + 12]);
	setrdccc(s[2 * CELL_NCH + 13]);
	setDelta_pl(s[2 * CELL_NCH + 14]);
	setSoC(s[2 * CELL_NCH + 15], checkV, print);
	setT(s[2 * CELL_NCH + 16]);
	setI(s[2 * CELL_NCH + 17], checkV, print);
	Vcell_valid = false;

	// check the voltage if desired
	if (checkV)
	{
		double V;
		try
		{
			V = getV(verb); // throws error if SoC is outside of range, but that should have been checked in validState()

#if TIMING
			T_setStates += (std::clock() - tstart) / (double)CLOCKS_PER_SEC; // time in seconds
#endif

			if (V < getVmin() || V > getVmax())
				throw 2; // throws error if SoC is valid, but V is not
			return V;
		}
		catch (int e)
		{
			if (verb)
				std::cerr << "ERROR in Cell_SPM::setStates when getting the voltage, which is " << V << endl
						  << flush;
			throw 2;
		}
	}
	else
	{
#if TIMING
		T_setStates += (std::clock() - tstart) / (double)CLOCKS_PER_SEC; // time in seconds
#endif
		return 0;
	}
}

bool Cell_SPM::validStates(double s[], int nin, bool print)
{
/*
 * Check the array contains states which are valid for an SPM
 * This function checks
 * 		the length of the array
 * 		0 < SoC < 1, but does not check SoC is compatible with the concentration
 * 		//todo calculate SoC based on the uniform concentration
 * 		Tmin < T < Tmax
 * 		0 < surface concentration < maximum concentration
 * 		a = 3e/R
 * 		all other values (except I) > 0
 */
#if TIMING
	std::clock_t tstart = std::clock();
#endif
	double tol = pow(10, -5);

	if (nin < getNstates())
	{
		std::cerr << "ERROR in Cell_SPM::validStates, array has the wrong length,  should be at least " << 2 * CELL_NCH + 18 << "  but is " << nin << '\n';
		throw 10;
	}

	bool verb = print && verbose_gl >= v_crit_gl; // print if the (global) verbose-setting is above the threshold

	bool range = true; // are all in the allowed range?

	// check whether SoC in the allowed range
	if (s[2 * CELL_NCH + 15] < 0.0 || s[2 * CELL_NCH + 15] > 1.0)
	{
		if (verb)
			std::cerr << "ERROR in Cell_SPM::validState, SoC is outside of the range, 0 <= SoC <= 1: " << s[2 * CELL_NCH + 15] << '\n';
		range = false;
	}

	// check whether T is in the allowed range
	if (s[2 * CELL_NCH + 16] < Tmin || s[2 * CELL_NCH + 16] > Tmax)
	{
		if (verb)
			std::cerr << "ERROR in Cell_SPM::validState, T is outside of the range, " << Tmin << " <= T <= " << Tmax << ", value is " << s[2 * CELL_NCH + 16] << '\n';
		range = false;
	}

	// Check the surface concentration
	double zpo[CELL_NCH], zno[CELL_NCH];
	for (int i = 0; i < CELL_NCH; i++)
	{
		zpo[i] = s[i];
		zno[i] = s[CELL_NCH + i];
	}
	double cps = getCp_surface(zpo);
	double cns = getCp_surface(zno);
	if (cps <= 0 || cns <= 0 || cps >= Cmaxpos || cns >= Cmaxneg)
	{
		if (verb)
		{ // print error message unless you want to suppress the output
			std::cerr << "ERROR in Cell_SPM::validState: concentration out of bounds. the positive lithium fraction is " << cps / Cmaxpos << " and the negative lithium fraction is " << cns / Cmaxneg;
			std::cerr << " they should both be between 0 and 1" << '\n';
		}
		range = false;
	}

	// check a = 3*e/R
	double ep = s[2 * CELL_NCH + 4];
	double en = s[2 * CELL_NCH + 5];
	double ap = s[2 * CELL_NCH + 6];
	double an = s[2 * CELL_NCH + 7];
	double app = 3 * ep / Rp;
	double ann = 3 * en / Rn;
	if (abs(ap - app) / ap > tol)
	{
		if (verb)
			std::cerr << "ERROR in Cell_SPM::validState. The value of ap is " << ap << " which is different from 3 ep/Rp, which has value " << app << " giving a relative error of " << abs(ap - app) / ap << '\n';
		range = false;
	}
	if (abs(an - ann) / an > tol)
	{
		if (verb)
			std::cerr << "ERROR in Cell_SPM::validState. The value of an is " << an << " which is different from 3 en/Rn, which has value " << ann << " giving a relative error of " << abs(an - ann) / an << '\n';
		range = false;
	}

	// all other properties must have a positive value
	for (int i = 0; i <= 16; i++)
	{
		if (s[2 * CELL_NCH + i] < 0)
		{
			if (verb)
				std::cerr << "ERROR in Cell_SPM::validState, one of the geometric states has a negative value, its value is " << s[2 * CELL_NCH + i] << endl
						  << flush;
			range = false;
		}
	}
	/* s[] consists of
	 * [0 CELL_NCH-1]			zp 			the transformed li concentration at the positive inner nodes of the positive particle (CELL_NCH values)
	 * [CELL_NCH 2*CELL_NCH-1]	zn			the transformed li concentration at the positive inner nodes of the negative particle (CELL_NCH values)
	 * 2*CELL_NCH + 0			delta 		the thickness of the SEI layer [m]
	 * 2*CELL_NCH + 1			LLI 		the lost lithium [As]
	 * 2*CELL_NCH + 2			thickp 		the thickness of the cathode [m]
	 * 2*CELL_NCH + 3			thickn		the thickness of the anode [m]
	 * 2*CELL_NCH + 4			ep 			the volume fraction of active material in the cathode [-]
	 * 2*CELL_NCH + 5			en 			the volume fraction of active material in the anode [-]
	 * 2*CELL_NCH + 6			ap 			the effective surface area of the cathode [m2 m-3]
	 * 2*CELL_NCH + 7			an			the effective surface area of the anode [m2 m-3]
	 * 2*CELL_NCH + 8			CS			the surface area of the cracks at the surface of the negative particle [m2]
	 * 2*CELL_NCH + 9			Dp			the diffusion constant at reference temperature of the cathode [m s-1]
	 * 2*CELL_NCH + 10			Dn			the diffusion constant at reference temperature of the anode [m s-1]
	 * 2*CELL_NCH + 11			rp 			the specific resistance of the cathode [Ohm m2]
	 * 2*CELL_NCH + 12			rn 			the specific resistance of the anode
	 * 2*CELL_NCH + 13			rcc 		the specific resistance of the separator
	 * 2*CELL_NCH + 14			delta_pl 	the thickness of the plated lithium layer [m]
	 * 2*CELL_NCH + 15			SoC 		the state of charge [-]
	 * 2*CELL_NCH + 16			T 			the temperature [K]
	 * 2*CELL_NCH + 17			I 			the current [A], negative for charging
	 */

#if TIMING
	T_validStates += (std::clock() - tstart) / (double)CLOCKS_PER_SEC; // time in seconds
#endif

	return range;
}

void Cell_SPM::SEI(double OCVnt, double etan, double *isei, double *den)
{
	/*
	 * Function to calculate the degradation effects of growth of the SEI layer
	 *
	 * IN
	 * OCVnt 	the OCV of the negative electrode at the battery temperature [V]
	 * etan 	the overpotential at the negative electrode [V]
	 *
	 * OUT
	 * isei 	current density for the SEI side-reaction [A m-2]
	 * den 		decrease in the active volume fraction as a result of SEI growth [sec-1]
	 *
	 * THROWS
	 * 106 		illegal value in id or por
	 * 107		too many degradation models
	 */

	// variables
	double kseit, Dseit;		// temperature-dependent SEI parameters, using Arrhenius' relation
	double isei1, isei2, isei3; // temporary parameters to calculate the SEI growth
	double is = 0;				// SEI side reaction current density of all models combined

	// Loop for each model to use
	for (int i = 0; i < deg_id.SEI_n; i++)
	{

		// Use a switch to calculate the magnitude of the SEI growth according to this degradation model
		switch (deg_id.SEI_id[i])
		{
		case 0: // no SEI growth
			is += 0;
			break;
		case 1:																														   // Kinetic model according to Ning & Popov, Journal of the Electrochemical Society 151 (10), 2004
			kseit = seiparam.sei1k * exp(seiparam.sei1k_T / Rg * (1 / T_ref - 1 / getT()));											   // Arrhenius relation for the rate parameter at the cell temperature
			is += nsei * F * kseit * exp(-nsei * F / (Rg * getT()) * alphasei * (OCVnt + etan - OCVsei + rsei * getDelta() * getI())); // Add the effect of this model
																																	   // eta_sei = OCVneg + etaneg - OCVsei + rsei*I
																																	   // isei = nFk exp(-nF/RT alpha eta_sei)
																																	   // on charge, I < 0 and etan < 0.
																																	   // so higher charging current -> more negative term in exponential -> larger isei
			break;
		case 2:																				// kinetics and diffusion according to Pinson & Bazant, Journal of the Electrochemical society 160 (2), 2013
			kseit = seiparam.sei2k * exp(seiparam.sei2k_T / Rg * (1 / T_ref - 1 / getT())); // Arrhenius relation for the rate parameter at the cell temperature
			Dseit = seiparam.sei2D * exp(seiparam.sei2D_T / Rg * (1 / T_ref - 1 / getT())); // Arrhenius relation for the diffusion constant at the cell temperature

			// derivation of the formula:
			// start equations (use the symbols from Yang, Leng, Zhang, Ge, Wang, Journal of Power Sources 360, 2017
			// but with opposite sign for j
			// j = nFk c exp(..)
			// j/nF = - D/delta (c - c0)
			// j/nF = D/delta (- j/(nFk exp(..)) + c0)
			// j * ( 1/(nFk exp(..)) + delta/DnF ) = c0
			// j = c0 / ( 1/(nFk exp(..)) + delta/DnF )
			isei2 = nsei * F * kseit * exp(-nsei * F / (Rg * getT()) * alphasei * (OCVnt + etan - OCVsei + rsei * getDelta() * getI()));
			isei3 = getDelta() / (nsei * F * Dseit);
			is += c_elec0 / (1 / isei2 + isei3); // Add the effects of this model
			break;
		case 3:																				// model from Christensen & Newmann, Journal of the Electrochemical Society 152 (4), 2005
			kseit = seiparam.sei3k * exp(seiparam.sei3k_T / Rg * (1 / T_ref - 1 / getT())); // Arrhenius relation for the rate parameter at the cell temperature
			Dseit = seiparam.sei3D * exp(seiparam.sei3D_T / Rg * (1 / T_ref - 1 / getT())); // Arrhenius relation for the diffusion constant at the cell temperature

			// Use equation [22] from the paper
			isei1 = 0.134461 * exp(-nsei * F * (etan + rsei * getDelta() * getI()) / (Rg * getT())); // the parameter a_L_K is set to 0.134461 but this constant can be lumped into the rate- and diffusion constants
			isei2 = nsei * F * kseit * exp(-nsei * F / (Rg * getT()) * alphasei * (OCVnt - OCVsei));
			isei3 = getDelta() / (nsei * F * Dseit);
			is += isei1 / (1 / isei2 + isei3); // Add the effects of this model
			break;
		case 4:																				// model from the optimisation in the paper
			kseit = seiparam.sei4k * exp(seiparam.sei4k_T / Rg * (1 / T_ref - 1 / getT())); // Arrhenius relation for the rate parameter at the cell temperature
			Dseit = seiparam.sei4D * exp(seiparam.sei4D_T / Rg * (1 / T_ref - 1 / getT())); // Arrhenius relation for the diffusion constant at the cell temperature

			isei1 = 0.134461 * exp(-nsei * F * (etan) / (Rg * getT())); // note: model used in optimisation had kpt/knt or vice versa, here a fixed value
			isei2 = nsei * F * kseit * exp(-nsei * F / (Rg * getT()) * alphasei * (OCVnt - OCVsei));
			isei3 = getDelta() / (nsei * F * Dseit);
			is += isei1 / (1 / isei2 + isei3);
			break;
		default: // unknown degradation model
			std::cerr << "ERROR in Cell::SEI, unknown SEI degradation model with identifier " << deg_id.SEI_id[i] << ". Only values 0 to 3 are allowed. Throw an error" << '\n';
			throw 106;
			break;
		}
	} // end loop for all the models you want to use

	// Make the output for the SEI side reaction current density
	*isei = is;

	// Calculate how much we decrease the volume fraction due to SEI growth
	if (deg_id.SEI_porosity == 0) // don't decrease volume fraction
		*den = 0;
	else if (deg_id.SEI_porosity == 1)
	{																	  // decrease volume fraction according to Ashwin, Chung, Wang, Journal of Power Sources 328, 2016
		double jn = getI() / elec_surf / (getAn() * n * F * getThickn()); // molar flux on the negative particle
		*den = -seiparam.sei_porosity * (jn * Vmain + *isei * Vsei);
		// note: they use J = volumetric current [A/m3] -> they multiply with 'an' but we already have density [A/m2]
		// - because they use the porosity while we use the volume fraction
	}
	else
	{ // unknown degradation model
		std::cerr << "ERROR in Cell::SEI, unknown value for decreasing the volume fraction " << deg_id.SEI_porosity << ". Only values 0 or 1 are allowed. Throw an error" << '\n';
		throw 106;
	}
}

void Cell_SPM::CS(double OCVnt, double etan, double *isei_multiplyer, double *dCS, double *dDn)
{
	/*
	 * function to calculate the degradation effect of surface cracking due to fatigue
	 *
	 * IN
	 * OCVnt 			the OCV of the negative electrode at the battery temperature [V]
	 * etan 			the overpotential at the negative electrode [V]
	 *
	 * OUT
	 * isei_multiplyer 	Extra SEI side reaction due to crack growth as a fraction of the original SEI side reaction current [-].
	 * 						i.e. total SEI growth = (1+isei_multiplyer)*isei
	 * dCS 				increase in surface area due to fatigue in this time step [m2 sec-1]
	 * dDn 				decrease in the diffusion constant of the graphite due to surface cracks [m s-1 s-1]
	 *
	 * THROWS
	 * 106 				illegal value in id or d
	 * 107				too many degradation models
	 * 108 				the stress values are not up to date
	 */

	// output parameters
	double ism = 0; // isei multiplier from all models to be considered
	double dcs = 0; // increase in surface area from all models to be considered

	double ASn = getAn() * getThickn() * elec_surf; // active surface area of the anode [m2]
													// this active area is used to translate absolute values (such as currents) to relative values (such as current densities)

	// Loop for each model we want to use
	for (int i = 0; i < deg_id.CS_n; i++)
	{

		// a switch to calculate the effect according to model i
		switch (deg_id.CS_id[i])
		{
		case 0: // no surface cracks
			dcs += 0;
			break;
		case 1: // Laresgoiti's stress and crack growth model (Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015)
				// this model calculates crack growth due to temporal variations in the li-concentration
			// check the calculated stress values are up to date
			if (!s_lares_update)
			{
				std::cerr << "ERROR in Cell::CS. The stress values for Laresgoiti's stress model are not updated. Throwing an error" << endl
						  << flush;
				// if you see this error, you have to call Cell::updateLaresgoitiStress(), which calculates the stress and stores the values, before you call this function
				throw 108;
			}

			// Implement the equation from the paper
			// capacity loss is m-power of the so-called stress amplitude (sigma_max - sigma_min)/2
			// sigma_max and sigma_min are the max and min stresses 'of the cyclic signal' i.e. within one charge/discharge
			// assume m = 1, then (max - min) = (max - t1) + (t1-t2) + (t2-t3) + ... + (tn - min)
			// so stress amplitude can be substituted by (the stress in the previous time step) - (the stress in this time step)
			dcs += csparam.CS1alpha * pow(abs(s_lares_n - s_lares_n_prev) / s_dt, 0.5); // equations (22)+ (27) from the paper. for s_dt, see LAM function with Dai lam
			break;
		case 2: // Laresgoiti's crack growth model (Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015)
				// but with stress model from Dai, Cai, White, Journal of Power sources 247, 2014
				// instead of Laresgoiti's stress from figure 5
				// this model calculates crack growth due to spatial (Dai) and temporal (Laresgoiti) variations in the li-concentration
			// ensure the stress values are up to date
			if (!s_dai_update)
			{
				std::cerr << "ERROR in Cell::CS. The stress values for Dai's stress model are not updated. Throwing an error" << endl
						  << flush;
				// if you see this error, you have to call Cell::updateDaiStress(), which calculates the stress and stores the values before you call this function
				throw 108;
			}

			// Add the effects of this model
			dcs += csparam.CS2alpha * pow(abs(s_dai_n - s_dai_n_prev) / s_dt, 0.5); // equations (22)+ (27) from the paper but with Dai's, for s_dt see Dai-stress in LAM function
			break;
		case 3: // model by Deshpande & Bernardi,Journal of the Electrochemical Society 164 (2), 2017
				// this model is adapted to calculate crack growth due to spatial variations in the li-concentration
			// get concentrations
			double cp[CELL_NCH + 2], cn[CELL_NCH + 2];
			getConcentration(CELL_NCH + 2, cp, cn);

			// Add the effects of this model
			dcs += csparam.CS3alpha * pow((cn[0] - cn[CELL_NCH + 1]) / Cmaxneg, 2); // equations (8) + (21)
																					// Note that eqn (8) refers to the change with respect to time while here we use the spatial variation
																					// This makes the model capture more of the spatial variation in stress
																					// Laresgoiti's model already accounted for temporal variation, so simply use that if you are interested in temporal rather than spatial variation
																					//	ism += s.getCS()/ASn; 											// increase SEI growth proportionally the crack surface
			break;
		case 4: // model from Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
			// equation (1a): CS = Amax(1 - exp(-m * Ah)) with m a fitting parameter and Ah the charge throughput up to now
			// 	dCS/dAh = m Amax exp(-m Ah) = m Amax - m Amax + m Amax exp(-m Ah) = m Amax - m CS = m(Amax - CS)
			// 	dCS/dt = dCS/dAh * dAh/dt = dCS/dAh * abs(I) = m(Amax - CS)*abs(I)
			// 		where we use the absolute value of I because 'Ah' is the total charge throughput, i.e. int ( abs(I) dt )

			double Amax;						  // 'maximum crack surface area', a fitting parameters
			Amax = max(csparam.CS4Amax, getCS()); // avoid negative crack growth if the crack surface becomes slightly larger than Amax
												  // this is possible due to discrete time steps: CS(t) is just smaller, but CS (t+1) = CS(t) + dCS*dt is just larger

			// Add the effects of this model
			dcs += csparam.CS4alpha * (Amax - getCS()) * abs(getI()); // see above, with m = csparam.CS4
			break;
		case 5:			   // model from Ekstrom and Lindbergh, Journal of the Electrochemical Society 162 (6), 2015
			double etasei; // overpotential for the crack-side-reaction = overpotential for the SEI reaction
			double kcr;	   // rate constant for the side reaction

			etasei = (OCVnt + etan - OCVsei + rsei * getDelta() * getI()); // overpotential [V], equation (6)

			// get surface concentration
			double cns;
			cns = getCn_surface(getZn());

			// Calculate the rate constant, equation (11) with an Arrhenius relation for the temperature (which wasn't considered by Ekstrom)
			if (getI() > 0)
				kcr = 0;
			else if (cns / Cmaxneg < 0.3)
				kcr = 2 * csparam.CS5k * exp(csparam.CS5k_T / Rg * (1 / T_ref - 1 / getT()));
			else if (cns / Cmaxneg < 0.7)
				kcr = 0;
			else
				kcr = csparam.CS5k * exp(csparam.CS5k_T / Rg * (1 / T_ref - 1 / getT()));

			// Add the effects of this model
			dcs += nsei * F * kcr * exp(-alphasei * nsei * F / (Rg * getT()) * etasei); // equation (9)
			break;
		default: // unknown degradation model
			std::cerr << "ERROR in Cell::CS, unknown crack growth model with identifier " << deg_id.CS_id[i] << ". Only values 0 to 5 are allowed. Throw an error" << '\n';
			throw 106;
			break;
		}
	}

	// increase SEI growth proportionally the crack surface
	ism += getCS() / ASn;

	// Make the output variables
	*isei_multiplyer = ism;
	*dCS = dcs;

	// Decrease the negative diffusion constant if needed
	if (deg_id.CS_diffusion == 0) // don't decrease the negative diffusion constant
		*dDn = 0;
	else if (deg_id.CS_diffusion == 1)
	{ // decrease it according to Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
		// equation (2) D(t) = D0 (1 - CS)^gamma
		// 	but this can become negative is CS is larger than 1, we assume there should be a term /Amax in eqn (2): D(t) = D0 (1 - (CS/Amax))^gamma, which becomes 0 if CS = Amax
		// so dD/dt = - gamma D0 (1-CS/Amax)^(gamma-1) 1/Amax dCS/dt

		double Dnmax;						  // cap the decrease rate at a maximum value of 2e-7 (2e-5% = kill the battery in  about 1000h)
		double Amax;						  // 'maximum crack surface area', a fitting parameters
		Amax = max(csparam.CS4Amax, getCS()); // avoid increasing diffusion coefficient if the crack surface becomes larger than Amax
											  // this is possible if the user chooses a different CS growth model, which does give larger crack surfaces (i.e. not CS4 which is Barai's crack growth model)
		Dnmax = csparam.CS_diffusion * pow(1 - getCS() / Amax, csparam.CS_diffusion - 1) / Amax * dcs;
		Dnmax = min(Dnmax, 2 * pow(10, -7));
		*dDn = -Dnmax * getDn();
	}
	else
	{ // unknown degradation model
		std::cerr << "ERROR in Cell::CS, unknown value for decreasing the diffusion constant " << deg_id.CS_diffusion << ". Only values 0 or 1 are allowed. Throw an error" << '\n';
		throw 106;
	}
}

void Cell_SPM::LAM(bool print, double zp_surf, double etap,
				   double *dthickp, double *dthickn, double *dap, double *dan, double *dep, double *den)
{
	/*
	 * Function to calculate the effect of loss of active material (LAM).
	 *
	 * LAM is simulated by decreasing the amount of active material.
	 * This will increase the current density on the particle for the same overall cell current (because there is less material to 'spread' it over).
	 * This will mean that for the same overall cell current,
	 * 		there will be a larger change in lithium concentration,
	 * 			so there is a larger change in open circuit voltage,
	 * 			so a smaller capacity before a voltage limit is reached
	 * 			so the capacity decreases,
	 * 		and there will be a larger resistive voltage drop, so the 'effective' resistance increases
	 *
	 * There are a couple of variables describing the amount of active material:
	 * 	elec_surf 	the (geometric) surface of the electrode, i.e. the product of the height and length of the electrode [m2]
	 * 	thick		the thickness of the electrode material [m]
	 * 	R 			radius of the particle of the single particle model [m]
	 * 	e			the volume fraction of active material	[-]
	 * 	a 			the effective surface, i.e. the surface per unit of electrode volume [m2 m-3]
	 * 				a = 3*e/R
	 *
	 * 	The current density is calculated as:
	 * 	i = I / (a * thick * elec_surf) = I / (3* e / R * thick * elec_surf)
	 *
	 * 	A decrease in any of these geometric parameters will have exactly the same effect on the battery:
	 * 	you can double i by halving a, or by halving thick, or by halving elec_surf.
	 * 	and you can halve a by halving e or doubling R.
	 * 	So it doesn't really matter which of the geometric parameters you decrease to account for LAM.
	 * 	However, this model doesn't allow to change elec_surf and R because these values are needed on multiple locations in the code.
	 * 	E.g. the value of R is needed to calculate values of the matrices used for the diffusion state space model (because the spatial discretisation depends on R).
	 * 	So if you were to change R, you have to re-calculate the matrices, which would needlessly complicate the model.
	 * 	Therefore, degradation can only decrease the values of 'thick', 'a' and 'e'
	 *
	 * IN
	 * print 	boolean indicating if we want to print error messages or not
	 * 				if true, error messages are printed
	 * 				if false no error messages are printed (but the error will still be thrown)
	 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
	 * zp_surf 		li-fraction at the surface of the positive particle [-]
	 * etap 		overpotential at the positive electrode [V]
	 *
	 * OUT
	 * dthickp 		change in electrode thickness of the positive electrode [m s-1]
	 * dthickn		change in electrode thickness of the negative electrode [m s-1]
	 * dap			change in effective electrode surface of the positive electrode [m2 m-3 s-1]
	 * dan			change in effective electrode surface of the negative electrode [m2 m-3 s-1]
	 * dep			change in volume fraction of active material in the positive electrode [s-1]
	 * den			change in volume fraction of active material in the negative electrode [s-1]
	 *
	 * THROWS
	 * 106 			illegal value in id
	 * 107			too many degradation models
	 * 108 			the stress values are not updated
	 */

	// output parameters
	double dthickpp = 0;
	double dthicknn = 0;
	double dapp = 0;
	double dann = 0;
	double depp = 0;
	double denn = 0;

	// loop for each model to use
	for (int i = 0; i < deg_id.LAM_n; i++)
	{

		// calculate the effect of this model
		switch (deg_id.LAM_id[i])
		{
		case 0: // no LAM
			dthickpp += 0;
			dthicknn += 0;
			dapp += 0;
			dann += 0;
			depp += 0;
			denn += 0;
			break;
		case 1: // Stress model from Dai, Cai, White, Journal of Power sources 247, 2014
				// LAM equation similar to CS equation from Laresgoiti
				// (Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015)
			// ensure the stress values are up to date
			if (!s_dai_update)
			{
				std::cerr << "ERROR in Cell::LAM. The stress values for Dai's stress model are not updated. Throwing an error" << endl
						  << flush;
				// if you see this error, you have to call Cell::updateDaiStress(), which calculates the stress and stores the values before calling this function
				throw 108;
			}

			// Laresgoiti's equation to link stress to LAM
			dthickpp += -lamparam.lam1p * pow(abs(s_dai_p - s_dai_p_prev) / s_dt, 1.0); // with Dai's stress model (values stored in s_dai_p)
			dthicknn += -lamparam.lam1n * pow(abs(s_dai_n - s_dai_n_prev) / s_dt, 1.0); // ageing fit
																						// you need to divide by the time step to counter the effect of the time step
																						// 	larger dt has double effect: increase the stress difference, and time integrate by larger time period
																						// 	assuming stress changes are constant at ds per second, it would be ds (time_now - time_prev), or ds * dt
																						// 	and to cover a period of T (so we need T/dt steps of dt each), the total effect of stress is (ds*dt) * dt * T/dt = ds*dt*T
																						// 	so the effect of stress increases linearly with the size of the time setp
																						// 	so divide by total time (dt*nstep) to cancel this out, i.e. ds is per second (and not per total time period)
			// assume the other effects are 0
			dapp += 0;
			dann += 0;
			depp += 0;
			denn += 0;
			break;
		case 2: // Model by Delacourt & Safari, Journal of the Electrochemical Society 159 (8), 2012
			// Get the molar flux on each particle
			double i_app, jp, jn;
			i_app = getI() / elec_surf;					   // current density on the electrode [A m-2]
			jp = -i_app / (getAp() * n * F * getThickp()); // molar flux on the positive particle [mol m-2 s-1]
			jn = i_app / (getAn() * n * F * getThickn());  // molar flux on the negative particle [mol m-2 s-1]

			// Use Arrhenius relations to update the fitting parameters for the cell temperature
			double ap, bp, an, bn;
			ap = lamparam.lam2ap * exp(lamparam.lam2t / Rg * (1 / T_ref - 1 / getT()));
			an = lamparam.lam2an * exp(lamparam.lam2t / Rg * (1 / T_ref - 1 / getT()));
			bp = lamparam.lam2bp * exp(lamparam.lam2t / Rg * (1 / T_ref - 1 / getT()));
			bn = lamparam.lam2bn * exp(lamparam.lam2t / Rg * (1 / T_ref - 1 / getT()));

			// Add the effects of this model
			depp += ap * abs(jp) + bp * sqrt(abs(jp)); // equation (5) from the paper
			denn += an * abs(jn) + bn * sqrt(abs(jn));
			// assume the other effects are 0
			dthickpp += 0;
			dthicknn += 0;
			dapp += 0;
			dann += 0;
			break;
		case 3:				 // Model by Kindermann, Keil, Frank, Jossen, Journal of the Electrochemical Society 164 (12), 2017
			double OCVpt;	 // cathode potential
			double etap_LAM; // overpotential for the NMC dissolution reaction
			double kt;		 // temperature dependent rate constant
			double idiss;	 // current density of the NMC dissolution reaction

			try
			{
				OCVpt = linInt(print, true, OCV_pos_x, OCV_pos_y, OCV_pos_n, zp_surf); // get OCV of positive electrode, throw error if out of bounds
																					   // this should be updated for the cell's temperature using the entropic coefficient of the cathode
																					   // but I couldn't find any data on this, so I have ignored the effect
			}
			catch (int e)
			{
				cout << "Error in Cell::LAM when calculating the cathode potential for LAM: " << e << ". Throwing it on" << endl
					 << flush;
				throw e;
			}

			etap_LAM = OCVpt + etap - OCVnmc;											 // equation (9) from the paper
			kt = lamparam.lam3k * exp(lamparam.lam3k_T / Rg * (1 / T_ref - 1 / getT())); // Arrhenius law
			idiss = -kt * exp(n * F / Rg / getT() * etap_LAM) / (n * F);				 // equation (8) from the paper
			idiss = max(idiss, -5 * pow(10, -6));										 // cap the effect at 5e-6 to avoid a very fast drop of capacity (which could  cause an error)
																						 // a value of 5e-6 gives dap = -3.5. The initial value is about 17000, so the cell is dead in 5,000 seconds
																						 // so this cap is quite high

			// Add the effects of this model
			depp += idiss;
			// assume the other effects are 0
			denn += 0;
			dthickpp += 0;
			dthicknn += 0;
			dapp += 0;
			dann += 0;
			break;
		case 4: // Model by Narayanrao, Joglekar, Inguva, Journal of the Electrochemical Society 160 (1), 2012
			// Add the effects of this model
			dapp += -lamparam.lam4p * getAp(); // equation (7) from the paper
			dann += -lamparam.lam4n * getAn();
			// assume the other effects are 0
			dthickpp += 0;
			dthicknn += 0;
			depp += 0;
			denn += 0;
			// todo: this will lead to errors in the code. ValidStates says a = 3e/R
			//  so if da/dt != 0 but de/dt = 0, this constraint will no longer be valid
			break;
		default: // unknown degradation model
			std::cerr << "ERROR in Cell::LAM, unknown LAM degradation model with identifier " << deg_id.LAM_id[i] << ". Only values 0 to 4 are allowed. Throw an error" << '\n';
			throw 106;
			break;
		}
	}

	// Make the output variables
	*dthickp = dthickpp;
	*dthickn = dthicknn;
	*dap = dapp;
	*dan = dann;
	*dep = depp;
	*den = denn;
}

void Cell_SPM::LiPlating(double OCVnt, double etan, double *ipl)
{
	/*
	 * Function to simulate the effect of lithium plating
	 *
	 * IN
	 * OCVnt 	the OCV of the negative electrode at the current battery temperature [V]
	 * etan 	the overpotential at the negative electrode [V]
	 *
	 * OUT
	 * ipl 		current density for the plating side-reaction [A m-2]
	 *
	 * THROWS
	 * 106 		illegal value in id
	 */

	// Arrhenius relation for temperature-dependent plating parameters
	double kplt = plparam.pl1k * exp(plparam.pl1k_T / Rg * (1 / T_ref - 1 / getT())); // Rate constant

	if (deg_id.pl_id == 0) // no plating
		*ipl = 0;
	else if (deg_id.pl_id == 1) // Yang, Leng, Zhang, Ge, Wang, Journal of Power Sources 360, 2017
		*ipl = npl * F * kplt * exp(-n * F / (Rg * getT()) * alphapl * (OCVnt + etan - OCVpl + rsei * getDelta() * getI()));
	else
	{
		std::cerr << "ERROR in Cell::LiPlating, illegal degradation model identifier " << deg_id.pl_id << ", only values 0 and 1 are allowed. Throwing an error" << '\n';
		throw 106;
	}
}

void Cell_SPM::dState_diffusion(bool print, double dzp[], double dzn[], double &dSoC)
{
	/*
	 * Calculate just the diffusion PDE
	 * IN
	 * print 	boolean indicating if we want to print error messages or not
	 * 				if true, error messages are printed
	 * 				if false no error messages are printed (but the error will still be thrown)
	 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
	 * blockDegradation 	if true, battery degradation is ignored (i.e. the time derivatives of those states are 0)
	 *
	 * OUT
	 * dstates	change in the states
	 * 		dzp			time derivative of the transformed concentration at the positive inner nodes of the positive electrode (dzp/dt)
	 * 		dzn			time derivative of the transformed concentration at the positive inner nodes of the negative electrode (dzn/dt)
	 *
	 */

#if TIMING
	std::clock_t tstart = std::clock();
#endif

	// current density
	double i_app = getI() / elec_surf;					  // current density on the electrode surfaces [A m-2]
	double jp = -i_app / (getAp() * n * F * getThickp()); // molar flux on the positive particle [mol m-2 s-1]
	double jn = i_app / (getAn() * n * F * getThickn());  // molar flux on the negative particle [mol m-2 s-1]

	// Arrhenius relation for temperature-dependent parameters
	double Dpt = getDp() * exp(Dp_T / Rg * (1 / T_ref - 1 / getT())); // Diffusion constant at the positive electrode at the cell's temperature [m s-1]
	double Dnt = getDn() * exp(Dn_T / Rg * (1 / T_ref - 1 / getT())); // Diffusion constant at the negative electrode at the cell's temperature [m s-1]

	// Calculate the effect of the main li-reaction on the (transformed) concentration
	double ctep, cten;
	for (int j = 0; j < CELL_NCH; j++)
	{								 // loop for each row of the matrix-vector product A * z
		ctep = M.Ap[j] * getZp()[j]; // A is diagonal, so the array M.A has only the diagonal elements
		cten = M.An[j] * getZn()[j];
		dzp[j] = (Dpt * ctep + M.Bp[j] * jp); // dz/dt = D * A * z + B * j
		dzn[j] = (Dnt * cten + M.Bn[j] * jn);
	}

	// SoC
	dSoC = -getI() / (getCap() * 3600); // dSoC	 	state of charge

#if TIMING
	T_dstate += (std::clock() - tstart) / (double)CLOCKS_PER_SEC; // time in seconds
#endif
}

void Cell_SPM::dState_thermal(bool print, double &dQgen)
{
/*
 * Calculate the internal heat generation
 * IN
 * print 	boolean indicating if we want to print error messages or not
 * 				if true, error messages are printed
 * 				if false no error messages are printed (but the error will still be thrown)
 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
 * blockDegradation 	if true, battery degradation is ignored (i.e. the time derivatives of those states are 0)
 *
 * OUT
 * dQgen 	change in generated heat
 */
#if TIMING
	std::clock_t tstart = std::clock();
#endif

	// Calculcate the lithium fractions at the surface of the particles
	double cps = getCp_surface(getZp());
	double cns = getCn_surface(getZn());

	// Calculate the overpotentials if needed
	if (!etacell_valid)
	{

		// Calculate the rate constants at the cell's temperature using an Arrhenius relation
		double kpt = kp * exp(kp_T / Rg * (1 / T_ref - 1 / getT()));
		double knt = kn * exp(kn_T / Rg * (1 / T_ref - 1 / getT()));

		// Calculate the overpotential using the Bulter-Volmer equation
		// 		if alpha is 0.5, the Bulter-Volmer relation can be inverted to eta = 2RT / (nF) asinh(x)
		double i_app = getI() / elec_surf;															  // current density on the electrodes [I m-2]
		double i0p = kpt * n * F * C_elec_sqrt * sqrt(cps) * sqrt(Cmaxpos - cps);					  // exchange current density of the positive electrode
		double i0n = knt * n * F * C_elec_sqrt * sqrt(cns) * sqrt(Cmaxneg - cns);					  // exchange current density of the negative electrode
		etapcell = (2 * Rg * getT()) / (n * F) * asinh(-0.5 * i_app / (getAp() * getThickp()) / i0p); // cathode overpotential [V], < 0 on discharge
		etancell = (2 * Rg * getT()) / (n * F) * asinh(0.5 * i_app / (getAn() * getThickn()) / i0n);  // anode overpotential [V],  > 0 on discharge
		etacell_valid = true;
	}

	// Calculate the entropic coefficient
	double dOCV;
	bool bound = true; // in linear interpolation, throw an error if you are outside of the allowed range of the data
	try
	{
		dOCV = linInt(print, bound, dOCV_tot_x, dOCV_tot_y, dOCV_tot_n, cps / Cmaxpos); // entropic coefficient of the entire cell OCV [V K-1]
	}
	catch (int e)
	{
		if (print)
			cout << "Error in Cell_SPM::dState when calculating the entropic coefficient " << e << ". Throwing it on" << endl
				 << flush;
		throw e;
	}

	// temperature model
	// Calculate the thermal sources/sinks/transfers per unit of volume of the battery
	// The battery volume is given by the product of the cell thickness and the electrode surface
	double Qrev = -getI() * getT() * dOCV;		  // reversible heat due to entropy changes [W]
	double Qrea = getI() * (etancell - etapcell); // reaction heat due to the kinetics [W]
	double Qohm = pow(getI(), 2) * getRdc();	  // Ohmic heat due to electrode resistance [W]
	// double Qc = -Qch*SAV*(getT() - T_env);					// cooling with the environment [W m-2]

	// total heat generation and cooling in W
	dQgen = (Qrev + Qrea + Qohm);
	// dQcool = (Qc) * (L*elec_surf);

	// dT = 1/(rho*Cp)*(Qrev+Qrea+Qohm+Qc);					// dT		cell temperature

#if TIMING
	T_dstate += (std::clock() - tstart) / (double)CLOCKS_PER_SEC; // time in seconds
#endif
}

void Cell_SPM::dState_degradation(bool print, double dstates[])
{
/*
 * calculate the effects of degradation
 *
 * IN
 * print 	boolean indicating if we want to print error messages or not
 * 				if true, error messages are printed
 * 				if false no error messages are printed (but the error will still be thrown)
 * 			we need this input from higher level functions because at this point we cannot know if this will be a critical error or not
 *
 * OUT
 * dstates	change in the states
 * 		dzp			time derivative of the transformed concentration at the positive inner nodes of the positive electrode (dzp/dt)
 * 		dzn			time derivative of the transformed concentration at the positive inner nodes of the negative electrode (dzn/dt)
 * 		ddelta 		time derivative of the SEI thickness [m s-1] (ddelta/dt)
 * 		dLLI 		time derivative of the lost lithium inventory [C s-1] (dLLI/dt)
 * 		dthickp 	time derivative of the thickness of the positive electrode [m s-1] (dthickp/dt), <0 (dthickp/dt)
 * 		dthickn		time derivative of the thickness of the negative electrode [m s-1] (dthickn/dt), <0 (dthickn/dt)
 * 		dep			time derivative of the volume fraction in the positive electrode [s-1] (dep/dt)
 * 		den			time derivative of the volume fraction in the negative electrode [s-1] (den/dt)
 * 		dap			time derivative of the effective surface area of the positive electrode [m2 m-3 s-1] (dap/dt)
 * 		dan			time derivative of the effective surface area of the negative electrode [m2 m-3 s-1] (dan/dt)
 * 		dCS 		time derivative of the crack surface [m2 s-1], dCS/st > 0 (dCS/dt)
 * 		dDp 		time derivative of the diffusion constant at reference temperature of the positive electrode [m s-1 s-1] (dDp/dt)
 * 		dDn			time derivative of the diffusion constant at reference temperature of the negative electrode [m s-1 s-1] (dDn/dt)
 * 		drdc_p		time derivative of the electrode resistance [Ohm m2 s-1] (dR/dt)
 * 		drdc_n		time derivative of the electrode resistance [Ohm m2 s-1] (dR/dt)
 * 		drdc_cc		time derivative of the electrode resistance [Ohm m2 s-1] (dR/dt)
 * 		ddelta_pl 	time derivative of the thickness of the plated lithium layer [m s-1] (ddelta_pl/dt)
 * 		dSoC
 * 		dT			time derivative of the battery temperature [K s-1] (dT/dt)
 * 		dI
 */
#if TIMING
	std::clock_t tstart = std::clock();
#endif

	// If we block degradation, just fill with all 0
	if (bockDegAndTherm)
	{
		for (int i = 0; i < 2 * CELL_NCH + 18; i++)
			dstates[i] = 0;
		return;
	}

	// Calculcate the lithium fractions at the surface of the particles
	double cps = getCp_surface(getZp());
	double cns = getCn_surface(getZn());

	// Calculate the overpotentials if needed
	if (!etacell_valid)
	{

		// Calculate the rate constants at the cell's temperature using an Arrhenius relation
		double kpt = kp * exp(kp_T / Rg * (1 / T_ref - 1 / getT()));
		double knt = kn * exp(kn_T / Rg * (1 / T_ref - 1 / getT()));

		// Calculate the overpotential using the Bulter-Volmer equation
		// 		if alpha is 0.5, the Bulter-Volmer relation can be inverted to eta = 2RT / (nF) asinh(x)
		double i_app = getI() / elec_surf;															  // current density on the electrodes [I m-2]
		double i0p = kpt * n * F * C_elec_sqrt * sqrt(cps) * sqrt(Cmaxpos - cps);					  // exchange current density of the positive electrode
		double i0n = knt * n * F * C_elec_sqrt * sqrt(cns) * sqrt(Cmaxneg - cns);					  // exchange current density of the negative electrode
		etapcell = (2 * Rg * getT()) / (n * F) * asinh(-0.5 * i_app / (getAp() * getThickp()) / i0p); // cathode overpotential [V], < 0 on discharge
		etancell = (2 * Rg * getT()) / (n * F) * asinh(0.5 * i_app / (getAn() * getThickn()) / i0n);  // anode overpotential [V],  > 0 on discharge
		etacell_valid = true;
	}

	// calculate the anode potential (needed for various degradation models)
	double OCV_n, dOCVn; // anode potential at reference temperature and entropic coefficient
	bool bound = true;	 // in linear interpolation, throw an error if you are outside of the allowed range of the data
	try
	{
		dOCVn = linInt(print, bound, dOCV_neg_x, dOCV_neg_y, dOCV_neg_n, cns / Cmaxneg); // entropic coefficient of the anode potential [V K-1]
		OCV_n = linInt(print, bound, OCV_neg_x, OCV_neg_y, OCV_neg_n, cps / Cmaxpos);	 // anode potential [V]
	}
	catch (int e)
	{
		if (print)
			cout << "Error in Cell_SPM::dState when calculating the anode potential " << e << ". Throwing it on" << endl
				 << flush;
		throw e;
	}
	double OCVnt = OCV_n + (getT() - T_ref) * dOCVn; // anode potential at the cell's temperature [V]

	// SEI growth
	double isei;			 // current density of the SEI growth side reaction [A m-2]
	double den_sei;			 // decrease in volume fraction due to SEI growth [s-1]
	double dznsei[CELL_NCH]; // additional diffusion in the anode due to isei
	try
	{
		SEI(OCVnt, etancell, &isei, &den_sei);
	}
	catch (int e)
	{
		if (print)
			cout << "Error in Cell::dState when calculating the effect of SEI growth: " << e << ". Throwing it on" << endl
				 << flush;
		throw e;
	}

	// Subtract Li from negative electrode (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + isei/nF)
	for (int j = 0; j < CELL_NCH; j++)
		dznsei[j] = (M.Bn[j] * isei / (nsei * F));

	// crack growth leading to additional exposed surface area
	double isei_multiplyer;		// relative increase in isei due to additional SEI growth on the extra exposed surface area [-]
	double dCS;					// increase in crack surface area [m2 s-1]
	double dDn;					// change in negative diffusion constant [m s-1 s-1]
	double dznsei_CS[CELL_NCH]; // additional diffusion in the anode due to extra SEI growth on the crack surface
	try
	{
		CS(OCVnt, etancell, &isei_multiplyer, &dCS, &dDn);
	}
	catch (int e)
	{
		if (print)
			cout << "Error in Cell::dState when calculating the effect of crack growth: " << e << ". Throwing it on" << endl
				 << flush;
		throw e;
	}

	// crack surface leads to extra SEI growth because the exposed surface area increases.
	// (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + isei/nF + isei_CS/nF)
	double isei_CS = isei * isei_multiplyer; // extra SEI side reaction current density due to the crack surface area [A m-2]
	for (int j = 0; j < CELL_NCH; j++)
		dznsei_CS[j] = (M.Bn[j] * isei_CS / (nsei * F));

	// loss of active material LAM
	double dthickp, dthickn, dap, dan, dep, den; // change in geometric parameters describing the amount of active material
	try
	{
		LAM(print, cps / Cmaxpos, etapcell, &dthickp, &dthickn, &dap, &dan, &dep, &den);
	}
	catch (int e)
	{
		if (print)
			cout << "Error in Cell::dState when calculating the LAM: " << e << ". Throwing it on" << endl
				 << flush;
		throw e;
	}

	// lithium plating
	double ipl;				 // current density of the plating side reaction [A m-2]
	double dzn_pl[CELL_NCH]; // additional diffusion in the anode due to ipl
	try
	{
		LiPlating(OCVnt, etancell, &ipl);
	}
	catch (int e)
	{
		if (print)
			cout << "Error in Cell::dState when calculating the lithium plating: " << e << ". Throwing it on" << endl
				 << flush;
		throw e;
	}

	// Subtract Li from negative electrode (like an extra current density -> add it in the boundary conditions: dCn/dx =  jn + ipl/nF)
	for (int j = 0; j < CELL_NCH; j++)
		dzn_pl[j] = (M.Bn[j] * ipl / (npl * F));

	// output
	for (int j = 0; j < CELL_NCH; j++)
	{
		dstates[j] = 0; // dzp		diffusion
		dstates[CELL_NCH + j] = (dznsei[j] + dznsei_CS[j] + dzn_pl[j]);
	}
	dstates[2 * CELL_NCH + 0] = isei / (nsei * F * rhosei);									// ddelta	SEI thickness
	dstates[2 * CELL_NCH + 1] = (isei + isei_CS + ipl) * elec_surf * getThickn() * getAn(); // dLLI		lost lithium
	dstates[2 * CELL_NCH + 2] = dthickp;													// dthickp 	electrode thickness
	dstates[2 * CELL_NCH + 3] = dthickn;													// dthickn
	dstates[2 * CELL_NCH + 4] = dep;														// dep		volume fraction of active material
	dstates[2 * CELL_NCH + 5] = den + den_sei;												// den
	dstates[2 * CELL_NCH + 6] = dap + 3 / Rp * dstates[2 * CELL_NCH + 4];					// dap		effective surface are, a = 3 e/R
	dstates[2 * CELL_NCH + 7] = dan + 3 / Rn * dstates[2 * CELL_NCH + 5];					// dan
	dstates[2 * CELL_NCH + 8] = dCS;														// dCS		surface area of the cracks
	dstates[2 * CELL_NCH + 9] = 0;															// dDp 		diffusion constant
	dstates[2 * CELL_NCH + 10] = dDn;														// dDn
	dstates[2 * CELL_NCH + 11] = 0;															// drdc_p 	cathode resistance
	dstates[2 * CELL_NCH + 12] = 0;															// drdc_n 	anode resistance
	dstates[2 * CELL_NCH + 13] = 0;															// drdc_cc 	current collector resistance
	dstates[2 * CELL_NCH + 14] = ipl / (npl * F * rhopl);									// ddelta_pl thickness of the plated lithium

#if TIMING
	T_dstate += (std::clock() - tstart) / (double)CLOCKS_PER_SEC; // time in seconds
#endif
}

double Cell_SPM::getThermalSurface()
{
	/*
	 * Assume heat transfer is on one side of the pouch (one of the faces).
	 * We assume that on the other side is another cell, with which this one will also exchange heat
	 *
	 * So even though the total surface area = 2*Acell (+small edges on top, bottom and sides), Atherm is just Acell
	 */
	return Acell;
}

double Cell_SPM::thermalModel_cell()
{
	/*
	 * Calculate the thermal model for this cell on its own.
	 * We do not account for neighbouring cells, and this cell is cooled convectively by the environment
	 */

	// total heat generation since last time the temperature was updated
	double Etot = Therm_Qgen;

	// cooling with the environment
	double Qc = Qch * getThermalSurface() * (T_ENV - getT()) * Therm_time; // cooling with the environment [W m-3]
	Etot += Qc;

	// Calculate the new temperature
	// rho * cp * dT/dt = Qtot / V
	// 		where 	Qtot is total power in W
	// 				V is the cell's volume L * elec_surf
	// so integrated over time this is
	// rho * cp * (Tnew - Told) = Etot / V
	double Tnew = getT() + Etot / (rho * Cp * L * elec_surf);

	// Check the new temperature is valid, and if so, set it
	if (Tnew > Tmax || Tnew < Tmin || isnan(Tnew))
	{
		if (verbose_gl >= v_crit_gl)
		{
			std::cerr << "ERROR in Cell_SPM::thermalModel, the new temperature of " << Tnew << " is outside the allowed range from " << Tmin << " to " << Tmax;
			std::cerr << ". The time since the last time this function was called is " << Therm_time << '\n';

			cout << "Total thermal energy " << Etot << ". internal heat generation " << Therm_Qgen << endl;
			cout << "giving change in temperature: " << Etot / (rho * Cp * L * elec_surf) << endl
				 << flush;
		}
		throw 9;
	}

	// setting the temperature is done by the parent module. else some cells will update their T before others, and we get inconsistencies
	// 		e.g. exchange from cell 2 to this cell will not be same as from this cell to cell 2 since T of this cell would have changed.

	// Reset the cumulative thermal variables
	Therm_Qgen = 0;
	Therm_time = 0;

	return Tnew;
}

double Cell_SPM::thermalModel_coupled(int Nneighbours, double Tneighbours[], double Kneighbours[], double Aneighb[], double tim)
{
	/*
	 * Calculate the thermal model for this cell and update its cell temperature.
	 * The heat exchange in [W] with element i is given by
	 * 		Qc = Kneighbours[i]*A*(Tneighbours[i] - getT());
	 * We assume all temperatures were constant for the past period, such that the thermal energy in [J] is given by
	 * 		Ec = Qc * Therm_time
	 * This is added up with the internal heath generated, and the cell's temperature is returned
	 *
	 * IN
	 * Nneigtbours 		the number of neighbouring cells / cooling systems etc.
	 * Tneighbours  	array with the temperature of the neighbouring  elements [K]
	 * Kneighbours		array with the heat transfer constants of the neighbouring elements, k or h
	 * Aneighb 			array with the surface area of the neigbouring cell
	 * 						the heat transfer happens over the smaller one of Aneighb[i] and the area of this cell
	 * tim 				the time since the last thermal model calculation [for verification]
	 *
	 * throws
	 * 8 	invalid time keeping
	 * 9 	invalid module temperature
	 */

	// check the time since the last checkup is not too large.
	// if the parent has not called the thermal model for a while, the equation becomes unstable
	// 		cause E = time * kA dT, so even a small dT will cause a huge E, and therefore a very large temperature swint
	//
	if (Therm_time > 15 * 60)
	{
		if (verbose_gl >= v_noncrit_gl)
			cout << "Warning in Cell_SPM::thermalModel, the time since this function was called last is very large, " << Therm_time << " which might lead to excessive temperature variations" << endl;
	}

	// then check whether our internal time keeping matches up with the external one
	if (abs(Therm_time - tim) > 1)
	{
		if (verbose_gl >= v_crit_gl)
		{
			std::cerr << "ERROR in Cell_SPM::thermalModel, according to the cell's internal timing, " << Therm_time << "s have passed since the last thermal model solution.";
			std::cerr << " The external time provided was " << tim << "s, which is more than 1s difference. Throwing an error" << endl;
		}
		throw 8;
	}

	// calculate the total thermal balance
	double Etot = Therm_Qgen;
	double Atherm;
	for (int i = 0; i < Nneighbours; i++)
	{
		Atherm = min(Aneighb[i], getThermalSurface());
		Etot += Kneighbours[i] * Atherm * (Tneighbours[i] - getT()) * Therm_time;
	}

	// Calculate the new temperature
	// rho * cp * dT/dt = Qtot / V
	// 		where 	Qtot is total power in W
	// 				V is the cell's volume L * elec_surf
	// so integrated over time this is
	// rho * cp * (Tnew - Told) = Etot / V
	double Tnew = getT() + Etot / (rho * Cp * L * elec_surf);

	// Check the new temperature is valid, and if so, set it
	if (Tnew > Tmax || Tnew < Tmin || isnan(Tnew))
	{
		if (verbose_gl >= v_crit_gl)
		{
			std::cerr << "ERROR in Cell_SPM::thermalModel, the new temperature of " << Tnew << " is outside the allowed range from " << Tmin << " to " << Tmax;
			std::cerr << ". The time since the last time this function was called is " << Therm_time << '\n';

			cout << "Total thermal energy " << Etot << ". internal heat generation " << Therm_Qgen << " and external contributions as below" << endl;
			for (int i = 0; i < Nneighbours; i++)
				cout << Aneighb[i] << ", " << Kneighbours[i] << "," << Tneighbours[i] << "," << Kneighbours[i] * Atherm * (Tneighbours[i] - getT()) * Therm_time << endl;
			cout << "giving change in temperature: " << Etot / (rho * Cp * L * elec_surf) << endl
				 << flush;
		}
		throw 9;
	}

	// setting the temperature is done by the parent module. else some cells will update their T before others, and we get inconsistencies
	// 		e.g. exchange from cell 2 to this cell will not be same as from this cell to cell 2 since T of this cell would have changed.

	// Reset the cumulative thermal variables
	Therm_Qgentot += Therm_Qgen;
	Therm_Qgen = 0;
	Therm_time = 0;

	return Tnew;
}

double Cell_SPM::thermal_getTotalHeat()
{
	return Therm_Qgentot; // total heat generation [J] since start of cell's life
}

double Cell_SPM::thermalModel(int Nneighbours, double Tneighbours[], double Kneighbours[], double Aneighb[], double tim)
{
	/*
	 * Calculate the thermal model for this cell and update its cell temperature.
	 * The heat exchange in [W] with element i is given by
	 * 		Qc = Kneighbours[i]*A*(Tneighbours[i] - getT());
	 * We assume all temperatures were constant for the past period, such that the thermal energy in [J] is given by
	 * 		Ec = Qc * Therm_time
	 * This is added up with the internal heath generated, and the cell's temperature is returned
	 *
	 * IN
	 * Nneigtbours 		the number of neighbouring cells / cooling systems etc.
	 * Tneighbours  	array with the temperature of the neighbouring  elements [K]
	 * Kneighbours		array with the heat transfer constants of the neighbouring elements, k or h
	 * Aneighb 			array with the surface area of the neigbouring cell
	 * 						the heat transfer happens over the smaller one of Aneighb[i] and the area of this cell
	 * tim 				the time since the last thermal model calculation [for verification]
	 *
	 * throws
	 * 8 	invalid time keeping
	 * 9 	invalid module temperature
	 */

	double Tnew;

	try
	{
		if (T_MODEL == 0)
			Tnew = getT();
		else if (T_MODEL == 1)
			Tnew = thermalModel_cell();
		if (T_MODEL == 2)
			Tnew = thermalModel_coupled(Nneighbours, Tneighbours, Kneighbours, Aneighb, tim);
	}
	catch (int e)
	{

		// indicate we have ran the thermal model
		Therm_Qgen = 0;
		Therm_time = 0;
		throw e;
	}

	// setting the temperature is done by the parent module. else some cells will update their T before others, and we get inconsistencies
	// 		e.g. exchange from cell 2 to this cell will not be same as from this cell to cell 2 since T of this cell would have changed.

	// Reset the cumulative thermal variables
	Therm_Qgen = 0;
	Therm_time = 0;

	return Tnew;
}

void Cell_SPM::timeStep_CC(double dt, bool addData, int nstep)
{
	/*
	 * take a number of time steps with a constant current
	 *
	 * The diffusion model is resolved to every time step of dt
	 * The thermal and degradation models are resolved once, for dt*nstep seconds
	 *
	 */

	// check the time step is positive
	if (dt < 0)
	{
		if (verbose_gl >= v_crit_gl)
			std::cerr << "ERROR in Cell_ECM::timeStep_CC, the time step dt must be 0 or positive, but has value " << dt << '\n';
		throw 10;
	}

	bool print = true;

	// Update the stress values stored in the attributes with the stress of the previous time step
	s_dai_p_prev = s_dai_p;		// Dai's stress in the positive particle in the previous time step
	s_dai_n_prev = s_dai_n;		// Dai's stress in the negative particle in the previous time step
	s_lares_n_prev = s_lares_n; // Laresgoiti's stress in the negative particle in the previous time step
	s_dt = nstep * dt;

	// *********************************************  Resolve the diffusion model for every dt time step ****************************************************************************************
	for (int t = 0; t < nstep; t++)
	{

		// Calculate the time derivatives
		double dzp[CELL_NCH], dzn[CELL_NCH], dSoC;
		try
		{
			dState_diffusion(print, dzp, dzn, dSoC);
		}
		catch (int e)
		{
			if (verbose_gl >= v_crit_gl)
				cout << "error in SPM cell " << getFullID() << " when calculating the time derivatives of the diffusion PDE: " << e << ", throwing it on" << endl
					 << flush;
			throw e;
		}

		// forward Euler time integration: s(t+1) = s(t) + ds/dt * dt
		try
		{
			for (int i = 0; i < CELL_NCH; i++)
			{
				dzp[i] *= dt;
				dzn[i] *= dt;
			}
			setdZp(dzp);
			setdZn(dzn);
			setSoC(getSoC() + dt * dSoC, false, true);
			Vcell_valid = false;
			etacell_valid = false;
			s_dai_update = false;
			s_lares_update = false;
		}
		catch (int e)
		{
			if (verbose_gl >= v_crit_gl)
				cout << "error in SPM cell " << getFullID() << " when doing Euler time integration of the diffusion PDE: " << e << ", throwing it on" << endl
					 << flush;
			throw e;
		}

		// increase the cumulative variables of this cell
		if (addData)
		{
#if DATASTORE_CELL > 0
			timetot += dt;
			ahtot += abs(getI()) * dt / 3600.0;
			whtot += abs(getI()) * getV() * dt / 3600.0;
#endif
		}
	}

	// **************************************************** Calculate the thermal model once for the nstep * dt time period *****************************************************************
	if (!bockDegAndTherm)
	{

		// Calculate the internal heat generation of the cell
		double dQgen;
		try
		{
			dState_thermal(print, dQgen);
		}
		catch (int e)
		{
			if (verbose_gl >= v_crit_gl)
				cout << "error in SPM cell " << getFullID() << ", error " << e << " when calculating the time derivatives of the thermal model, throwing it on" << endl;
			throw e;
		}
		Therm_Qgen += dQgen * dt * nstep;
		Therm_time += dt * nstep;

		// If this cell has a parent module, this parent will call the thermal model with the correct parameters
		// which will include heat exchange with the cell's neighbours and cooling from the cooling system of the module.

		// if there is no parent, this cell is a stand-alone cell.
		// We assume convective cooling with an environment
		// If there is no parent, assume we update T every nstep*dt. So update the temperature now
		if (parent == NULL)
		{
			double Tneigh[1] = {T_ENV};
			double Kneigh[1] = {Qch};				  // so cooling Qc = Qch * SAV * dT / (rho*cp) = Qch * A * dT / (rho*cp)
			double Atherm[1] = {getThermalSurface()}; // calculate the surface of this cell
			setT(thermalModel(1, Tneigh, Kneigh, Atherm, Therm_time));
		}
	}

	// else it is the responsibility of the parent to call the thermal model function with the correct settings

	// ******************************************************* Calculate the degradation model once for the nstep * dt time period **************************************************************

	if (!bockDegAndTherm)
	{
		// Calculate the stress values stored in the attributes for the stress in this time step
		if (s_dai) // only a few degradation models actually need the stress according to Dai, so only calculate it if needed
			updateDaiStress();
		if (s_lares) // only a few degradation models need the stress according to Laresgoiti
			updateLaresgoitiStress(true);
		// Calculate the time derivatives
		double ds[CELL_NSTATE_MAX];
		try
		{
			dState_degradation(print, ds);
		}
		catch (int e)
		{
			if (verbose_gl >= v_crit_gl)
				cout << "error in SPM cell " << getFullID() << " when calculating the time derivatives of the degradation: " << e << ", throwing it on" << endl
					 << flush;
			throw e;
		}

		// forward Euler time integration: s(t+1) = s(t) + ds/dt * dt
		// degradation accumulation
		for (int i = 0; i < 2 * CELL_NCH + 15; i++)
			degState[i] += ds[i] * dt * nstep;

		// update State [for now every time step, in future only when degradation model is called]
		for (int i = 0; i < 2 * CELL_NCH + 15; i++) // loop for all degradation states (no SoC, I and T)
			state[i] += degState[i];
		for (int i = 0; i < 2 * CELL_NCH + 15; i++)
			degState[i] = 0;
	}
}

shared_ptr<StorageUnit> Cell_SPM::copy()
{
	/*
	 * Return a copy.
	 * guaranteed properties which are copied
	 * 		states
	 * 		degid
	 *
	 * not copied across
	 * 		parent
	 * 		cycling data
	 */

	// get the states of this cell
	int nin = CELL_NSTATE_MAX;
	double s[nin];
	int nout;
	getStates(s, nin, nout);

	// make a new cell and copy the states across
	shared_ptr<Cell_SPM> c(new Cell_SPM(getID(), deg_id, var_cap, var_R, var_degSEI, var_degLAM));
	c->setStates(s, nout, false, true); // this sets the resistance
	c->setThroughput(timetot, ahtot, whtot);

// copy accross the statistics if we are storing stats
#if DATASTORE_CELL == 1
	// double HI[DATASTORE_NHIST], HV[DATASTORE_NHIST], HT[DATASTORE_NHIST];

	// Get the pointers to the start of both arrays
	double *histI, *histV, *histT;
	double *HI, *HV, *HT;
	int nout1, nout2;
	getHist(&histI, &histV, &histT, nout1);
	c->getHist(&HI, &HV, &HT, nout2);
	assert(nout1 == nout2);

	// Copy across from this one to the copy
	// (Note that this is dirty code, we should not be accessing private object directly through their memory addresses
	// 	the alternative is to make a full array, copy the values from here to the array, pass the array to a new function setStats
	//  and then the setStats of the copy-cell can transfer these values to its own histogram
	for (int i = 0; i < nout1; i++)
	{
		HI[i] = histI[i];
		HV[i] = histV[i];
		HT[i] = histT[i];
	}
#endif

#if TIMING
	c->setTimings(T_getOCV, T_getV, T_dstate, T_validStates, T_setStates);
#endif

	return c;
}

void Cell_SPM::storeData()
{
	/*
	 * Add another data point in the array.
	 *
	 * THROWS
	 * 5 	the data arrays are full
	 */

#if DATASTORE_CELL == 0
	if (verbose_gl > v_noncrit_gl)
		std::cerr << "ERROR in Cell_SPM::storeData, the settings in global.h are forbidding from storing data. CELL_NDATA_MAX = " << CELL_NDATA_MAX << '\n';
#endif

#if DATASTORE_CELL == 1
	hist_I_ind = increaseBin(hist_I_x, hist_I, getI(), hist_I_ind);
	hist_V_ind = increaseBin(hist_V_x, hist_V, getV(), hist_V_ind);
	hist_T_ind = increaseBin(hist_T_x, hist_T, getT(), hist_T_ind);
	// cout<<"SPM cell storeData in bins "<<hist_I_ind<<", "<<hist_V_ind<<", "<<hist_T_ind<<", number of values in these bins: "<<hist_I[hist_I_ind]<<", "<<hist_V[hist_V_ind]<<", "<<hist_T[hist_T_ind]<<", "<<endl;
#endif

// Then add new data specifically for SPM cells
#if DATASTORE_CELL == 2

	// if the arrays are full, throw an error since the data needs to be written
	if (index >= CELL_NDATA_MAX)
	{
		if (verbose_gl > v_crit_gl)
			std::cerr << "ERROR in Cell_SPM::storeData, the arrays are full. Write the data to a file and reset the data buffers for this cell" << '\n';
		throw 5;
	}

	// add a new point in all the arrays
	LLI[index] = getLLI() / 3600.0; // convert As to Ah
	LAM_en[index] = getEn();
	LAM_thickn[index] = getThickn();
	Rtot[index] = getRdc();

	// index++ is done in Cell::storeData
#endif

	// Finally call the function from Cell to add those fields
	// it also updates index
	Cell::storeData();
}

void Cell_SPM::writeData(string prefix)
{
/*
 * Writes data to a csv file.
 * The name of the csv file starts with the value of prefix, after which the identification string of this cell is appended
 *
 * Depending on the value of DATASTORE_CELL, different things are written
 * 	0 	nothing
 * 	1 	general info about the cell and usage statistics in file xxx_cellStats.csv
 * 	2 	cycling data (I, V, T, LLI, en, thickn, Rdc at every time step) in file xxx_cellData.csv
 *
 */

// store nothing
#if DATASTORE_CELL == 0
	if (verbose_gl > v_noncrit_gl)
		std::cerr << "ERROR in Cell_SPM::writeData, the settings in global.h are forbidding from storing data. CELL_NDATA_MAX = " << CELL_NDATA_MAX << '\n';
#endif

// store histograms and degradation state of cell utilisation
#if DATASTORE_CELL == 1
	string name = prefix + "_" + getFullID() + "_cellStats.csv"; // name of the file, start with the full hierarchy-ID to identify this cell
																 /*
																  * first line gives the total utilistion of the cell (time, charge and energy throughput)
																  * second line gives the parameters cell-to-cell variation (i.e. what were the values for this cell)
																  * third line gives the full battery state (which for an SPM cell indicates which degradation mechanism was active)
																  * 4-5 is empty
																  * 6-106 give the histogram values (number of data points in each bin), for I, V and T
																  * 107-108 is empty
																  * 109-208 gives the edges of the bins of the histogram: edge(i-1) < bin(i) < edge(i)
																  */

	// Get the state of the battery
	double st[CELL_NSTATE_MAX], varstate[5];
	int nout, noutdeg;
	getStates(st, CELL_NSTATE_MAX, nout);
	getVariations(varstate, 5, noutdeg);

	// append the new data to the existing file
	ofstream file;
	file.open(name, std::ios_base::app);

	if (!file.is_open())
	{
		std::cerr << "ERROR in Cell::writeData, could not open file " << name << endl;
		throw 11;
	}

	// write throughput data, cell-to-cell variations and the battery state:
	file << timetot << "," << ahtot << "," << whtot << "\n";
	for (int i = 0; i < noutdeg; i++)
		file << varstate[i] << ",";
	file << "\n";
	for (int i = 0; i < nout; i++)
		file << st[i] << ",";
	file << "\n";
	file << "\n";
	file << "\n";

	// write bin values
	for (int i = 0; i < DATASTORE_NHIST; i++)
		file << hist_I[i] << "," << hist_V[i] << "," << hist_T[i] << "\n";
	file << "\n";
	file << "\n";

	// write bin edges
	for (int i = 0; i < DATASTORE_NHIST - 1; i++)
		file << hist_I_x[i] << "," << hist_V_x[i] << "," << hist_T_x[i] << "\n";
	// data is in bin i if edge[i-1] <= data < edge[i]
	// except bin[0] if data < edge[0]
	//    and bin[N] if data >= edge[N-1]
	file << "\n";
	file << "\n";
	file << "\n";

	file.close();
#endif

// store cycling data at every time step
#if DATASTORE_CELL == 2
	string name = prefix + "_" + getFullID() + "_cellData.csv"; // name of the file, start with the full hierarchy-ID to identify this cell

	// append the new data to the existing file
	ofstream file;
	file.open(name, std::ios_base::app);

	if (!file.is_open())
	{
		std::cerr << "ERROR in Cell::writeData, could not open file " << name << endl;
		throw 11;
	}

	for (int i = 0; i < index; i++)
	{
		file << Itot[i] << "," << Vtot[i] << "," << SoCtot[i] << "," << Ttot[i] << "," << Ahtot[i] << "," << Whtot[i] << "," << Timetot[i]; // states from Cell, keep in same order to read-scripts stay valid
		file << "," << LLI[i] << "," << LAM_en[i] << "," << LAM_thickn[i] << "," << Rtot[i] << "\n";										// additional states for SPM
	}

	file.close();

	// reset the index to 0 since we can overwrite the stored data
	index = 0;
#endif
}

void Cell_SPM::getTimings(double &gOCV, double &tV, double &tdstate, double &tvalidstate, double &tsetState)
{
	/*
	 * return the timings for the batteries
	 */

	gOCV = 0;
	tV = 0;
	tdstate = 0;
	tvalidstate = 0;
	tsetState = 0;

#if TIMING
	gOCV = T_getOCV;
	tV = T_getV;
	tdstate = T_dstate;
	tvalidstate = T_validStates;
	tsetState = T_setStates;
#endif
}
void Cell_SPM::setTimings(double gOCV, double tV, double tdstate, double tvalidstate, double tsetState)
{
	/*
	 * return the timings for the batteries
	 */

#if TIMING
	T_getOCV = gOCV;
	T_getV = tV;
	T_dstate = tdstate;
	T_validStates = tvalidstate;
	T_setStates = tsetState;
#endif
}
