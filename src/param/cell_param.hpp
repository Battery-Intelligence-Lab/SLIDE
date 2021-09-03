

#pragma once

#include <vector>
#include <string>

#include "../read_CSVfiles.h"

// Define a structure with the identifications of which degradation model(s) to use
struct DEG_ID
{
	// identifiers for degradation models
	// Each array is made with a length of 'len', which is the maximum number of models per mechanisms.
	// If the user wants to use more models, you only have to increase the value of 'len' to the number you want to use
	// and change the '10' in the definition of the arrays
	int len = 10;	  // length of the arrays with identifications of which models to use
	int SEI_id[10];	  // array with identifications to decide which SEI models to use. Max length 10
					  /* 				0	no SEI growth
							 * 				1 	kinetic model only (Tafel kinetics)
							 * 				2 	Pinson&Bazant model: linear diffusion + Tafel kinetics
							 * 				3	Christensen and Newman model
							 */
	int SEI_n;		  // SEI_N 	number of SEI models to use (length of SEI_ID)
	int SEI_porosity; // integer deciding whether we reduce the active volume fraction due to SEI growth
					  /* 				0	don't reduce it
							 * 				1	use correlation from Ashwin et al. 2016
							 */
	int CS_id[10];	  // array with identifications for which model to use for surface cracking. Max length 10
					  /* 				0 	no surface cracking
							 * 				1 	Laresgoiti's stress + crack growth model
							 * 				2 	Dai stress model + Laresgoiti crack growth
							 * 				3 	model based on Deshpande and Bernardi, 2017
							 * 				4 	model from Barai et al
							 * 				5 	model from Ekstrom et al
							 */
	int CS_n;		  // number of surface crack growth models to use (length of CS_ID)
	int CS_diffusion; // integer deciding whether we reduce the negative diffusion constant due to surface cracks
					  /* 				0 	don't decrease diffusion
							 * 				1	decrease according to Barai et al. 2015
							 */
	int LAM_id[10];	  // array with the integers deciding which models is to be used for loss of active material. Max length 10
					  /* 				0 	no LAM
							 * 				1	Dai's stress model and Laresgoiti's correlation to get LAM
							 * 				2	delacourt's	correlation between abs(j) and porosity
							 * 				3 	Kindermann's model for cathode dissolution: tafel kinetics for increased porosity
							 * 				4 	Narayanrao's correlation which decreases the effective surface area proportionally to itself and j
							 */
	int LAM_n;		  // number of LAM models to be used (length of LAM_id)
	int pl_id;		  // integer deciding which model is to be used for li-plating
					  /* 				0 	no plating
							 * 				1	Yang et al thermodynamic plating (Tafel kinetics)
							 */

	std::string print() const
	{
		/*
	 * Function to get a string representation of the struct with the degradation settings
	 * This string is part of the names of the subfolders in which results are written.
	 *
	 * IN
	 * degid 	struct with the degradation identifiers
	 *
	 * OUT
	 * string 	string of the degradation identifiers
	 * 			identifiers of the same mechanism are separated by -
	 * 			identifiers of different mechanisms are separated by _
	 * 			e.g. if we use SEI model 1, no SEI porosity effect, no surface cracks, LAM model 2 and LAM model 3 and lithium plating model 1:
	 * 				2-0_0-0_2-3_1
	 * 				2 		SEI model 1
	 * 				0		no SEI porosity
	 * 				0		no surface cracks
	 * 				0 		don't decrease the diffusion due to surface cracks
	 * 				2		LAM model 2
	 * 				3		LAM model 3
	 * 				1		lithium plating model 1
	 *
	 */

		std::string id = "";

		// print SEI models and SEI_porosity (decreasing the active volume fraction due to SEI growth), separated by -
		for (int i = 0; i < SEI_n; i++)
			id += std::to_string(SEI_id[i]) + "-";

		id += std::to_string(SEI_porosity);

		// mechanism separator
		id += "_";

		// print CS models (surface cracking) and CS_diffusion (reducing the diffusion constant because of cracks), separated by -
		for (int i = 0; i < CS_n; i++)
			id += std::to_string(CS_id[i]) + "-";

		id += std::to_string(CS_diffusion);

		// mechanism separator
		id += "_";

		// print LAM models separated by -
		for (int i = 0; i < LAM_n; i++)
		{
			id += std::to_string(LAM_id[i]);
			if (i < LAM_n - 1) // print separator - only between LAM models
				id += "-";	   // not after the last LAM model
		}

		// mechanism separator
		id += "_";

		// print plating model
		id += std::to_string(pl_id);

		// output
		return id;
	}
};

// Define a structure with the fitting parameters of the SEI growth models (SEI)
struct SEIparam
{

	double sei1k;	// rate parameter of the SEI side reaction in the 1st SEI model
	double sei1k_T; // activation energy of sei1k

	double sei2k;	// rate parameter of the SEI side reaction in the 2nd SEI model
	double sei2k_T; // activation energy of sei2k
	double sei2D;	// diffusion constant of the SEI layer in the 2nd SEI model
	double sei2D_T; // activation energy of sei2D

	double sei3k;	// rate parameter of the SEI side reaction in the 3rd SEI model
	double sei3k_T; // activation energy of sei3k
	double sei3D;	// diffusion constant of the SEI layer in the 3rd SEI model
	double sei3D_T; // activation energy of sei3D

	double sei_porosity; // proportionality constant between the SEI growth and the decrease in volume fraction of active material
						 // (because the SEI layer blocks the pores at the surface)
};

// Define a structure with the fitting parameters of the surface crack growth models (CS)
struct CSparam
{

	double CS1alpha; // fitting parameter of the 1st surface crack growth model

	double CS2alpha; // fitting parameter of the 2nd surface crack growth model

	double CS3alpha; // fitting parameter of the 3rd surface crack growth model

	double CS4Amax;	 // maximum crack growth surface for the 4th surface crack growth model
	double CS4alpha; // fitting parameter of the 4th surface crack growth model

	double CS5k;   // rate parameter of the 5th surface crack growth model at reference temperature
	double CS5k_T; // activation energy of CS5k

	double CS_diffusion; // fitting parameter to decrease the diffusion constant due to surface cracks
};

// Define a structure with the fitting parameters of the models for loss of active material (LAM)
struct LAMparam
{
	double lam1p; // fitting parameter for the positive electrode for the 1st LAM model
	double lam1n; // fitting parameter for the negative electrode for the 1st LAM model

	double lam2ap; // fitting parameter 1 at reference temperature for the positive electrode for the 2nd LAM model
	double lam2bp; // fitting parameter 2 at reference temperature for the positive electrode for the 2nd LAM model
	double lam2an; // fitting parameter 1 at reference temperature for the negative electrode for the 2nd LAM model
	double lam2bn; // fitting parameter 2 at reference temperature for the negative electrode for the 2nd LAM model
	double lam2t;  // activation energy for all the parameters of the 2nd LAM model

	double lam3k;	// rate constant at reference temperature for the cathode dissolution side reaction
	double lam3k_T; // activation energy for lam3k

	double lam4p; // fitting parameter for the positive electrode for the 4th LAM model
	double lam4n; // fitting parameter for the negative electrode for the 4th LAM model
};

// Define a structure with the fitting parameters of the li-plating models (PL)
struct PLparam
{
	double pl1k;   // rate constant of the li-plating side reaction at reference temperature in the 1st model
	double pl1k_T; // activation energy of pl1k
};

struct StressParam
{
	// Constants for the stress model
	double omegap; // partial molar volume of positive electrode [m3 mol-1]
	double omegan; // partial molar volume of negative electrode [m3 mol-1]
	double Ep;	   // Young's modulus of positive electrode [GPa]
	double En;	   // Young's modulus of negative electrode [GPa]
	double nup;	   // Poisson's ratio of positive electrode [-]
	double nun;	   // Poisson's ratio of negative electrode [-]

	// values of the stress are often needed. Because it takes very long to calculate them, we calculate them once and store them so we don't need to repeat the same calculation twice
	bool s_dai;			   // do we need to calculate the stress according to Dai's model?
	bool s_lares;		   // do we need to calculate the stress according to Laresgoiti's model?
	bool s_dai_update;	   // boolean to indicate if Dai's stress are up to date with the battery state at this time step
	bool s_lares_update;   // boolean to indicate if Dai's stress are up to date with the battery state at this time step
	double s_dai_p;		   // maximum hydrostatic stress in the positive particle according to Dai's stress model
	double s_dai_n;		   // maximum hydrostatic stress in the negative particle according to Dai's stress model
	double s_lares_n;	   // stress in the negative particle according to Laresgoiti's stress model
	double s_dai_p_prev;   // maximum hydrostatic stress in the previous time step in the positive particle according to Dai's stress model
	double s_dai_n_prev;   // maximum hydrostatic stress in the previous time step in the negative particle according to Dai's stress model
	double s_lares_n_prev; // stress in the previous time step in the negative particle according to Laresgoiti's stress model
};

struct OCVcurves
{
	int OCV_pos_n{0}, OCV_neg_n{0};	  // number of data points in the OCV curve for the cathode/anode
	int dOCV_neg_n{0}, dOCV_tot_n{0}; // number of data points in the entropic coefficient for the anode/entire cell

	bool is_OCV_pos_fixed{false}, is_OCV_neg_fixed{false}, is_dOCV_neg_fixed{false}, is_dOCV_tot_fixed{false};

	std::string namepos{}, nameneg{}, nameentropicC{}, nameentropicCell{};

	std::vector<double> OCV_pos_x, OCV_neg_x;	// lithium fractions of the points of the cathode/anode OCV curve
	std::vector<double> OCV_pos_y, OCV_neg_y;	// voltage vs li/li+ of the points of the cathode/anode OCV curve [V]
	std::vector<double> dOCV_neg_x, dOCV_tot_x; // lithium fractions of the points of the anode/entire cell entropic coefficient curve
	std::vector<double> dOCV_neg_y, dOCV_tot_y; // entropic coefficient curve / the entire cell's entropic coefficient [V K-1]

	OCVcurves(const std::string &_namepos, int _OCV_pos_n,
			  const std::string &_nameneg, int _OCV_neg_n,
			  const std::string &_nameentropicC, int _dOCV_neg_n,
			  const std::string &_nameentropicCell, int _dOCV_tot_n)
		: OCV_pos_n(_OCV_pos_n), OCV_neg_n(_OCV_neg_n), dOCV_neg_n(_dOCV_neg_n), dOCV_tot_n(_dOCV_tot_n),
		  namepos(_namepos), nameneg(_nameneg), nameentropicC(_nameentropicC), nameentropicCell(_nameentropicCell),
		  OCV_pos_x(_OCV_pos_n), OCV_neg_x(_OCV_neg_n), OCV_pos_y(_OCV_pos_n), OCV_neg_y(_OCV_neg_n),
		  dOCV_neg_x(_dOCV_neg_n), dOCV_tot_x(_dOCV_tot_n), dOCV_neg_y(_dOCV_neg_n), dOCV_tot_y(_dOCV_tot_n)
	{

		loadCSV_2col(PathVar::data + nameneg, OCV_neg_x, OCV_neg_y);			// the OCV curve of the anode, the first column gives the lithium fractions (increasing), the 2nd column gives the OCV vs li/li+
		loadCSV_2col(PathVar::data + namepos, OCV_pos_x, OCV_pos_y);			// the OCV curve of the cathode, the first column gives the lithium fractions (increasing), the 2nd column gives the OCV vs li/li+
		loadCSV_2col(PathVar::data + nameentropicC, dOCV_neg_x, dOCV_neg_y);	// the entropic coefficient of the anode, the first column gives the lithium fractions (increasing), the 2nd column gives the entropic coefficient [V K-1]
		loadCSV_2col(PathVar::data + nameentropicCell, dOCV_tot_x, dOCV_tot_y); // the entropic coefficient of the entire cell, the first column gives the lithium fractions (increasing), the 2nd column gives the entropic coefficient [V K-1]

		// Check if it is fixed step.
		is_OCV_neg_fixed = check_is_fixed(OCV_neg_x);
		is_OCV_pos_fixed = check_is_fixed(OCV_pos_x);
		is_dOCV_neg_fixed = check_is_fixed(dOCV_neg_x);
		is_dOCV_tot_fixed = check_is_fixed(dOCV_tot_x);
	}

	double linInt_OCV_pos(double x, bool print = false, bool bound = true)
	{
		return ::linInt(print, bound, OCV_pos_x, OCV_pos_y, OCV_pos_n, x, is_OCV_pos_fixed);
	}
	double linInt_OCV_neg(double x, bool print = false, bool bound = true)
	{
		return ::linInt(print, bound, OCV_neg_x, OCV_neg_y, OCV_neg_n, x, is_OCV_neg_fixed);
	}
	double linInt_dOCV_neg(double x, bool print = false, bool bound = true)
	{
		return ::linInt(print, bound, dOCV_neg_x, dOCV_neg_y, dOCV_neg_n, x, is_dOCV_neg_fixed);
	}
	double linInt_dOCV_tot(double x, bool print = false, bool bound = true)
	{
		return ::linInt(print, bound, dOCV_tot_x, dOCV_tot_y, dOCV_tot_n, x, is_dOCV_tot_fixed);
	}

	OCVcurves() = default;
};

struct OCVcurves_array
{
	int OCV_pos_n, OCV_neg_n;	// number of data points in the OCV curve for the cathode/anode
	int dOCV_neg_n, dOCV_tot_n; // number of data points in the entropic coefficient for the anode/entire cell

	bool is_OCV_pos_fixed{false}, is_OCV_neg_fixed{false}, is_dOCV_neg_fixed{false}, is_dOCV_tot_fixed{false};

	std::string namepos, nameneg, nameentropicC, nameentropicCell;

	std::vector<double> OCV_pos_x, OCV_neg_x;	// lithium fractions of the points of the cathode/anode OCV curve
	std::vector<double> OCV_pos_y, OCV_neg_y;	// voltage vs li/li+ of the points of the cathode/anode OCV curve [V]
	std::vector<double> dOCV_neg_x, dOCV_tot_x; // lithium fractions of the points of the anode/entire cell entropic coefficient curve
	std::vector<double> dOCV_neg_y, dOCV_tot_y; // entropic coefficient curve / the entire cell's entropic coefficient [V K-1]

	OCVcurves_array(const std::string &_namepos, int _OCV_pos_n,
					const std::string &_nameneg, int _OCV_neg_n,
					const std::string &_nameentropicC, int _dOCV_neg_n,
					const std::string &_nameentropicCell, int _dOCV_tot_n)
		: OCV_pos_n(_OCV_pos_n), OCV_neg_n(_OCV_neg_n), dOCV_neg_n(_dOCV_neg_n), dOCV_tot_n(_dOCV_tot_n),
		  namepos(_namepos), nameneg(_nameneg), nameentropicC(_nameentropicC), nameentropicCell(_nameentropicCell),
		  OCV_pos_x(_OCV_pos_n), OCV_neg_x(_OCV_neg_n), OCV_pos_y(_OCV_pos_n), OCV_neg_y(_OCV_neg_n),
		  dOCV_neg_x(_dOCV_neg_n), dOCV_tot_x(_dOCV_tot_n), dOCV_neg_y(_dOCV_neg_n), dOCV_tot_y(_dOCV_tot_n)
	{

		loadCSV_2col(PathVar::data + nameneg, OCV_neg_x, OCV_neg_y);			// the OCV curve of the anode, the first column gives the lithium fractions (increasing), the 2nd column gives the OCV vs li/li+
		loadCSV_2col(PathVar::data + namepos, OCV_pos_x, OCV_pos_y);			// the OCV curve of the cathode, the first column gives the lithium fractions (increasing), the 2nd column gives the OCV vs li/li+
		loadCSV_2col(PathVar::data + nameentropicC, dOCV_neg_x, dOCV_neg_y);	// the entropic coefficient of the anode, the first column gives the lithium fractions (increasing), the 2nd column gives the entropic coefficient [V K-1]
		loadCSV_2col(PathVar::data + nameentropicCell, dOCV_tot_x, dOCV_tot_y); // the entropic coefficient of the entire cell, the first column gives the lithium fractions (increasing), the 2nd column gives the entropic coefficient [V K-1]

		// Check if it is fixed step.
		is_OCV_neg_fixed = check_is_fixed(OCV_neg_x);
		is_OCV_pos_fixed = check_is_fixed(OCV_pos_x);
		is_dOCV_neg_fixed = check_is_fixed(dOCV_neg_x);
		is_dOCV_tot_fixed = check_is_fixed(dOCV_tot_x);
	}

	double linInt_OCV_pos(double x, bool print = false, bool bound = true)
	{
		return ::linInt(print, bound, OCV_pos_x, OCV_pos_y, OCV_pos_n, x, is_OCV_pos_fixed);
	}
	double linInt_OCV_neg(double x, bool print = false, bool bound = true)
	{
		return ::linInt(print, bound, OCV_neg_x, OCV_neg_y, OCV_neg_n, x, is_OCV_neg_fixed);
	}
	double linInt_dOCV_neg(double x, bool print = false, bool bound = true)
	{
		return ::linInt(print, bound, dOCV_neg_x, dOCV_neg_y, dOCV_neg_n, x, is_dOCV_neg_fixed);
	}
	double linInt_dOCV_tot(double x, bool print = false, bool bound = true)
	{
		return ::linInt(print, bound, dOCV_tot_x, dOCV_tot_y, dOCV_tot_n, x, is_dOCV_tot_fixed);
	}

	OCVcurves_array()
	{
		std::cout << "Default constructor is called" << std::endl;
	};
};