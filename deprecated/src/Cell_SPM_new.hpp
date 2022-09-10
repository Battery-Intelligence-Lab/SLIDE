/*
 * CellSPM.hpp
 *
 *  Created on: 8 Feb 2020
 *   Author(s): Jorn Reniers
 */

#ifndef SRC_CELL_SPM_HPP_
#define SRC_CELL_SPM_HPP_

#include "Cell.hpp"
#include "SPMModel.h"
#include "global.h"

#include <ctime>

using namespace std;
//!< Define a structure with the identifications of which degradation model(s) to use
struct DEG_ID
{
  //!< identifiers for degradation models
  //!< Each array is made with a length of 'len', which is the maximum number of models per mechanisms.
  //!< If the user wants to use more models, you only have to increase the value of 'len' to the number you want to use
  //!< and change the '10' in the definition of the arrays
  int SEI_id[CELL_NDEG]; //!< array with identifications to decide which SEI models to use. Max length 10
                         /* 				0	no SEI growth
                          * 				1 	kinetic model only (Tafel kinetics)
                          * 				2 	Pinson&Bazant model: linear diffusion + Tafel kinetics
                          * 				3	Christensen and Newman model
                          * 				4 	model from paper on optimisation (Christ & Newman without Rdrop and with knt/kpt)
                          */
  int SEI_n;             //!< SEI_N 	number of SEI models to use (length of SEI_ID)
  int SEI_porosity;      //!< integer deciding whether we reduce the active volume fraction due to SEI growth
                         /* 				0	don't reduce it
                          * 				1	use correlation from Ashwin et al. 2016
                          */
  int CS_id[CELL_NDEG];  //!< array with identifications for which model to use for surface cracking. Max length 10
                         /* 				0 	no surface cracking
                          * 				1 	Laresgoiti's stress + crack growth model
                          * 				2 	Dai stress model + Laresgoiti crack growth
                          * 				3 	model based on Deshpande and Bernardi, 2017
                          * 				4 	model from Barai et al
                          * 				5 	model from Ekstrom et al
                          */
  int CS_n;              //!< number of surface crack growth models to use (length of CS_ID)
  int CS_diffusion;      //!< integer deciding whether we reduce the negative diffusion constant due to surface cracks
                         /* 				0 	don't decrease diffusion
                          * 				1	decrease according to Barai et al. 2015
                          */
  int LAM_id[CELL_NDEG]; //!< array with the integers deciding which models is to be used for loss of active material. Max length 10
                         /* 				0 	no LAM
                          * 				1	Dai's stress model and Laresgoiti's correlation to get LAM
                          * 				2	delacourt's	correlation between abs(j) and porosity
                          * 				3 	Kindermann's model for cathode dissolution: tafel kinetics for increased porosity
                          * 				4 	Narayanrao's correlation which decreases the effective surface area proportionally to itself and j
                          */
  int LAM_n;             //!< number of LAM models to be used (length of LAM_id)
  int pl_id;             //!< integer deciding which model is to be used for li-plating
                         /* 				0 	no plating
                          * 				1	Yang et al thermodynamic plating (Tafel kinetics)
                          */
};

//!< Define a structure with the fitting parameters of the SEI growth models (SEI)
struct SEIparam
{

  double sei1k;   //!< rate parameter of the SEI side reaction in the 1st SEI model
  double sei1k_T; //!< activation energy of sei1k

  double sei2k;   //!< rate parameter of the SEI side reaction in the 2nd SEI model
  double sei2k_T; //!< activation energy of sei2k
  double sei2D;   //!< diffusion constant of the SEI layer in the 2nd SEI model
  double sei2D_T; //!< activation energy of sei2D

  double sei3k;   //!< rate parameter of the SEI side reaction in the 3rd SEI model
  double sei3k_T; //!< activation energy of sei3k
  double sei3D;   //!< diffusion constant of the SEI layer in the 3rd SEI model
  double sei3D_T; //!< activation energy of sei3D

  double sei4k;   //!< rate parameter of the SEI side reaction in the 3rd SEI model
  double sei4k_T; //!< activation energy of sei3k
  double sei4D;   //!< diffusion constant of the SEI layer in the 3rd SEI model
  double sei4D_T; //!< activation energy of sei3D

  double sei_porosity; //!< proportionality constant between the SEI growth and the decrease in volume fraction of active material
                       //!< (because the SEI layer blocks the pores at the surface)
};

//!< Define a structure with the fitting parameters of the surface crack growth models (CS)
struct CSparam
{

  double CS1alpha; //!< fitting parameter of the 1st surface crack growth model

  double CS2alpha; //!< fitting parameter of the 2nd surface crack growth model

  double CS3alpha; //!< fitting parameter of the 3rd surface crack growth model

  double CS4Amax;  //!< maximum crack growth surface for the 4th surface crack growth model
  double CS4alpha; //!< fitting parameter of the 4th surface crack growth model

  double CS5k;   //!< rate parameter of the 5th surface crack growth model at reference temperature
  double CS5k_T; //!< activation energy of CS5k

  double CS_diffusion; //!< fitting parameter to decrease the diffusion constant due to surface cracks
};

//!< Define a structure with the fitting parameters of the models for loss of active material (LAM)
struct LAMparam
{
  double lam1p; //!< fitting parameter for the positive electrode for the 1st LAM model
  double lam1n; //!< fitting parameter for the negative electrode for the 1st LAM model

  double lam2ap; //!< fitting parameter 1 at reference temperature for the positive electrode for the 2nd LAM model
  double lam2bp; //!< fitting parameter 2 at reference temperature for the positive electrode for the 2nd LAM model
  double lam2an; //!< fitting parameter 1 at reference temperature for the negative electrode for the 2nd LAM model
  double lam2bn; //!< fitting parameter 2 at reference temperature for the negative electrode for the 2nd LAM model
  double lam2t;  //!< activation energy for all the parameters of the 2nd LAM model

  double lam3k;   //!< rate constant at reference temperature for the cathode dissolution side reaction
  double lam3k_T; //!< activation energy for lam3k

  double lam4p; //!< fitting parameter for the positive electrode for the 4th LAM model
  double lam4n; //!< fitting parameter for the negative electrode for the 4th LAM model
};

//!< Define a structure with the fitting parameters of the li-plating models (PL)
struct PLparam
{
  double pl1k;   //!< rate constant of the li-plating side reaction at reference temperature in the 1st model
  double pl1k_T; //!< activation energy of pl1k
};

class Cell_SPM : public Cell
{
protected:
  //!< state, inherits
  //!< SoC
  //!< I
  //!< T
  //!< nstates
  double state[2 * CELL_NCH + 15]; //!< array with all states
  /*
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
   */
  struct DEG_ID deg_id;             //!< structure with the identification of which degradation model(s) to use
  double Therm_Qgen;                //!< total heat generation since the last update [J]
  double Therm_Qgentot;             //!< variable for unit testing, total heat generation since the beginning of this cell's life [J]
  double Therm_time;                //!< time since the last update of the thermal model
  double degState[CELL_NSTATE_MAX]; //!< array with the degradation damage of all states
  double Vcell;                     //!< cell voltage at this point in time
  double etapcell;
  double etancell;
  bool etacell_valid;
  bool Vcell_valid; //!< is the value stored in Vcell valid at this point in time?

  //!< cell to cell variations
  double var_cap;    //!< relative factor increasing the capacity of the cell
  double var_R;      //!< relative factor increasing the DC resistance
  double var_degSEI; //!< relative factor to speed up or slow down the rate of SEI growth
  double var_degLAM; //!< relative factor to speed up or slow down the rate of LAM

  //!< Battery model constants
  double Cmaxpos;     //!< maximum lithium concentration in the cathode [mol m-3]
  double Cmaxneg;     //!< maximum lithium concentration in the anode [mol m-3]
  double C_elec_sqrt; //!< square root of the Li- concentration in electrolyte [mol m-3]
  double n;           //!< number of electrons involved in the main reaction [-]
  double F;           //!< Faraday's constant
  double Rg;          //!< ideal gas constant
  double T_ref;       //!< reference temperature [K]

  //!< Thermal model parameters
  double T_env; //!< environment temperature [K]
  double Qch;   //!< convective heat transfer coefficient for a stand-alone cell with convective cooling [W K-1 m-3]
  double rho;   //!< density of the battery
  double Cp;    //!< thermal capacity of the battery

  //!< Geometric parameters
  double L;         //!< thickness of the cell [m]
  double Acell;     //!< geometric surface area of the cell [m2]
  double elec_surf; //!< geometric surface area of the electrodes (electrode height * electrode width) [m2]
  double Rp;        //!< radius of the positive sphere of the Single Particle model [m]
  double Rn;        //!< radius of the negative sphere of the Single Particle model [m]
  //!< other geometric parameters are part of State because they can change over the battery's lifetime

  //!< parameters of the main li-insertion reaction
  double kp;   //!< rate constant of main reaction at positive electrode at reference temperature
  double kp_T; //!< activation energy for the Arrhenius relation of kp
  double kn;   //!< rate constant of main reaction at negative electrode at reference temperature
  double kn_T; //!< activation energy for the Arrhenius relation of kn
  //!< The diffusion constants at reference temperature are part of State because they can change over the battery's lifetime
  double Dp_T;       //!< activation energy for the Arrhenius relation of Dp
  double Dn_T;       //!< activation energy for the Arrhenius relation of Dn
  struct SPMModel M; //!< Matrices for spatial discretisation of the solid diffusion model

  //!< OCV curves
  int OCV_pos_n;                   //!< number of data points in the OCV curve for the cathode
  double OCV_pos_x[CELL_NOCV_MAX]; //!< lithium fractions of the points of the cathode OCV curve
  double OCV_pos_y[CELL_NOCV_MAX]; //!< voltage vs li/li+ of the points of the cathode OCV curve [V]

  int OCV_neg_n;                   //!< number of data points in the OCV curve for the anode
  double OCV_neg_x[CELL_NOCV_MAX]; //!< lithium fractions of the points of the anode OCV curve
  double OCV_neg_y[CELL_NOCV_MAX]; //!< voltage vs li/li+ of the points of the anode OCV curve [V]

  int dOCV_neg_n;                   //!< number of data points in the entropic coefficient for the anode
  double dOCV_neg_x[CELL_NOCV_MAX]; //!< lithium fractions of the points of the anode entropic coefficient curve
  double dOCV_neg_y[CELL_NOCV_MAX]; //!< entropic coefficient curve [V K-1]

  int dOCV_tot_n;                   //!< number of data points in the entropic coefficient for the entire cell OCV curve
  double dOCV_tot_x[CELL_NOCV_MAX]; //!< cathodic lithium fractions of the points of the entire cell's entropic coefficient
  double dOCV_tot_y[CELL_NOCV_MAX]; //!< the entire cell's entropic coefficient [V K-1]

  //!< SEI parameters
  double rsei;              //!< specific resistance times real surface area of the SEI film [Ohm m]
  double nsei;              //!< number of electrons involved in the SEI reaction [-]
  double alphasei;          //!< charge transfer coefficient of the SEI reaction [-]
  double OCVsei;            //!< equilibrium potential of the SEI side reaction [V]
  double rhosei;            //!< partial molar volume of the SEI layer [m3 mol-1]
  double c_elec0;           //!< bulk concentration of the electrolyte molecule participating in the SEI growth (e.g. EC) [mol m-3]
  double Vmain;             //!< partial molar volume of the main reaction, see Ashwin et al, 2016
  double Vsei;              //!< partial molar volume of the SEI side reaction, see Ashwin et al., 2016
  struct SEIparam seiparam; //!< structure with the fitting parameters of the different SEI growth models

  //!< surface crack parameters & constants
  double omegap; //!< partial molar volume of positive electrode [m3 mol-1]
  double omegan; //!< partial molar volume of negative electrode [m3 mol-1]
  double Ep;     //!< Young's modulus of positive electrode [GPa]
  double En;     //!< Young's modulus of negative electrode [GPa]
  double nup;    //!< Poisson's ratio of positive electrode [-]
  double nun;    //!< Poisson's ratio of negative electrode [-]
  //!< values of the stress are often needed. Because it takes very long to calculate them, we calculate them once and store them so we don't need to repeat the same calculation twice
  bool s_dai;             //!< do we need to calculate the stress according to Dai's model?
  bool s_lares;           //!< do we need to calculate the stress according to Laresgoiti's model?
  bool s_dai_update;      //!< boolean to indicate if Dai's stress are up to date with the battery state at this time step
  bool s_lares_update;    //!< boolean to indicate if Dai's stress are up to date with the battery state at this time step
  double s_dai_p;         //!< maximum hydrostatic stress in the positive particle according to Dai's stress model
  double s_dai_n;         //!< maximum hydrostatic stress in the negative particle according to Dai's stress model
  double s_lares_n;       //!< stress in the negative particle according to Laresgoiti's stress model
  double s_dai_p_prev;    //!< maximum hydrostatic stress in the previous time step in the positive particle according to Dai's stress model
  double s_dai_n_prev;    //!< maximum hydrostatic stress in the previous time step in the negative particle according to Dai's stress model
  double s_lares_n_prev;  //!< stress in the previous time step in the negative particle according to Laresgoiti's stress model
  double s_dt;            //!< time period between the 'previous' and 'current' stress [s]
  struct CSparam csparam; //!< structure with the fitting parameters of the different crack growth models

  //!< LAM parameters & constants
  double OCVnmc;            //!< equilibrium potential of the NMC dissolution side reaction [V]
  struct LAMparam lamparam; //!< structure with the fitting parameters of the different LAM models

  //!< Li-plating parameters & constants
  double npl;             //!< number of electrons involved in the plating reaction [-]
  double alphapl;         //!< charge transfer constant for the plating reaction [-]
  double OCVpl;           //!< OCV of the plating reaction [V]
  double rhopl;           //!< density of the plated lithium layer
  struct PLparam plparam; //!< structure with the fitting parameters of the different plating models

  //!< timing
#if TIMING
                          //!< std::clock_t tstart;
  double T_dstate;
  double T_getV;
  double T_getOCV;
  double T_validStates;
  double T_setStates;
#endif

//!< data storage
#if DATASTORE_CELL > 1
  double LLI[CELL_NDATA_MAX];        //!< total lost lithium inventory [Ah]
  double LAM_en[CELL_NDATA_MAX];     //!< volume fraction of active material on the anode [-], lost due to SEI pore clogging & Delacourt LAM
  double LAM_thickn[CELL_NDATA_MAX]; //!< thickness of the anode [m], lost due to Dai LAM
  double Rtot[CELL_NDATA_MAX];       //!< total resistance [Ohm]
#endif

  void getConcentration(int nin, double cp[], double cn[]);                                                                                                                                   //!< get the concentration at every node
  void setC(double fp0, double fn0);                                                                                                                                                          //!< set the concentrations to the given (uniform) concentration
  virtual void getDaiStress(int n, double *sigma_p, double *sigma_n, double sigma_r_p[], double sigma_r_n[], double sigma_t_p[], double sigma_t_n[], double sigma_h_p[], double sigma_h_n[]); //!< get the stresses at all nodes according to Dai's stress model
  virtual void updateDaiStress();                                                                                                                                                             //!< updated the stored stress values for Dai's stress model
  virtual void getLaresgoitiStress(bool print, double *sigma_n);                                                                                                                              //!< get the stresses at all nodes according to Laresgoiti's stress model
  virtual void updateLaresgoitiStress(bool print);                                                                                                                                            //!< update the stored stress values for Laresgoiti's stress model

  //!< thermal model
  double thermalModel_cell();
  double thermalModel_coupled(int Nneighb, double Tneighb[], double Kneighb[], double Aneighb[], double tim);

  //!< degradation functions
  virtual void SEI(double OCVnt, double etan, double *isei, double *den);                                                                             //!< calculate the effect of SEI growth
  virtual void CS(double OCVnt, double etan, double *isei_multiplyer, double *dCS, double *dDn);                                                      //!< calculate the effect of surface crack growth
  virtual void LAM(bool critical, double zp_surf, double etap, double *dthickp, double *dthickn, double *dap, double *dan, double *dep, double *den); //!< calculate the effect of LAM
  virtual void LiPlating(double OCVnt, double etan, double *isei);                                                                                    //!< calculate the effect of lithium plating

  //!< state space model
  void dState_diffusion(bool print, double dzp[], double dzn[], double &dSoC); //!< just diffusion PDE
  void dState_thermal(bool print, double &dQgen);                              //!< calculate the heat generation
  void dState_degradation(bool print, double dstates[]);                       //!< degradation effects

public:
  Cell_SPM();
  Cell_SPM(string IDi, const DEG_ID &degid, double capf, double resf, double degfsei, double degflam);
  virtual ~Cell_SPM();

  //!< basic getters and setters

  //!< Basic getters and setters for the states (public for unit test)
  double *getZp();
  double *getZn();
  double getCp_surface(double *zpi);
  double getCn_surface(double *zni);
  double getDelta();
  double getCS();
  double getDelta_pl();
  double getLLI();
  double getThickp();
  double getThickn();
  double getEp();
  double getEn();
  double getAp();
  double getAn();
  double getDp();
  double getDn();
  double getrdcp();
  double getrdcn();
  double getrdccc();
  double getRdc();
  void getVariations(double var[], int nin, int &nout);

  void setZp(double in[]);
  void setZn(double in[]);
  void setdZp(double din[]);
  void setdZn(double din[]);
  void setDelta(double in);
  void setCS(double in);
  void setDelta_pl(double in);
  void setLLI(double in);
  void setThickp(double in);
  void setThickn(double in);
  void setEp(double in);
  void setEn(double in);
  void setAp(double in);
  void setAn(double in);
  void setDp(double in);
  void setDn(double in);
  void setrdcp(double in);
  void setrdcn(double in);
  void setrdccc(double in);
  virtual double setSoC(double SoCnew, bool checkV = true, bool print = true);
  double setI(double Inew, bool checkV = true, bool print = true);

  //!< overwrite from Cell
  double getRtot();
  void getStates(double s[], int nin, int &nout);
  double getOCV(bool print = true);                         //!< print is an optional argument
  double getOCV(double fps, double fns, bool print = true); //!< print is an optional argument
  double getV(bool print = true);

  double setStates(double s[], int nin, bool checkV = true, bool print = true);
  bool validStates(double s[], int nin, bool print = true);
  void timeStep_CC(double dt, int steps = 1);

  //!< data storage
  virtual void storeData();
  virtual void writeData(string prefix);

  //!< thermal model
  double getThermalSurface();
  double thermalModel(int Nneighb, double Tneighb[], double Kneighb[], double Aneighb[], double tim);
  double thermal_getTotalHeat(); //!< function for unit testing

  //!< copy this object
  virtual shared_ptr<StorageUnit> copy();

  //!< keep timing record of how long all functions take
  void getTimings(double &gOCV, double &tV, double &tdstate, double &tvalidstate, double &tsetState);
  void setTimings(double gOCV, double tV, double tdstate, double tvalidstate, double tsetState);
};

#endif /* SRC_CELL_SPM_HPP_ */
