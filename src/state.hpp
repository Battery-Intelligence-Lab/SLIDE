/*
 * State.hpp
 *
 * Defines a class State which defines the state-variables of a cell for the state-space model formulation
 *
 * Copyright (c) 2019, The Chancellor, Masters and Scholars of the University
 * of Oxford, VITO nv, and the 'Slide' Developers.
 * See the licence file LICENCE.txt for more information.
 */

#pragma once

#include <array>
#include <utility>

#include "constants.hpp"

namespace slide
{
	using z_type = std::array<double, settings::nch>;
	using states_type = std::array<double, settings::ns>;
	using settings::nch;

	class State
	{

	public:
		State() = default; // Default constructor which DOESN'T initialise the states. All states are set to 0
		void initialise(z_type &zpi, z_type &zni, double Ti, double deltai, double LLIi,
						double thickpi, double thickni, double epi, double eni, double api, double ani,
						double CSi, double Dpi, double Dni, double Ri, double deltalii); // initialises all state variables to the given values

		auto const &getStates_arr() { return x; }

		auto const &get_zp(int i) { return x[i]; }		 // z_type, get the transformed concentrations
		auto const &get_zn(int i) { return x[nch + i]; } // z_type, get the transformed concentrations

		double &get_T() { return x[2 * nch + 0]; }		   // get the temperature [K]
		double &get_delta() { return x[2 * nch + 1]; }	   // get the SEI thickness
		double &get_LLI() { return x[2 * nch + 2]; }	   // get the lost lithium
		double &get_thickp() { return x[2 * nch + 3]; }	   // get the thickness of the cathode
		double &get_thickn() { return x[2 * nch + 4]; }	   // get the thickness of the anode
		double &get_ep() { return x[2 * nch + 5]; }		   // get the volume fraction of active material of the cathode
		double &get_en() { return x[2 * nch + 6]; }		   // get the volume fraction of active material of the anode
		double &get_ap() { return x[2 * nch + 7]; }		   // get the effective surface of the cathode
		double &get_an() { return x[2 * nch + 8]; }		   // get the effective surface of the anode
		double &get_CS() { return x[2 * nch + 9]; }		   // get the crack surface area
		double &get_Dp() { return x[2 * nch + 10]; }	   // get the diffusion constant of the cathode at reference temperature [m/s]
		double &get_Dn() { return x[2 * nch + 11]; }	   // get the diffusion constant of the anode at reference temperature [m/s]
		double &get_r() { return x[2 * nch + 12]; }		   // get the specific resistance (resistance times real surface area of the combined electrodes) [Ohm m2]
		double &get_delta_pl() { return x[2 * nch + 13]; } // get the thickness of the plated lithium layer

		//bool is_initialised() { return sini_ptr != nullptr; }

		double getR(double elec_surf) // get the total DC resistance of the electrodes
		{
			/*
			* IN
	 		* elec_surf geometric surface area of the electrodes [m2]
	 		*/
			const double surfp = get_ap() * get_thickp() * elec_surf; // real surface area of the positive electrode
			const double surfn = get_an() * get_thickn() * elec_surf; // real surface area of the negative electrode
			const double surf = (surfp + surfn) / 2;				  // mean surface area
			const double r_elec = get_r() / surf;					  // resistance of the electrodes [Ohm]
			return r_elec;
		}

		void setT(double Ti);				// set the temperature
		void setZ(z_type &zpi, z_type &zni) // set the transformed concentrations
		{

			std::copy(zpi.begin(), zpi.end(), x.begin());
			std::copy(zni.begin(), zni.end(), x.begin() + nch);
		}																											   // set the transformed concentration
		void setStates(slide::states_type &&states);																   // set the states to the values in the array
		void setIniStates(const slide::states_type &si);															   // set the initial states to the values in the array
		void overwriteGeometricStates(double thickpi, double thickni, double epi, double eni, double api, double ani); // overwrite the states related to the geometry of a cell
		void overwriteCharacterisationStates(double Dpi, double Dni, double ri);									   // overwrite the states related to the characterisation of a cell

	private:
		// battery states
		// zp[nch] zn[nch] T delta LLI thickp thickn ep en ap an CS Dp Dn R delta_liPlating
		slide::states_type x{}; // Array to hold all states.
								//	slide::State *sini_ptr{nullptr}; // array ptr with the initial battery states
	};
} // namespace slide