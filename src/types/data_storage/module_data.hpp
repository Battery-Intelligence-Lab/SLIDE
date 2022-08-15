/*
 * module_data.hpp
 *
 *  Created on: 16 May 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

namespace slide {
//!< struct ModuleData
//!< {
//!< 	double Ahtot;	//!< total charge throughput
//!< 	double Whtot;	//!< total energy throughput
//!< 	double Timetot; //!< total time [s]
//!< 	double Itot;	//!< total current [A]
//!< 	double Vtot;	//!< total voltage [V]
//!< 	double Ttot;	//!< temperature [K]

//!< 	//!< double I;						  //!< current of the cell at every step [A]
//!< 	//!< double V;						  //!< voltage of the cell at every step [V]
//!< 	//!< double OCVp;					      //!< cathode potential of the cell at every step [V]
//!< 	//!< double OCVn;					      //!< anode potential of the cell at every step [V]
//!< 	//!< double T;						  //!< temperature of the cell at every step [K]
//!< 	//!< double timeCha;					  //!< cumulative time spent on charging since the start at every step [s]
//!< 	//!< double AhCha;					  //!< cumulative charged throughput since the start at every step [A]
//!< 	//!< double WhCha;					  //!< cumulative charged energy throughput since the start at every step [Wh]
//!< 	//!< double timeDis;					  //!< cumulative time spent on discharging since the start at every step [s]
//!< 	//!< double AhDis;					  //!< cumulative discharged charge throughput since the start at every step [A]
//!< 	//!< double WhDis;					  //!< cumulative discharged energy throughput since the start at every step [Wh]
//!< 	//!< double timeRes;					  //!< cumulative time spent on rest since the start at every step [s]
//!< };

//!< template <int N>
//!< struct ModuleDataStorage
//!< {
//!< };

//!< template <>
//!< struct ModuleDataStorage<1>
//!< {
//!< 	double ahtot{};	  //!< charge throughput so far
//!< 	double whtot{};	  //!< energy throughput so far
//!< 	double timetot{}; //!< time so far [s]
//!< };

//!< template <>
//!< struct ModuleDataStorage<2>
//!< {
//!< 	double ahtot;	//!< charge throughput so far
//!< 	double whtot;	//!< energy throughput so far
//!< 	double timetot; //!< time so far [s]
//!< 	std::vector<ModuleData> moduleData;
//!< };
} // namespace slide