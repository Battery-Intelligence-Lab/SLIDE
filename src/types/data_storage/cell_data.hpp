/*
 * cell_data.hpp
 *
 *  Created on: 07 Feb 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include <iostream>

#include "../Histogram.hpp"

namespace slide
{
    struct CellCumulativeData
    {
        double Ah{0}, Wh{0}, Time{0}; // charge throughput, energy throughput, time so far. [s]
    };

    struct ProcedureThroughputData
    {
        int ID;                // the ID of the (dis)charge
        double charge;         // the charge throughput of every (dis)charge
        double energy;         // the energy throughput of every (dis)charge
        double coolSystemLoad; // the energy [Wh] required to run the cooling system of the battery in every cycle
        double convloss;       // the energy losses in the power electronic converter [Wh] in every cycle
    };

    struct BasicData
    {
        using value_type = float;
        value_type I;   // current [A]
        value_type V;   // voltage [V]
        value_type SOC; // State of charge [0-1]
        value_type T;   // Temperature [K]
    };

    struct CommonData : public BasicData, public CellCumulativeData
    {
    };

    struct CellCommonHist
    {
        Histogram<> I, V, T; // histograms for current, voltage, temperature
    };

    struct BatteryData
    {
        double Timetot;  // total time [s]
        double Icells;   // total current to the battery compartment [A]
        double Vcells;   // total voltage of the battery compartment [V]
        double Tbatt;    // temperature of the battery container [K]
        double convloss; // losses in the power electronic converter [W]
    };

    struct TimingData_Module_s
    {
        double setCurrent{}, timeStep{}, timeStepi{}, validSUs{};
    };

    struct TimingData_Module_p
    {
        double redistributeCurrent{}, setCurrent{}, timeStep{}, timeStepi{}, validSUs{};
    };
}