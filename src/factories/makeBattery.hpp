/**
 * @file makeBattery.hpp
 * @brief Battery factory function
 * @author Jorn Reniers
 * @author Volkan Kumtepeli
 * @date 02 Mar 2020
 */

#pragma once

#include "../StorageUnit.hpp"

#include <string>
#include <memory>

namespace slide {
Deep_ptr<StorageUnit> makeBattery(bool balance, bool capSpread, bool RcellSpread, bool degrateSpread, bool contactR, int coolControl);
Deep_ptr<StorageUnit> makeBattery2(bool balance, bool capSpread, bool RcellSpread, bool degrateSpread, bool contactR, int coolControl);
Deep_ptr<StorageUnit> makeBattery_EPFL(bool capSpread, bool RcellSpread, bool degrateSpread, bool contactR, int coolControl, std::string IDadditions, double RM);
Deep_ptr<StorageUnit> makeBattery_Test(bool capSpread, bool RcellSpread, bool degrateSpread, bool contactR, int coolControl, std::string IDadditions, double RM);
Deep_ptr<StorageUnit> makeBattery_TestParallel(bool capSpread, bool RcellSpread, bool degrateSpread, bool contactR, int coolControl, std::string IDadditions, double RM);
Deep_ptr<StorageUnit> makeBattery_EPFL_smaller(bool capSpread, bool RcellSpread, bool degrateSpread, bool contactR, int coolControl, std::string IDadditions, double RM);
} // namespace slide