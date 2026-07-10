/**
 * @file core_CellDesign_test.cpp
 * @brief Phase-1 tests for the physical description hierarchy (PLAN.md §3.3).
 */

#include "../../src/core/CellDesign.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace slide::core;

TEST_CASE("v4 Domain has explicit negative-first order and physical flux signs", "[core][design]")
{
  STATIC_REQUIRE(domain_index(Domain::neg) == 0);
  STATIC_REQUIRE(domain_index(Domain::pos) == 1);
  STATIC_REQUIRE(opposite(Domain::neg) == Domain::pos);
  STATIC_REQUIRE(opposite(Domain::pos) == Domain::neg);
  STATIC_REQUIRE(molar_flux_sign(Domain::neg) == 1);
  STATIC_REQUIRE(molar_flux_sign(Domain::pos) == -1);
}

TEST_CASE("CellDesign owns independent value-semantic electrodes", "[core][design]")
{
  CellDesign design;
  domain_value(design.electrode, Domain::neg).active_material.cs_max = 30'555.0;
  domain_value(design.electrode, Domain::pos).active_material.cs_max = 51'385.0;
  domain_value(design.electrode, Domain::neg).aging.push_back({ AgingMechanismKind::sei_kinetic, "SEI kinetic", { 1.0, 2.0 } });

  CellDesign copy = design;
  domain_value(copy.electrode, Domain::neg).active_material.cs_max = 31'000.0;
  domain_value(copy.electrode, Domain::neg).aging[0].coefficients[0] = 7.0;

  REQUIRE(domain_value(design.electrode, Domain::neg).active_material.cs_max == 30'555.0);
  REQUIRE(domain_value(design.electrode, Domain::pos).active_material.cs_max == 51'385.0);
  REQUIRE(domain_value(design.electrode, Domain::neg).aging[0].coefficients[0] == 1.0);
  REQUIRE(domain_value(copy.electrode, Domain::neg).active_material.cs_max == 31'000.0);
}
