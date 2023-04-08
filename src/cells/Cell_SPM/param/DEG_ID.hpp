/** \file DEG_ID.hpp
 *  \brief Defines a structure for handling degradation models.
 *  \author Jorn Reniers, Volkan Kumtepeli
 *  \date 12 May 2022
 */

#pragma once

#include <string>
#include <cstdint>
#include <array>
#include <iostream>

#include "../../../types/SmallVector.hpp"

//! Slide namespace contains all the types, classes, and functions for the simulation framework.
namespace slide {

/**
 * @brief DEG_ID structure handles the identifications of which degradation model(s) to use.
 */
struct DEG_ID
{
  using data_t = uint_fast8_t;      //!< Alias for uint_fast8_t as data_t.
  constexpr static data_t len = 10; //!< Length of the arrays with identifications of which models to use.

  /**
   * @brief DegArray class derived from SmallVector to handle degradation model arrays.
   */
  struct DegArray : public SmallVector<data_t, len>
  {
    /**
     * @brief Adds a model to the array.
     * @param elem Model identification number.
     */
    void add_model(data_t elem)
    {
      if (elem != 0)                //!< If element is 0 then there is no need to add.
        push_back(std::move(elem)); //!< #TODO better push_back needed.
    }
  };

  DegArray SEI_id{}; //!< Array with identifications to decide which SEI models to use.
  /**
   * - 0:	no SEI growth
   * - 1:	kinetic model only (Tafel kinetics), ref: Ning & Popov, Journal of the Electrochemical Society 151 (10), 2004
   * - 2:	Pinson&Bazant model: linear diffusion + Tafel kinetics, ref: Pinson & Bazant, Journal of the Electrochemical society 160 (2), 2013
   * - 3:	Christensen and Newman model, ref: Christensen & Newmann, Journal of the Electrochemical Society 152 (4), 2005
   */

  data_t SEI_porosity{ 0 }; //!< Integer deciding whether we reduce the active volume fraction due to SEI growth.
  /**
   * - 0:	don't reduce it
   * - 1:	use correlation from Ashwin et al. 2016, ref: Ashwin, Chung, Wang, Journal of Power Sources 328, 2016
   */

  DegArray CS_id; //!< Array with identifications for which model to use for surface cracking. Max length 10.
  /**
   * - 0:	no surface cracking
   * - 1:	Laresgoiti's stress + crack growth model, ref: Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015
   * - 2:	Dai stress model + Laresgoiti crack growth, ref: Laresgoiti, Kablitz, Ecker, Sauer, Journal of Power Sources 300, 2015; Dai, Cai, White, Journal of Power sources 247, 2014
   * - 3:	model based on Deshpande and Bernardi, ref: Deshpande & Bernardi, Journal of the Electrochemical Society 164 (2), 2017
   * - 4:	model from Barai et al, ref: Barai, Smith, Chen, Kim, Mukherjee, Journal of the Electrochemical Society 162 (9), 2015
   * - 5:	model from Ekstrom et al, ref: Ekstrom and Lindbergh, Journal of the Electrochemical Society 162 (6), 2015
   */

  data_t CS_diffusion{ 0 }; //!< Integer deciding whether we reduce the negative diffusion constant due to surface cracks.
  /**
   * - 0:	don't decrease diffusion
   * - 1:	decrease according to Barai et al. 2015
   */

  DegArray LAM_id; //!< Array with the integers deciding which models is to be used for loss of active material. Max length 10.
  /**
   * - 0:	no LAM
   * - 1:	Dai's stress model and Laresgoiti's correlation to get LAM
   * - 2:	delacourt's	correlation between abs(j) and porosity
   * - 3:	Kindermann's model for cathode dissolution: tafel kinetics for increased porosity
   * - 4:	Narayanrao's correlation which decreases the effective surface area proportionally to itself and j
   */

  data_t pl_id{ 0 }; //!< Integer deciding which model is to be used for li-plating.
  /**
   * - 0:	no plating
   * - 1:	Yang et al thermodynamic plating (Tafel kinetics)
   */

  auto print()
  {
    /**
     * Function to get a string representation of the struct with the degradation settings.
     * This string is part of the names of the subfolders in which results are written.
     *
     * @return std::string representation of the degradation identifiers.
     *         Identifiers of the same mechanism are separated by -.
     *         Identifiers of different mechanisms are separated by _.
     *         e.g., if we use SEI model 1, no SEI porosity effect, no surface cracks, LAM model 2 and LAM model 3 and lithium plating model 1:
     *         2-0_0-0_2-3_1
     *         2:		SEI model 1
     *         0:		no SEI porosity
     *         0:		no surface cracks
     *         0: 	don't decrease the diffusion due to surface cracks
     *         2:		LAM model 2
     *         3:		LAM model 3
     *         1:		lithium plating model 1
     */

    std::string id = "";
    id.reserve(30);

    // Print SEI models and SEI_porosity (decreasing the active volume fraction due to SEI growth), separated by -

    if (SEI_id.empty()) // If zero we do not store so we need to print even not storing.
      id += "0-";
    else
      for (const auto sei_id : SEI_id)
        id += std::to_string(sei_id) + '-';

    id += std::to_string(SEI_porosity);

    // Mechanism separator
    id += '_';

    // Print CS models (surface cracking) and CS_diffusion (reducing the diffusion constant because of cracks), separated by -

    if (CS_id.empty())
      id += "0-";
    else
      for (const auto cs_id : CS_id)
        id += std::to_string(cs_id) + '-';

    id += std::to_string(CS_diffusion);

    // Mechanism separator
    id += '_';

    // Print LAM models separated by -
    if (LAM_id.empty())
      id += "0-";
    else
      for (const auto lam_id : LAM_id)
        id += std::to_string(lam_id) + '-';

    // Mechanism separator
    if (LAM_id.empty())
      id += '_';
    else
      id.back() = '_'; // If LAM_id is not empty then it will print an extra '-' therefore we need to replace it with mechanism separator.

    // Print plating model
    id += std::to_string(pl_id);

    // Output
    return id;
  }
};
} // namespace slide
