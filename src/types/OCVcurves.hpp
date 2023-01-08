/*
 * OCV_curves.hpp
 *
 *  Created on: 27 Jun 2022
 *   Author(s): Jorn Reniers, Volkan Kumtepeli
 */

#pragma once

#include "XYdata.hpp"
#include "../utility/slide_aux.hpp"
#include "../settings/settings.hpp"

namespace slide {
struct OCVcurves
{
protected:
public:
  slide::XYdata_ss OCV_pos, OCV_neg;   //!< voltage vs lithium fractions li/li+ of the points of the cathode/anode OCV curve [V]
  slide::XYdata_ss dOCV_neg, dOCV_tot; //!< entropic coefficient curve / the entire cell's entropic coefficient vs li fractions  [V K-1]

  OCVcurves(const std::string &_namepos, const std::string &_nameneg,
            const std::string &_nameentropicC, const std::string &_nameentropicCell)
  {

    OCV_neg.setCurve(PathVar::data / _nameneg);           //!< the OCV curve of the anode, the first column gives the lithium fractions (increasing), the 2nd column gives the OCV vs li/li+
    OCV_pos.setCurve(PathVar::data / _namepos);           //!< the OCV curve of the cathode, the first column gives the lithium fractions (increasing), the 2nd column gives the OCV vs li/li+
    dOCV_neg.setCurve(PathVar::data / _nameentropicC);    //!< the entropic coefficient of the anode, the first column gives the lithium fractions (increasing), the 2nd column gives the entropic coefficient [V K-1]
    dOCV_tot.setCurve(PathVar::data / _nameentropicCell); //!< the entropic coefficient of the entire cell, the first column gives the lithium fractions (increasing), the 2nd column gives the entropic coefficient [V K-1]
  }

  OCVcurves() = default;

  static OCVcurves makeOCVcurves(cellType tp)
  {
    if (tp != cellType::KokamNMC) {
      std::cerr << "NOT IMPLEMENTED cellType\n";
      throw "NOT IMPLEMENTED cellType";
    }

    OCVcurves M(settings::path::Kokam::namepos, settings::path::Kokam::nameneg, settings::path::Kokam::nameentropicC, settings::path::Kokam::nameentropicCell);
    return M;
  }

  //!< static OCVcurves *makeOCVcurves(std::string _namepos = settings::path::Kokam::namepos,
  //!<                                 std::string _nameneg = settings::path::Kokam::namepos,
  //!<                                 std::string _nameentropicC = settings::path::Kokam::nameentropicC,
  //!<                                 std::string _nameentropicCell = settings::path::Kokam::nameentropicCell)

  //!< {
  //!< }
};

} // namespace slide