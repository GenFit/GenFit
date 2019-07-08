/* Copyright 2019, Belle II Collaboration
   Authors: Dmitrii Neverov

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "RKTrackRep.h"

namespace genfit {
  /**
   * Monopole track representation.
   * It is a minimal modification of RKTrackRep
   * with a different equations of motion for magnetic charges.
   *
   * In the current implementation the states on plane are 5-d:
   * u, v, u', v', q/p
   * except that q in this case is magnetic, and the monopole
   * has no electic charge.
   */
  class MplTrackRep : public RKTrackRep {

  public:

    MplTrackRep() : m_magCharge(0), m_mass(0) {};
    MplTrackRep(int pdgCode, float magCharge, char propDir = 0);
    ~MplTrackRep();

    // Hopefully this is the only function that is vastly different for monopoles
    double RKPropagate(M1x7& state7,
                       M7x7* jacobian,
                       M1x3& SA,
                       double S,
                       bool varField = true,
                       bool calcOnlyLastRowOfJ = false) const override;

    // Returns the magnetic charge instead of electric as in the base class.
    double getCharge(const StateOnPlane& state) const override; 

  private:

    const double m_magCharge; // the magnitude of magnetic charge in units of e+
    const double m_mass; // the mass of the monopole in units of GeV/c^2


  public:

    ClassDefOverride(MplTrackRep, 1)

  };
} //end genfit namespace
