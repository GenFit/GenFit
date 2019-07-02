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
   * A prototype for monopole track representation.
   * For now it would be a minimal modification of RKTrackRep
   * with a different equation of motion for monopoles.
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

    double getCharge(const StateOnPlane& state) const override; 

  private:

    const double m_magCharge;
    const double m_mass;

  public:

    ClassDef(MplTrackRep, 1)

  };
} //end genfit namespace
