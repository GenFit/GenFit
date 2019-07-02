/**************************************************************************
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2019 - Belle II Collaboration                             *
 *                                                                        *
 * Author: The Belle II Collaboration                                     *
 * Contributors: Dmitrii Neverov                                          *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

#include <genfit/RKTrackRep.h>

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
