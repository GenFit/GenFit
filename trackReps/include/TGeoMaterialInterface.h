/* Copyright 2008-2014, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

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

/** @addtogroup RKTrackRep
 * @{
 */

#ifndef genfit_TGeoMaterialInterface_h
#define genfit_TGeoMaterialInterface_h

#include "AbsMaterialInterface.h"


namespace genfit {

/**
 * @brief AbsMaterialInterface implementation for use with ROOT's TGeoManager.
 */
class TGeoMaterialInterface : public AbsMaterialInterface {

 public:

  TGeoMaterialInterface() {};
  ~TGeoMaterialInterface(){;};

  /** @brief Initialize the navigator at given position and with given
      direction.  Returns true if the volume changed.
   */
  bool initTrack(double posX, double posY, double posZ,
                 double dirX, double dirY, double dirZ);

  Material getMaterialParameters() override;

  /** @brief Make a step (following the curvature) until step length
   * sMax or the next boundary is reached.  After making a step to a
   * boundary, the position has to be beyond the boundary, i.e. the
   * current material has to be that beyond the boundary.  The actual
   * step made is returned.
   */
  double findNextBoundary(const RKTrackRep* rep,
                          const M1x7& state7,
                          double sMax,
                          bool varField = true);

  // ClassDef(TGeoMaterialInterface, 1);

 private:
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_TGeoMaterialInterface_h
