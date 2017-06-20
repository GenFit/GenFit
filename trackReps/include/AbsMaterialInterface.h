/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

#ifndef genfit_AbsMaterialInterface_h
#define genfit_AbsMaterialInterface_h

#include "RKTrackRep.h"
#include "Material.h"

#include <TObject.h>
#include <TVector3.h>


namespace genfit {

/**
 * @brief Abstract base class for geometry interfacing
 */
class AbsMaterialInterface : public TObject {

 public:

  AbsMaterialInterface() : debugLvl_(0) {;};
  virtual ~AbsMaterialInterface(){;};

  /** @brief Initialize the navigator at given position and with given direction.  Return true if volume changed.
   */
  virtual bool initTrack(double posX, double posY, double posZ,
                         double dirX, double dirY, double dirZ) = 0;

  /***
   * Get the material paramaters in the current material.
   * @return
   */
  virtual Material getMaterialParameters() = 0;

  /** @brief Make a step until maxStep or the next boundary is reached.
   *
   * After making a step to a boundary, the position has to be beyond the boundary,
   * i.e. in the current material has to be that beyond the boundary.
   * The actual step made is returned.
   */
  virtual double findNextBoundary(const RKTrackRep* rep,
                                  const M1x7& state7,
                                  double sMax,
                                  bool varField = true) = 0;

  virtual void setDebugLvl(unsigned int lvl = 1) {debugLvl_ = lvl;}

 protected:
  unsigned int debugLvl_;

  //ClassDef(AbsMaterialInterface, 1);

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_AbsMaterialInterface_h
