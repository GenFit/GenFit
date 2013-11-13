/* Copyright 2008-2010, Technische Universitaet Muenchen,
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

#ifndef genfit_mySpacepointMeasurement_h
#define genfit_mySpacepointMeasurement_h

#include "SpacepointMeasurement.h"
#include "TrackCandHit.h"
#include "mySpacepointDetectorHit.h"


namespace genfit {

/** @brief Example class for a spacepoint measurement which can be created
 * from mySpacepointDetectorHit via the MeasurementFactory.
 *
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 */
class mySpacepointMeasurement : public SpacepointMeasurement {

 public:

  /** Default constructor for ROOT IO. */
  mySpacepointMeasurement() :
     SpacepointMeasurement() {;}

  mySpacepointMeasurement(const mySpacepointDetectorHit* detHit, const TrackCandHit* hit) :
    SpacepointMeasurement()
  {
    rawHitCoords_(0) = detHit->getPos()(0);
    rawHitCoords_(1) = detHit->getPos()(1);
    rawHitCoords_(2) = detHit->getPos()(2);
    rawHitCov_ = detHit->getCov();
    detId_ = hit->getDetId();
    hitId_ = hit->getHitId();
  }

  virtual AbsMeasurement* clone() const {return new mySpacepointMeasurement(*this);}

  ClassDef(mySpacepointMeasurement,1)
};
/** @} */

} /* End of namespace genfit */

#endif // genfit_mySpacepointMeasurement_h
