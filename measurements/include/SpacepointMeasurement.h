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
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_SpacepointMeasurement_h
#define genfit_SpacepointMeasurement_h

#include "AbsMeasurement.h"
#include "AbsHMatrix.h"


namespace genfit {

/** @brief Class for measurements implementing a space point hit geometry.
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * For a space point the detector plane has to be defined with respect to
 * a track representation. SpacepointMeasurement implements a scheme where the
 * detectorplane is chosen perpendicular to the track.
 * In a track fit, only two of the three coordinates of a space point are
 * independent (the track is a one-dimensional object). Therefore the 3D
 * data of the hit is used to define a proper detector plane into which the
 * hit coordinates are then projected.
 */
class SpacepointMeasurement : public AbsMeasurement {

 public:
  SpacepointMeasurement(int nDim = 3);
  SpacepointMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint);

  virtual ~SpacepointMeasurement() {;}

  virtual AbsMeasurement* clone() const {return new SpacepointMeasurement(*this);}

  virtual SharedPlanePtr constructPlane(const StateOnPlane& state) const;

  virtual std::vector<MeasurementOnPlane*> constructMeasurementsOnPlane(const StateOnPlane& state) const;

  virtual const AbsHMatrix* constructHMatrix(const AbsTrackRep*) const;

  ClassDef(SpacepointMeasurement,1)
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_SpacepointMeasurement_h
