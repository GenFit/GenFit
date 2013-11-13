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

#ifndef genfit_ProlateSpacepointMeasurement_h
#define genfit_ProlateSpacepointMeasurement_h

#include "SpacepointMeasurement.h"


namespace genfit {

/** @brief Class for measurements implementing a space point hit geometry with a very prolate
 * form of the covariance matrix.
 *
 *  @author Johannes Rauch (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * Measurements from detectors measuring 3D space points with errors in one direction
 * much larger than the errors perpendicular should use this class.
 *
 * For these hits, a virtual detector plane lying in the POCA and
 * perpendicular to the track yields wrong results. Instead, the plane should contain the
 * direction of the largest error.
 *
 * The largest error direction can be set. Standard is in z.
 *
 */
class ProlateSpacepointMeasurement : public SpacepointMeasurement {

 public:
  ProlateSpacepointMeasurement(int nDim = 3);
  ProlateSpacepointMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint);

  virtual ~ProlateSpacepointMeasurement() {;}

  virtual AbsMeasurement* clone() const {return new ProlateSpacepointMeasurement(*this);}

  virtual SharedPlanePtr constructPlane(const StateOnPlane& state) const;


  const TVector3& getLargestErrorDirection(){return largestErrorDirection_;}
  void setLargestErrorDirection(const TVector3& dir){largestErrorDirection_ = dir.Unit();}

 protected:
  TVector3 largestErrorDirection_; // direction of largest error

 public:

  ClassDef(ProlateSpacepointMeasurement,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_ProlateSpacepointMeasurement_h
