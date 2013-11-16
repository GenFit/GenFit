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

#include "ProlateSpacepointMeasurement.h"

#include <cmath>

#include "Exception.h"
#include "RKTrackRep.h"


namespace genfit {

ProlateSpacepointMeasurement::ProlateSpacepointMeasurement(int nDim)
  : SpacepointMeasurement(nDim), largestErrorDirection_(0,0,1)
{
  ;
}

ProlateSpacepointMeasurement::ProlateSpacepointMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint)
  : SpacepointMeasurement(rawHitCoords, rawHitCov, detId, hitId, trackPoint), largestErrorDirection_(0,0,1)
{
  ;
}


SharedPlanePtr ProlateSpacepointMeasurement::constructPlane(const StateOnPlane& state) const {

  // copy state. Neglect covariance.
  StateOnPlane st(state);


  const TVector3 wire1(rawHitCoords_(0), rawHitCoords_(1), rawHitCoords_(2));

  const AbsTrackRep* rep = state.getRep();
  rep->extrapolateToLine(st, wire1, largestErrorDirection_);

  TVector3 dirInPoca = rep->getMom(st);
  dirInPoca.SetMag(1.);

  // check if direction is parallel to wire
  if (fabs(largestErrorDirection_.Angle(dirInPoca)) < 0.01){
    Exception exc("ProlateSpacepointMeasurement::constructPlane(): Cannot construct detector plane, track direction is parallel to largest error direction", __LINE__,__FILE__);
    throw exc;
  }

  // construct orthogonal vector
  TVector3 U = largestErrorDirection_.Cross(dirInPoca);

  return SharedPlanePtr(new DetPlane(wire1, U, largestErrorDirection_));
}


} /* End of namespace genfit */
