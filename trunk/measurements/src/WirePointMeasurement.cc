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

#include "WirePointMeasurement.h"

#include <Exception.h>
#include <RKTrackRep.h>
#include <HMatrixUV.h>

#include <cassert>
#include <algorithm>


namespace genfit {

WirePointMeasurement::WirePointMeasurement(int nDim)
  : WireMeasurement(nDim)
{
  assert(nDim >= 8);
}

WirePointMeasurement::WirePointMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint)
  : WireMeasurement(rawHitCoords, rawHitCov, detId, hitId, trackPoint)
{
  assert(rawHitCoords_.GetNrows() >= 8);
}


SharedPlanePtr WirePointMeasurement::constructPlane(const StateOnPlane& state) const {

  // copy state. Neglect covariance.
  StateOnPlane st(state);

  TVector3 wire1(rawHitCoords_(0), rawHitCoords_(1), rawHitCoords_(2));
  TVector3 wire2(rawHitCoords_(3), rawHitCoords_(4), rawHitCoords_(5));

  //std::cout << " wire1(" << rawHitCoords_(0) << ", " << rawHitCoords_(1) << ", " << rawHitCoords_(2) << ")" << std::endl;
  //std::cout << " wire2(" << rawHitCoords_(3) << ", " << rawHitCoords_(4) << ", " << rawHitCoords_(5) << ")" << std::endl;

  // unit vector along the wire (V)
  TVector3 wireDirection = wire2 - wire1;
  wireDirection.SetMag(1.);

  //std::cout << " wireDirection(" << wireDirection.X() << ", " << wireDirection.Y() << ", " << wireDirection.Z() << ")" << std::endl;

  // point of closest approach
  const AbsTrackRep* rep = state.getRep();
  rep->extrapolateToLine(st, wire1, wireDirection);
  //const TVector3& poca = rep->getPos(st);
  TVector3 dirInPoca = rep->getMom(st);
  dirInPoca.SetMag(1.);
  //const TVector3& pocaOnWire = wire1 + wireDirection.Dot(poca - wire1)*wireDirection;

  // check if direction is parallel to wire
  if (fabs(wireDirection.Angle(dirInPoca)) < 0.01){
    Exception exc("WireMeasurement::detPlane(): Cannot construct detector plane, direction is parallel to wire", __LINE__,__FILE__);
    throw exc;
  }

  // construct orthogonal vector
  TVector3 U = dirInPoca.Cross(wireDirection);
  // U.SetMag(1.); automatically assured

  return SharedPlanePtr(new DetPlane(wire1, U, wireDirection));
}


std::vector<MeasurementOnPlane*> WirePointMeasurement::constructMeasurementsOnPlane(const StateOnPlane& state) const
{
  MeasurementOnPlane* mopR = new MeasurementOnPlane(TVectorD(2),
       TMatrixDSym(2),
       state.getPlane(), state.getRep(), constructHMatrix(state.getRep()));

  mopR->getState()(0) = rawHitCoords_(6);
  mopR->getState()(1) = rawHitCoords_(7);

  mopR->getCov()(0,0) = rawHitCov_(6,6);
  mopR->getCov()(1,0) = rawHitCov_(7,6);
  mopR->getCov()(0,1) = rawHitCov_(6,7);
  mopR->getCov()(1,1) = rawHitCov_(7,7);


  MeasurementOnPlane* mopL = new MeasurementOnPlane(*mopR);
  mopL->getState()(0) *= -1;

  // set left/right weights
  if (leftRight_ < 0) {
    mopL->setWeight(1);
    mopR->setWeight(0);
  }
  else if (leftRight_ > 0) {
    mopL->setWeight(0);
    mopR->setWeight(1);
  }
  else {
    double val = 0.5 * pow(std::max(0., 1 - rawHitCoords_(6)/maxDistance_), 2.);
    mopL->setWeight(val);
    mopR->setWeight(val);
  }

  std::vector<MeasurementOnPlane*> retVal;
  retVal.push_back(mopL);
  retVal.push_back(mopR);
  return retVal;
}

const AbsHMatrix* WirePointMeasurement::constructHMatrix(const AbsTrackRep* rep) const {
  if (dynamic_cast<const RKTrackRep*>(rep) == NULL) {
    Exception exc("WirePointMeasurement default implementation can only handle state vectors of type RKTrackRep!", __LINE__,__FILE__);
    throw exc;
  }

  return new HMatrixUV();
}

} /* End of namespace genfit */
