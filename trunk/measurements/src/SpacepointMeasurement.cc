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

#include "SpacepointMeasurement.h"

#include "Exception.h"
#include "RKTrackRep.h"
#include "Tools.h"
#include "HMatrixUV.h"
#include "MeasurementOnPlane.h"

#include <cassert>


namespace genfit {

SpacepointMeasurement::SpacepointMeasurement(int nDim)
  : AbsMeasurement(nDim), weightedPlaneContruction_(true), G_(3), cutCov_(true)
{
  assert(nDim >= 3);

  //this->initG();
}

SpacepointMeasurement::SpacepointMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint,
    bool weightedPlaneContruction, bool cutCov)
  : AbsMeasurement(rawHitCoords, rawHitCov, detId, hitId, trackPoint),
    weightedPlaneContruction_(weightedPlaneContruction), cutCov_(cutCov)
{
  assert(rawHitCoords_.GetNrows() >= 3);

  if (weightedPlaneContruction_)
    this->initG();
}


SharedPlanePtr SpacepointMeasurement::constructPlane(const StateOnPlane& state) const {

  // copy state. Neglect covariance.
  StateOnPlane st(state);

  const TVector3 point(rawHitCoords_(0), rawHitCoords_(1), rawHitCoords_(2));

  if (weightedPlaneContruction_)
    st.extrapolateToPoint(point, G_);
  else
    st.extrapolateToPoint(point);

  return st.getPlane();
}


std::vector<MeasurementOnPlane*> SpacepointMeasurement::constructMeasurementsOnPlane(const StateOnPlane& state) const
{
  MeasurementOnPlane* mop = new MeasurementOnPlane(TVectorD(2),
       TMatrixDSym(3), // will be resized to 2x2 by Similarity later
       state.getPlane(), state.getRep(), constructHMatrix(state.getRep()));

  TVectorD& m = mop->getState();
  TMatrixDSym& V = mop->getCov();

  const TVector3& o(state.getPlane()->getO());
  const TVector3& u(state.getPlane()->getU());
  const TVector3& v(state.getPlane()->getV());

  // m
  m(0) = (rawHitCoords_(0)-o.X()) * u.X() +
         (rawHitCoords_(1)-o.Y()) * u.Y() +
         (rawHitCoords_(2)-o.Z()) * u.Z();

  m(1) = (rawHitCoords_(0)-o.X()) * v.X() +
         (rawHitCoords_(1)-o.Y()) * v.Y() +
         (rawHitCoords_(2)-o.Z()) * v.Z();

  //
  // V
  //
  V = rawHitCov_;

  // jac = dF_i/dx_j = s_unitvec * t_untivec, with s=u,v and t=x,y,z
  TMatrixD jac(3,2);
  jac(0,0) = u.X();
  jac(1,0) = u.Y();
  jac(2,0) = u.Z();
  jac(0,1) = v.X();
  jac(1,1) = v.Y();
  jac(2,1) = v.Z();

  if (cutCov_) { // cut (default)
    tools::invertMatrix(V);
    V.SimilarityT(jac);
    tools::invertMatrix(V);
  }
  else { // projection
    V.SimilarityT(jac);
  }

  std::vector<MeasurementOnPlane*> retVal;
  retVal.push_back(mop);
  return retVal;
}


const AbsHMatrix* SpacepointMeasurement::constructHMatrix(const AbsTrackRep* rep) const {
  if (dynamic_cast<const RKTrackRep*>(rep) == NULL) {
    Exception exc("SpacepointMeasurement default implementation can only handle state vectors of type RKTrackRep!", __LINE__,__FILE__);
    throw exc;
  }

  return new HMatrixUV();
}


void SpacepointMeasurement::initG() {
  rawHitCov_.GetSub(0, 2, G_);
  tools::invertMatrix(G_);
}


// Modified from auto-generated Streamer to account for non-persistent G_
void SpacepointMeasurement::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::SpacepointMeasurement.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(genfit::SpacepointMeasurement::Class(),this);

      if (weightedPlaneContruction_)
        this->initG();
   } else {
      R__b.WriteClassBuffer(genfit::SpacepointMeasurement::Class(),this);
   }
}


} /* End of namespace genfit */
