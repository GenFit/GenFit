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

#include "PlanarMeasurement.h"

#include <Exception.h>
#include <RKTrackRep.h>
#include <HMatrixU.h>
#include <HMatrixV.h>
#include <HMatrixUV.h>

#include <cassert>


namespace genfit {

PlanarMeasurement::PlanarMeasurement(int nDim)
  : AbsMeasurement(nDim), physicalPlane_(), planeId_(-1), stripV_(false)
{
  assert(nDim >= 1);
}

PlanarMeasurement::PlanarMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint)
  : AbsMeasurement(rawHitCoords, rawHitCov, detId, hitId, trackPoint), physicalPlane_(), planeId_(-1), stripV_(false)
{
  assert(rawHitCoords_.GetNrows() >= 1);
}


SharedPlanePtr PlanarMeasurement::constructPlane(const StateOnPlane&) const {
  if (!physicalPlane_) {
    Exception exc("PlanarMeasurement::constructPlane(): No plane has been set!", __LINE__,__FILE__);
    throw exc;
  }
  return physicalPlane_;
}


std::vector<MeasurementOnPlane*> PlanarMeasurement::constructMeasurementsOnPlane(const StateOnPlane& state) const {

  MeasurementOnPlane* mop = new MeasurementOnPlane(rawHitCoords_,
       rawHitCov_,
       state.getPlane(), state.getRep(), constructHMatrix(state.getRep()));

  std::vector<MeasurementOnPlane*> retVal;
  retVal.push_back(mop);
  return retVal;
}


const AbsHMatrix* PlanarMeasurement::constructHMatrix(const AbsTrackRep* rep) const {

  if (dynamic_cast<const RKTrackRep*>(rep) == NULL) {
    Exception exc("SpacepointMeasurement default implementation can only handle state vectors of type RKTrackRep!", __LINE__,__FILE__);
    throw exc;
  }

  switch(rawHitCoords_.GetNrows()) {
  case 1:
    if (stripV_)
      return new HMatrixV();
    return new HMatrixU();

  case 2:
    return new HMatrixUV();

  default:
    Exception exc("PlanarMeasurement default implementation can only handle 1D (strip) or 2D (pixel) measurements!", __LINE__,__FILE__);
    throw exc;
  }

}

void PlanarMeasurement::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::PlanarMeasurement.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::PlanarMeasurement thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsMeasurement baseClass0;
      baseClass0::Streamer(R__b);
      char flag;
      R__b >> flag;
      physicalPlane_.reset();
      if (flag) {
        physicalPlane_.reset(new DetPlane());
        physicalPlane_->Streamer(R__b);
      }
      R__b >> planeId_;
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsMeasurement baseClass0;
      baseClass0::Streamer(R__b);
      if (physicalPlane_) {
        R__b << (char)1;
        physicalPlane_->Streamer(R__b);
      } else {
        R__b << (char)0;
      }
      R__b << planeId_;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

} /* End of namespace genfit */
