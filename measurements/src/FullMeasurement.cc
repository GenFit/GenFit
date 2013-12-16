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

#include "FullMeasurement.h"

#include <Exception.h>
#include <RKTrackRep.h>
#include <HMatrixUnit.h>

#include <cassert>


namespace genfit {

FullMeasurement::FullMeasurement(int nDim)
  : AbsMeasurement(nDim), plane_()
{
  assert(nDim >= 1);
}


FullMeasurement::FullMeasurement(const MeasuredStateOnPlane& state, int detId, int hitId, TrackPoint* trackPoint)
  : AbsMeasurement(state.getState(), state.getCov(), detId, hitId, trackPoint), plane_(state.getPlane())
{
  assert(rawHitCoords_.GetNrows() == (int)state.getRep()->getDim());
}


SharedPlanePtr FullMeasurement::constructPlane(const StateOnPlane& state) const {
  if (!plane_) {
    Exception exc("FullMeasurement::constructPlane(): No plane has been set!", __LINE__,__FILE__);
    throw exc;
  }
  return plane_;
}


std::vector<MeasurementOnPlane*> FullMeasurement::constructMeasurementsOnPlane(const StateOnPlane& state) const {

  MeasurementOnPlane* mop = new MeasurementOnPlane(rawHitCoords_,
       rawHitCov_,
       state.getPlane(), state.getRep(), constructHMatrix(state.getRep()));

  std::vector<MeasurementOnPlane*> retVal;
  retVal.push_back(mop);
  return retVal;
}


const AbsHMatrix* FullMeasurement::constructHMatrix(const AbsTrackRep* rep) const {

  if (dynamic_cast<const RKTrackRep*>(rep) == NULL) {
    Exception exc("SpacepointMeasurement default implementation can only handle state vectors of type RKTrackRep!", __LINE__,__FILE__);
    throw exc;
  }

  return new HMatrixUnit();
}


void FullMeasurement::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::FullMeasurement.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::FullMeasurement thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsMeasurement baseClass0;
      baseClass0::Streamer(R__b);
      char flag;
      R__b >> flag;
      plane_.reset();
      if (flag) {
        plane_.reset(new DetPlane());
        plane_->Streamer(R__b);
      }
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      //This works around a msvc bug and should be harmless on other platforms
      typedef genfit::AbsMeasurement baseClass0;
      baseClass0::Streamer(R__b);
      if (plane_) {
        R__b << (char)1;
        plane_->Streamer(R__b);
      } else {
        R__b << (char)0;
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

} /* End of namespace genfit */
