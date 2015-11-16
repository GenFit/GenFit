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

#include "AbsFitterInfo.h"

#include <TClass.h>
#include <TBuffer.h>

namespace genfit {

AbsFitterInfo::AbsFitterInfo() :
  trackPoint_(NULL),
  rep_(NULL)
{
  ;
}

AbsFitterInfo::AbsFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep) :
  trackPoint_(trackPoint),
  rep_(rep)
{
  ;
}

void AbsFitterInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::AbsFitterInfo.
   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::AbsFitterInfo thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      //TObject::Streamer(R__b);
      // See the long comment in AbsFinitePlane::Streamer.  This
      // creates a duplicate of the DetPlane.
      TClass* cl = TClass::Load(R__b);
      DetPlane *p = (DetPlane*)(cl->New());
      cl->Streamer(p, R__b);
      sharedPlane_.reset(p);
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      //TObject::Streamer(R__b);
      // See the long comment in AbsFinitePlane::Streamer.  This
      // stores a duplicate of the DetPlane, allowing for typesafe
      // reading.
      sharedPlane_->IsA()->Store(R__b);
      sharedPlane_->Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

} /* End of namespace genfit */
