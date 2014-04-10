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

#include "StateOnPlane.h"
#include "AbsTrackRep.h"

#include <cassert>
#include <iostream>

namespace genfit {


void StateOnPlane::Print(Option_t*) const {
  std::cout << "genfit::StateOnPlane ";
  std::cout << " state vector: "; state_.Print();
  if (sharedPlane_ != NULL) {
    std::cout << " defined in plane "; sharedPlane_->Print();
    TVector3 pos, mom;
    getRep()->getPosMom(*this, pos, mom);
    std::cout << " 3D position: "; pos.Print();
    std::cout << " 3D momentum: "; mom.Print();
  }
}


// Modified from auto-generated Streamer to account for sharedPlane_
// also ignore rep_, this has to be set when loading the object owning
// this state.
void StateOnPlane::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::StateOnPlane.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::StateOnPlane thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      //TObject::Streamer(R__b);
      state_.Streamer(R__b);
      auxInfo_.Streamer(R__b);
      sharedPlane_.reset();  // needs to be set by owner;
      rep_ = NULL;  // needs to be set by owner
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
   } else {
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      //TObject::Streamer(R__b);
      state_.Streamer(R__b);
      auxInfo_.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}


} /* End of namespace genfit */
