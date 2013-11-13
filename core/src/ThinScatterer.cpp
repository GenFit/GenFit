/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

#include "ThinScatterer.h"

#include <iostream>


namespace genfit {


void ThinScatterer::Print(const Option_t*) const {
  std::cout << "ThinScatterer, defined in plane: ";
  sharedPlane_->Print();
  std::cout << "Material properties: ";
  material_.Print();
}


void ThinScatterer::Streamer(TBuffer &R__b)
{
  // Stream an object of class genfit::ThinScatterer.

  // TODO: test

  //This works around a msvc bug and should be harmless on other platforms
  typedef ::genfit::ThinScatterer thisClass;
  UInt_t R__s, R__c;
  if (R__b.IsReading()) {
     Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
     //TObject::Streamer(R__b);
     sharedPlane_.reset(new DetPlane());
     sharedPlane_->Streamer(R__b);
     material_.Streamer(R__b);
     R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
  } else {
     R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
     //TObject::Streamer(R__b);
     sharedPlane_->Streamer(R__b);
     material_.Streamer(R__b);
     R__b.SetByteCount(R__c, kTRUE);
  }
}



} /* End of namespace genfit */
