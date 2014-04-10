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

#include <iostream>
#include <TClass.h>

#include "MeasurementOnPlane.h"

namespace genfit {

MeasurementOnPlane::MeasurementOnPlane(const MeasurementOnPlane& other) :
    MeasuredStateOnPlane(other),
    weight_(other.weight_)
{
  hMatrix_.reset(other.hMatrix_->clone());
}


MeasurementOnPlane& MeasurementOnPlane::operator=(MeasurementOnPlane other) {
  swap(other);
  return *this;
}


void MeasurementOnPlane::swap(MeasurementOnPlane& other) {
  MeasuredStateOnPlane::swap(other);
  this->hMatrix_.swap(other.hMatrix_);
  std::swap(this->weight_, other.weight_);
}


void MeasurementOnPlane::Print(Option_t*) const
{
  std::cout << "genfit::MeasurementOnPlane, weight = " << weight_ << "\n";
  std::cout << " state vector: "; state_.Print();
  std::cout << " covariance matrix: "; cov_.Print();
  if (sharedPlane_ != NULL)
    std::cout << " defined in plane "; sharedPlane_->Print();
  std::cout << " hMatrix: "; hMatrix_->Print();

}


void MeasurementOnPlane::Streamer(TBuffer &R__b)
{
  // Stream an object of class genfit::MeasurementOnPlane.

  //This works around a msvc bug and should be harmless on other platforms
  typedef ::genfit::MeasurementOnPlane thisClass;
  UInt_t R__s, R__c;
  if (R__b.IsReading()) {
     Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
     MeasuredStateOnPlane::Streamer(R__b);
     hMatrix_.reset();
     char flag;
     R__b.ReadChar(flag);
     if (flag) {
       AbsHMatrix *h = 0;
       R__b >> h;
       hMatrix_.reset(h);
     }
     R__b >> weight_;
     R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
  } else {
     R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
     MeasuredStateOnPlane::Streamer(R__b);
     if (hMatrix_) {
       R__b.WriteChar(1);
       R__b << hMatrix_.get();
     } else {
       R__b.WriteChar(0);
     }
     R__b << weight_;
     R__b.SetByteCount(R__c, kTRUE);
  }
}

} /* End of namespace genfit */
