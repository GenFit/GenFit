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

#include "AbsMeasurement.h"

#include <cassert>
#include <iostream>


namespace genfit {

AbsMeasurement::AbsMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint)
  : rawHitCoords_(rawHitCoords), rawHitCov_(rawHitCov), detId_(detId), hitId_(hitId), trackPoint_(trackPoint)
{
  assert(rawHitCov_.GetNrows() == rawHitCoords_.GetNrows());
}


AbsMeasurement::AbsMeasurement(const AbsMeasurement& o)
  : TObject(o),
    rawHitCoords_(o.rawHitCoords_),
    rawHitCov_(o.rawHitCov_),
    detId_(o.detId_),
    hitId_(o.hitId_),
    trackPoint_(o.trackPoint_)
{
  ;
}


AbsMeasurement::~AbsMeasurement()
{
  ;
}


AbsMeasurement& AbsMeasurement::operator=(const AbsMeasurement&) {
  fputs ("must not call AbsMeasurement::operator=\n",stderr);
  abort();
  return *this;
}


void AbsMeasurement::Print(const Option_t*) const {
  std::cout << "genfit::AbsMeasurement, detId = " << detId_ << ". hitId = " << hitId_ << "\n";
  std::cout << "Raw hit coordinates: "; rawHitCoords_.Print();
  std::cout << "Raw hit covariance: "; rawHitCov_.Print();
}


} /* End of namespace genfit */
