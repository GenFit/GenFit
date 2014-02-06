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

#include "TrackCandHit.h"

#include <iostream>

namespace genfit {

TrackCandHit::TrackCandHit(int detId,
                               int hitId,
                               int planeId,
                               double sortingParameter)
  : detId_(detId),
    hitId_(hitId),
    planeId_(planeId),
    sortingParameter_(sortingParameter)
{
  ;
}


void TrackCandHit::Print(Option_t*) const {
  std::cout << "  TrackCandHit. DetId = " << detId_
            << " \t HitId = " << hitId_
            << " \t PlaneId = " << planeId_
            << " \t SortingParameter = " << sortingParameter_ << "\n";
}


bool operator== (const TrackCandHit& lhs, const TrackCandHit& rhs){
  if(lhs.detId_ == rhs.detId_ &&
     lhs.hitId_ == rhs.hitId_ &&
     lhs.planeId_ == rhs.planeId_)
    return true;
  return false;
}

} /* End of namespace genfit */
