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

#include "WireTrackCandHit.h"

#include <iostream>

namespace genfit {

WireTrackCandHit::WireTrackCandHit(int detId,
                               int hitId,
                               int planeId,
                               double sortingParameter,
                               char leftRight)
  : TrackCandHit(detId, hitId, planeId, sortingParameter),
    leftRight_(leftRight)
{
  ;
}


void WireTrackCandHit::Print(Option_t* option) const {
  std::cout << "  WireTrackCandHit. DetId = " << detId_
            << " \t HitId = " << hitId_
            << " \t PlaneId = " << planeId_
            << " \t SortingParameter = " << sortingParameter_
            << " \t leftRight = " << (int)leftRight_ << "\n";
}

} /* End of namespace genfit */
