/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

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

#include "GFTrackCandHit.h"

#include <iostream>


ClassImp(GFTrackCandHit)


GFTrackCandHit::GFTrackCandHit(int detId,
                               int hitId,
                               int planeId,
                               double rho)
  : fDetId(detId),
    fHitId(hitId),
    fPlaneId(planeId),
    fRho(rho)
{
  ;
}


GFTrackCandHit::~GFTrackCandHit(){
  ;
}


void GFTrackCandHit::Print(Option_t* option) const {
  std::cout << "  GFTrackCandHit. DetId = " << fDetId
            << " \t HitId = " << fHitId
            << " \t PlaneId = " << fPlaneId
            << " \t Rho = " << fRho << "\n";
}


bool operator== (const GFTrackCandHit& lhs, const GFTrackCandHit& rhs){
  if(lhs.fDetId == rhs.fDetId &&
     lhs.fHitId == rhs.fHitId &&
     lhs.fPlaneId == rhs.fPlaneId)
    return true;
  return false;
}

