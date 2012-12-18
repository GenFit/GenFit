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

#include <math.h>

#include "GFAbsProlateSpacepointHit.h"
#include "GFAbsTrackRep.h"
#include <GFException.h>


GFAbsProlateSpacepointHit::GFAbsProlateSpacepointHit(unsigned int dim) :
  GFAbsSpacepointHit(dim),
  fLargestErrorDirection(0, 0, 1)
{
  assert(dim >= 3);;
}


const GFDetPlane&
GFAbsProlateSpacepointHit::getDetPlane(GFAbsTrackRep* rep) {

  TVector3 wire1(fHitCoord(0), fHitCoord(1), fHitCoord(2));
  TVector3 wire2(wire1);
  wire2 += fLargestErrorDirection;

  // point of closest approach
  TVector3 poca, poca_onwire, dirInPoca;
  rep->extrapolateToLine(wire1, wire2, poca, dirInPoca, poca_onwire);

  // check if direction is parallel to wire
  if (fabs(fLargestErrorDirection.Angle(dirInPoca)) < 0.01){
    GFException exc("GFAbsProlateSpacepointHit::getDetPlane(): Cannot construct detector plane, track direction is parallel to largest error direction", __LINE__,__FILE__);
    throw exc;
  }

  // construct orthogonal vector
  TVector3 U = fLargestErrorDirection.Cross(dirInPoca);

  fDetPlane.set(poca_onwire, U, fLargestErrorDirection);

  return fDetPlane;
}

ClassImp(GFAbsProlateSpacepointHit)
