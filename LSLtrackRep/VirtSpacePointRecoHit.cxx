/* Copyright 2008-2009, Technische Universitaet Muenchen,
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
#include "VirtSpacePointRecoHit.h"

#include "LSLTrackRep.h"

ClassImp(VirtSpacePointRecoHit)

VirtSpacePointRecoHit::~VirtSpacePointRecoHit()
{}

VirtSpacePointRecoHit::VirtSpacePointRecoHit()
: SpacepointRecoHit(NparHitRep)
{}

VirtSpacePointRecoHit::VirtSpacePointRecoHit(double x, double y, double z)
  : SpacepointRecoHit(NparHitRep)
{
  fHitCoord[0][0] = x;
  fHitCoord[1][0] = y;
  fHitCoord[2][0] = z;
  fHitCov.UnitMatrix();
}

VirtSpacePointRecoHit::VirtSpacePointRecoHit(const TVector3& pos)
  : SpacepointRecoHit(NparHitRep)
{
  fHitCoord[0][0] = pos.X();
  fHitCoord[1][0] = pos.Y();
  fHitCoord[2][0] = pos.Z();
  fHitCov.UnitMatrix();
}

GFAbsRecoHit*
VirtSpacePointRecoHit::clone()
{ 
  return new VirtSpacePointRecoHit(*this);
}


TMatrixT<double>
VirtSpacePointRecoHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if (dynamic_cast<const LSLTrackRep*>(stateVector) != NULL) {
    TMatrixT<double> HMatrix (NparHitRep,5);

    HMatrix[0][0] = 1.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 0.;
    HMatrix[0][4] = 0.;

    HMatrix[1][0] = 0.;
    HMatrix[1][1] = 1.;
    HMatrix[1][2] = 0.;
    HMatrix[1][3] = 0.;
    HMatrix[1][4] = 0.;

    HMatrix[2][0] = 0.;
    HMatrix[2][1] = 0.;
    HMatrix[2][2] = 0.;
    HMatrix[2][3] = 0.;
    HMatrix[2][4] = 0.;
    return HMatrix;
  }
  else {
    std::cerr << "VirtSpacePointRecoHit can only handle state"
              << " vectors of type LSLTrackRep -> abort" << std::endl;
    throw;
  }

}


