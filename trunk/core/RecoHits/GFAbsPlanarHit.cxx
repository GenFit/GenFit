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

#include "GFAbsPlanarHit.h"
#include <GFException.h>

void
GFAbsPlanarHit::getMeasurement(const GFAbsTrackRep* rep,
                               const GFDetPlane& pl,
                               const TVectorD& statePred,
                               const TMatrixDSym& covPred,
                               TVectorD& m,
                               TMatrixDSym& V) {

  if (pl != fDetPlane){
    GFException exc("GFAbsPlanarHit::getMeasurement(): Trying to get measurement in a plane that does not match the physical detector plane!", __LINE__,__FILE__);
    throw exc;
  }

  static_cast<void>(rep);
  static_cast<void>(statePred);
  static_cast<void>(covPred);
  m.ResizeTo(fHitCoord);
  V.ResizeTo(fHitCov);
  m = fHitCoord;
  V = fHitCov;
}


ClassImp(GFAbsPlanarHit)
