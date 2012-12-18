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

#include "GFAbsWirepointHit.h"


void
GFAbsWirepointHit::getMeasurement(const GFAbsTrackRep* rep,
                                  const GFDetPlane& pl,
                                  const TVectorD& statePred,
                                  const TMatrixDSym& covPred,
                                  TVectorD& m,
                                  TMatrixDSym& V) {

  static_cast<void>(rep);
  static_cast<void>(statePred);
  static_cast<void>(covPred);

  checkPlane(pl);

  // m
  m.ResizeTo(2);
  m(0) = fHitCoord(6);
  m(1) = fHitCoord(7);

  
  // V
  V.ResizeTo(2,2);
  V(0,0) = fHitCov(6,6);
  V(1,0) = fHitCov(7,6);
  V(0,1) = fHitCov(6,7);
  V(1,1) = fHitCov(7,7);

}


ClassImp(GFAbsWirepointHit)
