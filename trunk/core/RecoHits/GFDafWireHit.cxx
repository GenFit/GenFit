/* Copyright 2011, Technische Universitaet Muenchen,
Authors: Karl Bicker

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


#include "GFDafWireHit.h"
#include <GFTools.h>
#include <GFException.h>
#include <cmath>


GFDafWireHit::GFDafWireHit(GFAbsWireHit* hit)
  : GFDafHit()
{
  // fill fRawHits with the hit
  fRawHits.push_back(hit);

  fHitUpd = false;

  // initialise the two weigths for left and right hit
	// set initial weights according to l/r resolution (default weights are left:1 right:1)
  fWeights.assign(2,1.);
	int lr = static_cast<GFAbsWireHit*>(fRawHits[0])->getLeftRightResolution();

	if (lr<0) { // left
	  fWeights[1] = 0.; // set right to 0
	}
	else if (lr>0){ // right
    fWeights[0] = 0.; // set left to 0
	}

	// set l/r resolution so that the plane is fixed
	static_cast<GFAbsWireHit*>(fRawHits[0])->setLeftRightResolution(1);
}


GFDafWireHit::~GFDafWireHit(){
  // set the l/r resolution of the wire hit according to what the DAF has determined
  if (fWeights[0] > fWeights[1]) static_cast<GFAbsWireHit*>(fRawHits[0])->setLeftRightResolution(-1);
  else static_cast<GFAbsWireHit*>(fRawHits[0])->setLeftRightResolution(1);
}


void GFDafWireHit::getMeasurement(const GFAbsTrackRep* rep,const GFDetPlane& pl,const TVectorD& statePred,const TMatrixDSym& covPred,TVectorD& m, TMatrixDSym& V) {

  if(fHitUpd && fDetPlane == pl) {
    m.ResizeTo(fHitCoord);
    V.ResizeTo(fHitCov);
    m = fHitCoord;
    V = fHitCov;
    return;
  }

  TMatrixDSym cov;
  TVectorD coord;

  fRawHits[0]->getMeasurement(rep, pl, statePred, covPred, coord, cov);

  // make sure fHitCoord and fHitCov have right dimensionality and set them to 0
  fHitCoord.ResizeTo(coord);
  fHitCov.ResizeTo(cov);

  GFTools::invertMatrix(cov);

  fHitCov = (fWeights[0] + fWeights[1]) * cov; // cov is already inverted

  // invert fHitCov
  GFTools::invertMatrix(fHitCov);

  //set the weighted-mean coord
  fHitCoord = fWeights[1] * cov * coord; // right side  // cov is already inverted
  coord(0) *= -1.; // invert the sign of the drift radius of the first (left) hit
  fHitCoord += fWeights[0] * cov * coord; // left side  // cov is already inverted


  fHitCoord *= fHitCov;

  //return by refernce
  m.ResizeTo(fHitCoord);
  V.ResizeTo(fHitCov);
  m = fHitCoord;
  V = fHitCov;
  fDetPlane = pl;
  fHitUpd = true;
}


void GFDafWireHit::getMeasurement(const GFAbsTrackRep* rep,const GFDetPlane& pl,const TVectorD& statePred,const TMatrixDSym& covPred,TVectorD& m, TMatrixDSym& V, unsigned int iHit){
  assert(iHit<2);

  fRawHits[0]->getMeasurement(rep, pl, statePred, covPred, m, V);

  if (iHit == 0) m[0] *= -1.; // invert the sign of the drift radius if left hit
}

ClassImp(GFDafWireHit)
