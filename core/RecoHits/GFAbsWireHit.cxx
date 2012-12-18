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
// ----------------------------------------------------------
// Please see GFAbsWireHit.h  before using this class.
// ----------------------------------------------------------

#include <math.h>

#include "GFAbsWireHit.h"

#include <assert.h>
#include "GFAbsTrackRep.h"
#include <GFException.h>


GFAbsWireHit::GFAbsWireHit(unsigned int dim) :
  GFAbsRecoHit(dim),
  fMaxdistance(1.E50),
  fLeftRight(0)
{
  assert(dim >= 7);
}


void
GFAbsWireHit::getMeasurement(const GFAbsTrackRep* rep,
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
  m.ResizeTo(1);
  m(0) = fHitCoord(6);

  // V
  V.ResizeTo(1,1);
  V(0,0) = fHitCov(6,6);
}


const GFDetPlane&
GFAbsWireHit::getDetPlane(GFAbsTrackRep* rep)
{
  assert(fHitCoord.GetNrows()>6);
  TVector3 wire1(fHitCoord(0), fHitCoord(1), fHitCoord(2));
  TVector3 wire2(fHitCoord(3), fHitCoord(4), fHitCoord(5));

  //  distance of one (the first) of the wire extremities from the plane
  //Double_t d_from_refDetPlane =  fDetPlane.dist(wire1).Mag();
  //if(d_from_refDetPlane < 1e-5) return fDetPlane;
  
  // point of closest approach
  TVector3 poca, poca_onwire, dirInPoca;
  rep->extrapolateToLine(wire1, wire2, poca, dirInPoca, poca_onwire);
  
  // check distance of poca to wire
  if((poca - poca_onwire).Mag() > fMaxdistance) {
    GFException exc("GFAbsWireHit::detPlane(): distance poca-wire > maxdistance", __LINE__,__FILE__);
    throw exc;    
  }

  // unitary vector along the wire (V)
  TVector3 wiredirection = wire2 - wire1; 
  wiredirection.SetMag(1.);
  
  // check if direction is parallel to wire
  if (fabs(wiredirection.Angle(dirInPoca)) < 0.01){
    GFException exc("GFAbsWireHit::detPlane(): Cannot construct detector plane, direction is parallel to wire", __LINE__,__FILE__);
    throw exc;
  }
  
  // construct orthogonal vector
  TVector3 U = dirInPoca.Cross(wiredirection);

  // check left/right ambiguity
  if (fLeftRight == 0){ // auto select
    if ((poca-poca_onwire)*U < 0) U *= -1.;
  }
  else if (fLeftRight < 0) U *= -1.;

  fDetPlane.set(wire1, U, wiredirection);
  
  return fDetPlane;
}


void
GFAbsWireHit::setLeftRightResolution(int lr){
  if (lr==0) fLeftRight = 0;
  else if (lr<0) fLeftRight = -1;
  else fLeftRight = 1;
}


void
GFAbsWireHit::checkPlane(const GFDetPlane& plane)
{
  assert(fHitCoord.GetNrows()>6);

  TVector3 wiredirection(fHitCoord(3) - fHitCoord(0),
                         fHitCoord(4) - fHitCoord(1),
                         fHitCoord(5) - fHitCoord(2));

  wiredirection.SetMag(1.);

  if(fabs(fabs(wiredirection.Dot(plane.getV())) - 1) > 1e-3) {
    GFException exc("GFAbsWireHit: plane not valid!!", __LINE__,__FILE__);
    throw exc;
  }
}


ClassImp(GFAbsWireHit)
