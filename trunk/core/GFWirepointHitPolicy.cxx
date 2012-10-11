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
// Please see GFWirepointHitPolicy.h  before using this class.
// ----------------------------------------------------------

#include "GFWirepointHitPolicy.h"

#include "assert.h"
#include <cmath>

#include "TMath.h"
#include "TVector3.h"

#include "GFAbsRecoHit.h"
#include "GFException.h"

const std::string GFWirepointHitPolicy::fPolicyName = "GFWirepointHitPolicy";

GFWirepointHitPolicy::GFWirepointHitPolicy() :
    fMaxdistance(1.E50),
    fLeftRight(0)
{
  ;
}

TMatrixT<double> 
GFWirepointHitPolicy::hitCoord(GFAbsRecoHit* hit,const GFDetPlane& plane)
{
  TMatrixT<double> returnMat(2,1);

  checkPlane(hit, plane);

  // raw x1, y1, z1, x2, y2, z2, rdrift, zreco
  TMatrixT<double> rC = hit->getRawHitCoord();

  returnMat(0,0) = rC(6,0);
  returnMat(1,0) = rC(7,0);
  return returnMat;
}

TMatrixT<double> 
GFWirepointHitPolicy::hitCov(GFAbsRecoHit* hit,const GFDetPlane& plane)
{
  checkPlane(hit, plane);

  TMatrixT<double> returnCov(2,2);
  TMatrixT<double> rawCov = hit->getRawHitCov();

  returnCov(0,0) = rawCov(6,6);
  returnCov(1,0) = rawCov(7,6);
  returnCov(0,1) = rawCov(6,7);
  returnCov(1,1) = rawCov(7,7);

  return  returnCov;
}



void GFWirepointHitPolicy::checkPlane(GFAbsRecoHit* hit,const GFDetPlane& plane)
{
  // raw x1, y1, z1, x2, y2, z2, rdrift, zreco
  TMatrixT<double> rC = hit->getRawHitCoord();

  assert(rC.GetNrows()==8);
  
  TVector3 wire1(rC(0,0), rC(1,0), rC(2,0));
  TVector3 wire2(rC(3,0), rC(4,0), rC(5,0));
  TVector3 wiredirection = wire1 - wire2;
  
  TVector3 vaxis = plane.getV();
  wiredirection.SetMag(1.);
  vaxis.SetMag(1.);

  if(fabs(TMath::Abs(wiredirection.Dot(vaxis)) - 1) > 1e-3) {
    std::cout << "GFWireHitPolicy: plane not valid!!" << std::endl;
  }
}


const GFDetPlane& 
GFWirepointHitPolicy::detPlane(GFAbsRecoHit* hit, GFAbsTrackRep* rep)
{
   //  distance of one (the first) of the wire extremities from the plane
   //Double_t d_from_refplane =  fDetPlane.dist(wire1).Mag();
   //if(d_from_refplane < 1e-5) return fDetPlane;

   TMatrixT<double> x = hit->getRawHitCoord();
   assert(x.GetNrows()==8);
   TVector3 wire1(x(0,0), x(1,0), x(2,0));
   TVector3 wire2(x(3,0), x(4,0), x(5,0));

   // point of closest approach
   TVector3 poca, poca_onwire, dirInPoca;
   rep->extrapolateToLine(wire1, wire2, poca, dirInPoca, poca_onwire);

   // check distance of poca to wire
   if((poca - poca_onwire).Mag() > fMaxdistance) {
     GFException exc("GFWireHitPolicy::detPlane(): distance poca-wire > maxdistance", __LINE__,__FILE__);
     throw exc;
   }

   // unitary vector along the wire (V)
   TVector3 wiredirection = wire2 - wire1;
   wiredirection.SetMag(1.);

   // check if direction is parallel to wire
   if (fabs(wiredirection.Angle(dirInPoca)) < 0.01){
     GFException exc("GFWireHitPolicy::detPlane(): Cannot construct detector plane, direction is parallel to wire", __LINE__,__FILE__);
     throw exc;
   }

   // construct orthogonal vector
   TVector3 U = wiredirection.Cross(dirInPoca);
   U.SetMag(1.);

   // check left/right ambiguity
   if (fLeftRight == 0){ // auto select
     if ((poca-poca_onwire)*U < 0) U *= -1.;
   }
   else if (fLeftRight < 0) U *= -1.;

   fDetPlane = GFDetPlane(poca_onwire, U, wiredirection);

   return fDetPlane;
}

void
GFWirepointHitPolicy::setLeftRightResolution(int lr){
  if (lr==0) fLeftRight = 0;
  else if (lr<0) fLeftRight = -1;
  else fLeftRight = 1;
}

ClassImp(GFWirepointHitPolicy)
