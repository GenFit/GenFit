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
#include "GFPseudoSpacepointWireHitPolicy.h"

#include "assert.h"

#include "TMath.h"

#include "GFAbsRecoHit.h"
#include "GFException.h"

const std::string GFPseudoSpacepointWireHitPolicy::fPolicyName = "GFPseudoSpacepointWireHitPolicy";


GFPseudoSpacepointWireHitPolicy::GFPseudoSpacepointWireHitPolicy() :
    fMaxdistance(1.E50),
    fWireDirection(0, 0, 1)
{
  ;
}


TMatrixT<double> 
GFPseudoSpacepointWireHitPolicy::hitCoord(GFAbsRecoHit* hit,const GFDetPlane& plane)
{
  TMatrixT<double> returnMat(2,1);

  TMatrixT<double> _D(3,1);
  TVector3 _U;
  TVector3 _V;

  _D(0,0) = plane.getO().X();
  _D(1,0) = plane.getO().Y();
  _D(2,0) = plane.getO().Z();

  _D *= -1.; 
  _D += hit->getRawHitCoord();
  //now the vector _D points from the origin of the plane to the hit point

  _U = plane.getU();
  _V = plane.getV();

  returnMat(0,0) = _D(0,0) * _U.X() + _D(1,0) * _U.Y() + _D(2,0) * _U.Z();
  returnMat(1,0) = _D(0,0) * _V.X() + _D(1,0) * _V.Y() + _D(2,0) * _V.Z();

  return returnMat;
}

TMatrixT<double> 
GFPseudoSpacepointWireHitPolicy::hitCov(GFAbsRecoHit* hit,const GFDetPlane& plane)
{
  TVector3 _U;
  TVector3 _V;

  _U = plane.getU();
  _V = plane.getV();

  TMatrixT<double> rawCov = hit->getRawHitCov();

  TMatrixT<double> jac(3,2);
  
  // jac = dF_i/dx_j = s_unitvec * t_untivec, with s=u,v and t=x,y,z
  jac(0,0) = _U.X();
  jac(1,0) = _U.Y();
  jac(2,0) = _U.Z();
  jac(0,1) = _V.X();
  jac(1,1) = _V.Y();
  jac(2,1) = _V.Z();

  TMatrixT<double> jac_orig = jac;
  TMatrixT<double> jac_t = jac.T();

  TMatrixT<double> result=jac_t * (rawCov * jac_orig);

  return  result;
}

const GFDetPlane&
GFPseudoSpacepointWireHitPolicy::detPlane(GFAbsRecoHit* hit, GFAbsTrackRep* rep)
{
  //  distance of one (the first) of the wire extremities from the plane
  //Double_t d_from_refplane =  fDetPlane.dist(wire1).Mag();
  //if(d_from_refplane < 1e-5) return fDetPlane;

  TMatrixT<double> rawcoord = hit->getRawHitCoord();
  TVector3 point(rawcoord(0,0), rawcoord(1,0), rawcoord(2,0));

  TVector3 wire1(point);
  TVector3 wire2(point);
  wire2 += fWireDirection;

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

  fDetPlane = GFDetPlane(poca_onwire, U, wiredirection);

  return fDetPlane;
}

ClassImp(GFPseudoSpacepointWireHitPolicy)
