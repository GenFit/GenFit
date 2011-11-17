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

  _D[0][0] = (plane.getO())[0];
  _D[1][0] = (plane.getO())[1];
  _D[2][0] = (plane.getO())[2];

  _D *= -1.; 
  _D += hit->getRawHitCoord();
  //now the vector _D points from the origin of the plane to the hit point


  _U = plane.getU();
  _V = plane.getV();


  returnMat[0][0] = _D[0][0] * _U[0] + _D[1][0] * _U[1] + _D[2][0] * _U[2];
  returnMat[1][0] = _D[0][0] * _V[0] + _D[1][0] * _V[1] + _D[2][0] * _V[2];
  //std::cout << "hitCoord="<<std::endl;
  //returnMat.Print();
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
  jac[0][0] = _U[0];
  jac[1][0] = _U[1];
  jac[2][0] = _U[2];
  jac[0][1] = _V[0];
  jac[1][1] = _V[1];
  jac[2][1] = _V[2];

  TMatrixT<double> jac_orig = jac;
  TMatrixT<double> jac_t = jac.T();

  TMatrixT<double> result=jac_t * (rawCov * jac_orig);
  //std::cout << "hitCov="<<std::endl;
  //result.Print();
  return  result;
}

const GFDetPlane&
GFPseudoSpacepointWireHitPolicy::detPlane(GFAbsRecoHit* hit, GFAbsTrackRep* rep)
{

  TMatrixT<double> rawcoord = hit->getRawHitCoord();
  TVector3 point(rawcoord[0][0],rawcoord[1][0],rawcoord[2][0]);

  TVector3 wire1(point);
  TVector3 wire2(point);
  wire2 += fWireDirection;

  //  distance of one (the first) of the wire extremities from the plane
  Double_t d_from_refplane =  fDetPlane.dist(wire1).Mag();
  if(d_from_refplane < 1e-5) return fDetPlane;


  // point of closest approach
  TVector3 poca, poca_onwire, dirInPoca;

  rep->extrapolateToLine(wire1, wire2, poca, dirInPoca, poca_onwire);


  Double_t distance;
  distance = TMath::Sqrt(fabs(((wire1-poca).Mag2()*(wire2-wire1).Mag2()-pow((wire1-poca).Dot(wire2-wire1),2))/(wire2-wire1).Mag2()));

  // check poca inside tube
  if(distance > fMaxdistance) {
    GFException exc("distance poca-wire > maxdistance", __LINE__,__FILE__);
    throw exc;
  }

  // find plane
  // unitary vector along distance
  // poca (on track), poca_onwire (on wire)
  TVector3 fromwiretoextr = poca - poca_onwire;
  fromwiretoextr.SetMag(1.);
  // unitary vector along the wire
  TVector3 wiredirection = wire2 - wire1;
  wiredirection.SetMag(1.);

  // check orthogonality
  if(fabs(fromwiretoextr * wiredirection) > 1e-3) {
    GFException exc("fromwiretoextr*wiredirection > 1e-3", __LINE__,__FILE__);
    throw exc;
  }

  TVector3 U;
  U = fromwiretoextr;
  TVector3 V;
  V = wiredirection;
  U.SetMag(1.);
  V.SetMag(1.);

  TVector3 O = poca_onwire;

  fDetPlane = GFDetPlane(O, U, V);

  return fDetPlane;
}

ClassImp(GFPseudoSpacepointWireHitPolicy)
