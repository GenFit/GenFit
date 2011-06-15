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
#include "GFSpacepointHitPolicy.h"

#include "assert.h"

#include "TMath.h"

#include "GFAbsRecoHit.h"

const std::string GFSpacepointHitPolicy::fPolicyName = "GFSpacepointHitPolicy";

TMatrixT<double> 
GFSpacepointHitPolicy::hitCoord(GFAbsRecoHit* hit,const GFDetPlane& plane)
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
GFSpacepointHitPolicy::hitCov(GFAbsRecoHit* hit,const GFDetPlane& plane)
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
GFSpacepointHitPolicy::detPlane(GFAbsRecoHit* hit, GFAbsTrackRep* rep)
{
  TMatrixT<double> rawcoord = hit->getRawHitCoord();
  TVector3 point(rawcoord[0][0],rawcoord[1][0],rawcoord[2][0]);

  TVector3 poca,dirInPoca;
  rep->extrapolateToPoint(point,poca,dirInPoca);

  fPlane.setO(point);
  fPlane.setNormal(dirInPoca);

  return fPlane;
}

ClassImp(GFSpacepointHitPolicy)
