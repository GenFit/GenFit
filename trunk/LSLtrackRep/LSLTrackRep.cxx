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
#include "LSLTrackRep.h"

#include <iostream>

#include "TMath.h"
#include "TMatrixD.h"

#include "Nystrom.h"
#include "RecoHits/GFAbsRecoHit.h"
#include "GFAbsBField.h"
#include "LSLEQM.h"
#include "AbsNystromEQM.h"

LSLTrackRep::LSLTrackRep()
  : GFAbsTrackRep(5), s(0) ,_acc(1E-2), _adaptive(false)
{
  
  _eqm=new LSLEQM(0); // this will leak memory in ROOT!
  fRefPlane=GFDetPlane(TVector3(0,0,0),TVector3(1,0,0),TVector3(0,1,0));
  
}

LSLTrackRep::LSLTrackRep(double z, double x, double y, 
			 double dxdz, double dydz, double invp,
			 double sigx, double sigy, 
			 double sigdxdz, double sigdydz, 
			 double siginvp,
			 GFAbsBField* field) 
  : GFAbsTrackRep(5), s(0), _acc(1E-2), _adaptive(false)
{
  s=z;
  fState[0]=x;
  fState[1]=y;
  fState[2]=dxdz;
  fState[3]=dydz;
  fState[4]=invp;

  fCov[0][0]=sigx;
  fCov[1][1]=sigy;
  fCov[2][2]=sigdxdz;
  fCov[3][3]=sigdydz;
  fCov[4][4]=siginvp;
  
  _eqm=new LSLEQM(field);

  fRefPlane=GFDetPlane(TVector3(x,y,z),TVector3(1,0,0),TVector3(0,1,0));
}

LSLTrackRep::LSLTrackRep(const LSLTrackRep& rep) 
  : GFAbsTrackRep(rep)
{
  _eqm=new LSLEQM(0);
  _acc=rep._acc;
  _adaptive=rep._adaptive;
  s=rep.s;
}


LSLTrackRep::~LSLTrackRep()
{
  if(_eqm!=NULL)delete _eqm;
}

void
LSLTrackRep::init(const TVector3& pos,
		  double dxdz, double dydz, double invp,
		  double sigx, double sigy, 
		  double sigdxdz, double sigdydz, 
		  double siginvp, 
		  GFAbsBField* field) 
{
  s=pos.Z();
  fState[0]=pos.X();
  fState[1]=pos.Y();
  fState[2]=dxdz;
  fState[3]=dydz;
  fState[4]=invp;

  fCov[0][0]=sigx;
  fCov[1][1]=sigy;
  fCov[2][2]=sigdxdz;
  fCov[3][3]=sigdydz;
  fCov[4][4]=siginvp;

  if(_eqm!=NULL)delete _eqm;
  _eqm=new LSLEQM(field);

  fRefPlane=GFDetPlane(pos,TVector3(1,0,0),TVector3(0,1,0));
}


void 
LSLTrackRep::SetBField(GFAbsBField* b)
{
  if(_eqm!=NULL)delete _eqm;
  _eqm=new LSLEQM(b);
}

double
LSLTrackRep::extrapolate(const GFDetPlane& pl, 
			 TVectorD& statePred)
{
  statePred.ResizeTo(fDimension,1);
  double sExtrapolateTo=pl.getO().Z();
  if(sExtrapolateTo<-1000 || sExtrapolateTo>5000)return 0;
  Nystrom rungeKutta(_eqm);
  rungeKutta.setAccuracy(_acc);
  rungeKutta.setAdaptive(_adaptive);
  //prepare the vectors

  //std::cout<<"s before extrapolation: "<<s<<std::endl;
  //std::cout<<"s_to: "<<sExtrapolateTo<<std::endl;

  TVectorT<double> u(3);u[0]=fState[0];u[1]=fState[1];u[2]=s;
  TVectorT<double> uprim(3);uprim[0]=fState[2];uprim[1]=fState[3];uprim[2]=1.;
  TVectorT<double> par(1);par[0]=fState[4];

  TVectorT<double> unew(3);unew=u;
  TVectorT<double> uprimnew(3);uprimnew=uprim;
  double l=rungeKutta.propagate(s,sExtrapolateTo,
		       u, uprim, par,
		       unew, uprimnew);
  // write results into statePred
  statePred[0]=unew[0];
  statePred[1]=unew[1];
  statePred[2]=uprimnew[0];
  statePred[3]=uprimnew[1];
  statePred[4]=fState[4];
  //std::cout<<"unew[2]=z="<<unew[2];
  //std::cout<<"s after extrapolation: "<<s<<std::endl;
  return l;
}

/*
void LSLTrackRep::extrapolate(const GFDetPlane& pl,
			      const TVectorD& stateFrom,
			      TVectorD& stateResult) {
  stateResult.ResizeTo(5,1);
  double s=
  Nystrom rungeKutta(_eqm);
  //prepare the vectors
  TVectorT<double> u(3);
  u[0]=stateFrom[0][0];
  u[1]=stateFrom[1][0];
  u[2]=sExtrapolateFrom;
  TVectorT<double> uprim(3);uprim[0]=stateFrom[2][0];uprim[1]=stateFrom[3][0];uprim[2]=1;
  TVectorT<double> par(1);par[0]=stateFrom[4][0];

  TVectorT<double> unew(3);
  TVectorT<double> uprimnew(3);
  rungeKutta.propagate(sExtrapolateFrom,sExtrapolateTo,
					   u, uprim, par,
					   unew, uprimnew);
  // write results into statePred
  stateResult[0][0]=unew[0];
  stateResult[1][0]=unew[1];
  stateResult[2][0]=uprimnew[0];
  stateResult[3][0]=uprimnew[1];
  stateResult[4][0]=state[4][0];
  //std::cout<<"unew[2]=z="<<unew[2];

}
*/

double
LSLTrackRep::extrapolate(const GFDetPlane& pl, 
			 TVectorD& statePred,
			 TMatrixDSym& covPred)
{
  statePred.ResizeTo(fDimension,1);
  covPred.ResizeTo(fDimension,fDimension);
  TMatrixD jacobian;
  //std::cout << "Extr from To: " << s << " " << sExtrapolateTo << std::endl;
  double l=extrapolate(pl,statePred);
  // covPred=JCovJ^T with J being Jacobian
  jacobian.ResizeTo(5,5);
  Jacobian(pl,statePred,jacobian);
  covPred=fCov;
  covPred.Similarity(jacobian);
  return l;
}

double
LSLTrackRep::extrapolateToPoint(const TVector3& p,
			       TVector3& poca,
			       TVector3& normVec){
  GFDetPlane plane;
  int dim = getDim();
  TVectorD statePred(dim,1);
  TMatrixDSym covPred(dim);
  plane.setO(p);
  plane.setU(TVector3(1,0,0));
  plane.setV(TVector3(0,1,0));
  double length = extrapolate(plane,statePred,covPred);
  poca.SetXYZ(statePred[0],statePred[1],plane.getO().Z());
  normVec.SetXYZ(statePred[2],statePred[3],plane.getO().Z());
  return length;
}


void 
LSLTrackRep::stepalong(double h){
  // create new detplane:
  GFDetPlane newp(fRefPlane);
  newp.setO(newp.getO()+TVector3(0,0,h));
  GFAbsTrackRep::extrapolate(newp);
}

void
LSLTrackRep::Jacobian(const GFDetPlane& pl,
		      const TVectorD& statePred,
		      TMatrixD& jacResult){
  
  TVectorD difPred(statePred);
  // do column wise differntiation:
  for(int icol=0;icol<5;++icol){
    // choose step
    double h=TMath::Abs(fState[icol])*1.e-4;
    if(h<1e-13)h=1.e-13;
    // vary the state
    fState[icol]+=h;
    extrapolate(pl,difPred);
    // difference:
    difPred-=statePred;
    // remove variation from state
    fState[icol]-=h;
    // fill jacobian with difference quotient
    for(int irow=0;irow<5;++irow)jacResult[irow][icol]=difPred[irow]/h;
  }
}
		

TVector3 
LSLTrackRep::getPos(const GFDetPlane& pl)
{
  double z=pl.getO().Z();
  TVectorD statePred(fState);
  GFDetPlane p(TVector3(0,0,z),TVector3(1,0,0),TVector3(0,1,0));
  extrapolate(p,statePred);
  return TVector3(statePred[0],statePred[1],z);
}
 
TVector3 
LSLTrackRep::getMom(const GFDetPlane& pl)
{
  double z=pl.getO().Z();
  TVectorD statePred(fState);
  //statePred.Print();
  GFDetPlane p(TVector3(0,0,z),TVector3(1,0,0),TVector3(0,1,0));
  extrapolate(p,statePred);
  //statePred.Print();
  TVector3 result(statePred[2],statePred[3],1);
  if(TMath::Abs(statePred[4])!=0){
    result.SetMag(1./TMath::Abs(statePred[4]));
  }
  else result.SetMag(100);
  return result;
}

void
LSLTrackRep::getPosMom(const GFDetPlane& pl,TVector3& pos,TVector3& mom)
{
  double z=pl.getO().Z();
  TVectorD statePred(fState);
  GFDetPlane p(TVector3(0,0,z),TVector3(1,0,0),TVector3(0,1,0));
  extrapolate(p,statePred);
  pos.SetXYZ(statePred[0],statePred[1],z);
  mom.SetXYZ(statePred[2],statePred[3],1);
  if(TMath::Abs(statePred[4])!=0){
    mom.SetMag(1./TMath::Abs(statePred[4]));
  }
  else mom.SetMag(100);
}


TVectorT<double> 
LSLTrackRep::getGlobal() {// (x,y,z,px,py,pz)
  TVector3 pos=GFAbsTrackRep::getPos();
  TVector3 mom=GFAbsTrackRep::getMom();
  double par[6];
  par[0]=pos.X();par[1]=pos.Y();par[2]=pos.Z();
  par[3]=mom.X();par[4]=mom.Y();par[5]=mom.Z();
  return TVectorT<double>(6,par);
}


TMatrixDSym
LSLTrackRep::getGlobalCov(){ // covariances
  TMatrixD L(6,5);
  double xp=fState[2];
  double yp=fState[3];
  double no=xp*xp+yp*yp+1;
  double sq=sqrt(no);
  double sq3inv=1/(sq*sq*sq);
  double q=getCharge();
  double p=q/fState[4];
  L[0][0]=1;
  L[1][1]=1;
  L[3][2]=p*(1+yp*yp)*sq3inv;
  L[3][3]=-p*xp*yp*sq3inv;
  L[3][4]=-p/(2.*fState[4])*xp/sq;
  L[4][2]=p*(1+xp*xp)*sq3inv;
  L[4][3]=-p*xp*yp*sq3inv;
  L[4][4]=-p/(2.*fState[4])*yp/sq;
  L[5][2]=-p*xp*sq3inv;
  L[5][3]=-p*yp*sq3inv;
  L[5][4]=-p/(2.*fState[4])/sq;
  
  // calculate new cov;
  TMatrixDSym result(fCov);
  result.Similarity(L);

  // set sigma_z by hand: 
  result[2][2]=0.01;

  return result;
  
}
	      
ClassImp(LSLTrackRep)
