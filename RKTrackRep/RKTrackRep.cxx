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

/* The Runge Kutta implementation stems from GEANT3 originally (R. Brun et al.)
   Porting to C goes back to Igor Gavrilenko @ CERN
   The code was taken from the Phast analysis package of the COMPASS experiment
   (Sergei Gerrassimov @ CERN)
*/

#include"RKTrackRep.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include"assert.h"
#include "stdlib.h"
#include"math.h"
#include"TMath.h"
#include"TGeoManager.h"
#include"TDatabasePDG.h"
#include"GFException.h"
#include"GFFieldManager.h"
#include"GFMaterialEffects.h"

#define MINSTEP 0.001   // minimum step [cm] for Runge Kutta and iteration to POCA


void RKTrackRep::setData(const TMatrixT<double>& st, const GFDetPlane& pl, const TMatrixT<double>* cov, const TMatrixT<double>* aux){
  if(aux != NULL) {
    fCacheSpu = (*aux)(0,0);
  } else {
    if(pl!=fCachePlane){
      std::cerr << "RKTrackRep::setData() - a fatal error ocurred! It was called with a reference plane which is not the same as the one from the last extrapolate(plane,state,cov)-> abort in line " << __LINE__ << std::endl;
      throw;
	}
  }
  GFAbsTrackRep::setData(st,pl,cov);
  fSpu=fCacheSpu;
}

const TMatrixT<double>* RKTrackRep::getAuxInfo(const GFDetPlane& pl) {

  if(pl!=fCachePlane) {
    std::cerr << "RKTrackRep::getAuxInfo() - Fatal error: Trying to get auxillary information with planes mismatch (Information returned does not belong to requested plane)! -> abort in line " << __LINE__ << std::endl;
	throw;
  }
  fAuxInfo.ResizeTo(1,1);
  fAuxInfo(0,0) = fCacheSpu;
  return &fAuxInfo;

}

RKTrackRep::~RKTrackRep(){
  //GFMaterialEffects is now a Singleton
}


RKTrackRep::RKTrackRep() : GFAbsTrackRep(5), fDirection(true), fPdg(0), fMass(0.), fCharge(-1), fCachePlane(), fCacheSpu(1), fSpu(1), fAuxInfo(1,1) {
}


RKTrackRep::RKTrackRep(const TVector3& pos,
                       const TVector3& mom,
                       const TVector3& poserr,
                       const TVector3& momerr,
                       const int& PDGCode) :
                       GFAbsTrackRep(5), fDirection(true), fCachePlane(), fCacheSpu(1), fAuxInfo(1,1) {
  setPDG(PDGCode); // also sets charge and mass
                         

  fRefPlane.setO(pos);
  fRefPlane.setNormal(mom);

  fState[0][0]=fCharge/mom.Mag();
  //u' and v'
  fState[1][0]=0.0;
  fState[2][0]=0.0;
  //u and v
  fState[3][0]=0.0;
  fState[4][0]=0.0;
  //spu
  fSpu=1.;

  TVector3 o=fRefPlane.getO();
  TVector3 u=fRefPlane.getU();
  TVector3 v=fRefPlane.getV();
  TVector3 w=u.Cross(v);
  double pu=0.;
  double pv=0.;
  double pw=mom.Mag();

  fCov[3][3] = poserr.X()*poserr.X() * u.X()*u.X() +
               poserr.Y()*poserr.Y() * u.Y()*u.Y() +
               poserr.Z()*poserr.Z() * u.Z()*u.Z();
  fCov[4][4] = poserr.X()*poserr.X() * v.X()*v.X() +
               poserr.Y()*poserr.Y() * v.Y()*v.Y() +
               poserr.Z()*poserr.Z() * v.Z()*v.Z();
  fCov[0][0] = fCharge*fCharge/pow(mom.Mag(),6.) * 
               (mom.X()*mom.X() * momerr.X()*momerr.X()+
                mom.Y()*mom.Y() * momerr.Y()*momerr.Y()+
                mom.Z()*mom.Z() * momerr.Z()*momerr.Z());
  fCov[1][1] = pow((u.X()/pw - w.X()*pu/(pw*pw)),2.)*momerr.X()*momerr.X() +
               pow((u.Y()/pw - w.Y()*pu/(pw*pw)),2.)*momerr.Y()*momerr.Y() +
               pow((u.Z()/pw - w.Z()*pu/(pw*pw)),2.)*momerr.Z()*momerr.Z();
  fCov[2][2] = pow((v.X()/pw - w.X()*pv/(pw*pw)),2.)*momerr.X()*momerr.X() +
               pow((v.Y()/pw - w.Y()*pv/(pw*pw)),2.)*momerr.Y()*momerr.Y() +
               pow((v.Z()/pw - w.Z()*pv/(pw*pw)),2.)*momerr.Z()*momerr.Z();

}

RKTrackRep::RKTrackRep(const TVector3& pos,
                       const TVector3& mom,
                       const int& PDGCode) :
                       GFAbsTrackRep(5),fDirection(true), fCachePlane(), fCacheSpu(1), fAuxInfo(1,1){
  setPDG(PDGCode); // also sets charge and mass
                         

  fRefPlane.setO(pos);
  fRefPlane.setNormal(mom);

  fState[0][0]=fCharge/mom.Mag();
  //u' and v'
  fState[1][0]=0.0;
  fState[2][0]=0.0;
  //u and v
  fState[3][0]=0.0;
  fState[4][0]=0.0;
  //spu
  fSpu=1.;

  TVector3 o=fRefPlane.getO();
  TVector3 u=fRefPlane.getU();
  TVector3 v=fRefPlane.getV();
  TVector3 w=u.Cross(v);
  double pu=0.;
  double pv=0.;
  double pw=mom.Mag();

  static const TVector3 stdPosErr(1.,1.,1.);
  static const TVector3 stdMomErr(1.,1.,1.);

  fCov[3][3] = stdPosErr.X()*stdPosErr.X() * u.X()*u.X() +
               stdPosErr.Y()*stdPosErr.Y() * u.Y()*u.Y() +
               stdPosErr.Z()*stdPosErr.Z() * u.Z()*u.Z();
  fCov[4][4] = stdPosErr.X()*stdPosErr.X() * v.X()*v.X() +
               stdPosErr.Y()*stdPosErr.Y() * v.Y()*v.Y() +
               stdPosErr.Z()*stdPosErr.Z() * v.Z()*v.Z();
  fCov[0][0] = fCharge*fCharge/pow(mom.Mag(),6.) * 
               (mom.X()*mom.X() * stdMomErr.X()*stdMomErr.X()+
                mom.Y()*mom.Y() * stdMomErr.Y()*stdMomErr.Y()+
                mom.Z()*mom.Z() * stdMomErr.Z()*stdMomErr.Z());
  fCov[1][1] = pow((u.X()/pw - w.X()*pu/(pw*pw)),2.)*stdMomErr.X()*stdMomErr.X() +
               pow((u.Y()/pw - w.Y()*pu/(pw*pw)),2.)*stdMomErr.Y()*stdMomErr.Y() +
               pow((u.Z()/pw - w.Z()*pu/(pw*pw)),2.)*stdMomErr.Z()*stdMomErr.Z();
  fCov[2][2] = pow((v.X()/pw - w.X()*pv/(pw*pw)),2.)*stdMomErr.X()*stdMomErr.X() +
               pow((v.Y()/pw - w.Y()*pv/(pw*pw)),2.)*stdMomErr.Y()*stdMomErr.Y() +
               pow((v.Z()/pw - w.Z()*pv/(pw*pw)),2.)*stdMomErr.Z()*stdMomErr.Z();

}

RKTrackRep::RKTrackRep(const GFDetPlane& pl,
                       const TVector3& mom,
                       const int& PDGCode) :
                       GFAbsTrackRep(5),fDirection(true), fCachePlane(), fCacheSpu(1), fAuxInfo(1,1){
  setPDG(PDGCode); // also sets charge and mass
                         

  fRefPlane = pl;
  TVector3 o=fRefPlane.getO();
  TVector3 u=fRefPlane.getU();
  TVector3 v=fRefPlane.getV();
  TVector3 w=u.Cross(v);

  double pu=mom*u;
  double pv=mom*v;
  double pw=mom*w;

  fState[0][0]=fCharge/mom.Mag();
  //u' and v'
  fState[1][0]=pu/pw;
  fState[2][0]=pv/pw;
  //u and v
  fState[3][0]=0.0;
  fState[4][0]=0.0;
  //spu
  fSpu=1.;

  static const TVector3 stdPosErr2(1.,1.,1.);
  static const TVector3 stdMomErr2(1.,1.,1.);

  fCov[3][3] = stdPosErr2.X()*stdPosErr2.X() * u.X()*u.X() +
               stdPosErr2.Y()*stdPosErr2.Y() * u.Y()*u.Y() +
               stdPosErr2.Z()*stdPosErr2.Z() * u.Z()*u.Z();
  fCov[4][4] = stdPosErr2.X()*stdPosErr2.X() * v.X()*v.X() +
               stdPosErr2.Y()*stdPosErr2.Y() * v.Y()*v.Y() +
               stdPosErr2.Z()*stdPosErr2.Z() * v.Z()*v.Z();
  fCov[0][0] = fCharge*fCharge/pow(mom.Mag(),6.) * 
               (mom.X()*mom.X() * stdMomErr2.X()*stdMomErr2.X()+
                mom.Y()*mom.Y() * stdMomErr2.Y()*stdMomErr2.Y()+
                mom.Z()*mom.Z() * stdMomErr2.Z()*stdMomErr2.Z());
  fCov[1][1] = pow((u.X()/pw - w.X()*pu/(pw*pw)),2.)*stdMomErr2.X()*stdMomErr2.X() +
               pow((u.Y()/pw - w.Y()*pu/(pw*pw)),2.)*stdMomErr2.Y()*stdMomErr2.Y() +
               pow((u.Z()/pw - w.Z()*pu/(pw*pw)),2.)*stdMomErr2.Z()*stdMomErr2.Z();
  fCov[2][2] = pow((v.X()/pw - w.X()*pv/(pw*pw)),2.)*stdMomErr2.X()*stdMomErr2.X() +
               pow((v.Y()/pw - w.Y()*pv/(pw*pw)),2.)*stdMomErr2.Y()*stdMomErr2.Y() +
               pow((v.Z()/pw - w.Z()*pv/(pw*pw)),2.)*stdMomErr2.Z()*stdMomErr2.Z();

}


void RKTrackRep::setPDG(int i){
  fPdg = i;
  TParticlePDG * part = TDatabasePDG::Instance()->GetParticle(fPdg);
  if(part == 0){
    std::cerr << "RKTrackRep::setPDG particle " << i 
              << " not known to TDatabasePDG -> abort" << std::endl;
    exit(1);
  }
  fMass = part->Mass();
  fCharge = part->Charge()/(3.);
}


TVector3 RKTrackRep::getPos(const GFDetPlane& pl){
  if(pl!=fRefPlane){
    TMatrixT<double> s(5,1);
    extrapolate(pl,s);
    return pl.getO()+s[3][0]*pl.getU()+s[4][0]*pl.getV();
  }
  return fRefPlane.getO()+fState[3][0]*fRefPlane.getU()+fState[4][0]*fRefPlane.getV();
}


TVector3 RKTrackRep::getMom(const GFDetPlane& pl){
  TMatrixT<double> statePred(fState);
  TVector3 retmom;
  if(pl!=fRefPlane) {
    extrapolate(pl,statePred);
    retmom = fCacheSpu*(pl.getNormal()+statePred[1][0]*pl.getU()+statePred[2][0]*pl.getV());
  }
  else{
    retmom = fSpu*(pl.getNormal()+statePred[1][0]*pl.getU()+statePred[2][0]*pl.getV());
  }
  retmom.SetMag(1./fabs(statePred[0][0]));
  return retmom;
}


void RKTrackRep::getPosMom(const GFDetPlane& pl,TVector3& pos,
                           TVector3& mom){
  TMatrixT<double> statePred(fState);
  if(pl!=fRefPlane) {
    extrapolate(pl,statePred);
    mom = fCacheSpu*(pl.getNormal()+statePred[1][0]*pl.getU()+statePred[2][0]*pl.getV());
  }
  else{
    mom = fSpu*(pl.getNormal()+statePred[1][0]*pl.getU()+statePred[2][0]*pl.getV());
  }
  mom.SetMag(1./fabs(statePred[0][0]));
  pos = pl.getO()+(statePred[3][0]*pl.getU())+(statePred[4][0]*pl.getV());
}




void RKTrackRep::extrapolateToPoint(const TVector3& pos,
				   TVector3& poca,
				   TVector3& dirInPoca){

  static const int maxIt(30);

  TVector3 o=fRefPlane.getO();
  TVector3 u=fRefPlane.getU();
  TVector3 v=fRefPlane.getV();
  TVector3 w=u.Cross(v);

  TVector3 pTilde = fSpu* (w + fState[1][0] * u + fState[2][0] * v);
  pTilde.SetMag(1.);

  TVector3 point = o + fState[3][0]*u + fState[4][0]*v;

  TMatrixT<double> state7(7,1);
  state7[0][0] = point.X();
  state7[1][0] = point.Y();
  state7[2][0] = point.Z();
  state7[3][0] = pTilde.X();
  state7[4][0] = pTilde.Y();
  state7[5][0] = pTilde.Z();
  state7[6][0] = fState[0][0];

  double coveredDistance(0.);

  GFDetPlane pl;
  int iterations(0);

  while(true){
    pl.setON(pos,TVector3(state7[3][0],state7[4][0],state7[5][0]));
    coveredDistance =  this->Extrap(pl,&state7);

    if(fabs(coveredDistance)<MINSTEP) break;
    if(++iterations == maxIt) {
      GFException exc("RKTrackRep::extrapolateToPoint ==> extrapolation to point failed, maximum number of iterations reached",__LINE__,__FILE__);
      throw exc;
    }
  }
  poca.SetXYZ(state7[0][0],state7[1][0],state7[2][0]);
  dirInPoca.SetXYZ(state7[3][0],state7[4][0],state7[5][0]);
}




TVector3 RKTrackRep::poca2Line(const TVector3& extr1,const TVector3& extr2,const TVector3& point) const {
  
  TVector3 theWire = extr2-extr1;
  if(theWire.Mag()<1.E-8){
    GFException exc("RKTrackRep::poca2Line ==> try to find poca between line and point, but the line is really just a point",__LINE__,__FILE__);
    throw exc;
  }
  double t = 1./(theWire*theWire)*(point*theWire+extr1*extr1-extr1*extr2);
  return (extr1+t*theWire);
}




void RKTrackRep::extrapolateToLine(const TVector3& point1,
                                   const TVector3& point2,
                                   TVector3& poca,
                                   TVector3& dirInPoca,
                                   TVector3& poca_onwire){
  static const int maxIt(30);

  TVector3 o=fRefPlane.getO();
  TVector3 u=fRefPlane.getU();
  TVector3 v=fRefPlane.getV();
  TVector3 w=u.Cross(v);

  TVector3 pTilde = fSpu* (w + fState[1][0] * u + fState[2][0] * v);
  pTilde.SetMag(1.);

  TVector3 point = o + fState[3][0]*u + fState[4][0]*v;

  TMatrixT<double> state7(7,1);
  state7[0][0] = point.X();
  state7[1][0] = point.Y();
  state7[2][0] = point.Z();
  state7[3][0] = pTilde.X();
  state7[4][0] = pTilde.Y();
  state7[5][0] = pTilde.Z();
  state7[6][0] = fState[0][0];

  double coveredDistance(0.);

  GFDetPlane pl;
  int iterations(0);

  while(true){
    pl.setO(point1);
    TVector3 currentDir(state7[3][0],state7[4][0],state7[5][0]);
    pl.setU(currentDir.Cross(point2-point1));
    pl.setV(point2-point1);
    coveredDistance = this->Extrap(pl,&state7);

    if(fabs(coveredDistance)<MINSTEP) break;
    if(++iterations == maxIt) {
      GFException exc("RKTrackRep::extrapolateToLine ==> extrapolation to line failed, maximum number of iterations reached",__LINE__,__FILE__);
      throw exc;
    }
  }
  poca.SetXYZ(state7[0][0],state7[1][0],state7[2][0]);
  dirInPoca.SetXYZ(state7[3][0],state7[4][0],state7[5][0]);
  poca_onwire = poca2Line(point1,point2,poca);
}




double RKTrackRep::extrapolate(const GFDetPlane& pl, 
                               TMatrixT<double>& statePred,
                               TMatrixT<double>& covPred){
  
  TMatrixT<double> cov7x7(7,7);
  TMatrixT<double> J_pM(7,5);

  TVector3 o=fRefPlane.getO();
  TVector3 u=fRefPlane.getU();
  TVector3 v=fRefPlane.getV();
  TVector3 w=u.Cross(v);

  J_pM[0][3] = u.X();J_pM[0][4]=v.X(); // dx/du
  J_pM[1][3] = u.Y();J_pM[1][4]=v.Y();
  J_pM[2][3] = u.Z();J_pM[2][4]=v.Z();

  TVector3 pTilde = fSpu * (w + fState[1][0] * u + fState[2][0] * v);
  double pTildeMag = pTilde.Mag();

  //J_pM matrix is d(x,y,z,ax,ay,az,q/p) / d(q/p,u',v',u,v)

  // da_x/du'
  J_pM[3][1] = fSpu/pTildeMag*(u.X()-pTilde.X()/(pTildeMag*pTildeMag)*u*pTilde);
  J_pM[4][1] = fSpu/pTildeMag*(u.Y()-pTilde.Y()/(pTildeMag*pTildeMag)*u*pTilde);
  J_pM[5][1] = fSpu/pTildeMag*(u.Z()-pTilde.Z()/(pTildeMag*pTildeMag)*u*pTilde);
  // da_x/dv'
  J_pM[3][2] = fSpu/pTildeMag*(v.X()-pTilde.X()/(pTildeMag*pTildeMag)*v*pTilde);
  J_pM[4][2] = fSpu/pTildeMag*(v.Y()-pTilde.Y()/(pTildeMag*pTildeMag)*v*pTilde);
  J_pM[5][2] = fSpu/pTildeMag*(v.Z()-pTilde.Z()/(pTildeMag*pTildeMag)*v*pTilde);
  // dqOp/dqOp
  J_pM[6][0] = 1.;

  TMatrixT<double> J_pM_transp(J_pM);
  J_pM_transp.T();

  cov7x7 = J_pM*(fCov*J_pM_transp);

  TVector3 pos = o + fState[3][0]*u + fState[4][0]*v;
  TMatrixT<double> state7(7,1);
  state7[0][0] = pos.X();
  state7[1][0] = pos.Y();
  state7[2][0] = pos.Z();
  state7[3][0] = pTilde.X()/pTildeMag;;
  state7[4][0] = pTilde.Y()/pTildeMag;;
  state7[5][0] = pTilde.Z()/pTildeMag;;
  state7[6][0] = fState[0][0];

  double coveredDistance = this->Extrap(pl,&state7,&cov7x7);

  TVector3 O = pl.getO();
  TVector3 U = pl.getU();
  TVector3 V = pl.getV();
  TVector3 W = pl.getNormal();

  double X = state7[0][0];
  double Y = state7[1][0];
  double Z = state7[2][0];
  double AX = state7[3][0];
  double AY = state7[4][0];
  double AZ = state7[5][0];
  double QOP = state7[6][0];
  TVector3 A(AX,AY,AZ);
  TVector3 Point(X,Y,Z);
  TMatrixT<double> J_Mp(5,7);
  
  // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,ax,ay,az,q/p)
  J_Mp[0][6] = 1.;
  //du'/da_x
  double AtW = A*W;
  J_Mp[1][3] = (U.X()*(AtW)-W.X()*(A*U))/(AtW*AtW);
  J_Mp[1][4] = (U.Y()*(AtW)-W.Y()*(A*U))/(AtW*AtW);
  J_Mp[1][5] = (U.Z()*(AtW)-W.Z()*(A*U))/(AtW*AtW);
  //dv'/da_x
  J_Mp[2][3] = (V.X()*(AtW)-W.X()*(A*V))/(AtW*AtW);
  J_Mp[2][4] = (V.Y()*(AtW)-W.Y()*(A*V))/(AtW*AtW);
  J_Mp[2][5] = (V.Z()*(AtW)-W.Z()*(A*V))/(AtW*AtW);
  //du/dx
  J_Mp[3][0] = U.X();
  J_Mp[3][1] = U.Y();
  J_Mp[3][2] = U.Z();
  //dv/dx
  J_Mp[4][0] = V.X();
  J_Mp[4][1] = V.Y();
  J_Mp[4][2] = V.Z();
  
  TMatrixT<double> J_Mp_transp(J_Mp);
  J_Mp_transp.T();

  covPred.ResizeTo(5,5);

  covPred = J_Mp*(cov7x7*J_Mp_transp);

  statePred.ResizeTo(5,1);
  statePred[0][0] = QOP;
  statePred[1][0] = (A*U)/(A*W);
  statePred[2][0] = (A*V)/(A*W);
  statePred[3][0] = (Point-O)*U;
  statePred[4][0] = (Point-O)*V;
  
  fCachePlane = pl;
  fCacheSpu = (A*W)/fabs(A*W);

  return coveredDistance;
}




double RKTrackRep::extrapolate(const GFDetPlane& pl, 
                               TMatrixT<double>& statePred){

  TVector3 o=fRefPlane.getO();
  TVector3 u=fRefPlane.getU();
  TVector3 v=fRefPlane.getV();
  TVector3 w=u.Cross(v);


  TVector3 pTilde = fSpu* (w + fState[1][0] * u + fState[2][0] * v);
  double pTildeMag = pTilde.Mag();

  TVector3 pos = o + fState[3][0]*u + fState[4][0]*v;

  TMatrixT<double> state7(7,1);
  state7[0][0] = pos.X();
  state7[1][0] = pos.Y();
  state7[2][0] = pos.Z();
  state7[3][0] = pTilde.X()/pTildeMag;
  state7[4][0] = pTilde.Y()/pTildeMag;
  state7[5][0] = pTilde.Z()/pTildeMag;
  state7[6][0] = fState[0][0];

  TVector3 O = pl.getO();
  TVector3 U = pl.getU();
  TVector3 V = pl.getV();
  TVector3 W = pl.getNormal();

  double coveredDistance = this->Extrap(pl,&state7);

  double X = state7[0][0];
  double Y = state7[1][0];
  double Z = state7[2][0];
  double AX = state7[3][0];
  double AY = state7[4][0];
  double AZ = state7[5][0];
  double QOP = state7[6][0];
  TVector3 A(AX,AY,AZ);
  TVector3 Point(X,Y,Z);

  statePred.ResizeTo(5,1);
  statePred[0][0] = QOP;
  statePred[1][0] = (A*U)/(A*W);
  statePred[2][0] = (A*V)/(A*W);
  statePred[3][0] = (Point-O)*U;
  statePred[4][0] = (Point-O)*V;

  return coveredDistance;
}




//
// Runge-Kutta method for tracking a particles through a magnetic field.            
// Uses Nystroem algorithm (See Handbook Nat. Bur. of Standards, procedure 25.5.20)  
//                                                                                  
// Input parameters:                                                               
//    SU     - plane parameters                                                                                                                         
//    SU[0]  - direction cosines normal to surface Ex                               
//    SU[1]  -          -------                    Ey                               
//    SU[2]  -          -------                    Ez; Ex*Ex+Ey*Ey+Ez*Ez=1          
//    SU[3]  - distance to surface from (0,0,0) > 0 cm                                 
//
//    ND     - number of variables for derivatives calculation                         
//    P      - initial parameters (coordinates(cm), direction cosines,              
//             charge/momentum (Gev-1) and derivatives this parameters  (8x7)            
//         
//    X        	Y        	Z        	Ax       	Ay       	Az       	q/P                   
//    P[ 0]   	P[ 1]   	P[ 2]   	P[ 3]   	P[ 4]   	P[ 5]   	P[ 6]  
//
//    dX/dp    	dY/dp    	dZ/dp    	dAx/dp   	dAy/dp   	dAz/dp   	d(q/P)/dp*P[6]         
//    P[ 7]   	P[ 8]   	P[ 9]   	P[10]   	P[11]   	P[12]   	P[13]   			      d()/dp1  
//
//    P[14]   	P[15]   	P[16]   	P[17]   	P[18]   	P[19]   	P[20]   		      	d()/dp2        
//    ............................................................................		d()/dpND       
//                                                                                  
// Output parameters:                                                               
//                                                                                  
//    P    -  output parameters and derivatives after propagation in magnetic field 
//            defined by Mfield (KGauss)                                            
//    Where a Mfield(R,H) - is interface to magnetic field information              
//    input  	R[ 0],R[ 1],R[ 2] - X     , Y      and Z  of the track                 
//    output 	H[ 0],H[ 1],H[ 2] - Hx    , Hy     and Hz of the magnetic field        
//           	H[ 3],H[ 4],H[ 5] - dHx/dx, dHx/dy and dHx/dz          //                
//           	H[ 6],H[ 7],H[ 8] - dHy/dx, dHy/dy and dHy/dz          // (not used)
//           	H[ 9],H[10],H[11] - dHz/dx, dHz/dy and dHz/dz          //                
//                                                                                  
// Authors: R.Brun, M.Hansroul, V.Perevoztchikov (Geant3)                           
//  
bool RKTrackRep::RKutta (const GFDetPlane& plane,
                         double* P, 
                         double& coveredDistance, 
                         std::vector<TVector3>& points,
                         std::vector<double>& pointPaths, 
                         const double& maxLen,  // currently not used
                         bool calcCov) const {

  static const double EC     = .000149896229;   // c/(2*10^12) resp. c/2Tera
  static const double DLT    = .0002;           // max. deviation for approximation-quality test
  static const double DLT32  = DLT/32.;         //
  static const double P3     = 1./3.;           // 1/3
  static const double Smax   = 100.;            // max. step allowed > 0 
  static const double Wmax   = 2000.;           // max. way allowed
  static const double Pmin   = 4.E-3;           // minimum momentum for propagation [GeV]
  static const int    ND     = 56;              // number of variables for derivatives calculation
  static const int    ND1    = ND-7;            // = 49
  double* R           = &P[0];                  // Start coordinates  in cm 	( x,  y,  z)
  double* A           = &P[3];                  // Start directions 	      	(ax, ay, az); 	ax^2+ay^2+az^2=1
  double  SA[3]       = {0.,0.,0.};             // Start directions derivatives 
  double  Pinv        = P[6]*EC;                // P[6] is charge/momentum in e/(Gev/c)
  double  Way         = 0.;                     // Total way of the trajectory
  double  Way2        = 0.;                     // Total way of the trajectory with correct signs
  bool    error       = false;                  // Error of propogation
  bool    stopBecauseOfMaterial = false;        // does not go through main loop again when stepsize is reduced by stepper

  points.clear();
  pointPaths.clear();
  //std::cout<<"coords  "<<R[0]<<"  "<<R[1]<<"   "<<R[2]<<  std::endl;
  //std::cout<<"R       "<<sqrt(pow(R[0],2)+pow(R[1],2))<<  std::endl;  
  //std::cout<<"momentum "<<fabs(fCharge/P[6])<< std::endl;
  if(fabs(fCharge/P[6])<Pmin){
    std::cerr << "RKTrackRep::RKutta ==> momentum too low: " << fabs(fCharge/P[6])*1000. << " MeV" << std::endl;
    return (false);
  }
  
  double SU[4];
  TVector3 O = plane.getO();
  TVector3 W = plane.getNormal(); 
  if(W*O > 0){ 		// make SU vector point away from origin
    SU[0] = W.X();
    SU[1] = W.Y();
    SU[2] = W.Z();
  }
  else{
    SU[0] = -1.*W.X();
    SU[1] = -1.*W.Y();
    SU[2] = -1.*W.Z();
  }
  SU[3] = plane.distance(0.,0.,0.);

  //
  // Step estimation until surface
  //
  double Step,An,Dist,Dis,S,Sl=0;

  points.push_back(TVector3(R[0],R[1],R[2]));
  pointPaths.push_back(0.);

  An=A[0]*SU[0]+A[1]*SU[1]+A[2]*SU[2];		// An = A * N;  component of A normal to surface
  if(fabs(An) < 1.E-6) {
    std::cerr << "RKTrackRep::RKutta ==> cannot propagate perpendicular to plane " << std::endl;
    return false;		          	// no propagation if A perpendicular to surface
  }
  if( plane.inActive(TVector3(R[0],R[1],R[2]),TVector3(A[0],A[1],A[2]))) {  // if direction is not pointing to active part of surface
    Dist=SU[3]-R[0]*SU[0]-R[1]*SU[1]-R[2]*SU[2];				// Distance between start coordinates and surface
    Step=Dist/An;
  }
  else{									                  // make sure to extrapolate towards surface
    if( (O.X()-R[0])*A[0] + (O.Y()-R[1])*A[1] + (O.Z()-R[2])*A[2] >0 ){	  	// if direction A pointing from start coordinates R towards surface
      Dist = sqrt((R[0]-O.X())*(R[0]-O.X())+					 // |R-O|; Distance between start coordinates and origin of surface
                  (R[1]-O.Y())*(R[1]-O.Y())+
                  (R[2]-O.Z())*(R[2]-O.Z()));      
    }
    else{									                // if direction pointing away from surface
      Dist = -1.*sqrt((R[0]-O.X())*(R[0]-O.X())+
                      (R[1]-O.Y())*(R[1]-O.Y())+
                      (R[2]-O.Z())*(R[2]-O.Z()));            
    }
    Step=Dist;    
  }

  if(fabs(Step)>Wmax) {
    std::cerr<<"RKTrackRep::RKutta ==> Too long extrapolation requested : "<<Step<<" cm !"<<std::endl;
    std::cerr<<"X = "<<R[0]<<" Y = "<<R[1]<<" Z = "<<R[2]
             <<"  COSx = "<<A[0]<<"  COSy = "<<A[1]<<"  COSz = "<<A[2]<<std::endl;
    std::cout<<"Destination  X = "<<SU[0]*SU[3]<<std::endl;
    return(false);
  }

  // reduce maximum stepsize S to Smax
  Step>Smax ? S=Smax : Step<-Smax ? S=-Smax : S=Step;	
  
  //
  // Main cycle of Runge-Kutta method
  //

  //for saving points only when the direction didnt change
  int Ssign=1;
  if(S<0) Ssign = -1;

  while(fabs(Step)>MINSTEP && !stopBecauseOfMaterial) {
    
    // call stepper and reduce stepsize
    double stepperLen;
    stepperLen = GFMaterialEffects::getInstance()->stepper(fabs(S),
                                  R[0],R[1],R[2],
                                  Ssign*A[0],Ssign*A[1],Ssign*A[2],
                                  fabs(fCharge/P[6]),
                                  fPdg);
    if (stepperLen<MINSTEP) stepperLen=MINSTEP; // prevents tiny stepsizes that can slow down the program
    if (S > stepperLen) {
      S = stepperLen;
      stopBecauseOfMaterial = true;
    }
    else if (S < -stepperLen) {
      S = -stepperLen;	
      stopBecauseOfMaterial = true;
    }

    double H0[12],H1[12],H2[12],r[3];
    double S3=P3*S, S4=.25*S, PS2=Pinv*S; 
    
    //
    // First point
    //   
    r[0]=R[0]      ; r[1]=R[1]      ; r[2]=R[2]      ;  
    TVector3 pos(r[0],r[1],r[2]);						// vector of start coordinates R0	(x, y, z)
    TVector3 H0vect = GFFieldManager::getFieldVal(pos);				// magnetic field in 10^-4 T = kGauss
    H0[0]=PS2*H0vect.X(); H0[1]=PS2*H0vect.Y(); H0[2]=PS2*H0vect.Z(); 		// H0 is PS2*(Hx, Hy, Hz) @ R0
    double A0=A[1]*H0[2]-A[2]*H0[1], B0=A[2]*H0[0]-A[0]*H0[2], C0=A[0]*H0[1]-A[1]*H0[0]; // (ax, ay, az) x H0
    double A2=A[0]+A0              , B2=A[1]+B0              , C2=A[2]+C0              ; // (A0, B0, C0) + (ax, ay, az)
    double A1=A2+A[0]              , B1=B2+A[1]              , C1=C2+A[2]              ; // (A0, B0, C0) + 2*(ax, ay, az)
      
    //
    // Second point
    //
    r[0]+=A1*S4    ; r[1]+=B1*S4    ; r[2]+=C1*S4    ;   //setup.Field(r,H1);
    pos.SetXYZ(r[0],r[1],r[2]);
    TVector3 H1vect = GFFieldManager::getFieldVal(pos);
    H1[0]=H1vect.X()*PS2; H1[1]=H1vect.Y()*PS2;H1[2]=H1vect.Z()*PS2;	// H1 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * [(A0, B0, C0) + 2*(ax, ay, az)]
    double A3,B3,C3,A4,B4,C4,A5,B5,C5;
    A3 = B2*H1[2]-C2*H1[1]+A[0]; B3=C2*H1[0]-A2*H1[2]+A[1]; C3=A2*H1[1]-B2*H1[0]+A[2]; // (A2, B2, C2) x H1 + (ax, ay, az)
    A4 = B3*H1[2]-C3*H1[1]+A[0]; B4=C3*H1[0]-A3*H1[2]+A[1]; C4=A3*H1[1]-B3*H1[0]+A[2]; // (A3, B3, C3) x H1 + (ax, ay, az)
    A5 = A4-A[0]+A4            ; B5=B4-A[1]+B4            ; C5=C4-A[2]+C4            ; //    2*(A4, B4, C4) - (ax, ay, az)

    //
    // Last point
    //
    r[0]=R[0]+S*A4 ; r[1]=R[1]+S*B4 ; r[2]=R[2]+S*C4 ;  //setup.Field(r,H2);
    pos.SetXYZ(r[0],r[1],r[2]);
    TVector3 H2vect = GFFieldManager::getFieldVal(pos);
    H2[0]=H2vect.X()*PS2; H2[1]=H2vect.Y()*PS2;H2[2]=H2vect.Z()*PS2;	// H2 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * (A4, B4, C4)
    double A6=B5*H2[2]-C5*H2[1], B6=C5*H2[0]-A5*H2[2], C6=A5*H2[1]-B5*H2[0]; // (A5, B5, C5) x H2

    //
    // Test approximation quality on given step and possible step reduction
    //
    double EST = fabs((A1+A6)-(A3+A4))+fabs((B1+B6)-(B3+B4))+fabs((C1+C6)-(C3+C4));  // EST = ||(ABC1+ABC6)-(ABC3+ABC4)||_1  =  ||(axzy x H0 + ABC5 x H2) - (ABC2 x H1 + ABC3 x H1)||_1
    if(EST>DLT) {
      S*=0.5; 
      stopBecauseOfMaterial = false;
      continue;
    }
    
    //
    // Derivatives of track parameters in last point
    //
    if(calcCov){
      for(int i=7; i!=ND; i+=7) {				// i = 7, 14, 21, 28, 35, 42, 49;    ND = 56;	ND1 = 49; rows of Jacobian
	
        double* dR = &P[i];			            		// dR = (dX/dpN,  dY/dpN,  dZ/dpN)
        double* dA = &P[i+3];				           	// dA = (dAx/dpN, dAy/dpN, dAz/dpN); N = X,Y,Z,Ax,Ay,Az,q/p
        
        //first point
        double dA0   = H0[ 2]*dA[1]-H0[ 1]*dA[2];		// dA0/dp	}
        double dB0   = H0[ 0]*dA[2]-H0[ 2]*dA[0];		// dB0/dp	 } = dA x H0	
        double dC0   = H0[ 1]*dA[0]-H0[ 0]*dA[1];		// dC0/dp	}
        
        if(i==ND1) {dA0+=A0; dB0+=B0; dC0+=C0;}			// if last row: (dA0, dB0, dC0) := (dA0, dB0, dC0) + (A0, B0, C0)
        
        double dA2   = dA0+dA[0];				// }
        double dB2   = dB0+dA[1]; 			//  } = (dA0, dB0, dC0) + dA
        double dC2   = dC0+dA[2];				// }
         
        //second point
        double dA3   = dA[0]+dB2*H1[2]-dC2*H1[1];		// dA3/dp	}
        double dB3   = dA[1]+dC2*H1[0]-dA2*H1[2];		// dB3/dp	 } = dA + (dA2, dB2, dC2) x H1	
        double dC3   = dA[2]+dA2*H1[1]-dB2*H1[0];		// dC3/dp	}
        
        if(i==ND1) {dA3+=A3-A[0]; dB3+=B3-A[1]; dC3+=C3-A[2];} // if last row: (dA3, dB3, dC3) := (dA3, dB3, dC3) + (A3, B3, C3) - (ax, ay, az)

        double dA4   = dA[0]+dB3*H1[2]-dC3*H1[1];		// dA4/dp	}
        double dB4   = dA[1]+dC3*H1[0]-dA3*H1[2];		// dB4/dp	 } = dA + (dA3, dB3, dC3) x H1	
        double dC4   = dA[2]+dA3*H1[1]-dB3*H1[0];		// dC4/dp	}
        
        if(i==ND1) {dA4+=A4-A[0]; dB4+=B4-A[1]; dC4+=C4-A[2];} // if last row: (dA4, dB4, dC4) := (dA4, dB4, dC4) + (A4, B4, C4) - (ax, ay, az)
        
        //last point	
        double dA5   = dA4+dA4-dA[0];				// }
        double dB5   = dB4+dB4-dA[1];				//  } =  2*(dA4, dB4, dC4) - dA
        double dC5   = dC4+dC4-dA[2]; 			// }

        double dA6   = dB5*H2[2]-dC5*H2[1];			// dA6/dp	}
        double dB6   = dC5*H2[0]-dA5*H2[2];			// dB6/dp	 } = (dA5, dB5, dC5) x H2	
        double dC6   = dA5*H2[1]-dB5*H2[0];			// dC6/dp	}	

        if(i==ND1) {dA6+=A6; dB6+=B6; dC6+=C6;}			// if last row: (dA6, dB6, dC6) := (dA6, dB6, dC6) + (A6, B6, C6)                                    
        
        dR[0]+=(dA2+dA3+dA4)*S3; dA[0] = (dA0+dA3+dA3+dA5+dA6)*P3;	// dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]      
        dR[1]+=(dB2+dB3+dB4)*S3; dA[1] = (dB0+dB3+dB3+dB5+dB6)*P3;	// dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
        dR[2]+=(dC2+dC3+dC4)*S3; dA[2] = (dC0+dC3+dC3+dC5+dC6)*P3;
      }
    }
    
    Way2 += S;				// add stepsize to way (signed)
    if((Way+=fabs(S))>Wmax){ 
      std::cerr<<"RKTrackRep::RKutta ==> Trajectory is longer than length limit : "<<Way<<" cm !"
      << " p/q = "<<1./P[6]<< " GeV"<<std::endl;
      return(false);
    }
    
    //
    // Track parameters in last point
    //   
    R[0]+=(A2+A3+A4)*S3; A[0]+=(SA[0]=(A0+A3+A3+A5+A6)*P3-A[0]);  // R  = R0 + S3*[(A2, B2, C2) +   (A3, B3, C3) + (A4, B4, C4)] 
    R[1]+=(B2+B3+B4)*S3; A[1]+=(SA[1]=(B0+B3+B3+B5+B6)*P3-A[1]);  // A  =     1/3*[(A0, B0, C0) + 2*(A3, B3, C3) + (A5, B5, C5) + (A6, B6, C6)]
    R[2]+=(C2+C3+C4)*S3; A[2]+=(SA[2]=(C0+C3+C3+C5+C6)*P3-A[2]); 	// SA = A_new - A_old
    Sl=S;	// last S used
    
    // if extrapolation has changed direction, delete the last point, because it is
    // not a consecutive point to be used for material estimations
    if(Ssign*S<0.) {
      pointPaths.at(pointPaths.size()-1)+=S;
      points.erase(points.end());
    }
    else{
      pointPaths.push_back(S);
    }

    points.push_back(TVector3(R[0],R[1],R[2]));

    double CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);	// 1/|A|
    A[0]*=CBA; A[1]*=CBA; A[2]*=CBA;				// normalize A
  
    // Step estimation until surface and test conditions for stop of propogation
    if(fabs(Way2)>Wmax) {
      Dis=0.;
      Dist=0.;
      S=0;
      Step=0.;
      break;
    }
    

    An=A[0]*SU[0]+A[1]*SU[1]+A[2]*SU[2];

    if(fabs(An) < 1.E-6) {
      error=true; 
      Step=0; 
      break;
    }
    
    if( plane.inActive(TVector3(R[0],R[1],R[2]),TVector3(A[0],A[1],A[2]))) {
      Dis=SU[3]-R[0]*SU[0]-R[1]*SU[1]-R[2]*SU[2];
      Step=Dis/An; 
    }
    else{
      if( (O.X()-R[0])*A[0] + (O.Y()-R[1])*A[1] + (O.Z()-R[2])*A[2] >0 ){
        Dis = sqrt((R[0]-O.X())*(R[0]-O.X())+
                   (R[1]-O.Y())*(R[1]-O.Y())+
                   (R[2]-O.Z())*(R[2]-O.Z()));      
      }
      else{
        Dis = -1.*sqrt((R[0]-O.X())*(R[0]-O.X())+
                       (R[1]-O.Y())*(R[1]-O.Y())+
                       (R[2]-O.Z())*(R[2]-O.Z()));            
      }
      Step = Dis; // signed distance to surface
    }

    if (Dis*Dist>0 && fabs(Dis)>fabs(Dist)) { // did not get closer to surface
      error=true; 
      Step=0; 
      break;
    }
    Dist=Dis;

    //
    // reset & check step size
    //
    // reset S to Step if extrapolation too long or in wrong direction
    if (S*Step<0. || fabs(S)>fabs(Step)) S=Step;
    else if (EST<DLT32 && fabs(2.*S)<=Smax) S*=2.;     
    
  } //end of main loop
  
  //
  // Output information preparation for main track parameteres
  //
  
  if (!stopBecauseOfMaterial) { // linear extrapolation to surface
    if(fabs(Sl) > 1.E-12) Sl=1./Sl;	      // Sl = inverted last Stepsize Sl
    A [0]+=(SA[0]*=Sl)*Step; 	// Step  = distance to surface
    A [1]+=(SA[1]*=Sl)*Step; 	// SA*Sl = delta A / delta way; local derivative of A with respect to the length of the way
    A [2]+=(SA[2]*=Sl)*Step;	// A = A + Step * SA*Sl

    P[0]      = R[0]+Step*(A[0]-.5*Step*SA[0]);    // P = R + Step*(A - 1/2*Step*SA); approximation for final point on surface                           
    P[1]      = R[1]+Step*(A[1]-.5*Step*SA[1]);
    P[2]      = R[2]+Step*(A[2]-.5*Step*SA[2]);
      
    points.push_back(TVector3(P[0],P[1],P[2]));
    pointPaths.push_back(Step);
  }
  
  double CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
  
  P[3]      = A[0]*CBA;	// normalize A
  P[4]      = A[1]*CBA;
  P[5]      = A[2]*CBA;
   
  //
  // Output derivatives of track parameters preparation 
  //
  An = A[0]*SU[0]+A[1]*SU[1]+A[2]*SU[2]; 
  fabs(An) < 1.E-6 ? An=1./An : An = 0; // 1/A_normal
  
  if(calcCov && !stopBecauseOfMaterial){
    for(int i=7; i!=ND; i+=7) {
      double* dR = &P[i];  double* dA = &P[i+3];	
      S = (dR[0]*SU[0]+dR[1]*SU[1]+dR[2]*SU[2])*An;	// dR_normal / A_normal
      dR[0]-=S*A [0];  dR[1]-=S*A [1]; dR[2]-=S*A [2]; 
      dA[0]-=S*SA[0];  dA[1]-=S*SA[1]; dA[2]-=S*SA[2]; 
    }
  }
  
  if(error){
    std::cerr << "RKTrackRep::RKutta ==> Do not get closer. Path = " << Way << " cm" << "  p/q = " << 1./P[6] << " GeV" << std::endl;
    return(false);
  }

  // calculate covered distance
  if (!stopBecauseOfMaterial) coveredDistance=Way2+Step;
  else coveredDistance=Way2;

  return(true);
}




double RKTrackRep::Extrap( const GFDetPlane& plane, TMatrixT<double>* state, TMatrixT<double>* cov) const {

  static const int maxNumIt(2000);
  int numIt(0);
  bool calcCov(true);
  if(cov==NULL) calcCov=false;  

  double *P;
  if(calcCov) {P = new double[56]; memset(P,0x00,56*sizeof(double));}
  else {P = new double[7];} // not needed memset(P,0x00,7*sizeof(double));};

  for(int i=0;i<7;++i){
    P[i] = (*state)[i][0];
  }
  
  TMatrixT<double> jac(7,7);
  TMatrixT<double> jacT(7,7);
  TMatrixT<double> oldCov(7,7);
  if(calcCov) oldCov=(*cov);
  double coveredDistance(0.);
  double sumDistance(0.);

  while(true){
    if(numIt++ > maxNumIt){
      GFException exc("RKTrackRep::Extrap ==> maximum number of iterations exceeded",
		      __LINE__,__FILE__);
      exc.setFatal();
      delete[] P;
      throw exc;
    }

    if(calcCov){
      memset(&P[7],0x00,49*sizeof(double));
      for(int i=0; i<6; ++i){
        P[(i+1)*7+i] = 1.;
      }
      P[55] =  (*state)[6][0];
    }

    double dir(1.);
    {
      TVector3 Pvect(P[0],P[1],P[2]); //position
      TVector3 Avect(P[3],P[4],P[5]); //direction
      TVector3 dist = plane.dist(Pvect); //from point to plane
      if(dist*Avect<0.) dir=-1.;
    }

    TVector3 directionBefore(P[3],P[4],P[5]); // direction before propagation
    directionBefore.SetMag(1.);
    
    // propagation
    std::vector<TVector3> points;
    std::vector<double> pointPaths;
    if( ! this->RKutta(plane, P, coveredDistance, points, pointPaths, -1., calcCov) ) { // maxLen currently not used
      GFException exc("RKTrackRep::Extrap ==> Runge Kutta propagation failed",__LINE__,__FILE__);
      exc.setFatal(); // stops propagation; faster, but some hits will be lost
      delete[] P;
      throw exc;
    }

    TVector3 directionAfter(P[3],P[4],P[5]); // direction after propagation
    directionAfter.SetMag(1.);
      
    sumDistance+=coveredDistance;
     
    // filter Points
    std::vector<TVector3> pointsFilt(1, points.at(0));
    std::vector<double> pointPathsFilt(1, 0.);
    // only if in right direction
    for(unsigned int i=1;i<points.size();++i){ 
      if (pointPaths.at(i) * coveredDistance > 0.) { 
        pointsFilt.push_back(points.at(i));
        pointPathsFilt.push_back(pointPaths.at(i));
      }
      else {
        pointsFilt.back() = points.at(i);
        pointPathsFilt.back() += pointPaths.at(i);
      }
      // clean up tiny steps
      int position = pointsFilt.size()-1;  // position starts with 0
      if (fabs(pointPathsFilt.back()) < MINSTEP && position > 1) {
        pointsFilt.at(position-1) = pointsFilt.at(position);
        pointsFilt.pop_back();        
        pointPathsFilt.at(position-1) += pointPathsFilt.at(position);
        pointPathsFilt.pop_back();
      }
    }
    
    //consistency check
    double checkSum(0.);
    for(unsigned int i=0;i<pointPathsFilt.size();++i){
      checkSum+=pointPathsFilt.at(i);
    }
    if(fabs(checkSum-coveredDistance)>1.E-7){
      GFException exc("RKTrackRep::Extrap ==> fabs(checkSum-coveredDistance)>1.E-7",__LINE__,__FILE__);
      exc.setFatal();
      delete[] P;
      throw exc;
    }
    
    if(calcCov){ //calculate Jacobian jac
      for(int i=0;i<7;++i){
	      for(int j=0;j<7;++j){
	        if(i<6) jac[i][j] = P[ (i+1)*7+j ];
	        else jac[i][j] = P[ (i+1)*7+j ]/P[6];
	      }  
      }
      jacT = jac;
      jacT.T();
    }
    
    TMatrixT<double> noise(7,7); // zero everywhere by default
    
    // call MatEffects
    double momLoss; // momLoss has a sign - negative loss means momentum gain
    
    momLoss = GFMaterialEffects::getInstance()->effects(pointsFilt,
                               pointPathsFilt,
                               fabs(fCharge/P[6]), // momentum
                               fPdg,
                               calcCov,
                               &noise,
                               &jac,
                               &directionBefore,
                               &directionAfter);
  
    if(fabs(P[6])>1.E-10){ // do momLoss only for defined 1/momentum .ne.0      
	    P[6] = fCharge/(fabs(fCharge/P[6])-momLoss);
    }
    
    if(calcCov){ //propagate cov and add noise
      if(!(oldCov < 1.E300)){
        GFException exc("RKTrackRep::Extrap ==> covariance matrix exceeds numerical limits",__LINE__,__FILE__);
        exc.setFatal();
        delete[] P;
        throw exc;
      }
      oldCov = *cov;
      *cov = jacT*((oldCov)*jac)+noise;
    }
    
    
    //we arrived at the destination plane, if we point to the active area
    //of the plane (if it is finite), and the distance is below threshold
    if( plane.inActive(TVector3(P[0],P[1],P[2]),TVector3(P[3],P[4],P[5]))) {
      if(plane.distance(P[0],P[1],P[2])<MINSTEP) break;
    }
  }
  (*state)[0][0] = P[0];  (*state)[1][0] = P[1];
  (*state)[2][0] = P[2];  (*state)[3][0] = P[3];
  (*state)[4][0] = P[4];  (*state)[5][0] = P[5];
  (*state)[6][0] = P[6];
    
  delete[] P;
  return sumDistance;
}

ClassImp(RKTrackRep)
