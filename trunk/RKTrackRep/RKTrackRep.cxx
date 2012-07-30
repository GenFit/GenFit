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

#include "RKTrackRep.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include "assert.h"
#include "stdlib.h"
#include "math.h"
#include "TMath.h"
#include "TGeoManager.h"
#include "TDatabasePDG.h"
#include "GFException.h"
#include "GFFieldManager.h"
#include "GFMaterialEffects.h"

#define MINSTEP 0.001   // minimum step [cm] for Runge Kutta and iteration to POCA
//#define DEBUG


RKTrackRep::~RKTrackRep(){
}


RKTrackRep::RKTrackRep() : GFAbsTrackRep(5), fDirection(0), fNoMaterial(false), fPdg(0), fMass(0.), fCharge(-1), fCachePlane(), fCacheSpu(1), fSpu(1), fAuxInfo(1,2) {
}


RKTrackRep::RKTrackRep(const TVector3& pos,
                       const TVector3& mom,
                       const TVector3& poserr,
                       const TVector3& momerr,
                       const int& PDGCode) :
                       GFAbsTrackRep(5), fDirection(0), fNoMaterial(false), fCachePlane(), fCacheSpu(1), fAuxInfo(1,2) {

  setPDG(PDGCode); // also sets charge and mass
  calcStateCov(pos, mom, poserr, momerr);
}


RKTrackRep::RKTrackRep(const TVector3& pos,
                       const TVector3& mom,
                       const int& PDGCode) :
                       GFAbsTrackRep(5), fDirection(0), fNoMaterial(false), fCachePlane(), fCacheSpu(1), fAuxInfo(1,2) {

  setPDG(PDGCode); // also sets charge and mass
  calcState(pos, mom);

  // set covariance diagonal elements to large number
  static const double value(1.E4);

  fCov[0][0] = value;
  fCov[1][1] = value;
  fCov[2][2] = value;
  fCov[3][3] = value;
  fCov[4][4] = value;
}


RKTrackRep::RKTrackRep(const GFTrackCand* aGFTrackCandPtr) :
                       GFAbsTrackRep(5), fDirection(0), fNoMaterial(false), fCachePlane(), fCacheSpu(1), fAuxInfo(1,2) {

  setPDG(aGFTrackCandPtr->getPdgCode()); // also sets charge and mass

  TVector3 pos = aGFTrackCandPtr->getPosSeed();
  TVector3 mom = aGFTrackCandPtr->getDirSeed();

  TVector3 poserr = aGFTrackCandPtr->getPosError();
  TVector3 momerr = aGFTrackCandPtr->getDirError();

  calcStateCov(pos, mom, poserr, momerr);
}



void RKTrackRep::setData(const TMatrixD& st, const GFDetPlane& pl, const TMatrixD* cov, const TMatrixD* aux){
  if(aux != NULL) {
    fCacheSpu = (*aux)(0,0);
    fDirection = (*aux)(0,1);
  }
  else {
    if(pl!=fCachePlane){
      std::cerr << "RKTrackRep::setData() - a fatal error occurred! It was called with a reference plane which is not the same as the one from the last extrapolate(plane,state,cov)-> abort in line " << __LINE__ << std::endl;
      throw;
    }
  }
  GFAbsTrackRep::setData(st,pl,cov);
  if (fCharge*fState[0][0] < 0) fCharge *= -1; // set charge accordingly! (fState[0][0] = q/p)
  fSpu = fCacheSpu;
}


const TMatrixD* RKTrackRep::getAuxInfo(const GFDetPlane& pl) {

  if(pl!=fCachePlane) {
    std::cerr << "RKTrackRep::getAuxInfo() - Fatal error: Trying to get auxiliary information with planes mismatch (Information returned does not belong to requested plane)! -> abort in line " << __LINE__ << std::endl;
	  throw;
  }
  fAuxInfo.ResizeTo(1,2);
  fAuxInfo(0,0) = fCacheSpu;
  fAuxInfo(0,1) = fDirection;
  return &fAuxInfo;
}


void RKTrackRep::setPDG(int i){
  fPdg = i;
  TParticlePDG * part = TDatabasePDG::Instance()->GetParticle(fPdg);
  if(part == 0){
    GFException exc("RKTrackRep::setPDG ==> particle id not known to TDatabasePDG",__LINE__,__FILE__);
    throw exc;
  }
  fMass = part->Mass();
  fCharge = part->Charge()/(3.);
}



void RKTrackRep::calcStateCov(const TVector3& pos,
                              const TVector3& mom,
                              const TVector3& poserr,
                              const TVector3& momerr){

  calcState(pos, mom);

  double pw = mom.Mag();
  double pu = 0.;
  double pv = 0.;

  TVector3 u = fRefPlane.getU();
  TVector3 v = fRefPlane.getV();
  TVector3 w = fRefPlane.getNormal();


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

  fCov[1][1] = pow((u.X()/pw - w.X()*pu/(pw*pw)),2.) * momerr.X()*momerr.X() +
               pow((u.Y()/pw - w.Y()*pu/(pw*pw)),2.) * momerr.Y()*momerr.Y() +
               pow((u.Z()/pw - w.Z()*pu/(pw*pw)),2.) * momerr.Z()*momerr.Z();

  fCov[2][2] = pow((v.X()/pw - w.X()*pv/(pw*pw)),2.) * momerr.X()*momerr.X() +
               pow((v.Y()/pw - w.Y()*pv/(pw*pw)),2.) * momerr.Y()*momerr.Y() +
               pow((v.Z()/pw - w.Z()*pv/(pw*pw)),2.) * momerr.Z()*momerr.Z();

}


void RKTrackRep::calcState(const TVector3& pos,
                           const TVector3& mom){

  fRefPlane.setON(pos, mom);
  fSpu=1.;

  fState[0][0] = fCharge/mom.Mag();

  //u' and v'
  fState[1][0] = 0.;
  fState[2][0] = 0.;

  //u and v
  fState[3][0] = 0.;
  fState[4][0] = 0.;

}



TMatrixD RKTrackRep::getState7() const{
  return getState7(fState, fRefPlane, fSpu);
}


TMatrixD RKTrackRep::getState7(const TMatrixD& state5, const GFDetPlane& pl, const double& spu) const{

  const TVector3 U = pl.getU();
  const TVector3 V = pl.getV();

  const TVector3 point = pl.getO() + state5[3][0]*U + state5[4][0]*V;

  TVector3 dir = spu * (pl.getNormal() + state5[1][0]*U + state5[2][0]*V);
  dir.SetMag(1.);

  TMatrixD state7(7,1);

  state7[0][0] = point.X();
  state7[1][0] = point.Y();
  state7[2][0] = point.Z();
  state7[3][0] = dir.X();
  state7[4][0] = dir.Y();
  state7[5][0] = dir.Z();
  state7[6][0] = state5[0][0];

  return state7;
}


TMatrixD RKTrackRep::getState5(const TMatrixD& state7, const GFDetPlane& pl, double& spu) const {

  const TVector3 O = pl.getO();
  const TVector3 U = pl.getU();
  const TVector3 V = pl.getV();

  const TVector3 Point(state7[0][0], state7[1][0], state7[2][0]);
  TVector3 A(state7[3][0], state7[4][0], state7[5][0]);

  // force A to be in normal direction and set spu accordingly
  double AtW = A * pl.getNormal();
  spu = 1.;
  if (AtW < 0) {
    A *= -1.;
    AtW *= -1.;
    spu = -1.;
  }

  TMatrixD state5(5,1);
  state5[0][0] = state7[6][0];
  state5[1][0] = A*U / AtW;
  state5[2][0] = A*V / AtW;
  state5[3][0] = (Point-O)*U;
  state5[4][0] = (Point-O)*V;

  return state5;
}



void RKTrackRep::transformPM(const TMatrixD& in5x5, TMatrixD& out,
                             const GFDetPlane& pl, const TMatrixD& state5, const double&  spu, TMatrixD* Jac) const{

  // check if out is 6x6 or 7x7
  bool sixD(false);
  if (out.GetNcols() == 7 && out.GetNrows() == 7) {}
  else if (out.GetNcols() == 6 && out.GetNrows() == 6) sixD = true;
  else {
    GFException exc("RKTrackRep::transformPM ==> output matrix has to be 6x6 or 7x7",__LINE__,__FILE__);
    throw exc;
  }

  // get vectors and aux variables
  const TVector3 U = pl.getU();
  const TVector3 V = pl.getV();
  const TVector3 W = pl.getNormal();

  const TVector3 pTilde =  spu * (W + state5[1][0]*U + state5[2][0]*V);
  const double pTildeMag = pTilde.Mag();
  const double pTildeMag2 = pTildeMag*pTildeMag;

  const double utpTilde = U*pTilde;
  const double vtpTilde = V*pTilde;

  //J_pM matrix is d(x,y,z,ax,ay,az,q/p) / d(q/p,u',v',u,v)   (out is 7x7)
  //J_pM matrix is d(x,y,z,px,py,pz) / d(q/p,u',v',u,v)       (out is 6x6)
  TMatrixD J_pM(5,7);
  if (sixD) J_pM.ResizeTo(5,6);

   // d(x,y,z)/d(u)
  J_pM[3][0] = U.X();
  J_pM[3][1] = U.Y();
  J_pM[3][2] = U.Z();
  // d(x,y,z)/d(v)
  J_pM[4][0] = V.X();
  J_pM[4][1] = V.Y();
  J_pM[4][2] = V.Z();

  if (!sixD) { // 7D
    // d(q/p)/d(q/p)
    J_pM[0][6] = 1.;
    // d(ax,ay,az)/d(u')
    double fact = spu / pTildeMag;
    J_pM[1][3] = fact * ( U.X() - pTilde.X()*utpTilde/pTildeMag2 );
    J_pM[1][4] = fact * ( U.Y() - pTilde.Y()*utpTilde/pTildeMag2 );
    J_pM[1][5] = fact * ( U.Z() - pTilde.Z()*utpTilde/pTildeMag2 );
    // d(ax,ay,az)/d(v')
    J_pM[2][3] = fact * ( V.X() - pTilde.X()*vtpTilde/pTildeMag2 );
    J_pM[2][4] = fact * ( V.Y() - pTilde.Y()*vtpTilde/pTildeMag2 );
    J_pM[2][5] = fact * ( V.Z() - pTilde.Z()*vtpTilde/pTildeMag2 );
  }
  else { // 6D
    const double qop = state5[0][0];
    const double p = fCharge/qop; // momentum

    // d(px,py,pz)/d(q/p)
    double fact = -1. * p / (pTildeMag * qop);
    J_pM[0][3] = fact * pTilde.X();
    J_pM[0][4] = fact * pTilde.Y();
    J_pM[0][5] = fact * pTilde.Z();
    // d(px,py,pz)/d(u')
    fact = p * spu / pTildeMag;
    J_pM[1][3] = fact * ( U.X() - pTilde.X()*utpTilde/pTildeMag2 );
    J_pM[1][4] = fact * ( U.Y() - pTilde.Y()*utpTilde/pTildeMag2 );
    J_pM[1][5] = fact * ( U.Z() - pTilde.Z()*utpTilde/pTildeMag2 );
    // d(px,py,pz)/d(v')
    J_pM[2][3] = fact * ( V.X() - pTilde.X()*vtpTilde/pTildeMag2 );
    J_pM[2][4] = fact * ( V.Y() - pTilde.Y()*vtpTilde/pTildeMag2 );
    J_pM[2][5] = fact * ( V.Z() - pTilde.Z()*vtpTilde/pTildeMag2 );
  }


  TMatrixD J_pM_transp(J_pM);
  J_pM_transp.T();

  out = J_pM_transp*(in5x5*J_pM);

  if (Jac!=NULL){
    Jac->ResizeTo(J_pM);
    *Jac = J_pM;
  }

}


void RKTrackRep::transformMP(const TMatrixD& in, TMatrixD& out5x5,
                             const GFDetPlane& pl, const TMatrixD& state7, TMatrixD* Jac) const {

  // check if in is 6x6 or 7x7
  bool sixD(false);
  if (in.GetNcols() == 7 && in.GetNrows() == 7) {}
  else if (in.GetNcols() == 6 && in.GetNrows() == 6) sixD = true;
  else {
    GFException exc("RKTrackRep::transformMP ==> input matrix has to be 6x6 or 7x7",__LINE__,__FILE__);
    throw exc;
  }

  // get vectors and aux variables
  const TVector3 U = pl.getU();
  const TVector3 V = pl.getV();
  const TVector3 W = pl.getNormal();

  const TVector3 A(state7[3][0], state7[4][0], state7[5][0]);

  const double AtU = A*U;
  const double AtV = A*V;
  const double AtW = A*W;


  // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,ax,ay,az,q/p)   (in is 7x7)
  // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,px,py,pz)       (in is 6x6)
  TMatrixD J_Mp(7,5);
  if (sixD) J_Mp.ResizeTo(6,5);

  if (!sixD) { // 7D
    // d(u')/d(ax,ay,az)
    double fact = 1./(AtW*AtW);
    J_Mp[3][1] = fact * (U.X()*AtW - W.X()*AtU);
    J_Mp[4][1] = fact * (U.Y()*AtW - W.Y()*AtU);
    J_Mp[5][1] = fact * (U.Z()*AtW - W.Z()*AtU);
    // d(v')/d(ax,ay,az)
    J_Mp[3][2] = fact * (V.X()*AtW - W.X()*AtV);
    J_Mp[4][2] = fact * (V.Y()*AtW - W.Y()*AtV);
    J_Mp[5][2] = fact * (V.Z()*AtW - W.Z()*AtV);
    // d(q/p)/d(q/p)
    J_Mp[6][0] = 1.;
  }
  else { // 6D
    const double qop = state7[6][0];
    const double p = fCharge/qop; // momentum

    // d(q/p)/d(px,py,pz)
    double fact = (-1.) * qop / p;
    J_Mp[3][0] = fact * A.X();
    J_Mp[4][0] = fact * A.Y();
    J_Mp[5][0] = fact * A.Z();
    // d(u')/d(px,py,pz)
    fact = 1./(p*AtW*AtW);
    J_Mp[3][1] = fact * (U.X()*AtW - W.X()*AtU);
    J_Mp[4][1] = fact * (U.Y()*AtW - W.Y()*AtU);
    J_Mp[5][1] = fact * (U.Z()*AtW - W.Z()*AtU);
    // d(v')/d(px,py,pz)
    J_Mp[3][2] = fact * (V.X()*AtW - W.X()*AtV);
    J_Mp[4][2] = fact * (V.Y()*AtW - W.Y()*AtV);
    J_Mp[5][2] = fact * (V.Z()*AtW - W.Z()*AtV);
  }
  //d(u)/d(x,y,z)
  J_Mp[0][3] = U.X();
  J_Mp[1][3] = U.Y();
  J_Mp[2][3] = U.Z();
  //d(v)/d(x,y,z)
  J_Mp[0][4] = V.X();
  J_Mp[1][4] = V.Y();
  J_Mp[2][4] = V.Z();


  TMatrixD J_Mp_transp(J_Mp);
  J_Mp_transp.T();

  out5x5.ResizeTo(5, 5);

  out5x5 = J_Mp_transp*(in*J_Mp);

  if (Jac!=NULL){
    Jac->ResizeTo(J_Mp);
    *Jac = J_Mp;
  }
}



TVector3 RKTrackRep::getPos(const GFDetPlane& pl){
  TMatrixD state7(getState7());
  if(pl!=fRefPlane){
    Extrap(pl, &state7);
  }
  return TVector3(state7[0][0], state7[1][0], state7[2][0]);
}


TVector3 RKTrackRep::getMom(const GFDetPlane& pl){
  TMatrixD state7(getState7());
  if(pl!=fRefPlane){
    Extrap(pl, &state7);
  }

  TVector3 mom(state7[3][0], state7[4][0], state7[5][0]);
  mom.SetMag(fCharge/state7[6][0]);
  return mom;
}


void RKTrackRep::getPosMom(const GFDetPlane& pl,TVector3& pos, TVector3& mom){
  TMatrixD state7(getState7());
  if(pl!=fRefPlane){
    Extrap(pl, &state7);
  }

  pos.SetXYZ(state7[0][0], state7[1][0], state7[2][0]);
  mom.SetXYZ(state7[3][0], state7[4][0], state7[5][0]);
  mom.SetMag(fCharge/state7[6][0]);
}


void RKTrackRep::getPosMomCov(const GFDetPlane& pl, TVector3& pos, TVector3& mom, TMatrixD& cov6x6){

  TMatrixD statePred(fState);
  TMatrixD covPred(fCov);
  double spu(fSpu);

  if(pl != fRefPlane) {
    extrapolate(pl, statePred, covPred);
    spu = fCacheSpu;
  }

  TMatrixD state7(getState7(statePred, pl, spu));

  // cov
  cov6x6.ResizeTo(6, 6); // make sure cov has correct dimensions
  transformPM(covPred, cov6x6, pl, statePred, spu);

  pos.SetXYZ(state7[0][0], state7[1][0], state7[2][0]);
  mom.SetXYZ(state7[3][0], state7[4][0], state7[5][0]);
  mom.SetMag(fCharge/state7[6][0]);
}


void RKTrackRep::setPosMomCov(const TVector3& pos, const TVector3& mom, const TMatrixD& cov6x6){

  if (cov6x6.GetNcols()!=6 || cov6x6.GetNrows()!=6){
    GFException exc("RKTrackRep::setPosMomCov ==> cov has to be 6x6 (x, y, z, px, py, pz)",__LINE__,__FILE__);
    throw exc;
  }

  // fCharge does not change!
  calcState(pos, mom);

  fCachePlane = fRefPlane;
  fCacheSpu = 1.;

  TMatrixD state7(getState7());

  transformMP(cov6x6, fCov, fRefPlane, state7);
}



void RKTrackRep::extrapolateToPoint(const TVector3& pos, TVector3& poca, TVector3& dirInPoca){

  static const int maxIt(30);

  TMatrixD state7(getState7());

  double coveredDistance(0.);

  GFDetPlane pl;
  int iterations(0);

  while(true){
    pl.setON(pos, TVector3(state7[3][0], state7[4][0], state7[5][0]));
    coveredDistance =  this->Extrap(pl, &state7);

    if(fabs(coveredDistance) < MINSTEP) break;
    if(++iterations == maxIt) {
      GFException exc("RKTrackRep::extrapolateToPoint ==> extrapolation to point failed, maximum number of iterations reached",__LINE__,__FILE__);
      throw exc;
    }
  }
  poca.SetXYZ(state7[0][0], state7[1][0], state7[2][0]);
  dirInPoca.SetXYZ(state7[3][0], state7[4][0], state7[5][0]);
}


TVector3 RKTrackRep::poca2Line(const TVector3& extr1,const TVector3& extr2,const TVector3& point) const {
  
  TVector3 theWire = extr2-extr1;
  if(theWire.Mag()<1.E-8){
    GFException exc("RKTrackRep::poca2Line ==> try to find POCA between line and point, but the line is really just a point",__LINE__,__FILE__);
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

  TMatrixD state7(getState7());

  double coveredDistance(0.);

  GFDetPlane pl;
  int iterations(0);

  while(true){
    pl.setO(point1);
    TVector3 currentDir(state7[3][0], state7[4][0], state7[5][0]);
    pl.setU(currentDir.Cross(point2-point1));
    pl.setV(point2-point1);
    coveredDistance = this->Extrap(pl, &state7);

    if(fabs(coveredDistance) < MINSTEP) break;
    if(++iterations == maxIt) {
      GFException exc("RKTrackRep::extrapolateToLine ==> extrapolation to line failed, maximum number of iterations reached",__LINE__,__FILE__);
      throw exc;
    }
  }

  poca.SetXYZ(state7[0][0], state7[1][0], state7[2][0]);
  dirInPoca.SetXYZ(state7[3][0], state7[4][0], state7[5][0]);
  poca_onwire = poca2Line(point1,point2,poca);
}


double RKTrackRep::extrapolate(const GFDetPlane& pl, 
                               TMatrixD& statePred,
                               TMatrixD& covPred){
  
  TMatrixD state7(getState7());
  TMatrixD cov7x7(7,7);

  transformPM(fCov, cov7x7, fRefPlane, fState, fSpu);

  double coveredDistance = Extrap(pl, &state7, &cov7x7);
  
  statePred.ResizeTo(5,1);
  statePred = getState5(state7, pl, fCacheSpu);
  fCachePlane = pl;

  covPred.ResizeTo(5,5);
  transformMP(cov7x7, covPred, pl, state7);

  return coveredDistance;
}


double RKTrackRep::extrapolate(const GFDetPlane& pl, 
                               TMatrixD& statePred){

  TMatrixD state7(getState7());
  double coveredDistance = Extrap(pl, &state7);
  double spu;
  statePred.ResizeTo(5,1);
  statePred = getState5(state7, pl, spu);

  return coveredDistance;
}


double RKTrackRep::stepalong(double h, TVector3& pos, TVector3& dir){

  TVector3 dest;

  static const int maxIt(30);

  TMatrixD state7(getState7());

  double coveredDistance(0.);

  GFDetPlane pl;
  int iterations(0);

  while(true){
    pos.SetXYZ(state7[0][0], state7[1][0], state7[2][0]);
    dir.SetXYZ(state7[3][0], state7[4][0], state7[5][0]);
    dir.SetMag(1.);

    dest = pos + (h - coveredDistance) * dir;

    pl.setON(dest, dir);
    coveredDistance += this->Extrap(pl, &state7);

    if(fabs(h - coveredDistance)<MINSTEP) break;
    if(++iterations == maxIt) {
      GFException exc("RKTrackRep::stepalong ==> maximum number of iterations reached",__LINE__,__FILE__);
      throw exc;
    }
  }

  pos.SetXYZ(state7[0][0], state7[1][0], state7[2][0]);
  dir.SetXYZ(state7[3][0], state7[4][0], state7[5][0]);

  return coveredDistance;
}



double RKTrackRep::Extrap( const GFDetPlane& plane, TMatrixD* state, TMatrixD* cov) {

  //std::cerr<<"RKTrackRep::Extrap"<<std::endl;
  static const int maxNumIt(200);
  int numIt(0);

  bool calcCov(true);
  if(cov==NULL) calcCov=false;

  double *P;
  if(calcCov) {P = new double[56];}
  else {P = new double[7];}

  for(int i=0;i<7;++i){
    P[i] = (*state)[i][0];
  }

  TMatrixD jac(7,7);
  TMatrixD jacT(7,7);
  TMatrixD oldCov(7,7);
  if(calcCov) oldCov=(*cov);
  double coveredDistance(0.);
  double sumDistance(0.);

  while(true){
    if(numIt++ > maxNumIt){
      GFException exc("RKTrackRep::Extrap ==> maximum number of iterations exceeded",__LINE__,__FILE__);
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

    TVector3 directionBefore(P[3],P[4],P[5]); // direction before propagation
    directionBefore.SetMag(1.);

    // propagation
    std::vector<TVector3> points;
    std::vector<double> pointPaths;
    if( ! this->RKutta(plane, P, coveredDistance, points, pointPaths, calcCov) ) {
      GFException exc("RKTrackRep::Extrap ==> Runge Kutta propagation failed",__LINE__,__FILE__);
      exc.setFatal();
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

    TMatrixD noise(7,7); // zero everywhere by default

    // call MatEffects
    if (!fNoMaterial){
      double momLoss; // momLoss has a sign - negative loss means momentum gain

      momLoss = GFMaterialEffects::getInstance()->effects(pointsFilt,
                                 pointPathsFilt,
                                 fabs(fCharge/P[6]), // momentum
                                 fPdg,
                                 fXX0,
                                 calcCov,
                                 &noise,
                                 &jac,
                                 &directionBefore,
                                 &directionAfter);

      #ifdef DEBUG
        std::cout << "momLoss: " << momLoss << " GeV \n";
      #endif

      if(fabs(P[6])>1.E-10){ // do momLoss only for defined 1/momentum .ne.0
        P[6] = fCharge/(fabs(fCharge/P[6])-momLoss);
      }
    }

    if(calcCov){ //propagate cov and add noise
      TMatrixD absCov(oldCov);
      absCov.Abs();
      if(!(absCov < 1.E200)){
        GFException exc("RKTrackRep::Extrap ==> covariance matrix exceeds numerical limits",__LINE__,__FILE__);
        oldCov.Print();
        exc.setFatal();
        delete[] P;
        throw exc;
      }
      oldCov = *cov;
      *cov = jacT*((oldCov)*jac)+noise;
    }

    #ifdef DEBUG
      jacT.Print();
      if(calcCov) cov->Print();
    #endif

    //we arrived at the destination plane, if we point to the active area of the plane (if it is finite), and the distance is below threshold
    if( plane.inActive(TVector3(P[0],P[1],P[2]), TVector3(P[3],P[4],P[5])) &&
        plane.distance(P[0],P[1],P[2]) < MINSTEP) break;

  }

  (*state)[0][0] = P[0];  (*state)[1][0] = P[1];  (*state)[2][0] = P[2];
  (*state)[3][0] = P[3];  (*state)[4][0] = P[4];  (*state)[5][0] = P[5];
  (*state)[6][0] = P[6];

  delete[] P;
  return sumDistance;
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
                         bool calcCov) const {

  // important fixed numbers
  static const double EC     = .000149896229;   // c/(2*10^12) resp. c/2Tera
  static const double P3     = 1./3.;           // 1/3
  static const int    ND     = 56;              // number of variables for derivatives calculation
  static const int    ND1    = ND-7;            // = 49
  // limits, check-values, etc. Can be tuned!
  static const double DLT    = .0002;           // max. deviation for approximation-quality test
  static const double Wmax   = 3000.;           // max. way allowed [cm]
  static const double AngleMax = 6.3;           // max. total angle change of momentum. Prevents extrapolating a curler round and round if no active plane is found.
  static const double Pmin   = 4.E-3;           // minimum momentum for propagation [GeV]
  static const unsigned int maxNumIt = 1000;    // maximum number of iterations in main loop
  // Aux parameters
  double* R           = &P[0];                  // Start coordinates  [cm] 	(x,  y,  z)
  double* A           = &P[3];                  // Start directions 	      (ax, ay, az); 	ax^2+ay^2+az^2=1
  double  SA[3]       = {0.,0.,0.};             // Start directions derivatives dA/S
  double  Pinv        = P[6]*EC;                // P[6] is charge/momentum in e/(GeV/c)
  double  Way         = 0.;                     // Sum of absolute values of all extrapolation steps [cm]
  bool    stopBecauseOfMaterial = false;        // does not go through main loop again when stepsize is reduced by stepper
  TVector3 pos(R[0],R[1],R[2]);                 // position
  TVector3 dir(A[0],A[1],A[2]);                 // direction
  TVector3 Hvect(0.,0.,0);
  double   momentum   = fabs(fCharge/P[6]);     // momentum [GeV]
  double   relMomLoss = 0;                      // relative momentum loss in RKutta
  double   deltaAngle = 0.;                     // total angle by which the momentum has changed during extrapolation
  double   An(0), S(0), Sl(0);
  double   S3(0), S4(0), PS2(0), EST(0), CBA(0);
  // Variables for RKutta main loop
  double   SU[4]={0.,0.,0.,0.}, H0[3]={0.,0.,0.}, H1[3]={0.,0.,0.}, H2[3]={0.,0.,0.}, r[3]={0.,0.,0.};
  double   A0(0), A1(0), A2(0), A3(0), A4(0), A5(0), A6(0);
  double   B0(0), B1(0), B2(0), B3(0), B4(0), B5(0), B6(0);
  double   C0(0), C1(0), C2(0), C3(0), C4(0), C5(0), C6(0);
  double   dA0(0), dA2(0), dA3(0), dA4(0), dA5(0), dA6(0);
  double   dB0(0), dB2(0), dB3(0), dB4(0), dB5(0), dB6(0);
  double   dC0(0), dC2(0), dC3(0), dC4(0), dC5(0), dC6(0);

  #ifdef DEBUG
    std::cout << "RKTrackRep::RKutta \n";
    std::cout << "position: "; pos.Print();
    std::cout << "direction: "; dir.Print();
    std::cout << "momentum: " << momentum << " GeV\n";
    std::cout << "destination: "; plane.Print();
  #endif

  // check momentum
  if(momentum < Pmin){
    std::cerr << "RKTrackRep::RKutta ==> momentum too low: " << fabs(fCharge/P[6])*1000. << " MeV" << std::endl;
    return (false);
  }
  
  // clear points and store starting point
  coveredDistance = 0;
  points.clear();
  pointPaths.clear();
  points.push_back(pos);
  pointPaths.push_back(0.);


  // make SU vector point away from origin
  const TVector3 W = plane.getNormal();
  if(W*plane.getO() > 0){
    SU[0] = W.X();
    SU[1] = W.Y();
    SU[2] = W.Z();
  }
  else{
    SU[0] = -1.*W.X();
    SU[1] = -1.*W.Y();
    SU[2] = -1.*W.Z();
  }
  SU[3] = plane.distance(0., 0., 0.);


  // Step estimation (signed)
  S = estimateStep(pos, dir, SU, plane, momentum, relMomLoss, deltaAngle, stopBecauseOfMaterial);
  if (fabs(S) < MINSTEP) return true;

  unsigned int counter(0);

  //
  // Main cycle of Runge-Kutta method
  //
  while (fabs(S) >= MINSTEP) {

    if(++counter > maxNumIt){
      std::cerr << "RKTrackRep::RKutta ==> maximum number of iterations exceeded\n";
      return(false);
    }

    #ifdef DEBUG
      std::cout << " RKutta main loop nr. " << counter << "\n";
    #endif

    //
    // Runge Kutta Extrapolation
    //
    S3 = P3*S;
    S4 = 0.25*S;
    PS2 = Pinv*S;
    
    // First point
    r[0] = R[0];           r[1] = R[1];           r[2]=R[2];
    pos.SetXYZ(r[0], r[1], r[2]); // vector of start coordinates R0	(x, y, z)
    Hvect = GFFieldManager::getFieldVal(pos);				// magnetic field in 10^-4 T = kGauss
    H0[0] = PS2*Hvect.X(); H0[1] = PS2*Hvect.Y(); H0[2] = PS2*Hvect.Z(); 		// H0 is PS2*(Hx, Hy, Hz) @ R0
    A0 = A[1]*H0[2]-A[2]*H0[1]; B0 = A[2]*H0[0]-A[0]*H0[2]; C0 = A[0]*H0[1]-A[1]*H0[0]; // (ax, ay, az) x H0
    A2 = A[0]+A0              ; B2 = A[1]+B0              ; C2 = A[2]+C0              ; // (A0, B0, C0) + (ax, ay, az)
    A1 = A2+A[0]              ; B1 = B2+A[1]              ; C1 = C2+A[2]              ; // (A0, B0, C0) + 2*(ax, ay, az)
      
    // Second point
    r[0] += A1*S4;         r[1] += B1*S4;         r[2] += C1*S4;
    pos.SetXYZ(r[0], r[1], r[2]);
    Hvect = GFFieldManager::getFieldVal(pos);
    H1[0] = Hvect.X()*PS2; H1[1] = Hvect.Y()*PS2; H1[2] = Hvect.Z()*PS2;	// H1 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * [(A0, B0, C0) + 2*(ax, ay, az)]
    A3 = B2*H1[2]-C2*H1[1]+A[0]; B3 = C2*H1[0]-A2*H1[2]+A[1]; C3 = A2*H1[1]-B2*H1[0]+A[2]; // (A2, B2, C2) x H1 + (ax, ay, az)
    A4 = B3*H1[2]-C3*H1[1]+A[0]; B4 = C3*H1[0]-A3*H1[2]+A[1]; C4 = A3*H1[1]-B3*H1[0]+A[2]; // (A3, B3, C3) x H1 + (ax, ay, az)
    A5 = A4-A[0]+A4            ; B5 = B4-A[1]+B4            ; C5 = C4-A[2]+C4            ; //    2*(A4, B4, C4) - (ax, ay, az)

    // Last point
    r[0]=R[0]+S*A4;         r[1]=R[1]+S*B4;         r[2]=R[2]+S*C4;  //setup.Field(r,H2);
    pos.SetXYZ(r[0], r[1], r[2]);
    Hvect = GFFieldManager::getFieldVal(pos);
    H2[0] = Hvect.X()*PS2;  H2[1] = Hvect.Y()*PS2;  H2[2] = Hvect.Z()*PS2;	// H2 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * (A4, B4, C4)
    A6 = B5*H2[2]-C5*H2[1]; B6 = C5*H2[0]-A5*H2[2]; C6 = A5*H2[1]-B5*H2[0]; // (A5, B5, C5) x H2


    // Test approximation quality on given step and possible step reduction
    EST = fabs((A1+A6)-(A3+A4))+fabs((B1+B6)-(B3+B4))+fabs((C1+C6)-(C3+C4));  // EST = ||(ABC1+ABC6)-(ABC3+ABC4)||_1  =  ||(axzy x H0 + ABC5 x H2) - (ABC2 x H1 + ABC3 x H1)||_1
    if(EST > DLT) {
      S *= 0.5;
      stopBecauseOfMaterial = false;
      #ifdef DEBUG
        std::cout << " Stepsize was halfed to " << S << "\n";
      #endif
      continue;
    }
    
    // update paths
    coveredDistance += S;				// add stepsize to way (signed)
    Way  += fabs(S);

    // check way limit
    if(Way > Wmax){
      std::cerr << "RKTrackRep::RKutta ==> Total extrapolation length is longer than length limit : " << Way << " cm !" << std::endl;
      return(false);
    }


    //
    // Derivatives of track parameters in last point
    //
    if(calcCov){
      // d(x, y, z)/d(x, y, z) submatrix is unit matrix
      P[7] = 1;  P[15] = 1;  P[23] = 1;
      // d(ax, ay, az)/d(ax, ay, az) submatrix is 0
      // start with d(x, y, z)/d(ax, ay, az)
      for(int i=4*7; i!=ND; i+=7) {				// i = 7, 14, 21, 28, 35, 42, 49;    ND = 56;	ND1 = 49; rows of Jacobian
	
        double* dR = &P[i];			            		// dR = (dX/dpN,  dY/dpN,  dZ/dpN)
        double* dA = &P[i+3];				           	// dA = (dAx/dpN, dAy/dpN, dAz/dpN); N = X,Y,Z,Ax,Ay,Az,q/p
        
        //first point
        dA0 = H0[2]*dA[1]-H0[1]*dA[2];		// dA0/dp	}
        dB0 = H0[0]*dA[2]-H0[2]*dA[0];		// dB0/dp	 } = dA x H0
        dC0 = H0[1]*dA[0]-H0[0]*dA[1];		// dC0/dp	}
        
        if(i==ND1) {dA0+=A0; dB0+=B0; dC0+=C0;}			// if last row: (dA0, dB0, dC0) := (dA0, dB0, dC0) + (A0, B0, C0)
        
        dA2 = dA0+dA[0];				// }
        dB2 = dB0+dA[1]; 			  //  } = (dA0, dB0, dC0) + dA
        dC2 = dC0+dA[2];				// }
         
        //second point
        dA3 = dA[0]+dB2*H1[2]-dC2*H1[1];		// dA3/dp	}
        dB3 = dA[1]+dC2*H1[0]-dA2*H1[2];		// dB3/dp	 } = dA + (dA2, dB2, dC2) x H1
        dC3 = dA[2]+dA2*H1[1]-dB2*H1[0];		// dC3/dp	}
        
        if(i==ND1) {dA3+=A3-A[0]; dB3+=B3-A[1]; dC3+=C3-A[2];} // if last row: (dA3, dB3, dC3) := (dA3, dB3, dC3) + (A3, B3, C3) - (ax, ay, az)

        dA4 = dA[0]+dB3*H1[2]-dC3*H1[1];		// dA4/dp	}
        dB4 = dA[1]+dC3*H1[0]-dA3*H1[2];		// dB4/dp	 } = dA + (dA3, dB3, dC3) x H1
        dC4 = dA[2]+dA3*H1[1]-dB3*H1[0];		// dC4/dp	}
        
        if(i==ND1) {dA4+=A4-A[0]; dB4+=B4-A[1]; dC4+=C4-A[2];} // if last row: (dA4, dB4, dC4) := (dA4, dB4, dC4) + (A4, B4, C4) - (ax, ay, az)
        
        //last point	
        dA5 = dA4+dA4-dA[0];			// }
        dB5 = dB4+dB4-dA[1];	  	//  } =  2*(dA4, dB4, dC4) - dA
        dC5 = dC4+dC4-dA[2]; 			// }

        dA6 = dB5*H2[2]-dC5*H2[1];			// dA6/dp	}
        dB6 = dC5*H2[0]-dA5*H2[2];			// dB6/dp	 } = (dA5, dB5, dC5) x H2
        dC6 = dA5*H2[1]-dB5*H2[0];			// dC6/dp	}

        if(i==ND1) {dA6+=A6; dB6+=B6; dC6+=C6;}			// if last row: (dA6, dB6, dC6) := (dA6, dB6, dC6) + (A6, B6, C6)                                    
        
        dR[0] += (dA2+dA3+dA4)*S3;  dA[0] = (dA0+dA3+dA3+dA5+dA6)*P3;	// dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
        dR[1] += (dB2+dB3+dB4)*S3;  dA[1] = (dB0+dB3+dB3+dB5+dB6)*P3;	// dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
        dR[2] += (dC2+dC3+dC4)*S3;  dA[2] = (dC0+dC3+dC3+dC5+dC6)*P3;
      }
    }
    
    //
    // Track parameters in last point
    //   
    R[0] += (A2+A3+A4)*S3;   A[0] += (SA[0]=(A0+A3+A3+A5+A6)*P3-A[0]);  // R  = R0 + S3*[(A2, B2, C2) +   (A3, B3, C3) + (A4, B4, C4)]
    R[1] += (B2+B3+B4)*S3;   A[1] += (SA[1]=(B0+B3+B3+B5+B6)*P3-A[1]);  // A  =     1/3*[(A0, B0, C0) + 2*(A3, B3, C3) + (A5, B5, C5) + (A6, B6, C6)]
    R[2] += (C2+C3+C4)*S3;   A[2] += (SA[2]=(C0+C3+C3+C5+C6)*P3-A[2]); 	// SA = A_new - A_old
    pos.SetXYZ(R[0], R[1], R[2]);

    // normalize A
    CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);	// 1/|A|
    A[0] *= CBA; A[1] *= CBA; A[2] *= CBA;
    dir.SetXYZ(A[0], A[1], A[2]);

    // append point and steplength
    pointPaths.push_back(S);
    points.push_back(pos);
    
    // if stepsize has been limited by material, break the loop
    if (stopBecauseOfMaterial) break;

    // estimate Step for next loop or linear extrapolation
    Sl = S;	// last S used
    S = estimateStep(pos, dir, SU, plane, momentum, relMomLoss, deltaAngle, stopBecauseOfMaterial);

    // if this estimation has set stopBecauseOfMaterial to true and S<MINSTEP, then the loop should break and no linear extrapolation should be performed!
    if (stopBecauseOfMaterial && S < MINSTEP) break;

    // check if we went back and forth multiple times -> we don't come closer to the plane!
    if (counter > 15){
      if (S                    *pointPaths[counter]   < 0 &&
          pointPaths[counter  ]*pointPaths[counter-1] < 0 &&
          pointPaths[counter-1]*pointPaths[counter-2] < 0 &&
          pointPaths[counter-2]*pointPaths[counter-3] < 0){
        std::cerr << "RKTrackRep::RKutta ==> Do not get closer to plane!"<<std::endl;
        return(false);
      }
    }

    // check if total angle is bigger than AngleMax. Can happen if a curler should be fitted and it does not hit the active area of the next plane.
    if (fabs(deltaAngle) > AngleMax){
      std::cerr << "RKTrackRep::RKutta ==> Do not get to an active plane! Already extrapolated " << deltaAngle * 180 / TMath::Pi() << "°." <<std::endl;
      return(false);
    }

  } //end of main loop
  

  //
  // linear extrapolation to surface
  //
  if (!stopBecauseOfMaterial) {
    bool bugfix (true);
    if( bugfix ) {
      if (fabs(Sl) > 1.E-12) Sl = 1./Sl;        // Sl = inverted last Stepsize Sl
      else{ 
	GFException exc("RKTrackRep::Extrap ==> last stepsize too small -> can't do linear extrapolation!",__LINE__,__FILE__);
	throw exc;
      }
    }
    // normalize SA
    SA[0]*=Sl;  SA[1]*=Sl;  SA[2]*=Sl;
    // calculate A
    A[0] += SA[0]*S;   	// S  = distance to surface
    A[1] += SA[1]*S;   	// SA*Sl = delta A / delta way; local derivative of A with respect to the length of the way
    A[2] += SA[2]*S;	  // A = A + S * SA*Sl

    // normalize A
    CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);  // 1/|A|
    A[0] *= CBA; A[1] *= CBA; A[2] *= CBA;

    R[0] = R[0]+S*(A[0]-0.5*S*SA[0]);    // P = R + S*(A - 1/2*S*SA); approximation for final point on surface
    R[1] = R[1]+S*(A[1]-0.5*S*SA[1]);
    R[2] = R[2]+S*(A[2]-0.5*S*SA[2]);

    //
    // Output derivatives of track parameters preparation
    //
    if(calcCov){
      An = A[0]*SU[0]+A[1]*SU[1]+A[2]*SU[2];
      fabs(An) > 1.E-7 ? An=1./An : An = 0; // 1/A_normal
      double norm;
      for(int i=7; i!=ND; i+=7) {
        double* dR = &P[i];
        double* dA = &P[i+3];
        norm = (dR[0]*SU[0] + dR[1]*SU[1] + dR[2]*SU[2])*An;	// dR_normal / A_normal
        dR[0] -= norm*A [0];   dR[1] -= norm*A [1];   dR[2] -= norm*A [2];
        dA[0] -= norm*SA[0];   dA[1] -= norm*SA[1];   dA[2] -= norm*SA[2];
      }
    }

    // update points and paths
    points.push_back(TVector3(R[0], R[1], R[2]));
    pointPaths.push_back(S);

    coveredDistance += S;
    Way  += fabs(S);
  }

  return(true);
}


double RKTrackRep::estimateStep(const TVector3& pos,
                                const TVector3& dir,
                                const double* SU,
                                const GFDetPlane& plane,
                                const double& momentum,
                                double& relMomLoss,
                                double& deltaAngle,
                                bool& stopBecauseOfMaterial) const{

  static const double Smax      = 100.;          // max. step allowed [cm]
  static const double dAngleMax = 0.1;           // max. deviation of angle between direction before and after the step [rad]
  double Step;

  #ifdef DEBUG
    std::cout << " RKTrackRep::estimateStep \n";
    std::cout << "  position: "; pos.Print();
    std::cout << "  direction: "; dir.Print();
  #endif


  // calculate distance to surface
  double Dist = SU[3] - (pos[0]*SU[0]+pos[1]*SU[1]+pos[2]*SU[2]);        // Distance between start coordinates and surface
  double An = dir[0]*SU[0] + dir[1]*SU[1] + dir[2]*SU[2];    // An = dir * N;  component of dir normal to surface

  if (fabs(An) > 1.E-10) Step = Dist/An;
  else {
    Step = Dist*1.E10;
    if (An<0) Step *= -1.;
  }

  // see if dir points towards surface (1) or not (-1)
  int dirSw(1);
  if (Step<0) dirSw = -1;

  Step = fabs(Step);

  #ifdef DEBUG
    std::cout << "  Distance to plane: " << Dist << "\n";
    std::cout << "  guess for Step (unsigned): " << Step << "\n";
    if (dirSw>0) std::cout << "  Direction is  pointing towards surface.\n";
    else  std::cout << "  Direction is pointing away from surface.\n";
  #endif

  // calculate way after which momentum angle has changed AngleMax
  TVector3 Hvect(GFFieldManager::getFieldVal(pos));       // magnetic field in 10^-4 T = kGauss
  double Hmag(Hvect.Mag());
  double SmaxAngle(Smax);
  double radius(0);
  if (Hmag > 1E-5){
    double p_perp = ( dir - Hvect*((dir*Hvect)/(Hmag*Hmag)) ).Mag() * momentum; // [GeV]
    radius = p_perp/(0.3E-3*Hmag); // [cm]
    double sinAngle = fabs(sin(dir.Angle(Hvect)));
    if (sinAngle > 1E-10){
      SmaxAngle = fabs(dAngleMax * radius / sinAngle); // [cm]
    }
  }


  //
  // Select direction
  //
  // auto select
  if (fDirection == 0){
    Step *= dirSw;
    #ifdef DEBUG
      std::cout << "  auto select direction. \n";
    #endif
  }
  // see if straight line approximation is ok
  else if ( Step < 0.2*SmaxAngle ){
    #ifdef DEBUG
      std::cout << "  straight line approximation is fine. Delta angle until surface is reached is approx " << Step/SmaxAngle * dAngleMax * 180 / TMath::Pi()  << " deg \n";
    #endif

    // if direction is pointing to active part of surface
    if( plane.inActive(pos, dir) ) {
      Step *= dirSw;
      #ifdef DEBUG
        std::cout << "  direction is pointing to active part of surface. \n";
      #endif
    }
    // if we are near the plane, but not pointing to the active area, make a big step!
    else {
      Step = fDirection*SmaxAngle;
      #ifdef DEBUG
        std::cout << "  we are near the plane, but not pointing to the active area. make a big step! \n";
      #endif
    }
  }
  // fDirection decides!
  else {
    Step *= fDirection;
    #ifdef DEBUG
      std::cout << "  select direction according to fDirection. \n";
    #endif
  }

  #ifdef DEBUG
    std::cout << "  guess for Step (signed): " << Step << "\n";
  #endif


  //
  // Limit stepsize
  //
  double Ssign = 1.;
  if (Step < 0) Ssign = -1;

  // reduce maximum stepsize Step to Smax
  if (Step > Smax) Step = Smax;
  else if (Step < -Smax) Step = -Smax;

  // also limit stepsize according to the change of the momentum direction!
  if (Step > SmaxAngle) Step = SmaxAngle;
  else if (Step < -SmaxAngle) Step = -SmaxAngle;

  #ifdef DEBUG
    std::cout << "  limit from maxangle: " << SmaxAngle << ", radius: " << radius << "\n";
  #endif

  // call stepper and reduce stepsize if step not too small
  if (fabs(Step) > MINSTEP){
    double stepperLen = GFMaterialEffects::getInstance()->stepper(fabs(Step),
                                                                  pos[0], pos[1], pos[2],
                                                                  Ssign*dir[0], Ssign*dir[1], Ssign*dir[2],
                                                                  momentum,
                                                                  relMomLoss,
                                                                  fPdg);

    if (Step > stepperLen) {
      Step = stepperLen;
      stopBecauseOfMaterial = true;
    }
    else if (Step < -stepperLen) {
      Step = -stepperLen;
      stopBecauseOfMaterial = true;
    }

    #ifdef DEBUG
      std::cout << "  limit from stepper: " << stepperLen << "\n";
    #endif
  }

  #ifdef DEBUG
    std::cout << "  --> Step used: " << Step << "\n";
  #endif

  deltaAngle += Step/SmaxAngle * dAngleMax;

  return Step;
}


void RKTrackRep::setPropDir(int dir){
  // make sure fDirection is -1, 0 or 1
  if (dir>0) fDirection = 1;
  else if (dir<0) fDirection = -1;
  else fDirection = 0;
}



ClassImp(RKTrackRep)
