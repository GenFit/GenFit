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
  initArrays();
}


RKTrackRep::RKTrackRep() : GFAbsTrackRep(5), fDirection(0), fNoMaterial(false), fPdg(0), fMass(0.), fCharge(-1), fCachePlane(), fCacheSpu(1), fSpu(1), fAuxInfo(1,2) {
  initArrays();
}


RKTrackRep::RKTrackRep(const TVector3& pos,
                       const TVector3& mom,
                       const TVector3& poserr,
                       const TVector3& momerr,
                       const int& PDGCode) :
                       GFAbsTrackRep(5), fDirection(0), fNoMaterial(false), fCachePlane(), fCacheSpu(1), fAuxInfo(1,2) {

  initArrays();
  setPDG(PDGCode); // also sets charge and mass
  calcStateCov(pos, mom, poserr, momerr);
}


RKTrackRep::RKTrackRep(const TVector3& pos,
                       const TVector3& mom,
                       const int& PDGCode) :
                       GFAbsTrackRep(5), fDirection(0), fNoMaterial(false), fCachePlane(), fCacheSpu(1), fAuxInfo(1,2) {

  initArrays();
  setPDG(PDGCode); // also sets charge and mass
  calcState(pos, mom);

  // set covariance diagonal elements to large number
  static const double value(1.E4);

  fCov(0,0) = value;
  fCov(1,1) = value;
  fCov(2,2) = value;
  fCov(3,3) = value;
  fCov(4,4) = value;
}


RKTrackRep::RKTrackRep(const GFTrackCand* aGFTrackCandPtr) :
                       GFAbsTrackRep(5), fDirection(0), fNoMaterial(false), fCachePlane(), fCacheSpu(1), fAuxInfo(1,2) {

  initArrays();
  setPDG(aGFTrackCandPtr->getPdgCode()); // also sets charge and mass

  double mom = aGFTrackCandPtr->getQoverPseed();
  if (fabs(mom) > 1.E-3) mom = fCharge/mom;
  else mom = 1.E3;

  calcStateCov(aGFTrackCandPtr->getPosSeed(),
               aGFTrackCandPtr->getDirSeed()*mom,
               aGFTrackCandPtr->getPosError(),
               aGFTrackCandPtr->getDirError()*mom);
}


void RKTrackRep::initArrays(){
  memset(fStateJac,0x00,(7+7*7)*sizeof(double));
  memset(fNoise,0x00,7*7*sizeof(double));
  memset(fStateJac,0x00,7*7*sizeof(double));
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
  if (fCharge*fState(0,0) < 0) fCharge *= -1; // set charge accordingly! (fState[0][0] = q/p)
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


  fCov(0,0) = fCharge*fCharge/pow(mom.Mag(),6.) *
              (mom.X()*mom.X() * momerr.X()*momerr.X()+
               mom.Y()*mom.Y() * momerr.Y()*momerr.Y()+
               mom.Z()*mom.Z() * momerr.Z()*momerr.Z());

  fCov(1,1) = pow((u.X()/pw - w.X()*pu/(pw*pw)),2.) * momerr.X()*momerr.X() +
              pow((u.Y()/pw - w.Y()*pu/(pw*pw)),2.) * momerr.Y()*momerr.Y() +
              pow((u.Z()/pw - w.Z()*pu/(pw*pw)),2.) * momerr.Z()*momerr.Z();

  fCov(2,2) = pow((v.X()/pw - w.X()*pv/(pw*pw)),2.) * momerr.X()*momerr.X() +
              pow((v.Y()/pw - w.Y()*pv/(pw*pw)),2.) * momerr.Y()*momerr.Y() +
              pow((v.Z()/pw - w.Z()*pv/(pw*pw)),2.) * momerr.Z()*momerr.Z();

  fCov(3,3) = poserr.X()*poserr.X() * u.X()*u.X() +
              poserr.Y()*poserr.Y() * u.Y()*u.Y() +
              poserr.Z()*poserr.Z() * u.Z()*u.Z();

  fCov(4,4) = poserr.X()*poserr.X() * v.X()*v.X() +
              poserr.Y()*poserr.Y() * v.Y()*v.Y() +
              poserr.Z()*poserr.Z() * v.Z()*v.Z();
}


void RKTrackRep::calcState(const TVector3& pos,
                           const TVector3& mom){

  fRefPlane.setON(pos, mom);
  fSpu=1.;

  fState(0,0) = fCharge/mom.Mag();

  //u' and v'
  fState(1,0) = 0.;
  fState(2,0) = 0.;

  //u and v
  fState(3,0) = 0.;
  fState(4,0) = 0.;
}



void RKTrackRep::getState7(M1x7& state7) const{
  getState7(state7, fState, fRefPlane, fSpu);
}


void RKTrackRep::getState7(M1x7& state7, const TMatrixD& state5, const GFDetPlane& pl, const double& spu) const{

  const TVector3 U = pl.getU();
  const TVector3 V = pl.getV();

  const TVector3 point = pl.getO() + state5(3,0)*U + state5(4,0)*V;

  TVector3 dir = spu * (pl.getNormal() + state5(1,0)*U + state5(2,0)*V);
  dir.SetMag(1.);

  for(unsigned int i=0; i<3; ++i){
    state7[i] = point(i); // (x, y, z
    state7[3+i] = dir(i); // a_x, a_y, a_z
  }
  state7[6] = state5(0,0);
}


TMatrixD RKTrackRep::getState5(const M1x7& state7, const GFDetPlane& pl, double& spu) const {

  const TVector3 O = pl.getO();
  const TVector3 U = pl.getU();
  const TVector3 V = pl.getV();

  const TVector3 Point(state7[0], state7[1], state7[2]);
  TVector3 A(state7[3], state7[4], state7[5]);

  // force A to be in normal direction and set spu accordingly
  double AtW = A * pl.getNormal();
  spu = 1.;
  if (AtW < 0) {
    A *= -1.;
    AtW *= -1.;
    spu = -1.;
  }

  TMatrixD state5(5,1);
  state5(0,0) = state7[6];
  state5(1,0) = A*U / AtW;
  state5(2,0) = A*V / AtW;
  state5(3,0) = (Point-O)*U;
  state5(4,0) = (Point-O)*V;

  return state5;
}



void RKTrackRep::transformPM7(const TMatrixD& in5x5, M7x7& out7x7,
                              const GFDetPlane& pl, const TMatrixD& state5, const double&  spu,
                              TMatrixD* Jac) const{

  // get vectors and aux variables
  const TVector3 U = pl.getU();
  const TVector3 V = pl.getV();

  const TVector3 pTilde =  spu * (pl.getNormal() + state5(1,0)*U + state5(2,0)*V);
  const double pTildeMag = pTilde.Mag();
  const double pTildeMag2 = pTildeMag*pTildeMag;

  const double utpTildeOverpTildeMag2 = U*pTilde / pTildeMag2;
  const double vtpTildeOverpTildeMag2 = V*pTilde / pTildeMag2;

  //J_pM matrix is d(x,y,z,ax,ay,az,q/p) / d(q/p,u',v',u,v)   (out is 7x7)
  M5x7 J_pM_;
  memset(J_pM_,0x00,5*7*sizeof(double));

   // d(x,y,z)/d(u)
  J_pM_[21] = U.X(); // [3][0]
  J_pM_[22] = U.Y(); // [3][1]
  J_pM_[23] = U.Z(); // [3][2]
  // d(x,y,z)/d(v)
  J_pM_[28] = V.X(); // [4][2]
  J_pM_[29] = V.Y(); // [4][2]
  J_pM_[30] = V.Z(); // [4][2]
  // d(q/p)/d(q/p)
  J_pM_[6] = 1.; // not needed for array matrix multiplication
  // d(ax,ay,az)/d(u')
  double fact = spu / pTildeMag;
  J_pM_[10] = fact * ( U.X() - pTilde.X()*utpTildeOverpTildeMag2 ); // [1][3]
  J_pM_[11] = fact * ( U.Y() - pTilde.Y()*utpTildeOverpTildeMag2 ); // [1][4]
  J_pM_[12] = fact * ( U.Z() - pTilde.Z()*utpTildeOverpTildeMag2 ); // [1][5]
  // d(ax,ay,az)/d(v')
  J_pM_[17] = fact * ( V.X() - pTilde.X()*vtpTildeOverpTildeMag2 ); // [2][3]
  J_pM_[18] = fact * ( V.Y() - pTilde.Y()*vtpTildeOverpTildeMag2 ); // [2][4]
  J_pM_[19] = fact * ( V.Z() - pTilde.Z()*vtpTildeOverpTildeMag2 ); // [2][5]


  // since the Jacobian contains a lot of zeros, and the resulting cov has to be symmetric,
  // the multiplication can be done much faster directly on array level
  // out = J_pM^T * in5x5 * J_pM
  const M5x5& in5x5_ = *((M5x5*) in5x5.GetMatrixArray());

  double JTC0  = J_pM_[21] * in5x5_[18] + J_pM_[28] * in5x5_[23];
  double JTC1  = J_pM_[21] * in5x5_[23] + J_pM_[28] * in5x5_[24];
  double JTC2  = J_pM_[21] * in5x5_[16] + J_pM_[28] * in5x5_[21];
  double JTC3  = J_pM_[21] * in5x5_[17] + J_pM_[28] * in5x5_[22];
  double JTC4  = J_pM_[22] * in5x5_[18] + J_pM_[29] * in5x5_[23];
  double JTC5  = J_pM_[22] * in5x5_[23] + J_pM_[29] * in5x5_[24];
  double JTC6  = J_pM_[22] * in5x5_[16] + J_pM_[29] * in5x5_[21];
  double JTC7  = J_pM_[22] * in5x5_[17] + J_pM_[29] * in5x5_[22];
  double JTC8  = J_pM_[23] * in5x5_[18] + J_pM_[30] * in5x5_[23];
  double JTC9  = J_pM_[23] * in5x5_[23] + J_pM_[30] * in5x5_[24];
  double JTC10 = J_pM_[23] * in5x5_[16] + J_pM_[30] * in5x5_[21];
  double JTC11 = J_pM_[23] * in5x5_[17] + J_pM_[30] * in5x5_[22];
  double JTC12 = J_pM_[10] * in5x5_[6]  + J_pM_[17] * in5x5_[11];
  double JTC13 = J_pM_[10] * in5x5_[11] + J_pM_[17] * in5x5_[12];
  double JTC14 = J_pM_[11] * in5x5_[6]  + J_pM_[18] * in5x5_[11];
  double JTC15 = J_pM_[11] * in5x5_[11] + J_pM_[18] * in5x5_[12];

  // loops are vectorizable by the compiler!
  for (int i=0; i<3; ++i) out7x7[i]  = JTC0 * J_pM_[21+i] + JTC1 * J_pM_[28+i];
  for (int i=0; i<3; ++i) out7x7[3+i]  = JTC2 * J_pM_[10+i] + JTC3 * J_pM_[17+i];
  out7x7[6]  =  J_pM_[21] * in5x5_[15] + J_pM_[28] * in5x5_[20];

  for (int i=0; i<2; ++i) out7x7[8+i]  = JTC4 * J_pM_[22+i] + JTC5 * J_pM_[29+i];
  for (int i=0; i<3; ++i) out7x7[10+i] = JTC6 * J_pM_[10+i] + JTC7 * J_pM_[17+i];
  out7x7[13] = J_pM_[22] * in5x5_[15] + J_pM_[29] * in5x5_[20];

  out7x7[16] = JTC8 * J_pM_[23] + JTC9 * J_pM_[30];
  for (int i=0; i<3; ++i) out7x7[17+i] = JTC10 * J_pM_[10+i] + JTC11 * J_pM_[17+i];
  out7x7[20] = J_pM_[23] * in5x5_[15] + J_pM_[30] * in5x5_[20];

  for (int i=0; i<3; ++i) out7x7[24+i] = JTC12 * J_pM_[10+i] + JTC13 * J_pM_[17+i];
  out7x7[27] = J_pM_[10] * in5x5_[5] + J_pM_[17] * in5x5_[10];

  for (int i=0; i<2; ++i) out7x7[32+i] = JTC14 * J_pM_[11+i] + JTC15 * J_pM_[18+i];
  out7x7[34] = J_pM_[11] * in5x5_[5] + J_pM_[18] * in5x5_[10];

  out7x7[40] = (J_pM_[12] * in5x5_[6] + J_pM_[19] * in5x5_[11]) * J_pM_[12] + (J_pM_[12] * in5x5_[11] + J_pM_[19] * in5x5_[12]) * J_pM_[19];
  out7x7[41] = J_pM_[12] * in5x5_[5] + J_pM_[19] * in5x5_[10];

  out7x7[48] = in5x5_[0];

  // symmetric part
  out7x7[7]  = out7x7[1];
  out7x7[14] = out7x7[2];  out7x7[15] = out7x7[9];
  out7x7[21] = out7x7[3];  out7x7[22] = out7x7[10];  out7x7[23] = out7x7[17];
  out7x7[28] = out7x7[4];  out7x7[29] = out7x7[11];  out7x7[30] = out7x7[18];  out7x7[31] = out7x7[25];
  out7x7[35] = out7x7[5];  out7x7[36] = out7x7[12];  out7x7[37] = out7x7[19];  out7x7[38] = out7x7[26];  out7x7[39] = out7x7[33];
  out7x7[42] = out7x7[6];  out7x7[43] = out7x7[13];  out7x7[44] = out7x7[20];  out7x7[45] = out7x7[27];  out7x7[46] = out7x7[34];  out7x7[47] = out7x7[41];

  if (Jac!=NULL){
    Jac->ResizeTo(5,7);
    *Jac = TMatrixD(5,7, &(J_pM_[0]));
  }
}


void RKTrackRep::transformPM6(const TMatrixD& in5x5, M6x6& out6x6,
                              const GFDetPlane& pl, const TMatrixD& state5, const double&  spu,
                              TMatrixD* Jac) const{

  // get vectors and aux variables
  const TVector3 U = pl.getU();
  const TVector3 V = pl.getV();

  const TVector3 pTilde =  spu * (pl.getNormal() + state5(1,0)*U + state5(2,0)*V);
  const double pTildeMag = pTilde.Mag();
  const double pTildeMag2 = pTildeMag*pTildeMag;

  const double utpTildeOverpTildeMag2 = U*pTilde / pTildeMag2;
  const double vtpTildeOverpTildeMag2 = V*pTilde / pTildeMag2;

  //J_pM matrix is d(x,y,z,px,py,pz) / d(q/p,u',v',u,v)       (out is 6x6)
  double J_pM_[5*6];
  memset(J_pM_,0x00,5*6*sizeof(double));

  const double qop = state5(0,0);
  const double p = fCharge/qop; // momentum

  // d(px,py,pz)/d(q/p)
  double fact = -1. * p / (pTildeMag * qop);
  J_pM_[3] = fact * pTilde.X(); // [0][3]
  J_pM_[4] = fact * pTilde.Y(); // [0][4]
  J_pM_[5] = fact * pTilde.Z(); // [0][5]
  // d(px,py,pz)/d(u')
  fact = p * spu / pTildeMag;
  J_pM_[9]  = fact * ( U.X() - pTilde.X()*utpTildeOverpTildeMag2 ); // [1][3]
  J_pM_[10] = fact * ( U.Y() - pTilde.Y()*utpTildeOverpTildeMag2 ); // [1][4]
  J_pM_[11] = fact * ( U.Z() - pTilde.Z()*utpTildeOverpTildeMag2 ); // [1][5]
  // d(px,py,pz)/d(v')
  J_pM_[15] = fact * ( V.X() - pTilde.X()*vtpTildeOverpTildeMag2 ); // [2][3]
  J_pM_[16] = fact * ( V.Y() - pTilde.Y()*vtpTildeOverpTildeMag2 ); // [2][4]
  J_pM_[17] = fact * ( V.Z() - pTilde.Z()*vtpTildeOverpTildeMag2 ); // [2][5]
  // d(x,y,z)/d(u)
  J_pM_[18] = U.X(); // [3][0]
  J_pM_[19] = U.Y(); // [3][1]
  J_pM_[20] = U.Z(); // [3][2]
  // d(x,y,z)/d(v)
  J_pM_[24] = V.X(); // [4][0]
  J_pM_[25] = V.Y(); // [4][1]
  J_pM_[26] = V.Z(); // [4][2]


  // do the transformation
  // out = J_pM^T * in5x5 * J_pM
  const M5x5& in5x5_ = *((M5x5*) in5x5.GetMatrixArray());

  double JTC0  = (J_pM_[3*6+0] * in5x5_[3*5+3] + J_pM_[4*6+0] * in5x5_[4*5+3]);
  double JTC1  = (J_pM_[3*6+0] * in5x5_[4*5+3] + J_pM_[4*6+0] * in5x5_[4*5+4]);
  double JTC2  = (J_pM_[3*6+0] * in5x5_[3*5+0] + J_pM_[4*6+0] * in5x5_[4*5+0]);
  double JTC3  = (J_pM_[3*6+0] * in5x5_[3*5+1] + J_pM_[4*6+0] * in5x5_[4*5+1]);
  double JTC4  = (J_pM_[3*6+0] * in5x5_[3*5+2] + J_pM_[4*6+0] * in5x5_[4*5+2]);
  double JTC5  = (J_pM_[3*6+1] * in5x5_[3*5+3] + J_pM_[4*6+1] * in5x5_[4*5+3]);
  double JTC6  = (J_pM_[3*6+1] * in5x5_[4*5+3] + J_pM_[4*6+1] * in5x5_[4*5+4]);
  double JTC7  = (J_pM_[3*6+1] * in5x5_[3*5+0] + J_pM_[4*6+1] * in5x5_[4*5+0]);
  double JTC8  = (J_pM_[3*6+1] * in5x5_[3*5+1] + J_pM_[4*6+1] * in5x5_[4*5+1]);
  double JTC9  = (J_pM_[3*6+1] * in5x5_[3*5+2] + J_pM_[4*6+1] * in5x5_[4*5+2]);
  double JTC10 = (J_pM_[3*6+2] * in5x5_[3*5+0] + J_pM_[4*6+2] * in5x5_[4*5+0]);
  double JTC11 = (J_pM_[3*6+2] * in5x5_[3*5+1] + J_pM_[4*6+2] * in5x5_[4*5+1]);
  double JTC12 = (J_pM_[3*6+2] * in5x5_[3*5+2] + J_pM_[4*6+2] * in5x5_[4*5+2]);
  double JTC13 = (J_pM_[0*6+3] * in5x5_[0*5+0] + J_pM_[1*6+3] * in5x5_[1*5+0] + J_pM_[2*6+3] * in5x5_[2*5+0]);
  double JTC14 = (J_pM_[0*6+3] * in5x5_[1*5+0] + J_pM_[1*6+3] * in5x5_[1*5+1] + J_pM_[2*6+3] * in5x5_[2*5+1]);
  double JTC15 = (J_pM_[0*6+3] * in5x5_[2*5+0] + J_pM_[1*6+3] * in5x5_[2*5+1] + J_pM_[2*6+3] * in5x5_[2*5+2]);
  double JTC16 = (J_pM_[0*6+4] * in5x5_[0*5+0] + J_pM_[1*6+4] * in5x5_[1*5+0] + J_pM_[2*6+4] * in5x5_[2*5+0]);
  double JTC17 = (J_pM_[0*6+4] * in5x5_[1*5+0] + J_pM_[1*6+4] * in5x5_[1*5+1] + J_pM_[2*6+4] * in5x5_[2*5+1]);
  double JTC18 = (J_pM_[0*6+4] * in5x5_[2*5+0] + J_pM_[1*6+4] * in5x5_[2*5+1] + J_pM_[2*6+4] * in5x5_[2*5+2]);

  // loops are vectorizable by the compiler!
  for (int i=0; i<3; ++i) out6x6[i] = JTC0 * J_pM_[18+i] + JTC1 * J_pM_[24+i];
  for (int i=0; i<3; ++i) out6x6[3+i] = JTC2 * J_pM_[3+i] + JTC3 * J_pM_[9+i] + JTC4 * J_pM_[15+i];

  for (int i=0; i<2; ++i) out6x6[7+i] = JTC5 * J_pM_[19+i] + JTC6 * J_pM_[25+i];
  for (int i=0; i<3; ++i) out6x6[9+i] = JTC7 * J_pM_[3+i] + JTC8 * J_pM_[9+i] + JTC9 * J_pM_[15+i];

  out6x6[12+2] = (J_pM_[18+2] * in5x5_[3*5+3] + J_pM_[24+2] * in5x5_[4*5+3]) * J_pM_[18+2] + (J_pM_[18+2] * in5x5_[4*5+3] + J_pM_[24+2] * in5x5_[4*5+4]) * J_pM_[24+2];
  for (int i=0; i<3; ++i) out6x6[15+i] = JTC10 * J_pM_[3+i] + JTC11 * J_pM_[9+i] + JTC12 * J_pM_[15+i];

  for (int i=0; i<3; ++i) out6x6[21+i] = JTC13 * J_pM_[3+i] + JTC14 * J_pM_[9+i] + JTC15 * J_pM_[15+i];

  for (int i=0; i<3; ++i) out6x6[28+i] = JTC16 * J_pM_[4+i] + JTC17 * J_pM_[10+i] + JTC18 * J_pM_[16+i];

  out6x6[30+5] = (J_pM_[5] * in5x5_[0*5+0] + J_pM_[6+5] * in5x5_[1*5+0] + J_pM_[12+5] * in5x5_[2*5+0]) * J_pM_[5] + (J_pM_[5] * in5x5_[1*5+0] + J_pM_[6+5] * in5x5_[1*5+1] + J_pM_[12+5] * in5x5_[2*5+1]) * J_pM_[6+5] + (J_pM_[5] * in5x5_[2*5+0] + J_pM_[6+5] * in5x5_[2*5+1] + J_pM_[12+5] * in5x5_[2*5+2]) * J_pM_[12+5];

  // symmetric part
  out6x6[6+0] = out6x6[1];
  out6x6[12+0] = out6x6[2];  out6x6[12+1] = out6x6[6+2];
  out6x6[18+0] = out6x6[3];  out6x6[18+1] = out6x6[6+3];  out6x6[18+2] = out6x6[12+3];
  out6x6[24+0] = out6x6[4];  out6x6[24+1] = out6x6[6+4];  out6x6[24+2] = out6x6[12+4];  out6x6[24+3] = out6x6[18+4];
  out6x6[30+0] = out6x6[5];  out6x6[30+1] = out6x6[6+5];  out6x6[30+2] = out6x6[12+5];  out6x6[30+3] = out6x6[18+5];  out6x6[30+4] = out6x6[24+5];

  if (Jac!=NULL){
    Jac->ResizeTo(5,6);
    *Jac = TMatrixD(5,6, &(J_pM_[0]));
  }
}


void RKTrackRep::transformM7P(const M7x7& in7x7, TMatrixD& out5x5,
                              const GFDetPlane& pl, const M1x7& state7,
                              TMatrixD* Jac) const {

  out5x5.ResizeTo(5, 5);

  // get vectors and aux variables
  const TVector3 U = pl.getU();
  const TVector3 V = pl.getV();
  const TVector3 W = pl.getNormal();

  const TVector3 A(state7[3], state7[4], state7[5]);

  const double AtU = A*U;
  const double AtV = A*V;
  const double AtW = A*W;

  // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,ax,ay,az,q/p)   (in is 7x7)
  double J_Mp_[7*5];
  memset(J_Mp_,0x00,7*5*sizeof(double));

  // d(u')/d(ax,ay,az)
  double fact = 1./(AtW*AtW);
  J_Mp_[16] = fact * (U.X()*AtW - W.X()*AtU); // [3][1]
  J_Mp_[21] = fact * (U.Y()*AtW - W.Y()*AtU); // [4][1]
  J_Mp_[26] = fact * (U.Z()*AtW - W.Z()*AtU); // [5][1]
  // d(v')/d(ax,ay,az)
  J_Mp_[17] = fact * (V.X()*AtW - W.X()*AtV); // [3][2]
  J_Mp_[22] = fact * (V.Y()*AtW - W.Y()*AtV); // [4][2]
  J_Mp_[27] = fact * (V.Z()*AtW - W.Z()*AtV); // [5][2]
  // d(q/p)/d(q/p)
  J_Mp_[30] = 1.; // [6][0]  - not needed for array matrix multiplication
  //d(u)/d(x,y,z)
  J_Mp_[3] = U.X(); // [0][3]
  J_Mp_[8] = U.Y(); // [1][3]
  J_Mp_[13] = U.Z(); // [2][3]
  //d(v)/d(x,y,z)
  J_Mp_[4] = V.X(); // [0][4]
  J_Mp_[9] = V.Y(); // [1][4]
  J_Mp_[14] = V.Z(); // [2][4]


  // since the Jacobian contains a lot of zeros, and the resulting cov has to be symmetric,
  // the multiplication can be done much faster directly on array level
  // out5x5 = J_Mp^T * in * J_Mp
  M5x5& out5x5_ = *((M5x5*) out5x5.GetMatrixArray());

  double JTC0  = (J_Mp_[16] * in7x7[24] + J_Mp_[21] * in7x7[31] + J_Mp_[26] * in7x7[38]);
  double JTC1  = (J_Mp_[16] * in7x7[31] + J_Mp_[21] * in7x7[32] + J_Mp_[26] * in7x7[39]);
  double JTC2  = (J_Mp_[16] * in7x7[38] + J_Mp_[21] * in7x7[39] + J_Mp_[26] * in7x7[40]);
  double JTC3  = (J_Mp_[16] * in7x7[21] + J_Mp_[21] * in7x7[28] + J_Mp_[26] * in7x7[35]);
  double JTC4  = (J_Mp_[16] * in7x7[22] + J_Mp_[21] * in7x7[29] + J_Mp_[26] * in7x7[36]);
  double JTC5  = (J_Mp_[16] * in7x7[23] + J_Mp_[21] * in7x7[30] + J_Mp_[26] * in7x7[37]);
  double JTC6  = (J_Mp_[17] * in7x7[21] + J_Mp_[22] * in7x7[28] + J_Mp_[27] * in7x7[35]);
  double JTC7  = (J_Mp_[17] * in7x7[22] + J_Mp_[22] * in7x7[29] + J_Mp_[27] * in7x7[36]);
  double JTC8  = (J_Mp_[17] * in7x7[23] + J_Mp_[22] * in7x7[30] + J_Mp_[27] * in7x7[37]);
  double JTC9  = (J_Mp_[3] * in7x7[0] + J_Mp_[8] * in7x7[7] + J_Mp_[13] * in7x7[14]);
  double JTC10 = (J_Mp_[3] * in7x7[7] + J_Mp_[8] * in7x7[8] + J_Mp_[13] * in7x7[15]);
  double JTC11 = (J_Mp_[3] * in7x7[14] + J_Mp_[8] * in7x7[15] + J_Mp_[13] * in7x7[16]);

  out5x5_[0] = in7x7[48];
  out5x5_[1] = J_Mp_[16] * in7x7[45] + J_Mp_[21] * in7x7[46] + J_Mp_[26] * in7x7[47];
  out5x5_[2] = J_Mp_[17] * in7x7[45] + J_Mp_[22] * in7x7[46] + J_Mp_[27] * in7x7[47];
  out5x5_[3] = J_Mp_[3] * in7x7[42] + J_Mp_[8] * in7x7[43] + J_Mp_[13] * in7x7[44];
  out5x5_[4] = J_Mp_[4] * in7x7[42] + J_Mp_[9] * in7x7[43] + J_Mp_[14] * in7x7[44];

  // loops are vectorizable by the compiler!
  for (int i=0; i<2; ++i) out5x5_[6+i] = JTC0 * J_Mp_[16+i] + JTC1 * J_Mp_[21+i] + JTC2 * J_Mp_[26+i];
  for (int i=0; i<2; ++i) out5x5_[8+i] = JTC3 * J_Mp_[3+i] + JTC4 * J_Mp_[8+i] + JTC5 * J_Mp_[13+i];

  out5x5_[12] = (J_Mp_[17] * in7x7[24] + J_Mp_[22] * in7x7[31] + J_Mp_[27] * in7x7[38]) * J_Mp_[17] + (J_Mp_[17] * in7x7[31] + J_Mp_[22] * in7x7[32] + J_Mp_[27] * in7x7[39]) * J_Mp_[22] + (J_Mp_[17] * in7x7[38] + J_Mp_[22] * in7x7[39] + J_Mp_[27] * in7x7[40]) * J_Mp_[27];
  for (int i=0; i<2; ++i) out5x5_[13+i] = JTC6 * J_Mp_[3+i] + JTC7 * J_Mp_[8+i] + JTC8 * J_Mp_[13+i];

  for (int i=0; i<2; ++i) out5x5_[18+i] = JTC9 * J_Mp_[3+i] + JTC10 * J_Mp_[8+i] + JTC11 * J_Mp_[13+i];

  out5x5_[24] = (J_Mp_[4] * in7x7[0] + J_Mp_[9] * in7x7[7] + J_Mp_[14] * in7x7[14]) * J_Mp_[4] + (J_Mp_[4] * in7x7[7] + J_Mp_[9] * in7x7[8] + J_Mp_[14] * in7x7[15]) * J_Mp_[9] + (J_Mp_[4] * in7x7[14] + J_Mp_[9] * in7x7[15] + J_Mp_[14] * in7x7[16]) * J_Mp_[14];

  // symmetric part
  out5x5_[5] = out5x5_[1];
  out5x5_[10] = out5x5_[2];  out5x5_[11] = out5x5_[7];
  out5x5_[15] = out5x5_[3];  out5x5_[16] = out5x5_[8];  out5x5_[17] = out5x5_[13];
  out5x5_[20] = out5x5_[4];  out5x5_[21] = out5x5_[9];  out5x5_[22] = out5x5_[14];  out5x5_[23] = out5x5_[19];

  if (Jac!=NULL){
    Jac->ResizeTo(7,5);
    *Jac = TMatrixD(7,5, &(J_Mp_[0]));
  }
}


void RKTrackRep::transformM6P(const M6x6& in6x6, TMatrixD& out5x5,
                              const GFDetPlane& pl, const M1x7& state7,
                              TMatrixD* Jac) const {

  out5x5.ResizeTo(5, 5);

  // get vectors and aux variables
  const TVector3 U = pl.getU();
  const TVector3 V = pl.getV();
  const TVector3 W = pl.getNormal();

  const TVector3 A(state7[3], state7[4], state7[5]);

  const double AtU = A*U;
  const double AtV = A*V;
  const double AtW = A*W;

  // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,px,py,pz)       (in is 6x6)
  double J_Mp_[6*5];
  memset(J_Mp_,0x00,6*5*sizeof(double));

  const double qop = state7[6];
  const double p = fCharge/qop; // momentum

  //d(u)/d(x,y,z)
  J_Mp_[3] = U.X();  // [0][3]
  J_Mp_[8] = U.Y();  // [1][3]
  J_Mp_[13] = U.Z(); // [2][3]
  //d(v)/d(x,y,z)
  J_Mp_[4] = V.X();  // [0][4]
  J_Mp_[9] = V.Y();  // [1][4]
  J_Mp_[14] = V.Z(); // [2][4]
  // d(q/p)/d(px,py,pz)
  double fact = (-1.) * qop / p;
  J_Mp_[15] = fact * A.X(); // [3][0]
  J_Mp_[20] = fact * A.Y(); // [4][0]
  J_Mp_[25] = fact * A.Z(); // [5][0]
  // d(u')/d(px,py,pz)
  fact = 1./(p*AtW*AtW);
  J_Mp_[16] = fact * (U.X()*AtW - W.X()*AtU); // [3][1]
  J_Mp_[21] = fact * (U.Y()*AtW - W.Y()*AtU); // [4][1]
  J_Mp_[26] = fact * (U.Z()*AtW - W.Z()*AtU); // [5][1]
  // d(v')/d(px,py,pz)
  J_Mp_[17] = fact * (V.X()*AtW - W.X()*AtV); // [3][2]
  J_Mp_[22] = fact * (V.Y()*AtW - W.Y()*AtV); // [4][2]
  J_Mp_[27] = fact * (V.Z()*AtW - W.Z()*AtV); // [5][2]

  // do the transformation
  // out5x5 = J_Mp^T * in * J_Mp
  M5x5& out5x5_ = *((M5x5*) out5x5.GetMatrixArray());

  double JTC0  = (J_Mp_[3*5+0] * in6x6[3*6+3] + J_Mp_[4*5+0] * in6x6[4*6+3] + J_Mp_[5*5+0] * in6x6[5*6+3]);
  double JTC1  = (J_Mp_[3*5+0] * in6x6[4*6+3] + J_Mp_[4*5+0] * in6x6[4*6+4] + J_Mp_[5*5+0] * in6x6[5*6+4]);
  double JTC2  = (J_Mp_[3*5+0] * in6x6[5*6+3] + J_Mp_[4*5+0] * in6x6[5*6+4] + J_Mp_[5*5+0] * in6x6[5*6+5]);
  double JTC3  = (J_Mp_[3*5+0] * in6x6[3*6+0] + J_Mp_[4*5+0] * in6x6[4*6+0] + J_Mp_[5*5+0] * in6x6[5*6+0]);
  double JTC4  = (J_Mp_[3*5+0] * in6x6[3*6+1] + J_Mp_[4*5+0] * in6x6[4*6+1] + J_Mp_[5*5+0] * in6x6[5*6+1]);
  double JTC5  = (J_Mp_[3*5+0] * in6x6[3*6+2] + J_Mp_[4*5+0] * in6x6[4*6+2] + J_Mp_[5*5+0] * in6x6[5*6+2]);
  double JTC6  = (J_Mp_[3*5+1] * in6x6[3*6+3] + J_Mp_[4*5+1] * in6x6[4*6+3] + J_Mp_[5*5+1] * in6x6[5*6+3]);
  double JTC7  = (J_Mp_[3*5+1] * in6x6[4*6+3] + J_Mp_[4*5+1] * in6x6[4*6+4] + J_Mp_[5*5+1] * in6x6[5*6+4]);
  double JTC8  = (J_Mp_[3*5+1] * in6x6[5*6+3] + J_Mp_[4*5+1] * in6x6[5*6+4] + J_Mp_[5*5+1] * in6x6[5*6+5]);
  double JTC9  = (J_Mp_[3*5+1] * in6x6[3*6+0] + J_Mp_[4*5+1] * in6x6[4*6+0] + J_Mp_[5*5+1] * in6x6[5*6+0]);
  double JTC10 = (J_Mp_[3*5+1] * in6x6[3*6+1] + J_Mp_[4*5+1] * in6x6[4*6+1] + J_Mp_[5*5+1] * in6x6[5*6+1]);
  double JTC11 = (J_Mp_[3*5+1] * in6x6[3*6+2] + J_Mp_[4*5+1] * in6x6[4*6+2] + J_Mp_[5*5+1] * in6x6[5*6+2]);
  double JTC12 = (J_Mp_[3*5+2] * in6x6[3*6+0] + J_Mp_[4*5+2] * in6x6[4*6+0] + J_Mp_[5*5+2] * in6x6[5*6+0]);
  double JTC13 = (J_Mp_[3*5+2] * in6x6[3*6+1] + J_Mp_[4*5+2] * in6x6[4*6+1] + J_Mp_[5*5+2] * in6x6[5*6+1]);
  double JTC14 = (J_Mp_[3*5+2] * in6x6[3*6+2] + J_Mp_[4*5+2] * in6x6[4*6+2] + J_Mp_[5*5+2] * in6x6[5*6+2]);
  double JTC15 = (J_Mp_[0*5+3] * in6x6[0*6+0] + J_Mp_[1*5+3] * in6x6[1*6+0] + J_Mp_[2*5+3] * in6x6[2*6+0]);
  double JTC16 = (J_Mp_[0*5+3] * in6x6[1*6+0] + J_Mp_[1*5+3] * in6x6[1*6+1] + J_Mp_[2*5+3] * in6x6[2*6+1]);
  double JTC17 = (J_Mp_[0*5+3] * in6x6[2*6+0] + J_Mp_[1*5+3] * in6x6[2*6+1] + J_Mp_[2*5+3] * in6x6[2*6+2]);

  // loops are vectorizable by the compiler!
  for (int i=0; i<3; ++i) out5x5_[i] = JTC0 * J_Mp_[15+i] + JTC1 * J_Mp_[20+i] + JTC2 * J_Mp_[25+i];
  for (int i=0; i<2; ++i) out5x5_[3+i] = JTC3 * J_Mp_[3+i] + JTC4 * J_Mp_[8+i] + JTC5 * J_Mp_[13+i];

  for (int i=0; i<2; ++i) out5x5_[6+i] = JTC6 * J_Mp_[16+i] + JTC7 * J_Mp_[21+i] + JTC8 * J_Mp_[26+i];
  for (int i=0; i<2; ++i) out5x5_[8+i] = JTC9 * J_Mp_[3+i] + JTC10 * J_Mp_[8+i] + JTC11 * J_Mp_[13+i];

  out5x5_[10+2] = (J_Mp_[15+2] * in6x6[3*6+3] + J_Mp_[20+2] * in6x6[4*6+3] + J_Mp_[5*5+2] * in6x6[5*6+3]) * J_Mp_[15+2] + (J_Mp_[15+2] * in6x6[4*6+3] + J_Mp_[20+2] * in6x6[4*6+4] + J_Mp_[5*5+2] * in6x6[5*6+4]) * J_Mp_[20+2] + (J_Mp_[15+2] * in6x6[5*6+3] + J_Mp_[20+2] * in6x6[5*6+4] + J_Mp_[5*5+2] * in6x6[5*6+5]) * J_Mp_[5*5+2];
  for (int i=0; i<2; ++i) out5x5_[13+i] = JTC12 * J_Mp_[3+i] + JTC13 * J_Mp_[8+i] + JTC14 * J_Mp_[13+i];

  for (int i=0; i<2; ++i) out5x5_[18+i] = JTC15 * J_Mp_[3+i] + JTC16 * J_Mp_[8+i] + JTC17 * J_Mp_[13+i];

  out5x5_[20+4] = (J_Mp_[4] * in6x6[0*6+0] + J_Mp_[5+4] * in6x6[1*6+0] + J_Mp_[10+4] * in6x6[2*6+0]) * J_Mp_[4] + (J_Mp_[4] * in6x6[1*6+0] + J_Mp_[5+4] * in6x6[1*6+1] + J_Mp_[10+4] * in6x6[2*6+1]) * J_Mp_[5+4] + (J_Mp_[4] * in6x6[2*6+0] + J_Mp_[5+4] * in6x6[2*6+1] + J_Mp_[10+4] * in6x6[2*6+2]) * J_Mp_[10+4];

  // symmetric part
  out5x5_[5+0] = out5x5_[1];
  out5x5_[10+0] = out5x5_[2];  out5x5_[10+1] = out5x5_[5+2];
  out5x5_[15+0] = out5x5_[3];  out5x5_[15+1] = out5x5_[5+3];  out5x5_[15+2] = out5x5_[10+3];
  out5x5_[20+0] = out5x5_[4];  out5x5_[20+1] = out5x5_[5+4];  out5x5_[20+2] = out5x5_[10+4];  out5x5_[20+3] = out5x5_[15+4];

  if (Jac!=NULL){
    Jac->ResizeTo(6,5);
    *Jac = TMatrixD(6,5, &(J_Mp_[0]));;
  }
}



TVector3 RKTrackRep::getPos(const GFDetPlane& pl){
  M1x7 state7;
  getState7(state7);
  if(pl!=fRefPlane) Extrap(pl, state7);
  return TVector3(state7[0], state7[1], state7[2]);
}


TVector3 RKTrackRep::getMom(const GFDetPlane& pl){
  M1x7 state7;
  getState7(state7);
  if(pl!=fRefPlane) Extrap(pl, state7);

  TVector3 mom(state7[3], state7[4], state7[5]);
  mom.SetMag(fCharge/state7[6]);
  return mom;
}


void RKTrackRep::getPosMom(const GFDetPlane& pl,TVector3& pos, TVector3& mom){
  M1x7 state7;
  getState7(state7);
  if(pl!=fRefPlane) Extrap(pl, state7);

  pos.SetXYZ(state7[0], state7[1], state7[2]);
  mom.SetXYZ(state7[3], state7[4], state7[5]);
  mom.SetMag(fCharge/state7[6]);
}


void RKTrackRep::getPosMomCov(const GFDetPlane& pl, TVector3& pos, TVector3& mom, TMatrixD& cov6x6){
  TMatrixD statePred(fState);
  TMatrixD covPred(fCov);
  double spu(fSpu);

  if(pl != fRefPlane) {
    extrapolate(pl, statePred, covPred);
    spu = fCacheSpu;
  }

  M1x7 state7;
  getState7(state7, statePred, pl, spu);

  // cov
  cov6x6.ResizeTo(6, 6); // make sure cov has correct dimensions
  M6x6& cov6x6_ = *((M6x6*) cov6x6.GetMatrixArray());
  transformPM6(covPred, cov6x6_, pl, statePred, spu);

  pos.SetXYZ(state7[0], state7[1], state7[2]);
  mom.SetXYZ(state7[3], state7[4], state7[5]);
  mom.SetMag(fCharge/state7[6]);
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

  M1x7 state7;
  getState7(state7);

  const M6x6& cov6x6_ = *((M6x6*) cov6x6.GetMatrixArray());

  transformM6P(cov6x6_, fCov, fRefPlane, state7);
}



void RKTrackRep::extrapolateToPoint(const TVector3& pos, TVector3& poca, TVector3& dirInPoca){

  static const unsigned int maxIt(30);

  M1x7 state7;
  getState7(state7);

  double coveredDistance(0.);

  GFDetPlane pl;
  unsigned int iterations(0);
  TVector3 currentDir;

  while(true){
    currentDir.SetXYZ(state7[3], state7[4], state7[5]);
    pl.setON(pos, currentDir);
    coveredDistance =  this->Extrap(pl, state7);

    if(fabs(coveredDistance) < MINSTEP) break;
    if(++iterations == maxIt) {
      GFException exc("RKTrackRep::extrapolateToPoint ==> extrapolation to point failed, maximum number of iterations reached",__LINE__,__FILE__);
      throw exc;
    }
  }
  poca.SetXYZ(state7[0], state7[1], state7[2]);
  dirInPoca.SetXYZ(state7[3], state7[4], state7[5]);
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
  static const unsigned int maxIt(30);

  M1x7 state7;
  getState7(state7);

  double coveredDistance(0.);

  GFDetPlane pl;
  unsigned int iterations(0);
  TVector3 currentDir;

  while(true){
    pl.setO(point1);
    currentDir.SetXYZ(state7[3], state7[4], state7[5]);
    pl.setU(currentDir.Cross(point2-point1));
    pl.setV(point2-point1);
    coveredDistance = this->Extrap(pl, state7);

    if(fabs(coveredDistance) < MINSTEP) break;
    if(++iterations == maxIt) {
      GFException exc("RKTrackRep::extrapolateToLine ==> extrapolation to line failed, maximum number of iterations reached",__LINE__,__FILE__);
      throw exc;
    }
  }

  poca.SetXYZ(state7[0], state7[1], state7[2]);
  dirInPoca.SetXYZ(state7[3], state7[4], state7[5]);
  poca_onwire = poca2Line(point1,point2,poca);
}


double RKTrackRep::extrapolate(const GFDetPlane& pl, 
                               TMatrixD& statePred,
                               TMatrixD& covPred){
  
  M1x7 state7;
  getState7(state7);
  M7x7 cov7x7;

  transformPM7(fCov, cov7x7, fRefPlane, fState, fSpu);

  double coveredDistance = Extrap(pl, state7, &cov7x7);
  
  statePred.ResizeTo(5,1);
  statePred = getState5(state7, pl, fCacheSpu);
  fCachePlane = pl;

  covPred.ResizeTo(5,5);
  transformM7P(cov7x7, covPred, pl, state7);

  return coveredDistance;
}


double RKTrackRep::extrapolate(const GFDetPlane& pl, 
                               TMatrixD& statePred){

  M1x7 state7;
  getState7(state7);
  double coveredDistance = Extrap(pl, state7);
  double spu;
  statePred.ResizeTo(5,1);
  statePred = getState5(state7, pl, spu);

  return coveredDistance;
}


double RKTrackRep::stepalong(double h, TVector3& pos, TVector3& dir){

  TVector3 dest;

  static const unsigned int maxIt(30);
  double coveredDistance(0.);

  M1x7 state7;
  getState7(state7);

  GFDetPlane pl;
  unsigned int iterations(0);

  while(true){
    pos.SetXYZ(state7[0], state7[1], state7[2]);
    dir.SetXYZ(state7[3], state7[4], state7[5]);
    dir.SetMag(1.);

    dest = pos + (h - coveredDistance) * dir;

    pl.setON(dest, dir);
    coveredDistance += this->Extrap(pl, state7);

    if(fabs(h - coveredDistance)<MINSTEP) break;
    if(++iterations == maxIt) {
      GFException exc("RKTrackRep::stepalong ==> maximum number of iterations reached",__LINE__,__FILE__);
      throw exc;
    }
  }

  pos.SetXYZ(state7[0], state7[1], state7[2]);
  dir.SetXYZ(state7[3], state7[4], state7[5]);

  return coveredDistance;
}



double RKTrackRep::Extrap( const GFDetPlane& plane, M1x7& state7, M7x7* cov) {

  static const unsigned int maxNumIt(200);
  unsigned int numIt(0);

  const bool calcCov(cov!=NULL);

  // set initial state
  memcpy(fStateJac, state7, 7*sizeof(double));

  double coveredDistance(0.);
  double sumDistance(0.);

  while(true){
    if(numIt++ > maxNumIt){
      GFException exc("RKTrackRep::Extrap ==> maximum number of iterations exceeded",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    // initialize fStateJac with unit matrix; last entry is not 1 but q/p
    if(calcCov){
      memset(&fStateJac[7],0x00,49*sizeof(double));
      for(int i=0; i<6; ++i){
        fStateJac[(i+1)*7+i] = 1.;
      }
      fStateJac[55] =  state7[6];
    }

    TVector3 directionBefore(fStateJac[3], fStateJac[4], fStateJac[5]); // direction before propagation

    // propagation
    std::vector<TVector3> points;
    points.reserve(10);
    std::vector<double> pointPaths;
    pointPaths.reserve(10);
    if( ! this->RKutta(plane, fStateJac, coveredDistance, points, pointPaths, calcCov) ) {
      GFException exc("RKTrackRep::Extrap ==> Runge Kutta propagation failed",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    TVector3 directionAfter(fStateJac[3], fStateJac[4], fStateJac[5]); // direction after propagation

    sumDistance+=coveredDistance;

    // filter Points
    unsigned int nPoints = points.size();
    std::vector<TVector3> pointsFilt(1, points.at(0));
    pointsFilt.reserve(nPoints);
    std::vector<double> pointPathsFilt(1, 0.);
    pointPathsFilt.reserve(nPoints);
    // only if in right direction
    for(unsigned int i=1; i<nPoints; ++i){
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

    if(calcCov){
      // calculate Jacobian -> divide last row of fStateJac by q/p
      // -> last column of jac is [0,0,0,0,0,0,1]
      for (unsigned int i=0; i<7; ++i) fStateJac[49+i] /= fStateJac[6];

      // set fNoise to 0
      memset(fNoise,0x00,7*7*sizeof(double));
    }


    // call MatFX
    if (!fNoMaterial){
      // momLoss has a sign - negative loss means momentum gain
      double momLoss = GFMaterialEffects::getInstance()->effects(pointsFilt,
                                 pointPathsFilt,
                                 fabs(fCharge/fStateJac[6]), // momentum
                                 fPdg,
                                 fXX0,
                                 calcCov,
                                 fNoise,
                                 &(fStateJac[7]),
                                 &directionBefore,
                                 &directionAfter);

      #ifdef DEBUG
        std::cout << "momLoss: " << momLoss << " GeV \n";
      #endif

      if(fabs(fStateJac[6])>1.E-10){ // do momLoss only for defined 1/momentum .ne.0
        fStateJac[6] = fCharge/(fabs(fCharge/fStateJac[6])-momLoss);
      }
    } // finished MatFX

    if(calcCov){ // propagate cov and add noise
      memcpy(fOldCov, (*cov), 7*7*sizeof(double));
      M7x7& cov_ = (*cov);

      // numerical check:
      for(unsigned int i=0; i<7*7; ++i){
        if(fabs((*cov)[i]) > 1.E100){
          GFException exc("RKTrackRep::Extrap ==> (*cov)ariance matrix exceeds numerical limits",__LINE__,__FILE__);
          exc.setFatal();
          throw exc;
        }
      }

      // cov = Jac^T * oldCov * Jac;
      // last column of jac is [0,0,0,0,0,0,1]
      // cov is symmetric
      double JTC0  = fStateJac[7] * fOldCov[0] + fStateJac[7+7] * fOldCov[7] + fStateJac[7+2*7] * fOldCov[2*7] + fStateJac[7+3*7] * fOldCov[3*7] + fStateJac[7+4*7] * fOldCov[4*7] + fStateJac[7+5*7] * fOldCov[5*7] + fStateJac[7+6*7] * fOldCov[6*7];
      double JTC1  = fStateJac[7] * fOldCov[7] + fStateJac[7+7] * fOldCov[7+1] + fStateJac[7+2*7] * fOldCov[2*7+1] + fStateJac[7+3*7] * fOldCov[3*7+1] + fStateJac[7+4*7] * fOldCov[4*7+1] + fStateJac[7+5*7] * fOldCov[5*7+1] + fStateJac[7+6*7] * fOldCov[6*7+1];
      double JTC2  = fStateJac[7] * fOldCov[2*7] + fStateJac[7+7] * fOldCov[2*7+1] + fStateJac[7+2*7] * fOldCov[2*7+2] + fStateJac[7+3*7] * fOldCov[3*7+2] + fStateJac[7+4*7] * fOldCov[4*7+2] + fStateJac[7+5*7] * fOldCov[5*7+2] + fStateJac[7+6*7] * fOldCov[6*7+2];
      double JTC3  = fStateJac[7] * fOldCov[3*7] + fStateJac[7+7] * fOldCov[3*7+1] + fStateJac[7+2*7] * fOldCov[3*7+2] + fStateJac[7+3*7] * fOldCov[3*7+3] + fStateJac[7+4*7] * fOldCov[4*7+3] + fStateJac[7+5*7] * fOldCov[5*7+3] + fStateJac[7+6*7] * fOldCov[6*7+3];
      double JTC4  = fStateJac[7] * fOldCov[4*7] + fStateJac[7+7] * fOldCov[4*7+1] + fStateJac[7+2*7] * fOldCov[4*7+2] + fStateJac[7+3*7] * fOldCov[4*7+3] + fStateJac[7+4*7] * fOldCov[4*7+4] + fStateJac[7+5*7] * fOldCov[5*7+4] + fStateJac[7+6*7] * fOldCov[6*7+4];
      double JTC5  = fStateJac[7] * fOldCov[5*7] + fStateJac[7+7] * fOldCov[5*7+1] + fStateJac[7+2*7] * fOldCov[5*7+2] + fStateJac[7+3*7] * fOldCov[5*7+3] + fStateJac[7+4*7] * fOldCov[5*7+4] + fStateJac[7+5*7] * fOldCov[5*7+5] + fStateJac[7+6*7] * fOldCov[6*7+5];
      double JTC6  = fStateJac[7] * fOldCov[6*7] + fStateJac[7+7] * fOldCov[6*7+1] + fStateJac[7+2*7] * fOldCov[6*7+2] + fStateJac[7+3*7] * fOldCov[6*7+3] + fStateJac[7+4*7] * fOldCov[6*7+4] + fStateJac[7+5*7] * fOldCov[6*7+5] + fStateJac[7+6*7] * fOldCov[6*7+6];

      double JTC7  = fStateJac[7+1] * fOldCov[0] + fStateJac[7+7+1] * fOldCov[7] + fStateJac[7+2*7+1] * fOldCov[2*7] + fStateJac[7+3*7+1] * fOldCov[3*7] + fStateJac[7+4*7+1] * fOldCov[4*7] + fStateJac[7+5*7+1] * fOldCov[5*7] + fStateJac[7+6*7+1] * fOldCov[6*7];
      double JTC8  = fStateJac[7+1] * fOldCov[7] + fStateJac[7+7+1] * fOldCov[7+1] + fStateJac[7+2*7+1] * fOldCov[2*7+1] + fStateJac[7+3*7+1] * fOldCov[3*7+1] + fStateJac[7+4*7+1] * fOldCov[4*7+1] + fStateJac[7+5*7+1] * fOldCov[5*7+1] + fStateJac[7+6*7+1] * fOldCov[6*7+1];
      double JTC9  = fStateJac[7+1] * fOldCov[2*7] + fStateJac[7+7+1] * fOldCov[2*7+1] + fStateJac[7+2*7+1] * fOldCov[2*7+2] + fStateJac[7+3*7+1] * fOldCov[3*7+2] + fStateJac[7+4*7+1] * fOldCov[4*7+2] + fStateJac[7+5*7+1] * fOldCov[5*7+2] + fStateJac[7+6*7+1] * fOldCov[6*7+2];
      double JTC10 = fStateJac[7+1] * fOldCov[3*7] + fStateJac[7+7+1] * fOldCov[3*7+1] + fStateJac[7+2*7+1] * fOldCov[3*7+2] + fStateJac[7+3*7+1] * fOldCov[3*7+3] + fStateJac[7+4*7+1] * fOldCov[4*7+3] + fStateJac[7+5*7+1] * fOldCov[5*7+3] + fStateJac[7+6*7+1] * fOldCov[6*7+3];
      double JTC11 = fStateJac[7+1] * fOldCov[4*7] + fStateJac[7+7+1] * fOldCov[4*7+1] + fStateJac[7+2*7+1] * fOldCov[4*7+2] + fStateJac[7+3*7+1] * fOldCov[4*7+3] + fStateJac[7+4*7+1] * fOldCov[4*7+4] + fStateJac[7+5*7+1] * fOldCov[5*7+4] + fStateJac[7+6*7+1] * fOldCov[6*7+4];
      double JTC12 = fStateJac[7+1] * fOldCov[5*7] + fStateJac[7+7+1] * fOldCov[5*7+1] + fStateJac[7+2*7+1] * fOldCov[5*7+2] + fStateJac[7+3*7+1] * fOldCov[5*7+3] + fStateJac[7+4*7+1] * fOldCov[5*7+4] + fStateJac[7+5*7+1] * fOldCov[5*7+5] + fStateJac[7+6*7+1] * fOldCov[6*7+5];
      double JTC13 = fStateJac[7+1] * fOldCov[6*7] + fStateJac[7+7+1] * fOldCov[6*7+1] + fStateJac[7+2*7+1] * fOldCov[6*7+2] + fStateJac[7+3*7+1] * fOldCov[6*7+3] + fStateJac[7+4*7+1] * fOldCov[6*7+4] + fStateJac[7+5*7+1] * fOldCov[6*7+5] + fStateJac[7+6*7+1] * fOldCov[6*7+6];

      double JTC14 = fStateJac[7+2] * fOldCov[0] + fStateJac[7+7+2] * fOldCov[7] + fStateJac[7+2*7+2] * fOldCov[2*7] + fStateJac[7+3*7+2] * fOldCov[3*7] + fStateJac[7+4*7+2] * fOldCov[4*7] + fStateJac[7+5*7+2] * fOldCov[5*7] + fStateJac[7+6*7+2] * fOldCov[6*7];
      double JTC15 = fStateJac[7+2] * fOldCov[7] + fStateJac[7+7+2] * fOldCov[7+1] + fStateJac[7+2*7+2] * fOldCov[2*7+1] + fStateJac[7+3*7+2] * fOldCov[3*7+1] + fStateJac[7+4*7+2] * fOldCov[4*7+1] + fStateJac[7+5*7+2] * fOldCov[5*7+1] + fStateJac[7+6*7+2] * fOldCov[6*7+1];
      double JTC16 = fStateJac[7+2] * fOldCov[2*7] + fStateJac[7+7+2] * fOldCov[2*7+1] + fStateJac[7+2*7+2] * fOldCov[2*7+2] + fStateJac[7+3*7+2] * fOldCov[3*7+2] + fStateJac[7+4*7+2] * fOldCov[4*7+2] + fStateJac[7+5*7+2] * fOldCov[5*7+2] + fStateJac[7+6*7+2] * fOldCov[6*7+2];
      double JTC17 = fStateJac[7+2] * fOldCov[3*7] + fStateJac[7+7+2] * fOldCov[3*7+1] + fStateJac[7+2*7+2] * fOldCov[3*7+2] + fStateJac[7+3*7+2] * fOldCov[3*7+3] + fStateJac[7+4*7+2] * fOldCov[4*7+3] + fStateJac[7+5*7+2] * fOldCov[5*7+3] + fStateJac[7+6*7+2] * fOldCov[6*7+3];
      double JTC18 = fStateJac[7+2] * fOldCov[4*7] + fStateJac[7+7+2] * fOldCov[4*7+1] + fStateJac[7+2*7+2] * fOldCov[4*7+2] + fStateJac[7+3*7+2] * fOldCov[4*7+3] + fStateJac[7+4*7+2] * fOldCov[4*7+4] + fStateJac[7+5*7+2] * fOldCov[5*7+4] + fStateJac[7+6*7+2] * fOldCov[6*7+4];
      double JTC19 = fStateJac[7+2] * fOldCov[5*7] + fStateJac[7+7+2] * fOldCov[5*7+1] + fStateJac[7+2*7+2] * fOldCov[5*7+2] + fStateJac[7+3*7+2] * fOldCov[5*7+3] + fStateJac[7+4*7+2] * fOldCov[5*7+4] + fStateJac[7+5*7+2] * fOldCov[5*7+5] + fStateJac[7+6*7+2] * fOldCov[6*7+5];
      double JTC20 = fStateJac[7+2] * fOldCov[6*7] + fStateJac[7+7+2] * fOldCov[6*7+1] + fStateJac[7+2*7+2] * fOldCov[6*7+2] + fStateJac[7+3*7+2] * fOldCov[6*7+3] + fStateJac[7+4*7+2] * fOldCov[6*7+4] + fStateJac[7+5*7+2] * fOldCov[6*7+5] + fStateJac[7+6*7+2] * fOldCov[6*7+6];

      double JTC21 = fStateJac[7+3] * fOldCov[0] + fStateJac[7+7+3] * fOldCov[7] + fStateJac[7+2*7+3] * fOldCov[2*7] + fStateJac[7+3*7+3] * fOldCov[3*7] + fStateJac[7+4*7+3] * fOldCov[4*7] + fStateJac[7+5*7+3] * fOldCov[5*7] + fStateJac[7+6*7+3] * fOldCov[6*7];
      double JTC22 = fStateJac[7+3] * fOldCov[7] + fStateJac[7+7+3] * fOldCov[7+1] + fStateJac[7+2*7+3] * fOldCov[2*7+1] + fStateJac[7+3*7+3] * fOldCov[3*7+1] + fStateJac[7+4*7+3] * fOldCov[4*7+1] + fStateJac[7+5*7+3] * fOldCov[5*7+1] + fStateJac[7+6*7+3] * fOldCov[6*7+1];
      double JTC23 = fStateJac[7+3] * fOldCov[2*7] + fStateJac[7+7+3] * fOldCov[2*7+1] + fStateJac[7+2*7+3] * fOldCov[2*7+2] + fStateJac[7+3*7+3] * fOldCov[3*7+2] + fStateJac[7+4*7+3] * fOldCov[4*7+2] + fStateJac[7+5*7+3] * fOldCov[5*7+2] + fStateJac[7+6*7+3] * fOldCov[6*7+2];
      double JTC24 = fStateJac[7+3] * fOldCov[3*7] + fStateJac[7+7+3] * fOldCov[3*7+1] + fStateJac[7+2*7+3] * fOldCov[3*7+2] + fStateJac[7+3*7+3] * fOldCov[3*7+3] + fStateJac[7+4*7+3] * fOldCov[4*7+3] + fStateJac[7+5*7+3] * fOldCov[5*7+3] + fStateJac[7+6*7+3] * fOldCov[6*7+3];
      double JTC25 = fStateJac[7+3] * fOldCov[4*7] + fStateJac[7+7+3] * fOldCov[4*7+1] + fStateJac[7+2*7+3] * fOldCov[4*7+2] + fStateJac[7+3*7+3] * fOldCov[4*7+3] + fStateJac[7+4*7+3] * fOldCov[4*7+4] + fStateJac[7+5*7+3] * fOldCov[5*7+4] + fStateJac[7+6*7+3] * fOldCov[6*7+4];
      double JTC26 = fStateJac[7+3] * fOldCov[5*7] + fStateJac[7+7+3] * fOldCov[5*7+1] + fStateJac[7+2*7+3] * fOldCov[5*7+2] + fStateJac[7+3*7+3] * fOldCov[5*7+3] + fStateJac[7+4*7+3] * fOldCov[5*7+4] + fStateJac[7+5*7+3] * fOldCov[5*7+5] + fStateJac[7+6*7+3] * fOldCov[6*7+5];
      double JTC27 = fStateJac[7+3] * fOldCov[6*7] + fStateJac[7+7+3] * fOldCov[6*7+1] + fStateJac[7+2*7+3] * fOldCov[6*7+2] + fStateJac[7+3*7+3] * fOldCov[6*7+3] + fStateJac[7+4*7+3] * fOldCov[6*7+4] + fStateJac[7+5*7+3] * fOldCov[6*7+5] + fStateJac[7+6*7+3] * fOldCov[6*7+6];

      double JTC28 = fStateJac[7+4] * fOldCov[0] + fStateJac[7+7+4] * fOldCov[7] + fStateJac[7+2*7+4] * fOldCov[2*7] + fStateJac[7+3*7+4] * fOldCov[3*7] + fStateJac[7+4*7+4] * fOldCov[4*7] + fStateJac[7+5*7+4] * fOldCov[5*7] + fStateJac[7+6*7+4] * fOldCov[6*7];
      double JTC29 = fStateJac[7+4] * fOldCov[7] + fStateJac[7+7+4] * fOldCov[7+1] + fStateJac[7+2*7+4] * fOldCov[2*7+1] + fStateJac[7+3*7+4] * fOldCov[3*7+1] + fStateJac[7+4*7+4] * fOldCov[4*7+1] + fStateJac[7+5*7+4] * fOldCov[5*7+1] + fStateJac[7+6*7+4] * fOldCov[6*7+1];
      double JTC30 = fStateJac[7+4] * fOldCov[2*7] + fStateJac[7+7+4] * fOldCov[2*7+1] + fStateJac[7+2*7+4] * fOldCov[2*7+2] + fStateJac[7+3*7+4] * fOldCov[3*7+2] + fStateJac[7+4*7+4] * fOldCov[4*7+2] + fStateJac[7+5*7+4] * fOldCov[5*7+2] + fStateJac[7+6*7+4] * fOldCov[6*7+2];
      double JTC31 = fStateJac[7+4] * fOldCov[3*7] + fStateJac[7+7+4] * fOldCov[3*7+1] + fStateJac[7+2*7+4] * fOldCov[3*7+2] + fStateJac[7+3*7+4] * fOldCov[3*7+3] + fStateJac[7+4*7+4] * fOldCov[4*7+3] + fStateJac[7+5*7+4] * fOldCov[5*7+3] + fStateJac[7+6*7+4] * fOldCov[6*7+3];
      double JTC32 = fStateJac[7+4] * fOldCov[4*7] + fStateJac[7+7+4] * fOldCov[4*7+1] + fStateJac[7+2*7+4] * fOldCov[4*7+2] + fStateJac[7+3*7+4] * fOldCov[4*7+3] + fStateJac[7+4*7+4] * fOldCov[4*7+4] + fStateJac[7+5*7+4] * fOldCov[5*7+4] + fStateJac[7+6*7+4] * fOldCov[6*7+4];
      double JTC33 = fStateJac[7+4] * fOldCov[5*7] + fStateJac[7+7+4] * fOldCov[5*7+1] + fStateJac[7+2*7+4] * fOldCov[5*7+2] + fStateJac[7+3*7+4] * fOldCov[5*7+3] + fStateJac[7+4*7+4] * fOldCov[5*7+4] + fStateJac[7+5*7+4] * fOldCov[5*7+5] + fStateJac[7+6*7+4] * fOldCov[6*7+5];
      double JTC34 = fStateJac[7+4] * fOldCov[6*7] + fStateJac[7+7+4] * fOldCov[6*7+1] + fStateJac[7+2*7+4] * fOldCov[6*7+2] + fStateJac[7+3*7+4] * fOldCov[6*7+3] + fStateJac[7+4*7+4] * fOldCov[6*7+4] + fStateJac[7+5*7+4] * fOldCov[6*7+5] + fStateJac[7+6*7+4] * fOldCov[6*7+6];

      // last row
      cov_[6] =  JTC6;
      cov_[7+6] =  JTC13;
      cov_[2*7+6] =  JTC20;
      cov_[3*7+6] =  JTC27;
      cov_[4*7+6] =  JTC34;
      cov_[5*7+6] =  fStateJac[7+5] * fOldCov[6*7] + fStateJac[7+7+5] * fOldCov[6*7+1] + fStateJac[7+2*7+5] * fOldCov[6*7+2] + fStateJac[7+3*7+5] * fOldCov[6*7+3] + fStateJac[7+4*7+5] * fOldCov[6*7+4] + fStateJac[7+5*7+5] * fOldCov[6*7+5] + fStateJac[7+6*7+5] * fOldCov[6*7+6];
      cov_[6*7+6] = fOldCov[6*7+6];

      // loops are vectorizable by the compiler!
      for (int i=0; i<6; ++i) cov_[i] = JTC0 * fStateJac[7+i] + JTC1 * fStateJac[14+i] + JTC2 * fStateJac[21+i] + JTC3  * fStateJac[28+i] + JTC4  * fStateJac[35+i] + JTC5 * fStateJac[42+i] + JTC6 * fStateJac[49+i];
      for (int i=0; i<5; ++i) cov_[8+i] = JTC7 * fStateJac[8+i] + JTC8 * fStateJac[15+i] + JTC9 * fStateJac[22+i] + JTC10 * fStateJac[29+i] + JTC11 * fStateJac[36+i] + JTC12 * fStateJac[43+i] + (JTC13) * fStateJac[50+i];
      for (int i=0; i<4; ++i) cov_[16+i] = JTC14 * fStateJac[9+i] + JTC15 * fStateJac[16+i] + JTC16 * fStateJac[23+i] + (JTC17) * fStateJac[30+i] + (JTC18) * fStateJac[37+i] + (JTC19) * fStateJac[44+i] + (JTC20) * fStateJac[51+i];
      for (int i=0; i<3; ++i) cov_[24+i] = JTC21 * fStateJac[10+i] + JTC22 * fStateJac[17+i] + JTC23 * fStateJac[24+i] + (JTC24) * fStateJac[31+i] + (JTC25) * fStateJac[38+i] + (JTC26) * fStateJac[45+i] + (JTC27) * fStateJac[52+i];
      for (int i=0; i<2; ++i) cov_[32+i] = JTC28 * fStateJac[11+i] + JTC29 * fStateJac[18+i] + JTC30 * fStateJac[25+i] + (JTC31) * fStateJac[32+i] + (JTC32) * fStateJac[39+i] + (JTC33) * fStateJac[46+i] + (JTC34) * fStateJac[53+i];
      cov_[5*7+5] = (fStateJac[7+5] * fOldCov[0] + fStateJac[7+7+5] * fOldCov[7] + fStateJac[7+2*7+5] * fOldCov[2*7] + fStateJac[7+3*7+5] * fOldCov[3*7] + fStateJac[7+4*7+5] * fOldCov[4*7] + fStateJac[7+5*7+5] * fOldCov[5*7] + fStateJac[7+6*7+5] * fOldCov[6*7]) * fStateJac[7+5] + (fStateJac[7+5] * fOldCov[7] + fStateJac[7+7+5] * fOldCov[7+1] + fStateJac[7+2*7+5] * fOldCov[2*7+1] + fStateJac[7+3*7+5] * fOldCov[3*7+1] + fStateJac[7+4*7+5] * fOldCov[4*7+1] + fStateJac[7+5*7+5] * fOldCov[5*7+1] + fStateJac[7+6*7+5] * fOldCov[6*7+1]) * fStateJac[7+7+5] + (fStateJac[7+5] * fOldCov[2*7] + fStateJac[7+7+5] * fOldCov[2*7+1] + fStateJac[7+2*7+5] * fOldCov[2*7+2] + fStateJac[7+3*7+5] * fOldCov[3*7+2] + fStateJac[7+4*7+5] * fOldCov[4*7+2] + fStateJac[7+5*7+5] * fOldCov[5*7+2] + fStateJac[7+6*7+5] * fOldCov[6*7+2]) * fStateJac[7+2*7+5] + (fStateJac[7+5] * fOldCov[3*7] + fStateJac[7+7+5] * fOldCov[3*7+1] + fStateJac[7+2*7+5] * fOldCov[3*7+2] + fStateJac[7+3*7+5] * fOldCov[3*7+3] + fStateJac[7+4*7+5] * fOldCov[4*7+3] + fStateJac[7+5*7+5] * fOldCov[5*7+3] + fStateJac[7+6*7+5] * fOldCov[6*7+3]) * fStateJac[7+3*7+5] + (fStateJac[7+5] * fOldCov[4*7] + fStateJac[7+7+5] * fOldCov[4*7+1] + fStateJac[7+2*7+5] * fOldCov[4*7+2] + fStateJac[7+3*7+5] * fOldCov[4*7+3] + fStateJac[7+4*7+5] * fOldCov[4*7+4] + fStateJac[7+5*7+5] * fOldCov[5*7+4] + fStateJac[7+6*7+5] * fOldCov[6*7+4]) * fStateJac[7+4*7+5] + (fStateJac[7+5] * fOldCov[5*7] + fStateJac[7+7+5] * fOldCov[5*7+1] + fStateJac[7+2*7+5] * fOldCov[5*7+2] + fStateJac[7+3*7+5] * fOldCov[5*7+3] + fStateJac[7+4*7+5] * fOldCov[5*7+4] + fStateJac[7+5*7+5] * fOldCov[5*7+5] + fStateJac[7+6*7+5] * fOldCov[6*7+5]) * fStateJac[7+5*7+5] + (fStateJac[7+5] * fOldCov[6*7] + fStateJac[7+7+5] * fOldCov[6*7+1] + fStateJac[7+2*7+5] * fOldCov[6*7+2] + fStateJac[7+3*7+5] * fOldCov[6*7+3] + fStateJac[7+4*7+5] * fOldCov[6*7+4] + fStateJac[7+5*7+5] * fOldCov[6*7+5] + fStateJac[7+6*7+5] * fOldCov[6*7+6]) * fStateJac[7+6*7+5];

      // symmetric part
      cov_[7]   = cov_[1];
      cov_[2*7] = cov_[2];  cov_[2*7+1] = cov_[9];
      cov_[3*7] = cov_[3];  cov_[3*7+1] = cov_[10];  cov_[3*7+2] = cov_[17];
      cov_[4*7] = cov_[4];  cov_[4*7+1] = cov_[11];  cov_[4*7+2] = cov_[18];  cov_[4*7+3] = cov_[25];
      cov_[5*7] = cov_[5];  cov_[5*7+1] = cov_[12];  cov_[5*7+2] = cov_[19];  cov_[5*7+3] = cov_[26];  cov_[5*7+4] = cov_[33];
      cov_[6*7] = cov_[6];  cov_[6*7+1] = cov_[13];  cov_[6*7+2] = cov_[20];  cov_[6*7+3] = cov_[27];  cov_[6*7+4] = cov_[34];  cov_[6*7+5] = cov_[41];

      for (int i=0; i<7*7; ++i) cov_[i] += fNoise[i];

    } // finished propagate cov and add noise

    #ifdef DEBUG
      jacT.Print();
      if(calcCov) cov->Print();
    #endif

    //we arrived at the destination plane, if we point to the active area of the plane (if it is finite), and the distance is below threshold
    TVector3 pos(fStateJac[0], fStateJac[1], fStateJac[2]);
    if( plane.inActive(pos, directionAfter) &&
        plane.distance(pos) < MINSTEP) break;

  }

  memcpy(state7, fStateJac, 7*sizeof(double));

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
                         M8x7& P,
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
  M1x3&   R           = *((M1x3*) &P[0]);       // Start coordinates  [cm] 	(x,  y,  z)
  M1x3&   A           = *((M1x3*) &P[3]);       // Start directions 	      (ax, ay, az); 	ax^2+ay^2+az^2=1
  M1x3    SA          = {0.,0.,0.};             // Start directions derivatives dA/S
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
  M1x4     SU = {0.,0.,0.,0.};
  M1x3     H0 = {0.,0.,0.}, H1 = {0.,0.,0.}, H2 = {0.,0.,0.}, r = {0.,0.,0.};
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
	
        M1x3& dR = *((M1x3*) &P[i]);			            		// dR = (dX/dpN,  dY/dpN,  dZ/dpN)
        M1x3& dA = *((M1x3*) &P[i+3]);				           	// dA = (dAx/dpN, dAy/dpN, dAz/dpN); N = X,Y,Z,Ax,Ay,Az,q/p
        
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
        M1x3& dR = *((M1x3*) &P[i]);
        M1x3& dA = *((M1x3*) &P[i+3]);
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
                                const M1x4& SU,
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
  if (fabs(Step) > MINSTEP && !fNoMaterial){
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
