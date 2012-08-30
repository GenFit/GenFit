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
#include "GFPointPath.h"

#define MINSTEP 0.001   // minimum step [cm] for Runge Kutta and iteration to POCA
//#define DEBUG


RKTrackRep::~RKTrackRep(){
}

void RKTrackRep::Streamer(TBuffer &R__b)
{
   // Stream an object of class RKTrackRep.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RKTrackRep::Class(),this);
      initArrays();
      fCachePlane = fRefPlane;
      fCacheSpu = fSpu;

   } else {
      R__b.WriteClassBuffer(RKTrackRep::Class(),this);
   }
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
  if (fabs(mom) > 1.E-3) mom = fabs(fCharge/mom);
  else mom = 1.E3;

  double fac(mom);

  // get momentum/direction vector
  TVector3 momVec(aGFTrackCandPtr->getDirSeed());
  // if the user has set the dirSeed magnitude already to the momentum, the dirError must not be multiplied with the momentum magnitude
  if (momVec.Mag()-mom < 0.01) fac = 1;
  // set magnitude
  momVec.SetMag(mom);

  calcStateCov(aGFTrackCandPtr->getPosSeed(),
               momVec,
               aGFTrackCandPtr->getPosError(),
               aGFTrackCandPtr->getDirError()*fac);
}


void RKTrackRep::initArrays(){
  memset(fStateJac,0x00,(7+7*7)*sizeof(double));
  memset(fNoise,0x00,7*7*sizeof(double));
  memset(fStateJac,0x00,7*7*sizeof(double));

  memset(fJ_pM_5x7,0x00,5*7*sizeof(double));
  memset(fJ_pM_5x6,0x00,5*6*sizeof(double));
  memset(fJ_Mp_7x5,0x00,7*5*sizeof(double));
  memset(fJ_Mp_6x5,0x00,6*5*sizeof(double));
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

  fU = fRefPlane.getU();
  fV = fRefPlane.getV();
  fW = fRefPlane.getNormal();


  fCov(0,0) = fCharge*fCharge/pow(mom.Mag(),6.) *
              (mom.X()*mom.X() * momerr.X()*momerr.X()+
               mom.Y()*mom.Y() * momerr.Y()*momerr.Y()+
               mom.Z()*mom.Z() * momerr.Z()*momerr.Z());

  fCov(1,1) = pow((fU.X()/pw - fW.X()*pu/(pw*pw)),2.) * momerr.X()*momerr.X() +
              pow((fU.Y()/pw - fW.Y()*pu/(pw*pw)),2.) * momerr.Y()*momerr.Y() +
              pow((fU.Z()/pw - fW.Z()*pu/(pw*pw)),2.) * momerr.Z()*momerr.Z();

  fCov(2,2) = pow((fV.X()/pw - fW.X()*pv/(pw*pw)),2.) * momerr.X()*momerr.X() +
              pow((fV.Y()/pw - fW.Y()*pv/(pw*pw)),2.) * momerr.Y()*momerr.Y() +
              pow((fV.Z()/pw - fW.Z()*pv/(pw*pw)),2.) * momerr.Z()*momerr.Z();

  fCov(3,3) = poserr.X()*poserr.X() * fU.X()*fU.X() +
              poserr.Y()*poserr.Y() * fU.Y()*fU.Y() +
              poserr.Z()*poserr.Z() * fU.Z()*fU.Z();

  fCov(4,4) = poserr.X()*poserr.X() * fV.X()*fV.X() +
              poserr.Y()*poserr.Y() * fV.Y()*fV.Y() +
              poserr.Z()*poserr.Z() * fV.Z()*fV.Z();
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



void RKTrackRep::getState7(M1x7& state7) {
  getState7(state7, fState, fRefPlane, fSpu);
}


void RKTrackRep::getState7(M1x7& state7, const TMatrixD& state5, const GFDetPlane& pl, const double& spu) {

  fU = pl.getU();
  fV = pl.getV();
  fO = pl.getO();
  fW = pl.getNormal();

  state7[0] = fO.X() + state5(3,0)*fU.X() + state5(4,0)*fV.X(); // x
  state7[1] = fO.Y() + state5(3,0)*fU.Y() + state5(4,0)*fV.Y(); // y
  state7[2] = fO.Z() + state5(3,0)*fU.Z() + state5(4,0)*fV.Z(); // z

  state7[3] = spu * (fW.X() + state5(1,0)*fU.X() + state5(2,0)*fV.X()); // a_x
  state7[4] = spu * (fW.Y() + state5(1,0)*fU.Y() + state5(2,0)*fV.Y()); // a_y
  state7[5] = spu * (fW.Z() + state5(1,0)*fU.Z() + state5(2,0)*fV.Z()); // a_z

  // normalize dir
  double norm = 1. / sqrt(state7[3]*state7[3] + state7[4]*state7[4] + state7[5]*state7[5]);
  for (unsigned int i=3; i<6; ++i) state7[i] *= norm;

  state7[6] = state5(0,0); // q/p
}


TMatrixD RKTrackRep::getState5(const M1x7& state7, const GFDetPlane& pl, double& spu) {

  fU = pl.getU();
  fV = pl.getV();

  fPos.SetXYZ(state7[0], state7[1], state7[2]);
  fPos -= pl.getO();

  fDir.SetXYZ(state7[3], state7[4], state7[5]);

  // force A to be in normal direction and set spu accordingly
  double AtW = fDir * pl.getNormal();
  spu = 1.;
  if (AtW < 0) {
    fDir *= -1.;
    AtW *= -1.;
    spu = -1.;
  }

  TMatrixD state5(5,1);
  state5(0,0) = state7[6];
  state5(1,0) = fDir*fU / AtW;
  state5(2,0) = fDir*fV / AtW;
  state5(3,0) = fPos*fU;
  state5(4,0) = fPos*fV;

  return state5;
}



void RKTrackRep::transformPM7(const TMatrixD& in5x5, M7x7& out7x7,
                              const GFDetPlane& pl, const TMatrixD& state5, const double&  spu,
                              TMatrixD* Jac) {

  // get vectors and aux variables
  fU = pl.getU();
  fV = pl.getV();
  fW = pl.getNormal();

  fpTilde.SetXYZ(spu * (fW.X() + state5(1,0)*fU.X() + state5(2,0)*fV.X()), // a_x
                 spu * (fW.Y() + state5(1,0)*fU.Y() + state5(2,0)*fV.Y()), // a_y
                 spu * (fW.Z() + state5(1,0)*fU.Z() + state5(2,0)*fV.Z()));// a_z


  const double pTildeMag = fpTilde.Mag();
  const double pTildeMag2 = pTildeMag*pTildeMag;

  const double utpTildeOverpTildeMag2 = fU*fpTilde / pTildeMag2;
  const double vtpTildeOverpTildeMag2 = fV*fpTilde / pTildeMag2;

  //J_pM matrix is d(x,y,z,ax,ay,az,q/p) / d(q/p,u',v',u,v)   (out is 7x7)

   // d(x,y,z)/d(u)
  fJ_pM_5x7[21] = fU.X(); // [3][0]
  fJ_pM_5x7[22] = fU.Y(); // [3][1]
  fJ_pM_5x7[23] = fU.Z(); // [3][2]
  // d(x,y,z)/d(v)
  fJ_pM_5x7[28] = fV.X(); // [4][2]
  fJ_pM_5x7[29] = fV.Y(); // [4][2]
  fJ_pM_5x7[30] = fV.Z(); // [4][2]
  // d(q/p)/d(q/p)
  fJ_pM_5x7[6] = 1.; // not needed for array matrix multiplication
  // d(ax,ay,az)/d(u')
  double fact = spu / pTildeMag;
  fJ_pM_5x7[10] = fact * ( fU.X() - fpTilde.X()*utpTildeOverpTildeMag2 ); // [1][3]
  fJ_pM_5x7[11] = fact * ( fU.Y() - fpTilde.Y()*utpTildeOverpTildeMag2 ); // [1][4]
  fJ_pM_5x7[12] = fact * ( fU.Z() - fpTilde.Z()*utpTildeOverpTildeMag2 ); // [1][5]
  // d(ax,ay,az)/d(v')
  fJ_pM_5x7[17] = fact * ( fV.X() - fpTilde.X()*vtpTildeOverpTildeMag2 ); // [2][3]
  fJ_pM_5x7[18] = fact * ( fV.Y() - fpTilde.Y()*vtpTildeOverpTildeMag2 ); // [2][4]
  fJ_pM_5x7[19] = fact * ( fV.Z() - fpTilde.Z()*vtpTildeOverpTildeMag2 ); // [2][5]


  // since the Jacobian contains a lot of zeros, and the resulting cov has to be symmetric,
  // the multiplication can be done much faster directly on array level
  // out = J_pM^T * in5x5 * J_pM
  const M5x5& in5x5_ = *((M5x5*) in5x5.GetMatrixArray());
  RKTools::J_pMTxcov5xJ_pM(fJ_pM_5x7, in5x5_, out7x7);

  if (Jac!=NULL){
    Jac->ResizeTo(5,7);
    *Jac = TMatrixD(5,7, &(fJ_pM_5x7[0]));
  }
}


void RKTrackRep::transformPM6(const TMatrixD& in5x5, M6x6& out6x6,
                              const GFDetPlane& pl, const TMatrixD& state5, const double&  spu,
                              TMatrixD* Jac) {

  // get vectors and aux variables
  fU = pl.getU();
  fV = pl.getV();

  fpTilde.SetXYZ(spu * (fW.X() + state5(1,0)*fU.X() + state5(2,0)*fV.X()), // a_x
                 spu * (fW.Y() + state5(1,0)*fU.Y() + state5(2,0)*fV.Y()), // a_y
                 spu * (fW.Z() + state5(1,0)*fU.Z() + state5(2,0)*fV.Z()));// a_z

  const double pTildeMag = fpTilde.Mag();
  const double pTildeMag2 = pTildeMag*pTildeMag;

  const double utpTildeOverpTildeMag2 = fU*fpTilde / pTildeMag2;
  const double vtpTildeOverpTildeMag2 = fV*fpTilde / pTildeMag2;

  //J_pM matrix is d(x,y,z,px,py,pz) / d(q/p,u',v',u,v)       (out is 6x6)

  const double qop = state5(0,0);
  const double p = fCharge/qop; // momentum

  // d(px,py,pz)/d(q/p)
  double fact = -1. * p / (pTildeMag * qop);
  fJ_pM_5x6[3] = fact * fpTilde.X(); // [0][3]
  fJ_pM_5x6[4] = fact * fpTilde.Y(); // [0][4]
  fJ_pM_5x6[5] = fact * fpTilde.Z(); // [0][5]
  // d(px,py,pz)/d(u')
  fact = p * spu / pTildeMag;
  fJ_pM_5x6[9]  = fact * ( fU.X() - fpTilde.X()*utpTildeOverpTildeMag2 ); // [1][3]
  fJ_pM_5x6[10] = fact * ( fU.Y() - fpTilde.Y()*utpTildeOverpTildeMag2 ); // [1][4]
  fJ_pM_5x6[11] = fact * ( fU.Z() - fpTilde.Z()*utpTildeOverpTildeMag2 ); // [1][5]
  // d(px,py,pz)/d(v')
  fJ_pM_5x6[15] = fact * ( fV.X() - fpTilde.X()*vtpTildeOverpTildeMag2 ); // [2][3]
  fJ_pM_5x6[16] = fact * ( fV.Y() - fpTilde.Y()*vtpTildeOverpTildeMag2 ); // [2][4]
  fJ_pM_5x6[17] = fact * ( fV.Z() - fpTilde.Z()*vtpTildeOverpTildeMag2 ); // [2][5]
  // d(x,y,z)/d(u)
  fJ_pM_5x6[18] = fU.X(); // [3][0]
  fJ_pM_5x6[19] = fU.Y(); // [3][1]
  fJ_pM_5x6[20] = fU.Z(); // [3][2]
  // d(x,y,z)/d(v)
  fJ_pM_5x6[24] = fV.X(); // [4][0]
  fJ_pM_5x6[25] = fV.Y(); // [4][1]
  fJ_pM_5x6[26] = fV.Z(); // [4][2]


  // do the transformation
  // out = J_pM^T * in5x5 * J_pM
  const M5x5& in5x5_ = *((M5x5*) in5x5.GetMatrixArray());
  RKTools::J_pMTxcov5xJ_pM(fJ_pM_5x6, in5x5_, out6x6);

  if (Jac!=NULL){
    Jac->ResizeTo(5,6);
    *Jac = TMatrixD(5,6, &(fJ_pM_5x6[0]));
  }
}


void RKTrackRep::transformM7P(const M7x7& in7x7, TMatrixD& out5x5,
                              const GFDetPlane& pl, const M1x7& state7,
                              TMatrixD* Jac) {

  out5x5.ResizeTo(5, 5);

  // get vectors and aux variables
  fU = pl.getU();
  fV = pl.getV();
  fW = pl.getNormal();

  fDir.SetXYZ(state7[3], state7[4], state7[5]);

  const double AtU = fDir*fU;
  const double AtV = fDir*fV;
  const double AtW = fDir*fW;

  // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,ax,ay,az,q/p)   (in is 7x7)

  // d(u')/d(ax,ay,az)
  double fact = 1./(AtW*AtW);
  fJ_Mp_7x5[16] = fact * (fU.X()*AtW - fW.X()*AtU); // [3][1]
  fJ_Mp_7x5[21] = fact * (fU.Y()*AtW - fW.Y()*AtU); // [4][1]
  fJ_Mp_7x5[26] = fact * (fU.Z()*AtW - fW.Z()*AtU); // [5][1]
  // d(v')/d(ax,ay,az)
  fJ_Mp_7x5[17] = fact * (fV.X()*AtW - fW.X()*AtV); // [3][2]
  fJ_Mp_7x5[22] = fact * (fV.Y()*AtW - fW.Y()*AtV); // [4][2]
  fJ_Mp_7x5[27] = fact * (fV.Z()*AtW - fW.Z()*AtV); // [5][2]
  // d(q/p)/d(q/p)
  fJ_Mp_7x5[30] = 1.; // [6][0]  - not needed for array matrix multiplication
  //d(u)/d(x,y,z)
  fJ_Mp_7x5[3]  = fU.X(); // [0][3]
  fJ_Mp_7x5[8]  = fU.Y(); // [1][3]
  fJ_Mp_7x5[13] = fU.Z(); // [2][3]
  //d(v)/d(x,y,z)
  fJ_Mp_7x5[4]  = fV.X(); // [0][4]
  fJ_Mp_7x5[9]  = fV.Y(); // [1][4]
  fJ_Mp_7x5[14] = fV.Z(); // [2][4]


  // since the Jacobian contains a lot of zeros, and the resulting cov has to be symmetric,
  // the multiplication can be done much faster directly on array level
  // out5x5 = J_Mp^T * in * J_Mp
  M5x5& out5x5_ = *((M5x5*) out5x5.GetMatrixArray());
  RKTools::J_MpTxcov7xJ_Mp(fJ_Mp_7x5, in7x7, out5x5_);

  if (Jac!=NULL){
    Jac->ResizeTo(7,5);
    *Jac = TMatrixD(7,5, &(fJ_Mp_7x5[0]));
  }
}


void RKTrackRep::transformM6P(const M6x6& in6x6, TMatrixD& out5x5,
                              const GFDetPlane& pl, const M1x7& state7,
                              TMatrixD* Jac) {

  out5x5.ResizeTo(5, 5);

  // get vectors and aux variables
  fU = pl.getU();
  fV = pl.getV();
  fW = pl.getNormal();

  fDir.SetXYZ(state7[3], state7[4], state7[5]);

  const double AtU = fDir*fU;
  const double AtV = fDir*fV;
  const double AtW = fDir*fW;

  // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,px,py,pz)       (in is 6x6)

  const double qop = state7[6];
  const double p = fCharge/qop; // momentum

  //d(u)/d(x,y,z)
  fJ_Mp_6x5[3]  = fU.X(); // [0][3]
  fJ_Mp_6x5[8]  = fU.Y(); // [1][3]
  fJ_Mp_6x5[13] = fU.Z(); // [2][3]
  //d(v)/d(x,y,z)
  fJ_Mp_6x5[4]  = fV.X(); // [0][4]
  fJ_Mp_6x5[9]  = fV.Y(); // [1][4]
  fJ_Mp_6x5[14] = fV.Z(); // [2][4]
  // d(q/p)/d(px,py,pz)
  double fact = (-1.) * qop / p;
  fJ_Mp_6x5[15] = fact * fDir.X(); // [3][0]
  fJ_Mp_6x5[20] = fact * fDir.Y(); // [4][0]
  fJ_Mp_6x5[25] = fact * fDir.Z(); // [5][0]
  // d(u')/d(px,py,pz)
  fact = 1./(p*AtW*AtW);
  fJ_Mp_6x5[16] = fact * (fU.X()*AtW - fW.X()*AtU); // [3][1]
  fJ_Mp_6x5[21] = fact * (fU.Y()*AtW - fW.Y()*AtU); // [4][1]
  fJ_Mp_6x5[26] = fact * (fU.Z()*AtW - fW.Z()*AtU); // [5][1]
  // d(v')/d(px,py,pz)
  fJ_Mp_6x5[17] = fact * (fV.X()*AtW - fW.X()*AtV); // [3][2]
  fJ_Mp_6x5[22] = fact * (fV.Y()*AtW - fW.Y()*AtV); // [4][2]
  fJ_Mp_6x5[27] = fact * (fV.Z()*AtW - fW.Z()*AtV); // [5][2]

  // do the transformation
  // out5x5 = J_Mp^T * in * J_Mp
  M5x5& out5x5_ = *((M5x5*) out5x5.GetMatrixArray());
  RKTools::J_MpTxcov6xJ_Mp(fJ_Mp_6x5, in6x6, out5x5_);

  if (Jac!=NULL){
    Jac->ResizeTo(6,5);
    *Jac = TMatrixD(6,5, &(fJ_Mp_6x5[0]));;
  }
}



TVector3 RKTrackRep::getPos(const GFDetPlane& pl){
#ifdef DEBUG
  std::cout << "RKTrackRep::getPos()\n";
#endif
  M1x7 state7;
  getState7(state7);
  if(pl!=fRefPlane) Extrap(pl, state7);
  return TVector3(state7[0], state7[1], state7[2]);
}


TVector3 RKTrackRep::getMom(const GFDetPlane& pl){
#ifdef DEBUG
  std::cout << "RKTrackRep::getMom()\n";
#endif
  M1x7 state7;
  getState7(state7);
  if(pl!=fRefPlane) Extrap(pl, state7);

  TVector3 mom(state7[3], state7[4], state7[5]);
  mom.SetMag(fCharge/state7[6]);
  return mom;
}


void RKTrackRep::getPosMom(const GFDetPlane& pl,TVector3& pos, TVector3& mom){
#ifdef DEBUG
  std::cout << "RKTrackRep::getPosMom()\n";
#endif
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

#ifdef DEBUG
  std::cout << "RKTrackRep::extrapolateToPoint()\n";
#endif

  static const unsigned int maxIt(300);

  M1x7 state7;
  getState7(state7);

  double coveredDistance(0.);

  GFDetPlane pl;
  unsigned int iterations(0);

  while(true){
    fDir.SetXYZ(state7[3], state7[4], state7[5]);
    pl.setON(pos, fDir);
    coveredDistance =  this->Extrap(pl, state7, NULL, true);

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

#ifdef DEBUG
  std::cout << "RKTrackRep::extrapolateToLine(), (x,y) = (" << point1.X() << ", " << point1.Y() << ")\n";
#endif

  static const unsigned int maxIt(300);

  M1x7 state7;
  getState7(state7);

  double coveredDistance(0.);

  GFDetPlane pl;
  unsigned int iterations(0);

  while(true){
    pl.setO(point1);
    fDir.SetXYZ(state7[3], state7[4], state7[5]);
    pl.setU(fDir.Cross(point2-point1));
    pl.setV(point2-point1);
    coveredDistance = this->Extrap(pl, state7, NULL, true);

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
  
#ifdef DEBUG
  std::cout << "RKTrackRep::extrapolate(pl, statePred, covPred)\n";
#endif

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

#ifdef DEBUG
  std::cout << "RKTrackRep::extrapolate(pl, statePred)\n";
#endif

  M1x7 state7;
  getState7(state7);
  double coveredDistance = Extrap(pl, state7);
  double spu;
  statePred.ResizeTo(5,1);
  statePred = getState5(state7, pl, spu);

  return coveredDistance;
}


double RKTrackRep::stepalong(double h, TVector3& pos, TVector3& dir){

#ifdef DEBUG
  std::cout << "RKTrackRep::stepalong()\n";
#endif

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



double RKTrackRep::Extrap( const GFDetPlane& plane, M1x7& state7, M7x7* cov, bool onlyOneStep) {

  static const unsigned int maxNumIt(200);
  unsigned int numIt(0);

  const bool calcCov(cov!=NULL);

  // set initial state
  memcpy(fStateJac, state7, 7*sizeof(double));

  double coveredDistance(0.);
  double sumDistance(0.);

  while(true){

    #ifdef DEBUG
      std::cout << "\n============ RKTrackRep::Extrap loop nr. " << numIt << " ============\n";
    #endif

    if(numIt++ > maxNumIt){
      GFException exc("RKTrackRep::Extrap ==> maximum number of iterations exceeded",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    // initialize fStateJac with unit matrix; last entry is not 1 but q/p
    if(calcCov){
      memset(&fStateJac[7],0x00,49*sizeof(double));
      for(int i=0; i<7; ++i){
        fStateJac[(i+1)*7+i] = 1.;
      }
    }

    fDirectionBefore.SetXYZ(fStateJac[3], fStateJac[4], fStateJac[5]); // direction before propagation

    // propagation
    std::vector<GFPointPath> points;

    bool checkJacProj = true;

    if( ! this->RKutta(plane, fStateJac, coveredDistance, points, checkJacProj, calcCov, onlyOneStep) ) {
      GFException exc("RKTrackRep::Extrap ==> Runge Kutta propagation failed",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    fPos.SetXYZ(fStateJac[0], fStateJac[1], fStateJac[2]);
    if (!fNoMaterial) points.push_back(GFPointPath(fPos, 0)); // add last point

    #ifdef DEBUG
      std::cout<<"Original points \n";
      for (unsigned int i=0; i<points.size(); ++i){
        points[i].Print();
      }
      std::cout<<"\n";
    #endif

    fDirectionAfter.SetXYZ(fStateJac[3], fStateJac[4], fStateJac[5]); // direction after propagation

    sumDistance+=coveredDistance;

    // filter points
    if (!fNoMaterial) { // points are only filled if mat fx are on
      if(points.size() > 2){ // check if there are at least three points
        double firstStep(points[0].getPath());
        for (unsigned int i=points.size()-2; i>0; --i){
          if (points[i].getPath() * firstStep < 0 || fabs(points[i].getPath()) < MINSTEP){
            points[i-1].addToPath(points[i].getPath());
            points.erase(points.begin()+i);
          }
        }
      }
      #ifdef DEBUG
        std::cout<<"Filtered points \n";
        for (unsigned int i=0; i<points.size(); ++i){
          points[i].Print();
        }
        std::cout<<"\n";
      #endif
    }


    if(calcCov) memset(fNoise,0x00,7*7*sizeof(double)); // set fNoise to 0


    // call MatFX
    unsigned int nPoints(points.size());
    if (!fNoMaterial && nPoints>0){
      // momLoss has a sign - negative loss means momentum gain
      double momLoss = GFMaterialEffects::getInstance()->effects(points,
                                                                 fabs(fCharge/fStateJac[6]), // momentum
                                                                 fPdg,
                                                                 fXX0,
                                                                 calcCov,
                                                                 fNoise,
                                                                 &(fStateJac[7]),
                                                                 &fDirectionBefore,
                                                                 &fDirectionAfter);

      #ifdef DEBUG
        std::cout << "momLoss: " << momLoss << " GeV; relative: " << momLoss/fabs(fCharge/fStateJac[6]) << "\n";
      #endif

      // do momLoss only for defined 1/momentum .ne.0
      if(fabs(fStateJac[6])>1.E-10) fStateJac[6] = fCharge/(fabs(fCharge/fStateJac[6])-momLoss);
    } // finished MatFX

    if(calcCov){ // propagate cov and add noise
      memcpy(fOldCov, (*cov), 7*7*sizeof(double));
      M7x7& cov_ = (*cov);

      // numerical check:
      for(unsigned int i=0; i<7*7; ++i){
        if(fabs((*cov)[i]) > 1.E100){
          GFException exc("RKTrackRep::Extrap ==> covariance matrix exceeds numerical limits",__LINE__,__FILE__);
          exc.setFatal();
          throw exc;
        }
      }

      // cov = Jac^T * oldCov * Jac;
      // last column of jac is [0,0,0,0,0,0,1]
      // cov is symmetric
      M7x7& J_MM = *((M7x7*) &(fStateJac[7]));
      RKTools::J_MMTxcov7xJ_MM(J_MM, fOldCov, cov_);

      // add noise to cov
      for (int i=0; i<7*7; ++i) cov_[i] += fNoise[i];

    } // finished propagate cov and add noise

    if (onlyOneStep) break;

    //we arrived at the destination plane, if we point to the active area of the plane (if it is finite), and the distance is below threshold
    if( plane.distance(fPos) < MINSTEP  &&  plane.inActive(fPos, fDirectionAfter)) {
      // check if Jacobian has been projected onto plane; Otherwise make another iteration
      if (calcCov && !checkJacProj && nPoints>0){
        #ifdef DEBUG
          std::cout << "Jacobian was not projected onto destination plane -> one more iteration. \n";
        #endif
        continue;
      }
      #ifdef DEBUG
        std::cout << "arrived at plane with a distance of  " << plane.distance(fPos) << " cm left. \n";
      #endif
      break;
    }

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
//    P      - initial parameters (coordinates(cm), direction,
//             charge/momentum (Gev-1) and derivatives this parameters  (8x7)            
//         
//    X        	Y        	Z        	Ax       	Ay       	Az       	q/P                   
//    P[ 0]   	P[ 1]   	P[ 2]   	P[ 3]   	P[ 4]   	P[ 5]   	P[ 6]  
//
//    dX/dp    	dY/dp    	dZ/dp    	dAx/dp   	dAy/dp   	dAz/dp   	d(q/P)/dp
//    P[ 7]   	P[ 8]   	P[ 9]   	P[10]   	P[11]   	P[12]   	P[13]   			      d()/dp1  
//
//    P[14]   	P[15]   	P[16]   	P[17]   	P[18]   	P[19]   	P[20]   		      	d()/dp2        
//    ............................................................................		d()/dpND       
//                                                                                  
// Authors: R.Brun, M.Hansroul, V.Perevoztchikov (Geant3)                           
//  
bool RKTrackRep::RKutta (const GFDetPlane& plane,
                         M8x7& P,
                         double& coveredDistance,
                         std::vector<GFPointPath>& points,
                         bool& checkJacProj,
                         bool calcCov,
                         bool onlyOneStep) {

  // important fixed numbers
  static const double EC     = 0.000149896229;  // c/(2*10^12) resp. c/2Tera
  static const double P3     = 1./3.;           // 1/3
  static const int    ND     = 56;              // number of variables for derivatives calculation
  static const int    ND1    = ND-7;            // = 49
  // limits, check-values, etc. Can be tuned!
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
  bool    atPlane = false;                      // stepper thinks that the plane will be reached in that step -> linear extrapolation and projection of jacobian
  bool    momLossExceeded = false;              // stepper had to limit stepsize due to momentum loss -> no next RKutta loop, no linear extrapolation
  fPos.SetXYZ(R[0],R[1],R[2]);                  // position
  fDir.SetXYZ(A[0],A[1],A[2]);                  // direction
  double   momentum   = fabs(fCharge/P[6]);     // momentum [GeV]
  double   relMomLoss = 0;                      // relative momentum loss in RKutta
  double   deltaAngle = 0.;                     // total angle by which the momentum has changed during extrapolation
  double   An(0), S(0), Sl(0), S3(0), S4(0), PS2(0), CBA(0);
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
    std::cout << "position: "; fPos.Print();
    std::cout << "direction: "; fDir.Print();
    std::cout << "momentum: " << momentum << " GeV\n";
    std::cout << "destination: "; plane.Print();
  #endif

  checkJacProj = false;

  // check momentum
  if(momentum < Pmin){
    std::ostringstream sstream;
    sstream << "RKTrackRep::RKutta ==> momentum too low: " << fabs(fCharge/P[6])*1000. << " MeV";
    GFException exc(sstream.str(),__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  

  // make SU vector point away from origin
  const TVector3 W = plane.getNormal();
  if(W*plane.getO() > 0){SU[0] =     W.X();  SU[1] =     W.Y();  SU[2] =     W.Z();}
  else                  {SU[0] = -1.*W.X();  SU[1] = -1.*W.Y();  SU[2] = -1.*W.Z();  }
  SU[3] = plane.distance(0., 0., 0.);

  unsigned int counter(0);

  // Step estimation (signed)
  fH = GFFieldManager::getFieldVal(fPos);
  S = estimateStep(points, fPos, fDir, SU, fH, plane, momentum, relMomLoss, deltaAngle, momLossExceeded, atPlane);
  if (fabs(S) < 0.001*MINSTEP) {
    #ifdef DEBUG
      std::cout << " RKutta - step too small -> break \n";
    #endif
    counter += 1; // skip the main loop, go to linear extrapolation step (will be skipped) and just project jacobian
  }

  //
  // Main loop of Runge-Kutta method
  //
  while (fabs(S) >= MINSTEP || counter == 0) {

    if(++counter > maxNumIt){
      GFException exc("RKTrackRep::RKutta ==> maximum number of iterations exceeded",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    #ifdef DEBUG
      std::cout << "------ RKutta main loop nr. " << counter-1 << " ------\n";
    #endif

    //
    // Runge Kutta Extrapolation
    //
    S3 = P3*S;
    S4 = 0.25*S;
    PS2 = Pinv*S;
    
    // First point
    r[0] = R[0];           r[1] = R[1];           r[2]=R[2];
    fPos.SetXYZ(r[0], r[1], r[2]); // vector of start coordinates R0	(x, y, z)
    fH = GFFieldManager::getFieldVal(fPos);				// magnetic field in 10^-4 T = kGauss
    H0[0] = PS2*fH.X(); H0[1] = PS2*fH.Y(); H0[2] = PS2*fH.Z(); 		// H0 is PS2*(Hx, Hy, Hz) @ R0
    A0 = A[1]*H0[2]-A[2]*H0[1]; B0 = A[2]*H0[0]-A[0]*H0[2]; C0 = A[0]*H0[1]-A[1]*H0[0]; // (ax, ay, az) x H0
    A2 = A[0]+A0              ; B2 = A[1]+B0              ; C2 = A[2]+C0              ; // (A0, B0, C0) + (ax, ay, az)
    A1 = A2+A[0]              ; B1 = B2+A[1]              ; C1 = C2+A[2]              ; // (A0, B0, C0) + 2*(ax, ay, az)
      
    // Second point
    r[0] += A1*S4;         r[1] += B1*S4;         r[2] += C1*S4;
    fPos.SetXYZ(r[0], r[1], r[2]);
    fH = GFFieldManager::getFieldVal(fPos);
    H1[0] = fH.X()*PS2; H1[1] = fH.Y()*PS2; H1[2] = fH.Z()*PS2;	// H1 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * [(A0, B0, C0) + 2*(ax, ay, az)]
    A3 = B2*H1[2]-C2*H1[1]+A[0]; B3 = C2*H1[0]-A2*H1[2]+A[1]; C3 = A2*H1[1]-B2*H1[0]+A[2]; // (A2, B2, C2) x H1 + (ax, ay, az)
    A4 = B3*H1[2]-C3*H1[1]+A[0]; B4 = C3*H1[0]-A3*H1[2]+A[1]; C4 = A3*H1[1]-B3*H1[0]+A[2]; // (A3, B3, C3) x H1 + (ax, ay, az)
    A5 = A4-A[0]+A4            ; B5 = B4-A[1]+B4            ; C5 = C4-A[2]+C4            ; //    2*(A4, B4, C4) - (ax, ay, az)

    // Last point
    r[0]=R[0]+S*A4;         r[1]=R[1]+S*B4;         r[2]=R[2]+S*C4;  //setup.Field(r,H2);
    fPos.SetXYZ(r[0], r[1], r[2]);
    fH = GFFieldManager::getFieldVal(fPos);
    H2[0] = fH.X()*PS2;  H2[1] = fH.Y()*PS2;  H2[2] = fH.Z()*PS2;	// H2 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * (A4, B4, C4)
    A6 = B5*H2[2]-C5*H2[1]; B6 = C5*H2[0]-A5*H2[2]; C6 = A5*H2[1]-B5*H2[0]; // (A5, B5, C5) x H2
    
    // update paths
    coveredDistance += S;				// add stepsize to way (signed)
    Way  += fabs(S);

    // check way limit
    if(Way > Wmax){
      std::ostringstream sstream;
      sstream << "RKTrackRep::RKutta ==> Total extrapolation length is longer than length limit : " << Way << " cm !";
      GFException exc(sstream.str(),__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    //
    // Derivatives of track parameters in last point
    //
    if(calcCov){
      // d(x, y, z)/d(x, y, z) submatrix is unit matrix
      P[7] = 1;  P[15] = 1;  P[23] = 1;
      // d(ax, ay, az)/d(ax, ay, az) submatrix is 0
      // start with d(x, y, z)/d(ax, ay, az)
      for(int i=4*7; i<ND; i+=7) {				// i = 28, 35, 42, 49;    ND = 56;	ND1 = 49; rows of Jacobian
	
        M1x3& dR = *((M1x3*) &P[i]);			            		// dR = (dX/dpN,  dY/dpN,  dZ/dpN)
        M1x3& dA = *((M1x3*) &P[i+3]);				           	// dA = (dAx/dpN, dAy/dpN, dAz/dpN); N = X,Y,Z,Ax,Ay,Az,q/p
        
        if(i==ND1) {dA[0]*=P[6]; dA[1]*=P[6]; dA[2]*=P[6];}

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
        
        if(i==ND1) {
          dR[0] += (dA2+dA3+dA4)*S3/P[6];  dA[0] = (dA0+dA3+dA3+dA5+dA6)*P3/P[6]; // dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
          dR[1] += (dB2+dB3+dB4)*S3/P[6];  dA[1] = (dB0+dB3+dB3+dB5+dB6)*P3/P[6]; // dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
          dR[2] += (dC2+dC3+dC4)*S3/P[6];  dA[2] = (dC0+dC3+dC3+dC5+dC6)*P3/P[6];
        }
        else {
          dR[0] += (dA2+dA3+dA4)*S3;  dA[0] = (dA0+dA3+dA3+dA5+dA6)*P3;	// dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
          dR[1] += (dB2+dB3+dB4)*S3;  dA[1] = (dB0+dB3+dB3+dB5+dB6)*P3;	// dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
          dR[2] += (dC2+dC3+dC4)*S3;  dA[2] = (dC0+dC3+dC3+dC5+dC6)*P3;
        }
      }
    }

    //
    // Track parameters in last point
    //   
    R[0] += (A2+A3+A4)*S3;   A[0] += (SA[0]=(A0+A3+A3+A5+A6)*P3-A[0]);  // R  = R0 + S3*[(A2, B2, C2) +   (A3, B3, C3) + (A4, B4, C4)]
    R[1] += (B2+B3+B4)*S3;   A[1] += (SA[1]=(B0+B3+B3+B5+B6)*P3-A[1]);  // A  =     1/3*[(A0, B0, C0) + 2*(A3, B3, C3) + (A5, B5, C5) + (A6, B6, C6)]
    R[2] += (C2+C3+C4)*S3;   A[2] += (SA[2]=(C0+C3+C3+C5+C6)*P3-A[2]); 	// SA = A_new - A_old
    fPos.SetXYZ(R[0], R[1], R[2]);

    // normalize A
    CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);	// 1/|A|
    A[0] *= CBA; A[1] *= CBA; A[2] *= CBA;
    fDir.SetXYZ(A[0], A[1], A[2]);
    
    if (onlyOneStep) return(true);

    // if stepsize has been limited by material, break the loop and return. No linear extrapolation!
    if (momLossExceeded) {
      #ifdef DEBUG
        std::cout<<" momLossExceeded -> return(true); \n";
      #endif
      return(true);
    }

    // estimate Step for next loop or linear extrapolation
    Sl = S;	// last S used
    S = estimateStep(points, fPos, fDir, SU, fH, plane, momentum, relMomLoss, deltaAngle, momLossExceeded, atPlane);

    if (atPlane && fabs(S) < MINSTEP) {
      #ifdef DEBUG
        std::cout<<" (atPlane && fabs(S) < MINSTEP) -> break; \n";
      #endif
      break; // else if at plane: do linear extrapolation
    }
    if (momLossExceeded && fabs(S) < MINSTEP) {
      #ifdef DEBUG
        std::cout<<" (momLossExceeded && fabs(S) < MINSTEP) -> return(true); \n";
      #endif
      return(true); // no linear extrapolation!
    }

    // check if total angle is bigger than AngleMax. Can happen if a curler should be fitted and it does not hit the active area of the next plane.
    if (fabs(deltaAngle) > AngleMax){
      std::ostringstream sstream;
      sstream << "RKTrackRep::RKutta ==> Do not get to an active plane! Already extrapolated " << deltaAngle * 180 / TMath::Pi() << "°.";
      GFException exc(sstream.str(),__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    // check if we went back and forth multiple times -> we don't come closer to the plane!
    if (counter > 3){
      if (S                          *points[counter-1].getPath() < 0 &&
          points[counter-1].getPath()*points[counter-2].getPath() < 0 &&
          points[counter-2].getPath()*points[counter-3].getPath() < 0){
        GFException exc("RKTrackRep::RKutta ==> Do not get closer to plane!",__LINE__,__FILE__);
        exc.setFatal();
        throw exc;
      }
    }

  } //end of main loop
  

  //
  // linear extrapolation to surface
  //
  if (atPlane) {

    if (fabs(Sl) > 0.001*MINSTEP){
      #ifdef DEBUG
        std::cout << " RKutta - linear extrapolation to surface\n";
      #endif
      Sl = 1./Sl;        // Sl = inverted last Stepsize Sl
      
      // normalize SA
      SA[0]*=Sl;  SA[1]*=Sl;  SA[2]*=Sl; // SA/Sl = delta A / delta way; local derivative of A with respect to the length of the way

      // calculate A
      A[0] += SA[0]*S;   	// S  = distance to surface
      A[1] += SA[1]*S;	  // A = A + S * SA*Sl
      A[2] += SA[2]*S;

      // normalize A
      CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);  // 1/|A|
      A[0] *= CBA; A[1] *= CBA; A[2] *= CBA;

      R[0] += S*(A[0]-0.5*S*SA[0]);    // P = R + S*(A - 0.5*S*SA); approximation for final point on surface
      R[1] += S*(A[1]-0.5*S*SA[1]);
      R[2] += S*(A[2]-0.5*S*SA[2]);
    }
#ifdef DEBUG
    else {
      std::cout << " RKutta - last stepsize too small -> can't do linear extrapolation! \n";
    }
#endif

    //
    // Project Jacobian of extrapolation onto destination plane
    //
    if(calcCov){
      if (checkJacProj && points.size()>0){
        GFException exc("RKTrackRep::Extrap ==> covariance is projected onto destination plane again",__LINE__,__FILE__);
        throw exc;
      }
      checkJacProj = true;
      #ifdef DEBUG
        std::cout << "  Project Jacobian of extrapolation onto destination plane\n";
      #endif
      An = A[0]*SU[0] + A[1]*SU[1] + A[2]*SU[2];
      fabs(An) > 1.E-7 ? An=1./An : An = 0; // 1/A_normal
      double norm;
      for(int i=7; i<ND; i+=7) {
        M1x3& dR = *((M1x3*) &P[i]);
        M1x3& dA = *((M1x3*) &P[i+3]);
        norm = (dR[0]*SU[0] + dR[1]*SU[1] + dR[2]*SU[2])*An;	// dR_normal / A_normal
        dR[0] -= norm*A [0];   dR[1] -= norm*A [1];   dR[2] -= norm*A [2];
        dA[0] -= norm*SA[0];   dA[1] -= norm*SA[1];   dA[2] -= norm*SA[2];
      }
    }

    coveredDistance += S;
    Way  += fabs(S);
  } // end of linear extrapolation to surface

  return(true);
}


double RKTrackRep::estimateStep(std::vector<GFPointPath>& points,
                                const TVector3& pos,
                                const TVector3& dir,
                                const M1x4& SU,
                                const TVector3& MagField,
                                const GFDetPlane& plane,
                                const double& momentum,
                                double& relMomLoss,
                                double& deltaAngle,
                                bool& momLossExceeded,
                                bool& atPlane) const {

  static const double Smax      = 10.;          // max. step allowed [cm]
  static const double dAngleMax = 0.05;         // max. deviation of angle between direction before and after the step [rad]
  double Step;

  #ifdef DEBUG
    std::cout << " RKTrackRep::estimateStep \n";
    std::cout << "  position: "; pos.Print();
    std::cout << "  direction: "; dir.Print();
  #endif


  // calculate distance to surface
  double Dist = SU[3] - (pos.X()*SU[0] + pos.Y()*SU[1] + pos.Z()*SU[2]);  // Distance between start coordinates and surface
  double An = dir.X()*SU[0] + dir.Y()*SU[1] + dir.Z()*SU[2];              // An = dir * N;  component of dir normal to surface

  if (fabs(An) > 1.E-10) Step = Dist/An;
  else {
    Step = Dist*1.E10;
    if (An<0) Step *= -1.;
  }

  // see if dir points towards surface (1) or not (-1)
  double StepSign(1);
  if (Step<0) StepSign = -1;

  #ifdef DEBUG
    std::cout << "  Distance to plane: " << Dist << "\n";
    std::cout << "  guess for Step: " << Step << "\n";
    if (StepSign>0) std::cout << "  Direction is  pointing towards surface.\n";
    else  std::cout << "  Direction is pointing away from surface.\n";
  #endif

  // calculate way SmaxAngle after which momentum angle has changed AngleMax
  double Hmag(MagField.Mag()), SmaxAngle(Smax), radius(0), p_perp(0);
  if (Hmag > 1E-5){
    p_perp = ( dir - MagField*((dir*MagField)/(Hmag*Hmag)) ).Mag() * momentum; // [GeV]
    radius = p_perp/(0.3E-3*Hmag); // [cm]
    double sinAngle = fabs(sin(dir.Angle(MagField)));
    if (sinAngle > 1E-10) SmaxAngle = fabs(dAngleMax * radius / sinAngle); // [cm]
  }


  //
  // Select direction
  //
  // auto select
  if (fDirection == 0){
    #ifdef DEBUG
      std::cout << "  auto select direction. \n";
    #endif
  }
  // see if straight line approximation is ok
  else if ( fabs(Step) < 0.2*SmaxAngle ){
    #ifdef DEBUG
      std::cout << "  straight line approximation is fine. Delta angle until surface is reached is approx " << Step/SmaxAngle * dAngleMax * 180 / TMath::Pi()  << " deg \n";
    #endif

    // if direction is pointing to active part of surface
    if( plane.inActive(pos, dir) ) {
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
    Step = fabs(Step)*fDirection;
    #ifdef DEBUG
      std::cout << "  select direction according to fDirection. \n";
    #endif
  }

  #ifdef DEBUG
    std::cout << "  guess for Step (signed): " << Step << "\n";
  #endif

  // re-check sign of Step
  if (Step>=0) StepSign = 1;
  else StepSign = -1;

  //
  // Limit stepsize
  //
  bool stopBecauseOfCurvature(false);

  // reduce maximum stepsize Step to Smax
  if (fabs(Step) > Smax) Step = StepSign*Smax;

  // also limit stepsize according to the change of the momentum direction!
  if (fabs(Step) > SmaxAngle) {
    Step = StepSign*SmaxAngle;
    stopBecauseOfCurvature = true;
  }

  #ifdef DEBUG
    std::cout << "  limit from maxangle: " << SmaxAngle << ", radius: " << radius << "\n";
  #endif


  // call stepper and reduce stepsize if step not too small
  if (!fNoMaterial){

    if(fabs(Step) > MINSTEP){ // only call stepper if step estimation big enough
      double StepMat = GFMaterialEffects::getInstance()->stepper(fabs(Step),
                                                                 SmaxAngle,
                                                                 pos.X(), pos.Y(), pos.Z(),
                                                                 StepSign*dir.X(), StepSign*dir.Y(), StepSign*dir.Z(),
                                                                 momentum,
                                                                 relMomLoss,
                                                                 fPdg);
      if (fabs(Step) > StepMat) {
        Step = StepSign*StepMat;
        momLossExceeded = true;
      }

      #ifdef DEBUG
        std::cout << "  limit from stepper: " << Step << "\n";
      #endif
    }
  }
  

  if (!stopBecauseOfCurvature && !momLossExceeded){
    atPlane = true;
    
    // improve step estimation to surface according to curvature
    if (Hmag > 1E-5 && fabs(Step) > 0.1*SmaxAngle){
      //
      // simplified Runge Kutta Extrapolation
      //
      double S3 = Step/3.;
      double PS2 = fCharge/momentum*0.000149896229 * Step;
      M1x3   H0 = {0.,0.,0.};
      double   A0(0), A2(0), A3(0), A4(0), A5(0), A6(0);
      double   B0(0), B2(0), B3(0), B4(0), B5(0), B6(0);
      double   C0(0), C2(0), C3(0), C4(0), C5(0), C6(0);

      // First point
      H0[0] = PS2*MagField.X(); H0[1] = PS2*MagField.Y(); H0[2] = PS2*MagField.Z();     // H0 is PS2*(Hx, Hy, Hz) @ R0
      A0 = dir.Y()*H0[2]-dir.Z()*H0[1]; B0 = dir.Z()*H0[0]-dir.X()*H0[2]; C0 = dir.X()*H0[1]-dir.Y()*H0[0]; // (ax, ay, az) x H0
      A2 = dir.X()+A0                 ; B2 = dir.Y()+B0                 ; C2 = dir.Z()+C0                 ; // (A0, B0, C0) + (ax, ay, az)

      // Second point
      A3 = B2*H0[2]-C2*H0[1]+dir.X(); B3 = C2*H0[0]-A2*H0[2]+dir.Y(); C3 = A2*H0[1]-B2*H0[0]+dir.Z(); // (A2, B2, C2) x H0 + (ax, ay, az)
      A4 = B3*H0[2]-C3*H0[1]+dir.X(); B4 = C3*H0[0]-A3*H0[2]+dir.Y(); C4 = A3*H0[1]-B3*H0[0]+dir.Z(); // (A3, B3, C3) x H0 + (ax, ay, az)
      A5 = A4-dir.X()+A4            ; B5 = B4-dir.Y()+B4            ; C5 = C4-dir.Z()+C4            ; //    2*(A4, B4, C4) - (ax, ay, az)

      // Last point
      A6 = B5*H0[2]-C5*H0[1]; B6 = C5*H0[0]-A5*H0[2]; C6 = A5*H0[1]-B5*H0[0]; // (A5, B5, C5) x H0

      // calculate distance to surface
      Dist = SU[3] - ((pos.X()+(A2+A3+A4)*S3) * SU[0] +
                      (pos.Y()+(B2+B3+B4)*S3) * SU[1] +
                      (pos.Z()+(C2+C3+C4)*S3) * SU[2]);        // Distance between start coordinates and surface

      An = (A0+A3+A3+A5+A6)/3. * SU[0] +
           (B0+B3+B3+B5+B6)/3. * SU[1] +
           (C0+C3+C3+C5+C6)/3. * SU[2];    // An = dir * N;  component of dir normal to surface

      Step += Dist/An;

      #ifdef DEBUG
        std::cout << "  Improved step estimation taking curvature into account: " << Step << "\n";
      #endif
    }

  }

  deltaAngle += Step/SmaxAngle * dAngleMax;
  points.push_back(GFPointPath(pos, Step));

  #ifdef DEBUG
    std::cout << "  --> Step to be used: " << Step << "\n";
  #endif

  return Step;
}


void RKTrackRep::setPropDir(int dir){
  // make sure fDirection is -1, 0 or 1
  if (dir>0) fDirection = 1;
  else if (dir<0) fDirection = -1;
  else fDirection = 0;
}



ClassImp(RKTrackRep)
