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

#include <TDatabasePDG.h>
#include <GFException.h>
#include <GFFieldManager.h>
#include <GFMaterialEffects.h>

//root stuff
#include <TMath.h>

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

      fAuxInfo.ResizeTo(1,2);
      fAuxInfo(0,0) = fCacheSpu;
      fAuxInfo(0,1) = fDirection;

   } else {
      R__b.WriteClassBuffer(RKTrackRep::Class(),this);
   }
}


RKTrackRep::RKTrackRep() : GFAbsTrackRep(5), fDirection(0), fNoMaterial(false), fPdg(0), fCharge(0), fSpu(1), fCachePlane(), fCacheSpu(1), fAuxInfo(1,2) {
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


RKTrackRep::RKTrackRep(const GFTrackCand* const aGFTrackCandPtr, int pdgCode) :
                       GFAbsTrackRep(5), fDirection(0), fNoMaterial(false), fCachePlane(), fCacheSpu(1), fAuxInfo(1,2) {

	if (pdgCode == 0){
		pdgCode = aGFTrackCandPtr->getPdgCode();
	}
  setPDG(pdgCode); // also sets charge and mass

  initArrays();

  TMatrixDSym cov6D = aGFTrackCandPtr->getCovSeed();
  setPosMomCov(aGFTrackCandPtr->getPosSeed(),
               aGFTrackCandPtr->getMomSeed(),
               cov6D);

  if( cov6D[0][0] < 0.0 ){ // no valid cov was set in the trackCand so just set a large one
    fCov.Zero();
    static const double value(1.E4);
    fCov(0,0) = value;
    fCov(1,1) = value;
    fCov(2,2) = value;
    fCov(3,3) = value;
    fCov(4,4) = value;
  }
}

RKTrackRep::RKTrackRep(const TVector3& pos, const TVector3& mom, const TMatrixDSym cov, const int& pdgCode) :
                       GFAbsTrackRep(5), fDirection(0), fNoMaterial(false), fCachePlane(), fCacheSpu(1), fAuxInfo(1,2) {

  setPDG(pdgCode); // also sets charge and mass
  initArrays();
  setPosMomCov(pos, mom, cov);
}


void RKTrackRep::initArrays(){
  memset(fNoise,0x00,7*7*sizeof(double));
  memset(fOldCov,0x00,7*7*sizeof(double));

  memset(fJ_pM_5x7,0x00,5*7*sizeof(double));
  memset(fJ_pM_5x6,0x00,5*6*sizeof(double));
  memset(fJ_Mp_7x5,0x00,7*5*sizeof(double));
  memset(fJ_Mp_6x5,0x00,6*5*sizeof(double));
}


void RKTrackRep::setData(const TVectorD& st, const GFDetPlane& pl, const TMatrixDSym* cov, const TMatrixD* aux){
  if(aux != NULL) {
    fCacheSpu = (*aux)(0,0);
    fDirection = (*aux)(0,1);
  }
  else {
    if(pl!=fCachePlane){
      GFException exc("RKTrackRep::setData() was called with a reference plane which is not the same as the one from the last extrapolate(plane,state,cov).",__LINE__,__FILE__);
      throw exc;
    }
  }
  GFAbsTrackRep::setData(st,pl,cov);
  if (fCharge*fState(0) < 0) fCharge *= -1; // set charge accordingly! (fState[0] = q/p)
  fSpu = fCacheSpu;
}


const TMatrixD* RKTrackRep::getAuxInfo(const GFDetPlane& pl) {

  if(pl!=fCachePlane) {
    GFException exc("RKTrackRep::getAuxInfo() - Trying to get auxiliary information with planes mismatch (Information returned does not belong to requested plane)!",__LINE__,__FILE__);
	  throw exc;
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

  const TVector3& U(fRefPlane.getU());
  const TVector3& V(fRefPlane.getV());
  TVector3 W(fRefPlane.getNormal());


  fCov(0,0) = fCharge*fCharge/pow(mom.Mag(),6.) *
              (mom.X()*mom.X() * momerr.X()*momerr.X()+
               mom.Y()*mom.Y() * momerr.Y()*momerr.Y()+
               mom.Z()*mom.Z() * momerr.Z()*momerr.Z());

  fCov(1,1) = pow((U.X()/pw - W.X()*pu/(pw*pw)),2.) * momerr.X()*momerr.X() +
              pow((U.Y()/pw - W.Y()*pu/(pw*pw)),2.) * momerr.Y()*momerr.Y() +
              pow((U.Z()/pw - W.Z()*pu/(pw*pw)),2.) * momerr.Z()*momerr.Z();

  fCov(2,2) = pow((V.X()/pw - W.X()*pv/(pw*pw)),2.) * momerr.X()*momerr.X() +
              pow((V.Y()/pw - W.Y()*pv/(pw*pw)),2.) * momerr.Y()*momerr.Y() +
              pow((V.Z()/pw - W.Z()*pv/(pw*pw)),2.) * momerr.Z()*momerr.Z();

  fCov(3,3) = poserr.X()*poserr.X() * U.X()*U.X() +
              poserr.Y()*poserr.Y() * U.Y()*U.Y() +
              poserr.Z()*poserr.Z() * U.Z()*U.Z();

  fCov(4,4) = poserr.X()*poserr.X() * V.X()*V.X() +
              poserr.Y()*poserr.Y() * V.Y()*V.Y() +
              poserr.Z()*poserr.Z() * V.Z()*V.Z();
}


void RKTrackRep::calcState(const TVector3& pos,
                           const TVector3& mom){

  fRefPlane.setON(pos, mom);
  fSpu=1.;

  fState(0) = fCharge/mom.Mag();

  //u' and v'
  fState(1) = 0.;
  fState(2) = 0.;

  //u and v
  fState(3) = 0.;
  fState(4) = 0.;
}



void RKTrackRep::getState7(M1x7& state7) {
  getState7(state7, fState, fRefPlane, fSpu);
}


void RKTrackRep::getState7(M1x7& state7, const TVectorD& state5, const GFDetPlane& pl, const double& spu) {

  const TVector3& U(pl.getU());
  const TVector3& V(pl.getV());
  const TVector3& O(pl.getO());
  TVector3 W(pl.getNormal());

  state7[0] = O.X() + state5(3)*U.X() + state5(4)*V.X(); // x
  state7[1] = O.Y() + state5(3)*U.Y() + state5(4)*V.Y(); // y
  state7[2] = O.Z() + state5(3)*U.Z() + state5(4)*V.Z(); // z

  state7[3] = spu * (W.X() + state5(1)*U.X() + state5(2)*V.X()); // a_x
  state7[4] = spu * (W.Y() + state5(1)*U.Y() + state5(2)*V.Y()); // a_y
  state7[5] = spu * (W.Z() + state5(1)*U.Z() + state5(2)*V.Z()); // a_z

  // normalize dir
  double norm = 1. / sqrt(state7[3]*state7[3] + state7[4]*state7[4] + state7[5]*state7[5]);
  for (unsigned int i=3; i<6; ++i) state7[i] *= norm;

  state7[6] = state5(0); // q/p
}


TVectorD RKTrackRep::getState5(const M1x7& state7, const GFDetPlane& pl, double& spu) {

  const TVector3& U(pl.getU());
  const TVector3& V(pl.getV());

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

  TVectorD state5(5);
  state5(0) = state7[6];
  state5(1) = fDir*U / AtW;
  state5(2) = fDir*V / AtW;
  state5(3) = fPos*U;
  state5(4) = fPos*V;

  return state5;
}



void RKTrackRep::transformPM7(const TMatrixD& in5x5, M7x7& out7x7,
                              const GFDetPlane& pl, const TVectorD& state5, const double&  spu,
                              TMatrixD* Jac) {

  // get vectors and aux variables
  const TVector3& U(pl.getU());
  const TVector3& V(pl.getV());
  TVector3 W(pl.getNormal());

  fpTilde.SetXYZ(spu * (W.X() + state5(1)*U.X() + state5(2)*V.X()), // a_x
                 spu * (W.Y() + state5(1)*U.Y() + state5(2)*V.Y()), // a_y
                 spu * (W.Z() + state5(1)*U.Z() + state5(2)*V.Z()));// a_z


  const double pTildeMag = fpTilde.Mag();
  const double pTildeMag2 = pTildeMag*pTildeMag;

  const double utpTildeOverpTildeMag2 = U*fpTilde / pTildeMag2;
  const double vtpTildeOverpTildeMag2 = V*fpTilde / pTildeMag2;

  //J_pM matrix is d(x,y,z,ax,ay,az,q/p) / d(q/p,u',v',u,v)   (out is 7x7)

   // d(x,y,z)/d(u)
  fJ_pM_5x7[21] = U.X(); // [3][0]
  fJ_pM_5x7[22] = U.Y(); // [3][1]
  fJ_pM_5x7[23] = U.Z(); // [3][2]
  // d(x,y,z)/d(v)
  fJ_pM_5x7[28] = V.X(); // [4][2]
  fJ_pM_5x7[29] = V.Y(); // [4][2]
  fJ_pM_5x7[30] = V.Z(); // [4][2]
  // d(q/p)/d(q/p)
  fJ_pM_5x7[6] = 1.; // not needed for array matrix multiplication
  // d(ax,ay,az)/d(u')
  double fact = spu / pTildeMag;
  fJ_pM_5x7[10] = fact * ( U.X() - fpTilde.X()*utpTildeOverpTildeMag2 ); // [1][3]
  fJ_pM_5x7[11] = fact * ( U.Y() - fpTilde.Y()*utpTildeOverpTildeMag2 ); // [1][4]
  fJ_pM_5x7[12] = fact * ( U.Z() - fpTilde.Z()*utpTildeOverpTildeMag2 ); // [1][5]
  // d(ax,ay,az)/d(v')
  fJ_pM_5x7[17] = fact * ( V.X() - fpTilde.X()*vtpTildeOverpTildeMag2 ); // [2][3]
  fJ_pM_5x7[18] = fact * ( V.Y() - fpTilde.Y()*vtpTildeOverpTildeMag2 ); // [2][4]
  fJ_pM_5x7[19] = fact * ( V.Z() - fpTilde.Z()*vtpTildeOverpTildeMag2 ); // [2][5]


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


void RKTrackRep::transformPM6(const TMatrixDSym& in5x5, M6x6& out6x6,
                              const GFDetPlane& pl, const TVectorD& state5, const double&  spu,
                              TMatrixD* Jac) {

  // get vectors and aux variables
  const TVector3& U(pl.getU());
  const TVector3& V(pl.getV());
  TVector3 W(pl.getNormal());

  fpTilde.SetXYZ(spu * (W.X() + state5(1)*U.X() + state5(2)*V.X()), // a_x
                 spu * (W.Y() + state5(1)*U.Y() + state5(2)*V.Y()), // a_y
                 spu * (W.Z() + state5(1)*U.Z() + state5(2)*V.Z()));// a_z

  const double pTildeMag = fpTilde.Mag();
  const double pTildeMag2 = pTildeMag*pTildeMag;

  const double utpTildeOverpTildeMag2 = U*fpTilde / pTildeMag2;
  const double vtpTildeOverpTildeMag2 = V*fpTilde / pTildeMag2;

  //J_pM matrix is d(x,y,z,px,py,pz) / d(q/p,u',v',u,v)       (out is 6x6)

  const double qop = state5(0);
  const double p = fCharge/qop; // momentum

  // d(px,py,pz)/d(q/p)
  double fact = -1. * p / (pTildeMag * qop);
  fJ_pM_5x6[3] = fact * fpTilde.X(); // [0][3]
  fJ_pM_5x6[4] = fact * fpTilde.Y(); // [0][4]
  fJ_pM_5x6[5] = fact * fpTilde.Z(); // [0][5]
  // d(px,py,pz)/d(u')
  fact = p * spu / pTildeMag;
  fJ_pM_5x6[9]  = fact * ( U.X() - fpTilde.X()*utpTildeOverpTildeMag2 ); // [1][3]
  fJ_pM_5x6[10] = fact * ( U.Y() - fpTilde.Y()*utpTildeOverpTildeMag2 ); // [1][4]
  fJ_pM_5x6[11] = fact * ( U.Z() - fpTilde.Z()*utpTildeOverpTildeMag2 ); // [1][5]
  // d(px,py,pz)/d(v')
  fJ_pM_5x6[15] = fact * ( V.X() - fpTilde.X()*vtpTildeOverpTildeMag2 ); // [2][3]
  fJ_pM_5x6[16] = fact * ( V.Y() - fpTilde.Y()*vtpTildeOverpTildeMag2 ); // [2][4]
  fJ_pM_5x6[17] = fact * ( V.Z() - fpTilde.Z()*vtpTildeOverpTildeMag2 ); // [2][5]
  // d(x,y,z)/d(u)
  fJ_pM_5x6[18] = U.X(); // [3][0]
  fJ_pM_5x6[19] = U.Y(); // [3][1]
  fJ_pM_5x6[20] = U.Z(); // [3][2]
  // d(x,y,z)/d(v)
  fJ_pM_5x6[24] = V.X(); // [4][0]
  fJ_pM_5x6[25] = V.Y(); // [4][1]
  fJ_pM_5x6[26] = V.Z(); // [4][2]


  // do the transformation
  // out = J_pM^T * in5x5 * J_pM
  const M5x5& in5x5_ = *((M5x5*) in5x5.GetMatrixArray());
  RKTools::J_pMTxcov5xJ_pM(fJ_pM_5x6, in5x5_, out6x6);

  if (Jac!=NULL){
    Jac->ResizeTo(5,6);
    *Jac = TMatrixD(5,6, &(fJ_pM_5x6[0]));
  }
}


void RKTrackRep::transformM7P(const M7x7& in7x7, TMatrixDSym& out5x5,
                              const GFDetPlane& pl, const M1x7& state7,
                              TMatrixD* Jac) {

  out5x5.ResizeTo(5, 5);

  // get vectors and aux variables
  const TVector3& U(pl.getU());
  const TVector3& V(pl.getV());
  TVector3 W(pl.getNormal());

  fDir.SetXYZ(state7[3], state7[4], state7[5]);

  const double AtU = fDir*U;
  const double AtV = fDir*V;
  const double AtW = fDir*W;

  // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,ax,ay,az,q/p)   (in is 7x7)

  // d(u')/d(ax,ay,az)
  double fact = 1./(AtW*AtW);
  fJ_Mp_7x5[16] = fact * (U.X()*AtW - W.X()*AtU); // [3][1]
  fJ_Mp_7x5[21] = fact * (U.Y()*AtW - W.Y()*AtU); // [4][1]
  fJ_Mp_7x5[26] = fact * (U.Z()*AtW - W.Z()*AtU); // [5][1]
  // d(v')/d(ax,ay,az)
  fJ_Mp_7x5[17] = fact * (V.X()*AtW - W.X()*AtV); // [3][2]
  fJ_Mp_7x5[22] = fact * (V.Y()*AtW - W.Y()*AtV); // [4][2]
  fJ_Mp_7x5[27] = fact * (V.Z()*AtW - W.Z()*AtV); // [5][2]
  // d(q/p)/d(q/p)
  fJ_Mp_7x5[30] = 1.; // [6][0]  - not needed for array matrix multiplication
  //d(u)/d(x,y,z)
  fJ_Mp_7x5[3]  = U.X(); // [0][3]
  fJ_Mp_7x5[8]  = U.Y(); // [1][3]
  fJ_Mp_7x5[13] = U.Z(); // [2][3]
  //d(v)/d(x,y,z)
  fJ_Mp_7x5[4]  = V.X(); // [0][4]
  fJ_Mp_7x5[9]  = V.Y(); // [1][4]
  fJ_Mp_7x5[14] = V.Z(); // [2][4]


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


void RKTrackRep::transformM6P(const M6x6& in6x6, TMatrixDSym& out5x5,
                              const GFDetPlane& pl, const M1x7& state7,
                              TMatrixD* Jac) {

  out5x5.ResizeTo(5, 5);

  // get vectors and aux variables
  const TVector3& U(pl.getU());
  const TVector3& V(pl.getV());
  TVector3 W(pl.getNormal());

  fDir.SetXYZ(state7[3], state7[4], state7[5]);

  const double AtU = fDir*U;
  const double AtV = fDir*V;
  const double AtW = fDir*W;

  // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,px,py,pz)       (in is 6x6)

  const double qop = state7[6];
  const double p = fCharge/qop; // momentum

  //d(u)/d(x,y,z)
  fJ_Mp_6x5[3]  = U.X(); // [0][3]
  fJ_Mp_6x5[8]  = U.Y(); // [1][3]
  fJ_Mp_6x5[13] = U.Z(); // [2][3]
  //d(v)/d(x,y,z)
  fJ_Mp_6x5[4]  = V.X(); // [0][4]
  fJ_Mp_6x5[9]  = V.Y(); // [1][4]
  fJ_Mp_6x5[14] = V.Z(); // [2][4]
  // d(q/p)/d(px,py,pz)
  double fact = (-1.) * qop / p;
  fJ_Mp_6x5[15] = fact * fDir.X(); // [3][0]
  fJ_Mp_6x5[20] = fact * fDir.Y(); // [4][0]
  fJ_Mp_6x5[25] = fact * fDir.Z(); // [5][0]
  // d(u')/d(px,py,pz)
  fact = 1./(p*AtW*AtW);
  fJ_Mp_6x5[16] = fact * (U.X()*AtW - W.X()*AtU); // [3][1]
  fJ_Mp_6x5[21] = fact * (U.Y()*AtW - W.Y()*AtU); // [4][1]
  fJ_Mp_6x5[26] = fact * (U.Z()*AtW - W.Z()*AtU); // [5][1]
  // d(v')/d(px,py,pz)
  fJ_Mp_6x5[17] = fact * (V.X()*AtW - W.X()*AtV); // [3][2]
  fJ_Mp_6x5[22] = fact * (V.Y()*AtW - W.Y()*AtV); // [4][2]
  fJ_Mp_6x5[27] = fact * (V.Z()*AtW - W.Z()*AtV); // [5][2]

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


void RKTrackRep::getPosMomCov(const GFDetPlane& pl, TVector3& pos, TVector3& mom, TMatrixDSym& cov6x6){
  TVectorD statePred(fState);
  TMatrixDSym covPred(fCov);
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


void RKTrackRep::setPosMomCov(const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov6x6){

  if (cov6x6.GetNcols()!=6 || cov6x6.GetNrows()!=6){
    GFException exc("RKTrackRep::setPosMomCov ==> cov has to be 6x6 (x, y, z, px, py, pz)",__LINE__,__FILE__);
    throw exc;
  }

  if (fCharge == 0){
    GFException exc("RKTrackRep::setPosMomCov ==> charge is 0. setPosMomCov cannot work with fCharge == 0 ",__LINE__,__FILE__);
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



double RKTrackRep::extrapolateToPoint(const TVector3& pos,
                                      TVector3& poca,
                                      TVector3& dirInPoca){

#ifdef DEBUG
  std::cout << "RKTrackRep::extrapolateToPoint()\n";
#endif

  static const unsigned int maxIt(1000);

  M1x7 state7;
  getState7(state7);
  fDir.SetXYZ(state7[3], state7[4], state7[5]);

  double step(0.), lastStep(0.), maxStep(1.E99), angle(0), distToPoca(0), tracklength(0);
  TVector3 lastDir(0,0,0);

  GFDetPlane pl(pos, fDir);
  unsigned int iterations(0);

  while(true){
    lastStep = step;
    lastDir = fDir;

    step =  this->Extrap(pl, state7, NULL, true, maxStep);
    tracklength += step;

    fDir.SetXYZ(state7[3], state7[4], state7[5]);
    poca.SetXYZ(state7[0], state7[1], state7[2]);

    // check break conditions
    angle = fabs(fDir.Angle((pos-poca))-TMath::PiOver2()); // angle between direction and connection to point - 90 deg
    distToPoca = (pos-poca).Mag();
    if (angle*distToPoca < 0.1*MINSTEP) break;
    if(++iterations == maxIt) {
      GFException exc("RKTrackRep::extrapolateToPoint ==> extrapolation to point failed, maximum number of iterations reached",__LINE__,__FILE__);
      throw exc;
    }

    // if lastStep and step have opposite sign, the real normal vector lies somewhere between the last two normal vectors (i.e. the directions)
    // -> try mean value of the two (normalization not needed)
    if (lastStep*step < 0){
      fDir += lastDir;
      maxStep = 0.5*fabs(lastStep); // make it converge!
    }

    pl.setNormal(fDir);
  }

  dirInPoca.SetXYZ(state7[3], state7[4], state7[5]);

#ifdef DEBUG
  std::cout << "RKTrackRep::extrapolateToPoint(): Reached POCA after " << iterations+1 << " iterations. Distance: " << (pos-poca).Mag() << " cm. Angle deviation: " << dirInPoca.Angle((pos-poca))-TMath::PiOver2() << " rad \n";
#endif

  return tracklength;
}


TVector3 RKTrackRep::poca2Line(const TVector3& extr1,const TVector3& extr2,const TVector3& point) const {
  
  TVector3 pocaOnLine(extr2);
  pocaOnLine -= extr1; // wireDir

  if(pocaOnLine.Mag()<1.E-8){
    GFException exc("RKTrackRep::poca2Line ==> try to find POCA between line and point, but the line is really just a point",__LINE__,__FILE__);
    throw exc;
  }

  double t = 1./(pocaOnLine.Mag2()) * ((point*pocaOnLine) + extr1.Mag2() - (extr1*extr2));
  pocaOnLine *= t;
  pocaOnLine += extr1;
  return pocaOnLine; // = extr1 + t*wireDir
}


double RKTrackRep::extrapolateToLine(const TVector3& point1,
                                     const TVector3& point2,
                                     TVector3& poca,
                                     TVector3& dirInPoca,
                                     TVector3& poca_onwire){

#ifdef DEBUG
  std::cout << "RKTrackRep::extrapolateToLine(), (x,y) = (" << point1.X() << ", " << point1.Y() << ")\n";
#endif

  static const unsigned int maxIt(1000);

  M1x7 state7;
  getState7(state7);
  fDir.SetXYZ(state7[3], state7[4], state7[5]);

  double step(0.), lastStep(0.), maxStep(1.E99), angle(0), distToPoca(0), tracklength(0);
  TVector3 wireDir(point2);
  wireDir -= point1;
  wireDir.SetMag(1.);
  TVector3 lastDir(0,0,0);

  GFDetPlane pl(point1, fDir.Cross(wireDir), wireDir);
  unsigned int iterations(0);

  while(true){
    lastStep = step;
    lastDir = fDir;

    step = this->Extrap(pl, state7, NULL, true, maxStep);
    tracklength += step;

    fDir.SetXYZ(state7[3], state7[4], state7[5]);
    poca.SetXYZ(state7[0], state7[1], state7[2]);
    poca_onwire = poca2Line(point1, point2, poca);

    // check break conditions
    angle = fabs(fDir.Angle((poca_onwire-poca))-TMath::PiOver2()); // angle between direction and connection to point - 90 deg
    distToPoca = (poca_onwire-poca).Mag();
    if (angle*distToPoca < 0.1*MINSTEP) break;
    if(++iterations == maxIt) {
      GFException exc("RKTrackRep::extrapolateToLine ==> extrapolation to line failed, maximum number of iterations reached",__LINE__,__FILE__);
      throw exc;
    }

    // if lastStep and step have opposite sign, the real normal vector lies somewhere between the last two normal vectors (i.e. the directions)
    // -> try mean value of the two (normalization not needed)
    if (lastStep*step < 0){
      fDir += lastDir;
      maxStep = 0.5*fabs(lastStep); // make it converge!
    }

    pl.setU(fDir.Cross(wireDir));
  }

  dirInPoca.SetXYZ(state7[3], state7[4], state7[5]);

#ifdef DEBUG
  std::cout << "RKTrackRep::extrapolateToLine(): Reached POCA after " << iterations+1 << " iterations. Distance: " << (poca_onwire-poca).Mag() << " cm. Angle deviation: " << dirInPoca.Angle((poca_onwire-poca))-TMath::PiOver2() << " rad \n";
#endif

  return tracklength;
}


double RKTrackRep::extrapolate(const GFDetPlane& pl, 
                               TVectorD& statePred,
                               TMatrixDSym& covPred){
  
#ifdef DEBUG
  std::cout << "RKTrackRep::extrapolate(pl, statePred, covPred)\n";
#endif

  M1x7 state7;
  getState7(state7);
  M7x7 cov7x7;

  transformPM7(fCov, cov7x7, fRefPlane, fState, fSpu);

  double coveredDistance = Extrap(pl, state7, &cov7x7);
  
  statePred.ResizeTo(5);
  statePred = getState5(state7, pl, fCacheSpu);
  fCachePlane = pl;

  covPred.ResizeTo(5, 5);
  transformM7P(cov7x7, covPred, pl, state7);

  return coveredDistance;
}


double RKTrackRep::extrapolate(const GFDetPlane& pl, 
                               TVectorD& statePred){

#ifdef DEBUG
  std::cout << "RKTrackRep::extrapolate(pl, statePred)\n";
#endif

  M1x7 state7;
  getState7(state7);
  double coveredDistance = Extrap(pl, state7);
  double spu;
  statePred.ResizeTo(5);
  statePred = getState5(state7, pl, spu);

  return coveredDistance;
}


double RKTrackRep::stepalong(double h, TVector3& pos, TVector3& dir){

#ifdef DEBUG
  std::cout << "RKTrackRep::stepalong()\n";
#endif

  TVector3 dest;

  static const unsigned int maxIt(1000);
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
    coveredDistance += this->Extrap(pl, state7, NULL, true);

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



double RKTrackRep::extrapolateToCylinder(double radius, TVector3& pos, TVector3& mom){

#ifdef DEBUG
  std::cout << "RKTrackRep::extrapolateToCylinder()\n";
#endif

  TVector3 dest, normal;

  static const unsigned int maxIt(1000);
  double coveredDistance(0.);

  M1x7 state7;
  getState7(state7);

  GFDetPlane pl;
  unsigned int iterations(0);

  while(true){
    pos.SetXYZ(state7[0], state7[1], state7[2]);
    mom.SetXYZ(state7[3], state7[4], state7[5]);
    mom.SetMag(1.);

    // solve quadratic equation
    double a = pow(mom.X(), 2) + pow(mom.Y(), 2);
    double b = 2. * (pos.X()*mom.X() + pos.Y()*mom.Y());
    double c = pow(pos.X(), 2) + pow(pos.Y(), 2) - radius*radius;
    double term = sqrt(b*b - 4.*a*c);
    double k1 = (-1.*b + term)/(2.*a);
    double k2 = (-1.*b - term)/(2.*a);

    // select smallest absolute solution -> closest cylinder surface
    double k = k1;
    if (fabs(k2)<fabs(k))
      k = k2;


    dest = pos + k * mom;
    normal = dest; // normal vector on cylinder surface
    normal.SetZ(0);
    normal.SetMag(1.);

    pl.setON(dest, normal);
    coveredDistance += this->Extrap(pl, state7, NULL, true);

    if(fabs(k)<MINSTEP) break;

    if(++iterations == maxIt) {
      GFException exc("RKTrackRep::extrapolateToCylinder ==> maximum number of iterations reached",__LINE__,__FILE__);
      throw exc;
    }
  }

  pos.SetXYZ(state7[0], state7[1], state7[2]);
  mom.SetXYZ(state7[3], state7[4], state7[5]);
  mom.SetMag(fCharge/state7[6]);

  return coveredDistance;

}



double RKTrackRep::Extrap( const GFDetPlane& plane, M1x7& state7, M7x7* cov, bool onlyOneStep, double maxStep) {

  static const unsigned int maxNumIt(500);
  unsigned int numIt(0);

  const bool calcCov(cov!=NULL);
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

    // initialize cov with unit matrix
    if(calcCov){
      memcpy(fOldCov, cov, 7*7*sizeof(double));
      memset(cov,0x00,49*sizeof(double));
      for(int i=0; i<7; ++i) (*cov)[8*i] = 1.;
    }

    fDirectionBefore.SetXYZ(state7[3], state7[4], state7[5]); // direction before propagation

    // propagation
    std::vector<GFPointPath> points;
    points.reserve(50);

    bool checkJacProj = true;

    if( ! this->RKutta(plane, state7, cov, coveredDistance, points, checkJacProj, onlyOneStep, maxStep) ) {
      GFException exc("RKTrackRep::Extrap ==> Runge Kutta propagation failed",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    fPos.SetXYZ(state7[0], state7[1], state7[2]);
    if (!fNoMaterial) points.push_back(GFPointPath(fPos, 0)); // add last point

    #ifdef DEBUG
      std::cout<<"Original points \n";
      for (unsigned int i=0; i<points.size(); ++i){
        points[i].Print();
      }
      std::cout<<"\n";
    #endif

    fDirectionAfter.SetXYZ(state7[3], state7[4], state7[5]); // direction after propagation

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
                                                                 fabs(fCharge/state7[6]), // momentum
                                                                 fPdg,
                                                                 fXX0,
                                                                 fNoise,
                                                                 (double *)cov,
                                                                 &fDirectionBefore,
                                                                 &fDirectionAfter);

      #ifdef DEBUG
        std::cout << "momLoss: " << momLoss << " GeV; relative: " << momLoss/fabs(fCharge/state7[6]) << "\n";
      #endif

      // do momLoss only for defined 1/momentum .ne.0
      if(fabs(state7[6])>1.E-10) state7[6] = fCharge/(fabs(fCharge/state7[6])-momLoss);
    } // finished MatFX

    if(calcCov){ // propagate cov and add noise
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
      RKTools::J_MMTxcov7xJ_MM(*cov, fOldCov);
      memcpy(cov, fOldCov, 7*7*sizeof(double));

      // add noise to cov
      for (int i=0; i<7*7; ++i) (*cov)[i] += fNoise[i];

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
//    state7 - initial parameters (coordinates(cm), direction,
//             charge/momentum (Gev-1) 
//    cov      and derivatives this parameters  (7x7)            
//         
//    X        	Y        	Z        	Ax       	Ay       	Az       	q/P                   
//    state7[0] state7[1] state7[2] state7[3] state7[4] state7[5] state7[6]  
//
//    dX/dp    	dY/dp    	dZ/dp    	dAx/dp   	dAy/dp   	dAz/dp   	d(q/P)/dp
//    cov[ 0]   cov[ 1]   cov[ 2]   cov[ 3]   cov[ 4]   cov[ 5]   cov[ 6]   			      d()/dp1  
//
//    cov[ 7]   cov[ 8]   cov[ 9]   cov[10]   cov[11]   cov[12]   cov[13]   		      	d()/dp2        
//    ............................................................................		d()/dpND       
//                                                                                  
// Authors: R.Brun, M.Hansroul, V.Perevoztchikov (Geant3)                           
//  
bool RKTrackRep::RKutta (const GFDetPlane& plane,
                         M1x7& state7,
                         M7x7* cov,
                         double& coveredDistance,
                         std::vector<GFPointPath>& points,
                         bool& checkJacProj,
                         bool onlyOneStep,
                         double maxStep) {

  // limits, check-values, etc. Can be tuned!
  static const double Wmax   = 3000.;           // max. way allowed [cm]
  static const double AngleMax = 6.3;           // max. total angle change of momentum. Prevents extrapolating a curler round and round if no active plane is found.
  static const double Pmin   = 4.E-3;           // minimum momentum for propagation [GeV]
  static const unsigned int maxNumIt = 1000;    // maximum number of iterations in main loop
  // Aux parameters
  M1x3&   R           = *((M1x3*) &state7[0]);  // Start coordinates  [cm] 	(x,  y,  z)
  M1x3&   A           = *((M1x3*) &state7[3]);  // Start directions 	      (ax, ay, az); 	ax^2+ay^2+az^2=1
  M1x3    SA          = {0.,0.,0.};             // Start directions derivatives dA/S
  double  Way         = 0.;                     // Sum of absolute values of all extrapolation steps [cm]
  bool    atPlane = false;                      // stepper thinks that the plane will be reached in that step -> linear extrapolation and projection of jacobian
  bool    momLossExceeded = false;              // stepper had to limit stepsize due to momentum loss -> no next RKutta loop, no linear extrapolation
  fPos.SetXYZ(R[0],R[1],R[2]);                  // position
  fDir.SetXYZ(A[0],A[1],A[2]);                  // direction
  double   momentum   = fabs(fCharge/state7[6]);// momentum [GeV]
  double   relMomLoss = 0;                      // relative momentum loss in RKutta
  double   deltaAngle = 0.;                     // total angle by which the momentum has changed during extrapolation
  double   An(0), S(0), Sl(0), CBA(0);
  M1x4     SU = {0.,0.,0.,0.};


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
    sstream << "RKTrackRep::RKutta ==> momentum too low: " << momentum*1000. << " MeV";
    GFException exc(sstream.str(),__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }
  

  // make SU vector point away from origin
  const TVector3 W(plane.getNormal());
  if(W*plane.getO() > 0){SU[0] =     W.X();  SU[1] =     W.Y();  SU[2] =     W.Z();}
  else                  {SU[0] = -1.*W.X();  SU[1] = -1.*W.Y();  SU[2] = -1.*W.Z();  }
  SU[3] = plane.distance(0., 0., 0.);

  unsigned int counter(0);

  // Step estimation (signed)
  S = estimateStep(points, fPos, fDir, SU, plane, momentum, relMomLoss, deltaAngle, momLossExceeded, atPlane, maxStep);
  if (fabs(S) < 0.001*MINSTEP) {
    #ifdef DEBUG
      std::cout << " RKutta - step too small -> break \n";
    #endif
    ++counter; // skip the main loop, go to linear extrapolation step (will be skipped) and just project jacobian
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

    RKPropagate(state7, cov, SA, S); // the actual Runkge Kutta propagation
    fPos.SetXYZ(R[0],R[1],R[2]);
    fDir.SetXYZ(A[0],A[1],A[2]);
    
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
    S = estimateStep(points, fPos, fDir, SU, plane, momentum, relMomLoss, deltaAngle, momLossExceeded, atPlane, maxStep);

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

      R[0] += S*(A[0]-0.5*S*SA[0]);    // R = R + S*(A - 0.5*S*SA); approximation for final point on surface
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
    if(cov != NULL){
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
      for(int i=0; i<49; i+=7) {
        norm = ((*cov)[i]*SU[0] + (*cov)[i+1]*SU[1] + (*cov)[i+2]*SU[2])*An;	// dR_normal / A_normal
        (*cov)[i]   -= norm*A [0];   (*cov)[i+1] -= norm*A [1];   (*cov)[i+2] -= norm*A [2];
        (*cov)[i+3] -= norm*SA[0];   (*cov)[i+4] -= norm*SA[1];   (*cov)[i+5] -= norm*SA[2];
      }
    }

    coveredDistance += S;
    Way  += fabs(S);
  } // end of linear extrapolation to surface

  return(true);
}


void
RKTrackRep::RKPropagate(M1x7& state7,
                        M7x7* cov,
                        M1x3& SA,
                        double S,
                        bool varField) const {

  // important fixed numbers
  static const double EC     = 0.000149896229;  // c/(2*10^12) resp. c/2Tera
  static const double P3     = 1./3.;           // 1/3
  // Aux parameters
  M1x3&   R           = *((M1x3*) &state7[0]);       // Start coordinates  [cm]  (x,  y,  z)
  M1x3&   A           = *((M1x3*) &state7[3]);       // Start directions         (ax, ay, az);   ax^2+ay^2+az^2=1
  double  S3(0), S4(0), PS2(0);
  M1x3     H0 = {0.,0.,0.}, H1 = {0.,0.,0.}, H2 = {0.,0.,0.}, r = {0.,0.,0.};
  // Variables for RKutta main loop
  double   A0(0), A1(0), A2(0), A3(0), A4(0), A5(0), A6(0);
  double   B0(0), B1(0), B2(0), B3(0), B4(0), B5(0), B6(0);
  double   C0(0), C1(0), C2(0), C3(0), C4(0), C5(0), C6(0);

  bool calcCov(cov != NULL);

  //
  // Runge Kutta Extrapolation
  //
  S3 = P3*S;
  S4 = 0.25*S;
  PS2 = state7[6]*EC * S;

  // First point
  r[0] = R[0];           r[1] = R[1];           r[2]=R[2];
  TVector3 pos(r[0], r[1], r[2]);// vector of start coordinates R0  (x, y, z)
  TVector3 field(GFFieldManager::getFieldVal(pos));       // magnetic field in 10^-4 T = kGauss
  H0[0] = PS2*field.X(); H0[1] = PS2*field.Y(); H0[2] = PS2*field.Z();     // H0 is PS2*(Hx, Hy, Hz) @ R0
  A0 = A[1]*H0[2]-A[2]*H0[1]; B0 = A[2]*H0[0]-A[0]*H0[2]; C0 = A[0]*H0[1]-A[1]*H0[0]; // (ax, ay, az) x H0
  A2 = A[0]+A0              ; B2 = A[1]+B0              ; C2 = A[2]+C0              ; // (A0, B0, C0) + (ax, ay, az)
  if (varField) {
    A1 = A2+A[0]            ; B1 = B2+A[1]              ; C1 = C2+A[2]              ; // (A0, B0, C0) + 2*(ax, ay, az)
  }

  // Second point
  if (varField) {
    r[0] += A1*S4;         r[1] += B1*S4;         r[2] += C1*S4;
    pos.SetXYZ(r[0], r[1], r[2]);
    field = GFFieldManager::getFieldVal(pos);
    H1[0] = field.X()*PS2; H1[1] = field.Y()*PS2; H1[2] = field.Z()*PS2; // H1 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * [(A0, B0, C0) + 2*(ax, ay, az)]
  }
  else if (calcCov) memcpy(H1, H0, 3*sizeof(double));
  A3 = B2*H1[2]-C2*H1[1]+A[0]; B3 = C2*H1[0]-A2*H1[2]+A[1]; C3 = A2*H1[1]-B2*H1[0]+A[2]; // (A2, B2, C2) x H1 + (ax, ay, az)
  A4 = B3*H1[2]-C3*H1[1]+A[0]; B4 = C3*H1[0]-A3*H1[2]+A[1]; C4 = A3*H1[1]-B3*H1[0]+A[2]; // (A3, B3, C3) x H1 + (ax, ay, az)
  A5 = A4-A[0]+A4            ; B5 = B4-A[1]+B4            ; C5 = C4-A[2]+C4            ; //    2*(A4, B4, C4) - (ax, ay, az)

  // Last point
  if (varField) {
    r[0]=R[0]+S*A4;         r[1]=R[1]+S*B4;         r[2]=R[2]+S*C4;  //setup.Field(r,H2);
    pos.SetXYZ(r[0], r[1], r[2]);
    field = GFFieldManager::getFieldVal(pos);
    H2[0] = field.X()*PS2;  H2[1] = field.Y()*PS2;  H2[2] = field.Z()*PS2; // H2 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * (A4, B4, C4)
  }
  else if (calcCov) memcpy(H2, H0, 3*sizeof(double));
  A6 = B5*H2[2]-C5*H2[1]; B6 = C5*H2[0]-A5*H2[2]; C6 = A5*H2[1]-B5*H2[0]; // (A5, B5, C5) x H2

  //
  // Derivatives of track parameters
  //
  if(calcCov){
    double   dA0(0), dA2(0), dA3(0), dA4(0), dA5(0), dA6(0);
    double   dB0(0), dB2(0), dB3(0), dB4(0), dB5(0), dB6(0);
    double   dC0(0), dC2(0), dC3(0), dC4(0), dC5(0), dC6(0);

    // d(x, y, z)/d(x, y, z) submatrix is unit matrix
    (*cov)[0] = 1;  (*cov)[8] = 1;  (*cov)[16] = 1;
    // d(ax, ay, az)/d(ax, ay, az) submatrix is 0
    // start with d(x, y, z)/d(ax, ay, az)
    for(int i=3*7; i<49; i+=7) {

      if(i==42) {(*cov)[i+3]*=state7[6]; (*cov)[i+4]*=state7[6]; (*cov)[i+5]*=state7[6];}

      //first point
      dA0 = H0[2]*(*cov)[i+4]-H0[1]*(*cov)[i+5];    // dA0/dp }
      dB0 = H0[0]*(*cov)[i+5]-H0[2]*(*cov)[i+3];    // dB0/dp  } = dA x H0
      dC0 = H0[1]*(*cov)[i+3]-H0[0]*(*cov)[i+4];    // dC0/dp }

      if(i==42) {dA0+=A0; dB0+=B0; dC0+=C0;}     // if last row: (dA0, dB0, dC0) := (dA0, dB0, dC0) + (A0, B0, C0)

      dA2 = dA0+(*cov)[i+3];        // }
      dB2 = dB0+(*cov)[i+4];        //  } = (dA0, dB0, dC0) + dA
      dC2 = dC0+(*cov)[i+5];        // }

      //second point
      dA3 = (*cov)[i+3]+dB2*H1[2]-dC2*H1[1];    // dA3/dp }
      dB3 = (*cov)[i+4]+dC2*H1[0]-dA2*H1[2];    // dB3/dp  } = dA + (dA2, dB2, dC2) x H1
      dC3 = (*cov)[i+5]+dA2*H1[1]-dB2*H1[0];    // dC3/dp }

      if(i==42) {dA3+=A3-A[0]; dB3+=B3-A[1]; dC3+=C3-A[2];} // if last row: (dA3, dB3, dC3) := (dA3, dB3, dC3) + (A3, B3, C3) - (ax, ay, az)

      dA4 = (*cov)[i+3]+dB3*H1[2]-dC3*H1[1];    // dA4/dp }
      dB4 = (*cov)[i+4]+dC3*H1[0]-dA3*H1[2];    // dB4/dp  } = dA + (dA3, dB3, dC3) x H1
      dC4 = (*cov)[i+5]+dA3*H1[1]-dB3*H1[0];    // dC4/dp }

      if(i==42) {dA4+=A4-A[0]; dB4+=B4-A[1]; dC4+=C4-A[2];} // if last row: (dA4, dB4, dC4) := (dA4, dB4, dC4) + (A4, B4, C4) - (ax, ay, az)

      //last point
      dA5 = dA4+dA4-(*cov)[i+3];      // }
      dB5 = dB4+dB4-(*cov)[i+4];      //  } =  2*(dA4, dB4, dC4) - dA
      dC5 = dC4+dC4-(*cov)[i+5];      // }

      dA6 = dB5*H2[2]-dC5*H2[1];      // dA6/dp }
      dB6 = dC5*H2[0]-dA5*H2[2];      // dB6/dp  } = (dA5, dB5, dC5) x H2
      dC6 = dA5*H2[1]-dB5*H2[0];      // dC6/dp }

      if(i==42) {dA6+=A6; dB6+=B6; dC6+=C6;}     // if last row: (dA6, dB6, dC6) := (dA6, dB6, dC6) + (A6, B6, C6)

      if(i==42) {
        (*cov)[i]   += (dA2+dA3+dA4)*S3/state7[6];  (*cov)[i+3] = (dA0+dA3+dA3+dA5+dA6)*P3/state7[6]; // dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
        (*cov)[i+1] += (dB2+dB3+dB4)*S3/state7[6];  (*cov)[i+4] = (dB0+dB3+dB3+dB5+dB6)*P3/state7[6]; // dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
        (*cov)[i+2] += (dC2+dC3+dC4)*S3/state7[6];  (*cov)[i+5] = (dC0+dC3+dC3+dC5+dC6)*P3/state7[6];
      }
      else {
        (*cov)[i]   += (dA2+dA3+dA4)*S3;  (*cov)[i+3] = (dA0+dA3+dA3+dA5+dA6)*P3; // dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
        (*cov)[i+1] += (dB2+dB3+dB4)*S3;  (*cov)[i+4] = (dB0+dB3+dB3+dB5+dB6)*P3; // dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
        (*cov)[i+2] += (dC2+dC3+dC4)*S3;  (*cov)[i+5] = (dC0+dC3+dC3+dC5+dC6)*P3;
      }
    }
  }

  //
  // Track parameters in last point
  //
  R[0] += (A2+A3+A4)*S3;   A[0] += (SA[0]=(A0+A3+A3+A5+A6)*P3-A[0]);  // R  = R0 + S3*[(A2, B2, C2) +   (A3, B3, C3) + (A4, B4, C4)]
  R[1] += (B2+B3+B4)*S3;   A[1] += (SA[1]=(B0+B3+B3+B5+B6)*P3-A[1]);  // A  =     1/3*[(A0, B0, C0) + 2*(A3, B3, C3) + (A5, B5, C5) + (A6, B6, C6)]
  R[2] += (C2+C3+C4)*S3;   A[2] += (SA[2]=(C0+C3+C3+C5+C6)*P3-A[2]);  // SA = A_new - A_old

  // normalize A
  double CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]); // 1/|A|
  A[0] *= CBA; A[1] *= CBA; A[2] *= CBA;
}


double RKTrackRep::estimateStep(std::vector<GFPointPath>& points,
                                const TVector3& pos,
                                const TVector3& dir,
                                const M1x4& SU,
                                const GFDetPlane& plane,
                                const double& momentum,
                                double& relMomLoss,
                                double& deltaAngle,
                                bool& momLossExceeded,
                                bool& atPlane,
                                double maxStep) const {

  static const double Smax      = 10.;          // max. step allowed [cm]
  static const double dAngleMax = 0.05;         // max. deviation of angle between direction before and after the step [rad]
  double Step;
  bool improveEstimation (true);

  momLossExceeded = false;
  atPlane = false;

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
  TVector3 MagField(GFFieldManager::getFieldVal(pos));
  double Hmag(MagField.Mag()), SmaxAngle(maxStep), radius(0);
  if (Hmag > 1E-5){
    double cosAngle = (dir.Dot(MagField))/Hmag;
    radius = momentum/(0.299792458E-3*Hmag) *
             sqrt( pow(dir.X() - cosAngle/Hmag * MagField.X(), 2) +
                   pow(dir.Y() - cosAngle/Hmag * MagField.Y(), 2) +
                   pow(dir.Z() - cosAngle/Hmag * MagField.Z(), 2)); // [cm]
    double sinAngle = sqrt(1 - cosAngle*cosAngle);
    if (sinAngle > 1E-10) SmaxAngle = fabs(dAngleMax * radius / sinAngle); // [cm]
  }


  //
  // Select direction
  //
  // auto select
  if (fDirection == 0 || !plane.isFinite()){
    #ifdef DEBUG
      std::cout << "  auto select direction";
      if (!plane.isFinite()) std::cout << ", plane is not finite";
      std::cout << ".\n";
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
      improveEstimation = false;
      #ifdef DEBUG
        std::cout << "  we are near the plane, but not pointing to the active area. make a big step! \n";
      #endif
    }
  }
  // fDirection decides!
  else {
    if (Step * fDirection < 0){
      Step = fDirection*SmaxAngle;
      improveEstimation = false;
      #ifdef DEBUG
        std::cout << "  invert Step according to fDirection and set Step to fDirection*SmaxAngle. \n";
      #endif
    }
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
  // reduce maximum stepsize Step to Smax and maxStep
  if (fabs(Step) > Smax) {
    Step = StepSign*Smax;
    improveEstimation = false;
  }
  if (fabs(Step) > maxStep) {
    Step = StepSign*maxStep;
    improveEstimation = false;
  }

  // also limit stepsize according to the change of the momentum direction!
  if (fabs(Step) > SmaxAngle) {
    Step = StepSign*SmaxAngle;
    improveEstimation = false;
  }

  #ifdef DEBUG
    std::cout << "  limit from maxangle: " << SmaxAngle << ", radius: " << radius << "\n";
  #endif


  // call stepper and reduce stepsize if step not too small
  if (!fNoMaterial){

    if(fabs(Step) > MINSTEP){ // only call stepper if step estimation big enough
      double StepMat = GFMaterialEffects::getInstance()->stepper(fabs(Step),
                                                                 SmaxAngle,
                                                                 pos,
                                                                 StepSign*dir,
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
  

  if (!momLossExceeded && improveEstimation){
    atPlane = true;
    
    // improve step estimation to surface according to curvature
    if (Hmag > 1E-5 && fabs(Step) > 0.1*SmaxAngle){

      M1x7 state7;
      state7[0] = pos.X();  state7[1] = pos.Y(); state7[2] = pos.Z();
      state7[3] = dir.X();  state7[4] = dir.Y(); state7[5] = dir.Z();
      state7[6] = fCharge/momentum;
      M1x3 SA;

      RKPropagate(state7, NULL, SA, Step, false);

      // calculate distance to surface
      Dist = SU[3] - (state7[0] * SU[0] +
                      state7[1] * SU[1] +
                      state7[2] * SU[2]); // Distance between position and surface

      An = state7[3] * SU[0] +
           state7[4] * SU[1] +
           state7[5] * SU[2];    // An = dir * N;  component of dir normal to surface

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