// This Class' Header ------------------
#include "PseudoSpacePointWireHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"
#include "TMatrixD.h"

// Class Member definitions -----------

ClassImp(PseudoSpacePointWireHit)


PseudoSpacePointWireHit::~PseudoSpacePointWireHit()
{}

PseudoSpacePointWireHit::PseudoSpacePointWireHit()
  : PseudoSpacePointWireRecoHit(NparHitRep)
{}

PseudoSpacePointWireHit::PseudoSpacePointWireHit(const TVector3& pos, const TVector3& wireDir,
                                                 double resPerp, double resWire, bool smear)
  : PseudoSpacePointWireRecoHit(NparHitRep){

  fHitCoord(0,0) = pos.X();
  fHitCoord(1,0) = pos.Y();
  fHitCoord(2,0) = pos.Z();

  fHitCov(0,0) = resPerp*resPerp;
  fHitCov(1,1) = resPerp*resPerp;
  fHitCov(2,2) = resWire*resWire;

  // rotation matrix
  TVector3 wDir(wireDir);
  wDir.SetMag(1);

  TVector3 xp = wDir.Orthogonal();
  xp.SetMag(1);
  TVector3 yp = wDir.Cross(xp);
  yp.SetMag(1);

  TMatrixD rot(3,3);

  rot(0,0) = xp.X();  rot(0,1) = yp.X();  rot(0,2) = wDir.X();
  rot(1,0) = xp.Y();  rot(1,1) = yp.Y();  rot(1,2) = wDir.Y();
  rot(2,0) = xp.Z();  rot(2,1) = yp.Z();  rot(2,2) = wDir.Z();


  if (smear) {
    TMatrixD smearVec(NparHitRep,1);
    TMatrixD smearVecRot(NparHitRep,1);
    smearVec(0,0) = resPerp;
    smearVec(1,0) = resPerp;
    smearVec(2,0) = resWire;
    smearVecRot.Mult(rot,smearVec);
    fHitCoord(0,0) += gRandom->Gaus(0, smearVecRot(0,0));
    fHitCoord(1,0) += gRandom->Gaus(0, smearVecRot(1,0));
    fHitCoord(2,0) += gRandom->Gaus(0, smearVecRot(2,0));
  }

  // rotate cov
  TMatrixD hitCovTemp(NparHitRep,NparHitRep);
  hitCovTemp.Mult(rot,fHitCov);
  fHitCov.MultT(hitCovTemp,rot);


  this->setWireDirection(wDir);
}


GFAbsRecoHit* 
PseudoSpacePointWireHit::clone(){
  return new PseudoSpacePointWireHit(*this);
}


TMatrixT<double>
PseudoSpacePointWireHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if ((dynamic_cast<const RKTrackRep*>(stateVector) != NULL)) {
    TMatrixT<double> HMatrix(2,5);

    HMatrix(0,0) = 0.;
    HMatrix(0,1) = 0.;
    HMatrix(0,2) = 0.;
    HMatrix(0,3) = 1.;
    HMatrix(0,4) = 0.;

    HMatrix(1,0) = 0.;
    HMatrix(1,1) = 0.;
    HMatrix(1,2) = 0.;
    HMatrix(1,3) = 0.;
    HMatrix(1,4) = 1.;

    return HMatrix;
  }
  else {
    std::cerr << "PseudoSpacePointWireHit can only handle state vectors of type RKTrackRep -> abort" << std::endl;
    throw;
  }
 
}


