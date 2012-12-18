// This Class' Header ------------------
#include "ProlateSpacepointHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"
#include "TMatrixD.h"

// Class Member definitions -----------

ClassImp(ProlateSpacepointHit)


ProlateSpacepointHit::~ProlateSpacepointHit()
{}

ProlateSpacepointHit::ProlateSpacepointHit()
  : GFAbsProlateSpacepointHit()
{}

ProlateSpacepointHit::ProlateSpacepointHit(const TVector3& pos, const TVector3& wireDir,
                                                 double resPerp, double resWire, bool smear)
  : GFAbsProlateSpacepointHit(){

  fHitCoord(0) = pos.X();
  fHitCoord(1) = pos.Y();
  fHitCoord(2) = pos.Z();

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
    TVectorD smearVec(getNparHit());
    smearVec(0) = resPerp;
    smearVec(1) = resPerp;
    smearVec(2) = resWire;
    smearVec *= rot;
    fHitCoord(0) += gRandom->Gaus(0, smearVec(0));
    fHitCoord(1) += gRandom->Gaus(0, smearVec(1));
    fHitCoord(2) += gRandom->Gaus(0, smearVec(2));
  }

  // rotate cov
  fHitCov.Similarity(rot);

  setLargestErrorDirection(wDir);
}


GFAbsRecoHit* 
ProlateSpacepointHit::clone(){
  return new ProlateSpacepointHit(*this);
}


const TMatrixT<double>&
ProlateSpacepointHit::getHMatrix(const GFAbsTrackRep* rep)
{
  if ((dynamic_cast<const RKTrackRep*>(rep) != NULL)) {
    static const double HMatrixContent[10] = {0, 0, 0, 1, 0,
                                              0, 0, 0, 0, 1};
    static const TMatrixT<double> HMatrix(2,5, HMatrixContent);

    return HMatrix;
  }
  else {
    std::cerr << "ProlateSpacepointHit can only handle state vectors of type RKTrackRep -> abort" << std::endl;
    throw;
  }
 
}


