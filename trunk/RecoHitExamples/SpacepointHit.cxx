// This Class' Header ------------------
#include "SpacepointHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"

// Class Member definitions -----------

ClassImp(SpacepointHit)


SpacepointHit::~SpacepointHit()
{}

SpacepointHit::SpacepointHit()
  : GFAbsSpacepointHit()
{}

SpacepointHit::SpacepointHit(const TVector3& point, const double& res, bool smear)
  : GFAbsSpacepointHit(){

  fHitCov(0,0) = res*res;
  fHitCov(1,1) = res*res;
  fHitCov(2,2) = res*res;

  if (smear){
    fHitCoord(0) = gRandom->Gaus(point.X(),res);
    fHitCoord(1) = gRandom->Gaus(point.Y(),res);
    fHitCoord(2) = gRandom->Gaus(point.Z(),res);
  }
  else {
    fHitCoord(0) = point.X();
    fHitCoord(1) = point.Y();
    fHitCoord(2) = point.Z();
  }
}

SpacepointHit::SpacepointHit(const TVector3& point, const TVector3& res, bool smear)
  : GFAbsSpacepointHit(){

  fHitCov(0,0) = res.X()*res.X();
  fHitCov(1,1) = res.Y()*res.Y();
  fHitCov(2,2) = res.Z()*res.Z();

  if (smear){
    fHitCoord(0) = gRandom->Gaus(point.X(),res.X());
    fHitCoord(1) = gRandom->Gaus(point.Y(),res.Y());
    fHitCoord(2) = gRandom->Gaus(point.Z(),res.Z());
  }
  else {
    fHitCoord(0) = point.X();
    fHitCoord(1) = point.Y();
    fHitCoord(2) = point.Z();
  }
}


GFAbsRecoHit* 
SpacepointHit::clone(){
  return new SpacepointHit(*this);
}


const TMatrixT<double>&
SpacepointHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if ((dynamic_cast<const RKTrackRep*>(stateVector) != NULL)) {
    static const double HMatrixContent[10] = {0, 0, 0, 1, 0,
                                              0, 0, 0, 0, 1};
    static const TMatrixT<double> HMatrix(2,5, HMatrixContent);

    return HMatrix;
  }
 else {
   std::cerr << "SpacepointHit can only handle state vectors of type RKTrackRep -> abort" << std::endl;
   throw;
 }
 
}


