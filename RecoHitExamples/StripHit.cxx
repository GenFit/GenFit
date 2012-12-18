#include "StripHit.h"

#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"
#include "TMath.h"

#include "math.h"

ClassImp(StripHit)


StripHit::~StripHit()
{}

StripHit::StripHit()
  : GFAbsPlanarHit(NparHitRep)
{}

StripHit::StripHit(const TVector3& point, const TVector3& norm,
		               const TVector3& u, double res,
		               bool smear)
  : GFAbsPlanarHit(NparHitRep){

  fHitCov(0,0) = res*res;

  assert(fabs(norm*u)<1.E-5);
  TVector3 v = u.Cross(norm);
  GFDetPlane plane(point,u,v);

  if(smear) fHitCoord(0) = gRandom->Gaus(0,res);
  else fHitCoord(0) = 0;

  setDetPlane(plane);
}


GFAbsRecoHit* 
StripHit::clone(){
  return new StripHit(*this);
}


const TMatrixT<double>&
StripHit::getHMatrix(const GFAbsTrackRep* rep) {

  if ( (dynamic_cast<const RKTrackRep*>(rep) != NULL)){
    static const double HMatrixContent[5] = {0, 0, 0, 1, 0};
    static const TMatrixT<double> HMatrix(1,5, HMatrixContent);

    return HMatrix;
  }
  else {
    std::cerr << "StripHit can only handle state vectors of type RKtrackRep -> abort" << std::endl;
    throw;
  }
 
}


