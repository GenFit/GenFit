#include "PixHit.h"

#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"
#include "TMath.h"

#include "math.h"

ClassImp(PixHit)


PixHit::~PixHit()
{}

PixHit::PixHit()
  : GFAbsPlanarHit(NparHitRep)
{}

PixHit::PixHit(const TVector3& point, const TVector3& norm,
		               const TVector3& u, double res,
		               bool smear)
  : GFAbsPlanarHit(NparHitRep){

  fHitCov(0,0) = res*res;
  fHitCov(1,1) = res*res;

  assert(fabs(norm*u)<1.E-5);
  TVector3 v = u.Cross(norm);
  GFDetPlane plane(point,u,v);

  if(smear) {
    fHitCoord(0) = gRandom->Gaus(0,res);
    fHitCoord(1) = gRandom->Gaus(0,res);
  }
  else {
    fHitCoord(0) = 0;
    fHitCoord(1) = 0;
  }

  setDetPlane(plane);
}


GFAbsRecoHit* 
PixHit::clone(){
  return new PixHit(*this);
}


const TMatrixT<double>&
PixHit::getHMatrix(const GFAbsTrackRep* rep) {

  if ( (dynamic_cast<const RKTrackRep*>(rep) != NULL)){
    static const double HMatrixContent[10] = {0, 0, 0, 1, 0,
                                              0, 0, 0, 0, 1};
    static const TMatrixT<double> HMatrix(2,5, HMatrixContent);

    return HMatrix;
  }
  else {
    std::cerr << "PixHit can only handle state vectors of type RKtrackRep -> abort" << std::endl;
    throw;
  }
 
}


