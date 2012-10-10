#include "PixHit.h"

#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"
#include "TMath.h"

#include"math.h"

ClassImp(PixHit)


PixHit::~PixHit()
{}

PixHit::PixHit()
  : PlanarRecoHit(NparHitRep)
{}

PixHit::PixHit(const TVector3& point, const TVector3& norm,
		               const TVector3& u, double res,
		               bool smear)
  : PlanarRecoHit(NparHitRep){

  fHitCov(0,0) = res*res;
  fHitCov(1,1) = res*res;

  assert(fabs(norm*u)<1.E-5);
  TVector3 v = u.Cross(norm);
  GFDetPlane d(point,u,v);

  if(smear) {
    fHitCoord(0,0) = gRandom->Gaus(0,res);
    fHitCoord(1,0) = gRandom->Gaus(0,res);
  }
  else {
    fHitCoord(0,0) = 0;
    fHitCoord(1,0) = 0;
  }

  fPolicy.setDetPlane(d);
}


GFAbsRecoHit* 
PixHit::clone(){
  return new PixHit(*this);
}


TMatrixT<double>
PixHit::getHMatrix(const GFAbsTrackRep* stateVector) {
  if ( (dynamic_cast<const RKTrackRep*>(stateVector) != NULL)){
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
    std::cerr << "PixHit can only handle state vectors of type RKtrackRep -> abort" << std::endl;
    throw;
  }
 
}


