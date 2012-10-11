// This Class' Header ------------------
#include "WireHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"

// Class Member definitions -----------

ClassImp(WireHit)


WireHit::~WireHit()
{}

WireHit::WireHit()
  : WireRecoHit(NparHitRep)
{}

WireHit::WireHit(const TVector3& wire1, const TVector3& wire2, double rdrift, double res, bool smear)
  : WireRecoHit(NparHitRep){

  fHitCov(6,6) = res*res;

  fHitCoord(0,0) = wire1.X();
  fHitCoord(1,0) = wire1.Y();
  fHitCoord(2,0) = wire1.Z();

  fHitCoord(3,0) = wire2.X();
  fHitCoord(4,0) = wire2.Y();
  fHitCoord(5,0) = wire2.Z();

  if (smear){
    fHitCoord(6,0) = gRandom->Gaus(rdrift,res);
  }
  else {
    fHitCoord(6,0) = rdrift;
  }
}


GFAbsRecoHit* 
WireHit::clone(){
  return new WireHit(*this);
}


TMatrixT<double>
WireHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if ((dynamic_cast<const RKTrackRep*>(stateVector) != NULL)) {
    TMatrixT<double> HMatrix(1,5);

    HMatrix(0,0) = 0.;
    HMatrix(0,1) = 0.;
    HMatrix(0,2) = 0.;
    HMatrix(0,3) = 1.;
    HMatrix(0,4) = 0.;

    return HMatrix;
  }
  else {
    std::cerr << "WireHit can only handle state vectors of type RKTrackRep -> abort" << std::endl;
    throw;
  }
 
}


