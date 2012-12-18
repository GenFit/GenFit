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
  : GFAbsWireHit()
{}

WireHit::WireHit(const TVector3& wire1, const TVector3& wire2, double rdrift, double res, bool smear)
  : GFAbsWireHit(){

  fHitCov(6,6) = res*res;

  fHitCoord(0) = wire1.X();
  fHitCoord(1) = wire1.Y();
  fHitCoord(2) = wire1.Z();

  fHitCoord(3) = wire2.X();
  fHitCoord(4) = wire2.Y();
  fHitCoord(5) = wire2.Z();

  if (smear){
    fHitCoord(6) = gRandom->Gaus(rdrift,res);
  }
  else {
    fHitCoord(6) = rdrift;
  }
}


GFAbsRecoHit* 
WireHit::clone(){
  return new WireHit(*this);
}


const TMatrixT<double>&
WireHit::getHMatrix(const GFAbsTrackRep* rep)
{
  if ((dynamic_cast<const RKTrackRep*>(rep) != NULL)) {
    static const double HMatrixContent[5] = {0, 0, 0, 1, 0};
    static const TMatrixT<double> HMatrix(1,5, HMatrixContent);

    return HMatrix;
  }
  else {
    std::cerr << "WireHit can only handle state vectors of type RKTrackRep -> abort" << std::endl;
    throw;
  }
 
}


