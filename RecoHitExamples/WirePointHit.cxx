// This Class' Header ------------------
#include "WirePointHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"

// Class Member definitions -----------

ClassImp(WirePointHit)


WirePointHit::~WirePointHit()
{}

WirePointHit::WirePointHit()
  : GFAbsWirepointHit()
{}

WirePointHit::WirePointHit(const TVector3& wire1, const TVector3& wire2, double rdrift, double z,
                           double resDrift, double resZ, bool smear)
  : GFAbsWirepointHit(){

  fHitCov(6,6) = resDrift*resDrift;
  fHitCov(7,7) = resZ*resZ;

  fHitCoord(0) = wire1.X();
  fHitCoord(1) = wire1.Y();
  fHitCoord(2) = wire1.Z();

  fHitCoord(3) = wire2.X();
  fHitCoord(4) = wire2.Y();
  fHitCoord(5) = wire2.Z();

  if (smear){
    fHitCoord(6) = gRandom->Gaus(rdrift,resDrift);
    fHitCoord(7) = gRandom->Gaus(z,resZ);
  }
  else {
    fHitCoord(6) = rdrift;
    fHitCoord(7) = z;
  }
}


GFAbsRecoHit* 
WirePointHit::clone(){
  return new WirePointHit(*this);
}


const TMatrixT<double>&
WirePointHit::getHMatrix(const GFAbsTrackRep* rep)
{
  if ((dynamic_cast<const RKTrackRep*>(rep) != NULL)) {
    static const double HMatrixContent[10] = {0, 0, 0, 1, 0,
                                              0, 0, 0, 0, 1};
    static const TMatrixT<double> HMatrix(2,5, HMatrixContent);

    return HMatrix;
  }
  else {
    std::cerr << "WirePointHit can only handle state vectors of type RKTrackRep -> abort" << std::endl;
    throw;
  }
 
}


