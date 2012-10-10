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
  : WirePointRecoHit(NparHitRep)
{}

WirePointHit::WirePointHit(const TVector3& wire1, const TVector3& wire2, double rdrift, double z,
                           double resDrift, double resZ, bool smear)
  : WirePointRecoHit(NparHitRep){

  fHitCov(6,6) = resDrift*resDrift;
  fHitCov(7,7) = resZ*resZ;

  fHitCoord(0,0) = wire1.X();
  fHitCoord(1,0) = wire1.Y();
  fHitCoord(2,0) = wire1.Z();

  fHitCoord(3,0) = wire2.X();
  fHitCoord(4,0) = wire2.Y();
  fHitCoord(5,0) = wire2.Z();

  if (smear){
    fHitCoord(6,0) = gRandom->Gaus(rdrift,resDrift);
    fHitCoord(7,0) = gRandom->Gaus(z,resZ);
  }
  else {
    fHitCoord(6,0) = rdrift;
    fHitCoord(7,0) = z;
  }
}


GFAbsRecoHit* 
WirePointHit::clone(){
  return new WirePointHit(*this);
}


TMatrixT<double>
WirePointHit::getHMatrix(const GFAbsTrackRep* stateVector)
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
    std::cerr << "WirePointHit can only handle state vectors of type RKTrackRep -> abort" << std::endl;
    throw;
  }
 
}


