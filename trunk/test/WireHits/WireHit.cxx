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

WireHit::WireHit(const TVector3& p1,const TVector3& p2,double driftLength,double res)
  : WireRecoHit(NparHitRep){

  fHitCoord[0][0] = p1.X();
  fHitCoord[1][0] = p1.Y();
  fHitCoord[2][0] = p1.Z();
  fHitCoord[3][0] = p2.X();
  fHitCoord[4][0] = p2.Y();
  fHitCoord[5][0] = p2.Z();
  fHitCoord[6][0] = gRandom->Gaus(driftLength,res);

  fHitCov[6][6] = res*res;

  fPolicy.setMaxDistance(2.);
}
GFAbsRecoHit* 
WireHit::clone(){
  return new WireHit(*this);
}


TMatrixT<double>
WireHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if ((dynamic_cast<const RKTrackRep*>(stateVector) != NULL)) {
    //I know, since this is the same everytime, it could be done in the
    //the constructor, but I do it here anyway, to make clear that in the
    //case of several track-reps per hit, it would have to be done here
   //    fHMatrix.ResizeTo(NparHitRep,5);
   TMatrixT<double> HMatrix(1,5);

    HMatrix[0][0] = 0.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 1.;
    HMatrix[0][4] = 0.;

    return HMatrix;
  }
  else {
    std::cerr << "WireHit can only handle state"
	      << " vectors of type RKTrackRep -> abort" << std::endl;
    throw;
  }
 
}


