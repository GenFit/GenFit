// This Class' Header ------------------
#include "PointHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"

// Class Member definitions -----------

ClassImp(PointHit)


PointHit::~PointHit()
{}

PointHit::PointHit()
  : SpacepointRecoHit(NparHitRep)
{}

PointHit::PointHit(TVector3 point,TVector3 res)
  : SpacepointRecoHit(NparHitRep){

  fHitCov[0][0] = res(0)*res(0);
  fHitCov[1][1] = res(1)*res(1);
  fHitCov[2][2] = res(2)*res(2);

  fHitCoord[0][0] = point(0);
  fHitCoord[1][0] = point(1);
  fHitCoord[2][0] = point(2);

}

GFAbsRecoHit* 
PointHit::clone(){
  return new PointHit(*this);
}

TMatrixT<double>
PointHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if ((dynamic_cast<const RKTrackRep*>(stateVector) != NULL)) {
    //I know, since this is the same everytime, it could be done in the
    //the constructor, but I do it here anyway, to make clear that in the
    //case of several track-reps per hit, it would have to be done here
   //    fHMatrix.ResizeTo(NparHitRep,5);
   TMatrixT<double> HMatrix(2,5);

    HMatrix[0][0] = 0.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 1.;
    HMatrix[0][4] = 0.;

    HMatrix[1][0] = 0.;
    HMatrix[1][1] = 0.;
    HMatrix[1][2] = 0.;
    HMatrix[1][3] = 0.;
    HMatrix[1][4] = 1.;

    return HMatrix;
  } else {
   std::cerr << "PointHit can only handle state"
	     << " vectors of type GeaneTrackRep2 -> abort" << std::endl;
   throw;
 }
 
}


