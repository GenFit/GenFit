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

PointHit::PointHit(const TVector3& point,const double& res)
  : SpacepointRecoHit(NparHitRep){

  fHitCov[0][0] = res*res;
  fHitCov[1][1] = res*res;
  fHitCov[2][2] = res*res;

  fHitCoord[0][0] = gRandom->Gaus(point.X(),res);
  fHitCoord[1][0] = gRandom->Gaus(point.Y(),res);
  fHitCoord[2][0] = gRandom->Gaus(point.Z(),res);
}

PointHit::PointHit(const TVector3& point,const TVector3& res)
  : SpacepointRecoHit(NparHitRep){

  fHitCov[0][0] = res.X()*res.X();
  fHitCov[1][1] = res.Y()*res.Y();
  fHitCov[2][2] = res.Z()*res.Z();

  fHitCoord[0][0] = gRandom->Gaus(point.X(),res.X());
  fHitCoord[1][0] = gRandom->Gaus(point.Y(),res.Y());
  fHitCoord[2][0] = gRandom->Gaus(point.Z(),res.Z());
}


GFAbsRecoHit* 
PointHit::clone(){
  return new PointHit(*this);
}


TMatrixT<double>
PointHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if ((dynamic_cast<const RKTrackRep*>(stateVector) != NULL)) {
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
  }
 else {
   std::cerr << "PointHit can only handle state"
	     << " vectors of type RKTrackRep -> abort" << std::endl;
   throw;
 }
 
}


