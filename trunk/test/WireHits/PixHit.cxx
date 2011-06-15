// This Class' Header ------------------
#include "PixHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "GFRectFinitePlane.h"
#include "TRandom.h"
#include "TMath.h"

// Class Member definitions -----------

ClassImp(PixHit)

PixHit::~PixHit()
{}

PixHit::PixHit()
  : PlanarRecoHit(NparHitRep)
{}

PixHit::PixHit(const GFDetPlane& pl,double res) : PlanarRecoHit(NparHitRep){

  fPolicy.setDetPlane(pl);
  fHitCoord[0][0] = gRandom->Gaus(0.,res);
  fHitCoord[1][0] = gRandom->Gaus(0.,res);
  fHitCov[0][0] = res*res;
  fHitCov[1][1] = res*res;

}

PixHit::PixHit(TVector3 point,double res)
  : PlanarRecoHit(NparHitRep){

  fHitCov[0][0] = res*res;
  fHitCov[1][1] = res*res;
  GFDetPlane d;

  d.setO(point);
  d.setNormal(TVector3(point.X(),point.Y(),0.));

  
  fHitCoord[0][0] = gRandom->Gaus(0,res);
  fHitCoord[1][0] = gRandom->Gaus(0,res);

  //  GFAbsFinitePlane * fin = new GFRectFinitePlane(-20,20,-3,3);
  //  d.setFinitePlane(fin);

  fPolicy.setDetPlane(d);
}

GFAbsRecoHit* 
PixHit::clone(){
  return new PixHit(*this);
}


TMatrixT<double>
PixHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if ( dynamic_cast<const RKTrackRep*>(stateVector) != NULL ) {
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
    std::cerr << "PixHit can only handle state"
	      << " vectors of type RKTrackRep -> abort" 
	      << std::endl;
    throw;
  }
}


