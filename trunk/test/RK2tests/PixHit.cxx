// This Class' Header ------------------
#include "PixHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "GeaneTrackRep2.h"
#include "RKtrackRep.h"
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

PixHit::PixHit(TVector3 point,TVector3 dir,double res)
  : PlanarRecoHit(NparHitRep){

  fHitCov[0][0] = res*res;
  fHitCov[1][1] = res*res;
  GFDetPlane d;

  //  res.Print();
  //  TVector3 W(dir);
  TVector3 temp(point);
  temp.SetZ(0.);
  temp.SetMag(1.);
  TVector3 W(temp);

  d.setO(point);
  TVector3 U(0.,0.,1.);
  TVector3 V = W.Cross(U);
  d.setUV(U,V);

  fHitCoord[0][0] = gRandom->Gaus(0,res);
  fHitCoord[1][0] = gRandom->Gaus(0,res);

  GFAbsFinitePlane * fin = new GFRectFinitePlane(-20,20,-3,3);
  d.setFinitePlane(fin);

  fPolicy.setDetPlane(d);
}

GFAbsRecoHit* 
PixHit::clone(){
  return new PixHit(*this);
}


TMatrixT<double>
PixHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if (dynamic_cast<const GeaneTrackRep2*>(stateVector) != NULL 
      ||
      dynamic_cast<const RKtrackRep*>(stateVector) != NULL ) {
    TMatrixT<double> HMatrix(2,6);
    
    HMatrix[0][0] = 0.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 1.;
    HMatrix[0][4] = 0.;
    HMatrix[0][5] = 0.;
    
    HMatrix[1][0] = 0.;
    HMatrix[1][1] = 0.;
    HMatrix[1][2] = 0.;
    HMatrix[1][3] = 0.;
    HMatrix[1][4] = 1.;
    HMatrix[1][5] = 0.;
    return HMatrix;
  }
  else {
    std::cerr << "PixHit can only handle state"
	      << " vectors of type GeaneTrackRep2 or RKtrackRep -> abort" 
	      << std::endl;
    throw;
  }
}


