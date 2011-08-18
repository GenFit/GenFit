// This Class' Header ------------------
#include "PixHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "GeaneTrackRep2.h"
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

PixHit::PixHit(TVector3 point,double res)
  : PlanarRecoHit(NparHitRep){

  fHitCov[0][0] = res*res;
  fHitCov[1][1] = res*res;
  GFDetPlane d;

  
  TVector3 N(point);
  N.SetZ(0.);
  //  N.SetXYZ(0.,0.,1.);
  d.setON(point,N);
  //  TVector3 O(point-0.1*d.getU()-0.1*d.getV());
  //  d.setO(point);
  //point.Print();
  //  point -= O;
  //point.Print();
  double u,v;
  //u = point*d.getU();
  //v = point*d.getV();
  //fHitCoord[0][0] = gRandom->Gaus(u,res);
  //fHitCoord[1][0] = gRandom->Gaus(v,res);
  fHitCoord[0][0] = gRandom->Gaus(0,res);
  fHitCoord[1][0] = gRandom->Gaus(0,res);

  GFAbsFinitePlane * fin = new GFRectFinitePlane(-6,6,-6,6);
  d.setFinitePlane(fin);

  //d.Print();
  //fHitCoord.Print();
  fPolicy.setDetPlane(d);
  //  exit(0);
}
GFAbsRecoHit* 
PixHit::clone(){
  return new PixHit(*this);
}


TMatrixT<double>
PixHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
 if (dynamic_cast<const GeaneTrackRep2*>(stateVector) != NULL) {
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
 else if(dynamic_cast<const RKTrackRep*>(stateVector) != NULL) {
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
	     << " vectors of type GeaneTrackRep2 or RKtrackRepXY -> abort" 
	     << std::endl;
   throw;
 }
 
}


