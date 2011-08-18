// This Class' Header ------------------
#include "StripHit.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "GeaneTrackRep2.h"
#include "RKTrackRep.h"
#include "GFDetPlane.h"
#include "TRandom.h"
#include "TMath.h"

// Class Member definitions -----------

ClassImp(StripHit)


StripHit::~StripHit()
{}

StripHit::StripHit()
  : PlanarRecoHit(NparHitRep)
{}

StripHit::StripHit(TVector3 point,double res,double sigma,int proj,int smearFlag)
  : PlanarRecoHit(NparHitRep){
  assert(smearFlag==1||smearFlag==2);
  fHitCov[0][0] = sigma*sigma;
  //  fHitCov[1][1] = res*res;
  GFDetPlane d;

  
  //std::cout << "proj " << proj << std::endl;
  TVector3 W(0.,0.,1.);
  TVector3 U(1.,0.,0.);
  //U.Print();
  if((proj % 2) == 0) {
    //std::cout << "even" << std::endl;
    U.SetXYZ(0.,1.,0.);
  }
  //U.Print();
  TVector3 V = W.Cross(U);

//   TVector3 O(gRandom->Uniform(10.,20.),
// 	     gRandom->Uniform(-20.,20.),
// 	     point.Z());
  TVector3 O(1.,1.,point.Z());
  d.setO(O);
  d.setU(U);
  d.setV(V);
  //point.Print();
  point -= O;
  //point.Print();
  double u,v;
  u = point*U;
  v = point*V;
  if(smearFlag==1) fHitCoord[0][0] = gRandom->Gaus(u,res);
  else fHitCoord[0][0] = gRandom->Uniform(u-res,u+res);
  //  fHitCoord[1][0] = gRandom->Gaus(v,res);
  //  std::cout << smearFlag << " " << res << " " << u << " " << fHitCoord[0][0] << std::endl;

  //  d.Print();
  //fHitCoord.Print();
  fPolicy.setDetPlane(d);
  //  exit(0);
}
GFAbsRecoHit* 
StripHit::clone(){
  return new StripHit(*this);
}


TMatrixT<double>
StripHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if ( (dynamic_cast<const RKTrackRep*>(stateVector) != NULL)){
   TMatrixT<double> HMatrix(1,5);

    HMatrix[0][0] = 0.;
    HMatrix[0][1] = 0.;
    HMatrix[0][2] = 0.;
    HMatrix[0][3] = 1.;
    HMatrix[0][4] = 0.;


    return HMatrix;
  }
  else   if ((dynamic_cast<const GeaneTrackRep2*>(stateVector) != NULL)){
   TMatrixT<double> HMatrix(1,6);

   HMatrix[0][0] = 0.;
   HMatrix[0][1] = 0.;
   HMatrix[0][2] = 0.;
   HMatrix[0][3] = 1.;
   HMatrix[0][4] = 0.;
   HMatrix[0][5] = 0.;
   
   
   return HMatrix;
   
  }
  else {
    std::cerr << "StripHit can only handle state"
	      << " vectors of type GeaneTrackRep2 or RKtrackRep -> abort" 
	      << std::endl;
    throw;
  }
 
}


