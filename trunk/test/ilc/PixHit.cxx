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

#include <cmath>

// Class Member definitions -----------

ClassImp(PixHit)


PixHit::~PixHit()
{}

PixHit::PixHit()
  : PlanarRecoHit(2)
{}

#define PLANAR true

PixHit::PixHit(const TVector3& point) : PlanarRecoHit(2){
  //determine to which pixel detector the hit belongs to:
  //0==innerSilicon (6 lay.), 1==outerSilicon (2 lay.), 2==ftd
  double R = sqrt(point.X()*point.X()+point.Y()*point.Y());
  double Z = point.Z();
  int iDet(-1);
#ifndef PLANAR
  if(fabs(Z)<6.3 && R<6.1){
    iDet=0;//inner silicon tracker (6 lay.)
  }
  else{
    if((fabs(Z)<37.131 && R>16.499 && R<16.53) ||
       (fabs(Z)<64.491 && R>30.899 && R<30.93) ){
      iDet=1;//outer silison tracker (2 lay.)
    }
    else{
      iDet=2;//can only be ftd now
    }
  }

  TVector3 fromZAxisToPoint = point;
  fromZAxisToPoint.SetZ(0.);
  fromZAxisToPoint.SetMag(1.);
  TVector3 U,V;
  double dU,dV;
  switch(iDet){
  case 0:
    U.SetXYZ(0.,0.,1.);
    V = fromZAxisToPoint.Cross(U);
    dU = 0.0028e-1;
    dV = 0.0028e-1;
    break;
  case 1:
    U.SetXYZ(0.,0.,1.);
    V = fromZAxisToPoint.Cross(U);
    dU = 0.05e-1;//z resolution
    dV = 0.007e-1;//rphi resolution
    break;
  case 2:
    U.SetXYZ(1.,0.,0.);
    V.SetXYZ(0.,1.,0.);
    dU = 0.007e-1;
    dV = 0.007e-1;
    break;
  default:
    throw;//wont happen
  }

  GFDetPlane d(point,U,V);

  //std::cout << "hit det " << iDet << std::endl;
  //point.Print();
  //d.getNormal().Print();
  //d.Print();

#else
  double X = point.X();
  if(fabs(Z)<6.3 && X<6.1){
    iDet=0;//inner silicon tracker (6 lay.)
  }
  else{
    if((fabs(Z)<37.131 && X>16.499 && X<16.53) ||
       (fabs(Z)<64.491 && X>30.899 && X<30.93) ){
      iDet=1;//outer silison tracker (2 lay.)
    }
    else{
      iDet=2;//can only be ftd now
    }
  }

  //  std::cout << "hit in det " << iDet << std::endl;

  TVector3 fromZAxisToPoint = point;
  fromZAxisToPoint.SetZ(0.);
  fromZAxisToPoint.SetMag(1.);
  TVector3 U,V;
  double dU,dV;
  switch(iDet){
  case 0:
    U.SetXYZ(0.,0.,1.);
    V.SetXYZ(0.,-1.,0.);
    dU = 0.0028e-1;
    dV = 0.0028e-1;
    break;
  case 1:
    U.SetXYZ(0.,0.,1.);
    V.SetXYZ(0.,-1.,0.);
    dU = 0.05e-1;//z resolution
    dV = 0.007e-1;//rphi resolution
    break;
  case 2:
    U.SetXYZ(1.,0.,0.);
    V.SetXYZ(0.,1.,0.);
    dU = 0.007e-1;
    dV = 0.007e-1;
    break;
  default:
    throw;//wont happen
  }

  GFDetPlane d(point,U,V);
  


#endif
  //std::cout << "@@@@@@@@@@@@@@@@@@@@" << std::endl;
  //point.Print();
  //std::cout << iDet << std::endl;
  //d.Print();
  //d.getNormal().Print();
  fPolicy.setDetPlane(d);

  fHitCoord[0][0] = gRandom->Gaus(0.,dU);
  fHitCoord[1][0] = gRandom->Gaus(0.,dV);

  fHitCov[0][0] = dU*dU;
  fHitCov[1][1] = dV*dV;
}


GFAbsRecoHit* 
PixHit::clone(){
  return new PixHit(*this);
}


TMatrixT<double>
PixHit::getHMatrix(const GFAbsTrackRep* stateVector)
{
  if (   dynamic_cast<const RKTrackRep*>(stateVector) != NULL ) {
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
  else   if (dynamic_cast<const GeaneTrackRep2*>(stateVector) != NULL  ) {
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
	      << " vectors of type GeaneTrackRep2 or RKTrackRep -> abort" 
	      << std::endl;
    throw;
  }
}


