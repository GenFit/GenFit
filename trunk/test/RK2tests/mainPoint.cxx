#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include"TFile.h"
#include"TTree.h"
#include"TROOT.h"
#include"TMath.h"
#include"TRandom.h"
#include"TGeoManager.h"

#include"GFException.h"
#include"GFAbsTrackRep.h"
#include"GeaneTrackRep2.h"
#include"RKtrackRep.h"

#include"GFConstField.h"
#include"GFFieldManager.h"

#include"PixHit.h"
#include"GFTrack.h"
#include"GFKalman.h"
#include"GFDaf.h"

#include"stdlib.h"
#include"math.h"

int main(){

  gRandom->SetSeed(5);
  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  GFFieldManager::getInstance()->init(new GFConstField(0.,20.,0.));

  TVector3 pos(0.0,0.0,0.0);
  TVector3 mom(0.,0.,1.);
  mom.SetMag(.5);
  //mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
  TVector3 posErr(1.,1.,1.);
  TVector3 momErr(1.,1.,1.);
  GFDetPlane plane(pos,mom);//TVector3(1.,0.,0.),TVector3(0.,1.,0.));
    
  GFAbsTrackRep *rephits = new RKtrackRep(pos,//plane,
					   mom,
					   posErr,
					   momErr,
					   211);

  TVector3 poca,dirInPoca,poca_onwire;
  TVector3 point(0.,-10.,12.);
  TVector3 point2(0.,10.,12.);

  GFDetPlane d(TVector3(0.,0.,1.03),TVector3(0.,0.0,1.));



  rephits->extrapolate(d);

  //  rephits->extrapolateToLine(point,point2,poca,dirInPoca,poca_onwire);  
  //poca.Print();
  //dirInPoca.Print();
  //poca_onwire.Print();

}

