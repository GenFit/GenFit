#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "TGeoManager.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "GFException.h"
#include "GFAbsTrackRep.h"
#include "RKTrackRep.h"
#include "GFConstField.h"
#include "GFFieldManager.h"

#include "PixHit.h"
#include "GFTrack.h"
#include "GFKalman.h"
#include "math.h"


#include "math.h"

TMatrixT<double> *stMCT;
TMatrixT<double> *covMCT;
TMatrixT<double> *stREC;
TMatrixT<double> *covREC;
double chi2;
int nfail;


int main(){


  stMCT  = new TMatrixT<double>;
  covMCT = new TMatrixT<double>;
  stREC  = new TMatrixT<double>;
  covREC = new TMatrixT<double>;

  TApplication app("app",NULL,NULL);

  TFile *file = TFile::Open("out.root","RECREATE");
  TTree *tree = new TTree("t","example output");
  
  tree->Branch("stMCT","TMatrixT<double>",&stMCT);
  tree->Branch("covMCT","TMatrixT<double>",&covMCT);
  tree->Branch("stREC","TMatrixT<double>",&stREC);
  tree->Branch("covREC","TMatrixT<double>",&covREC);
  tree->Branch("chi2",&chi2,"chi2/D");
  tree->Branch("nfail",&nfail,"nfail/I");

  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");

  TVector3 pos(0.01,0.01,0.01);
  TVector3 mom(1.,1.,1.);
  mom.SetMag(.05);
  mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
  TVector3 posErr(1.,1.,1.);
  TVector3 momErr(1.,1.,1.);


  GFFieldManager::getInstance()->init(new GFConstField(0.,0.,15));

  gRandom->SetSeed(5);

  const int nev = 1000;
  const double posSmear = 0.0001;
  const double momSmear = 0.0001;
  const int npoints = 800;
  const double pointDist = .1;
  const double resolution = 0.015;
  const bool eventDisplay = false;
  for(unsigned int iev = 0; iev<nev; ++iev){
    std::cerr << "@@@@@@@@@@@@@@@@ Doing event #" << iev << std::endl;

    TVector3 posM(pos);
    posM.SetX(gRandom->Gaus(posM.X(),posSmear*posM.X()));
    posM.SetY(gRandom->Gaus(posM.Y(),posSmear*posM.Y()));
    posM.SetZ(gRandom->Gaus(posM.Z(),posSmear*posM.Z()));

    TVector3 momM(mom);
    momM.SetX(gRandom->Gaus(momM.X(),momSmear*momM.X()));
    momM.SetY(gRandom->Gaus(momM.Y(),momSmear*momM.Y()));
    momM.SetZ(fabs(gRandom->Gaus(momM.Z(),momSmear*momM.Z())));



    GFAbsTrackRep *rephits = new RKTrackRep(pos,
					  mom,
					  posErr,
					  momErr,
					  211);
    GFAbsTrackRep *rep = new RKTrackRep(posM,
				      momM,
				      posErr,
				      momErr,
				      211);
    
    /*
    GFDetPlane d1(TVector3(1.,-1.,0.3232),
		TVector3(1,1.,0.),
		TVector3(-1.,1.,0.));

    TMatrixT<double> s(5,1);
    TMatrixT<double> c(5,5);
    s[0][0] = 2.;
    s[1][0] = -1.;
    s[2][0] = 1.;
    s[3][0] = 1;
    s[4][0] = .5;

    c[0][0] = 3.;
    //    c[0][1] = 1.;
    //    c[1][0] = 1.;
    c[1][1] = 2.;
    c[2][2] = 2.1;
    c[3][3] = 2.;
    c[4][4] = 1.;
    rep->setReferencePlane(d1);
    rep->setState(s);
    rep->setCov(c);

    TMatrixT<double> sP(5,1);
    TMatrixT<double> cP(5,5);
    GFDetPlane d2(TVector3(2.,2.,2.),
		TVector3(1.,1.,0.),
		TVector3(-1.,1.,0.));
    
    rep->extrapolate(d2,sP,cP);
    exit(0);
    */
    ///////////////////////////
    GFDetPlane plane = rephits->getReferencePlane();

    //for saving of MC truth    
    stMCT->ResizeTo(rephits->getState());
    *stMCT = rephits->getState();
    covMCT->ResizeTo(rephits->getCov());
    *covMCT = rephits->getCov();
    
    
    std::vector<TVector3> points_0;
    //    std::cout << "Creating point: ";
    
    for(unsigned int i=0;i<npoints;++i){
      TVector3 pos = rephits->getPos();
      //TAG
      //      GFDetPlane d(pos+TVector3(0.,0.,pointDist),TVector3(2,1.,0.),TVector3(-1,2,0.));
      GFDetPlane d(pos+TVector3(0.,0.,pointDist),TVector3(1.,0.,0.),TVector3(0,1,0.));

      TMatrixT<double> statePred(5,1);
      TMatrixT<double> covPred(5,5);
      try{
	rephits->extrapolate(d,statePred,covPred);
      }
      catch(GFException& e){
	e.what();
	std::cerr<<"Exceptions wont be further handled ->exit(1)"<<std::endl;
	exit(1);
      }
      rephits->setData(statePred,d,&covPred);
      TVector3 posR(rephits->getPos());
      points_0.push_back(posR);
      //      std::cout << i << " ";
    }//loop npoints
    //    std::cout << std::endl;
    delete rephits;
    std::vector<TVector3> points_1;
    std::vector<TVector3> points;
    int i=0;
    while(i<points_0.size()){
      TVector3 tmp(points_0.at(i));
      tmp.SetZ(0.);
      double r = tmp.Mag();
      if(r<4.||(r-((int)r) >0.12)) {
	++i;
	continue;
      }
      int intR = r;
      if((intR%5) == 0) {
	TVector3 T(points_0.at(i));
	T.SetZ(0.);
	std::cout << i << " r " << T.Mag() << std::endl;
	points.push_back(points_0.at(i));
	i+=10;
      }
      else{
	++i;
      }
    }
    GFTrack fitTrack(rep);//initialized with smeared rep
    for(unsigned int i=0;i<points.size();++i){
      fitTrack.addHit(new PixHit(points.at(i),resolution),
		      3,//dummy detector id
		      i);
    }

    if(eventDisplay){
      TPolyMarker3D *m = new TPolyMarker3D(points.size(),20);
      for(unsigned int i=0;i<points.size();++i){
	m->SetPoint(i,points.at(i).X(),points.at(i).Y(),points.at(i).Z());
      }
      new TCanvas;
      m->Draw();
      app.Run();
    }

    GFKalman k;
    k.setNumIterations(3);
    try{
      k.processTrack(&fitTrack);
    }
    catch(GFException& e){
      e.what();
      std::cerr<<"Exceptions wont be further handled ->exit(1)"<<std::endl;
      exit(1);
    }

    if(rep->getStatusFlag() != 0 ) continue;

    //extrapolate rep to start plane, so it at the same plane as the MC truth
    try{
      rep->extrapolate(plane);
    }
    catch(GFException& e){
      e.what();
      std::cerr<<"Exceptions wont be further handled ->exit(1)"<<std::endl;
      exit(1);
    }
    
    stREC->ResizeTo(rep->getState());
    *stREC = rep->getState();
    covREC->ResizeTo(rep->getCov());
    *covREC = rep->getCov();
    chi2 = rep->getRedChiSqu();
    nfail = fitTrack.getFailedHits();
    tree->Fill();
  }
  tree->Write();
  file->Close();

}
