#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TGeoManager.h"

#include "GFException.h"
#include "GFAbsTrackRep.h"
#include "GeaneTrackRep2.h"
#include "RKtrackRep.h"

#include "GFConstField.h"
#include "GFFieldManager.h"

#include "PixHit.h"
#include "GFTrack.h"
#include "GFKalman.h"
#include "GFDaf.h"

#include "stdlib.h"

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
  gROOT->Macro("../../config/Geane.C");


  TVector3 pos(0.01,0.01,0.01);
  TVector3 mom(0.,4.,1.);
  mom.SetMag(2.);
  mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
  TVector3 posErr(1.,1.,1.);
  TVector3 momErr(1.,1.,1.);
  GFDetPlane plane(pos,mom);//TVector3(1.,0.,0.),TVector3(0.,1.,0.));

  GFFieldManager::getInstance()->init(new GFConstField(0.,0.,15.));

  gRandom->SetSeed(5);
  //12MeV
  const int nev = 1000;
  const double posSmear = 0.03;
  const double momSmear = 0.03;
  const int npoints = 19;
  const double pointDist = 2.7;
  const double resolution = 0.03;
  const int numIT = 3;
  const bool eventDisplay = false;
  //  GFException::quiet();
  for(unsigned int iev = 0; iev<nev; ++iev){
    std::cerr << "@@@@@@@@@@@@@@@@ Doing event #" << iev << std::endl;

    TVector3 posM(pos);
    posM.SetX(gRandom->Gaus(posM.X(),posSmear*posM.X()));
    posM.SetY(gRandom->Gaus(posM.Y(),posSmear*posM.Y()));
    posM.SetZ(gRandom->Gaus(posM.Z(),posSmear*posM.Z()));

    TVector3 momM(mom);
    momM.SetX(gRandom->Gaus(momM.X(),momSmear*momM.X()));
    momM.SetY(gRandom->Gaus(momM.Y(),momSmear*momM.Y()));
    momM.SetZ(gRandom->Gaus(momM.Z(),momSmear*momM.Z()));


    GFDetPlane planeM(posM,TVector3(1.,0.,0.),TVector3(0.,1.,0.));
    GFDetPlane plane(pos,mom);//TVector3(1.,0.,0.),TVector3(0.,1.,0.));
    
    GFAbsTrackRep *rephits = new GeaneTrackRep2(plane,
					      mom,
					      posErr,
					      momErr,
					      -211);
    GFAbsTrackRep *rep = new RKtrackRep(posM,
					  momM,
					  posErr,
					  momErr,
					  -211);
    

    //for saving of MC truth    
    stMCT->ResizeTo(rephits->getState());
    *stMCT = rephits->getState();
    covMCT->ResizeTo(rephits->getCov());
    *covMCT = rephits->getCov();
    
    
    std::vector<TVector3> points;
    std::cout << "Creating point: ";
    for(unsigned int i=0;i<npoints;++i){
      TVector3 pos = rephits->getPos();
      TVector3 mom = rephits->getMom();
      mom.SetMag(pointDist);
      GFDetPlane d(pos+mom,mom);
      //      d.Print();
      TMatrixT<double> statePred(6,1);
      TMatrixT<double> covPred(6,6);
      try{
	rephits->extrapolate(d,statePred,covPred);
      }
      catch(GFException& e){
	e.what();
	std::cerr<<"Exceptions wont be further handled ->exit(1) "<<__LINE__<<std::endl;
	throw;
      }
      rephits->setState(statePred);
      rephits->setCov(covPred);
      rephits->setReferencePlane(d);
      TVector3 posR(rephits->getPos());
      points.push_back(posR);
      std::cout << i << " ";
    }//loop npoints
    std::cout << std::endl;
    delete rephits;

    GFTrack fitTrack(rep);//initialized with smeared rep
    int ihit = 0;
    for(unsigned int i=0;i<points.size();++i){
      GFAbsRecoHit *h;
      if(i==0){
	h = new PixHit(points.at(i),points.at(i),resolution);
      }
      else{
	h = new PixHit(points.at(i),points.at(i)-points.at(i-1),resolution);
      }
      //	h->Print();
      fitTrack.addHit(h,
		      4,//dummy detector id
		      ihit++);

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
    //GFDaf daf;
    k.setNumIterations(numIT);
    try{
      k.processTrack(&fitTrack);
    }
    catch(GFException& e){
      e.what();
      std::cerr<<"Exceptions wont be further handled ->exit(1) "<<__LINE__<<std::endl;
      continue;
    }
    if(rep->getStatusFlag()!=0) continue;//go to next event
    //extrapolate rep to start plane, so it at the same plane as the MC truth
    try{
      rep->extrapolate(plane);
    }
    catch(GFException& e){
      e.what();
      std::cerr<<"Exceptions wont be further handled ->exit(1) line "<<__LINE__<<std::endl;
      continue;
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

