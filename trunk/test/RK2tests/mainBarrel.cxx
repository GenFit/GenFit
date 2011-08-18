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
  TVector3 mom(0.,4.,0.8);
  mom.SetMag(.05);
  mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
  TVector3 posErr(1.,1.,1.);
  TVector3 momErr(1.,1.,1.);
  GFDetPlane plane(pos,mom);//TVector3(1.,0.,0.),TVector3(0.,1.,0.));

  GFFieldManager::getInstance()->init(new GFConstField(0.,0.,10.));

  gRandom->SetSeed(5);
  //12MeV
  const int nev = 1000;
  const double posSmear = 0.1;
  const double momSmear = 0.3;
  const int npoints = 1850;
  const double pointDist = .05;
  const double resolution = 0.03;
  const int numIT = 3;
  const bool eventDisplay = false;
  //  GFException::quiet();

  std::vector<double> rings;
  rings.push_back(2.);
  rings.push_back(5.);
  rings.push_back(11.);
  rings.push_back(18.);
  rings.push_back(25.);

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
    
    GFAbsTrackRep *rephits = new RKtrackRep(pos,//plane,
					     mom,
					     posErr,
					     momErr,
					     211);
    GFAbsTrackRep *rep = new RKtrackRep(posM,
					 momM,
					 posErr,
					 momErr,
					 211);
    

    //for saving of MC truth    
    stMCT->ResizeTo(rephits->getState());
    *stMCT = rephits->getState();
    covMCT->ResizeTo(rephits->getCov());
    *covMCT = rephits->getCov();
    
    
    std::vector<TVector3> points;
    std::vector<TVector3> pointsOther;
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
      for(unsigned int j=0;j<rings.size();++j){
	if(fabs(posR.Perp()-rings.at(j))<pointDist/2.){
	  points.push_back(posR);	  
	  std::cout << i << " ";
	  break;
	}	
	else{
	  pointsOther.push_back(posR);	  
	}
      }



    }//loop npoints
    std::cout << std::endl;
    delete rephits;

    std::cerr << "generated " << points.size() << " hits" << std::endl;

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
      TPolyMarker3D *mOther = new TPolyMarker3D(points.size(),20);
      for(unsigned int i=0;i<points.size();++i){
	m->SetPoint(i,points.at(i).X(),points.at(i).Y(),points.at(i).Z());
      }
      for(unsigned int i=0;i<pointsOther.size();++i){
	mOther->SetPoint(i,pointsOther.at(i).X(),pointsOther.at(i).Y(),pointsOther.at(i).Z());
      }
      m->SetMarkerSize(1.);
      m->SetMarkerColor(kRed);
      mOther->SetMarkerSize(0.2);
      mOther->SetMarkerColor(kBlue);
      new TCanvas;

      mOther->Draw();
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
    std::cout << "fit is done" << std::endl;
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

