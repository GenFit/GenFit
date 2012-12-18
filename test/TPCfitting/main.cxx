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
#include "RKTrackRep.h"

#include "GFConstField.h"
#include "GFFieldManager.h"

#include "PointHit.h"
#include "StripHit.h"
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

int main(int argc, char** argv){
  assert(argc==3);
  int nev(-1);
  std::istringstream str(argv[1]);
  str>>nev;
  assert(nev>0);

  

  stMCT  = new TMatrixT<double>;
  covMCT = new TMatrixT<double>;
  stREC  = new TMatrixT<double>;
  covREC = new TMatrixT<double>;

  TApplication app("app",NULL,NULL);

  TFile *file = TFile::Open(argv[2],"RECREATE");
  TTree *tree = new TTree("t","example output");
  
  tree->Branch("stMCT","TMatrixT<double>",&stMCT);
  tree->Branch("covMCT","TMatrixT<double>",&covMCT);
  tree->Branch("stREC","TMatrixT<double>",&stREC);
  tree->Branch("covREC","TMatrixT<double>",&covREC);
  tree->Branch("chi2",&chi2,"chi2/D");
  tree->Branch("nfail",&nfail,"nfail/I");

  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  gROOT->Macro("config/Geane.C");


  TVector3 pos(0.01,0.01,0.01);
  TVector3 mom(1.,0.,1.);
  mom.SetMag(.5);
  mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
  TVector3 posErr(1.,1.,1.);
  TVector3 momErr(1.,1.,1.);
  GFDetPlane plane(pos,TVector3(1.,0.,0.),TVector3(0.,1.,0.));

  GFFieldManager::getInstance()->init(new GFConstField(0.,0.,10.));

  gRandom->SetSeed(4);

  //  const int nev = 500;
  const double posSmear = 0.05;
  const double momSmear = 0.5;
  const int npoints = 10;
  const double pointDist = 3;
  const double resolution = 0.01;
  const int numIT = 1;
  const bool eventDisplay = false;
  //GFException::quiet();
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

    //dreggn    

    GFAbsTrackRep *rephits = new GeaneTrackRep2(plane,
					      mom,
					      posErr,
					      momErr,
					      211);
    GFAbsTrackRep *rep = new RKTrackRep(posM,
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
    std::cout << "Creating point: ";
    for(unsigned int i=0;i<npoints;++i){
      TVector3 pos = rephits->getPos();
      TVector3 mom = rephits->getMom();
      mom.SetMag(pointDist);
      GFDetPlane d(pos+mom,mom);

      TMatrixT<double> statePred(6,1);
      TMatrixT<double> covPred(6,6);
      try{
	rephits->extrapolate(d,statePred,covPred);
      }
      catch(GFException& e){
	e.what();
	std::cerr<<"Exceptions wont be further handled ->exit(1)"<<std::endl;
	//exit(1);
	continue;
      }
      //rephits->setState(statePred);
      //rephits->setCov(covPred);
      //rephits->setReferencePlane(d);
      rephits->setData(statePred,d,&covPred);
      TVector3 posR(rephits->getPos());
      points.push_back(posR);
      std::cout.flush();
      std::cout << i << " ";
    }//loop npoints
    std::cout << std::endl;
    delete rephits;

    GFTrack fitTrack(rep);//initialized with smeared rep

    int countHits(0);    
    for(unsigned int i=0;i<points.size();++i){
      if(i<10){
	fitTrack.addHit(new StripHit(points.at(i),resolution/4.,resolution/4.,i,1),
			3,//dummy detector id
			countHits++,points.at(i).Mag(),i);
	fitTrack.addHit(new StripHit(points.at(i),5*resolution/4.,resolution/4.,i,1),
				3,//dummy detector id
			countHits++,points.at(i).Mag(),i);
      }
      else{
	/*if(i>3) if(gRandom->Uniform()<.3) {
	    TVector3 addV(1.,1,0);
	    addV.SetPhi(gRandom->Uniform(0.,2.*TMath::Pi()));
	    points.at(i) = points.at(i)+addV;
	    }*/
	fitTrack.addHit(new PointHit(points.at(i),resolution),
			3,//dummy detector id
			countHits++,i);
      }
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

    //GFKalman k;
    GFDaf k;
    //dreggn
    //k.setNumIterations(numIT);
    //k.setNumIterations(6);
    //    k.setInitialDirection(-1);
    try{
      k.processTrack(&fitTrack);
      }
    catch(GFException& e){
      std::cerr << e.what();
      
      std::cerr<<"Exceptions wont be further handled ->exit(1)  line " << __LINE__<<std::endl;
      //exit(1);
      continue;
      //throw e;
    }
    //    fitTrack.printBookkeeping();
    std::cout << momM.Mag() << std::endl;
    std::cout << fitTrack.getMom().Mag() << std::endl;
    std::string dummy;
    //    std::cin >> dummy;
    std::cout << "ndf " << rep->getNDF() << std::endl;
    //extrapolate rep to start plane, so it at the same plane as the MC truth
    if(rep->getStatusFlag()==0){
      try{
	rep->extrapolate(plane);
      }
      catch(GFException& e){
	e.what();
	std::cerr<<"Exceptions wont be further handled ->exit(1)  line " << __LINE__<<std::endl;
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
  }
  tree->Write();
  file->Close();
}

