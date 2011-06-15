#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include"TFile.h"
#include"TTree.h"
#include"TROOT.h"
#include"TMath.h"
#include"TStopwatch.h"
#include"TRandom.h"
#include"TGeoManager.h"

#include"GFException.h"
#include"GFAbsTrackRep.h"
#include"GeaneTrackRep2.h"
#include"RKtrackRep.h"
#include"SlTrackRep.h"

#include"GFConstField.h"
#include"GFFieldManager.h"

#include"PointHit.h"
#include"StripHit.h"
#include"PixHit.h"
#include"GFTrack.h"
#include"GFKalman.h"
#include"GFDaf.h"

#include"stdlib.h"

#define PRINT 0

TMatrixT<double> *stMCT;
TMatrixT<double> *covMCT;
TMatrixT<double> *stREC;
TMatrixT<double> *covREC;
double chi2;
int nfail;
int ndf;

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
  tree->Branch("ndf",&ndf,"ndf/I");

  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
    gROOT->Macro("config/Geane.C");


  TVector3 pos(0.01,0.01,0.01);
  TVector3 mom(1.,1.,1.);
  mom.SetMag(1.);
  mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
  TVector3 posErr(1.,1.,1.);
  TVector3 momErr(1.,1.,1.);
  //GFDetPlane plane(pos,mom);

  GFFieldManager::getInstance()->init(new GFConstField(0.,0.,10.));
  if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
  gRandom->SetSeed(5);
  //12MeV
  const int nev = 10000;
  const double posSmear = 0.015;
  const double momSmear = 0.015;
  const int npoints = 30;
  const double pointDist = 1.;
  const double resolution = 0.015;
  const int numIT = 3;
  const bool eventDisplay = false;
  GFException::quiet();
  TStopwatch watch;
  watch.Start();
  for(unsigned int iev = 0; iev<nev; ++iev){
    if((iev%1) == 0) std::cerr << "@@@@@@@@@@@@@@@@ Doing event #" << iev << std::endl;

    TVector3 posM(pos);
    posM.SetX(gRandom->Gaus(posM.X(),posSmear*posM.X()));
    posM.SetY(gRandom->Gaus(posM.Y(),posSmear*posM.Y()));
    posM.SetZ(gRandom->Gaus(posM.Z(),posSmear*posM.Z()));

    TVector3 momM(mom);
    momM.SetX(gRandom->Gaus(momM.X(),momSmear*momM.X()));
    momM.SetY(gRandom->Gaus(momM.Y(),momSmear*momM.Y()));
    momM.SetZ(gRandom->Gaus(momM.Z(),momSmear*momM.Z()));

  if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
  GFDetPlane planeM(posM-0.5/momM.Mag()*momM,momM);
    GFDetPlane plane(pos,mom);
  GFAbsTrackRep *rephits = new RKtrackRep(pos,//plane,
					    mom,
					    posErr,
					    momErr,
					    211);
  if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
  GFAbsTrackRep *rep = new RKtrackRep(posM-.5/momM.Mag()*momM,//planeM,
					  momM,
					  posErr,
					  momErr,
					  211);
      if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
      //    GFDetPlane plane = rephits->getReferencePlane();
    //    plane.Print();
    //for saving of MC truth    
    stMCT->ResizeTo(rephits->getState());
    *stMCT = rephits->getState();
    covMCT->ResizeTo(rephits->getCov());
    *covMCT = rephits->getCov();
      if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
    
    std::vector<TVector3> points;
    //    std::cout << "Creating point: ";
    for(unsigned int i=0;i<npoints;++i){
      if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
      TVector3 pos = rephits->getPos();
      TVector3 mom = rephits->getMom();
      if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
      mom.SetMag(pointDist);
      if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
      GFDetPlane d(pos+mom,mom);
      if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
      if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
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
      if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
      rephits->setState(statePred);
      rephits->setCov(covPred);
      rephits->setReferencePlane(d);
      TVector3 posR(rephits->getPos());
      points.push_back(posR);
      if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
      //      std::cout << i << " ";
    }//loop npoints
    //    std::cout << std::endl;
    delete rephits;
  if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
    GFTrack fitTrack(rep);//initialized with smeared rep
    int ihit = 0;
    for(unsigned int i=0;i<points.size();++i){
      if(i==0){
	GFAbsRecoHit *h = new PixHit(plane,resolution);
	fitTrack.addHit(h,
			4,//dummy detector id
			ihit++);
      }
      /*
      else if(i<15){
	GFAbsRecoHit *h = new StripHit(points.at(i),resolution,resolution,i);
	fitTrack.addHit(h,
			2,//dummy detector id
			ihit++);
	
      }
      */
      else{
      	fitTrack.addHit(new PointHit(points.at(i),resolution),
       			3,//dummy detector id
       			ihit++);
      }
    }
  if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;

    if(eventDisplay){
      TPolyMarker3D *m = new TPolyMarker3D(points.size(),20);
      for(unsigned int i=0;i<points.size();++i){
	m->SetPoint(i,points.at(i).X(),points.at(i).Y(),points.at(i).Z());
      }
      new TCanvas;
      m->Draw();
      app.Run();
    }
  if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
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
  if(PRINT)std::cout << __FILE__ << " " << __LINE__ << std::endl;
    //plane.Print();
    //rep->getReferencePlane().Print();
  //    assert(plane==rep->getReferencePlane());
  if(plane!=rep->getReferencePlane()) continue;
    /*
    //extrapolate rep to start plane, so it at the same plane as the MC truth
    try{
      rep->extrapolate(plane);
    }
    catch(GFException& e){
      e.what();
      std::cerr<<"Exceptions wont be further handled ->exit(1) line "<<__LINE__<<std::endl;
      continue;
    }
    */
    stREC->ResizeTo(rep->getState());
    *stREC = rep->getState();
    covREC->ResizeTo(rep->getCov());
    *covREC = rep->getCov();
    chi2 = rep->getChiSqu();
    ndf = rep->getNDF();
    nfail = fitTrack.getFailedHits();
    tree->Fill();
  }
  watch.Stop();
  watch.Print();
  tree->Write();
  file->Close();
}

