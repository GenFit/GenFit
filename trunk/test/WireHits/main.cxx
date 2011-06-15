#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include"TFile.h"
#include"TTree.h"
#include"TROOT.h"
#include"TMath.h"
#include"TRandom.h"
#include"TGeoManager.h"

#include"GFException.h"
#include"GFAbsTrackRep.h"
#include"RKTrackRep.h"

#include"GFConstField.h"
#include"GFFieldManager.h"

#include"WireHit.h"
#include"PixHit.h"
#include"GFTrack.h"
#include"GFKalman.h"
#include"GFDaf.h"

#include"stdlib.h"

TMatrixT<double> *stMCT;
TMatrixT<double> *stREC;
TMatrixT<double> *covREC;
double chi2;
int ndf;
int nfail;

int main(int argc, char** argv){
  assert(argc==3);
  int nev(-1);
  std::istringstream str(argv[1]);
  str>>nev;
  assert(nev>0);

  

  stMCT  = new TMatrixT<double>;
  stREC  = new TMatrixT<double>;
  covREC = new TMatrixT<double>;

  TApplication app("app",NULL,NULL);

  TFile *file = TFile::Open(argv[2],"RECREATE");
  TTree *tree = new TTree("t","example output");
  
  tree->Branch("stMCT","TMatrixT<double>",&stMCT);
  tree->Branch("stREC","TMatrixT<double>",&stREC);
  tree->Branch("covREC","TMatrixT<double>",&covREC);
  tree->Branch("chi2",&chi2,"chi2/D");
  tree->Branch("ndf",&ndf,"ndf/I");
  tree->Branch("nfail",&nfail,"nfail/I");

  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");

  gRandom->SetSeed(4);

  TVector3 pos(0.01,0.01,0.0);
  pos.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
  TVector3 mom(1.,1.,0.1);
  mom.SetMag(1.0);
  mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));

  //this is the reference plane for comparison of parameters after fit
  GFDetPlane plane(pos,TVector3(pos.X(),pos.Y(),0.));

  GFFieldManager::getInstance()->init(new GFConstField(0.,0.,20.));



  //  const int nev = 500;
  const double posSmear = 0.1;
  const double momSmear = 0.1;
  const int numIT = 3;
  const bool eventDisplay = false;
  //GFException::quiet();
  for(unsigned int iev = 0; iev<nev; ++iev){
    std::cerr << "@@@@@@@@@@@@@@@@ Doing event #" << iev << std::endl;

    TVector3 posM(pos);
    posM.SetX(gRandom->Gaus(posM.X(),posSmear));
    posM.SetY(gRandom->Gaus(posM.Y(),posSmear));
    posM.SetZ(gRandom->Gaus(posM.Z(),posSmear));

    TVector3 momM(mom);
    momM.SetX(gRandom->Gaus(momM.X(),momSmear));
    momM.SetY(gRandom->Gaus(momM.Y(),momSmear));
    momM.SetZ(gRandom->Gaus(momM.Z(),momSmear));

    //used to make points along the trajectory
    GFAbsTrackRep *rephits = new RKTrackRep(plane,
					    mom,
					    211);
    

    GFAbsTrackRep *rep = new RKTrackRep(posM,
					momM,
					211);
    

    //for saving of MC truth    
    stMCT->ResizeTo(rephits->getState());
    *stMCT = rephits->getState();

    GFTrack fitTrack(rep);//initialized with smeared rep

    //make 5 points for planar hits
    std::vector<TVector3> planePoints;
    for(unsigned int i=0;i<5;++i){
      TVector3 pos = rephits->getPos();
      TVector3 mom = rephits->getMom();
      mom.SetMag(5.);
      GFDetPlane d(pos+mom,mom);

      TMatrixT<double> statePred(5,1);
      TMatrixT<double> covPred(5,5);
      try{
        rephits->extrapolate(d);
      }
      catch(GFException& e){
        e.what();
        std::cerr<<"Exceptions wont be further handled ->exit(1)"<<std::endl;
        //exit(1);
        continue;
      }

      TVector3 posR(rephits->getPos());
      planePoints.push_back(posR);
    } //done making 5 points for pixel hits

    //make the reco hits for pixel detectors
    int planeCountHits(0);
    for(unsigned int i=0;i<planePoints.size();++i){
      fitTrack.addHit(new PixHit(planePoints.at(i),0.005),
		      2,//dummy detector id
		      planeCountHits,planePoints.at(i).Mag(),planeCountHits);
      ++planeCountHits;
    }//done with planar hits


    //make 20 points for wire hits
    std::vector<TVector3> points;
    for(unsigned int i=0;i<20;++i){

      TVector3 pos = rephits->getPos();
      TVector3 dir = rephits->getPos();
      dir.SetZ(0.);
      dir.SetMag(2.);
      
      GFDetPlane d(pos+dir,dir);

      TMatrixT<double> statePred(5,1);
      TMatrixT<double> covPred(5,5);
      try{
	rephits->extrapolate(d,statePred,covPred);
      }
      catch(GFException& e){
	e.what();
	std::cerr<<"Exceptions wont be further handled ->exit(1)"<<std::endl;
	continue;
      }
      rephits->setData(statePred,d,&covPred);
      TVector3 posR(rephits->getPos());
      points.push_back(posR);
    }
    
    //now we have 20 points along the track in radial distaces of 2 cm

    //next make the wire hits 5 straight layers, 5 stereo one way,
    //5 stereo the other way, 5 straight
    int countHits(0);    

    std::vector<TVector3> wiresP1;
    std::vector<TVector3> wiresP2;
    for(unsigned int i=0;i<points.size();++i){
      TVector3 fromZAxisToPoint(points.at(i));
      fromZAxisToPoint.SetZ(0.);
      TVector3 U(0.,0.,1.);
      TVector3 V = fromZAxisToPoint.Cross(U);
      V.SetMag(1.);
      
      TVector3 P0 = points.at(i) + gRandom->Uniform(-0.5,0.5)*V;
      TVector3 wireDir;
      if(i>=5 && i<10){
	wireDir = U+0.5/5.71059*V;//make for 5 deg angle
      }
      else{
	if(i>=10 && i<15){
	  wireDir = U-0.5/5.71059*V;//make for 5 deg angle
	}
	else{
	  wireDir = U;
	}
      }
      
      wireDir.SetMag(1.);
      wiresP1.push_back( P0-20.*wireDir );
      wiresP2.push_back( P0+20.*wireDir );

      TVector3 poca,normVec,poca_onwire;
      try{
	rephits->extrapolateToLine(wiresP1.at(wiresP1.size()-1),
				    wiresP2.at(wiresP2.size()-1),
				    poca,normVec,poca_onwire);
      }
      catch(GFException& e){
	e.what();
	std::cerr<<"Exceptions wont be further handled ->exit(1)"<<std::endl;
	exit(1);
	//continue;
      }
      fitTrack.addHit(new WireHit(wiresP1.at(wiresP1.size()-1),wiresP2.at(wiresP2.size()-1),
				  (poca-poca_onwire).Mag(),0.01),
		      3,//dummy detector id
		      countHits,points.at(i).Mag(),countHits);
      ++countHits;
    }

    delete rephits;

    
    if(eventDisplay){
      TPolyMarker3D *mPlane = new TPolyMarker3D(planePoints.size(),21);
      for(unsigned int i=0;i<planePoints.size();++i){
	mPlane->SetPoint(i,planePoints.at(i).X(),planePoints.at(i).Y(),planePoints.at(i).Z());
      }

      TPolyMarker3D *m = new TPolyMarker3D(points.size(),20);
      TPolyLine3D *l[points.size()];
      for(unsigned int i=0;i<points.size();++i){
	m->SetPoint(i,points.at(i).X(),points.at(i).Y(),points.at(i).Z());
	l[i] = new TPolyLine3D(2);
	l[i]->SetPoint(0,wiresP1.at(i).X(),wiresP1.at(i).Y(),wiresP1.at(i).Z());
	l[i]->SetPoint(1,wiresP2.at(i).X(),wiresP2.at(i).Y(),wiresP2.at(i).Z());
      }
      new TCanvas;
      mPlane->SetMarkerColor(kBlue);
      mPlane->Draw();
      m->Draw();
      for(unsigned int i=0;i<points.size();++i) l[i]->Draw();
      gApplication->ReturnFromRun();
      gSystem->Run();//app.Run();
    }
    

    GFKalman k;
    //GFDaf k;
    k.setNumIterations(numIT);
    //k.setNumIterations(6);
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
    
    std::cout << "seed momentum: " << momM.Mag() << std::endl;
    std::cout << "recontructed momentum: " << fitTrack.getMom().Mag() << std::endl;
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
      chi2 = rep->getChiSqu();
      ndf = rep->getNDF();
      nfail = fitTrack.getFailedHits();
      tree->Fill();
    }
  }
  tree->Write();
  file->Close();
}

