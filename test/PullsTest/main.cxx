
#include <iostream>
#include <GenfitDisplay.h>
#include <GFConstField.h>
#include <GFException.h>
#include <GFFieldManager.h>
#include <GFKalman.h>
#include <GFTools.h>
#include <GFTrack.h>
#include <GFMaterialEffects.h>
#include <RKTrackRep.h>
#include <GFTGeoMaterialInterface.h>
#include <GeaneTrackRep2.h>
#include <GFRectFinitePlane.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TVector3.h>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include "TDatabasePDG.h"
#include <TMath.h>
#include <TString.h>

#include <TVirtualMC.h>
#include "TGeant3.h"

#include "PixHit.h"
#include "SpacepointHit.h"

//#define VALGRIND

int main() {
  std::cerr<<"main"<<std::endl;

  const unsigned int nEvents = 1000;
  const double BField = 15.;       // kGauss
  const double momentum = 0.1;     // GeV
  const double theta = 150;         // degree
  const double thetaDetPlane = 90;
  const unsigned int npoints = 7; // number of hits generated
  const double pointDist = 1;      // cm; approx. distance between hits generated
  const double resolution = 15.E-4;   // cm; resolution of generated hits
  const int pdg = 13;               // particle pdg code

  const bool smearPosMom = false;     // init the Reps with smeared pos and mom
  const double posSmear = 2*resolution;     // cm
  const double momSmear = 0.05*momentum;     // GeV
  
  const bool planarHits = true;    // true: create planar hits; false: create space point hits
  const bool GEANEhits = false;        // for creating the hits. true: use GeaneTrackRep2; false: use RKTrackRep
  const bool HelixTest = true;      // use helix for creating hits

  const int testRep = 0;          // 0=RKTrackRep, 1=RK7TrackRep, 2=GeaneTrackRep2, 3=RK+Geane
  const bool matFX = false;         // include material effects; can only be disabled for RKTrackRep!
  const bool smoothing = true;

  const bool debug = false;

  // init fitter
  GFKalman kalman;
  kalman.setNumIterations(3);

  // init mersenne twister
	gRandom->SetSeed(10);

#ifndef VALGRIND
  // init event display
	GenfitDisplay* display = GenfitDisplay::getInstance();
	display->reset();
#endif

  // init geometry and mag. field
  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  GFFieldManager::getInstance()->init(new GFConstField(0.,0.,BField));
  GFMaterialEffects::getInstance()->init(new GFTGeoMaterialInterface());

  // init Geane
  if (GEANEhits || testRep==2 || testRep==3) {
    if (debug) std::cerr<<"../../config/Geane.C"<<std::endl;
    gROOT->Macro("../../config/Geane.C");
  }
  
  const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);

  // init rootapp (for drawing histograms)
  TApplication* rootapp = new TApplication("rootapp", 0, 0);

  // prepare output tree for GFTracks
  GFTrack* trueTrack = new GFTrack();
  GFTrack* fitTrack = new GFTrack();
  TString outname = "out_Rep";
  outname += testRep;
  outname += "_degPlane";
  outname += thetaDetPlane;
  outname += ".root";
  TFile *file = TFile::Open(outname,"RECREATE");
  TTree *tree = new TTree("t","GFTracks");
  tree->Branch("trueTracks","GFTrack",&trueTrack);
  tree->Branch("fitTracks","GFTrack",&fitTrack);

  
  // create histograms
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptFit(1111);
  
  TH1D *hmomRes = new TH1D("hmomRes","mom res",2000,-0.02,0.02);
  TH1D *hupRes = new TH1D("hupRes","u' res",2000,-0.05,0.05);
  TH1D *hvpRes = new TH1D("hvpRes","v' res",2000,-0.05,0.05);
  TH1D *huRes = new TH1D("huRes","u res",2000,-0.05,0.05);
  TH1D *hvRes = new TH1D("hvRes","v res",2000,-0.05,0.05);

  TH1D *hqopPu = new TH1D("hqopPu","q/p pull",200,-6.,6.);
  TH1D *pVal = new TH1D("pVal","p-value",100,0.,1.);
  TH1D *hupPu = new TH1D("hupPu","u' pull",200,-6.,6.);
  TH1D *hvpPu = new TH1D("hvpPu","v' pull",200,-6.,6.);
  TH1D *huPu = new TH1D("huPu","u pull",200,-6.,6.);
  TH1D *hvPu = new TH1D("hvPu","v pull",200,-6.,6.);


  // main loop
  for (unsigned int iEvent=0; iEvent<nEvents; ++iEvent){
      
      if (debug || (iEvent+1)%10==0) std::cout << iEvent+1 << std::endl;

      // true start values
      TVector3 pos(0, 0, 0);
      TVector3 mom(1.,0,0);
      mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
      //mom.SetTheta(gRandom->Uniform(0.4*TMath::Pi(),0.6*TMath::Pi()));
      mom.SetTheta(theta*TMath::Pi()/180);
      mom.SetMag(momentum);

      // calc helix parameters
      TVector3 dir2D(mom);
      dir2D.SetZ(0);
      dir2D.SetMag(1.);
      double R = 100.*mom.Perp()/(0.03*BField);
      double sgn = 1;
      if (charge<0) sgn=-1.;
      TVector3 center = pos + sgn * R * dir2D.Orthogonal();
      double alpha0 = (pos-center).Phi();



      TVector3 posErr(1.,1.,1.);
      posErr *= posSmear;
      TVector3 momErr(1.,1.,1.);
      momErr *= momSmear;
      
      // smeared start values
      TVector3 posM(pos);
      TVector3 momM(mom);
      if (smearPosMom) {
        posM.SetX(gRandom->Gaus(posM.X(),posSmear));
        posM.SetY(gRandom->Gaus(posM.Y(),posSmear));
        posM.SetZ(gRandom->Gaus(posM.Z(),posSmear));

        momM.SetX(gRandom->Gaus(momM.X(),momSmear));
        momM.SetY(gRandom->Gaus(momM.Y(),momSmear));
        momM.SetZ(gRandom->Gaus(momM.Z(),momSmear));
      }

        
      // trackrep for creating hits
      GFAbsTrackRep* rephits;
      if (!GEANEhits) {
        rephits = new RKTrackRep(pos, mom, posErr, momErr, pdg);
        ((RKTrackRep*)rephits)->setPropDir(1);
      }
      else rephits = new GeaneTrackRep2(GFDetPlane(pos, mom), mom, posErr, momErr, pdg);

      if (!matFX) GFMaterialEffects::getInstance()->setNoEffects();
      					      
      // remember original initial plane and state					  
      GFDetPlane referencePlane;
      TVectorT<double> referenceState(rephits->getState());
      
      // create smeared hits
      std::vector<GFAbsRecoHit*> hits;
        
      TVector3 point, dir;
      if (debug) std::cerr << "Start creating hits ... \n";
      try{
        for (unsigned int i=0; i<npoints; ++i){
          // get current position and momentum
          if (!HelixTest) rephits->getPosMom(rephits->getReferencePlane(), point, dir);
          else{
            double angle = alpha0 - sgn * 0.1*i;
            TVector3 radius(R,0,0);
            radius.SetPhi(angle);
            point = center + radius;
            point.SetZ(pos.Z() + ((alpha0-angle)*R * TMath::Tan(mom.Theta()-TMath::Pi()*0.5) ));
            if (debug) point.Print();
            
            dir = mom;
            dir.SetPhi(mom.Phi()+(angle-alpha0));
            dir.SetMag(1);
          }
          // create hit
          if (i==0){ // get reference state
            TVector3 planeNorm(dir);
            planeNorm.SetTheta(thetaDetPlane*TMath::Pi()/180);
            TVector3 z(0,0,1);
            //z.SetTheta(thetaDetPlane*TMath::Pi()/180-TMath::PiOver2());
            referencePlane = GFDetPlane(point, planeNorm.Cross(z), (planeNorm.Cross(z)).Cross(planeNorm));
            if (!HelixTest) rephits->extrapolate(referencePlane, referenceState);
            else{
              double AtW = dir*referencePlane.getNormal();
              if (AtW<0) AtW *= -1.;
              referenceState[0] = charge/momentum;
              referenceState[1] = -1.*dir*referencePlane.getU()/AtW;
              referenceState[2] = -1.*dir*referencePlane.getV()/AtW;
              referenceState[3] = (point-referencePlane.getO())*referencePlane.getU();
              referenceState[4] = (point-referencePlane.getO())*referencePlane.getV();
              //std::cout<<"referenceState[2][0] "<<referenceState[2][0]<<"\n";
            }
          }
          //if (i>0){
            if (planarHits) {
              TVector3 planeNorm(dir);
              planeNorm.SetTheta(thetaDetPlane*TMath::Pi()/180);
              TVector3 z(0,0,1);
              //z.SetTheta(thetaDetPlane*TMath::Pi()/180+TMath::PiOver2());
              PixHit* hit = new PixHit(point,
                                       planeNorm,
                                       planeNorm.Cross(z),
                                       resolution,
                                       true);
              hits.push_back(hit);
            }
            else {
              SpacepointHit* hit = new SpacepointHit(point, resolution, true);
              hits.push_back(hit);
            }
          //}
          // stepalong (approximately)
          if (!HelixTest) {
            dir.SetMag(pointDist);
            GFDetPlane pl(point+dir, dir);
            rephits->extrapolate(pl);
          }
        }
      }
      catch(GFException& e){
	      std::cerr<<"Exception, next track"<<std::endl;
	      e.what();
	      continue; // here is a memleak!
      } 
      
      if (debug) std::cerr << "... done creating hits \n";
      
      
      
      // trackrep to be fitted and tested
      GFAbsTrackRep* rep;
      switch(testRep){
        case 0:
        case 3:
          rep = new RKTrackRep(posM, momM, posErr, momErr, pdg);
          //((RKTrackRep*)rep)->setPropDir(1);
          if (debug) {std::cout<<"init cov "; rep->getCov().Print();}
          break;
        case 1:
          std::cerr<<"Fieser Fehler!"<<std::endl;
          throw;
          //rep = new RK7TrackRep(posM, momM, /*posErr, momErr,*/ pdg);
          //((RKTrackRep*)rep)->setPropDir(1);
          break;
        case 2:
          rep = new GeaneTrackRep2(GFDetPlane(posM, momM), momM, posErr, momErr, pdg);
          break;
        default:
          std::cerr<<"Fieser Fehler!"<<std::endl;
          throw;
      }
      				          
		  // create track
      if (fitTrack != NULL) delete fitTrack;
		  fitTrack = new GFTrack(rep, smoothing); //initialized with smeared rep
		  
		  // optionally add GeaneTrackRep
		  if (testRep==3){
		    GFAbsTrackRep* rep2 = new GeaneTrackRep2(GFDetPlane(posM, momM), momM, posErr, momErr, pdg);
		    fitTrack->addTrackRep(rep2);
		  }
		  
		  // add hits
      for(unsigned int i=0; i<hits.size(); ++i){
        fitTrack->addHit(hits[i],
		                    3,//dummy detector id
		                    i);
      }
		  
      //test
      /*std::cout<<"pos and mom at last hit without fit\n";
      rep->getPos(hits.back()->getDetPlane(rep)).Print();
      rep->getMom(hits.back()->getDetPlane(rep)).Print();
      TMatrixD statePred(5,1);
      rep->extrapolate(hits.back()->getDetPlane(rep), statePred);
      statePred.Print();*/

		  
		  // do the fit
      try{
        if (debug) std::cerr<<"Starting the fitter"<<std::endl;
        kalman.processTrack(fitTrack);
        if (debug) std::cerr<<"fitter is finished!"<<std::endl;
      }
      catch(GFException& e){
        e.what();
	      std::cerr<<"Exception, next track"<<std::endl;
	      continue;
      }


      // check if fit was successfull
      if(rep->getStatusFlag() != 0 ) {
        continue;
      }

#ifndef VALGRIND
      // add track to event display
      std::vector<GFTrack*> event;
			event.push_back(fitTrack);
			display->addEvent(event);
#endif

			if (debug) {
			  std::cerr << "cov before extrapolating back to reference plane \n";
			  rep->getCov().Print();
			}
				          
		  // extrapolate back to reference plane. Not needed for planar hits!
		  if (!planarHits) rep->extrapolate(referencePlane);
		  

		  // check getPos etc methods
		  const double epsilon = 1E-3;
		  TVector3 getPos, getMom;

		  getPos = rep->getPos(referencePlane);
		  getMom = rep->getMom(referencePlane);
		  if ((pos-getPos).Mag() > 10*epsilon) {std::cout<<"pos mismatch\n"; pos.Print(); getPos.Print();}
		  if ((mom-getMom).Mag() > epsilon) {std::cout<<"mom mismatch\n"; mom.Print(); getMom.Print();}
		  rep->getPosMom(referencePlane, getPos, getMom);
		  if ((pos-getPos).Mag() > 10*epsilon) {std::cout<<"pos mismatch 2\n"; pos.Print(); getPos.Print();}
		  if ((mom-getMom).Mag() > epsilon) {std::cout<<"mom mismatch 2\n"; mom.Print(); getMom.Print();}

		  // calculate pulls  
		  TVectorT<double> state(rep->getState());
		  TMatrixTSym<double> cov(rep->getCov());
		  
		  if (debug) {
        state.Print();
        cov.Print();
		  }

      hmomRes->Fill( (charge/state[0]-momentum));
      hupRes->Fill(  (state[1]-referenceState[1]));
      hvpRes->Fill(  (state[2]-referenceState[2]));
      huRes->Fill(   (state[3]-referenceState[3]));
      hvRes->Fill(   (state[4]-referenceState[4]));

		  if (cov[0][0]>0) {
		    hqopPu->Fill( (state[0]-referenceState[0]) / sqrt(cov[0][0]) );
        pVal->Fill(   rep->getPVal());
        hupPu->Fill(  (state[1]-referenceState[1]) / sqrt(cov[1][1]) );
        hvpPu->Fill(  (state[2]-referenceState[2]) / sqrt(cov[2][2]) );
        huPu->Fill(   (state[3]-referenceState[3]) / sqrt(cov[3][3]) );
        hvPu->Fill(   (state[4]-referenceState[4]) / sqrt(cov[4][4]) );
		  }
		  
		  // print covariance
			if (debug) cov.Print();
			else if((iEvent)%100==0)  cov.Print();
				          
				          
		  // delete hits // todo: doesn't work; until then, we have a memleak
		  /*for(unsigned int i=0; i<npoints; ++i){
		    delete hits[i];
		  }*/
		  //delete rephits;
		  //delete rep;


			tree->Fill();

				          
				          
  }// end loop over events
#ifndef VALGRIND
  if (debug) std::cerr<<"Write Tree ...";
  tree->Write();
  if (debug) std::cerr<<"... done"<<std::endl;

  if (debug) std::cerr<<"Draw histograms ...";
	// fit and draw histograms
  TCanvas* c1 = new TCanvas();
  c1->Divide(2,3);

  c1->cd(1);
  hmomRes->Fit("gaus");
  hmomRes->Draw();

  c1->cd(3);
  hupRes->Fit("gaus");
  hupRes->Draw();

  c1->cd(4);
  hvpRes->Fit("gaus");
  hvpRes->Draw();

  c1->cd(5);
  huRes->Fit("gaus");
  huRes->Draw();

  c1->cd(6);
  hvRes->Fit("gaus");
  hvRes->Draw();

  c1->Write();

	TCanvas* c2 = new TCanvas();
	c2->Divide(2,3);
	
	c2->cd(1);
	hqopPu->Fit("gaus");
	hqopPu->Draw();
	
	c2->cd(2);
	pVal->Fit("pol1");
	pVal->Draw();
	c2->cd(3);
	hupPu->Fit("gaus");
	hupPu->Draw();
		
	c2->cd(4);
	hvpPu->Fit("gaus");
	hvpPu->Draw();
		
	c2->cd(5);
	huPu->Fit("gaus");
	huPu->Draw();	
		
	c2->cd(6);
	hvPu->Fit("gaus");
	hvPu->Draw();
	
	c2->Write();

	if (debug) std::cerr<<"... done"<<std::endl;
	
	// open event display
	display->setOptions("THDSPMA");
	display->open();
	
  rootapp->Run();		      

	file->Close();
  if (debug) std::cerr<<"... closed file"<<std::endl;
#endif
}

