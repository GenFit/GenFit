
#include <iostream>
#include <GenfitDisplay.h>
#include <GFConstField.h>
#include <GFException.h>
#include <GFFieldManager.h>
#include <GFKalman.h>
#include <GFTools.h>
#include <GFTrack.h>
#include <RKTrackRep.h>
#include <GeaneTrackRep2.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TH1D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TVector3.h>
#include <vector>

#include "PixHit.h"
#include "PointHit.h"

int main() {

  const unsigned int nEvents = 500; 
  const double momentum = 0.5;     // GeV
  const double posSmear = 0.1;     // cm
  const double momSmear = 0.1;     // GeV
  const unsigned int npoints = 30; // number of hits generated
  const double pointDist = 2.;     // cm; approx. distance between hits generated
  const double resolution = 0.06;  // cm; resolution of generated hits
  
  const bool planarHits = true;   // true: create planar hits; false: create space point hits
  const bool GEANE = false;       // true: use GeaneTrackRep2; false: use RKTrackRep

  // init fitter
  GFKalman kalman;
  kalman.setNumIterations(3);


  // init mersenne twister with TUUID
	TRandom3 rand(0);
	
  // init event display
	GenfitDisplay* display = GenfitDisplay::getInstance();
	display->reset();
  
  // init geometry and mag. field
  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  GFFieldManager::getInstance()->init(new GFConstField(0.,0.,15));
  
  // init rootapp (for drawing histograms)
  TApplication* rootapp = new TApplication("rootapp", 0, 0);
  
  // create histograms
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptFit(1111);
  
  TH1D *hqopPu = new TH1D("hqopPu","q/p pull",200,-6.,6.);
  TH1D *hupPu = new TH1D("hupPu","u' pull",200,-6.,6.);
  TH1D *hvpPu = new TH1D("hvpPu","v' pull",200,-6.,6.);
  TH1D *huPu = new TH1D("huPu","u pull",200,-6.,6.);
  TH1D *hvPu = new TH1D("hvPu","v pull",200,-6.,6.);

  // main loop
  for (unsigned int iEvent=0; iEvent<nEvents; ++iEvent){
      
      if((iEvent+1)%10==0) std::cout << iEvent+1 << std::endl;

      // true start values
      TVector3 pos(1.,0,0);
      TVector3 mom(1.,1.,1.);
      mom.SetMag(momentum);
      mom.SetPhi(rand.Uniform(0.,2*TMath::Pi()));
      mom.SetTheta(rand.Uniform(0.2*TMath::Pi(),0.8*TMath::Pi()));

      TVector3 posErr(1.,1.,1.);
      posErr *= 2.*posSmear;
      TVector3 momErr(1.,1.,1.);
      momErr *= 2.*momSmear;
      
      // smeared start values
      TVector3 posM(pos);
      posM.SetX(rand.Gaus(posM.X(),posSmear*posM.X()));
      posM.SetY(rand.Gaus(posM.Y(),posSmear*posM.Y()));
      posM.SetZ(rand.Gaus(posM.Z(),posSmear*posM.Z()));

      TVector3 momM(mom);
      momM.SetX(rand.Gaus(momM.X(),momSmear*momM.X()));
      momM.SetY(rand.Gaus(momM.Y(),momSmear*momM.Y()));
      momM.SetZ(rand.Gaus(momM.Z(),momSmear*momM.Z())); 
      

        
      // trackrep for creating hits
      GFAbsTrackRep* rephits;
      if (!GEANE) rephits = new RKTrackRep(pos, mom, posErr, momErr, 211);
      else rephits = new GeaneTrackRep2(GFDetPlane(pos, mom), mom, posErr, momErr, 211);
      					      
      // remember original initial plane and state					  
      GFDetPlane referencePlane(rephits->getReferencePlane());
      TMatrixT<double> referenceState(rephits->getState());
      
      // create smeared hits
      std::vector<GFAbsRecoHit*> hits;
        
      TVector3 point, dir;
      try{
        for (unsigned int i=0; i<npoints; ++i){
          // get current position and momentum
          rephits->getPosMom(rephits->getReferencePlane(), point, dir);
          // create hit
          if (planarHits) {
            PixHit* hit = new PixHit(GFDetPlane(point, dir), resolution);
            hits.push_back(hit);
          }
          else {
            TVector3 smearedPos(point);
            smearedPos.SetX(rand.Gaus(smearedPos.X(),resolution));
            smearedPos.SetY(rand.Gaus(smearedPos.Y(),resolution));
            smearedPos.SetZ(rand.Gaus(smearedPos.Z(),resolution)); 
            PointHit* hit = new PointHit(smearedPos, TVector3(resolution,resolution,resolution));
            hits.push_back(hit);
          }
          // stepalong (approximately)
          dir.SetMag(pointDist);
          GFDetPlane pl(point+dir, dir);
          rephits->extrapolate(pl);
        }
      }
      catch(GFException& e){
	      e.what();
	      std::cerr<<"Exception, next track"<<std::endl;
	      continue; // here is a memleak!
      } 
      
      
      
      
      // trackrep to be fitted and tested
      GFAbsTrackRep* rep;
      if (!GEANE) rep = new RKTrackRep(posM, momM, posErr, momErr, 211);
      else rep = new GeaneTrackRep2(GFDetPlane(posM, momM), momM, posErr, momErr, 211);
      				          
		  // create track, add hits 
		  GFTrack fitTrack(rep); //initialized with smeared rep
      for(unsigned int i=0; i<hits.size(); ++i){
        fitTrack.addHit(hits[i], 
		                    3,//dummy detector id
		                    i);
      }
		  
		  
		  // do the fit
      try{
        kalman.processTrack(&fitTrack);
      }
      catch(GFException& e){
        e.what();
	      std::cerr<<"Exception, next track"<<std::endl;
	      continue; // here is a memleak!
      }
      // check if fit was successfull
      if(rep->getStatusFlag() != 0 ) continue; // here is a memleak!


      // add track to event display
      std::vector<GFTrack*> event;
			event.push_back(&fitTrack);
			display->addEvent(event);

				          
		  // extrapolate back to reference plane
		  rep->extrapolate(referencePlane);
		  
		  // calculate pulls  
		  TMatrixT<double> state(rep->getState());
		  TMatrixT<double> cov(rep->getCov());
		  
		  hqopPu->Fill( (state[0][0]-referenceState[0][0]) / sqrt(cov[0][0]) );
			hupPu->Fill(  (state[1][0]-referenceState[1][0]) / sqrt(cov[1][1]) );
			hvpPu->Fill(  (state[2][0]-referenceState[2][0]) / sqrt(cov[2][2]) );
			huPu->Fill(   (state[3][0]-referenceState[3][0]) / sqrt(cov[3][3]) );
			hvPu->Fill(   (state[4][0]-referenceState[4][0]) / sqrt(cov[4][4]) );            
		  
		  // print covariance        
		  if((iEvent)%100==0)  cov.Print();
				          
				          
		  // delete hits // todo: doesn't work; until then, we have a memleak
		  /*for(unsigned int i=0; i<npoints; ++i){
		    delete hits[i];
		  }*/
		  //delete rephits;
		  //delete rep;

				          
				          
  }// end loop over events
  
				      
	// fit and draw histograms		
	TCanvas* c1 = new TCanvas();
	c1->Divide(2,3);
	
	c1->cd(1);
	hqopPu->Fit("gaus");
	hqopPu->Draw();
	
	c1->cd(3);
	hupPu->Fit("gaus");
	hupPu->Draw();
		
	c1->cd(4);
	hvpPu->Fit("gaus");
	hvpPu->Draw();
		
	c1->cd(5);
	huPu->Fit("gaus");
	huPu->Draw();	
		
	c1->cd(6);
	hvPu->Fit("gaus");
	hvPu->Draw();
	
	
	// open event display
	display->setOptions("THDSPM");
	display->open();
	
  rootapp->Run();		      
				      
}

