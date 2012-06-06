
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
#include <GFRectFinitePlane.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TVector3.h>
#include <vector>

#include <TROOT.h>

#include <TVirtualMC.h>

#include "PixHit.h"
#include "PointHit.h"

  typedef double M1x3[1*3];
  typedef double M1x4[1*4];
  typedef double M1x7[1*7];
  typedef double M5x5[5*5];
  typedef double M6x6[6*6];
  typedef double M7x7[7*7];
  typedef double M8x7[8*7];
  typedef double M6x5[6*5];
  typedef double M7x5[7*5];
  typedef double M5x6[5*6];
  typedef double M5x7[5*7];

int main() {

  const unsigned int nEvents = 100;
  const double BField = 20.;       // kGauss
  const double momentum = 0.6218547;     // GeV

  const bool debug = false;



  // init mersenne twister with TUUID
	TRandom3 rand(0);
	
  // init event display
	GenfitDisplay* display = GenfitDisplay::getInstance();
	display->reset();
  
  // init geometry and mag. field
  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  GFFieldManager::getInstance()->init(new GFConstField(0.,0.,BField));
  
  // init rootapp (for drawing histograms)
  TApplication* rootapp = new TApplication("rootapp", 0, 0);

  
  // create histograms
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptFit(1111);
  

  TH1D *diffs = new TH1D("diffs","diffs",2000,-1.,1.);
  TH1D *diffs2 = new TH1D("diffs2","diffs2",2000,-1.E15,1.E15);

  TH2D *phaseSpace1 = new TH2D("phaseSpace", "phaseSpace 7D direct -> ERROR", 360,-1.*TMath::Pi(),TMath::Pi(), 180,0,TMath::Pi() );
  TH2D *phaseSpaceBar1 = new TH2D("phaseSpaceBar", "phaseSpaceBar -> all ok", 360,-1.*TMath::Pi(),TMath::Pi(), 180,0,TMath::Pi() );

  TH2D *phaseSpace2 = new TH2D("phaseSpace", "phaseSpace 6D direct -> ERROR", 360,-1.*TMath::Pi(),TMath::Pi(), 180,0,TMath::Pi() );
  TH2D *phaseSpaceBar2 = new TH2D("phaseSpaceBar", "phaseSpaceBar -> all ok", 360,-1.*TMath::Pi(),TMath::Pi(), 180,0,TMath::Pi() );

  TH2D *phaseSpace3 = new TH2D("phaseSpace", "phaseSpace 6D get and set -> ERROR", 360,-1.*TMath::Pi(),TMath::Pi(), 180,0,TMath::Pi() );
  TH2D *phaseSpaceBar3 = new TH2D("phaseSpaceBar", "phaseSpaceBar -> all ok", 360,-1.*TMath::Pi(),TMath::Pi(), 180,0,TMath::Pi() );


  int method, dim, setting(0);
  bool sixD(false);

  // main loop
  for (unsigned int iEvent=0; iEvent<3*nEvents; ++iEvent){

    std::cout << "----------------------------------Event "<<iEvent<<"---------------------------------------------\n";

      if (iEvent==0){
        method = 0;
        dim = 7;
        sixD = false;
        std::cout << "                      =======================================================================================================\n\
                      ============================================== SETTING 1: 7D direct ===================================\n\
                      =======================================================================================================\n";
      }
      if (iEvent==nEvents){
        method = 0;
        dim = 6;
        sixD = true;
        ++setting;
        std::cout << "                      =======================================================================================================\n\
                      ============================================== SETTING 2: 6D direct ===================================\n\
                      =======================================================================================================\n";
      }
      if (iEvent==2*nEvents){
        method = 1;
        dim = 6;
        sixD = true;
        ++setting;
        std::cout << "                      =======================================================================================================\n\
                      ============================================== SETTING 3: 6D get and set ==============================\n\
                      =======================================================================================================\n";
      }

      if (debug || (iEvent+1)%10==0) std::cout << iEvent+1 << std::endl;

      // true start values
      TVector3 pos(0, 0, 0);
      TVector3 mom(1.,0,0);
      mom.SetMag(momentum);
      mom.SetPhi(rand.Uniform(0.,2*TMath::Pi()));
      mom.SetTheta(rand.Uniform(0*TMath::Pi(),1.0*TMath::Pi()));

      TVector3 posErr(1.,1.,1.);
      posErr *= 0.01;
      TVector3 momErr(1.,1.,1.);
      momErr *= 0.01;

      GFAbsTrackRep* rep = new RKTrackRep(pos, mom, /*posErr, momErr,*/ 211);
      				          
		  // extrapolate to reference plane.
      // ref origin
      TVector3 O(1,1,1);
      O.SetMag(rand.Uniform(5.,25.));
      O.SetPhi(rand.Uniform(0.,2*TMath::Pi()));
      O.SetTheta(rand.Uniform(0*TMath::Pi(),1.0*TMath::Pi()));
      // ref direction
      TVector3 N(1,0,0);
      N.SetPhi(rand.Uniform(0.,2*TMath::Pi()));
      N.SetTheta(rand.Uniform(0*TMath::Pi(),1.0*TMath::Pi()));

      GFDetPlane referencePlane(O, N);
      try{
        rep->extrapolate(referencePlane);
      }
      catch(GFException& e){
        std::cerr<<"Exception, next track"<<std::endl;
        e.what();
        continue; // here is a memleak!
      }


      TMatrixT<double> covBefore(5,5);
      covBefore = rep->getCov();
      TMatrixT<double> stateBefore(5,1);
      stateBefore = rep->getState();

      double state7[7];
      ((RKTrackRep*)rep)->getState7(state7);
      TMatrixT<double> state7Before(7,1, state7);
      M1x7& state7Before_ = *((M1x7*) state7Before.GetMatrixArray());
      GFDetPlane planeBefore(rep->getReferencePlane());
      TMatrixT<double> covTrans(dim, dim);
      TMatrixT<double> covAfter(5,5);

      if(method==0){
        // get covariance

        // transform
        TMatrixT<double> Jac_pM;
        double* covTrans_ = covTrans.GetMatrixArray();

        if(!sixD) ((RKTrackRep*)rep)->transformPM7(covBefore, *((M7x7*) covTrans_), planeBefore, stateBefore, ((RKTrackRep*)rep)->fSpu, &Jac_pM);
        else      ((RKTrackRep*)rep)->transformPM6(covBefore, *((M6x6*) covTrans_), planeBefore, stateBefore, ((RKTrackRep*)rep)->fSpu, &Jac_pM);

        // transform back
        TMatrixT<double> Jac_Mp;
        if(!sixD) ((RKTrackRep*)rep)->transformM7P(*((M7x7*) covTrans_), covAfter, planeBefore, state7Before_, &Jac_Mp);
        else      ((RKTrackRep*)rep)->transformM6P(*((M6x6*) covTrans_), covAfter, planeBefore, state7Before_, &Jac_Mp);

        std::cout << "fSpu = "<<((RKTrackRep*)rep)->fSpu<<"\n";

        // check if equal
        TMatrixT<double> covDiff(covAfter-covBefore);
        TMatrixT<double> JJ(Jac_Mp.T() * Jac_pM.T());
        //TMatrixT<double> JJ2(Jac_pM * Jac_Mp);

        double epsilon = 0.001;
        if (covDiff.Max() > epsilon || covDiff.Min() < -1.*epsilon){
          //std::cout << "product of jacobians Jac_pM * Jac_Mp";
          //JJ2.Print();

          std::cout << "product of jacobians Jac_Mp.T() * Jac_pM.T()";
          JJ.Print();

          /*std::cout << "J_pM^T ";
          Jac_pM.T().Print();
          std::cout << "J_Mp ";
          Jac_Mp.Print();*/
        }



        /*double epsilon = 1E-10;
        TMatrixT<double> unit2(JJ2);
        unit2.UnitMatrix();
        if ((JJ2-unit2).Max() > epsilon || (JJ2-unit2).Min() < -1.*epsilon){
        }

        // check Jacobians
        TMatrixT<double> unit(JJ);
        unit.UnitMatrix();
        if ((JJ-unit).Max() > epsilon || (JJ-unit).Min() < -1.*epsilon){

        }*/
      }
      else{
        // get pos mom
        ((RKTrackRep*)rep)->getPosMom(referencePlane, pos, mom);

        // extrapolate so that plane will be the same after setting posMomCov
        GFDetPlane pln(pos, mom);
        try{
          rep->extrapolate(pln);
        }
        catch(GFException& e){
          std::cerr<<"Exception, next track"<<std::endl;
          e.what();
          continue; // here is a memleak!
        }

        // get original Cov
        covBefore = rep->getCov();

        // get covariance
        ((RKTrackRep*)rep)->getPosMomCov(pln, pos, mom, covTrans);

        // set covariance
        ((RKTrackRep*)rep)->setPosMomCov(pos, mom, covTrans);

        // get cov after
        covAfter = rep->getCov();
      }
		  

		  
		  // check if equal
      TMatrixT<double> covDiff(covAfter-covBefore);

      double epsilon = 0.001;
      if (covDiff.Max() > epsilon || covDiff.Min() < -1.*epsilon){

        std::cout << "pos: "; pos.Print();
        std::cout << "mom: "; mom.Print();

        std::cout << "State";
        stateBefore.Print();
        state7Before.Print();
        std::cout << "Cov before: \n";
        covBefore.Print();
        std::cout << "Cov trans: \n";
        covTrans.Print();
        std::cout << "Cov after: \n";
        covAfter.Print();
        std::cout << "Cov diff: \n";
        covDiff.Print();

        switch (setting){
          case 0:
            phaseSpace1->Fill(mom.Phi(), mom.Theta());
            break;
          case 1:
            phaseSpace2->Fill(mom.Phi(), mom.Theta());
            break;
          case 2:
            phaseSpace3->Fill(mom.Phi(), mom.Theta());
        }
      }
      else {
        switch (setting){
          case 0:
            phaseSpaceBar1->Fill(mom.Phi(), mom.Theta());
            break;
          case 1:
            phaseSpaceBar2->Fill(mom.Phi(), mom.Theta());
            break;
          case 2:
            phaseSpaceBar3->Fill(mom.Phi(), mom.Theta());
        }
      }


				          
  }// end loop over events

  // fit and draw histograms
  TCanvas* c1 = new TCanvas();
  c1->Divide(2,3);

  int iCanv(0);

  c1->cd(++iCanv);
  phaseSpace1->Draw("colz");

  c1->cd(++iCanv);
  phaseSpaceBar1->Draw("colz");


  c1->cd(++iCanv);
  phaseSpace2->Draw("colz");

  c1->cd(++iCanv);
  phaseSpaceBar2->Draw("colz");


  c1->cd(++iCanv);
  phaseSpace3->Draw("colz");

  c1->cd(++iCanv);
  phaseSpaceBar3->Draw("colz");

  // open event display
  display->setOptions("THDSPM");
  display->open();

  rootapp->Run();
}

