#include <iostream>
#include <assert.h>
#include <math.h>
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "TRandom.h"
#include "TClonesArray.h"
#include "TGeoManager.h"

#include "RKTrackRep.h"
//#include "GeaneTrackRep2.h"
#include "GFTrack.h"
#include "GFKalman.h"
#include "PointHit.h"
#include "PixHit.h"

#include "GFConstField.h"
#include "GFFieldManager.h"
#include "GFException.h"

#define PDGID 13

#define DIMREP 5

int main(int argc, char** argv){



  assert(argc==4);

  TFile* geoFile = TFile::Open(argv[1]);
  assert(geoFile!=NULL);
  TGeoManager *geo = (TGeoManager*)(gROOT->FindObject("FAIRGeom"));
  //  gROOT->Macro("Geane.C");

  TFile* inFile = TFile::Open(argv[2]);
  assert(inFile!=NULL);
  TTree *inTree = (TTree*)(gROOT->FindObject("cbmsim"));

  TFile *outFile = new TFile(argv[3],"RECREATE");
  TTree* t = new TTree("t","ilc study output");
  
  float stMCT[DIMREP];
  float stREC[DIMREP];
  float covREC[DIMREP];

  char buf[50];
  sprintf(buf,"stMCT[%d]/F",DIMREP);
  t->Branch("stMCT",stMCT,buf);
  sprintf(buf,"stREC[%d]/F",DIMREP);
  t->Branch("stREC",stREC,buf);
  sprintf(buf,"covREC[%d]/F",DIMREP);
  t->Branch("covREC",covREC,buf);


  //initialize constant magnetic field
  //the real field map class (if present) would be instantiated here as well
  GFFieldManager::getInstance()->init(new GFConstField(0.,0.,35.));

  gRandom->SetSeed(1);

  TClonesArray *silPoints=0;
  TClonesArray *tpcPoints=0;
  TClonesArray *MCTpos=0;
  TClonesArray *MCTmom=0;
  inTree->SetBranchAddress("ilcPoints",&silPoints);
  inTree->SetBranchAddress("ilcTpcPoints",&tpcPoints);
  inTree->SetBranchAddress("ilcMCTpos",&MCTpos);
  inTree->SetBranchAddress("ilcMCTmom",&MCTmom);

  GFKalman kalman;

  TH1D *momRec = new TH1D("momRec","reconstructed momentum",200,49.,51.);

  int nEv = inTree->GetEntries();
  //nEv=1;
  for(int iev=0;iev<nEv;++iev){
    inTree->GetEntry(iev);
    std::cout << "########## processing event number " << iev << std::endl;

    //extract MC truth values for track pos and mom from tree
    assert(MCTpos->GetEntriesFast()==1);
    assert(MCTmom->GetEntriesFast()==1);
    TVector3 mctPos = *((TVector3*)MCTpos->At(0));
    TVector3 mctMom = *((TVector3*)MCTmom->At(0));

    std::cout << "MC truth posision and momentum:" << std::endl;
    mctPos.Print();
    mctMom.Print();
    
    //calculate seed values for position and momentum, i.e. smear them
    TVector3 seedPos = mctPos;
    seedPos.SetX( gRandom->Gaus(seedPos.X(),0.1) );
    seedPos.SetY( gRandom->Gaus(seedPos.Y(),0.1) );
    seedPos.SetZ( gRandom->Gaus(seedPos.Z(),0.1) );
    TVector3 seedMom = mctMom;
    seedMom.SetX( gRandom->Gaus(seedMom.X(),0.1*seedMom.X()) );
    seedMom.SetY( gRandom->Gaus(seedMom.Y(),0.1*seedMom.Y()) );
    seedMom.SetZ( gRandom->Gaus(seedMom.Z(),0.1*seedMom.Z()) );

    //TVector3 err(1.,1.,1.);
    //GFDetPlane startPlane(seedPos,seedMom);
    //GFAbsTrackRep* fitRep = new GeaneTrackRep2(startPlane,seedMom,err,err,PDGID);//pdg id #define from simulation
    RKTrackRep* fitRep = new RKTrackRep(seedPos,seedMom,PDGID);//pdg id #define from simulation
    //GFTrack takes ownership over trackRep
    GFTrack fitTrack(fitRep);

    int nSilPoints = silPoints->GetEntriesFast();
    int silHitCounter(0);
    std::cout << "number of hits in silicon detectors: " << nSilPoints << std::endl;
    for(int i=0;i<nSilPoints;++i){
      TVector3* P = (TVector3*)(silPoints->At(i));
      fitTrack.addHit(new PixHit(*P),
		      1,//some dummy det ID
		      P->Mag(),//ordering parameter (is not used in this case)
		      silHitCounter++);//plane ID (set dummy values for TPC >100)
      
    }

    int nTpcPoints = tpcPoints->GetEntriesFast();
    int tpcHitCounter(0);
    std::cout << "number of hits in tpc detector: " << nTpcPoints << std::endl;
    for(int i=0;i<nTpcPoints;++i){
      TVector3* P = (TVector3*)(tpcPoints->At(i));
      double R = sqrt(P->X()*P->X()+P->Y()*P->Y());
      if(R<39.5 || R>173.9) continue;//skip hits that are outside the pad plane radius
      
      fitTrack.addHit(new PointHit(*P),
		      2,//some dummy det ID
		      P->Mag(),//ordering parameter (is not used in this case)
		      (tpcHitCounter++) + 100);//plane ID (set dummy values for TPC >100)
      
    }

    kalman.processTrack(&fitTrack);
    std::cout << "status flag: " << fitRep->getStatusFlag() << std::endl;
    if(fitRep->getStatusFlag()!=0) continue;
    //fit was successful
    std::cout << "resulting momentum vector:" << std::endl;
    fitTrack.getMom().Print();
    std::cout << "number of degrees of freedom :" << fitRep->getNDF() << std::endl;
    momRec->Fill(    fitTrack.getMom().Mag() );

    //define the reference plane where MC truth state and reconstructed state are compared
    TVector3 refU(0.,0.,1.);
    TVector3 refNormal(mctMom.X(),mctMom.Y(),0.);//normal vector of refPlane
    refNormal.SetMag(1.);
    TVector3 refV = refNormal.Cross(refU);


    GFDetPlane refPlane(mctPos,refU,refV);
    TMatrixT<double> stREF(DIMREP,1);
    TMatrixT<double> covREF(DIMREP,DIMREP);
    try{
      //this operation does not modify the track representation, although
      //it would matter, because we are done with it anyways.
      fitRep->extrapolate(refPlane,stREF,covREF);
    }
    catch(GFException& e){
      std::cerr << e.what();
      std::cerr << "### excpetion occured in extrapolation to reference plane -> skip event" << std::endl;
      continue;
    }
    for(int i=0;i<DIMREP;++i){
      //save the output in standard arrays, because saving TMatrixT<double> to the tree
      //is a pain in the ass when it comes to TTree::Draw()
      stREC[i] = stREF[i][0];
      covREC[i] = sqrt(covREF[i][i]);
    }

    //calculate MC truth state vector for output tree

    //we define the reference plain to have origin in MC truth pos, so:
    stMCT[3]=0.;//u
    stMCT[4]=0.;//v
    //need the particle charge to calculate the MC truth value for charge
    TParticlePDG * part = TDatabasePDG::Instance()->GetParticle(PDGID);//define at beginning of file
    double charge = part->Charge()/(3.);
    stMCT[0] = charge/mctMom.Mag();
    //calculate du/dw and dv/dw
    double pu=mctMom*refU;
    double pv=mctMom*refV;
    double pw=mctMom*refNormal;
    assert(fabs(pw)>1.E-4); //this can never happen by construction
    stMCT[1] = pu/pw;
    stMCT[2] = pv/pw;

    //fill the output tree
    t->Fill();

  }
  inFile->Close();
  outFile->cd();
  t->Write();
  momRec->Write();
  outFile->Close();
}
