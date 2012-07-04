
#include <iostream>
#include <GenfitDisplay.h>
#include <GFConstField.h>
#include <GFException.h>
#include <GFFieldManager.h>
#include <GFTrack.h>
#include <RKTrackRep.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TH1D.h>
#include <TH3D.h>
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

#include "GFRaveVertexFactory.h"
#include "GFRaveVertex.h"


//#define VALGRIND

int main() {
  std::cerr<<"main"<<std::endl;

  const unsigned int nEvents = 1000;
  const double BField = 15.;       // kGauss
  const unsigned int nTracks = 3;
  const TVector3 vertexPos(0,0,0);
  const double vertexErr = 0.1;
  std::string fMethod("kalman-smoothing:1");
  const bool fUseBeamspot = false;

  const bool useReps = false; // false = use GFTracks as input; true = use GFTrackReps as input

  const double mom = 0.75;

  const bool debug = false;




  Int_t fVerbose = 0;
  if (debug) fVerbose=10;

  // init mersenne twister with TUUID
  TRandom3 rand(10);

  // init geometry and mag. field
  TGeoManager* geom = new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  GFFieldManager::getInstance()->init(new GFConstField(0.,0.,BField));

  TMatrixD fBeamCov(3,3);
  fBeamCov.UnitMatrix();
  fBeamCov *= vertexErr;

  // init vertex factory
  GFRaveVertexFactory* fVertexFactory = new GFRaveVertexFactory(fVerbose, false);
  fVertexFactory->setBeamspot(vertexPos, fBeamCov);
  fVertexFactory->setMethod(fMethod);

  std::vector < GFRaveVertex* >* fVertexBuffer = new std::vector < GFRaveVertex* >;

  // init rootapp (for drawing histograms)
  TApplication* rootapp = new TApplication("rootapp", 0, 0);

  // create histograms
  TH3D *vertexRes = new TH3D("vertexRes","vertexRes",100, vertexPos.X()-5.*vertexErr, vertexPos.X()+5.*vertexErr,
                                                     100, vertexPos.Y()-5.*vertexErr, vertexPos.Y()+5.*vertexErr,
                                                     100, vertexPos.Z()-5.*vertexErr, vertexPos.Z()+5.*vertexErr);

  // main loop
  for (unsigned int iEvent=0; iEvent<nEvents; ++iEvent){

    if (debug || (iEvent+1)%10==0) std::cout << iEvent+1 << std::endl;

    std::vector<GFTrack*> tracks;
    std::vector<GFAbsTrackRep*> reps;

    // create and fill track(rep)s
    for (unsigned int iTr=0; iTr<nTracks; ++iTr){

      TVector3 trackPos(vertexPos);
      trackPos.SetXYZ(rand.Gaus(trackPos.X(), vertexErr),
                      rand.Gaus(trackPos.Y(), vertexErr),
                      rand.Gaus(trackPos.Z(), vertexErr));

      TVector3 trackMom(mom,0,0);
      trackMom.SetPhi(rand.Uniform(0.,2*TMath::Pi()));
      trackMom.SetTheta(rand.Uniform(0.,TMath::Pi()));

      RKTrackRep* rep = new RKTrackRep(trackPos, trackMom, TVector3(vertexErr,vertexErr,vertexErr), TVector3(vertexErr,vertexErr,vertexErr), 211);
      GFTrack* track = new GFTrack(rep, false);

      tracks.push_back(track);
      reps.push_back(rep);

    }

    if (!useReps){
      // method A
      fVertexFactory->findVertices(fVertexBuffer, tracks, fUseBeamspot);
    }
    else {
      // method B
      fVertexFactory->findVertices(fVertexBuffer, reps, fUseBeamspot);
    }

    // test everything
    unsigned int nVert = fVertexBuffer->size();

    for (unsigned int i=0; i<nVert; ++i){
      if (debug) (*fVertexBuffer)[i]->Print();
      vertexRes->Fill((*fVertexBuffer)[i]->getPos().X(),
                      (*fVertexBuffer)[i]->getPos().Y(),
                      (*fVertexBuffer)[i]->getPos().Z());
    }


    // delete tracks, reps and vertices
    for (unsigned int i=0; i<nTracks; ++i) delete tracks[i];
    for (unsigned int i=0; i<nVert; ++i) delete (*fVertexBuffer)[i];
    fVertexBuffer->clear();

  }// end of main loop


  delete fVertexFactory;
  delete fVertexBuffer;


#ifndef VALGRIND

  if (debug) std::cerr<<"Draw histograms ...";
  // fit and draw histograms
  TCanvas* c1 = new TCanvas();


  c1->cd(1);
  vertexRes->Draw();


  c1->Write();



  if (debug) std::cerr<<"... done"<<std::endl;

  // open event display
  //display->setOptions("THDSPMAG"); // G show geometry
  //display->open();

  rootapp->Run();

#endif
}

