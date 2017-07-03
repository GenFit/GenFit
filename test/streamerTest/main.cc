//
// Write fit results to tree, read again, compare.
//

#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>

#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <TGeoMaterialInterface.h>

#include <EventDisplay.h>

#include <HelixTrackModel.h>
#include <MeasurementCreator.h>

#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TRandom.h>
#include <TVector3.h>
#include <vector>

#include "TDatabasePDG.h"
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>


#define FILENAME "/tmp/streamerTest.root"

constexpr bool verbose = false;

bool emptyTrackTest()
{
  TFile *f = TFile::Open(FILENAME, "RECREATE");
  f->cd();
  genfit::Track *t = new genfit::Track();
  try {
    t->checkConsistency();
  } catch (genfit::Exception) {
    return false;
  }

  t->Write("direct");
  f->Close();
  delete t;
  delete f;

  f = TFile::Open(FILENAME, "READ");
  t = (genfit::Track*)f->Get("direct");

  bool result = false;
  try {
    t->checkConsistency();
    result = true;
  } catch (genfit::Exception) {
    result = false;
  }
  delete t;
  delete f;
  return result;
}


int main() {
  if (!emptyTrackTest()) {
    std::cout << "emptyTrackTest failed." << std::endl;
    return 1;
  }


  // prepare output tree for Tracks 
  // std::unique_ptr<genfit::Track> fitTrack(new genfit::Track());
  genfit::Track* fitTrack = new genfit::Track();
  TVectorD stateFinal;
  TMatrixDSym covFinal;
  genfit::DetPlane planeFinal;

  TFile* fOut = new TFile(FILENAME, "RECREATE");
  fOut->cd();
  TTree* tResults = new TTree("tResults", "results from track fit");
  // it is important to set the splitLevel to -1, because some of the genfit classes use custom streamers
  tResults->Branch("gfTrack", "genfit::Track", &fitTrack, 32000, -1);
  tResults->Branch("stateFinal", &stateFinal);
  tResults->Branch("covFinal", &covFinal, 32000, -1);
  tResults->Branch("planeFinal", &planeFinal, 32000, -1);





  gRandom->SetSeed(14);

  // init MeasurementCreator
  genfit::MeasurementCreator measurementCreator;


  // init geometry and mag. field
  new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0., 15.)); // 15 kGauss
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());


  // init fitter
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

  unsigned int nEvents(100);

  // main loop
  for (unsigned int iEvent=0; iEvent<nEvents; ++iEvent){

    // true start values
    TVector3 pos(0, 0, 0);
    TVector3 mom(1.,0,0);
    mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
    mom.SetTheta(gRandom->Uniform(0.4*TMath::Pi(),0.6*TMath::Pi()));
    mom.SetMag(gRandom->Uniform(0.2, 1.));


    // helix track model
    const int pdg = 13;               // particle pdg code
    const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge()/(3.);
    genfit::HelixTrackModel* helix = new genfit::HelixTrackModel(pos, mom, charge);
    measurementCreator.setTrackModel(helix);


    unsigned int nMeasurements = gRandom->Uniform(5, 15);


    // smeared start values
    const bool smearPosMom = true;     // init the Reps with smeared pos and mom
    const double posSmear = 0.1;     // cm
    const double momSmear = 3. /180.*TMath::Pi();     // rad
    const double momMagSmear = 0.1;   // relative

    TVector3 posM(pos);
    TVector3 momM(mom);
    if (smearPosMom) {
      posM.SetX(gRandom->Gaus(posM.X(),posSmear));
      posM.SetY(gRandom->Gaus(posM.Y(),posSmear));
      posM.SetZ(gRandom->Gaus(posM.Z(),posSmear));

      momM.SetPhi(gRandom->Gaus(mom.Phi(),momSmear));
      momM.SetTheta(gRandom->Gaus(mom.Theta(),momSmear));
      momM.SetMag(gRandom->Gaus(mom.Mag(), momMagSmear*mom.Mag()));
    }
    // approximate covariance
    TMatrixDSym covM(6);
    double resolution = 0.01;
    for (int i = 0; i < 3; ++i)
      covM(i,i) = resolution*resolution;
    for (int i = 3; i < 6; ++i)
      covM(i,i) = pow(resolution / nMeasurements / sqrt(3), 2);


    // trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

    // smeared start state
    genfit::MeasuredStateOnPlane stateSmeared(rep);
    rep->setPosMomCov(stateSmeared, posM, momM, covM);


    // create track
    TVectorD seedState(6);
    TMatrixDSym seedCov(6);
    rep->get6DStateCov(stateSmeared, seedState, seedCov);
    fitTrack = new genfit::Track(rep, seedState, seedCov);


    // create random measurement types
    std::vector<genfit::eMeasurementType> measurementTypes;
    for (unsigned int i = 0; i < nMeasurements; ++i)
      measurementTypes.push_back(genfit::eMeasurementType(gRandom->Uniform(8)));


    // create smeared measurements and add to track
    try{
      for (unsigned int i=0; i<measurementTypes.size(); ++i){
        std::vector<genfit::AbsMeasurement*> measurements = measurementCreator.create(measurementTypes[i], i*5.);
        genfit::TrackPoint* tp = new genfit::TrackPoint(measurements, fitTrack);
        // test scatterers
        genfit::ThinScatterer* sc = new genfit::ThinScatterer(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(1,1,1), TVector3(1,1,1))),
                                                              genfit::Material(1,2,3,4,5));
        tp->setScatterer(sc);

        fitTrack->insertPoint(tp);
      }
    }
    catch(genfit::Exception& e){
      std::cerr<<"Exception, next track"<<std::endl;
      std::cerr << e.what();
      delete fitTrack;
      fitTrack = 0;
      continue;
    }

    fitTrack->checkConsistency();

    // do the fit
    fitter->processTrack(fitTrack);

    fitTrack->checkConsistency();


    stateFinal.ResizeTo(fitTrack->getFittedState().getState());
    stateFinal = fitTrack->getFittedState().getState();
    covFinal.ResizeTo(fitTrack->getFittedState().getCov());
    covFinal = fitTrack->getFittedState().getCov();
    planeFinal =  *(fitTrack->getFittedState().getPlane());

    switch (iEvent % 5) {
    case 0:
      break;
    case 1:
      fitTrack->prune("FL");
      break;
    case 2:
      fitTrack->prune("W");
      break;
    case 3:
      fitTrack->prune("RC");
      break;
    case 4:
      fitTrack->prune("U");
      break;
    }

    tResults->Fill();

    //fitTrack->Print();

    delete fitTrack;
    fitTrack = 0;

  }// end loop over events

  delete fitter;




  fOut->Write();
  delete fOut;

  fOut = TFile::Open(FILENAME, "READ");
  fOut->GetObject("tResults", tResults);
  TVectorD* pState = 0;
  tResults->SetBranchAddress("stateFinal", &pState);
  TMatrixDSym* pMatrix = 0;
  tResults->SetBranchAddress("covFinal", &pMatrix);
  genfit::DetPlane* plane = 0;
  tResults->SetBranchAddress("planeFinal", &plane);
  tResults->SetBranchAddress("gfTrack", &fitTrack);

  int fail(0);

  for (Long_t nEntry = 0; nEntry < tResults->GetEntries(); ++nEntry) {
    tResults->GetEntry(nEntry);

    try {
      fitTrack->checkConsistency();
    } catch (genfit::Exception& e) {
      std::cout << e.getExcString() << std::endl;
      return 1;
    }

    try {
      if (*pState == fitTrack->getFittedState().getState() &&
          *pMatrix == fitTrack->getFittedState().getCov() &&
          *plane ==  *(fitTrack->getFittedState().getPlane())) {
        // track ok
      }
      else {
        if (verbose) {
          std::cout << "stored track not equal, small differences can occur if some info has been pruned." << std::endl;
          pState->Print();
          fitTrack->getFittedState().getState().Print();
          pMatrix->Print();
          fitTrack->getFittedState().getCov().Print();
          plane->Print();
          fitTrack->getFittedState().getPlane()->Print();
          }
        ++fail;
        //return 1;
      }
    }
    catch (genfit::Exception& e) {
        if (verbose) {
            std::cerr << e.what();
        }
      return 1;
    }
  }
  std::cout << nEvents - fail << " stored tracks are identical to fitted tracks, as far as tested." << std::endl;
  std::cout << fail << " tracks were not identical to fitted tracks, as far as tested." << std::endl;
  delete fitTrack;
  std::cout << "deleteing didn't segfault" << std::endl;

  return 0;
}
