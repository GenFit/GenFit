#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackCand.h>
#include <TrackPoint.h>

#include <MeasurementProducer.h>
#include <MeasurementFactory.h>

#include "mySpacepointDetectorHit.h"
#include "mySpacepointMeasurement.h"

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




int main() {

  gRandom->SetSeed(14);

  // init MeasurementCreator
  genfit::MeasurementCreator measurementCreator;


  // init geometry and mag. field
  new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0.,0., 15.)); // 15 kGauss
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());


  // init event display
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();


  // init fitter
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();


  TClonesArray myDetectorHitArray("genfit::mySpacepointDetectorHit");

  // init the factory
  int myDetId(1);
  genfit::MeasurementFactory<genfit::AbsMeasurement> factory;
  genfit::MeasurementProducer<genfit::mySpacepointDetectorHit, genfit::mySpacepointMeasurement> myProducer(&myDetectorHitArray);
  factory.addProducer(myDetId, &myProducer);


  // main loop
  for (unsigned int iEvent=0; iEvent<100; ++iEvent){

    myDetectorHitArray.Clear();

    //TrackCand
    genfit::TrackCand myCand;

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

    // covariance
    double resolution = 0.01;
    TMatrixDSym cov(3);
    for (int i = 0; i < 3; ++i)
      cov(i,i) = resolution*resolution;

    for (unsigned int i=0; i<nMeasurements; ++i) {
      // "simulate" the detector
      TVector3 currentPos = helix->getPos(i*2.);
      currentPos.SetX(gRandom->Gaus(currentPos.X(), resolution));
      currentPos.SetY(gRandom->Gaus(currentPos.Y(), resolution));
      currentPos.SetZ(gRandom->Gaus(currentPos.Z(), resolution));

      // Fill the TClonesArray and the TrackCand
      // In a real experiment, you detector code would deliver mySpacepointDetectorHits and fill the TClonesArray.
      // The patternRecognition would create the TrackCand.
      new(myDetectorHitArray[i]) genfit::mySpacepointDetectorHit(currentPos, cov);
      myCand.addHit(myDetId, i);
    }


    // smeared start values (would come from the pattern recognition)
    const bool smearPosMom = true;   // init the Reps with smeared pos and mom
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


    // set start values and pdg to cand
    myCand.setPosMomSeedAndPdgCode(posM, momM, pdg);


    // create track
    genfit::Track fitTrack(myCand, factory, new genfit::RKTrackRep(pdg));


    // do the fit
    try{
      fitter->processTrack(&fitTrack);
    }
    catch(genfit::Exception& e){
      std::cerr << e.what();
      std::cerr << "Exception, next track" << std::endl;
      continue;
    }

    //check
    assert(fitTrack.checkConsistency());


    if (iEvent < 1000) {
      // add track to event display
      display->addEvent(&fitTrack);
    }


  }// end loop over events

  delete fitter;

  // open event display
  display->open();

}


