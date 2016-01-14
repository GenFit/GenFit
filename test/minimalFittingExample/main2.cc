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

#include <PlanarMeasurement.h>

#include <TEveManager.h>
#include <TGeoManager.h>
#include <TVector3.h>
#include <vector>

#include "TDatabasePDG.h"
#include <TMath.h>




int main() {

  // init geometry and mag. field
  new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0. ,10., 0.)); // 1 T


  // init event display
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();


  // init fitter
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();


  // particle pdg code; pion hypothesis
  const int pdg = 211;

  // start values for the fit, e.g. from pattern recognition
  TVector3 pos(0, 0, 0);
  TVector3 mom(0, 0, 3);


  // trackrep
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

  // create track
  genfit::Track fitTrack(rep, pos, mom);


  const int detId(0); // detector ID
  int planeId(0); // detector plane ID
  int hitId(0); // hit ID

  double detectorResolution(0.001); // resolution of planar detectors
  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
  hitCov *= detectorResolution*detectorResolution;


  // add some planar hits to track with coordinates I just made up
  TVectorD hitCoords(2);
  hitCoords[0] = 0;
  hitCoords[1] = 0;
  genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
  measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,0), TVector3(1,0,0), TVector3(0,1,0))), ++planeId);
  fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

  hitCoords[0] = -0.15;
  hitCoords[1] = 0;
  measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
  measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,10), TVector3(1,0,0), TVector3(0,1,0))), ++planeId);
  fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

  hitCoords[0] = -0.4;
  hitCoords[1] = 0;
  measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
  measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,20), TVector3(1,0,0), TVector3(0,1,0))), ++planeId);
  fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));



  //check
  assert(fitTrack.checkConsistency());

  // do the fit
  fitter->processTrack(&fitTrack);

  // print fit result
  fitTrack.getFittedState().Print();

  //check
  assert(fitTrack.checkConsistency());


  display->addEvent(&fitTrack);


  delete fitter;

  // open event display
  display->open();

}


