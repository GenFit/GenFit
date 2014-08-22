
#include "EventDisplay.h"

#include <assert.h>
#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <sys/time.h>

#include "AbsMeasurement.h"
#include "FullMeasurement.h"
#include "PlanarMeasurement.h"
#include "ProlateSpacepointMeasurement.h"
#include "SpacepointMeasurement.h"
#include "WireMeasurement.h"
#include "WirePointMeasurement.h"
#include "AbsTrackRep.h"
#include "ConstField.h"
#include "DetPlane.h"
#include "Exception.h"
#include "FieldManager.h"
#include "Tools.h"
#include "KalmanFitterInfo.h"
#include "KalmanFitter.h"
#include "DAF.h"
#include "KalmanFitterRefTrack.h"
#include "RKTrackRep.h"

#include <TApplication.h>
#include <TEveBrowser.h>
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveGeoNode.h>
#include <TEveGeoShape.h>
#include <TEveStraightLineSet.h>
#include <TEveTriangleSet.h>
#include <TDecompSVD.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TGeoEltu.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TGeoSphere.h>
#include <TGeoTube.h>
#include <TMath.h>
#include <TMatrixT.h>
#include <TMatrixTSym.h>
#include <TMatrixDSymEigen.h>
#include <TROOT.h>
#include <TVector2.h>
#include <TVectorD.h>
#include <TSystem.h>

#include "boost/scoped_ptr.hpp"


ClassImp(genfit::EventDisplay)

namespace genfit {


EventDisplay* EventDisplay::eventDisplay_ = NULL;

EventDisplay::EventDisplay() :
  errorScale_(1.),
  drawGeometry_(false),
  drawDetectors_(true),
  drawHits_(true),
  drawErrors_(true),
  drawPlanes_(true),
  drawTrackMarkers_(true),
  drawTrack_(true),
  drawRefTrack_(true),
  drawForward_(true),
  drawBackward_(true),
  drawAutoScale_(true),
  drawScaleMan_(false),
  drawSilent_(false),
  drawCardinalRep_(true),
  repId_(0),
  drawAllTracks_(true),
  trackId_(0),
  refit_(false),
  debugLvl_(0),
  fitterId_(SimpleKalman),
  mmHandling_(weightedAverage),
  squareRootFormalism_(false),
  dPVal_(1.E-3),
  dRelChi2_(0.2),
  dChi2Ref_(1.),
  nMinIter_(2),
  nMaxIter_(4),
  nMaxFailed_(-1),
  resort_(false)
{

  if((!gApplication) || (gApplication && gApplication->TestBit(TApplication::kDefaultApplication))) {
    std::cout << "In EventDisplay ctor: gApplication not found, creating..." << std::flush;
    new TApplication("ROOT_application", 0, 0);
    std::cout << "done!" << std::endl;
  }
  if(!gEve) {
    std::cout << "In EventDisplay ctor: gEve not found, creating..." << std::flush;
    TEveManager::Create();
    std::cout << "done!" << std::endl;
  }

  eventId_ = 0;

}

void EventDisplay::setOptions(std::string opts) {

  if(opts != "") {
    for(size_t i = 0; i < opts.length(); ++i) {
      if(opts[i] == 'A') drawAutoScale_ = true;
      if(opts[i] == 'B') drawBackward_ = true;
      if(opts[i] == 'D') drawDetectors_ = true;
      if(opts[i] == 'E') drawErrors_ = true;
      if(opts[i] == 'F') drawForward_ = true;
      if(opts[i] == 'H') drawHits_ = true;
      if(opts[i] == 'M') drawTrackMarkers_ = true;
      if(opts[i] == 'P') drawPlanes_ = true;
      if(opts[i] == 'S') drawScaleMan_ = true;
      if(opts[i] == 'T') drawTrack_ = true;
      if(opts[i] == 'X') drawSilent_ = true;
      if(opts[i] == 'G') drawGeometry_ = true;
    }
  }

}

void EventDisplay::setErrScale(double errScale) { errorScale_ = errScale; }

double EventDisplay::getErrScale() { return errorScale_; }

EventDisplay* EventDisplay::getInstance() {

  if(eventDisplay_ == NULL) {
    eventDisplay_ = new EventDisplay();
  }
  return eventDisplay_;

}

EventDisplay::~EventDisplay() { reset(); }

void EventDisplay::reset() {

  for(unsigned int i = 0; i < events_.size(); i++) {

    for(unsigned int j = 0; j < events_[i]->size(); j++) {

      delete events_[i]->at(j);

    }
    delete events_[i];
  }

  events_.clear();
}


void EventDisplay::addEvent(std::vector<Track*>& tracks) {

  std::vector<Track*>* vec = new std::vector<Track*>;

  for(unsigned int i = 0; i < tracks.size(); i++) {
    vec->push_back(new Track(*(tracks[i])));
  }

  events_.push_back(vec);
}


void EventDisplay::addEvent(std::vector<const Track*>& tracks) {

  std::vector<Track*>* vec = new std::vector<Track*>;

  for(unsigned int i = 0; i < tracks.size(); i++) {
    vec->push_back(new Track(*(tracks[i])));
  }

  events_.push_back(vec);
}


void EventDisplay::addEvent(const Track* tr) {

  std::vector<Track*>* vec = new std::vector<Track*>;
  vec->push_back(new Track(*tr));
  events_.push_back(vec);
}


void EventDisplay::next(unsigned int stp) {

  gotoEvent(eventId_ + stp);

}

void EventDisplay::prev(unsigned int stp) {

  if(events_.size() == 0) return;
  if(eventId_ < stp) {
    gotoEvent(0);
  } else {
    gotoEvent(eventId_ - stp);
  }

}

int EventDisplay::getNEvents() { return events_.size(); }


void EventDisplay::gotoEvent(unsigned int id) {

  if (events_.size() == 0)
    return;
  else if(id >= events_.size())
    id = events_.size() - 1;

  bool resetCam = true;

  if (id == eventId_)
    resetCam = false;

  eventId_ = id;

  std::cout << "At event " << id << std::endl;
  if (gEve->GetCurrentEvent()) {
    gEve->GetCurrentEvent()->DestroyElements();
  }
  double old_error_scale = errorScale_;
  drawEvent(eventId_, resetCam);
  if(old_error_scale != errorScale_) {
    if (gEve->GetCurrentEvent()) {
      gEve->GetCurrentEvent()->DestroyElements();
    }
    drawEvent(eventId_, resetCam); // if autoscaling changed the error, draw again.
  }
  errorScale_ = old_error_scale;

}

void EventDisplay::open() {

  std::cout << "EventDisplay::open(); " << getNEvents() << " events loaded" << std::endl;

  if(getNEvents() > 0) {
    double old_error_scale = errorScale_;
    drawEvent(0);
    if(old_error_scale != errorScale_) {
      std::cout << "autoscaling changed the error, draw again." << std::endl;
      gotoEvent(0); // if autoscaling changed the error, draw again.
    }
    errorScale_ = old_error_scale;
  }


  if(!drawSilent_) {
    makeGui();
    gApplication->Run(kTRUE);
  }

  std::cout << "opened" << std::endl;

}


void EventDisplay::drawEvent(unsigned int id, bool resetCam) {

  std::cout << "EventDisplay::drawEvent(" << id << ")" << std::endl;


  // draw the geometry, does not really work yet. If it's fixed, the docu in the header file should be changed.
  if(drawGeometry_) {
    TGeoNode* top_node = gGeoManager->GetTopNode();
    assert(top_node != NULL);

    //Set transparency & color of geometry
    TObjArray* volumes = gGeoManager->GetListOfVolumes();
    for(int i = 0; i < volumes->GetEntriesFast(); i++) {
      TGeoVolume* volume = dynamic_cast<TGeoVolume*>(volumes->At(i));
      assert(volume != NULL);
      volume->SetLineColor(12);
      volume->SetTransparency(50);
    }

    TEveGeoTopNode* eve_top_node = new TEveGeoTopNode(gGeoManager, top_node);
    eve_top_node->IncDenyDestroy();
    gEve->AddElement(eve_top_node);
  }


  for(unsigned int i = 0; i < events_.at(id)->size(); i++) { // loop over all tracks in an event

    if (!drawAllTracks_ && trackId_ != i)
      continue;

    Track* track = events_[id]->at(i);
    if (! track->checkConsistency()){
      std::cerr<<"track is not consistent"<<std::endl;
      continue;
    }


    boost::scoped_ptr<Track> refittedTrack(NULL);
    if (refit_) {

      std::cout << "Refit track:" << std::endl;

      boost::scoped_ptr<AbsKalmanFitter> fitter;
      switch (fitterId_) {
        case SimpleKalman:
          fitter.reset(new KalmanFitter(nMaxIter_, dPVal_));
          fitter->setMultipleMeasurementHandling(mmHandling_);
          (static_cast<KalmanFitter*>(fitter.get()))->useSquareRootFormalism(squareRootFormalism_);
          break;

        case RefKalman:
          fitter.reset(new KalmanFitterRefTrack(nMaxIter_, dPVal_));
          fitter->setMultipleMeasurementHandling(mmHandling_);
          static_cast<KalmanFitterRefTrack*>(fitter.get())->setDeltaChi2Ref(dChi2Ref_);
          break;

        case DafSimple:
          fitter.reset(new DAF(false));
          ( static_cast<KalmanFitter*>( (static_cast<DAF*>(fitter.get()))->getKalman() ) )->useSquareRootFormalism(squareRootFormalism_);
          break;
        case DafRef:
          fitter.reset(new DAF());
          ( static_cast<KalmanFitterRefTrack*>( (static_cast<DAF*>(fitter.get()))->getKalman() ) )->setDeltaChi2Ref(dChi2Ref_);
          break;

      }
      fitter->setDebugLvl(std::max(0, (int)debugLvl_-1));
      fitter->setMinIterations(nMinIter_);
      fitter->setMaxIterations(nMaxIter_);
      fitter->setRelChi2Change(dRelChi2_);
      fitter->setMaxFailedHits(nMaxFailed_);


      refittedTrack.reset(new Track(*track));
      refittedTrack->deleteFitterInfo();

      if (debugLvl_>0)
        refittedTrack->Print("C");

      timeval startcputime, endcputime;

      try{
        gettimeofday(&startcputime, NULL);
        fitter->processTrack(refittedTrack.get(), resort_);
        gettimeofday(&endcputime, NULL);
      }
      catch(genfit::Exception& e){
        std::cerr << e.what();
        std::cerr << "Exception, could not refit track" << std::endl;
        continue;
      }

      int microseconds = 1000000*(endcputime.tv_sec - startcputime.tv_sec) + (endcputime.tv_usec - startcputime.tv_usec);
      std::cout << "it took " << double(microseconds) /  1000 << " ms of CPU to fit the track\n";

      if (! refittedTrack->checkConsistency()){
        std::cerr<<"refittedTrack is not consistent"<<std::endl;
        continue;
      }

      track = refittedTrack.get();
    }




    AbsTrackRep* rep;

    if (drawCardinalRep_) {
      rep = track->getCardinalRep();
      std::cout << "Draw cardinal rep" << std::endl;
    }
    else {
      if (repId_ >= track->getNumReps())
        repId_ = track->getNumReps() - 1;
       rep = track->getTrackRep(repId_);
       std::cout << "Draw rep" << repId_ << std::endl;
    }

    if (debugLvl_>0) {
      //track->Print();
      track->Print("C");
      track->getFitStatus(rep)->Print();

      if (track->getFitStatus(rep)->isFitted()) {
        try {
          std::cout << "fitted state: \n";
          track->getFittedState().Print();
        }
        catch (Exception& e) {
          std::cerr << e.what();
        }
      }
    }



    rep->setPropDir(0);

    unsigned int numhits = track->getNumPointsWithMeasurement();

    KalmanFitterInfo* fi;
    KalmanFitterInfo* prevFi = 0;
    const MeasuredStateOnPlane* fittedState(NULL);
    const MeasuredStateOnPlane* prevFittedState(NULL);

    for(unsigned int j = 0; j < numhits; j++) { // loop over all hits in the track

      fittedState = NULL;

      TrackPoint* tp = track->getPointWithMeasurement(j);
      if (! tp->hasRawMeasurements()) {
        std::cerr<<"trackPoint has no raw measurements"<<std::endl;
        continue;
      }

      const AbsMeasurement* m = tp->getRawMeasurement();
      int hit_coords_dim = m->getDim();

      // check if multiple AbsMeasurements are of same type
      if (tp->getNumRawMeasurements() > 1) {
        bool sameTypes(true);
        for (unsigned int iM=1; iM<tp->getNumRawMeasurements(); ++iM) {
          if (typeid(*(tp->getRawMeasurement(iM))) != typeid(*m))
            sameTypes = false;
        }
        if (!sameTypes) {
          std::cerr<<"cannot draw trackpoint containing multiple Measurements of differend types"<<std::endl;
          continue;
        }
      }



      // get the fitter infos ------------------------------------------------------------------
      if (! tp->hasFitterInfo(rep)) {
        std::cerr<<"trackPoint has no fitterInfo for rep"<<std::endl;
        continue;
      }

      AbsFitterInfo* fitterInfo = tp->getFitterInfo(rep);

      fi = dynamic_cast<KalmanFitterInfo*>(fitterInfo);
      if(fi == NULL) {
        std::cerr<<"can only display KalmanFitterInfo"<<std::endl;
        continue;
      }
      if (! fi->hasPredictionsAndUpdates()) {
        std::cerr<<"KalmanFitterInfo does not have all predictions and updates"<<std::endl;
        //continue;
      }
      else {
        try {
          fittedState = &(fi->getFittedState(true));
        }
        catch (Exception& e) {
          std::cerr << e.what();
          std::cerr<<"can not get fitted state"<<std::endl;
          fittedState = NULL;
          prevFi = fi;
          prevFittedState = fittedState;
          continue;
        }
      }

      if (fittedState == NULL) {
        if (fi->hasForwardUpdate()) {
          fittedState = fi->getForwardUpdate();
        }
        else if (fi->hasBackwardUpdate()) {
          fittedState = fi->getBackwardUpdate();
        }
        else if (fi->hasForwardPrediction()) {
          fittedState = fi->getForwardPrediction();
        }
        else if (fi->hasBackwardPrediction()) {
          fittedState = fi->getBackwardPrediction();
        }
      }

      if (fittedState == NULL) {
        std::cout << "canot get any state from fitterInfo, continue.\n";
        prevFi = fi;
        prevFittedState = fittedState;
        continue;
      }

      TVector3 track_pos = fittedState->getPos();
      double charge = fittedState->getCharge();

      //std::cout << "trackPos: "; track_pos.Print();


      // determine measurement type
      bool full_hit = false;
      bool planar_hit = false;
      bool planar_pixel_hit = false;
      bool space_hit = false;
      bool wire_hit = false;
      bool wirepoint_hit = false;
      if (dynamic_cast<const FullMeasurement*>(m) != NULL) {
        full_hit = true;
      }
      else if(dynamic_cast<const PlanarMeasurement*>(m) != NULL) {
        planar_hit = true;
        if(hit_coords_dim == 2) {
          planar_pixel_hit = true;
        }
      } else if (dynamic_cast<const SpacepointMeasurement*>(m) != NULL) {
        space_hit = true;
      } else if (dynamic_cast<const WireMeasurement*>(m) != NULL) {
        wire_hit = true;
        if (dynamic_cast<const WirePointMeasurement*>(m) != NULL) {
          wirepoint_hit = true;
        }
      } else {
        std::cout << "Track " << i << ", Hit " << j << ": Unknown measurement type: skipping hit!" << std::endl;
        continue;
      }


      // loop over MeasurementOnPlanes
      unsigned int nMeas = fi->getNumMeasurements();
      for (unsigned int iMeas = 0; iMeas < nMeas; ++iMeas) {

        if (iMeas > 0 && wire_hit)
          break;

        const MeasurementOnPlane* mop = fi->getMeasurementOnPlane(iMeas);
        const TVectorT<double>& hit_coords = mop->getState();
        const TMatrixTSym<double>& hit_cov = mop->getCov();

        // finished getting the hit infos -----------------------------------------------------

        // sort hit infos into variables ------------------------------------------------------
        TVector3 o = fittedState->getPlane()->getO();
        TVector3 u = fittedState->getPlane()->getU();
        TVector3 v = fittedState->getPlane()->getV();

        double_t hit_u = 0;
        double_t hit_v = 0;
        double_t plane_size = 4;
        TVector2 stripDir(1,0);

        if(planar_hit) {
          if(!planar_pixel_hit) {
            if (dynamic_cast<RKTrackRep*>(rep) != NULL) {
              const TMatrixD& H = mop->getHMatrix()->getMatrix();
              stripDir.Set(H(0,3), H(0,4));
            }
            hit_u = hit_coords(0);
          } else {
            hit_u = hit_coords(0);
            hit_v = hit_coords(1);
          }
        } else if (wire_hit) {
          hit_u = fabs(hit_coords(0));
          hit_v = v*(track_pos-o); // move the covariance tube so that the track goes through it
          if (wirepoint_hit) {
            hit_v = hit_coords(1);
          }
        }

        if(plane_size < 4) plane_size = 4;
        // finished setting variables ---------------------------------------------------------

        // draw planes if corresponding option is set -----------------------------------------
        if(iMeas == 0 &&
           (drawPlanes_ || (drawDetectors_ && planar_hit))) {
          TVector3 move(0,0,0);
          if (planar_hit) move = track_pos-o;
          if (wire_hit) move = v*(v*(track_pos-o)); // move the plane along the wire until the track goes through it
          TEveBox* box = boxCreator(o + move, u, v, plane_size, plane_size, 0.01);
          if (drawDetectors_ && planar_hit) {
            box->SetMainColor(kCyan);
          } else {
            box->SetMainColor(kGray);
          }
          box->SetMainTransparency(50);
          gEve->AddElement(box);
        }
        // finished drawing planes ------------------------------------------------------------

        // draw track if corresponding option is set ------------------------------------------
        if (j == 0) {
          if (drawBackward_) {
              MeasuredStateOnPlane update ( *fi->getBackwardUpdate() );
              update.extrapolateBy(-3.);
              makeLines(&update, fi->getBackwardUpdate(), rep, kMagenta, 1, drawTrackMarkers_, drawErrors_, 1);
          }
        }
        if (j > 0 && prevFi != NULL) {
          if(drawTrack_) {
            makeLines(prevFittedState, fittedState, rep, charge > 0 ? kRed : kBlue, 1, drawTrackMarkers_, drawErrors_, 3);
            if (drawErrors_) { // make sure to draw errors in both directions
              makeLines(prevFittedState, fittedState, rep, charge > 0 ? kRed : kBlue, 1, false, drawErrors_, 0, 0);
            }
          }
          if (drawForward_) {
            makeLines(prevFi->getForwardUpdate(), fi->getForwardPrediction(), rep, kCyan, 1, drawTrackMarkers_, drawErrors_, 1, 0);
            if (j == numhits-1) {
              MeasuredStateOnPlane update ( *fi->getForwardUpdate() );
              update.extrapolateBy(3.);
              makeLines(fi->getForwardUpdate(), &update, rep, kCyan, 1, drawTrackMarkers_, drawErrors_, 1, 0);
            }
          }
          if (drawBackward_) {
            makeLines(prevFi->getBackwardPrediction(), fi->getBackwardUpdate(), rep, kMagenta, 1, drawTrackMarkers_, drawErrors_, 1);
          }
          // draw reference track if corresponding option is set ------------------------------------------
          if(drawRefTrack_ && fi->hasReferenceState() && prevFi->hasReferenceState())
            makeLines(prevFi->getReferenceState(), fi->getReferenceState(), rep, charge > 0 ? kRed + 2 : kBlue + 2, 2, drawTrackMarkers_, false, 3);
        }
        else if (j > 0 && prevFi == NULL) {
          std::cout << "previous FitterInfo == NULL \n";
        }

        // draw detectors if option is set, only important for wire hits ----------------------
        if(drawDetectors_) {

          if(wire_hit) {
            TEveGeoShape* det_shape = new TEveGeoShape("det_shape");
            det_shape->IncDenyDestroy();
            det_shape->SetShape(new TGeoTube(std::max(0., (double)(hit_u-0.0105/2.)), hit_u+0.0105/2., plane_size));

            TVector3 norm = u.Cross(v);
            TGeoRotation det_rot("det_rot", (u.Theta()*180)/TMath::Pi(), (u.Phi()*180)/TMath::Pi(),
                (norm.Theta()*180)/TMath::Pi(), (norm.Phi()*180)/TMath::Pi(),
                (v.Theta()*180)/TMath::Pi(), (v.Phi()*180)/TMath::Pi()); // move the tube to the right place and rotate it correctly
            TVector3 move = v*(v*(track_pos-o)); // move the tube along the wire until the track goes through it
            TGeoCombiTrans det_trans(o(0) + move.X(),
                                     o(1) + move.Y(),
                                     o(2) + move.Z(),
                                     &det_rot);
            det_shape->SetTransMatrix(det_trans);
            det_shape->SetMainColor(kCyan);
            det_shape->SetMainTransparency(25);
            if((drawHits_ && (hit_u+0.0105/2 > 0)) || !drawHits_) {
              gEve->AddElement(det_shape);
            }
          }

        }
        // finished drawing detectors ---------------------------------------------------------

        if(drawHits_) {

          // draw planar hits, with distinction between strip and pixel hits ----------------
          if (full_hit) {

            StateOnPlane dummy(rep);
            StateOnPlane dummy2(TVectorD(rep->getDim()), static_cast<const FullMeasurement*>(m)->constructPlane(dummy), rep);
            MeasuredStateOnPlane sop = *(static_cast<const FullMeasurement*>(m)->constructMeasurementsOnPlane(dummy2)[0]);
            sop.getCov()*=errorScale_;

            MeasuredStateOnPlane prevSop(sop);
            prevSop.extrapolateBy(-3);
            makeLines(&sop, &prevSop, rep, kYellow, 1, false, true, 0, 0);

            prevSop = sop;
            prevSop.extrapolateBy(3);
            makeLines(&sop, &prevSop, rep, kYellow, 1, false, true, 0, 0);
          }

          if(planar_hit) {
            if(!planar_pixel_hit) {
              TEveBox* hit_box;
              TVector3 stripDir3 = stripDir.X()*u + stripDir.Y()*v;
              TVector3 stripDir3perp = stripDir.Y()*u - stripDir.X()*v;
              TVector3 move = stripDir3perp*(stripDir3perp*(track_pos-o));
              hit_box = boxCreator((o + move + hit_u*stripDir3), stripDir3, stripDir3perp, errorScale_*std::sqrt(hit_cov(0,0)), plane_size, 0.0105);
              hit_box->SetMainColor(kYellow);
              hit_box->SetMainTransparency(0);
              gEve->AddElement(hit_box);
            } else {
              // calculate eigenvalues to draw error-ellipse ----------------------------
              TMatrixDSymEigen eigen_values(hit_cov);
              TEveGeoShape* cov_shape = new TEveGeoShape("cov_shape");
              cov_shape->IncDenyDestroy();
              TVectorT<double> ev = eigen_values.GetEigenValues();
              TMatrixT<double> eVec = eigen_values.GetEigenVectors();
              double pseudo_res_0 = errorScale_*std::sqrt(ev(0));
              double pseudo_res_1 = errorScale_*std::sqrt(ev(1));
              // finished calcluating, got the values -----------------------------------

              // do autoscaling if necessary --------------------------------------------
              if(drawAutoScale_) {
                double min_cov = std::min(pseudo_res_0,pseudo_res_1);
                if(min_cov < 1e-5) {
                  std::cout << "Track " << i << ", Hit " << j << ": Invalid covariance matrix (Eigenvalue < 1e-5), autoscaling not possible!" << std::endl;
                } else {
                  if(min_cov < 0.049) {
                    double cor = 0.05 / min_cov;
                    std::cout << "Track " << i << ", Hit " << j << ": Pixel covariance too small, rescaling by " << cor;
                    errorScale_ *= cor;
                    pseudo_res_0 *= cor;
                    pseudo_res_1 *= cor;
                    std::cout << " to " << errorScale_ << std::endl;
                  }
                }
              }
              // finished autoscaling ---------------------------------------------------

              // calculate the semiaxis of the error ellipse ----------------------------
              cov_shape->SetShape(new TGeoEltu(pseudo_res_0, pseudo_res_1, 0.0105));
              TVector3 pix_pos = o + hit_u*u + hit_v*v;
              TVector3 u_semiaxis = (pix_pos + eVec(0,0)*u + eVec(1,0)*v)-pix_pos;
              TVector3 v_semiaxis = (pix_pos + eVec(0,1)*u + eVec(1,1)*v)-pix_pos;
              TVector3 norm = u.Cross(v);
              // finished calculating ---------------------------------------------------

              // rotate and translate everything correctly ------------------------------
              TGeoRotation det_rot("det_rot", (u_semiaxis.Theta()*180)/TMath::Pi(), (u_semiaxis.Phi()*180)/TMath::Pi(),
                  (v_semiaxis.Theta()*180)/TMath::Pi(), (v_semiaxis.Phi()*180)/TMath::Pi(),
                  (norm.Theta()*180)/TMath::Pi(), (norm.Phi()*180)/TMath::Pi());
              TGeoCombiTrans det_trans(pix_pos(0),pix_pos(1),pix_pos(2), &det_rot);
              cov_shape->SetTransMatrix(det_trans);
              // finished rotating and translating --------------------------------------

              cov_shape->SetMainColor(kYellow);
              cov_shape->SetMainTransparency(0);
              gEve->AddElement(cov_shape);
            }
          }
          // finished drawing planar hits ---------------------------------------------------

          // draw spacepoint hits -----------------------------------------------------------
          if(space_hit) {
            {
              // get eigenvalues of covariance to know how to draw the ellipsoid ------------
              TMatrixDSymEigen eigen_values(m->getRawHitCov());
              TEveGeoShape* cov_shape = new TEveGeoShape("cov_shape");
              cov_shape->IncDenyDestroy();
              cov_shape->SetShape(new TGeoSphere(0.,1.));
              TVectorT<double> ev = eigen_values.GetEigenValues();
              TMatrixT<double> eVec = eigen_values.GetEigenVectors();
              TVector3 eVec1(eVec(0,0),eVec(1,0),eVec(2,0));
              TVector3 eVec2(eVec(0,1),eVec(1,1),eVec(2,1));
              TVector3 eVec3(eVec(0,2),eVec(1,2),eVec(2,2));
              const TVector3 norm = u.Cross(v);
              // got everything we need -----------------------------------------------------

              static const double radDeg(180./TMath::Pi());
              TGeoRotation det_rot("det_rot", eVec1.Theta()*radDeg, eVec1.Phi()*radDeg,
                  eVec2.Theta()*radDeg, eVec2.Phi()*radDeg,
                  eVec3.Theta()*radDeg, eVec3.Phi()*radDeg);

              if (! det_rot.IsValid()){
                // hackish fix if eigenvectors are not orthonogonal
                if (fabs(eVec2*eVec3) > 1.e-10)
                  eVec3 = eVec1.Cross(eVec2);

                det_rot.SetAngles(eVec1.Theta()*radDeg, eVec1.Phi()*radDeg,
                    eVec2.Theta()*radDeg, eVec2.Phi()*radDeg,
                    eVec3.Theta()*radDeg, eVec3.Phi()*radDeg);
              }

              // set the scaled eigenvalues -------------------------------------------------
              double pseudo_res_0 = errorScale_*std::sqrt(ev(0));
              double pseudo_res_1 = errorScale_*std::sqrt(ev(1));
              double pseudo_res_2 = errorScale_*std::sqrt(ev(2));
              if(drawScaleMan_) { // override again if necessary
                pseudo_res_0 = errorScale_*0.5;
                pseudo_res_1 = errorScale_*0.5;
                pseudo_res_2 = errorScale_*0.5;
              }
              // finished scaling -----------------------------------------------------------

              // autoscale if necessary -----------------------------------------------------
              if(drawAutoScale_) {
                double min_cov = std::min(pseudo_res_0,std::min(pseudo_res_1,pseudo_res_2));
                if(min_cov < 1e-5) {
                  std::cout << "Track " << i << ", Hit " << j << ": Invalid covariance matrix (Eigenvalue < 1e-5), autoscaling not possible!" << std::endl;
                } else {
                  if(min_cov <= 0.149) {
                    double cor = 0.15 / min_cov;
                    std::cout << "Track " << i << ", Hit " << j << ": Space hit covariance too small, rescaling by " << cor;
                    errorScale_ *= cor;
                    pseudo_res_0 *= cor;
                    pseudo_res_1 *= cor;
                    pseudo_res_2 *= cor;
                    std::cout << " to " << errorScale_ << std::endl;

                  }
                }
              }
              // finished autoscaling -------------------------------------------------------

              // rotate and translate -------------------------------------------------------
              TGeoGenTrans det_trans(o(0),o(1),o(2),
                                     //std::sqrt(pseudo_res_0/pseudo_res_1/pseudo_res_2), std::sqrt(pseudo_res_1/pseudo_res_0/pseudo_res_2), std::sqrt(pseudo_res_2/pseudo_res_0/pseudo_res_1), // this workaround is necessary due to the "normalization" performed in  TGeoGenTrans::SetScale
                                     //1/(pseudo_res_0),1/(pseudo_res_1),1/(pseudo_res_2),
                                     pseudo_res_0, pseudo_res_1, pseudo_res_2,
                                     &det_rot);
              cov_shape->SetTransMatrix(det_trans);
              // finished rotating and translating ------------------------------------------

              cov_shape->SetMainColor(kYellow);
              cov_shape->SetMainTransparency(10);
              gEve->AddElement(cov_shape);
            }


            {
              // calculate eigenvalues to draw error-ellipse ----------------------------
              TMatrixDSymEigen eigen_values(hit_cov);
              TEveGeoShape* cov_shape = new TEveGeoShape("cov_shape");
              cov_shape->IncDenyDestroy();
              TVectorT<double> ev = eigen_values.GetEigenValues();
              TMatrixT<double> eVec = eigen_values.GetEigenVectors();
              double pseudo_res_0 = errorScale_*std::sqrt(ev(0));
              double pseudo_res_1 = errorScale_*std::sqrt(ev(1));
              // finished calcluating, got the values -----------------------------------

              // do autoscaling if necessary --------------------------------------------
              if(drawAutoScale_) {
                double min_cov = std::min(pseudo_res_0,pseudo_res_1);
                if(min_cov < 1e-5) {
                  std::cout << "Track " << i << ", Hit " << j << ": Invalid covariance matrix (Eigenvalue < 1e-5), autoscaling not possible!" << std::endl;
                } else {
                  if(min_cov < 0.049) {
                    double cor = 0.05 / min_cov;
                    std::cout << "Track " << i << ", Hit " << j << ": Pixel covariance too small, rescaling by " << cor;
                    errorScale_ *= cor;
                    pseudo_res_0 *= cor;
                    pseudo_res_1 *= cor;
                    std::cout << " to " << errorScale_ << std::endl;
                  }
                }
              }
              // finished autoscaling ---------------------------------------------------

              // calculate the semiaxis of the error ellipse ----------------------------
              cov_shape->SetShape(new TGeoEltu(pseudo_res_0, pseudo_res_1, 0.0105));
              TVector3 pix_pos = o + hit_u*u + hit_v*v;
              TVector3 u_semiaxis = (pix_pos + eVec(0,0)*u + eVec(1,0)*v)-pix_pos;
              TVector3 v_semiaxis = (pix_pos + eVec(0,1)*u + eVec(1,1)*v)-pix_pos;
              TVector3 norm = u.Cross(v);
              // finished calculating ---------------------------------------------------

              // rotate and translate everything correctly ------------------------------
              static const double radDeg(180./TMath::Pi());
              TGeoRotation det_rot("det_rot", u_semiaxis.Theta()*radDeg, u_semiaxis.Phi()*radDeg,
                  v_semiaxis.Theta()*radDeg, v_semiaxis.Phi()*radDeg,
                  norm.Theta()*radDeg, norm.Phi()*radDeg);
              /*if (! det_rot.IsValid()){
                u_semiaxis.Print();
                v_semiaxis.Print();
                norm.Print();
              }*/
              TGeoCombiTrans det_trans(pix_pos(0),pix_pos(1),pix_pos(2), &det_rot);
              cov_shape->SetTransMatrix(det_trans);
              // finished rotating and translating --------------------------------------

              cov_shape->SetMainColor(kYellow);
              cov_shape->SetMainTransparency(0);
              gEve->AddElement(cov_shape);
            }
          }
          // finished drawing spacepoint hits -----------------------------------------------

          // draw wire hits -----------------------------------------------------------------
          if(wire_hit) {
            TEveGeoShape* cov_shape = new TEveGeoShape("cov_shape");
            cov_shape->IncDenyDestroy();
            double pseudo_res_0 = errorScale_*std::sqrt(hit_cov(0,0));
            double pseudo_res_1 = plane_size;
            if (wirepoint_hit) pseudo_res_1 = errorScale_*std::sqrt(hit_cov(1,1));

            // autoscale if necessary -----------------------------------------------------
            if(drawAutoScale_) {
              if(pseudo_res_0 < 1e-5) {
                std::cout << "Track " << i << ", Hit " << j << ": Invalid wire resolution (< 1e-5), autoscaling not possible!" << std::endl;
              } else {
                if(pseudo_res_0 < 0.0049) {
                  double cor = 0.005 / pseudo_res_0;
                  std::cout << "Track " << i << ", Hit " << j << ": Wire covariance too small, rescaling by " << cor;
                  errorScale_ *= cor;
                  pseudo_res_0 *= cor;
                  std::cout << " to " << errorScale_ << std::endl;
                }
              }

              if(wirepoint_hit && pseudo_res_1 < 1e-5) {
                std::cout << "Track " << i << ", Hit " << j << ": Invalid wire resolution (< 1e-5), autoscaling not possible!" << std::endl;
              } else {
                if(pseudo_res_1 < 0.0049) {
                  double cor = 0.005 / pseudo_res_1;
                  std::cout << "Track " << i << ", Hit " << j << ": Wire covariance too small, rescaling by " << cor;
                  errorScale_ *= cor;
                  pseudo_res_1 *= cor;
                  std::cout << " to " << errorScale_ << std::endl;
                }
              }
            }
            // finished autoscaling -------------------------------------------------------

            TEveBox* hit_box;
            TVector3 move = v*(v*(track_pos-o));
            hit_box = boxCreator((o + move + hit_u*u), u, v, errorScale_*std::sqrt(hit_cov(0,0)), pseudo_res_1, 0.0105);
            hit_box->SetMainColor(kYellow);
            hit_box->SetMainTransparency(0);
            gEve->AddElement(hit_box);

            hit_box = boxCreator((o + move - hit_u*u), u, v, errorScale_*std::sqrt(hit_cov(0,0)), pseudo_res_1, 0.0105);
            hit_box->SetMainColor(kYellow);
            hit_box->SetMainTransparency(0);
            gEve->AddElement(hit_box);
          }
          // finished drawing wire hits -----------------------------------------------------

        } // finished drawing hits

      } // finished looping over MeasurmentOnPlanes


      prevFi = fi;
      prevFittedState = fittedState;

    }

  }

  gEve->Redraw3D(resetCam);

}




TEveBox* EventDisplay::boxCreator(TVector3 o, TVector3 u, TVector3 v, float ud, float vd, float depth) {

  TEveBox* box = new TEveBox("detPlane_shape");
  float vertices[24];

  TVector3 norm = u.Cross(v);
  u *= (0.5*ud);
  v *= (0.5*vd);
  norm *= (0.5*depth);

  vertices[0] = o(0) - u(0) - v(0) - norm(0);
  vertices[1] = o(1) - u(1) - v(1) - norm(1);
  vertices[2] = o(2) - u(2) - v(2) - norm(2);

  vertices[3] = o(0) + u(0) - v(0) - norm(0);
  vertices[4] = o(1) + u(1) - v(1) - norm(1);
  vertices[5] = o(2) + u(2) - v(2) - norm(2);

  vertices[6] = o(0) + u(0) - v(0) + norm(0);
  vertices[7] = o(1) + u(1) - v(1) + norm(1);
  vertices[8] = o(2) + u(2) - v(2) + norm(2);

  vertices[9] = o(0) - u(0) - v(0) + norm(0);
  vertices[10] = o(1) - u(1) - v(1) + norm(1);
  vertices[11] = o(2) - u(2) - v(2) + norm(2);

  vertices[12] = o(0) - u(0) + v(0) - norm(0);
  vertices[13] = o(1) - u(1) + v(1) - norm(1);
  vertices[14] = o(2) - u(2) + v(2) - norm(2);

  vertices[15] = o(0) + u(0) + v(0) - norm(0);
  vertices[16] = o(1) + u(1) + v(1) - norm(1);
  vertices[17] = o(2) + u(2) + v(2) - norm(2);

  vertices[18] = o(0) + u(0) + v(0) + norm(0);
  vertices[19] = o(1) + u(1) + v(1) + norm(1);
  vertices[20] = o(2) + u(2) + v(2) + norm(2);

  vertices[21] = o(0) - u(0) + v(0) + norm(0);
  vertices[22] = o(1) - u(1) + v(1) + norm(1);
  vertices[23] = o(2) - u(2) + v(2) + norm(2);


  for(int k = 0; k < 24; k += 3) box->SetVertex((k/3), vertices[k], vertices[k+1], vertices[k+2]);

  return box;

}


void EventDisplay::makeLines(const StateOnPlane* prevState, const StateOnPlane* state, const AbsTrackRep* rep,
    const Color_t& color, const Style_t& style, bool drawMarkers, bool drawErrors, double lineWidth, int markerPos)
{
  if (prevState == NULL || state == NULL) {
    std::cerr << "prevState == NULL || state == NULL\n";
    return;
  }

  TVector3 pos, dir, oldPos, oldDir;
  rep->getPosDir(*state, pos, dir);
  rep->getPosDir(*prevState, oldPos, oldDir);

  double distA = (pos-oldPos).Mag();
  double distB = distA;
  if ((pos-oldPos)*oldDir < 0)
    distA *= -1.;
  if ((pos-oldPos)*dir < 0)
    distB *= -1.;
  TVector3 intermediate1 = oldPos + 0.3 * distA * oldDir;
  TVector3 intermediate2 = pos - 0.3 * distB * dir;
  TEveStraightLineSet* lineSet = new TEveStraightLineSet;
  lineSet->AddLine(oldPos(0), oldPos(1), oldPos(2), intermediate1(0), intermediate1(1), intermediate1(2));
  lineSet->AddLine(intermediate1(0), intermediate1(1), intermediate1(2), intermediate2(0), intermediate2(1), intermediate2(2));
  lineSet->AddLine(intermediate2(0), intermediate2(1), intermediate2(2), pos(0), pos(1), pos(2));
  lineSet->SetLineColor(color);
  lineSet->SetLineStyle(style);
  lineSet->SetLineWidth(lineWidth);
  if (drawMarkers) {
    if (markerPos == 0)
      lineSet->AddMarker(oldPos(0), oldPos(1), oldPos(2));
    else
      lineSet->AddMarker(pos(0), pos(1), pos(2));
  }

  if (lineWidth > 0)
    gEve->AddElement(lineSet);


  if (drawErrors) {
    const MeasuredStateOnPlane* measuredState;
    if (markerPos == 0)
      measuredState = dynamic_cast<const MeasuredStateOnPlane*>(prevState);
    else
      measuredState = dynamic_cast<const MeasuredStateOnPlane*>(state);

    if (measuredState != NULL) {

      // step for evaluate at a distance from the original plane
      TVector3 eval;
      if (markerPos == 0) {
        if (fabs(distA) < 1.) {
          distA < 0 ? distA = -1 : distA = 1;
        }
        eval = 0.2 * distA * oldDir;
      }
      else {
        if (fabs(distB) < 1.) {
          distB < 0 ? distB = -1 : distB = 1;
        }
        eval = -0.2 * distB * dir;
      }


      // get cov at first plane
      TMatrixDSym cov;
      TVector3 position, direction;
      rep->getPosMomCov(*measuredState, position, direction, cov);

      // get eigenvalues & -vectors
      TMatrixDSymEigen eigen_values(cov.GetSub(0,2, 0,2));
      TVectorT<double> ev = eigen_values.GetEigenValues();
      TMatrixT<double> eVec = eigen_values.GetEigenVectors();
      TVector3 eVec1, eVec2;
      // limit
      static const double maxErr = 1000.;
      double ev0 = std::min(ev(0), maxErr);
      double ev1 = std::min(ev(1), maxErr);
      double ev2 = std::min(ev(2), maxErr);

      // get two largest eigenvalues/-vectors
      if (ev0 < ev1 && ev0 < ev2) {
        eVec1.SetXYZ(eVec(0,1),eVec(1,1),eVec(2,1));
        eVec1 *= sqrt(ev1);
        eVec2.SetXYZ(eVec(0,2),eVec(1,2),eVec(2,2));
        eVec2 *= sqrt(ev2);
      }
      else if (ev1 < ev0 && ev1 < ev2) {
        eVec1.SetXYZ(eVec(0,0),eVec(1,0),eVec(2,0));
        eVec1 *= sqrt(ev0);
        eVec2.SetXYZ(eVec(0,2),eVec(1,2),eVec(2,2));
        eVec2 *= sqrt(ev2);
      }
      else {
        eVec1.SetXYZ(eVec(0,0),eVec(1,0),eVec(2,0));
        eVec1 *= sqrt(ev0);
        eVec2.SetXYZ(eVec(0,1),eVec(1,1),eVec(2,1));
        eVec2 *= sqrt(ev1);
      }

      if (eVec1.Cross(eVec2)*eval < 0)
        eVec2 *= -1;
      //assert(eVec1.Cross(eVec2)*eval > 0);

      const TVector3 oldEVec1(eVec1);
      const TVector3 oldEVec2(eVec2);

      const int nEdges = 24;
      std::vector<TVector3> vertices;

      vertices.push_back(position);

      // vertices at plane
      for (int i=0; i<nEdges; ++i) {
        const double angle = 2*TMath::Pi()/nEdges * i;
        vertices.push_back(position + cos(angle)*eVec1 + sin(angle)*eVec2);
      }



      DetPlane* newPlane = new DetPlane(*(measuredState->getPlane()));
      newPlane->setO(position + eval);

      MeasuredStateOnPlane stateCopy(*measuredState);
      try{
        rep->extrapolateToPlane(stateCopy, SharedPlanePtr(newPlane));
      }
      catch(Exception& e){
        std::cerr<<e.what();
        return;
      }

      // get cov at 2nd plane
      rep->getPosMomCov(stateCopy, position, direction, cov);

      // get eigenvalues & -vectors
      TMatrixDSymEigen eigen_values2(cov.GetSub(0,2, 0,2));
      ev = eigen_values2.GetEigenValues();
      eVec = eigen_values2.GetEigenVectors();
      // limit
      ev0 = std::min(ev(0), maxErr);
      ev1 = std::min(ev(1), maxErr);
      ev2 = std::min(ev(2), maxErr);

      // get two largest eigenvalues/-vectors
      if (ev0 < ev1 && ev0 < ev2) {
        eVec1.SetXYZ(eVec(0,1),eVec(1,1),eVec(2,1));
        eVec1 *= sqrt(ev1);
        eVec2.SetXYZ(eVec(0,2),eVec(1,2),eVec(2,2));
        eVec2 *= sqrt(ev2);
      }
      else if (ev1 < ev0 && ev1 < ev2) {
        eVec1.SetXYZ(eVec(0,0),eVec(1,0),eVec(2,0));
        eVec1 *= sqrt(ev0);
        eVec2.SetXYZ(eVec(0,2),eVec(1,2),eVec(2,2));
        eVec2 *= sqrt(ev2);
      }
      else {
        eVec1.SetXYZ(eVec(0,0),eVec(1,0),eVec(2,0));
        eVec1 *= sqrt(ev0);
        eVec2.SetXYZ(eVec(0,1),eVec(1,1),eVec(2,1));
        eVec2 *= sqrt(ev1);
      }

      if (eVec1.Cross(eVec2)*eval < 0)
        eVec2 *= -1;
      //assert(eVec1.Cross(eVec2)*eval > 0);

      if (oldEVec1*eVec1 < 0) {
        eVec1 *= -1;
        eVec2 *= -1;
      }

      // vertices at 2nd plane
      double angle0 = eVec1.Angle(oldEVec1);
      if (eVec1*(eval.Cross(oldEVec1)) < 0)
        angle0 *= -1;
      for (int i=0; i<nEdges; ++i) {
        const double angle = 2*TMath::Pi()/nEdges * i - angle0;
        vertices.push_back(position + cos(angle)*eVec1 + sin(angle)*eVec2);
      }

      vertices.push_back(position);


      TEveTriangleSet* error_shape = new TEveTriangleSet(vertices.size(), nEdges*2);
      for(unsigned int k = 0; k < vertices.size(); ++k) {
        error_shape->SetVertex(k, vertices[k].X(), vertices[k].Y(), vertices[k].Z());
      }

      assert(vertices.size() == 2*nEdges+2);

      int iTri(0);
      for (int i=0; i<nEdges; ++i) {
        //error_shape->SetTriangle(iTri++,  0,             i+1,        (i+1)%nEdges+1);
        error_shape->SetTriangle(iTri++,  i+1,           i+1+nEdges, (i+1)%nEdges+1);
        error_shape->SetTriangle(iTri++, (i+1)%nEdges+1, i+1+nEdges, (i+1)%nEdges+1+nEdges);
        //error_shape->SetTriangle(iTri++,  2*nEdges+1,    i+1+nEdges, (i+1)%nEdges+1+nEdges);
      }

      //assert(iTri == nEdges*4);

      error_shape->SetMainColor(color);
      error_shape->SetMainTransparency(25);
      gEve->AddElement(error_shape);
    }
  }
}


void EventDisplay::makeGui() {

  TEveBrowser* browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft);

  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain->SetWindowName("XX GUI");
  frmMain->SetCleanup(kDeepCleanup);

  TGLabel* lbl = 0;
  TGTextButton* tb = 0;
  EventDisplay*  fh = EventDisplay::getInstance();

  TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain); {
    // evt number entry
    lbl = new TGLabel(hf, "Go to event: ");
    hf->AddFrame(lbl);
    guiEvent = new TGNumberEntry(hf, 0, 9,999, TGNumberFormat::kNESInteger,
                          TGNumberFormat::kNEANonNegative,
                          TGNumberFormat::kNELLimitMinMax,
                          0, 99999);
    hf->AddFrame(guiEvent);
    guiEvent->Connect("ValueSet(Long_t)", "genfit::EventDisplay", fh, "guiGoto()");

    // redraw button
    tb = new TGTextButton(hf, "Redraw Event");
    hf->AddFrame(tb);
    tb->Connect("Clicked()", "genfit::EventDisplay", fh, "guiGoto()");
  }
  frmMain->AddFrame(hf);

  // draw options
  hf = new TGHorizontalFrame(frmMain); {
    lbl = new TGLabel(hf, "\n Draw Options");
    hf->AddFrame(lbl);
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiDrawGeometry_ =  new TGCheckButton(hf, "Draw geometry");
    if(drawGeometry_) guiDrawGeometry_->Toggle();
    hf->AddFrame(guiDrawGeometry_);
    guiDrawGeometry_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiDrawDetectors_ =  new TGCheckButton(hf, "Draw detectors");
    if(drawDetectors_) guiDrawDetectors_->Toggle();
    hf->AddFrame(guiDrawDetectors_);
    guiDrawDetectors_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiDrawHits_ =  new TGCheckButton(hf, "Draw hits");
    if(drawHits_) guiDrawHits_->Toggle();
    hf->AddFrame(guiDrawHits_);
    guiDrawHits_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);



  hf = new TGHorizontalFrame(frmMain); {
    guiDrawPlanes_ =  new TGCheckButton(hf, "Draw planes");
    if(drawPlanes_) guiDrawPlanes_->Toggle();
    hf->AddFrame(guiDrawPlanes_);
    guiDrawPlanes_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiDrawTrackMarkers_ =  new TGCheckButton(hf, "Draw track markers");
    if(drawTrackMarkers_) guiDrawTrackMarkers_->Toggle();
    hf->AddFrame(guiDrawTrackMarkers_);
    guiDrawTrackMarkers_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);


  hf = new TGHorizontalFrame(frmMain); {
    guiDrawTrack_ =  new TGCheckButton(hf, "Draw track");
    if(drawTrack_) guiDrawTrack_->Toggle();
    hf->AddFrame(guiDrawTrack_);
    guiDrawTrack_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiDrawRefTrack_ =  new TGCheckButton(hf, "Draw reference track");
    if(drawRefTrack_) guiDrawRefTrack_->Toggle();
    hf->AddFrame(guiDrawRefTrack_);
    guiDrawRefTrack_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiDrawErrors_ =  new TGCheckButton(hf, "Draw track errors");
    if(drawErrors_) guiDrawErrors_->Toggle();
    hf->AddFrame(guiDrawErrors_);
    guiDrawErrors_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiDrawForward_ =  new TGCheckButton(hf, "Draw forward fit");
    if(drawForward_) guiDrawForward_->Toggle();
    hf->AddFrame(guiDrawForward_);
    guiDrawForward_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiDrawBackward_ =  new TGCheckButton(hf, "Draw backward fit");
    if(drawBackward_) guiDrawBackward_->Toggle();
    hf->AddFrame(guiDrawBackward_);
    guiDrawBackward_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);


  hf = new TGHorizontalFrame(frmMain); {
    guiDrawAutoScale_ =  new TGCheckButton(hf, "Auto-scale errors");
    if(drawAutoScale_) guiDrawAutoScale_->Toggle();
    hf->AddFrame(guiDrawAutoScale_);
    guiDrawAutoScale_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiDrawScaleMan_ =  new TGCheckButton(hf, "Manually scale errors");
    if(drawScaleMan_) guiDrawScaleMan_->Toggle();
    hf->AddFrame(guiDrawScaleMan_);
    guiDrawScaleMan_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiErrorScale_ = new TGNumberEntry(hf, errorScale_, 6,999, TGNumberFormat::kNESReal,
                          TGNumberFormat::kNEANonNegative,
                          TGNumberFormat::kNELLimitMinMax,
                          1.E-4, 1.E5);
    hf->AddFrame(guiErrorScale_);
    guiErrorScale_->Connect("ValueSet(Long_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
    lbl = new TGLabel(hf, "Error scale");
    hf->AddFrame(lbl);
  }
  frmMain->AddFrame(hf);



  hf = new TGHorizontalFrame(frmMain); {
    lbl = new TGLabel(hf, "\n TrackRep options");
    hf->AddFrame(lbl);
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiDrawCardinalRep_ =  new TGCheckButton(hf, "Draw cardinal rep");
    if(drawCardinalRep_) guiDrawCardinalRep_->Toggle();
    hf->AddFrame(guiDrawCardinalRep_);
    guiDrawCardinalRep_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiRepId_ = new TGNumberEntry(hf, repId_, 6,999, TGNumberFormat::kNESInteger,
                          TGNumberFormat::kNEANonNegative,
                          TGNumberFormat::kNELLimitMinMax,
                          0, 99);
    hf->AddFrame(guiRepId_);
    guiRepId_->Connect("ValueSet(Long_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
    lbl = new TGLabel(hf, "Else draw rep with id");
    hf->AddFrame(lbl);
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiDrawAllTracks_ =  new TGCheckButton(hf, "Draw all tracks");
    if(drawAllTracks_) guiDrawAllTracks_->Toggle();
    hf->AddFrame(guiDrawAllTracks_);
    guiDrawAllTracks_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain); {
    guiTrackId_ = new TGNumberEntry(hf, trackId_, 6,999, TGNumberFormat::kNESInteger,
                          TGNumberFormat::kNEANonNegative,
                          TGNumberFormat::kNELLimitMinMax,
                          0, 99);
    hf->AddFrame(guiTrackId_);
    guiTrackId_->Connect("ValueSet(Long_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
    lbl = new TGLabel(hf, "Else draw track nr. ");
    hf->AddFrame(lbl);
  }
  frmMain->AddFrame(hf);



  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();

  browser->StopEmbedding();
  browser->SetTabTitle("Draw Control", 0);


  browser->StartEmbedding(TRootBrowser::kLeft);
  TGMainFrame* frmMain2 = new TGMainFrame(gClient->GetRoot(), 1000, 600);
  frmMain2->SetWindowName("XX GUI");
  frmMain2->SetCleanup(kDeepCleanup);

  hf = new TGHorizontalFrame(frmMain2); {
    // evt number entry
    lbl = new TGLabel(hf, "Go to event: ");
    hf->AddFrame(lbl);
    guiEvent2 = new TGNumberEntry(hf, 0, 9,999, TGNumberFormat::kNESInteger,
                          TGNumberFormat::kNEANonNegative,
                          TGNumberFormat::kNELLimitMinMax,
                          0, 99999);
    hf->AddFrame(guiEvent2);
    guiEvent2->Connect("ValueSet(Long_t)", "genfit::EventDisplay", fh, "guiGoto2()");

    // redraw button
    tb = new TGTextButton(hf, "Redraw Event");
    hf->AddFrame(tb);
    tb->Connect("Clicked()", "genfit::EventDisplay", fh, "guiGoto()");
  }
  frmMain2->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain2); {
    lbl = new TGLabel(hf, "\n Fitting options");
    hf->AddFrame(lbl);
  }
  frmMain2->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain2); {
    guiRefit_ =  new TGCheckButton(hf, "Refit");
    if(refit_) guiRefit_->Toggle();
    hf->AddFrame(guiRefit_);
    guiRefit_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain2->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain2); {
    guiDebugLvl_ = new TGNumberEntry(hf, debugLvl_, 6,999, TGNumberFormat::kNESInteger,
                          TGNumberFormat::kNEANonNegative,
                          TGNumberFormat::kNELLimitMinMax,
                          0, 999);
    hf->AddFrame(guiDebugLvl_);
    guiDebugLvl_->Connect("ValueSet(Long_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
    lbl = new TGLabel(hf, "debug level");
    hf->AddFrame(lbl);
  }
  frmMain2->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain2); {
    guiFitterId_ = new TGButtonGroup(hf,"Fitter type:");
    guiFitterId_->Connect("Clicked(Int_t)","genfit::EventDisplay", fh, "guiSelectFitterId(int)");
    hf->AddFrame(guiFitterId_, new TGLayoutHints(kLHintsTop));
      TGRadioButton* fitterId_button = new TGRadioButton(guiFitterId_, "Simple Kalman");
      new TGRadioButton(guiFitterId_, "Reference Kalman");
      new TGRadioButton(guiFitterId_, "DAF w/ simple Kalman");
      new TGRadioButton(guiFitterId_, "DAF w/ reference Kalman");
      fitterId_button->SetDown(true, false);
      guiFitterId_->Show();
  }
  frmMain2->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain2); {
    guiMmHandling_ = new TGButtonGroup(hf,"Multiple measurement handling in Kalman:");
    guiMmHandling_->Connect("Clicked(Int_t)","genfit::EventDisplay", fh, "guiSelectMmHandling(int)");
    hf->AddFrame(guiMmHandling_, new TGLayoutHints(kLHintsTop));
      TGRadioButton* mmHandling_button = new TGRadioButton(guiMmHandling_, "weighted average");
      new TGRadioButton(guiMmHandling_, "unweighted average");
      new TGRadioButton(guiMmHandling_, "weighted, closest to reference");
      new TGRadioButton(guiMmHandling_, "unweighted, closest to reference");
      new TGRadioButton(guiMmHandling_, "weighted, closest to prediction");
      new TGRadioButton(guiMmHandling_, "unweighted, closest to prediction");
      new TGRadioButton(guiMmHandling_, "weighted, closest to reference for WireMeasurements, weighted average else");
      new TGRadioButton(guiMmHandling_, "unweighted, closest to reference for WireMeasurements, unweighted average else");
      new TGRadioButton(guiMmHandling_, "weighted, closest to prediction for WireMeasurements, weighted average else");
      new TGRadioButton(guiMmHandling_, "unweighted, closest to prediction for WireMeasurements, unweighted average else");
      mmHandling_button->SetDown(true, false);
      guiMmHandling_->Show();
  }
  frmMain2->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain2); {
    guiSquareRootFormalism_ =  new TGCheckButton(hf, "Use square root formalism (simple Kalman/simple DAF)");
    if(squareRootFormalism_) guiSquareRootFormalism_->Toggle();
    hf->AddFrame(guiSquareRootFormalism_);
    guiSquareRootFormalism_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain2->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain2); {
    guiDPVal_ = new TGNumberEntry(hf, dPVal_, 6,9999, TGNumberFormat::kNESReal,
                          TGNumberFormat::kNEANonNegative,
                          TGNumberFormat::kNELLimitMinMax,
                          0, 999);
    hf->AddFrame(guiDPVal_);
    guiDPVal_->Connect("ValueSet(Long_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
    lbl = new TGLabel(hf, "delta pVal (convergence criterium)");
    hf->AddFrame(lbl);
  }
  frmMain2->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain2); {
    guiRelChi2_ = new TGNumberEntry(hf, dRelChi2_, 6,9999, TGNumberFormat::kNESReal,
                          TGNumberFormat::kNEANonNegative,
                          TGNumberFormat::kNELLimitMinMax,
                          0, 999);
    hf->AddFrame(guiRelChi2_);
    guiRelChi2_->Connect("ValueSet(Long_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
    lbl = new TGLabel(hf, "rel chi^2 change (non-convergence criterium)");
    hf->AddFrame(lbl);
  }
  frmMain2->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain2); {
    guiDChi2Ref_ = new TGNumberEntry(hf, dChi2Ref_, 6,9999, TGNumberFormat::kNESReal,
                          TGNumberFormat::kNEANonNegative,
                          TGNumberFormat::kNELLimitMinMax,
                          0, 999);
    hf->AddFrame(guiDChi2Ref_);
    guiDChi2Ref_->Connect("ValueSet(Long_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
    lbl = new TGLabel(hf, "min chi^2 change for re-calculating reference track (Ref Kalman)");
    hf->AddFrame(lbl);
  }
  frmMain2->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain2); {
    guiNMinIter_ = new TGNumberEntry(hf, nMinIter_, 6,999, TGNumberFormat::kNESInteger,
                          TGNumberFormat::kNEANonNegative,
                          TGNumberFormat::kNELLimitMinMax,
                          1, 100);
    hf->AddFrame(guiNMinIter_);
    guiNMinIter_->Connect("ValueSet(Long_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
    lbl = new TGLabel(hf, "Minimum nr of iterations");
    hf->AddFrame(lbl);
  }
  frmMain2->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain2); {
    guiNMaxIter_ = new TGNumberEntry(hf, nMaxIter_, 6,999, TGNumberFormat::kNESInteger,
                          TGNumberFormat::kNEANonNegative,
                          TGNumberFormat::kNELLimitMinMax,
                          1, 100);
    hf->AddFrame(guiNMaxIter_);
    guiNMaxIter_->Connect("ValueSet(Long_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
    lbl = new TGLabel(hf, "Maximum nr of iterations");
    hf->AddFrame(lbl);
  }
  frmMain2->AddFrame(hf);

  hf = new TGHorizontalFrame(frmMain2); {
    guiNMaxFailed_ = new TGNumberEntry(hf, nMaxFailed_, 6,999, TGNumberFormat::kNESInteger,
                          TGNumberFormat::kNEAAnyNumber,
                          TGNumberFormat::kNELLimitMinMax,
                          -1, 1000);
    hf->AddFrame(guiNMaxFailed_);
    guiNMaxFailed_->Connect("ValueSet(Long_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
    lbl = new TGLabel(hf, "Maximum nr of failed hits");
    hf->AddFrame(lbl);
  }
  frmMain2->AddFrame(hf);


  hf = new TGHorizontalFrame(frmMain2); {
    guiResort_ =  new TGCheckButton(hf, "Resort track");
    if(resort_) guiResort_->Toggle();
    hf->AddFrame(guiResort_);
    guiResort_->Connect("Toggled(Bool_t)", "genfit::EventDisplay", fh, "guiSetDrawParams()");
  }
  frmMain2->AddFrame(hf);




  frmMain2->MapSubwindows();
  frmMain2->Resize();
  frmMain2->MapWindow();

  browser->StopEmbedding();
  browser->SetTabTitle("Refit Control", 0);
}


void EventDisplay::guiGoto(){
  Long_t n = guiEvent->GetNumberEntry()->GetIntNumber();
  guiEvent2->SetIntNumber(n);
  gotoEvent(n);
}

void EventDisplay::guiGoto2(){
  Long_t n = guiEvent2->GetNumberEntry()->GetIntNumber();
  guiEvent->SetIntNumber(n);
  gotoEvent(n);
}


void EventDisplay::guiSetDrawParams(){

  drawGeometry_ = guiDrawGeometry_->IsOn();
  drawDetectors_ = guiDrawDetectors_->IsOn();
  drawHits_ = guiDrawHits_->IsOn();
  drawErrors_ = guiDrawErrors_->IsOn();

  drawPlanes_ = guiDrawPlanes_->IsOn();
  drawTrackMarkers_ = guiDrawTrackMarkers_->IsOn();
  drawTrack_ = guiDrawTrack_->IsOn();
  drawRefTrack_ = guiDrawRefTrack_->IsOn();
  drawForward_ = guiDrawForward_->IsOn();
  drawBackward_ = guiDrawBackward_->IsOn();

  drawAutoScale_ = guiDrawAutoScale_->IsOn();
  drawScaleMan_ = guiDrawScaleMan_->IsOn();

  errorScale_ = guiErrorScale_->GetNumberEntry()->GetNumber();

  drawCardinalRep_ = guiDrawCardinalRep_->IsOn();
  repId_ = guiRepId_->GetNumberEntry()->GetNumber();

  drawAllTracks_ = guiDrawAllTracks_->IsOn();
  trackId_ = guiTrackId_->GetNumberEntry()->GetNumber();


  refit_ = guiRefit_->IsOn();
  debugLvl_ = guiDebugLvl_->GetNumberEntry()->GetNumber();

  squareRootFormalism_ = guiSquareRootFormalism_->IsOn();
  dPVal_ = guiDPVal_->GetNumberEntry()->GetNumber();
  dRelChi2_ = guiRelChi2_->GetNumberEntry()->GetNumber();
  dChi2Ref_ = guiDChi2Ref_->GetNumberEntry()->GetNumber();
  nMinIter_ = guiNMinIter_->GetNumberEntry()->GetNumber();
  nMaxIter_ = guiNMaxIter_->GetNumberEntry()->GetNumber();
  nMaxFailed_ = guiNMaxFailed_->GetNumberEntry()->GetNumber();
  resort_ = guiResort_->IsOn();

  gotoEvent(eventId_);
}


void EventDisplay::guiSelectFitterId(int val){
  fitterId_ = eFitterType(val-1);
  gotoEvent(eventId_);
}

void EventDisplay::guiSelectMmHandling(int val){
  mmHandling_ = eMultipleMeasurementHandling(val-1);
  gotoEvent(eventId_);
}


} // end of namespace genfit
