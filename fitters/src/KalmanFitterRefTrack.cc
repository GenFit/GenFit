/* Copyright 2013, Ludwig-Maximilians Universität München, Technische Universität München
   Authors: Tobias Schlüter & Johannes Rauch

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/

/* This implements the Kalman fitter with reference track.  */

#include "Tools.h"
#include "Track.h"
#include "TrackPoint.h"
#include "Exception.h"
#include "IO.h"

#include "KalmanFitterRefTrack.h"
#include "KalmanFitterInfo.h"
#include "KalmanFitStatus.h"

#include <TDecompChol.h>
#include <Math/ProbFunc.h>


using namespace genfit;


TrackPoint* KalmanFitterRefTrack::fitTrack(Track* tr, const AbsTrackRep* rep, double& chi2, double& ndf, int direction)
{

  //if (!isTrackPrepared(tr, rep)) {
  //  Exception exc("KalmanFitterRefTrack::fitTrack ==> track is not properly prepared.",__LINE__,__FILE__);
  //  throw exc;
  //}

  unsigned int dim = rep->getDim();

  chi2 = 0;
  ndf = -1. * dim;
  KalmanFitterInfo* prevFi(nullptr);

  TrackPoint* retVal(nullptr);

  if (debugLvl_ > 0) {
    debugOut << tr->getNumPoints() << " TrackPoints with measurements in this track." << std::endl;
  }

  bool alreadyFitted(!refitAll_);

  p_.ResizeTo(dim);
  C_.ResizeTo(dim, dim);

  for (size_t i = 0; i < tr->getNumPointsWithMeasurement(); ++i) {
      TrackPoint *tp = 0;
      assert(direction == +1 || direction == -1);
      if (direction == +1)
        tp = tr->getPointWithMeasurement(i);
      else if (direction == -1)
        tp = tr->getPointWithMeasurement(-i-1);

      if (! tp->hasFitterInfo(rep)) {
        if (debugLvl_ > 0) {
          debugOut << "TrackPoint " << i << " has no fitterInfo -> continue. \n";
        }
        continue;
      }

      KalmanFitterInfo* fi = static_cast<KalmanFitterInfo*>(tp->getFitterInfo(rep));

      if (alreadyFitted && fi->hasUpdate(direction)) {
        if (debugLvl_ > 0) {
          debugOut << "TrackPoint " << i << " is already fitted -> continue. \n";
        }
        prevFi = fi;
        chi2 += fi->getUpdate(direction)->getChiSquareIncrement();
        ndf += fi->getUpdate(direction)->getNdf();
        continue;
      }

      alreadyFitted = false;

      if (debugLvl_ > 0) {
        debugOut << " process TrackPoint nr. " << i << "\n";
      }
      processTrackPoint(fi, prevFi, tp, chi2, ndf, direction);
      retVal = tp;

      prevFi = fi;
  }

  return retVal;
}


void KalmanFitterRefTrack::processTrackWithRep(Track* tr, const AbsTrackRep* rep, bool resortHits)
{
  if (tr->hasFitStatus(rep) && tr->getFitStatus(rep)->isTrackPruned()) {
    Exception exc("KalmanFitterRefTrack::processTrack: Cannot process pruned track!", __LINE__,__FILE__);
    throw exc;
  }

  double oldChi2FW = 1e6;
  double oldPvalFW = 0.;
  double oldChi2BW = 1e6;
  double oldPvalBW = 0.;
  double chi2FW(0), ndfFW(0);
  double chi2BW(0), ndfBW(0);
  int nFailedHits(0);

  KalmanFitStatus* status = new KalmanFitStatus();
  tr->setFitStatus(status, rep);

  status->setIsFittedWithReferenceTrack(true);

  unsigned int nIt=0;
  for (;;) {

    try {
      if (debugLvl_ > 0) {
        debugOut << " KalmanFitterRefTrack::processTrack with rep " << rep
        << " (id == " << tr->getIdForRep(rep) << ")"<< ", iteration nr. " << nIt << "\n";
      }

      // prepare
      if (!prepareTrack(tr, rep, resortHits, nFailedHits) && !refitAll_) {
        if (debugLvl_ > 0) {
          debugOut << "KalmanFitterRefTrack::processTrack. Track preparation did not change anything!\n";
        }

        status->setIsFitted();

        status->setIsFitConvergedPartially();
        if (nFailedHits == 0)
          status->setIsFitConvergedFully();
        else
          status->setIsFitConvergedFully(false);

        status->setNFailedPoints(nFailedHits);

        status->setHasTrackChanged(false);
        status->setCharge(rep->getCharge(*static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(0)->getFitterInfo(rep))->getBackwardUpdate()));
        status->setNumIterations(nIt);
        status->setForwardChi2(chi2FW);
        status->setBackwardChi2(chi2BW);
        status->setForwardNdf(std::max(0., ndfFW));
        status->setBackwardNdf(std::max(0., ndfBW));
        if (debugLvl_ > 0) {
          status->Print();
        }
        return;
      }

      if (debugLvl_ > 0) {
        debugOut << "KalmanFitterRefTrack::processTrack. Prepared Track:";
        tr->Print("C");
        //tr->Print();
      }

      // resort
      if (resortHits) {
        if (tr->sort()) {
          if (debugLvl_ > 0) {
            debugOut << "KalmanFitterRefTrack::processTrack. Resorted Track:";
            tr->Print("C");
          }
          prepareTrack(tr, rep, resortHits, nFailedHits);// re-prepare if order of hits has changed!
          status->setNFailedPoints(nFailedHits);
          if (debugLvl_ > 0) {
            debugOut << "KalmanFitterRefTrack::processTrack. Prepared resorted Track:";
            tr->Print("C");
          }
        }
      }


      // fit forward
      if (debugLvl_ > 0)
        debugOut << "KalmanFitterRefTrack::forward fit\n";
      TrackPoint* lastProcessedPoint = fitTrack(tr, rep, chi2FW, ndfFW, +1);

      // fit backward
      if (debugLvl_ > 0) {
        debugOut << "KalmanFitterRefTrack::backward fit\n";
      }

      // backward fit must not necessarily start at last hit, set prediction = forward update and blow up cov
      if (lastProcessedPoint != nullptr) {
        KalmanFitterInfo* lastInfo = static_cast<KalmanFitterInfo*>(lastProcessedPoint->getFitterInfo(rep));
        if (! lastInfo->hasBackwardPrediction()) {
          lastInfo->setBackwardPrediction(new MeasuredStateOnPlane(*(lastInfo->getForwardUpdate())));
          lastInfo->getBackwardPrediction()->blowUpCov(blowUpFactor_, resetOffDiagonals_, blowUpMaxVal_);  // blow up cov
          if (debugLvl_ > 0) {
            debugOut << "blow up cov for backward fit at TrackPoint " << lastProcessedPoint << "\n";
          }
        }
      }

      fitTrack(tr, rep, chi2BW, ndfBW, -1);

      ++nIt;


      double PvalBW = std::max(0.,ROOT::Math::chisquared_cdf_c(chi2BW, ndfBW));
      double PvalFW = (debugLvl_ > 0) ? std::max(0.,ROOT::Math::chisquared_cdf_c(chi2FW, ndfFW)) : 0; // Don't calculate if not debugging as this function potentially takes a lot of time.

      if (debugLvl_ > 0) {
        debugOut << "KalmanFitterRefTrack::Track after fit:"; tr->Print("C");

        debugOut << "old chi2s: " << oldChi2BW << ", " << oldChi2FW
            << " new chi2s: " << chi2BW << ", " << chi2FW << std::endl;
        debugOut << "old pVals: " << oldPvalBW << ", " << oldPvalFW
            << " new pVals: " << PvalBW << ", " << PvalFW << std::endl;
      }

      // See if p-value only changed little.  If the initial
      // parameters are very far off, initial chi^2 and the chi^2
      // after the first iteration will be both very close to zero, so
      // we need to have at least two iterations here.  Convergence
      // doesn't make much sense before running twice anyway.
      bool converged(false);
      bool finished(false);
      if (nIt >= minIterations_ && fabs(oldPvalBW - PvalBW) < deltaPval_)  {
        // if pVal ~ 0, check if chi2 has changed significantly
        if (fabs(1 - fabs(oldChi2BW / chi2BW)) > relChi2Change_) {
          finished = false;
        }
        else {
          finished = true;
          converged = true;
        }

        if (PvalBW == 0.)
          converged = false;
      }

      if (finished) {
        if (debugLvl_ > 0) {
          debugOut << "Fit is finished! ";
          if(converged)
            debugOut << "Fit is converged! ";
          debugOut << "\n";
        }

        if (nFailedHits == 0)
          status->setIsFitConvergedFully(converged);
        else
          status->setIsFitConvergedFully(false);

        status->setIsFitConvergedPartially(converged);
        status->setNFailedPoints(nFailedHits);

        break;
      }
      else {
        oldPvalBW = PvalBW;
        oldChi2BW = chi2BW;
        if (debugLvl_ > 0) {
          oldPvalFW = PvalFW;
          oldChi2FW = chi2FW;
        }
      }

      if (nIt >= maxIterations_) {
        if (debugLvl_ > 0) {
          debugOut << "KalmanFitterRefTrack::number of max iterations reached!\n";
        }
        break;
      }
    }
    catch(Exception& e) {
      errorOut << e.what();
      status->setIsFitted(false);
      status->setIsFitConvergedFully(false);
      status->setIsFitConvergedPartially(false);
      status->setNFailedPoints(nFailedHits);
      if (debugLvl_ > 0) {
        status->Print();
      }
      return;
    }

  }


  TrackPoint* tp = tr->getPointWithMeasurementAndFitterInfo(0, rep);

  double charge(0);
  if (tp != nullptr) {
    if (static_cast<KalmanFitterInfo*>(tp->getFitterInfo(rep))->hasBackwardUpdate())
      charge = static_cast<KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()->getCharge();
  }
  status->setCharge(charge);

  if (tp != nullptr) {
    status->setIsFitted();
  }
  else { // none of the trackPoints has a fitterInfo
    status->setIsFitted(false);
    status->setIsFitConvergedFully(false);
    status->setIsFitConvergedPartially(false);
    status->setNFailedPoints(nFailedHits);
  }

  status->setHasTrackChanged(false);
  status->setNumIterations(nIt);
  status->setForwardChi2(chi2FW);
  status->setBackwardChi2(chi2BW);
  status->setForwardNdf(ndfFW);
  status->setBackwardNdf(ndfBW);

  if (debugLvl_ > 0) {
    status->Print();
  }
}


bool KalmanFitterRefTrack::prepareTrack(Track* tr, const AbsTrackRep* rep, bool setSortingParams, int& nFailedHits) {

  if (debugLvl_ > 0) {
    debugOut << "KalmanFitterRefTrack::prepareTrack \n";
  }

  int notChangedUntil, notChangedFrom;

  // remove outdated reference states
  bool changedSmthg = removeOutdated(tr, rep,  notChangedUntil, notChangedFrom);


  // declare matrices
  FTransportMatrix_.ResizeTo(rep->getDim(), rep->getDim());
  FTransportMatrix_.UnitMatrix();
  BTransportMatrix_.ResizeTo(rep->getDim(), rep->getDim());

  FNoiseMatrix_.ResizeTo(rep->getDim(), rep->getDim());
  BNoiseMatrix_.ResizeTo(rep->getDim(), rep->getDim());

  forwardDeltaState_.ResizeTo(rep->getDim());
  backwardDeltaState_.ResizeTo(rep->getDim());

  // declare stuff
  KalmanFitterInfo* prevFitterInfo(nullptr);
  std::unique_ptr<MeasuredStateOnPlane> firstBackwardUpdate;

  ReferenceStateOnPlane* referenceState(nullptr);
  ReferenceStateOnPlane* prevReferenceState(nullptr);

  const MeasuredStateOnPlane* smoothedState(nullptr);
  const MeasuredStateOnPlane* prevSmoothedState(nullptr);

  double trackLen(0);

  bool newRefState(false); // has the current Point a new reference state?
  bool prevNewRefState(false); // has the last successfull point a new reference state?

  unsigned int nPoints = tr->getNumPoints();


  unsigned int i=0;
  nFailedHits = 0;


  // loop over TrackPoints
  for (; i<nPoints; ++i) {

    try {

      if (debugLvl_ > 0) {
        debugOut << "Prepare TrackPoint " << i << "\n";
      }

      TrackPoint* trackPoint = tr->getPoint(i);

      // check if we have a measurement
      if (!trackPoint->hasRawMeasurements()) {
        if (debugLvl_ > 0) {
          debugOut << "TrackPoint has no rawMeasurements -> continue \n";
        }
        continue;
      }

      newRefState = false;


      // get fitterInfo
      KalmanFitterInfo* fitterInfo(nullptr);
      if (trackPoint->hasFitterInfo(rep))
        fitterInfo = dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(rep));

      // create new fitter info if none available
      if (fitterInfo == nullptr) {
        if (debugLvl_ > 0) {
          debugOut << "create new KalmanFitterInfo \n";
        }
        changedSmthg = true;
        fitterInfo = new KalmanFitterInfo(trackPoint, rep);
        trackPoint->setFitterInfo(fitterInfo);
      }
      else {
        if (debugLvl_ > 0) {
          debugOut << "TrackPoint " << i << " (" << trackPoint << ") already has KalmanFitterInfo \n";
        }

        if (prevFitterInfo == nullptr) {
          if (fitterInfo->hasBackwardUpdate())
            firstBackwardUpdate.reset(new MeasuredStateOnPlane(*(fitterInfo->getBackwardUpdate())));
        }
      }

      // get smoothedState if available
      if (fitterInfo->hasPredictionsAndUpdates()) {
        smoothedState = &(fitterInfo->getFittedState(true));
        if (debugLvl_ > 0) {
          debugOut << "got smoothed state \n";
          //smoothedState->Print();
        }
      }
      else {
        smoothedState = nullptr;
      }


      if (fitterInfo->hasReferenceState()) {

        referenceState = fitterInfo->getReferenceState();


        if (!prevNewRefState) {
          if (debugLvl_ > 0) {
            debugOut << "TrackPoint already has referenceState and previous referenceState has not been altered -> continue \n";
          }
          trackLen += referenceState->getForwardSegmentLength();
          if (setSortingParams)
            trackPoint->setSortingParameter(trackLen);

          prevNewRefState = newRefState;
          prevReferenceState = referenceState;
          prevFitterInfo = fitterInfo;
          prevSmoothedState = smoothedState;

          continue;
        }


        if (prevReferenceState == nullptr) {
          if (debugLvl_ > 0) {
            debugOut << "TrackPoint already has referenceState but previous referenceState is nullptr -> reset forward info of current reference state and continue \n";
          }

          referenceState->resetForward();

          if (setSortingParams)
            trackPoint->setSortingParameter(trackLen);

          prevNewRefState = newRefState;
          prevReferenceState = referenceState;
          prevFitterInfo = fitterInfo;
          prevSmoothedState = smoothedState;

          continue;
        }

        // previous refState has been altered ->need to update transport matrices
        if (debugLvl_ > 0) {
          debugOut << "TrackPoint already has referenceState but previous referenceState has been altered -> update transport matrices and continue \n";
        }
        StateOnPlane stateToExtrapolate(*prevReferenceState);

        // make sure track is consistent if extrapolation fails
        prevReferenceState->resetBackward();
        referenceState->resetForward();

        double segmentLen = rep->extrapolateToPlane(stateToExtrapolate, fitterInfo->getReferenceState()->getPlane(), false, true);
        if (debugLvl_ > 0) {
          debugOut << "extrapolated stateToExtrapolate (prevReferenceState) by " << segmentLen << " cm.\n";
        }
        trackLen += segmentLen;

        if (segmentLen == 0) {
          FTransportMatrix_.UnitMatrix();
          FNoiseMatrix_.Zero();
          forwardDeltaState_.Zero();
          BTransportMatrix_.UnitMatrix();
          BNoiseMatrix_.Zero();
          backwardDeltaState_.Zero();
        }
        else {
          rep->getForwardJacobianAndNoise(FTransportMatrix_, FNoiseMatrix_, forwardDeltaState_);
          rep->getBackwardJacobianAndNoise(BTransportMatrix_, BNoiseMatrix_, backwardDeltaState_);
        }

        prevReferenceState->setBackwardSegmentLength(-segmentLen);
        prevReferenceState->setBackwardTransportMatrix(BTransportMatrix_);
        prevReferenceState->setBackwardNoiseMatrix(BNoiseMatrix_);
        prevReferenceState->setBackwardDeltaState(backwardDeltaState_);

        referenceState->setForwardSegmentLength(segmentLen);
        referenceState->setForwardTransportMatrix(FTransportMatrix_);
        referenceState->setForwardNoiseMatrix(FNoiseMatrix_);
        referenceState->setForwardDeltaState(forwardDeltaState_);

        newRefState = true;

        if (setSortingParams)
          trackPoint->setSortingParameter(trackLen);

        prevNewRefState = newRefState;
        prevReferenceState = referenceState;
        prevFitterInfo = fitterInfo;
        prevSmoothedState = smoothedState;

        continue;
      }


      // Construct plane
      SharedPlanePtr plane;
      if (smoothedState != nullptr) {
        if (debugLvl_ > 0)
          debugOut << "construct plane with smoothedState \n";
        plane = trackPoint->getRawMeasurement(0)->constructPlane(*smoothedState);
      }
      else if (prevSmoothedState != nullptr) {
        if (debugLvl_ > 0) {
          debugOut << "construct plane with prevSmoothedState \n";
          //prevSmoothedState->Print();
        }
        plane = trackPoint->getRawMeasurement(0)->constructPlane(*prevSmoothedState);
      }
      else if (prevReferenceState != nullptr) {
        if (debugLvl_ > 0) {
          debugOut << "construct plane with prevReferenceState \n";
        }
        plane = trackPoint->getRawMeasurement(0)->constructPlane(*prevReferenceState);
      }
      else if (rep != tr->getCardinalRep() &&
                trackPoint->hasFitterInfo(tr->getCardinalRep()) &&
                dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep())) != nullptr &&
                static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->hasPredictionsAndUpdates() ) {
        if (debugLvl_ > 0) {
          debugOut << "construct plane with smoothed state of cardinal rep fit \n";
        }
        TVector3 pos, mom;
        const MeasuredStateOnPlane& fittedState = static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->getFittedState(true);
        tr->getCardinalRep()->getPosMom(fittedState, pos, mom);
        StateOnPlane cardinalState(rep);
        rep->setPosMom(cardinalState, pos, mom);
        rep->setQop(cardinalState, tr->getCardinalRep()->getQop(fittedState));
        plane = trackPoint->getRawMeasurement(0)->constructPlane(cardinalState);
      }
      else {
        if (debugLvl_ > 0) {
          debugOut << "construct plane with state from track \n";
        }
        StateOnPlane seedFromTrack(rep);
        rep->setPosMom(seedFromTrack, tr->getStateSeed()); // also fills auxInfo
        plane = trackPoint->getRawMeasurement(0)->constructPlane(seedFromTrack);
      }

      if (plane.get() == nullptr) {
        Exception exc("KalmanFitterRefTrack::prepareTrack ==> construced plane is nullptr!",__LINE__,__FILE__);
        exc.setFatal();
        throw exc;
      }



      // do extrapolation and set reference state infos
      std::unique_ptr<StateOnPlane> stateToExtrapolate(nullptr);
      if (prevFitterInfo == nullptr) { // first measurement
        if (debugLvl_ > 0) {
          debugOut << "prevFitterInfo == nullptr \n";
        }
        if (smoothedState != nullptr) {
          if (debugLvl_ > 0) {
            debugOut << "extrapolate smoothedState to plane\n";
          }
          stateToExtrapolate.reset(new StateOnPlane(*smoothedState));
        }
        else if (referenceState != nullptr) {
          if (debugLvl_ > 0) {
            debugOut << "extrapolate referenceState to plane\n";
          }
          stateToExtrapolate.reset(new StateOnPlane(*referenceState));
        }
        else if (rep != tr->getCardinalRep() &&
                  trackPoint->hasFitterInfo(tr->getCardinalRep()) &&
                  dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep())) != nullptr &&
                  static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->hasPredictionsAndUpdates() ) {
          if (debugLvl_ > 0) {
            debugOut << "extrapolate smoothed state of cardinal rep fit to plane\n";
          }
          TVector3 pos, mom;
          const MeasuredStateOnPlane& fittedState = static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->getFittedState(true);
          stateToExtrapolate.reset(new StateOnPlane(fittedState));
          stateToExtrapolate->setRep(rep);
        }
        else {
          if (debugLvl_ > 0) {
            debugOut << "extrapolate seed from track to plane\n";
          }
          stateToExtrapolate.reset(new StateOnPlane(rep));
          rep->setPosMom(*stateToExtrapolate, tr->getStateSeed());
      	  rep->setTime(*stateToExtrapolate, tr->getTimeSeed());
        }
      } // end if (prevFitterInfo == nullptr)
      else {
        if (prevSmoothedState != nullptr) {
          if (debugLvl_ > 0) {
            debugOut << "extrapolate prevSmoothedState to plane \n";
          }
          stateToExtrapolate.reset(new StateOnPlane(*prevSmoothedState));
        }
        else {
          assert (prevReferenceState != nullptr);
          if (debugLvl_ > 0) {
            debugOut << "extrapolate prevReferenceState to plane \n";
          }
          stateToExtrapolate.reset(new StateOnPlane(*prevReferenceState));
        }
      }

      // make sure track is consistent if extrapolation fails
      if (prevReferenceState != nullptr)
        prevReferenceState->resetBackward();
      fitterInfo->deleteReferenceInfo();

      if (prevFitterInfo != nullptr) {
        rep->extrapolateToPlane(*stateToExtrapolate, prevFitterInfo->getPlane());
        if (debugLvl_ > 0) {
          debugOut << "extrapolated stateToExtrapolate to plane of prevFitterInfo (plane could have changed!) \n";
        }
      }

      double segmentLen = rep->extrapolateToPlane(*stateToExtrapolate, plane, false, true);
      trackLen += segmentLen;
      if (debugLvl_ > 0) {
        debugOut << "extrapolated stateToExtrapolate by " << segmentLen << " cm.\t";
        debugOut << "charge of stateToExtrapolate: " << rep->getCharge(*stateToExtrapolate) << " \n";
      }

      // get jacobians and noise matrices
      if (segmentLen == 0) {
        FTransportMatrix_.UnitMatrix();
        FNoiseMatrix_.Zero();
        forwardDeltaState_.Zero();
        BTransportMatrix_.UnitMatrix();
        BNoiseMatrix_.Zero();
        backwardDeltaState_.Zero();
      }
      else {
        if (i>0)
          rep->getForwardJacobianAndNoise(FTransportMatrix_, FNoiseMatrix_, forwardDeltaState_);
        rep->getBackwardJacobianAndNoise(BTransportMatrix_, BNoiseMatrix_, backwardDeltaState_);
      }


      if (i==0) {
        // if we are at first measurement and seed state is defined somewhere else
        segmentLen = 0;
        trackLen = 0;
      }
      if (setSortingParams)
        trackPoint->setSortingParameter(trackLen);


      // set backward matrices for previous reference state
      if (prevReferenceState != nullptr) {
        prevReferenceState->setBackwardSegmentLength(-segmentLen);
        prevReferenceState->setBackwardTransportMatrix(BTransportMatrix_);
        prevReferenceState->setBackwardNoiseMatrix(BNoiseMatrix_);
        prevReferenceState->setBackwardDeltaState(backwardDeltaState_);
      }


      // create new reference state
      newRefState = true;
      changedSmthg = true;
      referenceState = new ReferenceStateOnPlane(stateToExtrapolate->getState(),
             stateToExtrapolate->getPlane(),
             stateToExtrapolate->getRep(),
             stateToExtrapolate->getAuxInfo());
      referenceState->setForwardSegmentLength(segmentLen);
      referenceState->setForwardTransportMatrix(FTransportMatrix_);
      referenceState->setForwardNoiseMatrix(FNoiseMatrix_);
      referenceState->setForwardDeltaState(forwardDeltaState_);

      referenceState->resetBackward();

      fitterInfo->setReferenceState(referenceState);


      // get MeasurementsOnPlane
      std::vector<double> oldWeights = fitterInfo->getWeights();
      bool oldWeightsFixed = fitterInfo->areWeightsFixed();
      fitterInfo->deleteMeasurementInfo();
      const std::vector<AbsMeasurement*>& rawMeasurements = trackPoint->getRawMeasurements();
      for ( std::vector< genfit::AbsMeasurement* >::const_iterator measurement = rawMeasurements.begin(), lastMeasurement = rawMeasurements.end(); measurement != lastMeasurement; ++measurement) {
        assert((*measurement) != nullptr);
        fitterInfo->addMeasurementsOnPlane((*measurement)->constructMeasurementsOnPlane(*referenceState));
      }
      if (oldWeights.size() == fitterInfo->getNumMeasurements()) {
        fitterInfo->setWeights(oldWeights);
        fitterInfo->fixWeights(oldWeightsFixed);
      }


      // if we made it here, no Exceptions were thrown and the TrackPoint could successfully be processed
      prevNewRefState = newRefState;
      prevReferenceState = referenceState;
      prevFitterInfo = fitterInfo;
      prevSmoothedState = smoothedState;

    }
    catch (Exception& e) {

      if (debugLvl_ > 0) {
        errorOut << "exception at hit " << i << "\n";
        debugOut << e.what();
      }


      ++nFailedHits;
      if (maxFailedHits_<0 || nFailedHits <= maxFailedHits_) {
        prevNewRefState = true;
        referenceState = nullptr;
        smoothedState = nullptr;
        tr->getPoint(i)->deleteFitterInfo(rep);

        if (setSortingParams)
          tr->getPoint(i)->setSortingParameter(trackLen);

        if (debugLvl_ > 0) {
          debugOut << "There was an exception, try to continue with next TrackPoint " << i+1 << " \n";
        }

        continue;
      }


      // clean up
      removeForwardBackwardInfo(tr, rep, notChangedUntil, notChangedFrom);

      // set sorting parameters of rest of TrackPoints and remove FitterInfos
      for (; i<nPoints; ++i) {
        TrackPoint* trackPoint = tr->getPoint(i);

        if (setSortingParams)
          trackPoint->setSortingParameter(trackLen);

        trackPoint->deleteFitterInfo(rep);
      }
      return true;

    }

  } // end loop over TrackPoints




  removeForwardBackwardInfo(tr, rep, notChangedUntil, notChangedFrom);

  if (firstBackwardUpdate && tr->getPointWithMeasurementAndFitterInfo(0, rep)) {
    KalmanFitterInfo* fi = static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurementAndFitterInfo(0, rep)->getFitterInfo(rep));
    if (fi && ! fi->hasForwardPrediction()) {
      if (debugLvl_ > 0) {
        debugOut << "set backwards update of first point as forward prediction (with blown up cov) \n";
      }
      if (fi->getPlane() != firstBackwardUpdate->getPlane()) {
        rep->extrapolateToPlane(*firstBackwardUpdate, fi->getPlane());
      }
      firstBackwardUpdate->blowUpCov(blowUpFactor_, resetOffDiagonals_, blowUpMaxVal_);
      fi->setForwardPrediction(new MeasuredStateOnPlane(*firstBackwardUpdate));
    }
  }

  KalmanFitStatus* fitStatus = dynamic_cast<KalmanFitStatus*>(tr->getFitStatus(rep));
  if (fitStatus != nullptr)
    fitStatus->setTrackLen(trackLen);

  if (debugLvl_ > 0) {
    debugOut << "trackLen of reference track = " << trackLen << "\n";
  }

  // self check
  //assert(tr->checkConsistency());
  assert(isTrackPrepared(tr, rep));

  return changedSmthg;
}


bool
KalmanFitterRefTrack::removeOutdated(Track* tr, const AbsTrackRep* rep, int& notChangedUntil, int& notChangedFrom) {

  if (debugLvl_ > 0) {
    debugOut << "KalmanFitterRefTrack::removeOutdated \n";
  }

  bool changedSmthg(false);

  unsigned int nPoints = tr->getNumPoints();
  notChangedUntil = nPoints-1;
  notChangedFrom = 0;

  // loop over TrackPoints
  for (unsigned int i=0; i<nPoints; ++i) {

    TrackPoint* trackPoint = tr->getPoint(i);

    // check if we have a measurement
    if (!trackPoint->hasRawMeasurements()) {
      if (debugLvl_ > 0) {
        debugOut << "TrackPoint has no rawMeasurements -> continue \n";
      }
      continue;
    }

    // get fitterInfo
    KalmanFitterInfo* fitterInfo(nullptr);
    if (trackPoint->hasFitterInfo(rep))
      fitterInfo = dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(rep));

    if (fitterInfo == nullptr)
      continue;


    // check if we need to calculate or update reference state
    if (fitterInfo->hasReferenceState()) {

      if (! fitterInfo->hasPredictionsAndUpdates()) {
        if (debugLvl_ > 0) {
          debugOut << "reference state but not all predictions & updates -> do not touch reference state. \n";
        }
        continue;
      }


      const MeasuredStateOnPlane& smoothedState = fitterInfo->getFittedState(true);
      resM_.ResizeTo(smoothedState.getState());
      resM_ = smoothedState.getState();
      resM_ -= fitterInfo->getReferenceState()->getState();
      double chi2(0);

      // calculate chi2, ignore off diagonals
      double* resArray = resM_.GetMatrixArray();
      for (int j=0; j<smoothedState.getCov().GetNcols(); ++j)
        chi2 += resArray[j]*resArray[j] / smoothedState.getCov()(j,j);

      if (chi2 < deltaChi2Ref_) {
        // reference state is near smoothed state ->  do not update reference state
        if (debugLvl_ > 0) {
          debugOut << "reference state is near smoothed state ->  do not update reference state, chi2 = " << chi2 << "\n";
        }
        continue;
      } else {
        if (debugLvl_ > 0) {
          debugOut << "reference state is not close to smoothed state, chi2 = " << chi2 << "\n";
        }
      }
    }

    if (debugLvl_ > 0) {
      debugOut << "remove reference info \n";
    }

    fitterInfo->deleteReferenceInfo();
    changedSmthg = true;

    // decided to update reference state -> set notChangedUntil (only once)
    if (notChangedUntil == (int)nPoints-1)
      notChangedUntil = i-1;

    notChangedFrom = i+1;

  } // end loop over TrackPoints


  if (debugLvl_ > 0) {
    tr->Print("C");
  }

  return changedSmthg;
}


void
KalmanFitterRefTrack::removeForwardBackwardInfo(Track* tr, const AbsTrackRep* rep, int notChangedUntil, int notChangedFrom) const {

  unsigned int nPoints = tr->getNumPoints();

  if (refitAll_) {
    tr->deleteForwardInfo(0, -1, rep);
    tr->deleteBackwardInfo(0, -1, rep);
    return;
  }

  // delete forward/backward info from/to points where reference states have changed
  if (notChangedUntil != (int)nPoints-1) {
    tr->deleteForwardInfo(notChangedUntil+1, -1, rep);
  }
  if (notChangedFrom != 0)
    tr->deleteBackwardInfo(0, notChangedFrom-1, rep);

}


void
KalmanFitterRefTrack::processTrackPoint(KalmanFitterInfo* fi, const KalmanFitterInfo* prevFi, const TrackPoint* tp, double& chi2, double& ndf, int direction)
{
  if(squareRootFormalism_) {
    processTrackPointSqrt(fi, prevFi, tp, chi2, ndf, direction);
    return;
  }

  if (debugLvl_ > 0) {
    debugOut << " KalmanFitterRefTrack::processTrackPoint " << fi->getTrackPoint() << "\n";
  }

  unsigned int dim = fi->getRep()->getDim();

  p_.Zero(); // p_{k|k-1}
  C_.Zero(); // C_{k|k-1}

  // predict
  if (prevFi != nullptr) {
    const TMatrixD& F = fi->getReferenceState()->getTransportMatrix(direction); // Transport matrix
    assert(F.GetNcols() == (int)dim);
    const TMatrixDSym& N = fi->getReferenceState()->getNoiseMatrix(direction); // Noise matrix
    //p_ = ( F * prevFi->getUpdate(direction)->getState() ) + fi->getReferenceState()->getDeltaState(direction);
    p_ = prevFi->getUpdate(direction)->getState();
    p_ *= F;
    p_ += fi->getReferenceState()->getDeltaState(direction);

    C_ = prevFi->getUpdate(direction)->getCov();
    C_.Similarity(F);
    C_ += N;
    fi->setPrediction(new MeasuredStateOnPlane(p_, C_, fi->getReferenceState()->getPlane(), fi->getReferenceState()->getRep(), fi->getReferenceState()->getAuxInfo()), direction);
    if (debugLvl_ > 1) {
      debugOut << "\033[31m";
      debugOut << "F (Transport Matrix) "; F.Print();
      debugOut << "p_{k,r} (reference state) "; fi->getReferenceState()->getState().Print();
      debugOut << "c (delta state) "; fi->getReferenceState()->getDeltaState(direction).Print();
      debugOut << "F*p_{k-1,r} + c "; (F *prevFi->getReferenceState()->getState() + fi->getReferenceState()->getDeltaState(direction)).Print();
    }
  }
  else {
    if (fi->hasPrediction(direction)) {
      if (debugLvl_ > 0) {
        debugOut << "  Use prediction as start \n";
      }
      p_ = fi->getPrediction(direction)->getState();
      C_ = fi->getPrediction(direction)->getCov();
    }
    else {
      if (debugLvl_ > 0) {
        debugOut << "  Use reference state and seed cov as start \n";
      }
      const AbsTrackRep *rep = fi->getReferenceState()->getRep();
      p_ = fi->getReferenceState()->getState();

      // Convert from 6D covariance of the seed to whatever the trackRep wants.
      TMatrixDSym dummy(p_.GetNrows());
      MeasuredStateOnPlane mop(p_, dummy, fi->getReferenceState()->getPlane(), rep, fi->getReferenceState()->getAuxInfo());
      TVector3 pos, mom;
      rep->getPosMom(mop, pos, mom);
      rep->setPosMomCov(mop, pos, mom, fi->getTrackPoint()->getTrack()->getCovSeed());
      // Blow up, set.
      mop.blowUpCov(blowUpFactor_, resetOffDiagonals_, blowUpMaxVal_);
      fi->setPrediction(new MeasuredStateOnPlane(mop), direction);
      C_ = mop.getCov();
    }
    if (debugLvl_ > 1) {
      debugOut << "\033[31m";
      debugOut << "p_{k,r} (reference state)"; fi->getReferenceState()->getState().Print();
    }
  }

  if (debugLvl_ > 1) {
    debugOut << " p_{k|k-1} (state prediction)"; p_.Print();
    debugOut << " C_{k|k-1} (covariance prediction)"; C_.Print();
    debugOut << "\033[0m";
  }

  // update(s)
  double chi2inc = 0;
  double ndfInc = 0;
  const std::vector<MeasurementOnPlane *> measurements = getMeasurements(fi, tp, direction);
  for (std::vector<MeasurementOnPlane *>::const_iterator it = measurements.begin(); it != measurements.end(); ++it) {
    const MeasurementOnPlane& m = **it;

    if (!canIgnoreWeights() && m.getWeight() <= 1.01E-10) {
      if (debugLvl_ > 1) {
        debugOut << "Weight of measurement is almost 0, continue ... /n";
      }
      continue;
    }

    const AbsHMatrix* H(m.getHMatrix());
    // (weighted) cov
    const TMatrixDSym& V((!canIgnoreWeights() && m.getWeight() < 0.99999) ?
                          1./m.getWeight() * m.getCov() :
                          m.getCov());

    covSumInv_.ResizeTo(C_);
    covSumInv_ = C_; // (V_k + H_k C_{k|k-1} H_k^T)^(-1)
    H->HMHt(covSumInv_);
    covSumInv_ += V;

    tools::invertMatrix(covSumInv_);

    const TMatrixD& CHt(H->MHt(C_));

    res_.ResizeTo(m.getState());
    res_ = m.getState();
    res_ -= H->Hv(p_); // residual
    if (debugLvl_ > 1) {
      debugOut << "\033[34m";
      debugOut << "m (measurement) "; m.getState().Print();
      debugOut << "V ((weighted) measurement covariance) "; (1./m.getWeight() * m.getCov()).Print();
      debugOut << "residual        "; res_.Print();
      debugOut << "\033[0m";
    }
    p_ +=  TMatrixD(CHt, TMatrixD::kMult, covSumInv_) * res_; // updated state
    if (debugLvl_ > 1) {
      debugOut << "\033[32m";
      debugOut << " update"; (TMatrixD(CHt, TMatrixD::kMult, covSumInv_) * res_).Print();
      debugOut << "\033[0m";
    }

    covSumInv_.Similarity(CHt); // with (C H^T)^T = H C^T = H C  (C is symmetric)
    C_ -= covSumInv_; // updated Cov

    if (debugLvl_ > 1) {
      //debugOut << " C update "; covSumInv_.Print();
      debugOut << "\033[32m";
      debugOut << " p_{k|k} (updated state)"; p_.Print();
      debugOut << " C_{k|k} (updated covariance)"; C_.Print();
      debugOut << "\033[0m";
    }

    // Calculate chi² increment.  At the first point chi2inc == 0 and
    // the matrix will not be invertible.
    res_ = m.getState();
    res_ -= H->Hv(p_); // new residual
    if (debugLvl_ > 1) {
      debugOut << " resNew ";
      res_.Print();
    }

    // only calculate chi2inc if res != nullptr.
    // If matrix inversion fails, chi2inc = 0
    if (res_ != 0) {
      Rinv_.ResizeTo(C_);
      Rinv_ = C_;
      H->HMHt(Rinv_);
      Rinv_ -= V;
      Rinv_ *= -1;

      bool couldInvert(true);
      try {
        tools::invertMatrix(Rinv_);
      }
      catch (Exception& e) {
        if (debugLvl_ > 1) {
          debugOut << e.what();
        }
        couldInvert = false;
      }

      if (couldInvert) {
        if (debugLvl_ > 1) {
          debugOut << " Rinv ";
          Rinv_.Print();
        }
        chi2inc += Rinv_.Similarity(res_);
      }
    }


    if (!canIgnoreWeights()) {
      ndfInc += m.getWeight() * m.getState().GetNrows();
    }
    else
      ndfInc += m.getState().GetNrows();


  } // end loop over measurements

  chi2 += chi2inc;
  ndf += ndfInc;


  KalmanFittedStateOnPlane* upState = new KalmanFittedStateOnPlane(p_, C_, fi->getReferenceState()->getPlane(), fi->getReferenceState()->getRep(), fi->getReferenceState()->getAuxInfo(), chi2inc, ndfInc);
  upState->setAuxInfo(fi->getReferenceState()->getAuxInfo());
  fi->setUpdate(upState, direction);


  if (debugLvl_ > 0) {
    debugOut << " chi² inc " << chi2inc << "\t";
    debugOut << " ndf inc  " << ndfInc << "\t";
    debugOut << " charge of update  " << fi->getRep()->getCharge(*upState) << "\n";
  }

  // check
  if (not fi->checkConsistency()) {
    throw genfit::Exception("Consistency check failed ", __LINE__, __FILE__);
  }

}




void
KalmanFitterRefTrack::processTrackPointSqrt(KalmanFitterInfo* fi, const KalmanFitterInfo* prevFi,
					    const TrackPoint* tp, double& chi2, double& ndf, int direction)
{
  if (debugLvl_ > 0) {
    debugOut << " KalmanFitterRefTrack::processTrackPointSqrt " << fi->getTrackPoint() << "\n";
  }

  unsigned int dim = fi->getRep()->getDim();

  p_.Zero(); // p_{k|k-1}
  C_.Zero(); // C_{k|k-1}

  TMatrixD S(dim, dim); // sqrt(C_);

  // predict
  if (prevFi != nullptr) {
    const TMatrixD& F = fi->getReferenceState()->getTransportMatrix(direction); // Transport matrix
    assert(F.GetNcols() == (int)dim);
    const TMatrixDSym& N = fi->getReferenceState()->getNoiseMatrix(direction); // Noise matrix
    //N = 0;

    //p_ = ( F * prevFi->getUpdate(direction)->getState() ) + fi->getReferenceState()->getDeltaState(direction);
    p_ = prevFi->getUpdate(direction)->getState();
    p_ *= F;
    p_ += fi->getReferenceState()->getDeltaState(direction);


    TDecompChol decompS(prevFi->getUpdate(direction)->getCov());
    decompS.Decompose();
    TMatrixD Q;
    tools::noiseMatrixSqrt(N, Q);
    tools::kalmanPredictionCovSqrt(decompS.GetU(), F, Q, S);

    fi->setPrediction(new MeasuredStateOnPlane(p_, TMatrixDSym(TMatrixDSym::kAtA, S), fi->getReferenceState()->getPlane(), fi->getReferenceState()->getRep(), fi->getReferenceState()->getAuxInfo()), direction);
    if (debugLvl_ > 1) {
      debugOut << "\033[31m";
      debugOut << "F (Transport Matrix) "; F.Print();
      debugOut << "p_{k,r} (reference state) "; fi->getReferenceState()->getState().Print();
      debugOut << "c (delta state) "; fi->getReferenceState()->getDeltaState(direction).Print();
      debugOut << "F*p_{k-1,r} + c "; (F *prevFi->getReferenceState()->getState() + fi->getReferenceState()->getDeltaState(direction)).Print();
    }
  }
  else {
    if (fi->hasPrediction(direction)) {
      if (debugLvl_ > 0) {
        debugOut << "  Use prediction as start \n";
      }
      p_ = fi->getPrediction(direction)->getState();
      TDecompChol decompS(fi->getPrediction(direction)->getCov());
      decompS.Decompose();
      S = decompS.GetU();
    }
    else {
      if (debugLvl_ > 0) {
        debugOut << "  Use reference state and seed cov as start \n";
      }
      const AbsTrackRep *rep = fi->getReferenceState()->getRep();
      p_ = fi->getReferenceState()->getState();

      // Convert from 6D covariance of the seed to whatever the trackRep wants.
      TMatrixDSym dummy(p_.GetNrows());
      MeasuredStateOnPlane mop(p_, dummy, fi->getReferenceState()->getPlane(), rep, fi->getReferenceState()->getAuxInfo());
      TVector3 pos, mom;
      rep->getPosMom(mop, pos, mom);
      rep->setPosMomCov(mop, pos, mom, fi->getTrackPoint()->getTrack()->getCovSeed());
      // Blow up, set.
      mop.blowUpCov(blowUpFactor_, resetOffDiagonals_, blowUpMaxVal_);
      fi->setPrediction(new MeasuredStateOnPlane(mop), direction);
      TDecompChol decompS(mop.getCov());
      decompS.Decompose();
      S = decompS.GetU();
    }
    if (debugLvl_ > 1) {
      debugOut << "\033[31m";
      debugOut << "p_{k,r} (reference state)"; fi->getReferenceState()->getState().Print();
    }
  }

  if (debugLvl_ > 1) {
    debugOut << " p_{k|k-1} (state prediction)"; p_.Print();
    debugOut << " C_{k|k-1} (covariance prediction)"; C_.Print();
    debugOut << "\033[0m";
  }

  // update(s)
  double chi2inc = 0;
  double ndfInc = 0;

  const std::vector<MeasurementOnPlane *> measurements = getMeasurements(fi, tp, direction);
  for (std::vector<MeasurementOnPlane *>::const_iterator it = measurements.begin(); it != measurements.end(); ++it) {
    const MeasurementOnPlane& m = **it;

    if (!canIgnoreWeights() && m.getWeight() <= 1.01E-10) {
      if (debugLvl_ > 1) {
        debugOut << "Weight of measurement is almost 0, continue ... /n";
      }
      continue;
    }

    const AbsHMatrix* H(m.getHMatrix());
    // (weighted) cov
    const TMatrixDSym& V((!canIgnoreWeights() && m.getWeight() < 0.99999) ?
                          1./m.getWeight() * m.getCov() :
                          m.getCov());
    TDecompChol decompR(V);
    decompR.Decompose();
    const TMatrixD& R(decompR.GetU());

    res_.ResizeTo(m.getState());
    res_ = m.getState();
    res_ -= H->Hv(p_); // residual

    TVectorD update;
    tools::kalmanUpdateSqrt(S, res_, R, H,
			    update, S);

    if (debugLvl_ > 1) {
      debugOut << "\033[34m";
      debugOut << "m (measurement) "; m.getState().Print();
      debugOut << "V ((weighted) measurement covariance) "; (1./m.getWeight() * m.getCov()).Print();
      debugOut << "residual        "; res_.Print();
      debugOut << "\033[0m";
    }

    p_ += update;
    if (debugLvl_ > 1) {
      debugOut << "\033[32m";
      debugOut << " update"; update.Print();
      debugOut << "\033[0m";
    }

    if (debugLvl_ > 1) {
      //debugOut << " C update "; covSumInv_.Print();
      debugOut << "\033[32m";
      debugOut << " p_{k|k} (updated state)"; p_.Print();
      debugOut << " C_{k|k} (updated covariance)"; C_.Print();
      debugOut << "\033[0m";
    }

    // Calculate chi² increment.  At the first point chi2inc == 0 and
    // the matrix will not be invertible.
    res_ -= H->Hv(update); // new residual
    if (debugLvl_ > 1) {
      debugOut << " resNew ";
      res_.Print();
    }

    // only calculate chi2inc if res != 0.
    // If matrix inversion fails, chi2inc = 0
    if (res_ != 0) {
      Rinv_.ResizeTo(V);
      Rinv_ = V - TMatrixDSym(TMatrixDSym::kAtA, H->MHt(S));

      bool couldInvert(true);
      try {
        tools::invertMatrix(Rinv_);
      }
      catch (Exception& e) {
        if (debugLvl_ > 1) {
          debugOut << e.what();
        }
        couldInvert = false;
      }

      if (couldInvert) {
        if (debugLvl_ > 1) {
          debugOut << " Rinv ";
          Rinv_.Print();
        }
        chi2inc += Rinv_.Similarity(res_);
      }
    }

    if (!canIgnoreWeights()) {
      ndfInc += m.getWeight() * m.getState().GetNrows();
    }
    else
      ndfInc += m.getState().GetNrows();


  } // end loop over measurements

  C_ = TMatrixDSym(TMatrixDSym::kAtA, S);

  chi2 += chi2inc;
  ndf += ndfInc;


  KalmanFittedStateOnPlane* upState = new KalmanFittedStateOnPlane(p_, C_, fi->getReferenceState()->getPlane(), fi->getReferenceState()->getRep(), fi->getReferenceState()->getAuxInfo(), chi2inc, ndfInc);
  upState->setAuxInfo(fi->getReferenceState()->getAuxInfo());
  fi->setUpdate(upState, direction);


  if (debugLvl_ > 0) {
    debugOut << " chi² inc " << chi2inc << "\t";
    debugOut << " ndf inc  " << ndfInc << "\t";
    debugOut << " charge of update  " << fi->getRep()->getCharge(*upState) << "\n";
  }

  // check
  if (not fi->checkConsistency()) {
    throw genfit::Exception("Consistency check failed ", __LINE__, __FILE__);
  }

}
