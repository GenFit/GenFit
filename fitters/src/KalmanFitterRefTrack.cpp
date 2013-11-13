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

#include "KalmanFitterRefTrack.h"
#include "KalmanFitterInfo.h"
#include "KalmanFitStatus.h"

#include "boost/scoped_ptr.hpp"

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
  KalmanFitterInfo* prevFi(NULL);

  TrackPoint* retVal(NULL);

  if (debugLvl_ > 0)
    std::cout << tr->getNumPoints() << " TrackPoints with measurements in this track." << std::endl;

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
        if (debugLvl_ > 0)
          std::cout << "TrackPoint " << i << " has no fitterInfo -> continue. \n";
        continue;
      }

      KalmanFitterInfo* fi = static_cast<KalmanFitterInfo*>(tp->getFitterInfo(rep));

      if (alreadyFitted && fi->hasUpdate(direction)) {
        if (debugLvl_ > 0)
          std::cout << "TrackPoint " << i << " is already fitted -> continue. \n";
        prevFi = fi;
        chi2 += fi->getUpdate(direction)->getChiSquareIncrement();
        ndf += fi->getUpdate(direction)->getNdf();
        continue;
      }

      alreadyFitted = false;

      if (debugLvl_ > 0)
        std::cout << " process TrackPoint nr. " << i << "\n";
      processTrackPoint(fi, prevFi, tp, chi2, ndf, direction);
      retVal = tp;

      prevFi = fi;
  }

  return retVal;
}


void KalmanFitterRefTrack::processTrackWithRep(Track* tr, const AbsTrackRep* rep, bool resortHits)
{

  if (tr->getFitStatus(rep) != NULL && tr->getFitStatus(rep)->isTrackPruned()) {
    Exception exc("KalmanFitterRefTrack::processTrack: Cannot process pruned track!", __LINE__,__FILE__);
    throw exc;
  }

  double oldChi2FW = 1e6;
  double oldPvalFW = 0.;
  double oldChi2BW = 1e6;
  double oldPvalBW = 0.;
  double chi2FW(0), ndfFW(0);
  double chi2BW(0), ndfBW(0);

  KalmanFitStatus* status = new KalmanFitStatus();
  tr->setFitStatus(status, rep);

  status->setIsFittedWithReferenceTrack(true);

  unsigned int nIt=0;
  for (;;) {

    try {
      if (debugLvl_ > 0)
        std::cout << " KalmanFitterRefTrack::processTrack with rep " << rep
        << " (id == " << tr->getIdForRep(rep) << ")"<< ", iteration nr. " << nIt << "\n";

      // prepare
      if (!prepareTrack(tr, rep, resortHits) && !refitAll_) {
        if (debugLvl_ > 0)
          std::cout << "KalmanFitterRefTrack::processTrack. Track preparation did not change anything!\n";
        status->setIsFitted();
        status->setIsFitConverged();
        status->setHasTrackChanged(false);
        status->setCharge(rep->getCharge(*static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(0)->getFitterInfo(rep))->getBackwardUpdate()));
        status->setNumIterations(nIt);
        status->setForwardChi2(chi2FW);
        status->setBackwardChi2(chi2BW);
        status->setForwardNdf(std::max(0., ndfFW));
        status->setBackwardNdf(std::max(0., ndfBW));
        if (debugLvl_ > 0)
          status->Print();
        return;
      }

      if (debugLvl_ > 0) {
        std::cout << "KalmanFitterRefTrack::processTrack. Prepared Track:";
        tr->Print("C");
        //tr->Print();
      }

      // resort
      if (resortHits) {
        if (tr->sort()) {
          if (debugLvl_ > 0) {
            std::cout << "KalmanFitterRefTrack::processTrack. Resorted Track:";
            tr->Print("C");
          }
          prepareTrack(tr, rep, resortHits);// re-prepare if order of hits has changed!
          if (debugLvl_ > 0) {
            std::cout << "KalmanFitterRefTrack::processTrack. Prepared resorted Track:";
            tr->Print("C");
          }
        }
      }


      // fit forward
      if (debugLvl_ > 0)
        std::cout << "KalmanFitterRefTrack::forward fit\n";
      TrackPoint* lastProcessedPoint = fitTrack(tr, rep, chi2FW, ndfFW, +1);

      // fit backward
      if (debugLvl_ > 0)
        std::cout << "KalmanFitterRefTrack::backward fit\n";

      // backward fit must not necessarily start at last hit, set prediction = forward update and blow up cov
      if (lastProcessedPoint != NULL) {
        KalmanFitterInfo* lastInfo = static_cast<KalmanFitterInfo*>(lastProcessedPoint->getFitterInfo(rep));
        if (! lastInfo->hasBackwardPrediction()) {
          lastInfo->setBackwardPrediction(new MeasuredStateOnPlane(*(lastInfo->getForwardUpdate())));
          lastInfo->getBackwardPrediction()->blowUpCov(blowUpFactor_);  // blow up cov
          if (debugLvl_ > 0)
            std::cout << "blow up cov for backward fit at TrackPoint " << lastProcessedPoint << "\n";
        }
      }

      fitTrack(tr, rep, chi2BW, ndfBW, -1);

      ++nIt;


      double PvalBW = ROOT::Math::chisquared_cdf_c(chi2BW, ndfBW);
      double PvalFW = (debugLvl_ > 0) ? ROOT::Math::chisquared_cdf_c(chi2FW, ndfFW) : 0; // Don't calculate if not debugging as this function potentially takes a lot of time.

      if (debugLvl_ > 0) {
        std::cout << "KalmanFitterRefTrack::Track after fit:"; tr->Print("C");

        std::cout << "old chi2s: " << oldChi2BW << ", " << oldChi2FW
            << " new chi2s: " << chi2BW << ", " << chi2FW << std::endl;
        std::cout << "old pVals: " << oldPvalBW << ", " << oldPvalFW
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
          std::cout << "Fit is finished! ";
          if(converged)
            std::cout << "Fit is converged! ";
          std::cout << "\n";
        }
        status->setIsFitConverged(converged);
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
        if (debugLvl_ > 0)
          std::cout << "KalmanFitterRefTrack::number of max iterations reached!\n";
        break;
      }
    }
    catch(Exception& e) {
      std::cerr << e.what();
      status->setIsFitted(false);
      status->setIsFitConverged(false);
      if (debugLvl_ > 0)
        status->Print();
      return;
    }

  }


  TrackPoint* tp = tr->getPointWithMeasurementAndFitterInfo(0, rep);
  if (tp != NULL &&
      static_cast<KalmanFitterInfo*>(tp->getFitterInfo(rep))->hasBackwardUpdate())
    status->setCharge(rep->getCharge(*static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(0)->getFitterInfo(rep))->getBackwardUpdate()));

  if (tp != NULL) {
    status->setIsFitted();
  }
  else { // none of the trackPoints has a fitterInfo
    status->setIsFitted(false);
    status->setIsFitConverged(false);
  }

  status->setHasTrackChanged(false);
  status->setNumIterations(nIt);
  status->setForwardChi2(chi2FW);
  status->setBackwardChi2(chi2BW);
  status->setForwardNdf(ndfFW);
  status->setBackwardNdf(ndfBW);

  if (debugLvl_ > 0)
    status->Print();
}


bool KalmanFitterRefTrack::prepareTrack(Track* tr, const AbsTrackRep* rep, bool setSortingParams) {

  if (debugLvl_ > 0)
    std::cout << "KalmanFitterRefTrack::prepareTrack \n";

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
  KalmanFitterInfo* prevFitterInfo(NULL);
  boost::scoped_ptr<MeasuredStateOnPlane> firstBackwardUpdate;

  ReferenceStateOnPlane* referenceState(NULL);
  ReferenceStateOnPlane* prevReferenceState(NULL);

  const MeasuredStateOnPlane* smoothedState(NULL);
  const MeasuredStateOnPlane* prevSmoothedState(NULL);

  double trackLen(0);

  bool newRefState(false);

  unsigned int nPoints = tr->getNumPoints();


  unsigned int i=0;

  try {

    // loop over TrackPoints
    for (; i<nPoints; ++i) {

      TrackPoint* trackPoint = tr->getPoint(i);

      // check if we have a measurement
      if (!trackPoint->hasRawMeasurements()) {
        if (debugLvl_ > 0)
          std::cout << "TrackPoint has no rawMeasurements -> continue \n";
        continue;
      }


      // get fitterInfo
      KalmanFitterInfo* fitterInfo(NULL);
      if (trackPoint->hasFitterInfo(rep))
        fitterInfo = dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(rep));

      // create new fitter info if none available
      if (fitterInfo == NULL) {
        if (debugLvl_ > 0)
          std::cout << "create new KalmanFitterInfo \n";
        fitterInfo = new KalmanFitterInfo(trackPoint, rep);
        trackPoint->setFitterInfo(fitterInfo);
        changedSmthg = true;
      }
      else {
        if (debugLvl_ > 0)
          std::cout << "TrackPoint " << i << " (" << trackPoint << ") already has KalmanFitterInfo \n";

        if (prevFitterInfo == NULL) {
          if (fitterInfo->hasBackwardUpdate())
            firstBackwardUpdate.reset(new MeasuredStateOnPlane(*(fitterInfo->getBackwardUpdate())));
        }
      }

      // get smoothedState if available
      if (fitterInfo->hasPredictionsAndUpdates()) {
        smoothedState = &(fitterInfo->getFittedState(true));
      }
      else {
        smoothedState = NULL;
      }


      if (fitterInfo->hasReferenceState()) {

        referenceState = fitterInfo->getReferenceState();
        prevFitterInfo = fitterInfo;
        prevSmoothedState = smoothedState;

        if (!newRefState) {
          if (debugLvl_ > 0)
            std::cout << "TrackPoint already has referenceState and previous referenceState has not been altered -> continue \n";
          prevReferenceState = referenceState;
          trackLen += referenceState->getForwardSegmentLength();
          if (setSortingParams)
            trackPoint->setSortingParameter(trackLen);
          continue;
        }

        // previous refState has been altered ->need to update transport matrices
        if (debugLvl_ > 0)
          std::cout << "TrackPoint already has referenceState but previous referenceState has been altered -> update transport matrices and continue \n";
        StateOnPlane stateToExtrapolate(*prevReferenceState);

        // make sure track is consistent if extrapolation fails
        prevReferenceState->resetBackward();
        referenceState->resetForward();

        double segmentLen = rep->extrapolateToPlane(stateToExtrapolate, fitterInfo->getReferenceState()->getPlane(), false, true);
        if (debugLvl_ > 0)
          std::cout << "extrapolated stateToExtrapolate (prevReferenceState) by " << segmentLen << " cm.\n";
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

        if (setSortingParams)
          trackPoint->setSortingParameter(trackLen);


        prevReferenceState = referenceState;
        newRefState = false;

        continue;
      }

      newRefState = false;


      // Construct plane
      SharedPlanePtr plane;
      if (smoothedState != NULL) {
        if (debugLvl_ > 0)
          std::cout << "construct plane with smoothedState \n";
        plane = trackPoint->getRawMeasurement(0)->constructPlane(*smoothedState);
      }
      else if (prevSmoothedState != NULL) {
        if (debugLvl_ > 0)
          std::cout << "construct plane with prevSmoothedState \n";
        plane = trackPoint->getRawMeasurement(0)->constructPlane(*prevSmoothedState);
      }
      else if (prevReferenceState != NULL) {
        if (debugLvl_ > 0)
          std::cout << "construct plane with prevReferenceState \n";
        plane = trackPoint->getRawMeasurement(0)->constructPlane(*prevReferenceState);
      }
      else if (rep != tr->getCardinalRep() &&
                trackPoint->hasFitterInfo(tr->getCardinalRep()) &&
                dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep())) != NULL &&
                static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->hasPredictionsAndUpdates() ) {
        if (debugLvl_ > 0)
          std::cout << "construct plane with smoothed state of cardinal rep fit \n";
        TVector3 pos, mom;
        const MeasuredStateOnPlane& fittedState = static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->getFittedState(true);
        tr->getCardinalRep()->getPosMom(fittedState, pos, mom);
        StateOnPlane cardinalState(rep);
        rep->setPosMom(cardinalState, pos, mom);
        rep->setQop(cardinalState, tr->getCardinalRep()->getQop(fittedState));
        plane = trackPoint->getRawMeasurement(0)->constructPlane(cardinalState);
      }
      else {
        if (debugLvl_ > 0)
          std::cout << "construct plane with state from track \n";
        StateOnPlane seedFromTrack(rep);
        rep->setPosMom(seedFromTrack, tr->getStateSeed()); // also fills auxInfo
        plane = trackPoint->getRawMeasurement(0)->constructPlane(seedFromTrack);
      }

      assert (plane.get() != NULL);



      // do extrapolation and set reference state infos
      boost::scoped_ptr<StateOnPlane> stateToExtrapolate(NULL);
      if (prevFitterInfo == NULL) { // first measurement
        if (debugLvl_ > 0)
          std::cout << "prevFitterInfo == NULL \n";
        if (smoothedState != NULL) {
          if (debugLvl_ > 0)
            std::cout << "extrapolate smoothedState to plane\n";
          stateToExtrapolate.reset(new StateOnPlane(*smoothedState));
        }
        else if (referenceState != NULL) {
          if (debugLvl_ > 0)
            std::cout << "extrapolate referenceState to plane\n";
          stateToExtrapolate.reset(new StateOnPlane(*referenceState));
        }
        else if (rep != tr->getCardinalRep() &&
                  trackPoint->hasFitterInfo(tr->getCardinalRep()) &&
                  dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep())) != NULL &&
                  static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->hasPredictionsAndUpdates() ) {
          if (debugLvl_ > 0)
            std::cout << "extrapolate smoothed state of cardinal rep fit to plane\n";
          TVector3 pos, mom;
          const MeasuredStateOnPlane& fittedState = static_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(tr->getCardinalRep()))->getFittedState(true);
          tr->getCardinalRep()->getPosMom(fittedState, pos, mom);
          stateToExtrapolate.reset(new StateOnPlane(rep));
          rep->setPosMom(*stateToExtrapolate, pos, mom);
          rep->setQop(*stateToExtrapolate, tr->getCardinalRep()->getQop(fittedState));
        }
        else {
          if (debugLvl_ > 0)
            std::cout << "extrapolate seed from track to plane\n";
          stateToExtrapolate.reset(new StateOnPlane(rep));
          rep->setPosMom(*stateToExtrapolate, tr->getStateSeed());
        }
      } // end if (prevFitterInfo == NULL)
      else {
        if (prevSmoothedState != NULL) {
          if (debugLvl_ > 0)
            std::cout << "extrapolate prevSmoothedState to plane \n";
          stateToExtrapolate.reset(new StateOnPlane(*prevSmoothedState));
        }
        else {
          assert (prevReferenceState != NULL);
          if (debugLvl_ > 0)
            std::cout << "extrapolate prevReferenceState to plane \n";
          stateToExtrapolate.reset(new StateOnPlane(*prevReferenceState));
        }
      }

      // make sure track is consistent if extrapolation fails
      if (prevReferenceState != NULL)
        prevReferenceState->resetBackward();
      fitterInfo->deleteReferenceInfo();

      if (prevFitterInfo != NULL) {
        rep->extrapolateToPlane(*stateToExtrapolate, prevFitterInfo->getPlane());
        if (debugLvl_ > 0)
          std::cout << "extrapolated stateToExtrapolate to plane of prevFitterInfo (plane could have changed!) \n";
      }

      double segmentLen = rep->extrapolateToPlane(*stateToExtrapolate, plane, false, true);
      trackLen += segmentLen;
      if (debugLvl_ > 0) {
        std::cout << "extrapolated stateToExtrapolate by " << segmentLen << " cm.\t";
        std::cout << "charge of stateToExtrapolate: " << rep->getCharge(*stateToExtrapolate) << " \n";
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
      if (prevReferenceState != NULL) {
        prevReferenceState->setBackwardSegmentLength(-segmentLen);
        prevReferenceState->setBackwardTransportMatrix(BTransportMatrix_);
        prevReferenceState->setBackwardNoiseMatrix(BNoiseMatrix_);
        prevReferenceState->setBackwardDeltaState(backwardDeltaState_);
      }


      // create new reference state
      newRefState = true;
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
        assert((*measurement) != NULL);
        fitterInfo->addMeasurementsOnPlane((*measurement)->constructMeasurementsOnPlane(rep, plane));
      }
      if (oldWeights.size() == fitterInfo->getNumMeasurements()) {
        fitterInfo->setWeights(oldWeights);
        fitterInfo->fixWeights(oldWeightsFixed);
      }

      changedSmthg = true;

      prevReferenceState = referenceState;
      prevFitterInfo = fitterInfo;
      prevSmoothedState = smoothedState;

    } // end loop over TrackPoints

  }
  catch (Exception& e) {

    if (debugLvl_ > 0) {
      std::cout << "exception at hit " << i << "\n";
      std::cerr << e.what();
    }

    // clean up
    removeForwardBackwardInfo(tr, rep, notChangedUntil, notChangedFrom);

    // set sorting parameters of TrackPoints where no reference state could be calculated
    for (; i<nPoints; ++i) {
      TrackPoint* trackPoint = tr->getPoint(i);

      if (setSortingParams)
        trackPoint->setSortingParameter(trackLen);

      trackPoint->deleteFitterInfo(rep);
    }

    //prevReferenceState->resetForward();
    //referenceState->resetBackward();

    //Exception exc("KalmanFitterRefTrack::prepareTrack: got an exception.",__LINE__,__FILE__);
    //exc.setFatal();
    //throw exc;

    return true;
  }




  removeForwardBackwardInfo(tr, rep, notChangedUntil, notChangedFrom);

  if (firstBackwardUpdate) {
    KalmanFitterInfo* fi = static_cast<KalmanFitterInfo*>(tr->getPointWithMeasurement(0)->getFitterInfo(rep));
    if (! fi->hasForwardPrediction()) {
      if (debugLvl_ > 0)
        std::cout << "set backwards update of first point as forward prediction (with blown up cov) \n";
      if (fi->getPlane() != firstBackwardUpdate->getPlane()) {
        rep->extrapolateToPlane(*firstBackwardUpdate, fi->getPlane());
      }
      firstBackwardUpdate->blowUpCov(blowUpFactor_);
      fi->setForwardPrediction(new MeasuredStateOnPlane(*firstBackwardUpdate));
    }
  }

  KalmanFitStatus* fitStatus = dynamic_cast<KalmanFitStatus*>(tr->getFitStatus(rep));
  if (fitStatus != NULL)
    fitStatus->setTrackLen(trackLen);

  if (debugLvl_ > 0)
    std::cout << "trackLen of reference track = " << trackLen << "\n";

  // self check
  //assert(tr->checkConsistency());
  assert(isTrackPrepared(tr, rep));

  return changedSmthg;
}


bool
KalmanFitterRefTrack::removeOutdated(Track* tr, const AbsTrackRep* rep, int& notChangedUntil, int& notChangedFrom) {

  if (debugLvl_ > 0)
    std::cout << "KalmanFitterRefTrack::removeOutdated \n";

  bool changedSmthg(false);

  unsigned int nPoints = tr->getNumPoints();
  notChangedUntil = nPoints-1;
  notChangedFrom = 0;

  // loop over TrackPoints
  for (unsigned int i=0; i<nPoints; ++i) {

    TrackPoint* trackPoint = tr->getPoint(i);

    // check if we have a measurement
    if (!trackPoint->hasRawMeasurements()) {
      if (debugLvl_ > 0)
        std::cout << "TrackPoint has no rawMeasurements -> continue \n";
      continue;
    }

    // get fitterInfo
    KalmanFitterInfo* fitterInfo(NULL);
    if (trackPoint->hasFitterInfo(rep))
      fitterInfo = dynamic_cast<KalmanFitterInfo*>(trackPoint->getFitterInfo(rep));

    if (fitterInfo == NULL)
      continue;


    // check if we need to calculate or update reference state
    if (fitterInfo->hasReferenceState()) {

      if (! fitterInfo->hasPredictionsAndUpdates()) {
        if (debugLvl_ > 0)
          std::cout << "reference state but not all predictions & updates -> do not touch reference state. \n";
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
        if (debugLvl_ > 0)
          std::cout << "reference state is near smoothed state ->  do not update reference state, chi2 = " << chi2 << "\n";
        continue;
      } else {
        if (debugLvl_ > 0)
          std::cout << "reference state is not close to smoothed state, chi2 = " << chi2 << "\n";
      }
    }

    if (debugLvl_ > 0)
      std::cout << "remove reference info \n";

    fitterInfo->deleteReferenceInfo();
    changedSmthg = true;

    // decided to update reference state -> set notChangedUntil (only once)
    if (notChangedUntil == (int)nPoints-1)
      notChangedUntil = i-1;

    notChangedFrom = i+1;

  } // end loop over TrackPoints


  if (debugLvl_ > 0)
    tr->Print("C");

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
  if (debugLvl_ > 0)
    std::cout << " KalmanFitterRefTrack::processTrackPoint " << fi->getTrackPoint() << "\n";

  unsigned int dim = fi->getRep()->getDim();

  p_.Zero(); // p_{k|k-1}
  C_.Zero(); // C_{k|k-1}

  // predict
  if (prevFi != NULL) {
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
      std::cout << "\033[31m";
      std::cout << "F (Transport Matrix) "; F.Print();
      std::cout << "p_{k,r} (reference state) "; fi->getReferenceState()->getState().Print();
      std::cout << "c (delta state) "; fi->getReferenceState()->getDeltaState(direction).Print();
      std::cout << "F*p_{k-1,r} + c "; (F *prevFi->getReferenceState()->getState() + fi->getReferenceState()->getDeltaState(direction)).Print();
    }
  }
  else {
    if (fi->hasPrediction(direction)) {
      if (debugLvl_ > 0)
        std::cout << "  Use prediction as start \n";
      p_ = fi->getPrediction(direction)->getState();
      C_ = fi->getPrediction(direction)->getCov();
    }
    else {
      if (debugLvl_ > 0)
        std::cout << "  Use reference state and seed cov as start \n";
      const AbsTrackRep *rep = fi->getReferenceState()->getRep();
      p_ = fi->getReferenceState()->getState();

      // Convert from 6D covariance of the seed to whatever the trackRep wants.
      TMatrixDSym dummy(p_.GetNrows());
      MeasuredStateOnPlane mop(p_, dummy, fi->getReferenceState()->getPlane(), rep, fi->getReferenceState()->getAuxInfo());
      TVector3 pos, mom;
      rep->getPosMom(mop, pos, mom);
      rep->setPosMomCov(mop, pos, mom, fi->getTrackPoint()->getTrack()->getCovSeed());
      // Blow up, set.
      mop.blowUpCov(blowUpFactor_);
      fi->setPrediction(new MeasuredStateOnPlane(mop), direction);
      C_ = mop.getCov();
    }
    if (debugLvl_ > 1) {
      std::cout << "\033[31m";
      std::cout << "p_{k,r} (reference state)"; fi->getReferenceState()->getState().Print();
    }
  }

  if (debugLvl_ > 1) {
    std::cout << " p_{k|k-1} (state prediction)"; p_.Print();
    std::cout << " C_{k|k-1} (covariance prediction)"; C_.Print();
    std::cout << "\033[0m";
  }

  // update(s)
  double chi2inc = 0;
  double ndfInc = 0;
  const std::vector<MeasurementOnPlane *> measurements = getMeasurements(fi, tp, direction);
  for (std::vector<MeasurementOnPlane *>::const_iterator it = measurements.begin(); it != measurements.end(); ++it) {
    const MeasurementOnPlane& m = **it;

    if (!canIgnoreWeights() && m.getWeight() <= 1.01E-10) {
      if (debugLvl_ > 1) {
        std::cout << "Weight of measurement is almost 0, continue ... /n";
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
      std::cout << "\033[34m";
      std::cout << "m (measurement) "; m.getState().Print();
      std::cout << "V ((weighted) measurement covariance) "; (1./m.getWeight() * m.getCov()).Print();
      std::cout << "residual        "; res_.Print();
      std::cout << "\033[0m";
    }
    p_ +=  TMatrixD(CHt, TMatrixD::kMult, covSumInv_) * res_; // updated state
    if (debugLvl_ > 1) {
      std::cout << "\033[32m";
      std::cout << " update"; (TMatrixD(CHt, TMatrixD::kMult, covSumInv_) * res_).Print();
      std::cout << "\033[0m";
    }

    covSumInv_.Similarity(CHt); // with (C H^T)^T = H C^T = H C  (C is symmetric)
    C_ -= covSumInv_; // updated Cov

    if (debugLvl_ > 1) {
      //std::cout << " C update "; covSumInv_.Print();
      std::cout << "\033[32m";
      std::cout << " p_{k|k} (updated state)"; p_.Print();
      std::cout << " C_{k|k} (updated covariance)"; C_.Print();
      std::cout << "\033[0m";
    }

    // Calculate chi² increment.  At the first point chi2inc == 0 and
    // the matrix will not be invertible.
    res_ = m.getState();
    res_ -= H->Hv(p_); // new residual
    if (debugLvl_ > 1) {
      std::cout << " resNew ";
      res_.Print();
    }

    // only calculate chi2inc if res != NULL.
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
        if (debugLvl_ > 1)
          std::cerr << e.what();
        couldInvert = false;
      }

      if (couldInvert) {
        if (debugLvl_ > 1) {
          std::cout << " Rinv ";
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
    std::cout << " chi² inc " << chi2inc << "\t";
    std::cout << " ndf inc  " << ndfInc << "\t";
    std::cout << " charge of update  " << fi->getRep()->getCharge(*upState) << "\n";
  }

  // check
  assert(fi->checkConsistency());

}
