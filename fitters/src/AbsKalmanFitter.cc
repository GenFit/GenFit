/* Copyright 2013, Ludwig-Maximilians Universität München,
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

/* This implements the simple Kalman fitter with no reference track
   that uses the stateSeed only until it forgets about it after the
   first few hits.  */

#include "Track.h"
#include "TrackPoint.h"
#include "Exception.h"
#include "KalmanFitterInfo.h"
#include "AbsMeasurement.h"

#include "AbsKalmanFitter.h"
#include <Math/ProbFunc.h>


//#define DEBUG

using namespace genfit;



void AbsKalmanFitter::getChiSquNdf(const Track* tr, const AbsTrackRep* rep,
    double& bChi2, double& fChi2,
    double& bNdf,  double& fNdf) const {

  bChi2 = 0;
  fChi2 = 0;
  bNdf = -1. * rep->getDim();
  fNdf = -1. * rep->getDim();

  const std::vector<TrackPoint*>& pointsWM = tr->getPointsWithMeasurement();
  for (std::vector<TrackPoint*>::const_iterator tpIter = pointsWM.begin(), endIter = pointsWM.end(); tpIter != endIter; ++tpIter) {
    if (! (*tpIter)->hasFitterInfo(rep))
      continue;

    AbsFitterInfo* afi = (*tpIter)->getFitterInfo(rep);
    KalmanFitterInfo* fi = dynamic_cast<KalmanFitterInfo*>(afi);
    if (!fi) {
      Exception exc("AbsKalmanFitter::getChiSqu(): fitterInfo is not a KalmanFitterInfo", __LINE__,__FILE__);
      throw exc;
    }

    KalmanFittedStateOnPlane* fup = fi->getForwardUpdate();
    KalmanFittedStateOnPlane* bup = fi->getBackwardUpdate();

    if (fup == NULL || bup == NULL) {
      Exception exc("AbsKalmanFitter::getChiSqu(): fup == NULL || bup == NULL", __LINE__,__FILE__);
      throw exc;
    }

    bChi2 += bup->getChiSquareIncrement();
    fChi2 += fup->getChiSquareIncrement();

    bNdf += bup->getNdf();
    fNdf += fup->getNdf();
  }

  if (bNdf < 0)
    bNdf = 0;

  if (fNdf < 0)
    fNdf = 0;
}


double AbsKalmanFitter::getChiSqu(const Track* tr, const AbsTrackRep* rep, int direction) const {
  double bChi2(0), fChi2(0), bNdf(0), fNdf(0);

  getChiSquNdf(tr, rep, bChi2, fChi2, bNdf, fNdf);

  if (direction < 0)
    return bChi2;
  return fChi2;
}

double AbsKalmanFitter::getNdf(const Track* tr, const AbsTrackRep* rep, int direction) const {
  double bChi2(0), fChi2(0), bNdf(0), fNdf(0);

  getChiSquNdf(tr, rep, bChi2, fChi2, bNdf, fNdf);

  if (direction < 0)
    return bNdf;
  return fNdf;
}

double AbsKalmanFitter::getRedChiSqu(const Track* tr, const AbsTrackRep* rep, int direction) const {
  double bChi2(0), fChi2(0), bNdf(0), fNdf(0);

  getChiSquNdf(tr, rep, bChi2, fChi2, bNdf, fNdf);

  if (direction < 0)
    return bChi2/bNdf;
  return fChi2/fNdf;
}

double AbsKalmanFitter::getPVal(const Track* tr, const AbsTrackRep* rep, int direction) const {
  double bChi2(0), fChi2(0), bNdf(0), fNdf(0);

  getChiSquNdf(tr, rep, bChi2, fChi2, bNdf, fNdf);

  if (direction < 0)
    return std::max(0.,ROOT::Math::chisquared_cdf_c(bChi2, bNdf));
  return std::max(0.,ROOT::Math::chisquared_cdf_c(fChi2, fNdf));
}


bool AbsKalmanFitter::isTrackPrepared(const Track* tr, const AbsTrackRep* rep) const {
  const std::vector<TrackPoint*>& points = tr->getPointsWithMeasurement();

  if (points.size() == 0)
    return true;

  for (std::vector<TrackPoint*>::const_iterator pIt = points.begin(), pEnd = points.end(); pIt != pEnd; ++pIt) {
    KalmanFitterInfo* fi = dynamic_cast<KalmanFitterInfo*>((*pIt)->getFitterInfo(rep));

    if (!fi)
      continue;

    if (!(fi->checkConsistency()))
      return false;

    if (!(fi->hasReferenceState()))
      return false;
  }

  return true;
}

bool AbsKalmanFitter::isTrackFitted(const Track* tr, const AbsTrackRep* rep) const {
  if (! tr->getFitStatus(rep)->isFitted())
    return false;

  // in depth check

  const std::vector< TrackPoint* >& points = tr->getPointsWithMeasurement();

  if (points.size() == 0)
    return true;

  for (std::vector<TrackPoint*>::const_iterator pIt = points.begin(), pEnd = points.end(); pIt != pEnd; ++pIt) {
    KalmanFitterInfo* fi = dynamic_cast<KalmanFitterInfo*>((*pIt)->getFitterInfo(rep));
    if (!fi)
      return false;

    if (!(fi->checkConsistency()))
      return false;

    if (!(fi->hasForwardUpdate()))
      return false;

    if (!(fi->hasBackwardUpdate()))
      return false;
  }

  return true;
}


const std::vector<MeasurementOnPlane *> AbsKalmanFitter::getMeasurements(const KalmanFitterInfo* fi, const TrackPoint* tp, int direction) const {

  switch (multipleMeasurementHandling_) {
    case weightedAverage :
    case unweightedAverage :
      return fi->getMeasurementsOnPlane();

    case weightedClosestToReference :
    case unweightedClosestToReference :
    {
      if (!fi->hasReferenceState()) {
        Exception e("AbsKalmanFitter::getMeasurement: no ReferenceState.", __LINE__,__FILE__);
        e.setFatal();
        throw e;
      }
      std::vector<MeasurementOnPlane *> retVal;
      retVal.push_back(fi->getClosestMeasurementOnPlane(fi->getReferenceState()));
      return retVal;
    }

    case weightedClosestToPrediction :
    case unweightedClosestToPrediction :
    {
      if (!fi->hasPrediction(direction)) {
        Exception e("AbsKalmanFitter::getMeasurement: no prediction.", __LINE__,__FILE__);
        e.setFatal();
        throw e;
      }
      std::vector<MeasurementOnPlane *> retVal;
      retVal.push_back(fi->getClosestMeasurementOnPlane(fi->getPrediction(direction)));
      return retVal;
    }


    case weightedClosestToReferenceWire :
    case unweightedClosestToReferenceWire :
    {
      if (tp->getNumRawMeasurements() == 1 && tp->getRawMeasurement()->isLeftRightMeasurement()) {
        if (!fi->hasReferenceState()) {
          Exception e("AbsKalmanFitter::getMeasurement: no ReferenceState.", __LINE__,__FILE__);
          e.setFatal();
          throw e;
        }
        std::vector<MeasurementOnPlane *> retVal;
        retVal.push_back(fi->getClosestMeasurementOnPlane(fi->getReferenceState()));
        return retVal;
      }
      else
        return fi->getMeasurementsOnPlane();
    }

    case weightedClosestToPredictionWire :
    case unweightedClosestToPredictionWire :
    {
      if (tp->getNumRawMeasurements() == 1 && tp->getRawMeasurement()->isLeftRightMeasurement()) {
        if (!fi->hasPrediction(direction)) {
          Exception e("AbsKalmanFitter::getMeasurement: no prediction.", __LINE__,__FILE__);
          e.setFatal();
          throw e;
        }
        std::vector<MeasurementOnPlane *> retVal;
        retVal.push_back(fi->getClosestMeasurementOnPlane(fi->getPrediction(direction)));
        return retVal;
      }
      else
        return fi->getMeasurementsOnPlane();
    }


    default:
    {
      Exception e("AbsKalmanFitter::getMeasurement: choice not valid.", __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }
  }
}


bool AbsKalmanFitter::canIgnoreWeights() const {
  switch (multipleMeasurementHandling_) {
    case unweightedAverage :
    case unweightedClosestToReference :
    case unweightedClosestToPrediction :
    case unweightedClosestToReferenceWire :
    case unweightedClosestToPredictionWire :
      return true;

    case weightedAverage :
    case weightedClosestToReference :
    case weightedClosestToPrediction :
    case weightedClosestToReferenceWire :
    case weightedClosestToPredictionWire :
      return false;

    default:
    {
      Exception e("AbsKalmanFitter::canIgnoreWeights: choice not valid.", __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }
  }
}
