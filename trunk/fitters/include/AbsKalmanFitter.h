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
/** @addtogroup genfit
 * @{
 */

#ifndef genfit_AbsKalmanFitter_h
#define genfit_AbsKalmanFitter_h

#include "AbsFitter.h"
#include "MeasurementOnPlane.h"
#include "TrackPoint.h"


namespace genfit {

class KalmanFitterInfo;

enum eMultipleMeasurementHandling {
  weightedAverage, /**<  weighted average between measurements; used by DAF */
  unweightedAverage, /**<  average between measurements, all weighted with 1 */
  weightedClosestToReference, /**<  closest to reference, weighted with its weight_ */
  unweightedClosestToReference, /**<  closest to reference, weighted with 1 */
  weightedClosestToPrediction, /**<  closest to prediction, weighted with its weight_ */
  unweightedClosestToPrediction, /**<  closest to prediction, weighted with 1 */
  weightedClosestToReferenceWire, /**<  if corresponding TrackPoint has one WireMeasurement, select closest to reference, weighted with its weight_. Otherwise use weightedAverage. */
  unweightedClosestToReferenceWire, /**<  if corresponding TrackPoint has one WireMeasurement, select closest to reference, weighted with 1. Otherwise use unweightedAverage. */
  weightedClosestToPredictionWire, /**<  if corresponding TrackPoint has one WireMeasurement, select closest to prediction, weighted with its weight_. Otherwise use weightedAverage. */
  unweightedClosestToPredictionWire /**<  if corresponding TrackPoint has one WireMeasurement, select closest to prediction, weighted with 1. Otherwise use unweightedAverage. Recommended for KalmanFitter to 'resolve' l/r ambiguities */
};

/**
 * @brief Abstract base class for Kalman fitter and derived fitting algorithms
 */
class AbsKalmanFitter : public AbsFitter {

 public:

  AbsKalmanFitter(unsigned int maxIterations = 4, double deltaPval = 1e-3, double blowUpFactor = 1e3)
    : AbsFitter(), minIterations_(2), maxIterations_(maxIterations), deltaPval_(deltaPval), relChi2Change_(0.2),
      blowUpFactor_(blowUpFactor), resetOffDiagonals_(true), blowUpMaxVal_(1.E6),
      multipleMeasurementHandling_(unweightedClosestToPredictionWire),
      maxFailedHits_(-1) {
    if (minIterations_ > maxIterations_)
      minIterations_ = maxIterations_;
  }

  virtual ~AbsKalmanFitter() {;}

  //virtual void fitTrack(Track* tr, const AbsTrackRep* rep, double& chi2, double& ndf, int direction) = 0;

  void getChiSquNdf(const Track* tr, const AbsTrackRep* rep, double& bChi2, double& fChi2, double& bNdf,  double& fNdf) const;
  double getChiSqu(const Track* tr, const AbsTrackRep* rep, int direction = -1) const;
  double getNdf(const Track* tr, const AbsTrackRep* rep, int direction = -1) const;
  double getRedChiSqu(const Track* tr, const AbsTrackRep* rep, int direction = -1) const;
  double getPVal(const Track* tr, const AbsTrackRep* rep, int direction = -1) const;

  unsigned int getMinIterations() const {return minIterations_;}
  unsigned int getMaxIterations() const {return maxIterations_;}
  double getDeltaPval() const {return deltaPval_;}
  double getRelChi2Change() const {return relChi2Change_;}
  double getBlowUpFactor() const {return blowUpFactor_;}
  bool getResetOffDiagonals() const {return resetOffDiagonals_;}
  double getBlowUpMaxVal() const {return blowUpMaxVal_;}
  eMultipleMeasurementHandling getMultipleMeasurementHandling() const {return multipleMeasurementHandling_;}
  int getMaxFailedHits() const {return maxFailedHits_;}

  //! Set the minimum number of iterations
  virtual void setMinIterations(unsigned int n) {minIterations_ = std::max((unsigned int)1,n); if (maxIterations_ < minIterations_) maxIterations_ = minIterations_;}
  //! Set the maximum number of iterations
  virtual void setMaxIterations(unsigned int n) {maxIterations_ = n; if (minIterations_ > maxIterations_) minIterations_ = maxIterations_;}

  /**
   * @brief Set Convergence criterion
   *
   * if track total P-value changes less than this between consecutive iterations, consider the track converged.
   * chi² from the backwards fit is used.
   */
  void setDeltaPval(double deltaPval) {deltaPval_ = deltaPval;}

  /**
   * @ brief Set Non-convergence criterion
   *
   * if the relative chi^2 between two iterations is larger than relChi2Change_, the fit is NOT converged and
   * further iterations will be done, even if the deltaPval_ convergence criterium is met.
   * This is especially useful for fits which have a p-value of almost 0 (possibly due to bad start values),
   * where the p-value from one iteration to the next might not change much. However, a significant change in
   * chi^2 tells us, that the fit might not yet be converged.
   */
  void setRelChi2Change(double relChi2Change) {relChi2Change_ = relChi2Change;}

  void setBlowUpFactor(double blowUpFactor) {blowUpFactor_ = blowUpFactor;}
  void setResetOffDiagonals(bool resetOffDiagonals) {resetOffDiagonals_ = resetOffDiagonals;}
  //! Limit the cov entries to this maximum value when blowing up the cov. Set to negative value to disable. Default is 1.E6.
  /**
   * This is especially useful for fits where the measurements do not constrain one direction,
   * e.g. parallel wire measurements. The fit will not be constrained along the wire direction.
   * This also means that the covariance along the wire direction will not get smaller during the fit.
   * However, it will be blown up before the next iteration, leading to an exponential growth of
   * the covariance element(s) corresponding to the wire direction. This can then lead to numerical problems.
   * To prevent this, the maximum value of the covariance elements after blowing up can be limited.
   */
  void setBlowUpMaxVal(double blowUpMaxVal) {blowUpMaxVal_= blowUpMaxVal;}

  //! How should multiple measurements be handled?
  void setMultipleMeasurementHandling(eMultipleMeasurementHandling mmh) {multipleMeasurementHandling_ = mmh;}

  virtual void setMaxFailedHits(int val) {maxFailedHits_ = val;}

  bool isTrackPrepared(const Track* tr, const AbsTrackRep* rep) const;
  bool isTrackFitted(const Track* tr, const AbsTrackRep* rep) const;

  //! returns if the fitter can ignore the weights and handle the MeasurementOnPlanes as if they had weight 1.
  bool canIgnoreWeights() const;

 protected:

  //! get the measurementsOnPlane taking the multipleMeasurementHandling_ into account
  const std::vector<MeasurementOnPlane *> getMeasurements(const KalmanFitterInfo* fi, const TrackPoint* tp, int direction) const;

  //! Minimum number of iterations to attempt.  Forward and backward are counted as one iteration.
  unsigned int minIterations_;

  //! Maximum number of iterations to attempt.  Forward and backward are counted as one iteration.
  unsigned int maxIterations_;
  /**
   * @brief Convergence criterion
   *
   * if track total P-value changes less than this between consecutive iterations, consider the track converged.
   * chi² from the backwards fit is used.
   */
  double deltaPval_;
  /**
   * @ brief Non-convergence criterion
   *
   * if the relative chi^2 between two iterations is larger than relChi2Change_, the fit is NOT converged and
   * further iterations will be done, even if the deltaPval_ convergence criterium is met.
   * This is especially useful for fits which have a p-value of almost 0 (possibly due to bad start values),
   * where the p-value from one iteration to the next might not change much. However, a significant change in
   * chi^2 tells us, that the fit might not yet be converged.
   */
  double relChi2Change_;

  //! Blow up the covariance of the forward (backward) fit by this factor before seeding the backward (forward) fit.
  double blowUpFactor_;
  //! Reset the off-diagonals to 0 when blowing up the cov
  bool resetOffDiagonals_;
  //! Limit the cov entries to this maxuimum value when blowing up the cov.
  double blowUpMaxVal_;

  //! How to handle if there are multiple MeasurementsOnPlane
  eMultipleMeasurementHandling multipleMeasurementHandling_;

  //! after how many failed hits (exception during construction of plane, extrapolation etc.) should the fit be cancelled.
  //! -1 means don't cancel
  int maxFailedHits_;

 public:

  ClassDef(AbsKalmanFitter, 2)

};

} /* End of namespace genfit */
/** @} */

#endif //genfit_AbsKalmanFitter_h
