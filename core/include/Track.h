/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch & Tobias Schlüter

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

#ifndef genfit_Track_h
#define genfit_Track_h

#include "AbsTrackRep.h"
#include "FitStatus.h"
#include "MeasurementFactory.h"
#include "TrackCand.h"
#include "TrackPoint.h"

#include <vector>
#include <TObject.h>
#include <TVectorD.h>


namespace genfit {

class KalmanFitStatus;

/**
 * @brief Helper class for TrackPoint sorting, used in Track::sort().
 */
class TrackPointComparator {
 public:
  /**
   * Comparison operator used in Track::sort(). Compares sorting parameter.
   */
  bool operator() (const TrackPoint* lhs, const TrackPoint* rhs) const {
    return lhs->getSortingParameter() < rhs->getSortingParameter();
  }
};


/**
 * @brief Collection of TrackPoint objects, AbsTrackRep objects and FitStatus objects
 *
 *  Holds a number of AbsTrackRep objects, which correspond to the
 *  different particle hypotheses or track models which should be fitted.
 *  A 6D seed #stateSeed_ (x,y,z,p_x,p_y,p_z) and 6x6 #covSeed_ should be provided as start values for fitting.
 *  When fitting the Track with a AbsFitter,
 *  a FitStatus object will be created, containing information about the fit.
 *  The fitted states will be stored in AbsFitterInfo objects in every TrackPoints.
 *
 *  The fit will be performed for every AbsTrackRep,
 *  so after the fit there will be one AbsFitterInfo for each AbsTrackRep
 *  in every TrackPoint, as well as one FitStatus for every AbsTrackRep.
 *
 */
class Track : public TObject {

 public:

  Track();

  /**
   * @ brief Construct Track from TrackCand, using a MeasurementFactory
   *
   * The MeasurementFactory will be used to create AbsMeasuremen objects.
   * TrackPoints will be created.
   * If two or more consecutive PlanarMeasurement objects with the same detector- and planeId
   * are created by the factory, they will be put into the same TrackPoint.
   *
   * Optionally, a AbsTrackRep can be provided.
   *
   * The stateSeed_ of the Track will be filled with the seed of the TrackCand.
   * A guess for covSeed_ will be made using the largest entry of the cov of the first measurement
   * and the number of measurements (For the covSeed_, it is just important that it will be
   * big enough not to bias the fit too much, but not too big in order to avoid
   * numerical problems).
   */
  Track(const TrackCand& trackCand, const MeasurementFactory<genfit::AbsMeasurement>& factory, AbsTrackRep* rep = nullptr);

  Track(AbsTrackRep* trackRep, const TVectorD& stateSeed);
  Track(AbsTrackRep* trackRep, const TVector3& posSeed, const TVector3& momSeed);
  Track(AbsTrackRep* trackRep, const TVectorD& stateSeed, const TMatrixDSym& covSeed);

  Track(const Track&); // copy constructor
  Track& operator=(Track); // assignment operator
  void swap(Track& other); // nothrow

  virtual ~Track();
  virtual void Clear(Option_t* = "");

  void createMeasurements(const TrackCand& trackCand, const MeasurementFactory<genfit::AbsMeasurement>& factory);

  TrackPoint* getPoint(int id) const;
  const std::vector< genfit::TrackPoint* > & getPoints() const {return trackPoints_;}
  unsigned int getNumPoints() const {return trackPoints_.size();}

  TrackPoint* getPointWithMeasurement(int id) const;
  const std::vector< genfit::TrackPoint* > & getPointsWithMeasurement() const  {return trackPointsWithMeasurement_;}
  unsigned int getNumPointsWithMeasurement() const {return trackPointsWithMeasurement_.size();}

  TrackPoint* getPointWithMeasurementAndFitterInfo(int id, const AbsTrackRep* rep = nullptr) const;
  TrackPoint* getPointWithFitterInfo(int id, const AbsTrackRep* rep = nullptr) const;

  /**
   * @brief Shortcut to get FittedStates.
   *
   * Uses getPointWithFitterInfo(id, rep).
   * Gets the fitted state at trackpoint id for the track representation rep.
   * Per default, the fitted state of the fitterInfo of the first TrackPoint
   * with one or more AbsFitterInfo objects
   * is returned. If no AbsTrackRep is specified, the AbsFitterInfo of the cardinal rep will be used.
   */
  const MeasuredStateOnPlane& getFittedState(int id = 0, const AbsTrackRep* rep = nullptr, bool biased = true) const;

  AbsTrackRep* getTrackRep(int id) const {return trackReps_.at(id);}
  /// Return the track representations as a list of pointers.
  const std::vector<genfit::AbsTrackRep*>& getTrackReps() const {return trackReps_;}
  unsigned int getNumReps() const {return trackReps_.size();}

  //! This is used when streaming TrackPoints.
  int getIdForRep(const AbsTrackRep* rep) const;

  /** @brief Get cardinal track representation
   *
   * The user has to choose which AbsTrackRep should be considered the
   * best one after the fit. E.g. the track representation giving the
   * smallest chi2 could be chosen. By default the first in the list is returned.
   * @sa #determineCardinalRep()
   */
  AbsTrackRep* getCardinalRep() const {return trackReps_.at(cardinalRep_);}
  unsigned int getCardinalRepId() const {return cardinalRep_;}

  //! Get the MCT track id, for MC simulations - default value = -1
  int getMcTrackId() const {return mcTrackId_;}

  //! Check if track has a FitStatus for given AbsTrackRep. Per default, check for cardinal rep.
  bool hasFitStatus(const AbsTrackRep* rep = nullptr) const;
  //! Get FitStatus for a AbsTrackRep. Per default, return FitStatus for cardinalRep.
  FitStatus* getFitStatus(const AbsTrackRep* rep = nullptr) const {if (rep == nullptr) rep = getCardinalRep(); return fitStatuses_.at(rep);}

  //! Check if track has a KalmanFitStatus for given AbsTrackRep. Per default, check for cardinal rep.
  bool hasKalmanFitStatus(const AbsTrackRep* rep = nullptr) const;
  //! If FitStatus is a KalmanFitStatus, return it. Otherwise return nullptr
  KalmanFitStatus* getKalmanFitStatus(const AbsTrackRep* rep = nullptr) const;

  void setFitStatus(FitStatus* fitStatus, const AbsTrackRep* rep);

  double getTimeSeed() const {return timeSeed_;}
  void setTimeSeed(double time) {timeSeed_ = time;}

  const TVectorD& getStateSeed() const {return stateSeed_;}
  void setStateSeed(const TVectorD& s) {stateSeed_.ResizeTo(s); stateSeed_ = s;}
  void setStateSeed(const TVector3& pos, const TVector3& mom);

  const TMatrixDSym& getCovSeed() const {return covSeed_;}
  void setCovSeed(const TMatrixDSym& c) {covSeed_.ResizeTo(c); covSeed_ = c;}

  //! Set the MCT track id, for MC simulations
  void setMcTrackId(int i) {mcTrackId_ = i;}

  /**
   * @brief Insert TrackPoint BEFORE TrackPoint with position id, if id >= 0.
   *
   * Id -1 means after last TrackPoint. Id -2 means before last TrackPoint. ...
   * Also deletes backwardInfos before new point and forwardInfos after new point.
   * Also sets Track backpointer of point accordingly.
   */
  void insertPoint(TrackPoint* point, int id = -1);

  /**
   * @brief Insert TrackPoints BEFORE TrackPoint with position id, if id >= 0.
   *
   * Id -1 means after last TrackPoint. Id -2 means before last TrackPoint. ...
   * Also deletes backwardInfos before and for new points and forwardInfos after and for new points.
   * Also sets Track backpointers of points accordingly.
   */
  void insertPoints(std::vector<genfit::TrackPoint*> points, int id = -1);

  void deletePoint(int id);

  //! Creates a new TrackPoint containing the measurement, and adds it to the track
  void insertMeasurement(AbsMeasurement* measurement, int id = -1);

  //! Delete all measurement information and the track points of the track. Does not delete track representations.
  void deleteTrackPointsAndFitStatus();
  /**
   * @brief Merge two tracks.
   *
   * The TrackPoint objects of other will be cloned and inserted
   * after id (per default, they will be appended at the end).
   * The other Track will not be altered, the TrackPoint objects will be (deep) copied.
   * Only copies the TrackPoint objects, NOT the AbsTrackRep, FitStatus, seed state and other objects of the other track.
   */
  void mergeTrack(const Track* other, int id = -1);

  void addTrackRep(AbsTrackRep* trackRep);

  //! Delete a AbsTrackRep and all corresponding AbsFitterInfo objects in every TrackPoint.
  void deleteTrackRep(int id);

  void setCardinalRep(int id);
  //! See with which AbsTrackRep the track was fitted best (converged fit w/ smallest chi2) and set the cardinal rep accordingly.
  void determineCardinalRep();

  /**
   * @brief Sort TrackPoint and according to their sorting parameters.
   *
   * Returns if the order of the TrackPoint has actually changed.
   */
  bool sort();

  //! Try to set the fitted state as seed. Return if it was successfull.
  //! Adapt the sign of all TrackReps' pdg to the actual fitted charge.
  bool udpateSeed(int id = 0, AbsTrackRep* rep = nullptr, bool biased = true);

  //! Flip the ordering of the TrackPoints
  void reverseTrackPoints();

  //! Flip direction of momentum seed
  void reverseMomSeed() {
    stateSeed_(3) *= -1; stateSeed_(4) *= -1; stateSeed_(5) *= -1;
  }

  //! Switch the pdg signs of specified rep (of all reps if rep == nullptr).
  void switchPDGSigns(AbsTrackRep* rep = nullptr);

  //! Make track ready to be fitted in reverse direction
  /**
   * Flip the order of TrackPoints and the momentum direction of the seed state.
   * If possible, take the smoothed state of the last hit as new seed state.
   * Flip charge of the TrackReps.
   */
  void reverseTrack();


  void deleteForwardInfo(int startId = 0, int endId = -1, const AbsTrackRep* rep = nullptr); // delete in range [startId, endId]. If rep == nullptr, delete for ALL reps, otherwise only for rep.
  void deleteBackwardInfo(int startId = 0, int endId = -1, const AbsTrackRep* rep = nullptr); // delete in range [startId, endId]. If rep == nullptr, delete for ALL reps, otherwise only for rep.
  void deleteReferenceInfo(int startId = 0, int endId = -1, const AbsTrackRep* rep = nullptr); // delete in range [startId, endId]. If rep == nullptr, delete for ALL reps, otherwise only for rep.
  void deleteMeasurementInfo(int startId = 0, int endId = -1, const AbsTrackRep* rep = nullptr); // delete in range [startId, endId]. If rep == nullptr, delete for ALL reps, otherwise only for rep.
  void deleteFitterInfo(int startId = 0, int endId = -1, const AbsTrackRep* rep = nullptr); // delete in range [startId, endId]. If rep == nullptr, delete for ALL reps, otherwise only for rep.

  //! get TrackLength between to trackPoints (if nullptr, for cardinal rep)
  double getTrackLen(AbsTrackRep* rep = nullptr, int startId = 0, int endId = -1) const;
  //! get time of flight in ns between to trackPoints (if nullptr, for cardinal rep)
  double getTOF(AbsTrackRep* rep = nullptr, int startId = 0, int endId = -1) const;

  /**
   * Delete the fit status and all the FitStates of the TrackPoints
   * for the given hypothesis.
   * This is equal to resetting the track for the rep, so another fit
   * can start from scratch.
   * Useful if you have changed some seeds.
   */
  void deleteFittedState(const genfit::AbsTrackRep* rep); 

  //! Construct a new TrackCand containing the hit IDs of the measurements
  /**
   * The idea is hat you can get a TrackCand for storing the hit IDs after a track has been fitted.
   * His could have been reordered, added or removed, so that the original TackCand no longer
   * represents the Track correctly.
   * You might want to call determineCardinalRep() and/or udpateSeed() before.
   */
  TrackCand* constructTrackCand() const;

  //! Helper function: For all KalmanFitterInfos belonging to rep (if nullptr, for all reps),
  //! call the fixWeights() function, so that e.g. the DAF will not alter weights anymore.
  void fixWeights(AbsTrackRep* rep = nullptr, int startId = 0, int endId = -1);

  /**
   * @brief Delete unneeded information from the Track.
   *
   * Possible options: (see also PruneFlags defined in FitStatus.h)
   * C:  prune all reps except cardinalRep
   * F:  prune all points except first point (also prune referenceInfo from fitterInfos)
   * L:  prune all points except last point (also prune referenceInfo from fitterInfos)
   * FL: prune all points except first and last point (also prune referenceInfo from fitterInfos)
   * W:  prune rawMeasurements from TrackPoints
   * R:  prune referenceInfo from fitterInfos
   * M:  prune measurementInfo from fitterInfos
   * I:  if F, L, or FL is set, prune forward (backward) info of first (last) point
   * U:  if fitterInfo is a KalmanFitterInfo, prune predictions and keep updates
   */
  void prune(const Option_t* = "CFLWRMIU");

  void Print(const Option_t* = "") const;

  void checkConsistency() const;

 private:

  void trackHasChanged();

  void fillPointsWithMeasurement();

  std::vector<AbsTrackRep*> trackReps_; // Ownership
  unsigned int cardinalRep_; // THE selected rep, default = 0;

  std::vector<TrackPoint*> trackPoints_; // Ownership
  std::vector<TrackPoint*> trackPointsWithMeasurement_; //! helper

  std::map< const AbsTrackRep*, FitStatus* > fitStatuses_; // Ownership over FitStatus*

  int mcTrackId_; /**< if MC simulation, store the mc track id here */
  double timeSeed_;
  TVectorD stateSeed_; // 6D: position, momentum
  TMatrixDSym covSeed_; // 6D


 public:
  ClassDef(Track,3)
  // Class version history:
  //  ver 3: introduces timeSeed_
};

} /* End of namespace genfit */
/** @} */

#endif // genfit_Track_h
