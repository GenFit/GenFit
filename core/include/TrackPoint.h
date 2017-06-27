/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

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

#ifndef genfit_TrackPoint_h
#define genfit_TrackPoint_h

#include "AbsMeasurement.h"
#include "AbsFitterInfo.h"
#include "ThinScatterer.h"

#include <TObject.h>

#include <map>
#include <vector>
#include <memory>


namespace genfit {

class Track;
class KalmanFitterInfo;

/**
 * @brief Object containing AbsMeasurement and AbsFitterInfo objects.
 *
 */
class TrackPoint : public TObject {

 public:

  TrackPoint();
  TrackPoint(Track* track);

  /**
   * @brief Contructor taking list of measurements.
   *
   * AbsMeasurement::setTrackPoint() of each measurement will be called.
   * TrackPoint takes ownership over rawMeasurements.
   */
  TrackPoint(const std::vector< genfit::AbsMeasurement* >& rawMeasurements, Track* track);

  /**
   * @brief Contructor taking one measurement.
   *
   * AbsMeasurement::setTrackPoint() of the measurement will be called.
   * TrackPoint takes ownership over the rawMeasurement.
   */
  TrackPoint(genfit::AbsMeasurement* rawMeasurement, Track* track);

  TrackPoint(const TrackPoint&); // copy constructor
  TrackPoint& operator=(TrackPoint); // assignment operator
  void swap(TrackPoint& other);

  /**
   * custom copy constructor where all TrackRep pointers are exchanged according to the map.
   * FitterInfos with a rep in repsToIgnore will NOT be copied.
   */
  TrackPoint(const TrackPoint& rhs,
      const std::map<const genfit::AbsTrackRep*, genfit::AbsTrackRep*>& map,
      const std::vector<const genfit::AbsTrackRep*> * repsToIgnore = nullptr);

  virtual ~TrackPoint();


  double getSortingParameter() const {return sortingParameter_;}

  Track* getTrack() const {return track_;}
  void setTrack(Track* track) {track_ = track;}

  const std::vector< genfit::AbsMeasurement* >& getRawMeasurements() const {return rawMeasurements_;}
  AbsMeasurement* getRawMeasurement(int i = 0) const;
  unsigned int getNumRawMeasurements() const {return rawMeasurements_.size();}
  bool hasRawMeasurements() const {return (rawMeasurements_.size() != 0);}
  //! Get list of all fitterInfos
  std::vector< genfit::AbsFitterInfo* > getFitterInfos() const;
  //! Get fitterInfo for rep. Per default, use cardinal rep
  AbsFitterInfo* getFitterInfo(const AbsTrackRep* rep = nullptr) const;
  //! Helper to avoid casting
  KalmanFitterInfo* getKalmanFitterInfo(const AbsTrackRep* rep = nullptr) const;
  bool hasFitterInfo(const AbsTrackRep* rep) const {
    return (fitterInfos_.find(rep) != fitterInfos_.end());
  }

  ThinScatterer* getMaterialInfo() const {return thinScatterer_.get();}
  bool hasThinScatterer() const {return thinScatterer_.get() != nullptr;}


  void setSortingParameter(double sortingParameter) {sortingParameter_ = sortingParameter;}
  //! Takes ownership and sets this as measurement's trackPoint
  void addRawMeasurement(genfit::AbsMeasurement* rawMeasurement) {assert(rawMeasurement!=nullptr); rawMeasurement->setTrackPoint(this); rawMeasurements_.push_back(rawMeasurement);}
  void deleteRawMeasurements();
  //! Takes Ownership
  void setFitterInfo(genfit::AbsFitterInfo* fitterInfo);
  void deleteFitterInfo(const AbsTrackRep* rep) {delete fitterInfos_[rep]; fitterInfos_.erase(rep);}

  void setScatterer(ThinScatterer* scatterer) {thinScatterer_.reset(scatterer);}

  void Print(const Option_t* = "") const;

  /**
   * This function is used when reading the TrackPoint and is called
   * by the owner in order to build fitterInfos_ from vFitterInfos_.
   * This requires that the track_ be set.  It also empties
   * vFitterInfos_ which has served its purpose after this function is
   * called.
   */
  void fixupRepsForReading();

 private:
  double sortingParameter_;

  //! Pointer to Track where TrackPoint belongs to
  Track* track_; //! No ownership

  //! Can be more than one, e.g. multiple measurements in the same Si detector, left and right measurements of a wire detector etc.
  std::vector<AbsMeasurement*> rawMeasurements_; // Ownership

  std::map< const AbsTrackRep*, AbsFitterInfo* > fitterInfos_; //! Ownership over FitterInfos

  /**
   * The following map is read while streaming.  After reading the
   * TrackPoint, the Track's streamer will call fixupRepsForReading,
   * and this map will be translated into the map fitterInfos. The
   * map is indexed by the ids of the corresponding TrackReps.
   */
  std::map<unsigned int, AbsFitterInfo*> vFitterInfos_; //!

  std::unique_ptr<ThinScatterer> thinScatterer_; // Ownership

 public:

  ClassDef(TrackPoint,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_TrackPoint_h
