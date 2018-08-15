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

#ifndef genfit_KalmanFitterInfo_h
#define genfit_KalmanFitterInfo_h

#include "AbsFitterInfo.h"
#include "KalmanFittedStateOnPlane.h"
#include "MeasuredStateOnPlane.h"
#include "MeasurementOnPlane.h"
#include "ReferenceStateOnPlane.h"
#include "StateOnPlane.h"

#include <vector>

#include <memory>


namespace genfit {


/**
 *  @brief Collects information needed and produced by a AbsKalmanFitter implementations and is specific to one AbsTrackRep of the Track.
 */
class KalmanFitterInfo : public AbsFitterInfo {

 public:

  KalmanFitterInfo();
  KalmanFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep);
  virtual ~KalmanFitterInfo();

  virtual KalmanFitterInfo* clone() const;

  ReferenceStateOnPlane* getReferenceState() const {return referenceState_.get();}
  MeasuredStateOnPlane* getForwardPrediction() const {return forwardPrediction_.get();}
  MeasuredStateOnPlane* getBackwardPrediction() const {return backwardPrediction_.get();}
  MeasuredStateOnPlane* getPrediction(int direction) const {if (direction >=0) return forwardPrediction_.get(); return backwardPrediction_.get();}
  KalmanFittedStateOnPlane* getForwardUpdate() const {return forwardUpdate_.get();}
  KalmanFittedStateOnPlane* getBackwardUpdate() const {return backwardUpdate_.get();}
  KalmanFittedStateOnPlane* getUpdate(int direction) const {if (direction >=0) return forwardUpdate_.get(); return backwardUpdate_.get();}
  const std::vector< genfit::MeasurementOnPlane* >& getMeasurementsOnPlane() const {return measurementsOnPlane_;}
  MeasurementOnPlane* getMeasurementOnPlane(int i = 0) const {if (i<0) i += measurementsOnPlane_.size(); return measurementsOnPlane_.at(i);}
  //! Get weighted mean of all measurements.
  //! @param ignoreWeights If set, the weights of the individual measurements will be ignored (they will be treated as if they all had weight 1)
  MeasurementOnPlane getAvgWeightedMeasurementOnPlane(bool ignoreWeights = false) const;
  //! Get measurements which is closest to state.
  MeasurementOnPlane* getClosestMeasurementOnPlane(const StateOnPlane*) const;
  unsigned int getNumMeasurements() const {return measurementsOnPlane_.size();}
  //! Get weights of measurements.
  std::vector<double> getWeights() const;
  //! Are the weights fixed?
  bool areWeightsFixed() const {return fixWeights_;}
  //! Get unbiased or biased (default) smoothed state.
  const MeasuredStateOnPlane& getFittedState(bool biased = true) const;
  //! Get unbiased (default) or biased residual from ith measurement.
  MeasurementOnPlane getResidual(unsigned int iMeasurement = 0, bool biased = false, bool onlyMeasurementErrors = true) const; // calculate residual, track and measurement errors are added if onlyMeasurementErrors is false
  double getSmoothedChi2(unsigned int iMeasurement = 0) const;

  bool hasMeasurements() const {return getNumMeasurements() > 0;}
  bool hasReferenceState() const {return (referenceState_.get() != nullptr);}
  bool hasForwardPrediction() const {return (forwardPrediction_.get()  != nullptr);}
  bool hasBackwardPrediction() const {return (backwardPrediction_.get() != nullptr);}
  bool hasForwardUpdate() const {return (forwardUpdate_.get() != nullptr);}
  bool hasBackwardUpdate() const {return (backwardUpdate_.get() != nullptr);}
  bool hasUpdate(int direction) const {if (direction < 0) return hasBackwardUpdate(); return hasForwardUpdate();}
  bool hasPredictionsAndUpdates() const {return (hasForwardPrediction() && hasBackwardPrediction() && hasForwardUpdate() && hasBackwardUpdate());}

  void setReferenceState(ReferenceStateOnPlane* referenceState);
  void setForwardPrediction(MeasuredStateOnPlane* forwardPrediction);
  void setBackwardPrediction(MeasuredStateOnPlane* backwardPrediction);
  void setPrediction(MeasuredStateOnPlane* prediction, int direction)  {if (direction >=0) setForwardPrediction(prediction); else setBackwardPrediction(prediction);}
  void setForwardUpdate(KalmanFittedStateOnPlane* forwardUpdate);
  void setBackwardUpdate(KalmanFittedStateOnPlane* backwardUpdate);
  void setUpdate(KalmanFittedStateOnPlane* update, int direction)  {if (direction >=0) setForwardUpdate(update); else setBackwardUpdate(update);}
  void setMeasurementsOnPlane(const std::vector< genfit::MeasurementOnPlane* >& measurementsOnPlane);
  void addMeasurementOnPlane(MeasurementOnPlane* measurementOnPlane);
  void addMeasurementsOnPlane(const std::vector< genfit::MeasurementOnPlane* >& measurementsOnPlane);
  //! Set weights of measurements.
  void setWeights(const std::vector<double>&);
  void fixWeights(bool arg = true) {fixWeights_ = arg;}
  void setRep(const AbsTrackRep* rep);

  void deleteForwardInfo();
  void deleteBackwardInfo();
  void deletePredictions();
  void deleteReferenceInfo() {setReferenceState(nullptr);}
  void deleteMeasurementInfo();

  virtual void Print(const Option_t* = "") const;

  virtual bool checkConsistency(const genfit::PruneFlags* = nullptr) const;

 private:

  //! Reference state. Used by KalmanFitterRefTrack.
  std::unique_ptr<ReferenceStateOnPlane> referenceState_; // Ownership
  std::unique_ptr<MeasuredStateOnPlane> forwardPrediction_; // Ownership
  std::unique_ptr<KalmanFittedStateOnPlane> forwardUpdate_; // Ownership
  std::unique_ptr<MeasuredStateOnPlane> backwardPrediction_; // Ownership
  std::unique_ptr<KalmanFittedStateOnPlane> backwardUpdate_; // Ownership
  mutable std::unique_ptr<MeasuredStateOnPlane> fittedStateUnbiased_; //!  cache
  mutable std::unique_ptr<MeasuredStateOnPlane> fittedStateBiased_; //!  cache

 //> TODO ! ptr implement: to the special ownership version
  /* class owned_pointer_vector : private std::vector<MeasuredStateOnPlane*> {
   public: 
    ~owned_pointer_vector() { for (size_t i = 0; i < this->size(); ++i)
                         delete this[i]; }
    size_t size() const { return this->size(); };
    void push_back(MeasuredStateOnPlane* measuredState) { this->push_back(measuredState); };
    const  MeasuredStateOnPlane* at(size_t i)  const { return this->at(i); }; 
	//owned_pointer_vector::iterator erase(owned_pointer_vector::iterator position) ;
	//owned_pointer_vector::iterator erase(owned_pointer_vector::iterator first, owned_pointer_vector::iterator last);
};
	*/

  std::vector<MeasurementOnPlane*> measurementsOnPlane_; // Ownership
  bool fixWeights_; // weights should not be altered by fitters anymore

 public:

  ClassDef(KalmanFitterInfo,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_KalmanFitterInfo_h
