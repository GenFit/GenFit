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

#ifndef genfit_MeasurementProducer_h
#define genfit_MeasurementProducer_h

#include "Exception.h"
#include "TrackCand.h"

#include <assert.h>
#include <TClonesArray.h>


namespace genfit {

class AbsMeasurement;

/** @brief Abstract interface class for MeasurementProducer
 *
 * Defines the very basic interface of a producer.
 */
template <class measurement_T>
class AbsMeasurementProducer {
public:
  /** @brief Virtual abstract method to produce a Measurement.
   * Implemented in MeasurementProducer
   */
  virtual measurement_T* produce(int index, const TrackCandHit* hit) = 0;
  virtual ~AbsMeasurementProducer() {};
};


/** @brief Template class for a measurement producer module
 *
 * A MeasurementProducer module is used by MeasurementFactory to create Measurements for
 * one specific detector type.
 *
 * It is assumed that each detector has as output of its digitization /
 * clustering some sort of hit or cluster class which stores all information that
 * corresponds to a measured hit in that detector. The MeasurementProducer
 * converts this information into a class that can be handled by genfit.
 * This class is realized as a Measurement (a class inheriting from AbsMeasurement).
 *
 * In order to use the MeasurementProducer facility, a
 * Measurement has to implement a constructor which takes as an argument
 * a pointer to the cluster class and a TrackCandHit. This constructor serves as the initializing
 * constructor for the Measurement.
 *
 * The MeasurementProducer will fetch the cluster objects from a TClonesArray and
 * use the initializing constructor to build the corresponding Measurement.
 *
 * @param hit_t template parameter specifying hit/cluster class
 * @param measurement_T template parameter specifying Measurement
 */
template <class hit_T, class measurement_T>
class MeasurementProducer : public AbsMeasurementProducer<genfit::AbsMeasurement> {
 private:
  /** @brief pointer to array with cluster data */
  TClonesArray* hitArrayTClones_;

 public:
  /** @brief Constructor takes pointer to the hit array */
  MeasurementProducer(TClonesArray*);
  virtual ~MeasurementProducer();

  /** @brief Create a Measurement from the cluster at position index
   * in TClonesArray
   */
  virtual AbsMeasurement* produce(int index, const TrackCandHit* hit);
};


template <class hit_T, class measurement_T>
  MeasurementProducer<hit_T, measurement_T>::MeasurementProducer(TClonesArray* theArr) {
  hitArrayTClones_ = theArr;
  //std::cout << "hit array with " << hitArrayTClones_->GetEntries() << " entries." << std::endl;
}

template <class hit_T, class measurement_T>
MeasurementProducer<hit_T, measurement_T>::~MeasurementProducer() {
  // we don't assume ownership over the hit arrays
}

template <class hit_T, class measurement_T>
AbsMeasurement* MeasurementProducer<hit_T, measurement_T>::produce(int index, const TrackCandHit* hit) {
  assert(hitArrayTClones_ != NULL);
  //std::cout << "hit array with " << hitArrayTClones_->GetEntries() << " entries, looking for entry " << index << "." << std::endl;
  if(hitArrayTClones_->At(index) == 0) {
    Exception e("In MeasurementProducer: index for hit in TClonesArray out of bounds",__LINE__,__FILE__);
    e.setFatal();
    throw e;
  }
  return ( new measurement_T( (hit_T*) hitArrayTClones_->At(index), hit ) );
}


} /* End of namespace genfit */
/** @} */

#endif // genfit_MeasurementProducer_h
