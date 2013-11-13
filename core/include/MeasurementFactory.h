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


#ifndef genfit_MeasurementFactory_h
#define genfit_MeasurementFactory_h

#include "MeasurementProducer.h"
#include "TrackCand.h"

#include <vector>
#include <map>


namespace genfit {

class AbsMeasurement;

/** @brief Factory object to create AbsMeasurement objects from digitized and clustered data
 *
 * The MeasurementFactory is used to automatically fill Track objects with
 * hit data. For each detector type that is used, one AbsMeasurementProducer
 * has to be registered in the factory. The factory can the use the index
 * information from a TrackCand object to load the indexed hits into
 * the Track.
 *
 * @sa AbsMeasurementProducer
 * @sa TrackCand
 */
template <class measurement_T>
class MeasurementFactory {
 private:
  std::map<int, AbsMeasurementProducer<measurement_T>*> hitProdMap_;


 public:
  MeasurementFactory() {};
  virtual ~MeasurementFactory() { clear(); }

  /** @brief Register a producer module to the factory
   *
   * For each type of hit a separate producer is needed. The type of hit
   * is identified by the detector ID (detID). This index corresponds to the
   * detector ID that is stored in the TrackCand.
   */
  void addProducer(int detID, AbsMeasurementProducer<measurement_T>* hitProd);

  /** @brief Clear all hit producers
   */
  void clear();

  /** @brief Create a Measurement
   *
   * Measurements have to implement a Constructor which takes the cluster object
   * from which the Measurement is build as the only parameter.
   * @sa AbsMeasurementProducer
   */
  measurement_T* createOne (int detID, int index, const TrackCandHit* hit) const;

  /** @brief Create a collection of Measurements
   *
   * This is the standard way to prepare the hit collection for a Track. The
   * resulting collection can contain hits from several detectors. The order
   * of the hits is the same as in the TrackCand. It is assumed that this order
   * is already along the Track.
   *
   * Measurements have to implement a constructor which takes the cluster object
   * from which the Measurement is build as the only parameter.
   * @sa AbsMeasurementProducer
   */
  std::vector<measurement_T*> createMany(const TrackCand& cand) const;

};


template <class measurement_T>
void MeasurementFactory<measurement_T>::addProducer(int detID, AbsMeasurementProducer<measurement_T>* hitProd) {
  typename std::map<int, AbsMeasurementProducer<measurement_T>*>::iterator it = hitProdMap_.find(detID);
  if(it == hitProdMap_.end()) {
    hitProdMap_[detID] = hitProd;
  } else {
    Exception exc("MeasurementFactory: detID already in use",__LINE__,__FILE__);
    exc.setFatal();
    std::vector<double> numbers;
    numbers.push_back(detID);
    exc.setNumbers("detID",numbers);
    throw exc;
  }
}

template <class measurement_T>
void MeasurementFactory<measurement_T>::clear(){
  typename std::map<int, AbsMeasurementProducer<measurement_T>*>::iterator it=hitProdMap_.begin();
  while(it!=hitProdMap_.end()){
    delete it->second;
    ++it;
  }
  hitProdMap_.clear();
}

template <class measurement_T>
measurement_T* MeasurementFactory<measurement_T>::createOne(int detID, int index, const TrackCandHit* hit) const {
  typename std::map<int, AbsMeasurementProducer<measurement_T>*>::const_iterator it = hitProdMap_.find(detID);

  if(it != hitProdMap_.end()) {
    return it->second->produce(index, hit);
  } else {
    Exception exc("MeasurementFactory: no hitProducer for this detID available",__LINE__,__FILE__);
    exc.setFatal();
    std::vector<double> numbers;
    numbers.push_back(detID);
    exc.setNumbers("detID", numbers);
    throw exc;
  }
}

template <class measurement_T>
typename std::vector<measurement_T*> MeasurementFactory<measurement_T>::createMany(const TrackCand& cand) const {
  typename std::vector<measurement_T*> hitVec;
  unsigned int nHits=cand.getNHits();
  for(unsigned int i=0;i<nHits;i++) {
    int detID, index;
    const TrackCandHit* hit = cand.getHit(i);
    cand.getHit(i, detID, index);
    hitVec.push_back( MeasurementFactory<measurement_T>::createOne(hit->getDetId(), hit->getHitId(), hit) );
  }
  return hitVec;
}


} /* End of namespace genfit */
/** @} */

#endif // genfit_MeasurementFactory_h
