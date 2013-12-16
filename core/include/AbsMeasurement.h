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

#ifndef genfit_AbsMeasurement_h
#define genfit_AbsMeasurement_h

#include "MeasurementOnPlane.h"
#include "AbsHMatrix.h"

#include <TObject.h>


namespace genfit {

class AbsTrackRep;
class TrackPoint;

/**
 *  @brief Contains the measurement and covariance in raw detector coordinates.
 *
 *  Detector and hit ids can be used to point back to the original detector hits (clusters etc.).
 */
class AbsMeasurement : public TObject {

 public:

  AbsMeasurement() : rawHitCoords_(), rawHitCov_(), detId_(-1), hitId_(-1) {;}
  AbsMeasurement(int nDims) : rawHitCoords_(nDims), rawHitCov_(nDims), detId_(-1), hitId_(-1) {;}
  AbsMeasurement(const TVectorD& rawHitCoords, const TMatrixDSym& rawHitCov, int detId, int hitId, TrackPoint* trackPoint);

  virtual ~AbsMeasurement();

  //! Deep copy ctor for polymorphic class.
  virtual AbsMeasurement* clone() const = 0;

  TrackPoint* getTrackPoint() const {return trackPoint_;}
  void setTrackPoint(TrackPoint* tp) {trackPoint_ = tp;}

  const TVectorD& getRawHitCoords() const {return rawHitCoords_;}
  const TMatrixDSym& getRawHitCov() const {return rawHitCov_;}
  int getDetId() const {return detId_;}
  int getHitId() const {return hitId_;}

  unsigned int getDim() const {return rawHitCoords_.GetNrows();}

  void setDetId(int detId) {detId_ = detId;}
  void setHitId(int hitId) {hitId_ = hitId;}


  /**
   * Construct (virtual) detector plane (use state's AbsTrackRep).
   * It's possible to make corrections to the plane here.
   * The state should be defined somewhere near the measurement.
   * For virtual planes, the state will be extrapolated to the POCA to point (SpacepointMeasurement)
   * or line (WireMeasurement), and from this info the plane will be constructed.
   */
  virtual SharedPlanePtr constructPlane(const StateOnPlane& state) const = 0;

  /**
   * Construct MeasurementOnPlane on plane of the state
   * and wrt the states TrackRep.
   * The state will usually be the prediction or reference state,
   * and has to be defined AT the measurement.
   * The AbsMeasurement will be projected onto the plane.
   * It's possible to make corrections to the coordinates here (e.g. by using the state coordinates).
   * Usually the vector will contain only one element. But in the case of e.g. a WireMeasurement, it will be 2 (left and right).
   */
  virtual std::vector<genfit::MeasurementOnPlane*> constructMeasurementsOnPlane(const StateOnPlane& state) const = 0;

  /**
   * Returns a new AbsHMatrix object. Caller must take ownership.
   */
  virtual const AbsHMatrix* constructHMatrix(const AbsTrackRep*) const = 0;

  virtual void Print(const Option_t* = "") const;


 private:
  //! protect from calling assignment operator from outside the class. Use #clone() if you want a copy!
  virtual AbsMeasurement& operator=(const AbsMeasurement&); // default cannot work because TVector and TMatrix = operators don't do resizing

 protected:
  //! protect from calling copy c'tor from outside the class. Use #clone() if you want a copy!
  AbsMeasurement(const AbsMeasurement&);

  TVectorD rawHitCoords_;
  TMatrixDSym rawHitCov_;
  int detId_; // detId id is -1 per default
  int hitId_; // hitId id is -1 per default

  //! Pointer to TrackPoint where the measurement belongs to
  TrackPoint* trackPoint_; //! No ownership

 public:
  ClassDef(AbsMeasurement, 1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_AbsMeasurement_h
