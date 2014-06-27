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

#ifndef genfit_AbsFitterInfo_h
#define genfit_AbsFitterInfo_h

#include "MeasurementOnPlane.h"
#include "FitStatus.h"

#include <TObject.h>
#include <TVectorD.h>


namespace genfit {

class AbsTrackRep;
class TrackPoint;

/**
 *  @brief This class collects all information needed and produced by a specific  AbsFitter and is specific to one AbsTrackRep of the Track.
 */
class AbsFitterInfo : public TObject {

 public:

  AbsFitterInfo();
  AbsFitterInfo(const TrackPoint* trackPoint, const AbsTrackRep* rep);

  virtual ~AbsFitterInfo() {};

  //! Deep copy ctor for polymorphic class.
  virtual AbsFitterInfo* clone() const = 0;

  const TrackPoint* getTrackPoint() const {return trackPoint_;}
  const AbsTrackRep* getRep() const {return rep_;}

  void setTrackPoint(const TrackPoint *tp) {trackPoint_ = tp;}
  virtual void setRep(const AbsTrackRep* rep) {rep_ = rep;}

  virtual bool hasMeasurements() const = 0;
  virtual bool hasReferenceState() const = 0;
  virtual bool hasForwardPrediction() const = 0;
  virtual bool hasBackwardPrediction() const = 0;
  virtual bool hasPrediction(int direction) const {if (direction >=0) return hasForwardPrediction(); return hasBackwardPrediction();}
  virtual bool hasForwardUpdate() const = 0;
  virtual bool hasBackwardUpdate() const = 0;
  virtual bool hasUpdate(int direction) const {if (direction >=0) return hasForwardUpdate(); return hasBackwardUpdate();}

  virtual void deleteForwardInfo() = 0;
  virtual void deleteBackwardInfo() = 0;
  virtual void deleteReferenceInfo() = 0;
  virtual void deleteMeasurementInfo() = 0;

  const SharedPlanePtr& getPlane() const {return sharedPlane_;}
  virtual const MeasuredStateOnPlane& getFittedState(bool biased = true) const = 0;
  virtual MeasurementOnPlane getResidual(unsigned int iMeasurement = 0, bool biased = true, bool onlyMeasurementErrors = false) const = 0;

  void setPlane(const SharedPlanePtr& plane) {sharedPlane_ = plane;}

  virtual void Print(const Option_t* = "") const {;}

  virtual bool checkConsistency(const PruneFlags* = NULL) const = 0;

 protected:

  /** Pointer to TrackPoint where the FitterInfo belongs to
   */
  const TrackPoint* trackPoint_; //! No ownership

  /** Pointer to AbsTrackRep with respect to which the FitterInfo is defined
   */
  const AbsTrackRep* rep_; //! No ownership

  SharedPlanePtr sharedPlane_; //! Shared ownership.  '!' shuts up ROOT.


 private:
  AbsFitterInfo(const AbsFitterInfo&); // copy constructor
  AbsFitterInfo& operator=(const AbsFitterInfo&); // assignment operator


 public:
  ClassDef(AbsFitterInfo,1)

};

//! Needed for boost cloneability:
inline AbsFitterInfo* new_clone( const AbsFitterInfo & a)
{
  return a.clone();
}

} /* End of namespace genfit */
/** @} */

#endif // genfit_AbsFitterInfo_h
