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

#ifndef genfit_StateOnPlane_h
#define genfit_StateOnPlane_h

#include "SharedPlanePtr.h"
#include "AbsTrackRep.h"

#include <TObject.h>
#include <TVectorD.h>


namespace genfit {

/**
 * @brief A state with arbitrary dimension defined in a DetPlane.
 *
 * The dimension and meaning of the #state_ vector are defined by the track parameterization of the #rep_.
 * #sharedPlane_ is a shared_pointer, the ownership over that plane is shared between all StateOnPlane objects defined in that plane.
 * The definition of the state is bound to the TrackRep #rep_. Therefore, the StateOnPlane contains a pointer to a AbsTrackRep.
 * It will provide functionality to extrapolate it and translate the state it into cartesian coordinates.
 * Shortcuts to all functions of the AbsTrackRep which use this StateOnPlane are also provided here.
 */
class StateOnPlane : public TObject {

 public:


  StateOnPlane(const AbsTrackRep* rep = NULL);
  //! state is defined by the TrackReps parameterization
  StateOnPlane(const TVectorD& state, const SharedPlanePtr& plane, const AbsTrackRep* rep);
  StateOnPlane(const TVectorD& state, const SharedPlanePtr& plane, const AbsTrackRep* rep, const TVectorD& auxInfo);

  StateOnPlane& operator=(StateOnPlane other);
  void swap(StateOnPlane& other); // nothrow

  virtual ~StateOnPlane() {}

  const TVectorD& getState() const {return state_;}
  TVectorD& getState() {return state_;}
  const TVectorD& getAuxInfo() const {return auxInfo_;}
  TVectorD& getAuxInfo() {return auxInfo_;}
  const SharedPlanePtr& getPlane() const {return sharedPlane_;}
  const AbsTrackRep* getRep() const {return rep_;}

  void setState(const TVectorD& state) {if(state_.GetNrows() == 0) state_.ResizeTo(state); state_ = state;}
  void setPlane(const SharedPlanePtr& plane) {sharedPlane_ = plane;}
  void setStatePlane(const TVectorD& state, const SharedPlanePtr& plane) {state_ = state; sharedPlane_ = plane;}
  void setAuxInfo(const TVectorD& auxInfo) {if(auxInfo_.GetNrows() == 0) auxInfo_.ResizeTo(auxInfo); auxInfo_ = auxInfo;}
  void setRep(const AbsTrackRep* rep) {rep_ = rep;}

  // Shortcuts to TrackRep functions
  double extrapolateToPlane(const SharedPlanePtr& plane,
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToPlane(*this, plane, stopAtBoundary, calcJacobianNoise);}
  double extrapolateToLine(const TVector3& linePoint,
        const TVector3& lineDirection,
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToLine(*this, linePoint, lineDirection, stopAtBoundary, calcJacobianNoise);}
  double extrapolateToPoint(const TVector3& point,
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToPoint(*this, point, stopAtBoundary, calcJacobianNoise);}
  double extrapolateToPoint(const TVector3& point,
        const TMatrixDSym& G, // weight matrix (metric)
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToPoint(*this, point, G, stopAtBoundary, calcJacobianNoise);}
  double extrapolateToCylinder(double radius,
        const TVector3& linePoint = TVector3(0.,0.,0.),
        const TVector3& lineDirection = TVector3(0.,0.,1.),
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToCylinder(*this, radius, linePoint, lineDirection, stopAtBoundary, calcJacobianNoise);}
  double extrapolateToSphere(double radius,
        const TVector3& point = TVector3(0.,0.,0.),
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToSphere(*this, radius, point, stopAtBoundary, calcJacobianNoise);}
  double extrapolateBy(double step,
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateBy(*this, step, stopAtBoundary, calcJacobianNoise);}
  double extrapolateToMeasurement(const AbsMeasurement* measurement,
        bool stopAtBoundary = false,
        bool calcJacobianNoise = false) {return rep_->extrapolateToMeasurement(*this, measurement, stopAtBoundary, calcJacobianNoise);}


  TVector3 getPos() const {return rep_->getPos(*this);}
  TVector3 getMom() const {return rep_->getMom(*this);}
  TVector3 getDir() const {return rep_->getDir(*this);}
  void getPosMom(TVector3& pos, TVector3& mom) const {rep_->getPosMom(*this, pos, mom);}
  void getPosDir(TVector3& pos, TVector3& dir) const {rep_->getPosDir(*this, pos, dir);}
  TVectorD get6DState() const {return rep_->get6DState(*this);}
  double getMomMag() const {return rep_->getMomMag(*this);}
  int getPDG() const {return rep_->getPDG();}
  double getCharge() const {return rep_->getCharge(*this);}
  double getQop() const {return rep_->getQop(*this);}
  double getMass() const {return rep_->getMass(*this);}

  void setPosMom(const TVector3& pos, const TVector3& mom) {rep_->setPosMom(*this, pos, mom);}
  void setPosMom(const TVectorD& state6) {rep_->setPosMom(*this, state6);}
  void setChargeSign(double charge) {rep_->setChargeSign(*this, charge);}
  void setQop(double qop) {rep_->setQop(*this, qop);}


  virtual void Print(Option_t* option = "") const;

 protected:

  TVectorD state_; // state vector
  TVectorD auxInfo_; // auxiliary information (e.g. charge, flight direction etc.)
  SharedPlanePtr sharedPlane_; //! Shared ownership.  '!' in order to silence ROOT, custom streamer writes and reads this.

 private:

  /** Pointer to TrackRep with respect to which StateOnPlane is defined
   */
  const AbsTrackRep* rep_; //! No ownership

 public:
  ClassDef(StateOnPlane,1)

};


inline StateOnPlane::StateOnPlane(const AbsTrackRep* rep) :
  state_(0), auxInfo_(0), sharedPlane_(), rep_(rep)
{
  if (rep != NULL) {
    state_.ResizeTo(rep->getDim());
  }
}

inline StateOnPlane::StateOnPlane(const TVectorD& state, const SharedPlanePtr& plane, const AbsTrackRep* rep) :
  state_(state), auxInfo_(0), sharedPlane_(plane), rep_(rep)
{
  assert(rep != NULL);
  assert(sharedPlane_.get() != NULL);
}

inline StateOnPlane::StateOnPlane(const TVectorD& state, const SharedPlanePtr& plane, const AbsTrackRep* rep, const TVectorD& auxInfo) :
  state_(state), auxInfo_(auxInfo), sharedPlane_(plane), rep_(rep)
{
  assert(rep != NULL);
  assert(sharedPlane_.get() != NULL);
}

inline StateOnPlane& StateOnPlane::operator=(StateOnPlane other) {
  swap(other);
  return *this;
}

inline void StateOnPlane::swap(StateOnPlane& other) {
  this->state_.ResizeTo(other.state_);
  std::swap(this->state_, other.state_);
  this->auxInfo_.ResizeTo(other.auxInfo_);
  std::swap(this->auxInfo_, other.auxInfo_);
  this->sharedPlane_.swap(other.sharedPlane_);
  std::swap(this->rep_, other.rep_);
}

} /* End of namespace genfit */
/** @} */

#endif // genfit_StateOnPlane_h
