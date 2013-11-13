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

#ifndef genfit_ReferenceStateOnPlane_h
#define genfit_ReferenceStateOnPlane_h

#include "StateOnPlane.h"


namespace genfit {

/**
 * @brief #StateOnPlane with linearized transport to that #ReferenceStateOnPlane from previous and next #ReferenceStateOnPlane
 *
 * Transport matrices describe transport TO that plane.
 * We have transport matrix F, noise matrix N and delta state c.
 * Now, state p and covariance C follow this mathematics:
 *
 * p = F * p_old + c
 * C = F * C_old * F^T + N
 *
 */
class ReferenceStateOnPlane : public StateOnPlane {

 public:

  ReferenceStateOnPlane();
  ReferenceStateOnPlane(const TVectorD& state,
      const SharedPlanePtr& plane,
      const AbsTrackRep* rep);
  ReferenceStateOnPlane(const TVectorD& state,
      const SharedPlanePtr& plane,
      const AbsTrackRep* rep,
      const TVectorD& auxInfo);
  ReferenceStateOnPlane(const StateOnPlane& state);

  StateOnPlane& operator=(ReferenceStateOnPlane other);
  void swap(ReferenceStateOnPlane& other); // nothrow

  virtual ~ReferenceStateOnPlane() {}

  void setForwardSegmentLength(double len) {forwardSegmentLength_ = len;}
  void setBackwardSegmentLength(double len) {backwardSegmentLength_ = len;}
  void setForwardTransportMatrix(const TMatrixD& mat) {forwardTransportMatrix_.ResizeTo(mat); forwardTransportMatrix_=mat;}
  void setBackwardTransportMatrix(const TMatrixD& mat) {backwardTransportMatrix_.ResizeTo(mat); backwardTransportMatrix_=mat;}
  void setTransportMatrix(const TMatrixD& mat, int direction) {if (direction >= 0) setForwardTransportMatrix(mat); else setBackwardTransportMatrix(mat);}
  void setForwardNoiseMatrix(const TMatrixDSym& mat) {forwardNoiseMatrix_.ResizeTo(mat); forwardNoiseMatrix_=mat;}
  void setBackwardNoiseMatrix(const TMatrixDSym& mat) {backwardNoiseMatrix_.ResizeTo(mat); backwardNoiseMatrix_=mat;}
  void setNoiseMatrix(const TMatrixDSym& mat, int direction) {if (direction >= 0) setForwardNoiseMatrix(mat); else setBackwardNoiseMatrix(mat);}
  void setForwardDeltaState(const TVectorD& mat) {forwardDeltaState_.ResizeTo(mat); forwardDeltaState_=mat;}
  void setBackwardDeltaState(const TVectorD& mat) {backwardDeltaState_.ResizeTo(mat); backwardDeltaState_=mat;}
  void setDeltaState(const TVectorD& mat, int direction) {if (direction >= 0) setForwardDeltaState(mat); else setBackwardDeltaState(mat);}


  double getForwardSegmentLength() const {return forwardSegmentLength_;}
  double getBackwardSegmentLength() const {return backwardSegmentLength_;}
  const TMatrixD& getForwardTransportMatrix() const {return forwardTransportMatrix_;}
  const TMatrixD& getBackwardTransportMatrix() const {return backwardTransportMatrix_;}
  const TMatrixD& getTransportMatrix(int direction) const {if (direction >= 0) return forwardTransportMatrix_; return backwardTransportMatrix_;}
  const TMatrixDSym& getForwardNoiseMatrix() const {return forwardNoiseMatrix_;}
  const TMatrixDSym& getBackwardNoiseMatrix() const {return backwardNoiseMatrix_;}
  const TMatrixDSym& getNoiseMatrix(int direction) const {if (direction >= 0) return forwardNoiseMatrix_; return backwardNoiseMatrix_;}
  const TVectorD& getForwardDeltaState() const {return forwardDeltaState_;}
  const TVectorD& getBackwardDeltaState() const {return backwardDeltaState_;}
  const TVectorD& getDeltaState(int direction) const {if (direction >= 0) return forwardDeltaState_; return backwardDeltaState_;}

  void resetForward();
  void resetBackward();

  virtual void Print(Option_t* option = "") const;

 protected:

  double forwardSegmentLength_; /**< Segment length from previous referenceState */
  double backwardSegmentLength_; /**< Segment length from next referenceState */
  TMatrixD forwardTransportMatrix_; /**< transport matrix F from previous referenceState */
  TMatrixD backwardTransportMatrix_; /**< transport matrix F from next referenceState */
  TMatrixDSym forwardNoiseMatrix_; /**< noise matrix N for transport from previous referenceState */
  TMatrixDSym backwardNoiseMatrix_; /**< noise matrix N for transport from next referenceState */
  TVectorD forwardDeltaState_; /**< c */
  TVectorD backwardDeltaState_; /**< c */


 public:

  ClassDef(ReferenceStateOnPlane,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_ReferenceStateOnPlane_h
