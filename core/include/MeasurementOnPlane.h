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

#ifndef genfit_MeasurementOnPlane_h
#define genfit_MeasurementOnPlane_h

#include "MeasuredStateOnPlane.h"
#include "AbsHMatrix.h"

#include <TMatrixD.h>

#include <cmath>
#include <memory>


namespace genfit {

/**
 * @brief Measured coordinates on a plane.
 *
 * The dimensionality will usually be 1 or 2.
 * The HMatrix corresponds to a projection matrix, which is used to
 * project the track parameters with the original dimensionality down
 * to the measured dimensionality.
 */
class MeasurementOnPlane : public MeasuredStateOnPlane {

 public:

  MeasurementOnPlane(const AbsTrackRep* rep = nullptr) :
    MeasuredStateOnPlane(rep), hMatrix_(nullptr), weight_(0) {}
  MeasurementOnPlane(const TVectorD& state, const TMatrixDSym& cov, SharedPlanePtr plane, const AbsTrackRep* rep, const AbsHMatrix* hMatrix, double weight = 1.) :
    MeasuredStateOnPlane(state, cov, plane, rep), hMatrix_(hMatrix), weight_(weight) {}

  //! copy constructor
  MeasurementOnPlane(const MeasurementOnPlane& other);
  //! assignment operator
  MeasurementOnPlane& operator=(MeasurementOnPlane other);
  void swap(MeasurementOnPlane& other);

  virtual ~MeasurementOnPlane() {}

  const AbsHMatrix* getHMatrix() const {return hMatrix_.get();}
  double getWeight() const {return weight_;}

  TMatrixDSym getWeightedCov() {return weight_*cov_;}

  void setHMatrix(const AbsHMatrix* hMatrix) {hMatrix_.reset(hMatrix);}
  void setWeight(double weight) {weight_ = fmax(weight, 1.E-10);}

  void Print(Option_t* option = "") const ;

 private:
  TVector3 getPos() const;
  TVector3 getMom() const;
  TVector3 getDir() const;
  void getPosMom(TVector3& pos, TVector3& mom) const;
  void getPosDir(TVector3& pos, TVector3& dir) const;
  TVectorD get6DState() const;
  double getMomMag() const;
  int getPDG() const;
  double getCharge() const;
  double getQop() const;
  double getMass() const;
  double getTime() const;

  void setPosMom(const TVector3& pos, const TVector3& mom);
  void setPosMom(const TVectorD& state6);
  void setChargeSign(double charge);
  void setQop(double qop);
  void setTime(double time);


 protected:

  std::unique_ptr<const AbsHMatrix> hMatrix_; // Ownership
  double weight_;

 public:
  ClassDef(MeasurementOnPlane,1)

};

} /* End of namespace   */
/** @} */

#endif //  _MeasurementOnPlane_h
