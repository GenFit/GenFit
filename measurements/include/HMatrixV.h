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

#ifndef genfit_HMatrixV_h
#define genfit_HMatrixV_h

#include "AbsHMatrix.h"


namespace genfit {

/**
 * @brief AbsHMatrix implementation for one-dimensional MeasurementOnPlane and RKTrackRep parameterization.
 *
 * This projects out v.
 * H = (0, 0, 0, 0, 1)
 */
class HMatrixV : public AbsHMatrix {

 public:

  HMatrixV() {;}

  const TMatrixD& getMatrix() const override;

  TVectorD Hv(const TVectorD& v) const override;

  TMatrixD MHt(const TMatrixDSym& M) const override;
  TMatrixD MHt(const TMatrixD& M) const override;

  void HMHt(TMatrixDSym& M) const override;

  virtual HMatrixV* clone() const override {return new HMatrixV(*this);}

  virtual bool isEqual(const AbsHMatrix& other) const override {return (dynamic_cast<const HMatrixV*>(&other) != nullptr);}

  virtual void Print(const Option_t* = "") const override;

  ClassDefOverride(HMatrixV,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_HMatrixV_h
