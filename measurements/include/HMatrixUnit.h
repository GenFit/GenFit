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

#ifndef genfit_HMatrixUnit_h
#define genfit_HMatrixUnit_h

#include "AbsHMatrix.h"


namespace genfit {

/**
 * @brief AbsHMatrix implementation for 5-dimensional MeasurementOnPlane and RKTrackRep parameterization.
 *
 * H = (1, 0, 0, 0, 0)
 *     (0, 1, 0, 0, 0)
 *     (0, 0, 1, 0, 0)
 *     (0, 0, 0, 1, 0)
 *     (0, 0, 0, 0, 1)
 */
class HMatrixUnit : public AbsHMatrix {

 public:

  HMatrixUnit() {;}

  const TMatrixD& getMatrix() const;

  TVectorD Hv(const TVectorD& v) const {return v;}

  TMatrixD MHt(const TMatrixDSym& M) const {return TMatrixD(M);}
  TMatrixD MHt(const TMatrixD& M) const {return M;}

  void HMHt(TMatrixDSym& M) const {return;}

  virtual AbsHMatrix* clone() const {return new HMatrixUnit(*this);}

  virtual bool isEqual(const AbsHMatrix& other) const {return (dynamic_cast<const HMatrixUnit*>(&other) != NULL);}

  ClassDef(HMatrixUnit,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_HMatrixUnit_h
