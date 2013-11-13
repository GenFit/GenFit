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

#ifndef genfit_AbsHMatrix_h
#define genfit_AbsHMatrix_h

#include <TObject.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>


namespace genfit {

/**
 * @brief HMatrix for projecting from AbsTrackRep parameters to measured parameters in a DetPlane.
 *
 */
class AbsHMatrix : public TObject {

 public:

  AbsHMatrix() {;}

  virtual ~AbsHMatrix() {;}

  //! Get the actual matrix representation
  virtual const TMatrixD& getMatrix() const = 0;

  //! H*v
  virtual TVectorD Hv(const TVectorD& v) const {return getMatrix()*v;}

  //! M*H^t
  virtual TMatrixD MHt(const TMatrixDSym& M) const {return TMatrixD(M, TMatrixD::kMultTranspose, getMatrix());}
  virtual TMatrixD MHt(const TMatrixD& M) const {return TMatrixD(M, TMatrixD::kMultTranspose, getMatrix());}

  //! similarity: H*M*H^t
  virtual void HMHt(TMatrixDSym& M) const {M.Similarity(getMatrix());}

  virtual AbsHMatrix* clone() const = 0;

  bool operator==(const AbsHMatrix& other) const {return this->isEqual(other);}
  bool operator!=(const AbsHMatrix& other) const {return !(this->isEqual(other));}
  virtual bool isEqual(const AbsHMatrix& other) const = 0;

  virtual void Print(const Option_t* = "") const {;}

 protected:
  // protect from calling copy c'tor or assignment operator from outside the class. Use #clone() if you want a copy!
  AbsHMatrix(const AbsHMatrix& o) : TObject(o) {;}
  AbsHMatrix& operator=(const AbsHMatrix&);

 public:
  ClassDef(AbsHMatrix,1)

};

} /* End of namespace genfit */
/** @} */

#endif // genfit_AbsHMatrix_h
