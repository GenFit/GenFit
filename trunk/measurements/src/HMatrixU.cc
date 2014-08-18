/* Copyright 2013, Technische Universitaet Muenchen,
   Authors: Johannes Rauch

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

#include "HMatrixU.h"
#include <cassert>
#include <alloca.h>
#include <iostream>


namespace genfit {


// 0, 0, 0, 1, 0

const TMatrixD& HMatrixU::getMatrix() const {
  static const double HMatrixContent[5] = {0, 0, 0, 1, 0};

  static const TMatrixD HMatrix(1,5, HMatrixContent);

  return HMatrix;
}


TVectorD HMatrixU::Hv(const TVectorD& v) const {
  assert (v.GetNrows() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * 1);

  retValArray[0] = v(3); // u

  return TVectorD(1, retValArray);
}


TMatrixD HMatrixU::MHt(const TMatrixDSym& M) const {
  assert (M.GetNcols() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * 5);
  const double* MatArray = M.GetMatrixArray();

  for (unsigned int i=0; i<5; ++i) {
    retValArray[i] = MatArray[i*5 + 3];
  }

  return TMatrixD(5,1, retValArray);
}


TMatrixD HMatrixU::MHt(const TMatrixD& M) const {
  assert (M.GetNcols() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * M.GetNrows());
  const double* MatArray = M.GetMatrixArray();

  for (int i = 0; i < M.GetNrows(); ++i) {
    retValArray[i] = MatArray[i*5 + 3];
  }

  return TMatrixD(M.GetNrows(),1, retValArray);
}


void HMatrixU::HMHt(TMatrixDSym& M) const {
  assert (M.GetNrows() == 5);

  M(0,0) = M(3,3);

  M.ResizeTo(1,1);
}


void HMatrixU::Print(const Option_t*) const {
  std::cout << "U" << std::endl;
}


} /* End of namespace genfit */
