/* Copyright 2013, Technische Universitaet Muenchen, Ludwig-Maximilians-Universität München
   Authors: Johannes Rauch, Tobias Schlüter

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

#include "HMatrixV.h"

#include "IO.h"

#include <cassert>
#include <alloca.h>

namespace genfit {


// 0, 0, 0, 0, 1

const TMatrixD& HMatrixV::getMatrix() const {
  static const double HMatrixContent[5] = {0, 0, 0, 0, 1};

  static const TMatrixD HMatrix(1,5, HMatrixContent);

  return HMatrix;
}


TVectorD HMatrixV::Hv(const TVectorD& v) const {
  assert (v.GetNrows() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * 1);

  retValArray[0] = v(4); // v

  return TVectorD(1, retValArray);
}


TMatrixD HMatrixV::MHt(const TMatrixDSym& M) const {
  assert (M.GetNcols() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * 5);
  const double* MatArray = M.GetMatrixArray();

  for (unsigned int i=0; i<5; ++i) {
    retValArray[i] = MatArray[i*5 + 4];
  }

  return TMatrixD(5,1, retValArray);
}


TMatrixD HMatrixV::MHt(const TMatrixD& M) const {
  assert (M.GetNcols() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * M.GetNrows());
  const double* MatArray = M.GetMatrixArray();

  for (int i = 0; i < M.GetNrows(); ++i) {
    retValArray[i] = MatArray[i*5 + 4];
  }

  return TMatrixD(M.GetNrows(),1, retValArray);
}


void HMatrixV::HMHt(TMatrixDSym& M) const {
  assert (M.GetNrows() == 5);

  M(0,0) = M(4,4);

  M.ResizeTo(1,1);
}


void HMatrixV::Print(const Option_t*) const {
  printOut << "V" << std::endl;
}

} /* End of namespace genfit */
