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

#include "HMatrixUV.h"

#include "IO.h"

#include <cassert>
#include <alloca.h>

namespace genfit {


// 0, 0, 0, 1, 0
// 0, 0, 0, 0, 1

const TMatrixD& HMatrixUV::getMatrix() const {
  static const double HMatrixContent[2*5] = {0, 0, 0, 1, 0,
                                             0, 0, 0, 0, 1};

  static const TMatrixD HMatrix(2,5, HMatrixContent);

  return HMatrix;
}


TVectorD HMatrixUV::Hv(const TVectorD& v) const {
  assert (v.GetNrows() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * 2);
  const double* VecArray = v.GetMatrixArray();

  retValArray[0] = VecArray[3]; // u
  retValArray[1] = VecArray[4]; // v

  return TVectorD(2, retValArray);
}


TMatrixD HMatrixUV::MHt(const TMatrixDSym& M) const {
  assert (M.GetNcols() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * 5*2);
  const double* MatArray = M.GetMatrixArray();

  for (unsigned int i=0; i<5; ++i) {
    retValArray[i*2] = MatArray[i*5 + 3];
    retValArray[i*2 + 1] = MatArray[i*5 + 4];
  }

  return TMatrixD(5,2, retValArray);
}


TMatrixD HMatrixUV::MHt(const TMatrixD& M) const {
  assert (M.GetNcols() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * M.GetNrows()*2);
  const double* MatArray = M.GetMatrixArray();

  for (int i = 0; i < M.GetNrows(); ++i) {
    retValArray[i*2] = MatArray[i*5 + 3];
    retValArray[i*2 + 1] = MatArray[i*5 + 4];
  }

  return TMatrixD(M.GetNrows(),2, retValArray);
}


void HMatrixUV::HMHt(TMatrixDSym& M) const {
  assert (M.GetNrows() == 5);
  double* MatArray = M.GetMatrixArray();

  //
  //  HMH^t = ( M_33  M_34 ) where M_34 == M_43
  //          ( M_43  M_44 )
  //
  double uu = MatArray[3*5 + 3];
  double uv = MatArray[3*5 + 4];
  double vv = MatArray[4*5 + 4];

  M.ResizeTo(2,2);
  MatArray = M.GetMatrixArray();
  MatArray[0] = uu; MatArray[1] = uv;
  MatArray[2] = uv; MatArray[3] = vv;
}


void HMatrixUV::Print(const Option_t*) const {
  printOut << "UV" << std::endl;
}


} /* End of namespace genfit */
