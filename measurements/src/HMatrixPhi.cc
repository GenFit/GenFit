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

#include "HMatrixPhi.h"

#include "IO.h"

#include <TBuffer.h>

#include <cassert>
#include <alloca.h>
#include <math.h>
#include <TBuffer.h>

namespace genfit {


// 0, 0, 0, cos(phi), sin(phi)


HMatrixPhi::HMatrixPhi(double phi) :
  phi_(phi),
  cosPhi_(cos(phi)),
  sinPhi_(sin(phi))
{
  ;
}

const TMatrixD& HMatrixPhi::getMatrix() const {
  static const double HMatrixContent[5] = {0, 0, 0, cosPhi_, sinPhi_};

  static const TMatrixD HMatrix(1,5, HMatrixContent);

  return HMatrix;
}


TVectorD HMatrixPhi::Hv(const TVectorD& v) const {
  assert (v.GetNrows() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * 1);

  retValArray[0] = cosPhi_*v(3) + sinPhi_*v(4);

  return TVectorD(1, retValArray);
}


TMatrixD HMatrixPhi::MHt(const TMatrixDSym& M) const {
  assert (M.GetNcols() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * 5);
  const double* MatArray = M.GetMatrixArray();

  for (unsigned int i=0; i<5; ++i) {
    retValArray[i] = cosPhi_*MatArray[i*5 + 3] + sinPhi_*MatArray[i*5 + 4];
  }

  return TMatrixD(5,1, retValArray);
}


TMatrixD HMatrixPhi::MHt(const TMatrixD& M) const {
  assert (M.GetNcols() == 5);

  double* retValArray =(double *)alloca(sizeof(double) * M.GetNrows());
  const double* MatArray = M.GetMatrixArray();

  for (int i = 0; i < M.GetNrows(); ++i) {
    retValArray[i] = cosPhi_*MatArray[i*5 + 3] + sinPhi_*MatArray[i*5 + 4];
  }

  return TMatrixD(M.GetNrows(),1, retValArray);
}


void HMatrixPhi::HMHt(TMatrixDSym& M) const {
  assert (M.GetNrows() == 5);

  M(0,0) =   cosPhi_ * (cosPhi_*M(3,3) + sinPhi_*M(3,4))
           + sinPhi_ * (cosPhi_*M(4,3) + sinPhi_*M(4,4));

  M.ResizeTo(1,1);
}


bool HMatrixPhi::isEqual(const AbsHMatrix& other) const {
  if (dynamic_cast<const HMatrixPhi*>(&other) == NULL)
    return false;

  return (phi_ == static_cast<const HMatrixPhi*>(&other)->phi_);
}

void HMatrixPhi::Print(const Option_t*) const
{
  printOut << "phi = " << phi_ << std::endl;
}

void HMatrixPhi::Streamer(TBuffer &R__b) {
  // Stream an object of class genfit::HMatrixPhi.

  // Modified from auto-generated streamer to set non-persistent members after reading

  if (R__b.IsReading()) {
    R__b.ReadClassBuffer(genfit::HMatrixPhi::Class(),this);
    cosPhi_ = cos(phi_);
    sinPhi_ = sin(phi_);
  } else {
    R__b.WriteClassBuffer(genfit::HMatrixPhi::Class(),this);
  }
}


} /* End of namespace genfit */
