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

#include "HMatrixUnit.h"
#include <cassert>
#include <alloca.h>


namespace genfit {


// 1, 0, 0, 0, 0,
// 0, 1, 0, 0, 0,
// 0, 0, 1, 0, 0,
// 0, 0, 0, 1, 0,
// 0, 0, 0, 0, 1

const TMatrixD& HMatrixUnit::getMatrix() const {
  static const double HMatrixContent[5*5] = {1, 0, 0, 0, 0,
                                               0, 1, 0, 0, 0,
                                               0, 0, 1, 0, 0,
                                               0, 0, 0, 1, 0,
                                               0, 0, 0, 0, 1};

  static const TMatrixD HMatrix(5,5, HMatrixContent);

  return HMatrix;
}


} /* End of namespace genfit */
