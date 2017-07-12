/* Copyright 2008-2009, Technische Universitaet Muenchen,
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

/** @addtogroup RKTrackRep
 * @{
 */

#ifndef genfit_RKTools_h
#define genfit_RKTools_h

#include <stddef.h>
#include <algorithm>

namespace genfit {

template <size_t nRows, size_t nCols>
struct RKMatrix {
  double vals[nRows * nCols];

  double& operator()(size_t iRow, size_t iCol) {
    return vals[nCols*iRow + iCol];
  }
  double& operator[](size_t n) {
    return vals[n];
  }
  const double& operator[](size_t n) const {
    return vals[n];
  }
  double* begin() { return vals; }
  double* end() { return vals + nRows * nCols; }
  const double* begin() const { return vals; }
  const double* end() const { return vals + nRows * nCols; }
  RKMatrix<nRows, nCols>& operator=(const RKMatrix<nRows, nCols>& o) {
    std::copy(o.begin(), o.end(), this->begin());
    return *this;
  }

};

typedef RKMatrix<1, 3> M1x3;
typedef RKMatrix<1, 4> M1x4;
typedef RKMatrix<1, 7> M1x7;
typedef RKMatrix<5, 5> M5x5;
typedef RKMatrix<6, 6> M6x6;
typedef RKMatrix<7, 7> M7x7;
typedef RKMatrix<6, 5> M6x5;
typedef RKMatrix<7, 5> M7x5;
typedef RKMatrix<5, 6> M5x6;
typedef RKMatrix<5, 7> M5x7;

} /* End of namespace genfit */
/** @} */

#endif // genfit_RKTools_h

