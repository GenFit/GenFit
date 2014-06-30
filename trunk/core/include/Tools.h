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

#ifndef genfit_Tools_h
#define genfit_Tools_h

#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>

/**
 * @brief Matrix inversion tools.
 */
namespace genfit {
namespace tools {

  /** @brief Invert a matrix, throwing an Exception when inversion fails.
   * Optional calculation of determinant.
   */
  void invertMatrix(const TMatrixDSym& mat, TMatrixDSym& inv, double* determinant = NULL);
  /** @brief Same, replacing its argument.
   */
  void invertMatrix(TMatrixDSym& mat, double* determinant = NULL);

  /** @brief Solves R^t x = b, replacing b with the solution for x.  R is
   *  assumed to be upper diagonal.
   */
  bool transposedForwardSubstitution(const TMatrixD& R, TVectorD& b);
  /** @brief Same, for a column of the matrix b.  */
  bool transposedForwardSubstitution(const TMatrixD& R, TMatrixD& b, int nCol);
  /** @brief Inverts the transpose of the upper right matrix R into inv.  */
  bool transposedInvert(const TMatrixD& R, TMatrixD& inv);

  /** @brief Replaces A with an upper right matrix connected to A by
   *  an orthongonal transformation.  I.e., it computes R from a QR
   *  decomposition of A = QR, replacing A.
   */
  void QR(TMatrixD& A);

  /** @brief Replaces A with an upper right matrix connected to A by
   *  an orthongonal transformation.  I.e., it computes R from a QR
   *  decomposition of A = QR, replacing A.  Also replaces b by Q'b
   *  where Q' is the transposed of Q.
   */
  void QR(TMatrixD& A, TVectorD& b);

  /** @brief This averages the covariance matrices C1, C2 in a
   *  numerically stable way by using matrix square roots.  This code
   *  is in no way optimized so use with care if speed is a concern.
   */
void safeAverage(const TMatrixDSym& C1, const TMatrixDSym& C2,
		 TMatrixDSym& result);

} /* End of namespace tools */
} /* End of namespace genfit */
/** @} */

#endif // genfit_Tools_h
