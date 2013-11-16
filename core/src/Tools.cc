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

#include "Tools.h"

#include <cmath>
#include <memory>
#include <typeinfo>
#include <cassert>

#include <TDecompChol.h>
#include <TMatrixTSymCramerInv.h>
#include <TMath.h>

#include "Exception.h"

// Use Cramer inversion for small matrices?
static const bool useCramer = false;

namespace genfit {

void tools::invertMatrix(const TMatrixDSym& mat, TMatrixDSym& inv, double* determinant){
  inv.ResizeTo(mat);

  // check if numerical limits are reached (i.e at least one entry < 1E-100 and/or at least one entry > 1E100)
  if (!(mat<1.E100) || !(mat>-1.E100)){
    Exception e("Tools::invertMatrix() - cannot invert matrix, entries too big (>1e100)",
        __LINE__,__FILE__);
    e.setFatal();
    throw e;
  }
  // do the trivial inversions for 1x1 and 2x2 matrices manually
  if (mat.GetNrows() == 1){
    if (determinant != NULL) *determinant = mat(0,0);
    inv(0,0) = 1./mat(0,0);
    return;
  }

  if (mat.GetNrows() == 2){
    double det = mat(0,0)*mat(1,1) - mat(1,0)*mat(1,0);
    if (determinant != NULL) *determinant = det;
    if(fabs(det) < 1E-50){
      Exception e("Tools::invertMatrix() - cannot invert matrix , determinant = 0",
          __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }
    det = 1./det;
    inv(0,0) =             det * mat(1,1);
    inv(0,1) = inv(1,0) = -det * mat(1,0);
    inv(1,1) =             det * mat(0,0);
    return;
  }


  if (useCramer && mat.GetNrows() <= 6){
    Bool_t (*inversion)(TMatrixDSym&, Double_t*) = 0;
    inv.ResizeTo(mat);
    inv = mat;
    switch (mat.GetNrows()) {
    case 3:
      inversion = TMatrixTSymCramerInv::Inv3x3; break;
    case 4:
      inversion = TMatrixTSymCramerInv::Inv4x4; break;
    case 5:
      inversion = TMatrixTSymCramerInv::Inv5x5; break;
    case 6:
      inversion = TMatrixTSymCramerInv::Inv6x6; break;
    }

    Bool_t success = inversion(inv, determinant);
    if (!success){
      Exception e("Tools::invertMatrix() - cannot invert matrix, determinant = 0",
          __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }
    return;
  }

  // else use TDecompChol
  bool status = 0;
  TDecompChol invertAlgo(mat, 1E-50);

  status = invertAlgo.Invert(inv);
  if(status == 0){
    Exception e("Tools::invertMatrix() - cannot invert matrix, status = 0",
        __LINE__,__FILE__);
    e.setFatal();
    throw e;
  }

  if (determinant != NULL) {
    double d1, d2;
    invertAlgo.Det(d1, d2);
    *determinant = ldexp(d1, d2);
  }
}

void tools::invertMatrix(TMatrixDSym& mat, double* determinant){
  // check if numerical limits are reached (i.e at least one entry < 1E-100 and/or at least one entry > 1E100)
  if (!(mat<1.E100) || !(mat>-1.E100)){
    Exception e("Tools::invertMatrix() - cannot invert matrix, entries too big (>1e100)",
        __LINE__,__FILE__);
    e.setFatal();
    throw e;
  }
  // do the trivial inversions for 1x1 and 2x2 matrices manually
  if (mat.GetNrows() == 1){
    if (determinant != NULL) *determinant = mat(0,0);
    mat(0,0) = 1./mat(0,0);
    return;
  }

  if (mat.GetNrows() == 2){
    double *arr = mat.GetMatrixArray();
    double det = arr[0]*arr[3] - arr[1]*arr[1];
    if (determinant != NULL) *determinant = det;
    if(fabs(det) < 1E-50){
      Exception e("Tools::invertMatrix() - cannot invert matrix, determinant = 0",
          __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }
    det = 1./det;
    double temp[3];
    temp[0] =  det * arr[3];
    temp[1] = -det * arr[1];
    temp[2] =  det * arr[0];
    //double *arr = mat.GetMatrixArray();
    arr[0] = temp[0];
    arr[1] = arr[2] = temp[1];
    arr[3] = temp[2];
    return;
  }

  if (useCramer && mat.GetNrows() <= 6){
    Bool_t (*inversion)(TMatrixDSym&, Double_t*) = 0;
    switch (mat.GetNrows()) {
    case 3:
      inversion = TMatrixTSymCramerInv::Inv3x3; break;
    case 4:
      inversion = TMatrixTSymCramerInv::Inv4x4; break;
    case 5:
      inversion = TMatrixTSymCramerInv::Inv5x5; break;
    case 6:
      inversion = TMatrixTSymCramerInv::Inv6x6; break;
    }

    Bool_t success = inversion(mat, determinant);
    if (!success){
      Exception e("Tools::invertMatrix() - cannot invert matrix, determinant = 0",
          __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }
    return;
  }

  // else use TDecompChol
  bool status = 0;
  TDecompChol invertAlgo(mat, 1E-50);

  status = invertAlgo.Invert(mat);
  if(status == 0){
    Exception e("Tools::invertMatrix() - cannot invert matrix, status = 0",
        __LINE__,__FILE__);
    e.setFatal();
    throw e;
  }

  if (determinant != NULL) {
    double d1, d2;
    invertAlgo.Det(d1, d2);
    *determinant = ldexp(d1, d2);
  }
}


// Solves R^T x = b, replaces b with the result x.  R is assumed
// to be upper-diagonal.  This is forward substitution, but with
// indices flipped.
bool tools::transposedForwardSubstitution(const TMatrixD& R, TVectorD& b)
{
  size_t n = R.GetNrows();
  for (unsigned int i = 0; i < n; ++i) {
    double sum = b(i);
    for (unsigned int j = 0; j < i; ++j) {
      sum -= b(j)*R(j,i);  // already replaced previous elements in b.
    }
    if (R(i,i) == 0)
      return false;
    b(i) = sum / R(i,i);
  }
  return true;
}

// Same, but for one column of the matrix b.  Used by transposedInvert below
bool tools::transposedForwardSubstitution(const TMatrixD& R, TMatrixD& b, int nCol)
{
  size_t n = R.GetNrows();
  for (unsigned int i = 0; i < n; ++i) {
    double sum = b(i, nCol);
    for (unsigned int j = 0; j < i; ++j) {
      sum -= b(j, nCol)*R(j,i);  // already replaced previous elements in b.
    }
    if (R(i,i) == 0)
      return false;
    b(i, nCol) = sum / R(i,i);
  }
  return true;
}

// inv will be the inverse of the transposed of the upper-right matrix R
bool tools::transposedInvert(const TMatrixD& R, TMatrixD& inv)
{
  bool result = true;

  inv.ResizeTo(R);
  for (int i = 0; i < inv.GetNrows(); ++i)
    for (int j = 0; j < inv.GetNcols(); ++j)
      inv(i, j) = (i == j);

  for (int i = 0; i < inv.GetNcols(); ++i)
    result = result && transposedForwardSubstitution(R, inv, i);

  return result;
}

// This replaces A with an upper right matrix connected to A by a
// orthogonal transformation.  I.e., it computes the R from a QR
// decomposition of A replacing A.
void tools::QR(TMatrixD& A)
{
  int nCols = A.GetNcols();
  int nRows = A.GetNrows();
  assert(nRows >= nCols);
  // This uses Businger and Golub's algorithm from Handbook for
  // Automatical Computation, Vol. 2, Chapter 8, but without
  // pivoting.  I.e., we stop at the middle of page 112.  We don't
  // explicitly calculate the orthogonal matrix.

  double *const ak = A.GetMatrixArray();
  // No variable-length arrays in C++, alloca does the exact same thing ...
  double *const u = (double *)alloca(sizeof(double)*nRows);

  // Main loop over matrix columns.
  for (int k = 0; k < nCols; ++k) {
    double akk = ak[k*nCols + k];

    double sum = akk*akk;
    // Put together a housholder transformation.
    for (int i = k + 1; i < nRows; ++i) {
      sum += ak[i*nCols + k]*ak[i*nCols + k];
      u[i] = ak[i*nCols + k];
    }
    double sigma = sqrt(sum);
    double beta = 1/(sum + sigma*fabs(akk));
    // The algorithm uses only the uk[i] for i >= k.
    u[k] = copysign(sigma + fabs(akk), akk);

    // Calculate y (again taking into account zero entries).  This
    // encodes how the (sub)matrix changes by the householder transformation.
    for (int i = k; i < nCols; ++i) {
      double y = 0;
      for (int j = k; j < nRows; ++j)
	y += u[j]*ak[j*nCols + i];
      y *= beta;
      // ... and apply the changes.
      for (int j = k; j < nRows; ++j)
	ak[j*nCols + i] -= u[j]*y; //y[j];
    }
  }

  // Zero below diagonal
  for (int i = 1; i < nCols; ++i)
    for (int j = 0; j < i; ++j)
      ak[i*nCols + j] = 0.;
  for (int i = nCols; i < nRows; ++i)
    for (int j = 0; j < nCols; ++j)
      ak[i*nCols + j] = 0.;
}

// This averages the covariance matrices C1, C2 in a numerically
// stable way by using matrix square roots.  No optimizations
// performed, so use with care.
void tools::safeAverage(const TMatrixDSym& C1, const TMatrixDSym& C2,
			TMatrixDSym& result)
{
  /*
    The algorithm proceeds as follows:
    write C1 = S1 S1' (prime for transpose),
          C2 = S2 S2'
    Then the inverse of the average can be written as ("." for matrix
    multiplication)
         C^-1 = ((S1'^-1, S2'^-1) . (S1'^-1) )
                (                   (S2'^-1) )
    Inserting an orthogonal matrix T in the middle:
         C^-1 = ((S1'^-1, S2'^-1) . T . T' . (S1'^-1) )
                (                            (S2'^-1) )
    doesn't change this because T.T' = 1.
    Now choose T s.t. T'.(S1'^-1, S2'^-1)' is an upper right matrix.  We
    use Tools::QR for the purpose, as we don't actually need T.

    Then the inverse needed to obtain the covariance matrix can be
    obtained by inverting the upper right matrix, which is squared to
    obtained the new covariance matrix.  */
  TDecompChol dec1(C1);
  dec1.Decompose();
  TDecompChol dec2(C2);
  dec2.Decompose();

  const TMatrixD& S1 = dec1.GetU();
  const TMatrixD& S2 = dec2.GetU();

  TMatrixD S1inv, S2inv;
  transposedInvert(S1, S1inv);
  transposedInvert(S2, S2inv);

  TMatrixD A(2 * S1.GetNrows(), S1.GetNcols());
  for (int i = 0; i < S1.GetNrows(); ++i) {
    for (int j = 0; j < S2.GetNcols(); ++j) {
      A(i, j) = S1inv(i, j);
      A(i + S1.GetNrows(), j) = S2inv(i, j);
    }
  }

  QR(A);
  A.ResizeTo(S1.GetNrows(), S1.GetNrows());

  TMatrixD inv;
  transposedInvert(A, inv);

  result.ResizeTo(inv.GetNcols(), inv.GetNcols());
  result = TMatrixDSym(TMatrixDSym::kAtA, inv);
}
} /* End of namespace genfit */
