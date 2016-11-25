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
#include <TMatrixDSymEigen.h>
#include <TMatrixTSymCramerInv.h>
#include <TMath.h>

#include "AbsHMatrix.h"
#include "Exception.h"


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
    if (determinant != nullptr) *determinant = mat(0,0);
    inv(0,0) = 1./mat(0,0);
    return;
  }

  if (mat.GetNrows() == 2){
    double det = mat(0,0)*mat(1,1) - mat(1,0)*mat(1,0);
    if (determinant != nullptr) *determinant = det;
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

  if (determinant != nullptr) {
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
    if (determinant != nullptr) *determinant = mat(0,0);
    mat(0,0) = 1./mat(0,0);
    return;
  }

  if (mat.GetNrows() == 2){
    double *arr = mat.GetMatrixArray();
    double det = arr[0]*arr[3] - arr[1]*arr[1];
    if (determinant != nullptr) *determinant = det;
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

  if (determinant != nullptr) {
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
  double *const bk = b.GetMatrixArray();
  const double *const Rk = R.GetMatrixArray();
  for (unsigned int i = 0; i < n; ++i) {
    double sum = bk[i];
    for (unsigned int j = 0; j < i; ++j) {
      sum -= bk[j]*Rk[j*n + i];  // already replaced previous elements in b.
    }
    if (Rk[i*n+i] == 0)
      return false;
    bk[i] = sum / Rk[i*n + i];
  }
  return true;
}


// Same, but for one column of the matrix b.  Used by transposedInvert below
// assumes b(i,j) == (i == j)
bool tools::transposedForwardSubstitution(const TMatrixD& R, TMatrixD& b, int nCol)
{
  size_t n = R.GetNrows();
  double *const bk = b.GetMatrixArray() + nCol;
  const double *const Rk = R.GetMatrixArray();
  for (unsigned int i = nCol; i < n; ++i) {
    double sum = (i == (size_t)nCol);
    for (unsigned int j = 0; j < i; ++j) {
      sum -= bk[j*n]*Rk[j*n + i];  // already replaced previous elements in b.
    }
    if (Rk[i*n+i] == 0)
      return false;
    bk[i*n] = sum / Rk[i*n + i];
  }
  return true;
}


// inv will be the inverse of the transposed of the upper-right matrix R
bool tools::transposedInvert(const TMatrixD& R, TMatrixD& inv)
{
  bool result = true;

  inv.ResizeTo(R);
  double *const invk = inv.GetMatrixArray();
  int nRows = inv.GetNrows();
  for (int i = 0; i < nRows; ++i)
    for (int j = 0; j < nRows; ++j)
      invk[i*nRows + j] = (i == j);

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

// This replaces A with an upper right matrix connected to A by a
// orthogonal transformation.  I.e., it computes the R from a QR
// decomposition of A replacing A.  Simultaneously it transforms b by
// the inverse orthogonal transformation.
// 
// The purpose is this: the least-squared problem
//   ||Ax - b|| = min
// is equivalent to
//   ||QRx - b|| = ||Rx - Q'b|| = min
// where Q' denotes the transposed (i.e. inverse).
void tools::QR(TMatrixD& A, TVectorD& b)
{
  int nCols = A.GetNcols();
  int nRows = A.GetNrows();
  assert(nRows >= nCols);
  assert(b.GetNrows() == nRows);
  // This uses Businger and Golub's algorithm from Handbook for
  // Automatic Computation, Vol. 2, Chapter 8, but without pivoting.
  // I.e., we stop at the middle of page 112.  We don't explicitly
  // calculate the orthogonal matrix, but Q'b which is not done
  // explicitly in Businger et al.
  // Also in Numer. Math. 7, 269-276 (1965)

  double *const ak = A.GetMatrixArray();
  double *const bk = b.GetMatrixArray();
  // No variable-length arrays in C++, alloca does the exact same thing ...
  //double * u = (double *)alloca(sizeof(double)*nRows);
  double u[500];

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

    // Calculate b (again taking into account zero entries).  This
    // encodes how the (sub)vector changes by the householder transformation.
    double yb = 0;
    for (int j = k; j < nRows; ++j)
      yb += u[j]*bk[j];
    yb *= beta;
    // ... and apply the changes.
    for (int j = k; j < nRows; ++j)
      bk[j] -= u[j]*yb;

    // Calculate y (again taking into account zero entries).  This
    // encodes how the (sub)matrix changes by the householder transformation.
    for (int i = k; i < nCols; ++i) {
      double y = 0;
      for (int j = k; j < nRows; ++j)
	y += u[j]*ak[j*nCols + i];
      y *= beta;
      // ... and apply the changes.
      for (int j = k; j < nRows; ++j)
	ak[j*nCols + i] -= u[j]*y;
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


void
tools::noiseMatrixSqrt(const TMatrixDSym& noise,
		       TMatrixD& noiseSqrt)
{
  // This is the slowest part of the whole Sqrt Kalman.  Using an LDLt
  // transform is probably the easiest way of remedying this.
  TMatrixDSymEigen eig(noise);
  noiseSqrt.ResizeTo(noise);
  noiseSqrt = eig.GetEigenVectors();
  double* pNoiseSqrt = noiseSqrt.GetMatrixArray();
  const TVectorD& evs(eig.GetEigenValues());
  const double* pEvs = evs.GetMatrixArray();
  // GetEigenVectors is such that noise = noiseSqrt * evs * noiseSqrt'
  // We're evaluating the first product with the eigenvalues replaced
  // by their square roots, so we're multiplying with a diagonal
  // matrix from the right.
  int iCol = 0;
  for (; iCol < noiseSqrt.GetNrows(); ++iCol) {
    double ev = pEvs[iCol] > 0 ? sqrt(pEvs[iCol]) : 0;
    // if (ev == 0)
    //  break;
    for (int j = 0; j < noiseSqrt.GetNrows(); ++j) {
      pNoiseSqrt[j*noiseSqrt.GetNcols() + iCol] *= ev;
    }
  }
  if (iCol < noiseSqrt.GetNcols()) {
    // Hit zero eigenvalue, resize matrix
    noiseSqrt.ResizeTo(noiseSqrt.GetNrows(), iCol);
  }

  // noiseSqrt * noiseSqrt' = noise
}


// Transports the square root of the covariance matrix using a
// square-root formalism
//
// With covariance square root S, transport matrix F and noise matrix
// square root Q.
void
tools::kalmanPredictionCovSqrt(const TMatrixD& S,
			       const TMatrixD& F, const TMatrixD& Q,
			       TMatrixD& Snew)
{
  Snew.ResizeTo(S.GetNrows() + Q.GetNcols(),
		S.GetNcols());

  // This overwrites all elements, no precautions necessary
  Snew.SetSub(0, 0, TMatrixD(S, TMatrixD::kMultTranspose, F));
  if (Q.GetNcols() != 0)
    Snew.SetSub(S.GetNrows(), 0, TMatrixD(TMatrixD::kTransposed, Q));

  tools::QR(Snew);

  // The result is in the upper right corner of the matrix.
  Snew.ResizeTo(S.GetNrows(), S.GetNrows());
}


// Kalman measurement update (no transport)
// x, S : state prediction, covariance square root
// res, R, H : residual, measurement covariance square root, H matrix of the measurement
// gives the update (new state = x + update) and the updated covariance square root.
// S and Snew are allowed to refer to the same object.
void
tools::kalmanUpdateSqrt(const TMatrixD& S,
			const TVectorD& res, const TMatrixD& R,
			const AbsHMatrix* H,
			TVectorD& update, TMatrixD& SNew)
{
  TMatrixD pre(S.GetNrows() + R.GetNrows(),
	       S.GetNcols() + R.GetNcols());
  pre.SetSub(0,            0,         R); /* Zeros in upper right block */
  pre.SetSub(R.GetNrows(), 0, H->MHt(S)); pre.SetSub(R.GetNrows(), R.GetNcols(), S);

  tools::QR(pre);
  const TMatrixD& r = pre;

  const TMatrixD& a(r.GetSub(0, R.GetNrows()-1,
			     0, R.GetNcols()-1));
  TMatrixD K(TMatrixD::kTransposed, r.GetSub(0, R.GetNrows()-1, R.GetNcols(), pre.GetNcols()-1));
  SNew = r.GetSub(R.GetNrows(), pre.GetNrows()-1, R.GetNcols(), pre.GetNcols()-1);

  update.ResizeTo(res);
  update = res;
  tools::transposedForwardSubstitution(a, update);
  update *= K;
}

} /* End of namespace genfit */
