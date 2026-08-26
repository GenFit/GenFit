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
#include <cstring>
#include <algorithm>
#include <vector>

#include <TDecompChol.h>
#include <TMatrixDSymEigen.h>
#include <TMatrixTSymCramerInv.h>
#include <TMath.h>

#include "AbsHMatrix.h"
#include "Exception.h"


namespace {

  /** Largest dimension for which a similarity result is assembled on the stack. */
  const int c_maxSimilarityDim = 7;

  /** Similarity transform out = b * sym * b^T for dimensions known at compile
   * time, with b of size N x M and sym of size M x M.
   *
   * The two matrix products are accumulated in the order TMatrixDSym::Similarity()
   * uses, and as it does only the upper triangle of the result is computed and
   * then mirrored, so the result is bit for bit the one ROOT produces.
   */
  template<int N, int M>
  void similarityFixed(const double* b, const double* sym, double* out)
  {
    double ba[N * M];
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < M; ++j) {
        double element = 0;
        for (int k = 0; k < M; ++k) element += b[i * M + k] * sym[k * M + j];
        ba[i * M + j] = element;
      }
    }
    for (int i = 0; i < N; ++i) {
      for (int j = i; j < N; ++j) {
        double element = 0;
        for (int k = 0; k < M; ++k) element += ba[i * M + k] * b[j * M + k];
        out[i * N + j] = element;
        out[j * N + i] = element;
      }
    }
  }

  /** Same, for dimensions only known at run time. */
  void similarityDynamic(const double* b, int nRows, int nCols, const double* sym, double* out)
  {
    std::vector<double> ba(nRows * nCols);
    for (int i = 0; i < nRows; ++i) {
      for (int j = 0; j < nCols; ++j) {
        double element = 0;
        for (int k = 0; k < nCols; ++k) element += b[i * nCols + k] * sym[k * nCols + j];
        ba[i * nCols + j] = element;
      }
    }
    for (int i = 0; i < nRows; ++i) {
      for (int j = i; j < nRows; ++j) {
        double element = 0;
        for (int k = 0; k < nCols; ++k) element += ba[i * nCols + k] * b[j * nCols + k];
        out[i * nRows + j] = element;
        out[j * nRows + i] = element;
      }
    }
  }

  /** Pick the implementation matching the given dimensions.  The shapes that
   * occur in the Kalman fitters and the track representations are worth
   * unrolling; everything else falls back to the run-time loops.
   */
  void similarityInto(const double* b, int nRows, int nCols, const double* sym, double* out)
  {
    if (nRows == 5 && nCols == 5) similarityFixed<5, 5>(b, sym, out);
    else if (nRows == 7 && nCols == 7) similarityFixed<7, 7>(b, sym, out);
    else if (nRows == 5 && nCols == 1) similarityFixed<5, 1>(b, sym, out);
    else if (nRows == 5 && nCols == 2) similarityFixed<5, 2>(b, sym, out);
    else similarityDynamic(b, nRows, nCols, sym, out);
  }

} // anonymous namespace

namespace genfit {

void tools::similarity(const TMatrixD& b, const TMatrixDSym& sym, TMatrixDSym& out)
{
  const int nRows = b.GetNrows();
  const int nCols = b.GetNcols();

  if (nCols != sym.GetNrows()) {
    Exception e("Tools::similarity() - matrix dimensions do not match", __LINE__, __FILE__);
    e.setFatal();
    throw e;
  }

  out.ResizeTo(nRows, nRows);
  similarityInto(b.GetMatrixArray(), nRows, nCols, sym.GetMatrixArray(), out.GetMatrixArray());
}


void tools::similarity(const TMatrixD& b, TMatrixDSym& sym)
{
  const int nRows = b.GetNrows();
  const int nCols = b.GetNcols();

  if (nCols != sym.GetNrows()) {
    Exception e("Tools::similarity() - matrix dimensions do not match", __LINE__, __FILE__);
    e.setFatal();
    throw e;
  }

  // The result is n x n while sym is m x m, so it is assembled elsewhere before
  // sym can be overwritten.
  if (nRows <= c_maxSimilarityDim) {
    double result[c_maxSimilarityDim * c_maxSimilarityDim];
    similarityInto(b.GetMatrixArray(), nRows, nCols, sym.GetMatrixArray(), result);
    sym.ResizeTo(nRows, nRows);
    std::copy(result, result + nRows * nRows, sym.GetMatrixArray());
  } else {
    TMatrixDSym result;
    tools::similarity(b, sym, result);
    sym.ResizeTo(result);
    sym = result;
  }
}


double tools::similarity(const TVectorD& v, const TMatrixDSym& sym)
{
  const int n = v.GetNrows();

  if (n != sym.GetNrows()) {
    Exception e("Tools::similarity() - vector and matrix dimensions do not match", __LINE__, __FILE__);
    e.setFatal();
    throw e;
  }

  const double* const vArray = v.GetMatrixArray();
  const double* const symArray = sym.GetMatrixArray();

  // Order of operations as in TMatrixDSym::Similarity(const TVectorD&): first
  // sym * v, then the scalar product with v.
  std::vector<double> heapWork;
  double stackWork[c_maxSimilarityDim];
  double* work = stackWork;
  if (n > c_maxSimilarityDim) {
    heapWork.resize(n);
    work = heapWork.data();
  }
  for (int i = 0; i < n; ++i) {
    double element = 0;
    for (int k = 0; k < n; ++k) element += symArray[i * n + k] * vArray[k];
    work[i] = element;
  }

  double sum = 0;
  for (int i = 0; i < n; ++i) sum += work[i] * vArray[i];
  return sum;
}


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
  inv.ResizeTo(R);

  const int n = R.GetNrows();
  const double* const Rk = R.GetMatrixArray();
  double* const invk = inv.GetMatrixArray();

  // inv is lower triangular (it is the transpose of the inverse of the
  // upper-right R).  Start from the identity: the substitution below overwrites
  // the lower triangle including the diagonal, so the strict upper triangle
  // keeps its zeros, and on a singular R the entries never reached stay as the
  // identity.
  std::memset(invk, 0, sizeof(double) * n * n);
  for (int i = 0; i < n; ++i)
    invk[i*n + i] = 1.;

  // Forward substitution solving R^T inv = I, one column at a time.  Column c
  // vanishes above row c, and within row i only j in [c, i) contributes, so
  // both loops start at c and no structural zero is ever multiplied in.
  bool result = true;
  for (int c = 0; c < n && result; ++c) {
    const double* const invc = invk + c;        // inv(j,c) == invc[j*n]
    for (int i = c; i < n; ++i) {
      const double* const Ri = Rk + i;          // R(j,i) == Ri[j*n]
      double sum = (i == c);
      for (int j = c; j < i; ++j)
        sum -= invc[j*n] * Ri[j*n];
      const double Rii = Rk[i*n + i];
      if (Rii == 0) {
        result = false;
        break;
      }
      invk[i*n + c] = sum / Rii;
    }
  }

  return result;
}

namespace {

// Upper-triangular Cholesky factorisation of a symmetric positive-definite
// matrix C (n x n, full row-major storage as kept by TMatrixDSym), writing the
// upper-triangular factor U (with C = U^T U) into the upper triangle of the
// n x n buffer U.  Only the upper triangle and the diagonal of U are written;
// the strictly-lower triangle is left as-is (the sole consumer here,
// transposedInvert(), reads only the upper triangle).
//
// The summation order, and the per-row reciprocal 1/ujj that scales the
// off-diagonal entries, are those of ROOT's TDecompChol::Decompose(), so U is
// bit-for-bit the factor TDecompChol produces.  Returns false if C is not
// positive definite.
bool choleskyUpper(const double* C, double* U, int n)
{
  for (int icol = 0; icol < n; ++icol) {
    const int rowOff = icol*n;

    // Diagonal element U(icol,icol); test for non-positive-definiteness.
    double ujj = C[rowOff + icol];
    for (int irow = 0; irow < icol; ++irow) {
      const double u = U[irow*n + icol];
      ujj -= u*u;
    }
    if (ujj <= 0)
      return false;
    ujj = sqrt(ujj);
    U[rowOff + icol] = ujj;
    const double inv_ujj = 1.0 / ujj;

    // Off-diagonal elements of this row.
    for (int j = icol + 1; j < n; ++j) {
      double s = C[rowOff + j];
      for (int i = 0; i < icol; ++i)
        s -= U[i*n + j]*U[i*n + icol];
      U[rowOff + j] = s * inv_ujj;
    }
  }
  return true;
}

} // anonymous namespace

bool tools::averageState(const TVectorD& state1, const TMatrixDSym& cov1,
                         const TVectorD& state2, const TMatrixDSym& cov2,
                         TVectorD& avgState, TMatrixD& avgCovFactor)
{
  // See genfit::calcAverageState() for the derivation.  In short: with the
  // upper Cholesky factors S1, S2 (cov1 = S1' S1, cov2 = S2' S2) the combined
  // information is (S1inv', S2inv').(S1inv; S2inv) where Sinv = (S')^-1.  A QR
  // decomposition of A = (S1inv; S2inv) gives an upper triangular R with
  // R'R = cov1^-1 + cov2^-1, hence avgCov = R^-1 R'^-1 and the averaged state
  // follows from the same orthogonal transformation applied to (S1inv.x1;
  // S2inv.x2).
  const int nRows = cov1.GetNrows();
  assert(cov2.GetNrows() == nRows);

  // Upper Cholesky factors S1, S2 of the two covariance matrices.
  TMatrixD S1(nRows, nRows), S2(nRows, nRows);
  if (!choleskyUpper(cov1.GetMatrixArray(), S1.GetMatrixArray(), nRows))
    return false;
  if (!choleskyUpper(cov2.GetMatrixArray(), S2.GetMatrixArray(), nRows))
    return false;

  // S1inv = (S1')^-1, S2inv = (S2')^-1 -- both lower triangular.
  TMatrixD S1inv, S2inv;
  transposedInvert(S1, S1inv);
  transposedInvert(S2, S2inv);

  // Assemble A = (S1inv; S2inv) and b = (S1inv.x1; S2inv.x2).  A is zero-filled
  // by its constructor, so the strictly-upper-triangular halves of the two
  // lower-triangular blocks stay zero and only j <= i is written.
  TMatrixD A(2*nRows, nRows);
  TVectorD b(2*nRows);
  double *const Ak = A.GetMatrixArray();
  double *const bk = b.GetMatrixArray();
  const double* const S1invk = S1inv.GetMatrixArray();
  const double* const S2invk = S2inv.GetMatrixArray();
  const double* const x1 = state1.GetMatrixArray();
  const double* const x2 = state2.GetMatrixArray();
  for (int i = 0; i < nRows; ++i) {
    double sum1 = 0;
    double sum2 = 0;
    for (int j = 0; j <= i; ++j) {
      const double s1 = S1invk[i*nRows + j];
      const double s2 = S2invk[i*nRows + j];
      Ak[i*nRows + j]           = s1;
      Ak[(i + nRows)*nRows + j] = s2;
      sum1 += s1*x1[j];
      sum2 += s2*x2[j];
    }
    bk[i]         = sum1;
    bk[i + nRows] = sum2;
  }

  // QR decomposition: R (upper triangular) ends up in the top nRows rows of A,
  // and Q'b in b (only its first nRows entries are needed below).
  QR(A, b);
  A.ResizeTo(nRows, nRows);

  // avgCovFactor = (R')^-1 (lower triangular); avgCov = avgCovFactor' avgCovFactor.
  transposedInvert(A, avgCovFactor);

  // Averaged state = R^-1 . (top of Q'b), i.e. avgCovFactor' . (Q'b).
  const double* const invk = avgCovFactor.GetMatrixArray();
  avgState.ResizeTo(nRows);
  double* const sk = avgState.GetMatrixArray();
  for (int i = 0; i < nRows; ++i) {
    double sum = 0;
    for (int j = i; j < nRows; ++j)
      sum += invk[j*nRows + i] * bk[j];
    sk[i] = sum;
  }

  return true;
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
  // No variable-length arrays in C++; alloca gives the nRows-long scratch
  // vector holding the Householder reflector.
  double *const u = (double *)alloca(sizeof(double)*nRows);

  // The Householder update of the trailing submatrix, (I - beta u u') A, is
  // applied in one of two loop orders.  The rank-1 form runs its inner loop
  // contiguously across the remaining columns, which is cache friendly and
  // vectorises, but it needs enough columns to pay off; the column-by-column
  // form strides down the rows and pipelines better on tall-skinny shapes.
  // The two orders accumulate in the same sequence and give identical results.
  const bool wide = (nRows <= 2*nCols);
  double *const w = wide ? (double *)alloca(sizeof(double)*nCols) : nullptr;

  // Main loop over matrix columns.
  for (int k = 0; k < nCols; ++k) {
    double akk = ak[k*nCols + k];

    double sum = akk*akk;
    // Put together a housholder transformation.
    for (int i = k + 1; i < nRows; ++i) {
      double aik = ak[i*nCols + k];
      sum += aik*aik;
      u[i] = aik;
    }
    double sigma = sqrt(sum);
    double beta = 1/(sum + sigma*fabs(akk));
    // The algorithm uses only the uk[i] for i >= k.
    u[k] = copysign(sigma + fabs(akk), akk);

    if (wide) {
      // First pass: accumulate w[i] = sum_j u[j] A[j][i] for the trailing
      // columns (contiguous inner loop) and yb = sum_j u[j] b[j].
      for (int i = k; i < nCols; ++i) w[i] = 0;
      double yb = 0;
      for (int j = k; j < nRows; ++j) {
	const double uj = u[j];
	const double* arow = ak + j*nCols;
	for (int i = k; i < nCols; ++i)
	  w[i] += uj*arow[i];
	yb += uj*bk[j];
      }
      for (int i = k; i < nCols; ++i) w[i] *= beta;
      yb *= beta;
      // Second pass: rank-1 update A[j][i] -= u[j] w[i] and b[j] -= u[j] yb.
      for (int j = k; j < nRows; ++j) {
	const double uj = u[j];
	double* arow = ak + j*nCols;
	for (int i = k; i < nCols; ++i)
	  arow[i] -= uj*w[i];
	bk[j] -= uj*yb;
      }
    } else {
      // Tall-skinny: apply the transformation to b, then column by column to A.
      double yb = 0;
      for (int j = k; j < nRows; ++j)
	yb += u[j]*bk[j];
      yb *= beta;
      for (int j = k; j < nRows; ++j)
	bk[j] -= u[j]*yb;

      for (int i = k; i < nCols; ++i) {
	double y = 0;
	for (int j = k; j < nRows; ++j)
	  y += u[j]*ak[j*nCols + i];
	y *= beta;
	for (int j = k; j < nRows; ++j)
	  ak[j*nCols + i] -= u[j]*y;
      }
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
