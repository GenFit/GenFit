// Tests for genfit::tools::averageState, the numerical kernel of
// genfit::calcAverageState().
//
// The checks are that, across a range of covariance dimensions, averageState
//   * agrees with a reference built on ROOT's TDecompChol,
//   * reproduces the textbook inverse-variance weighted average,
// and that a non-positive-definite covariance is rejected.

#include <gtest/gtest.h>

#include <Tools.h>

#include <TDecompChol.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>

#include <cmath>
#include <random>

namespace genfit {

namespace {

// Reference implementation of the same average, driven by ROOT's TDecompChol.
// Cross-check for tools::averageState().
bool refROOT(const TVectorD& x1, const TMatrixDSym& C1,
             const TVectorD& x2, const TMatrixDSym& C2,
             TVectorD& avgState, TMatrixD& inv)
{
  TDecompChol d1(C1);
  bool success = d1.Decompose();
  TDecompChol d2(C2);
  success &= d2.Decompose();
  if (!success)
    return false;

  const int nRows = d1.GetU().GetNrows();
  TMatrixD S1inv, S2inv;
  tools::transposedInvert(d1.GetU(), S1inv);
  tools::transposedInvert(d2.GetU(), S2inv);

  TMatrixD A(2*nRows, nRows);
  TVectorD b(2*nRows);
  double *const bData = b.GetMatrixArray();
  double *const AData = A.GetMatrixArray();
  const double* S1invData = S1inv.GetMatrixArray();
  const double* S2invData = S2inv.GetMatrixArray();
  for (int i = 0; i < nRows; ++i) {
    double sum1 = 0;
    double sum2 = 0;
    for (int j = 0; j <= i; ++j) {
      AData[i*nRows + j] = S1invData[i*nRows + j];
      AData[(i + nRows)*nRows + j] = S2invData[i*nRows + j];
      sum1 += S1invData[i*nRows + j]*x1.GetMatrixArray()[j];
      sum2 += S2invData[i*nRows + j]*x2.GetMatrixArray()[j];
    }
    bData[i] = sum1;
    bData[i + nRows] = sum2;
  }

  tools::QR(A, b);
  A.ResizeTo(nRows, nRows);

  tools::transposedInvert(A, inv);
  b.ResizeTo(nRows);
  for (int i = 0; i < nRows; ++i) {
    double sum = 0;
    for (int j = i; j < nRows; ++j)
      sum += inv.GetMatrixArray()[j*nRows+i] * b[j];
    b.GetMatrixArray()[i] = sum;
  }

  avgState.ResizeTo(nRows);
  avgState = b;
  return true;
}

// Random symmetric positive definite covariance C = M^T M + (n+1) I, and a
// random state x.
void makeSymPosDef(TMatrixDSym& C, TVectorD& x, int n, std::mt19937_64& rng)
{
  std::uniform_real_distribution<double> d(-1.0, 1.0);
  TMatrixD M(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      M(i, j) = d(rng);
  C.ResizeTo(n, n);
  C = TMatrixDSym(TMatrixDSym::kAtA, M);
  for (int i = 0; i < n; ++i)
    C(i, i) += (n + 1.0);
  x.ResizeTo(n);
  for (int i = 0; i < n; ++i)
    x(i) = d(rng);
}

} // anonymous namespace

// Agreement with a reference that Cholesky-factors the covariances with ROOT's
// TDecompChol.  Both paths share tools::QR and tools::transposedInvert, which
// TestTools.cpp pins separately, so this isolates the inline factorisation.
// It follows TDecompChol's operation order, so the averaged state and the
// covariance factor have to come out exactly equal, not merely close.
TEST(AverageState, AgreesWithRootTDecompChol)
{
  std::mt19937_64 rng(20240629);

  for (int n : {1, 2, 3, 4, 5, 6, 7, 10, 15, 25}) {
    for (int t = 0; t < 200; ++t) {
      TMatrixDSym C1, C2; TVectorD x1, x2;
      makeSymPosDef(C1, x1, n, rng);
      makeSymPosDef(C2, x2, n, rng);

      TVectorD sRef, sOpt;
      TMatrixD factorRef, factorOpt;
      ASSERT_TRUE(refROOT(x1, C1, x2, C2, sRef, factorRef));
      ASSERT_TRUE(tools::averageState(x1, C1, x2, C2, sOpt, factorOpt));

      for (int i = 0; i < n; ++i)
        ASSERT_EQ(sOpt(i), sRef(i)) << "state entry " << i << ", n = " << n;

      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          ASSERT_EQ(factorOpt(i, j), factorRef(i, j))
              << "covariance factor entry (" << i << ", " << j << "), n = " << n;
    }
  }
}

// Reproduces the textbook inverse-variance weighted average, checked with
// ROOT's dense inversion:  avgCov == (C1^-1 + C2^-1)^-1  and
//                          avgState == avgCov (C1^-1 x1 + C2^-1 x2).
TEST(AverageState, MatchesTextbookAverage)
{
  std::mt19937_64 rng(987654321);
  const double tol = 1e-8;
  double worst = 0.0;

  for (int n : {1, 2, 3, 5, 6, 7, 12}) {
    for (int t = 0; t < 200; ++t) {
      TMatrixDSym C1, C2; TVectorD x1, x2;
      makeSymPosDef(C1, x1, n, rng);
      makeSymPosDef(C2, x2, n, rng);

      TVectorD sOpt; TMatrixD covFactor;
      ASSERT_TRUE(tools::averageState(x1, C1, x2, C2, sOpt, covFactor));

      TMatrixDSym C1inv(C1); C1inv.Invert();
      TMatrixDSym C2inv(C2); C2inv.Invert();
      TMatrixDSym sumInv(C1inv); sumInv += C2inv;
      TMatrixDSym avgCovExp(sumInv); avgCovExp.Invert();
      TVectorD stateExp = avgCovExp * ((C1inv * x1) + (C2inv * x2));

      TMatrixDSym avgCov(TMatrixDSym::kAtA, covFactor);

      double scaleC = 1e-300;
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          scaleC = std::max(scaleC, std::fabs(avgCovExp(i, j)));
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          worst = std::max(worst, std::fabs(avgCov(i, j) - avgCovExp(i, j)) / scaleC);

      double scaleS = 1e-300;
      for (int i = 0; i < n; ++i) scaleS = std::max(scaleS, std::fabs(stateExp(i)));
      for (int i = 0; i < n; ++i)
        worst = std::max(worst, std::fabs(sOpt(i) - stateExp(i)) / scaleS);
    }
  }

  EXPECT_LT(worst, tol) << "worst relative averaging residual";
}

// A covariance that is not positive definite must be rejected (return false).
TEST(AverageState, RejectsNonPositiveDefinite)
{
  std::mt19937_64 rng(13);
  const int n = 5;
  TMatrixDSym C1, C2; TVectorD x1, x2;
  makeSymPosDef(C1, x1, n, rng);
  makeSymPosDef(C2, x2, n, rng);
  C1(2, 2) = -1.0;   // break positive-definiteness

  TVectorD s; TMatrixD cov;
  EXPECT_FALSE(tools::averageState(x1, C1, x2, C2, s, cov));
}

} // namespace genfit
