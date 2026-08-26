// Tests for genfit::tools::QR and genfit::tools::transposedInvert.
//
// Each case runs the function on a fixed input and compares against hard-coded
// expected values, so any change in the numerical result shows up here.  The QR
// cases cover both of its internal loop orders, nRows <= 2*nCols and
// nRows > 2*nCols.

#include <gtest/gtest.h>

#include <Tools.h>

#include <TMatrixD.h>
#include <TVectorD.h>


namespace genfit {

namespace {

// Tolerance for the comparison against the hard-coded expected values.  The
// slack only absorbs last-bit differences between compilers; a correct
// implementation matches far more closely.
const double kTol = 1e-12;

TMatrixD makeMatrix(int nRows, int nCols, const double* v)
{
  TMatrixD M(nRows, nCols);
  for (int i = 0; i < nRows * nCols; ++i)
    M.GetMatrixArray()[i] = v[i];
  return M;
}

TVectorD makeVector(int n, const double* v)
{
  TVectorD x(n);
  for (int i = 0; i < n; ++i)
    x(i) = v[i];
  return x;
}

void expectMatrixNear(const TMatrixD& got, const double* expected, double tol)
{
  const int N = got.GetNrows() * got.GetNcols();
  for (int i = 0; i < N; ++i)
    EXPECT_NEAR(got.GetMatrixArray()[i], expected[i], tol) << "element " << i;
}

void expectVectorNear(const TVectorD& got, const double* expected, double tol)
{
  for (int i = 0; i < got.GetNrows(); ++i)
    EXPECT_NEAR(got(i), expected[i], tol) << "element " << i;
}

} // anonymous namespace


// ============================================================================
// transposedInvert -- inv is the inverse of the transpose of upper-right R.
// ============================================================================

TEST(ToolsTransposedInvert, Matches3x3Reference)
{
  const double in[9] = { 2, -1, 0.5,
                         0,  3, 1,
                         0,  0, 4 };
  // expected
  const double kInvExpected3[9] = {
     0.5,                   0,                    0,
     0.16666666666666666,   0.33333333333333331,  0,
    -0.10416666666666666,  -0.083333333333333329, 0.25,
  };

  TMatrixD R = makeMatrix(3, 3, in);
  TMatrixD inv;
  EXPECT_TRUE(tools::transposedInvert(R, inv));
  ASSERT_EQ(inv.GetNrows(), 3);
  ASSERT_EQ(inv.GetNcols(), 3);
  expectMatrixNear(inv, kInvExpected3, kTol);
}

TEST(ToolsTransposedInvert, Matches4x4Reference)
{
  const double in[16] = { 1.5, 0.2, -0.3,  0.1,
                          0,   2.0,  0.5, -0.4,
                          0,   0,    0.8,  0.25,
                          0,   0,    0,    1.25 };
  // expected
  const double kInvExpected4[16] = {
     0.66666666666666663,   0,       0,     0,
    -0.066666666666666666,  0.5,     0,     0,
     0.29166666666666663,  -0.3125,  1.25,  0,
    -0.13300000000000001,   0.2225, -0.25,  0.80000000000000004,
  };

  TMatrixD R = makeMatrix(4, 4, in);
  TMatrixD inv;
  EXPECT_TRUE(tools::transposedInvert(R, inv));
  ASSERT_EQ(inv.GetNrows(), 4);
  ASSERT_EQ(inv.GetNcols(), 4);
  expectMatrixNear(inv, kInvExpected4, kTol);
}

// A singular upper-triangular R (zero on the diagonal) must be reported as a
// failure.
TEST(ToolsTransposedInvert, RejectsSingular)
{
  const double in[9] = { 2, 1, 1,
                         0, 0, 1,    // zero pivot
                         0, 0, 3 };
  TMatrixD R = makeMatrix(3, 3, in);
  TMatrixD inv;
  EXPECT_FALSE(tools::transposedInvert(R, inv));
}


// ============================================================================
// QR -- replaces A by R (upper triangular) and b by Q'b.
// Two shapes exercise both internal loop orders of QR
// (nRows <= 2*nCols "wide", and nRows > 2*nCols "tall").
// ============================================================================

TEST(ToolsQR, MatchesWideReference)
{
  const double inA[12] = { 1, 2, 3,
                           4, 5, 6,
                           7, 8, 10,
                           1, 0, 1 };
  const double inB[4]  = { 1, 2, 3, 4 };
  // expected
  const double kRwide[12] = {
    -8.1853527718724521, -9.5292166597918087, -11.972605546917913,
     0,                   1.4812257933030573,   1.2897748404271523,
     0,                   0,                    0.99659283506934926,
     0,                   0,                    0,
  };
  const double kQtbWide[4] = {
    -4.1537611081143782, -2.4183278258009064, 2.3959183911598769, -1.0758876551830792,
  };

  TMatrixD A = makeMatrix(4, 3, inA);
  TVectorD b = makeVector(4, inB);
  tools::QR(A, b);
  expectMatrixNear(A, kRwide, kTol);
  expectVectorNear(b, kQtbWide, kTol);

  // R must be upper triangular (exact zeros below the diagonal).
  for (int i = 0; i < A.GetNrows(); ++i)
    for (int j = 0; j < A.GetNcols(); ++j)
      if (i > j) {
        EXPECT_DOUBLE_EQ(A(i, j), 0.0);
      }
}

TEST(ToolsQR, MatchesTallReference)
{
  const double inA[14] = { 1, 2,
                           3, 4,
                           5, 6,
                           7, 8,
                           2, 1,
                           0, 3,
                           1, 1 };
  const double inB[7]  = { 1, 2, 3, 4, 5, 6, 7 };
  // expected
  const double kRtall[14] = {
    -9.4339811320566032, -10.917978164065509,
     0,                   -3.4347857005916342,
     0,                    0,
     0,                    0,
     0,                    0,
     0,                    0,
     0,                    0,
  };
  const double kQtbTall[7] = {
    -7.1019857960426123, -3.627787944720116, -0.71376370727510574,
    -0.84467000412574733, 4.7555916682978587, 3.3405061045544975, 6.4345468515746784,
  };

  TMatrixD A = makeMatrix(7, 2, inA);
  TVectorD b = makeVector(7, inB);
  tools::QR(A, b);
  expectMatrixNear(A, kRtall, kTol);
  expectVectorNear(b, kQtbTall, kTol);

  for (int i = 0; i < A.GetNrows(); ++i)
    for (int j = 0; j < A.GetNcols(); ++j)
      if (i > j) {
        EXPECT_DOUBLE_EQ(A(i, j), 0.0);
      }
}

} // namespace genfit
