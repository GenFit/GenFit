// Tests for the genfit::tools::similarity overloads.
//
// They exist to give the same result as TMatrixDSym::Similarity() with less
// work, so every case here compares against ROOT and requires the results to be
// exactly equal, not merely close.  The shapes cover the unrolled dimensions as
// well as the run-time fallback, and for the in-place overload both the stack
// and the heap path.

#include <gtest/gtest.h>

#include <Exception.h>
#include <Tools.h>

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>

#include <random>

namespace genfit {

namespace {

// Shapes the unrolled implementations are selected for, plus shapes that fall
// through to the run-time loops, including nRows < nCols and a dimension above
// the stack-buffer limit.
const int kShapes[][2] = {
  {5, 5}, {7, 7}, {5, 1}, {5, 2},                    // unrolled
  {1, 1}, {2, 3}, {3, 2}, {4, 4}, {6, 6}, {8, 3}, {3, 8}, {10, 10}  // dynamic
};

TMatrixD randomMatrix(int nRows, int nCols, std::mt19937_64& rng)
{
  std::uniform_real_distribution<double> d(-3., 3.);
  TMatrixD M(nRows, nCols);
  for (int i = 0; i < nRows; ++i)
    for (int j = 0; j < nCols; ++j)
      M(i, j) = d(rng);
  return M;
}

TMatrixDSym randomSym(int n, std::mt19937_64& rng)
{
  std::uniform_real_distribution<double> d(-3., 3.);
  TMatrixDSym S(n);
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      const double v = d(rng);
      S(i, j) = v;
      S(j, i) = v;
    }
  }
  return S;
}

TVectorD randomVector(int n, std::mt19937_64& rng)
{
  std::uniform_real_distribution<double> d(-3., 3.);
  TVectorD v(n);
  for (int i = 0; i < n; ++i)
    v(i) = d(rng);
  return v;
}

void expectEqualToRoot(const TMatrixDSym& got, const TMatrixD& b, const TMatrixDSym& sym)
{
  TMatrixDSym expected(sym);
  expected.Similarity(b);

  ASSERT_EQ(got.GetNrows(), expected.GetNrows());
  ASSERT_EQ(got.GetNcols(), expected.GetNcols());
  for (int i = 0; i < expected.GetNrows(); ++i) {
    for (int j = 0; j < expected.GetNcols(); ++j) {
      ASSERT_EQ(got(i, j), expected(i, j))
          << "entry (" << i << ", " << j << ") for "
          << b.GetNrows() << "x" << b.GetNcols();
    }
  }
}

} // anonymous namespace


// ============================================================================
// similarity(b, sym, out) -- out = b * sym * b^T.
// ============================================================================

TEST(ToolsSimilarity, MatchesRootForEveryShape)
{
  std::mt19937_64 rng(4711);

  for (const auto& shape : kShapes) {
    const int nRows = shape[0];
    const int nCols = shape[1];
    for (int t = 0; t < 50; ++t) {
      const TMatrixD b = randomMatrix(nRows, nCols, rng);
      const TMatrixDSym sym = randomSym(nCols, rng);

      TMatrixDSym out;
      tools::similarity(b, sym, out);
      expectEqualToRoot(out, b, sym);
    }
  }
}

// The result is symmetric by construction: only the upper triangle is computed
// and then mirrored.
TEST(ToolsSimilarity, ResultIsExactlySymmetric)
{
  std::mt19937_64 rng(1234);

  for (const auto& shape : kShapes) {
    const TMatrixD b = randomMatrix(shape[0], shape[1], rng);
    const TMatrixDSym sym = randomSym(shape[1], rng);

    TMatrixDSym out;
    tools::similarity(b, sym, out);

    for (int i = 0; i < out.GetNrows(); ++i)
      for (int j = 0; j < out.GetNcols(); ++j)
        ASSERT_EQ(out(i, j), out(j, i)) << "entry (" << i << ", " << j << ")";
  }
}

// An output matrix that already has the wrong size is resized.
TEST(ToolsSimilarity, ResizesOutput)
{
  std::mt19937_64 rng(99);
  const TMatrixD b = randomMatrix(5, 2, rng);
  const TMatrixDSym sym = randomSym(2, rng);

  TMatrixDSym out(3);
  tools::similarity(b, sym, out);

  EXPECT_EQ(out.GetNrows(), 5);
  expectEqualToRoot(out, b, sym);
}

TEST(ToolsSimilarity, RejectsMismatchedDimensions)
{
  std::mt19937_64 rng(7);
  const TMatrixD b = randomMatrix(5, 3, rng);
  const TMatrixDSym sym = randomSym(4, rng);

  TMatrixDSym out;
  EXPECT_THROW(tools::similarity(b, sym, out), Exception);
}


// ============================================================================
// similarity(b, sym) -- in place, sym = b * sym * b^T.
// ============================================================================

TEST(ToolsSimilarityInPlace, MatchesRootForEveryShape)
{
  std::mt19937_64 rng(20250826);

  for (const auto& shape : kShapes) {
    const int nRows = shape[0];
    const int nCols = shape[1];
    for (int t = 0; t < 50; ++t) {
      const TMatrixD b = randomMatrix(nRows, nCols, rng);
      const TMatrixDSym sym = randomSym(nCols, rng);

      TMatrixDSym inPlace(sym);
      tools::similarity(b, inPlace);

      ASSERT_EQ(inPlace.GetNrows(), nRows)
          << "in-place result not resized for " << nRows << "x" << nCols;
      expectEqualToRoot(inPlace, b, sym);
    }
  }
}

TEST(ToolsSimilarityInPlace, RejectsMismatchedDimensions)
{
  std::mt19937_64 rng(21);
  const TMatrixD b = randomMatrix(5, 3, rng);
  TMatrixDSym sym = randomSym(4, rng);

  EXPECT_THROW(tools::similarity(b, sym), Exception);
}


// ============================================================================
// similarity(v, sym) -- the scalar v^T * sym * v.
// ============================================================================

TEST(ToolsSimilarityVector, MatchesRoot)
{
  std::mt19937_64 rng(555);

  // Covers both the stack buffer and, for n above its limit, the heap one.
  for (int n : {1, 2, 3, 5, 6, 7, 8, 12, 25}) {
    for (int t = 0; t < 50; ++t) {
      const TVectorD v = randomVector(n, rng);
      const TMatrixDSym sym = randomSym(n, rng);

      const double expected = TMatrixDSym(sym).Similarity(v);
      ASSERT_EQ(tools::similarity(v, sym), expected) << "n = " << n;
    }
  }
}

TEST(ToolsSimilarityVector, RejectsMismatchedDimensions)
{
  std::mt19937_64 rng(31);
  const TVectorD v = randomVector(3, rng);
  const TMatrixDSym sym = randomSym(5, rng);

  EXPECT_THROW(tools::similarity(v, sym), Exception);
}

} // namespace genfit
