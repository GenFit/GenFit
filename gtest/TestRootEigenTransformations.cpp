#include <gtest/gtest.h>

#include <TVectorD.h>
#include <Eigen/Dense>

#include <RootEigenTransformations.h>

namespace genfit {

    class RootEigenTransformations : public ::testing::Test {
    protected:
    };

    TEST_F(RootEigenTransformations, VectorTransformation) {
        TVectorD myTVectorD(5);
        myTVectorD(0) = 0.5;
        myTVectorD(0) = 0.4;
        myTVectorD(0) = 0.3;
        myTVectorD(0) = 0.2;
        myTVectorD(0) = 0.1;

        const auto myEigenVector(rootVectorToEigenVector<5>(myTVectorD));
        EXPECT_EQ(myTVectorD(0), myEigenVector(0));
        EXPECT_EQ(myTVectorD(1), myEigenVector(1));
        EXPECT_EQ(myTVectorD(2), myEigenVector(2));
        EXPECT_EQ(myTVectorD(3), myEigenVector(3));
        EXPECT_EQ(myTVectorD(4), myEigenVector(4));

        const auto myRootVector(eigenVectorToRootVector<5>(myEigenVector));
        EXPECT_EQ(myTVectorD(0), myRootVector(0));
        EXPECT_EQ(myTVectorD(1), myRootVector(1));
        EXPECT_EQ(myTVectorD(2), myRootVector(2));
        EXPECT_EQ(myTVectorD(3), myRootVector(3));
        EXPECT_EQ(myTVectorD(4), myRootVector(4));
    }

    TEST_F(RootEigenTransformations, SymmetricMatrixTransformation) {
        TMatrixDSym myTMatrixDSym(5);
        for (unsigned int row=0; row<5; ++row) {
            for (unsigned int col=0; col<5; ++col) {
                myTMatrixDSym(row, col) = row*col;  // Put some non-zero value here.
            }
        }

        const auto myEigenMatrix(rootMatrixSymToEigenMatrix<5>(myTMatrixDSym));
        EXPECT_EQ(myTMatrixDSym(0, 0), myEigenMatrix(0, 0));
        EXPECT_EQ(myTMatrixDSym(0, 1), myEigenMatrix(0, 1));
        EXPECT_EQ(myTMatrixDSym(0, 2), myEigenMatrix(0, 2));
        EXPECT_EQ(myTMatrixDSym(0, 3), myEigenMatrix(0, 3));
        EXPECT_EQ(myTMatrixDSym(0, 4), myEigenMatrix(0, 4));
        EXPECT_EQ(myTMatrixDSym(1, 0), myEigenMatrix(1, 0));
        EXPECT_EQ(myTMatrixDSym(1, 1), myEigenMatrix(1, 1));
        EXPECT_EQ(myTMatrixDSym(1, 2), myEigenMatrix(1, 2));
        EXPECT_EQ(myTMatrixDSym(1, 3), myEigenMatrix(1, 3));
        EXPECT_EQ(myTMatrixDSym(1, 4), myEigenMatrix(1, 4));
        EXPECT_EQ(myTMatrixDSym(2, 0), myEigenMatrix(2, 0));
        EXPECT_EQ(myTMatrixDSym(2, 1), myEigenMatrix(2, 1));
        EXPECT_EQ(myTMatrixDSym(2, 2), myEigenMatrix(2, 2));
        EXPECT_EQ(myTMatrixDSym(2, 3), myEigenMatrix(2, 3));
        EXPECT_EQ(myTMatrixDSym(2, 4), myEigenMatrix(2, 4));
        EXPECT_EQ(myTMatrixDSym(3, 0), myEigenMatrix(3, 0));
        EXPECT_EQ(myTMatrixDSym(3, 1), myEigenMatrix(3, 1));
        EXPECT_EQ(myTMatrixDSym(3, 2), myEigenMatrix(3, 2));
        EXPECT_EQ(myTMatrixDSym(3, 3), myEigenMatrix(3, 3));
        EXPECT_EQ(myTMatrixDSym(3, 4), myEigenMatrix(3, 4));
        EXPECT_EQ(myTMatrixDSym(4, 0), myEigenMatrix(4, 0));
        EXPECT_EQ(myTMatrixDSym(4, 1), myEigenMatrix(4, 1));
        EXPECT_EQ(myTMatrixDSym(4, 2), myEigenMatrix(4, 2));
        EXPECT_EQ(myTMatrixDSym(4, 3), myEigenMatrix(4, 3));
        EXPECT_EQ(myTMatrixDSym(4, 4), myEigenMatrix(4, 4));

        const auto myRootMatrix(eigenMatrixToRootMatrixSym<5>(myEigenMatrix));
        EXPECT_EQ(myTMatrixDSym(0, 0), myRootMatrix(0, 0));
        EXPECT_EQ(myTMatrixDSym(0, 1), myRootMatrix(0, 1));
        EXPECT_EQ(myTMatrixDSym(0, 2), myRootMatrix(0, 2));
        EXPECT_EQ(myTMatrixDSym(0, 3), myRootMatrix(0, 3));
        EXPECT_EQ(myTMatrixDSym(0, 4), myRootMatrix(0, 4));
        EXPECT_EQ(myTMatrixDSym(1, 0), myRootMatrix(1, 0));
        EXPECT_EQ(myTMatrixDSym(1, 1), myRootMatrix(1, 1));
        EXPECT_EQ(myTMatrixDSym(1, 2), myRootMatrix(1, 2));
        EXPECT_EQ(myTMatrixDSym(1, 3), myRootMatrix(1, 3));
        EXPECT_EQ(myTMatrixDSym(1, 4), myRootMatrix(1, 4));
        EXPECT_EQ(myTMatrixDSym(2, 0), myRootMatrix(2, 0));
        EXPECT_EQ(myTMatrixDSym(2, 1), myRootMatrix(2, 1));
        EXPECT_EQ(myTMatrixDSym(2, 2), myRootMatrix(2, 2));
        EXPECT_EQ(myTMatrixDSym(2, 3), myRootMatrix(2, 3));
        EXPECT_EQ(myTMatrixDSym(2, 4), myRootMatrix(2, 4));
        EXPECT_EQ(myTMatrixDSym(3, 0), myRootMatrix(3, 0));
        EXPECT_EQ(myTMatrixDSym(3, 1), myRootMatrix(3, 1));
        EXPECT_EQ(myTMatrixDSym(3, 2), myRootMatrix(3, 2));
        EXPECT_EQ(myTMatrixDSym(3, 3), myRootMatrix(3, 3));
        EXPECT_EQ(myTMatrixDSym(3, 4), myRootMatrix(3, 4));
        EXPECT_EQ(myTMatrixDSym(4, 0), myRootMatrix(4, 0));
        EXPECT_EQ(myTMatrixDSym(4, 1), myRootMatrix(4, 1));
        EXPECT_EQ(myTMatrixDSym(4, 2), myRootMatrix(4, 2));
        EXPECT_EQ(myTMatrixDSym(4, 3), myRootMatrix(4, 3));
        EXPECT_EQ(myTMatrixDSym(4, 4), myRootMatrix(4, 4));

    }

    TEST_F(RootEigenTransformations, MatrixTransformation) {
        const unsigned int rows = 7;
        const unsigned int cols = 5;
        TMatrixD myTMatrix(rows, cols);
        for (unsigned int row=0; row<rows; ++row) {
            for (unsigned int col=0; col<cols; ++col) {
                myTMatrix(row, col) = (row+1)*(col*col+1);  // Put some non-zero value here.
            }
        }

        const auto myEigenMatrix(rootMatrixToEigenMatrix<7, 5>(myTMatrix));
        EXPECT_EQ(myTMatrix(0, 0), myEigenMatrix(0, 0));
        EXPECT_EQ(myTMatrix(0, 1), myEigenMatrix(0, 1));
        EXPECT_EQ(myTMatrix(0, 2), myEigenMatrix(0, 2));
        EXPECT_EQ(myTMatrix(0, 3), myEigenMatrix(0, 3));
        EXPECT_EQ(myTMatrix(0, 4), myEigenMatrix(0, 4));
        EXPECT_EQ(myTMatrix(1, 0), myEigenMatrix(1, 0));
        EXPECT_EQ(myTMatrix(1, 1), myEigenMatrix(1, 1));
        EXPECT_EQ(myTMatrix(1, 2), myEigenMatrix(1, 2));
        EXPECT_EQ(myTMatrix(1, 3), myEigenMatrix(1, 3));
        EXPECT_EQ(myTMatrix(1, 4), myEigenMatrix(1, 4));
        EXPECT_EQ(myTMatrix(2, 0), myEigenMatrix(2, 0));
        EXPECT_EQ(myTMatrix(2, 1), myEigenMatrix(2, 1));
        EXPECT_EQ(myTMatrix(2, 2), myEigenMatrix(2, 2));
        EXPECT_EQ(myTMatrix(2, 3), myEigenMatrix(2, 3));
        EXPECT_EQ(myTMatrix(2, 4), myEigenMatrix(2, 4));
        EXPECT_EQ(myTMatrix(3, 0), myEigenMatrix(3, 0));
        EXPECT_EQ(myTMatrix(3, 1), myEigenMatrix(3, 1));
        EXPECT_EQ(myTMatrix(3, 2), myEigenMatrix(3, 2));
        EXPECT_EQ(myTMatrix(3, 3), myEigenMatrix(3, 3));
        EXPECT_EQ(myTMatrix(3, 4), myEigenMatrix(3, 4));
        EXPECT_EQ(myTMatrix(4, 0), myEigenMatrix(4, 0));
        EXPECT_EQ(myTMatrix(4, 1), myEigenMatrix(4, 1));
        EXPECT_EQ(myTMatrix(4, 2), myEigenMatrix(4, 2));
        EXPECT_EQ(myTMatrix(4, 3), myEigenMatrix(4, 3));
        EXPECT_EQ(myTMatrix(4, 4), myEigenMatrix(4, 4));
        EXPECT_EQ(myTMatrix(5, 0), myEigenMatrix(5, 0));
        EXPECT_EQ(myTMatrix(5, 1), myEigenMatrix(5, 1));
        EXPECT_EQ(myTMatrix(5, 2), myEigenMatrix(5, 2));
        EXPECT_EQ(myTMatrix(5, 3), myEigenMatrix(5, 3));
        EXPECT_EQ(myTMatrix(5, 4), myEigenMatrix(5, 4));
        EXPECT_EQ(myTMatrix(6, 0), myEigenMatrix(6, 0));
        EXPECT_EQ(myTMatrix(6, 1), myEigenMatrix(6, 1));
        EXPECT_EQ(myTMatrix(6, 2), myEigenMatrix(6, 2));
        EXPECT_EQ(myTMatrix(6, 3), myEigenMatrix(6, 3));
        EXPECT_EQ(myTMatrix(6, 4), myEigenMatrix(6, 4));

        const auto myRootMatrix(eigenMatrixToRootMatrix<7, 5>(myEigenMatrix));
        EXPECT_EQ(myTMatrix(0, 0), myRootMatrix(0, 0));
        EXPECT_EQ(myTMatrix(0, 1), myRootMatrix(0, 1));
        EXPECT_EQ(myTMatrix(0, 2), myRootMatrix(0, 2));
        EXPECT_EQ(myTMatrix(0, 3), myRootMatrix(0, 3));
        EXPECT_EQ(myTMatrix(0, 4), myRootMatrix(0, 4));
        EXPECT_EQ(myTMatrix(1, 0), myRootMatrix(1, 0));
        EXPECT_EQ(myTMatrix(1, 1), myRootMatrix(1, 1));
        EXPECT_EQ(myTMatrix(1, 2), myRootMatrix(1, 2));
        EXPECT_EQ(myTMatrix(1, 3), myRootMatrix(1, 3));
        EXPECT_EQ(myTMatrix(1, 4), myRootMatrix(1, 4));
        EXPECT_EQ(myTMatrix(2, 0), myRootMatrix(2, 0));
        EXPECT_EQ(myTMatrix(2, 1), myRootMatrix(2, 1));
        EXPECT_EQ(myTMatrix(2, 2), myRootMatrix(2, 2));
        EXPECT_EQ(myTMatrix(2, 3), myRootMatrix(2, 3));
        EXPECT_EQ(myTMatrix(2, 4), myRootMatrix(2, 4));
        EXPECT_EQ(myTMatrix(3, 0), myRootMatrix(3, 0));
        EXPECT_EQ(myTMatrix(3, 1), myRootMatrix(3, 1));
        EXPECT_EQ(myTMatrix(3, 2), myRootMatrix(3, 2));
        EXPECT_EQ(myTMatrix(3, 3), myRootMatrix(3, 3));
        EXPECT_EQ(myTMatrix(3, 4), myRootMatrix(3, 4));
        EXPECT_EQ(myTMatrix(4, 0), myRootMatrix(4, 0));
        EXPECT_EQ(myTMatrix(4, 1), myRootMatrix(4, 1));
        EXPECT_EQ(myTMatrix(4, 2), myRootMatrix(4, 2));
        EXPECT_EQ(myTMatrix(4, 3), myRootMatrix(4, 3));
        EXPECT_EQ(myTMatrix(4, 4), myRootMatrix(4, 4));
        EXPECT_EQ(myTMatrix(5, 0), myRootMatrix(5, 0));
        EXPECT_EQ(myTMatrix(5, 1), myRootMatrix(5, 1));
        EXPECT_EQ(myTMatrix(5, 2), myRootMatrix(5, 2));
        EXPECT_EQ(myTMatrix(5, 3), myRootMatrix(5, 3));
        EXPECT_EQ(myTMatrix(5, 4), myRootMatrix(5, 4));
        EXPECT_EQ(myTMatrix(6, 0), myRootMatrix(6, 0));
        EXPECT_EQ(myTMatrix(6, 1), myRootMatrix(6, 1));
        EXPECT_EQ(myTMatrix(6, 2), myRootMatrix(6, 2));
        EXPECT_EQ(myTMatrix(6, 3), myRootMatrix(6, 3));
        EXPECT_EQ(myTMatrix(6, 4), myRootMatrix(6, 4));

    }

}