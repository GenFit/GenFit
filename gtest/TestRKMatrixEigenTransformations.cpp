#include <gtest/gtest.h>

#include <RKTools.h>
#include <EigenMatrixTypedefs.h>
#include <RKMatrixEigenTransformations.h>

namespace genfit {

    class RKMatrixEigenTransformations : public::testing::Test {
        protected:
    };

    TEST_F(RKMatrixEigenTransformations, VectorTransformation)
    {
        RKMatrix<1, 4> myVector4 = {{ 1, 2, 3, 4 }};

        auto myEigenMatrix(RKMatrixToEigenMatrix<1, 4>(myVector4));
        EXPECT_EQ(myVector4[0], myEigenMatrix[0]);
        EXPECT_EQ(myVector4[1], myEigenMatrix[1]);
        EXPECT_EQ(myVector4[2], myEigenMatrix[2]);
        EXPECT_EQ(myVector4[3], myEigenMatrix[3]);

        auto myRKMatrix(eigenMatrixToRKMatrix<1, 4>(myEigenMatrix));
        EXPECT_EQ(myVector4[0], myRKMatrix[0]);
        EXPECT_EQ(myVector4[1], myRKMatrix[1]);
        EXPECT_EQ(myVector4[2], myRKMatrix[2]);
        EXPECT_EQ(myVector4[3], myRKMatrix[3]);
    }

    TEST_F(RKMatrixEigenTransformations, MatrixTransformation)
    {
        RKMatrix<2, 3> myMatrix2x3 = {
                1, 2, 3,
                4, 5, 6,
        };

        auto myEigenMatrix(RKMatrixToEigenMatrix<2, 3>(myMatrix2x3));
        EXPECT_EQ(myMatrix2x3(0, 0), myEigenMatrix(0, 0));
        EXPECT_EQ(myMatrix2x3(0, 1), myEigenMatrix(0, 1));
        EXPECT_EQ(myMatrix2x3(0, 2), myEigenMatrix(0, 2));
        EXPECT_EQ(myMatrix2x3(1, 0), myEigenMatrix(1, 0));
        EXPECT_EQ(myMatrix2x3(1, 1), myEigenMatrix(1, 1));
        EXPECT_EQ(myMatrix2x3(1, 2), myEigenMatrix(1, 2));

        auto myRKMatrix(eigenMatrixToRKMatrix<2, 3>(myEigenMatrix));
        EXPECT_EQ(myMatrix2x3(0, 0), myRKMatrix(0, 0));
        EXPECT_EQ(myMatrix2x3(0, 1), myRKMatrix(0, 1));
        EXPECT_EQ(myMatrix2x3(0, 2), myRKMatrix(0, 2));
        EXPECT_EQ(myMatrix2x3(1, 0), myRKMatrix(1, 0));
        EXPECT_EQ(myMatrix2x3(1, 1), myRKMatrix(1, 1));
        EXPECT_EQ(myMatrix2x3(1, 2), myRKMatrix(1, 2));
    }

}