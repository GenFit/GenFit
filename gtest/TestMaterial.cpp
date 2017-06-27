#include <gtest/gtest.h>

#include <Material.h>

namespace genfit {

    class MaterialTests : public ::testing::Test {
    protected:
    };

    TEST_F (MaterialTests, Constructor01) {
        Material myMaterial;
        EXPECT_FLOAT_EQ(0, myMaterial.density);
        EXPECT_FLOAT_EQ(0, myMaterial.Z);
        EXPECT_FLOAT_EQ(0, myMaterial.A);
        EXPECT_FLOAT_EQ(0, myMaterial.radiationLength);
        EXPECT_FLOAT_EQ(0, myMaterial.mEE);
    }

    TEST_F (MaterialTests, Constructor02) {
        Material myMaterial(1, 2, 3, 4, 5);
        EXPECT_FLOAT_EQ(1, myMaterial.density);
        EXPECT_FLOAT_EQ(2, myMaterial.Z);
        EXPECT_FLOAT_EQ(3, myMaterial.A);
        EXPECT_FLOAT_EQ(4, myMaterial.radiationLength);
        EXPECT_FLOAT_EQ(5, myMaterial.mEE);
    }

    TEST_F (MaterialTests, ComparisonOperators) {
        Material myMaterial01(0, 0, 0, 0, 0);
        Material myMaterial02(1, 2, 3, 4, 5);
        Material myMaterial03(1, 2, 3, 4, 5);

        EXPECT_FALSE(myMaterial01 == myMaterial02);
        EXPECT_TRUE(myMaterial02 == myMaterial02);
        EXPECT_TRUE(myMaterial02 == myMaterial03);

        EXPECT_TRUE(myMaterial01 != myMaterial02);
        EXPECT_FALSE(myMaterial02 != myMaterial02);
        EXPECT_FALSE(myMaterial02 != myMaterial03);
    }

}
