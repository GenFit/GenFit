#include <gtest/gtest.h>

#include <Units.h>

namespace genfit {

    class UnitTests : public ::testing::Test {
    protected:
    };

    TEST_F (UnitTests, Energy) {
        EXPECT_FLOAT_EQ(1, 1 * Unit::GeV);
        EXPECT_FLOAT_EQ(0.001, 1 * Unit::MeV);
        EXPECT_FLOAT_EQ(0.000001, 1 * Unit::keV);
        EXPECT_FLOAT_EQ(0.000000001, 1 * Unit::eV);
    }

}