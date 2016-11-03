#include <gtest/gtest.h>

#include <AbsMaterialInterface.h>
#include <MaterialEffects.h>

namespace genfit {

    class MaterialEffectsTests : public ::testing::Test {
    protected:
        virtual void SetUp() {
            genfit::MaterialEffects::getInstance()->init(nullptr);
        }

        virtual void TearDown() {
            genfit::MaterialEffects::getInstance()->destruct();
        }
    };


    // TODO: Write a ConstMaterialInterface similiar to the ConstMagneticField for testing purposes.
    // TODO: We can easily check the formulas then! Yeah...
}