#include <gtest/gtest.h>

#include <EigenMatrixTypedefs.h>
#include <TGeoMaterialInterface.h>
#include <MaterialEffects.h>

#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>


namespace genfit {

    class MaterialEffectsTests : public ::testing::Test {
    protected:
        virtual void SetUp() {

            new TGeoManager("testGeometry", "Test geometry");

            const Scalar Z = 18;
            const Scalar rho = 1.390;  // g / cm^3
            const Scalar A = 39.95; // g/mol;

            TGeoMaterial* liquidArgonMaterial = new TGeoMaterial("liquidArgonMaterial", A, Z, rho);
            TGeoMedium *liquidArgonMedium = new TGeoMedium("liquidArgonMedium", 0, liquidArgonMaterial, 0);

            assert(gGeoManager->GetMaterial(0));
            assert(gGeoManager->GetMedium("liquidArgonMedium"));

            TGeoVolume *top = gGeoManager->MakeBox("Box", liquidArgonMedium, 1000., 1000., 1000.);
            gGeoManager->SetTopVolume(top);

            assert(gGeoManager->GetCurrentVolume());
            assert(gGeoManager->GetCurrentVolume()->GetMedium());
            assert(gGeoManager->GetCurrentVolume()->GetMedium()->GetMaterial());

            genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
        }

        virtual void TearDown() {
            genfit::MaterialEffects::getInstance()->destruct();
        }
    };


    TEST_F (MaterialEffectsTests, Initialization) {
        MaterialProperties myMaterialProperties;
        genfit::MaterialEffects::getInstance()->materialInterface_->getMaterialParameters(myMaterialProperties);
        EXPECT_FLOAT_EQ(18, myMaterialProperties.getZ());
        EXPECT_FLOAT_EQ(1.390, myMaterialProperties.getDensity());
        EXPECT_FLOAT_EQ(39.95, myMaterialProperties.getA());
        EXPECT_FLOAT_EQ(14.04828, myMaterialProperties.getRadLen());
        EXPECT_FLOAT_EQ(188, myMaterialProperties.getMEE());
    }


}