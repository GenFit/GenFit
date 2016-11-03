#include <gtest/gtest.h>

#include <TVector3.h>

#include <Exception.h>
#include <ConstField.h>
#include <FieldManager.h>

namespace genfit {

    class ConstFieldUninitializedTests : public ::testing::Test {
    protected:
    };

    TEST_F (ConstFieldUninitializedTests, UninitializedBehaviour) {
        EXPECT_EQ(false, genfit::FieldManager::getInstance()->isInitialized());
        EXPECT_THROW(genfit::FieldManager::getInstance()->getField(), std::runtime_error);
        EXPECT_THROW(genfit::FieldManager::getInstance()->checkInitialized(), std::runtime_error);
        EXPECT_THROW(genfit::FieldManager::getInstance()->getFieldVal(TVector3(1, 1, 1)), std::runtime_error);
    }


    class ConstFieldInitializedTests : public ::testing::Test {
    protected:
        virtual void SetUp() {
            const double bFieldZ = 20;  // kGauss. Is 1.5T
            m_constField = new genfit::ConstField(0., 0., bFieldZ);
            genfit::FieldManager::getInstance()->init(m_constField);
        }
        virtual void TearDown() {
            genfit::FieldManager::getInstance()->destruct();

        }
        genfit::ConstField* m_constField;
    };

    TEST_F (ConstFieldInitializedTests, InitializedBehaviour) {
        EXPECT_EQ(true, genfit::FieldManager::getInstance()->isInitialized());
        EXPECT_NO_THROW(genfit::FieldManager::getInstance()->checkInstanciated());
        EXPECT_NO_THROW(genfit::FieldManager::getInstance()->checkInitialized());
        EXPECT_EQ(TVector3(0, 0, 20), genfit::FieldManager::getInstance()->getFieldVal(TVector3(1, 1, 1)));
    }

}