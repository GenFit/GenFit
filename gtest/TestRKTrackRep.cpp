#include <gtest/gtest.h>

#include <TVector3.h>

#include <Exception.h>
#include <RKTrackRep.h>
#include <ConstField.h>
#include <FieldManager.h>
#include <StateOnPlane.h>
#include <SharedPlanePtr.h>
#include <MeasurementOnPlane.h>


namespace genfit {

    class RKTrackRepTests : public ::testing::Test {

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

    TEST_F (RKTrackRepTests, RKStep) {
        genfit::RKStep myRKStep;

        EXPECT_EQ(0, myRKStep.state7_[0]);
        EXPECT_EQ(0, myRKStep.state7_[1]);
        EXPECT_EQ(0, myRKStep.state7_[2]);
        EXPECT_EQ(0, myRKStep.state7_[3]);
        EXPECT_EQ(0, myRKStep.state7_[4]);
        EXPECT_EQ(0, myRKStep.state7_[5]);
        EXPECT_EQ(0, myRKStep.state7_[6]);
    }

    TEST_F (RKTrackRepTests, ExtrapStep) {
        genfit::ExtrapStep myExtrapStep;

        EXPECT_EQ(0, myExtrapStep.jac7_(0, 0));
        EXPECT_EQ(0, myExtrapStep.jac7_(0, 1));
        EXPECT_EQ(0, myExtrapStep.jac7_(0, 2));
        EXPECT_EQ(0, myExtrapStep.jac7_(0, 3));
        EXPECT_EQ(0, myExtrapStep.jac7_(0, 4));
        EXPECT_EQ(0, myExtrapStep.jac7_(0, 5));
        EXPECT_EQ(0, myExtrapStep.jac7_(0, 6));

        EXPECT_EQ(0, myExtrapStep.jac7_(1, 0));
        EXPECT_EQ(0, myExtrapStep.jac7_(1, 1));
        EXPECT_EQ(0, myExtrapStep.jac7_(1, 2));
        EXPECT_EQ(0, myExtrapStep.jac7_(1, 3));
        EXPECT_EQ(0, myExtrapStep.jac7_(1, 4));
        EXPECT_EQ(0, myExtrapStep.jac7_(1, 5));
        EXPECT_EQ(0, myExtrapStep.jac7_(1, 6));

        EXPECT_EQ(0, myExtrapStep.jac7_(2, 0));
        EXPECT_EQ(0, myExtrapStep.jac7_(2, 1));
        EXPECT_EQ(0, myExtrapStep.jac7_(2, 2));
        EXPECT_EQ(0, myExtrapStep.jac7_(2, 3));
        EXPECT_EQ(0, myExtrapStep.jac7_(2, 4));
        EXPECT_EQ(0, myExtrapStep.jac7_(2, 5));
        EXPECT_EQ(0, myExtrapStep.jac7_(2, 6));

        EXPECT_EQ(0, myExtrapStep.jac7_(3, 0));
        EXPECT_EQ(0, myExtrapStep.jac7_(3, 1));
        EXPECT_EQ(0, myExtrapStep.jac7_(3, 2));
        EXPECT_EQ(0, myExtrapStep.jac7_(3, 3));
        EXPECT_EQ(0, myExtrapStep.jac7_(3, 4));
        EXPECT_EQ(0, myExtrapStep.jac7_(3, 5));
        EXPECT_EQ(0, myExtrapStep.jac7_(3, 6));

        EXPECT_EQ(0, myExtrapStep.jac7_(4, 0));
        EXPECT_EQ(0, myExtrapStep.jac7_(4, 1));
        EXPECT_EQ(0, myExtrapStep.jac7_(4, 2));
        EXPECT_EQ(0, myExtrapStep.jac7_(4, 3));
        EXPECT_EQ(0, myExtrapStep.jac7_(4, 4));
        EXPECT_EQ(0, myExtrapStep.jac7_(4, 5));
        EXPECT_EQ(0, myExtrapStep.jac7_(4, 6));

        EXPECT_EQ(0, myExtrapStep.jac7_(5, 0));
        EXPECT_EQ(0, myExtrapStep.jac7_(5, 1));
        EXPECT_EQ(0, myExtrapStep.jac7_(5, 2));
        EXPECT_EQ(0, myExtrapStep.jac7_(5, 3));
        EXPECT_EQ(0, myExtrapStep.jac7_(5, 4));
        EXPECT_EQ(0, myExtrapStep.jac7_(5, 5));
        EXPECT_EQ(0, myExtrapStep.jac7_(5, 6));

        EXPECT_EQ(0, myExtrapStep.jac7_(6, 0));
        EXPECT_EQ(0, myExtrapStep.jac7_(6, 1));
        EXPECT_EQ(0, myExtrapStep.jac7_(6, 2));
        EXPECT_EQ(0, myExtrapStep.jac7_(6, 3));
        EXPECT_EQ(0, myExtrapStep.jac7_(6, 4));
        EXPECT_EQ(0, myExtrapStep.jac7_(6, 5));
        EXPECT_EQ(0, myExtrapStep.jac7_(6, 6));

        EXPECT_EQ(0, myExtrapStep.noise7_(0, 0));
        EXPECT_EQ(0, myExtrapStep.noise7_(0, 1));
        EXPECT_EQ(0, myExtrapStep.noise7_(0, 2));
        EXPECT_EQ(0, myExtrapStep.noise7_(0, 3));
        EXPECT_EQ(0, myExtrapStep.noise7_(0, 4));
        EXPECT_EQ(0, myExtrapStep.noise7_(0, 5));
        EXPECT_EQ(0, myExtrapStep.noise7_(0, 6));

        EXPECT_EQ(0, myExtrapStep.noise7_(1, 0));
        EXPECT_EQ(0, myExtrapStep.noise7_(1, 1));
        EXPECT_EQ(0, myExtrapStep.noise7_(1, 2));
        EXPECT_EQ(0, myExtrapStep.noise7_(1, 3));
        EXPECT_EQ(0, myExtrapStep.noise7_(1, 4));
        EXPECT_EQ(0, myExtrapStep.noise7_(1, 5));
        EXPECT_EQ(0, myExtrapStep.noise7_(1, 6));

        EXPECT_EQ(0, myExtrapStep.noise7_(2, 0));
        EXPECT_EQ(0, myExtrapStep.noise7_(2, 1));
        EXPECT_EQ(0, myExtrapStep.noise7_(2, 2));
        EXPECT_EQ(0, myExtrapStep.noise7_(2, 3));
        EXPECT_EQ(0, myExtrapStep.noise7_(2, 4));
        EXPECT_EQ(0, myExtrapStep.noise7_(2, 5));
        EXPECT_EQ(0, myExtrapStep.noise7_(2, 6));

        EXPECT_EQ(0, myExtrapStep.noise7_(3, 0));
        EXPECT_EQ(0, myExtrapStep.noise7_(3, 1));
        EXPECT_EQ(0, myExtrapStep.noise7_(3, 2));
        EXPECT_EQ(0, myExtrapStep.noise7_(3, 3));
        EXPECT_EQ(0, myExtrapStep.noise7_(3, 4));
        EXPECT_EQ(0, myExtrapStep.noise7_(3, 5));
        EXPECT_EQ(0, myExtrapStep.noise7_(3, 6));

        EXPECT_EQ(0, myExtrapStep.noise7_(4, 0));
        EXPECT_EQ(0, myExtrapStep.noise7_(4, 1));
        EXPECT_EQ(0, myExtrapStep.noise7_(4, 2));
        EXPECT_EQ(0, myExtrapStep.noise7_(4, 3));
        EXPECT_EQ(0, myExtrapStep.noise7_(4, 4));
        EXPECT_EQ(0, myExtrapStep.noise7_(4, 5));
        EXPECT_EQ(0, myExtrapStep.noise7_(4, 6));

        EXPECT_EQ(0, myExtrapStep.noise7_(5, 0));
        EXPECT_EQ(0, myExtrapStep.noise7_(5, 1));
        EXPECT_EQ(0, myExtrapStep.noise7_(5, 2));
        EXPECT_EQ(0, myExtrapStep.noise7_(5, 3));
        EXPECT_EQ(0, myExtrapStep.noise7_(5, 4));
        EXPECT_EQ(0, myExtrapStep.noise7_(5, 5));
        EXPECT_EQ(0, myExtrapStep.noise7_(5, 6));

        EXPECT_EQ(0, myExtrapStep.noise7_(6, 0));
        EXPECT_EQ(0, myExtrapStep.noise7_(6, 1));
        EXPECT_EQ(0, myExtrapStep.noise7_(6, 2));
        EXPECT_EQ(0, myExtrapStep.noise7_(6, 3));
        EXPECT_EQ(0, myExtrapStep.noise7_(6, 4));
        EXPECT_EQ(0, myExtrapStep.noise7_(6, 5));
        EXPECT_EQ(0, myExtrapStep.noise7_(6, 6));
    }

    /// Black-Box-Test
    TEST_F (RKTrackRepTests, RKPropagate01) {
        genfit::M1x7 myState7;
        myState7.vals[0] = 0;
        myState7.vals[1] = 0;
        myState7.vals[2] = 0;
        myState7.vals[3] = 0;
        myState7.vals[4] = 0;
        myState7.vals[5] = 1;
        myState7.vals[6] = 1;
        genfit::M7x7* myJacobian = nullptr;
        genfit::M1x3 mySA;
        double myS = 0.5;
        const bool varField = false;
        const bool calcOnlyLastRowOfJ = false;

        genfit::RKTrackRep myRKTrackRep;
        EXPECT_EQ(10, myRKTrackRep.RKPropagate(myState7, myJacobian, mySA, myS, varField, calcOnlyLastRowOfJ));
        EXPECT_EQ(0, myState7[0]);
        EXPECT_EQ(0, myState7[1]);
        EXPECT_EQ(0.5, myState7[2]);
        EXPECT_EQ(0, myState7[3]);
        EXPECT_EQ(0, myState7[4]);
        EXPECT_EQ(1, myState7[5]);
        EXPECT_EQ(1, myState7[6]);
        EXPECT_EQ(0, mySA[0]);
        EXPECT_EQ(0, mySA[1]);
        EXPECT_EQ(0, mySA[2]);
    }

    /// Black-Box-Test
    TEST_F (RKTrackRepTests, RKPropagate02) {
        genfit::M1x7 myState7;
        myState7.vals[0] = 1;
        myState7.vals[1] = 1;
        myState7.vals[2] = 1;
        myState7.vals[3] = 1;
        myState7.vals[4] = 1;
        myState7.vals[5] = 1;
        myState7.vals[6] = 1;
        genfit::M7x7* myJacobian = nullptr;
        genfit::M1x3 mySA;
        double myS = 0.1;
        const bool varField = false;
        const bool calcOnlyLastRowOfJ = false;

        genfit::RKTrackRep myRKTrackRep;
        EXPECT_FLOAT_EQ(10, myRKTrackRep.RKPropagate(myState7, myJacobian, mySA, myS, varField, calcOnlyLastRowOfJ));
        EXPECT_FLOAT_EQ(1.1000299732532006, myState7[0]);
        EXPECT_FLOAT_EQ(1.0999700147633968, myState7[1]);
        EXPECT_FLOAT_EQ(1.1, myState7[2]);
        EXPECT_FLOAT_EQ(0.5776963359022331, myState7[3]);
        EXPECT_FLOAT_EQ(0.5770039949184069, myState7[4]);
        EXPECT_FLOAT_EQ(0.5773502691896258, myState7[5]);
        EXPECT_FLOAT_EQ(1, myState7[6]);
        EXPECT_FLOAT_EQ(0.000599405129044106, mySA[0]);
        EXPECT_FLOAT_EQ(-0.0005997646311051152, mySA[1]);
        EXPECT_FLOAT_EQ(0., mySA[2]);
    }

    /// White-Box-Test
    TEST_F (RKTrackRepTests, getState7) {
        RKTrackRep myRKTrackRep;
        genfit::M1x7 myState7;
        myState7.vals[0] = 0;
        myState7.vals[1] = 0;
        myState7.vals[2] = 0;
        myState7.vals[3] = 0;
        myState7.vals[4] = 0;
        myState7.vals[5] = 0;
        myState7.vals[6] = 0;

        genfit::MeasurementOnPlane myMeasurementOnPlane;
        EXPECT_THROW(myRKTrackRep.getState7(myMeasurementOnPlane, myState7), genfit::Exception);
        EXPECT_EQ(0, myState7[0]);
        EXPECT_EQ(0, myState7[1]);
        EXPECT_EQ(0, myState7[2]);
        EXPECT_EQ(0, myState7[3]);
        EXPECT_EQ(0, myState7[4]);
        EXPECT_EQ(0, myState7[5]);
        EXPECT_EQ(0, myState7[6]);

        TVectorD myState5(5);
        myState5[0] = -1;
        myState5[1] = 1;
        myState5[2] = 1;
        myState5[3] = 0.1;
        myState5[4] = 0.1;
        genfit::SharedPlanePtr mySharedPlanePtr(new genfit::DetPlane(TVector3(0, 0, 1), TVector3(0, 0, 1), nullptr));
        genfit::StateOnPlane myStateOnPlane(myState5, mySharedPlanePtr, &myRKTrackRep);
        EXPECT_NO_THROW(myRKTrackRep.getState7(myStateOnPlane, myState7));
        EXPECT_FLOAT_EQ(-0.1, myState7[0]);
        EXPECT_FLOAT_EQ(-0.1, myState7[1]);
        EXPECT_FLOAT_EQ(1, myState7[2]);
        EXPECT_FLOAT_EQ(-0.57735, myState7[3]);
        EXPECT_FLOAT_EQ(-0.57735, myState7[4]);
        EXPECT_FLOAT_EQ(0.57735, myState7[5]);
        EXPECT_FLOAT_EQ(-1, myState7[6]);

    }

    /// White-Box-Test
    TEST_F (RKTrackRepTests, getState5) {
        RKTrackRep myRKTrackRep;
        genfit::M1x7 myState7;
        // state7 must already lie on plane of state
        // We take the values from the getState7 test, because we know them.
        myState7[0] = -0.1;
        myState7[1] = -0.1;
        myState7[2] = 1;
        myState7[3] = -0.57735;
        myState7[4] = -0.57735;
        myState7[5] = 0.57735;
        myState7[6] = -1;
        TVectorD myState5(5);
        genfit::SharedPlanePtr mySharedPlanePtr(new genfit::DetPlane(TVector3(0, 0, 1), TVector3(0, 0, 1), nullptr));
        genfit::StateOnPlane myStateOnPlane(myState5, mySharedPlanePtr, &myRKTrackRep);
        EXPECT_NO_THROW(myRKTrackRep.getState5(myStateOnPlane, myState7));
        EXPECT_FLOAT_EQ(-0.1, myState7[0]);
        EXPECT_FLOAT_EQ(-0.1, myState7[1]);
        EXPECT_FLOAT_EQ(1, myState7[2]);
        EXPECT_FLOAT_EQ(-0.57735, myState7[3]);
        EXPECT_FLOAT_EQ(-0.57735, myState7[4]);
        EXPECT_FLOAT_EQ(0.57735, myState7[5]);
        EXPECT_FLOAT_EQ(-1, myState7[6]);

        EXPECT_EQ(-1, myStateOnPlane.getState()[0]);
        EXPECT_EQ(1, myStateOnPlane.getState()[1]);
        EXPECT_EQ(1, myStateOnPlane.getState()[2]);
        EXPECT_EQ(0.1, myStateOnPlane.getState()[3]);
        EXPECT_EQ(0.1, myStateOnPlane.getState()[4]);
    }

    /// White-Box-Test
    TEST_F (RKTrackRepTests, calcJ_pM_5x7) {
        // TODO: Implement
    }

    /// White-Box-Test
    TEST_F (RKTrackRepTests, transformPM6) {
        // TODO: Implement
    }

    /// White-Box-Test
    TEST_F (RKTrackRepTests, calcJ_Mp_7x5) {
        // TODO: Implement
    }

    /// White-Box-Test
    TEST_F (RKTrackRepTests, calcForwardJacobianAndNoise) {
        const TVector3 startPlaneO(0, 0, 1);
        const TVector3 startPlaneN(0, 0, 1);
        const TVector3 destPlaneO(0, 0, -1);
        const TVector3 destPlaneN(0, 0, 1);

        genfit::RKTrackRep myRKTrackRep;
        genfit::DetPlane myStartPlane(startPlaneO, startPlaneN);
        genfit::DetPlane myDestPlane(destPlaneO, destPlaneN);
        genfit::M1x7 myStartState;
        genfit::M1x7 myDestState;

        EXPECT_EQ(0, myRKTrackRep.ExtrapSteps_.size());
        EXPECT_THROW(myRKTrackRep.calcForwardJacobianAndNoise(myStartState, myStartPlane, myDestState, myDestPlane),
                     genfit::Exception);
    }

    /// White-Box-Test
    TEST_F (RKTrackRepTests, transformM6P) {
        // TODO: Implement
    }

    /// White-Box-Test
    TEST_F (RKTrackRepTests, RKutta) {
        // TODO: Implement
    }

    /// White-Box-Test
    TEST_F (RKTrackRepTests, estimateStep) {
        // TODO: Implement
    }

    /// White-Box-Test
    TEST_F (RKTrackRepTests, Extrap) {
        // TODO: Implement
    }

    /// White-Box-Test
    TEST_F (RKTrackRepTests, momMag) {
        genfit::RKTrackRep myRKTrackRep;
        genfit::M1x7 myState7 = {0, 0, 0, 0, 0, 0, 0};
        EXPECT_EQ(std::numeric_limits<double>::infinity(), myRKTrackRep.momMag(myState7));
        myState7 = {0, 0, 0, 0, 0, 0, 1};
        EXPECT_EQ(1.0, myRKTrackRep.momMag(myState7));
        myState7 = {0, 0, 0, 0, 0, 0, -1};
        EXPECT_EQ(1.0, myRKTrackRep.momMag(myState7));
    }

}
