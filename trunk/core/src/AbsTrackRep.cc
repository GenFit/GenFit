/* Copyright 2008-2009, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "AbsTrackRep.h"
#include "StateOnPlane.h"
#include "AbsMeasurement.h"
#include "IO.h"

#include <TDatabasePDG.h>



namespace genfit {

AbsTrackRep::AbsTrackRep() :
  pdgCode_(0), propDir_(0), debugLvl_(0)
{
  ;
}

AbsTrackRep::AbsTrackRep(int pdgCode, char propDir) :
  pdgCode_(pdgCode), propDir_(propDir), debugLvl_(0)
{
  ;
}

AbsTrackRep::AbsTrackRep(const AbsTrackRep& rep) :
  TObject(rep), pdgCode_(rep.pdgCode_), propDir_(rep.propDir_), debugLvl_(rep.debugLvl_)
{
  ;
}


double AbsTrackRep::extrapolateToMeasurement(StateOnPlane& state,
    const AbsMeasurement* measurement,
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  return this->extrapolateToPlane(state, measurement->constructPlane(state), stopAtBoundary, calcJacobianNoise);
}


TVectorD AbsTrackRep::get6DState(const StateOnPlane& state) const {
  TVector3 pos, mom;
  getPosMom(state, pos, mom);

  TVectorD stateVec(6);

  stateVec(0) = pos.X();
  stateVec(1) = pos.Y();
  stateVec(2) = pos.Z();

  stateVec(3) = mom.X();
  stateVec(4) = mom.Y();
  stateVec(5) = mom.Z();

  return stateVec;
}


void AbsTrackRep::get6DStateCov(const MeasuredStateOnPlane& state, TVectorD& stateVec, TMatrixDSym& cov) const {
  TVector3 pos, mom;
  getPosMomCov(state, pos, mom, cov);

  stateVec.ResizeTo(6);

  stateVec(0) = pos.X();
  stateVec(1) = pos.Y();
  stateVec(2) = pos.Z();

  stateVec(3) = mom.X();
  stateVec(4) = mom.Y();
  stateVec(5) = mom.Z();
}


double AbsTrackRep::getPDGCharge() const {
  TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdgCode_);
  assert(particle != NULL);
  return particle->Charge()/(3.);
}


double AbsTrackRep::getMass(const StateOnPlane& /*state*/) const {
  return TDatabasePDG::Instance()->GetParticle(pdgCode_)->Mass();
}


void AbsTrackRep::calcJacobianNumerically(const genfit::StateOnPlane& origState,
                                               const genfit::SharedPlanePtr destPlane,
                                               TMatrixD& jacobian) const {

  // Find the transport matrix for track propagation from origState to destPlane
  // I.e. this finds
  //            d stateDestPlane / d origState |_origState

  jacobian.ResizeTo(getDim(), getDim());

  // no science behind these values, I verified that forward and
  // backward propagation yield inverse matrices to good
  // approximation.  In order to avoid bad roundoff errors, the actual
  // step taken is determined below, separately for each direction.
  const double defaultStepX = 1.E-5;
  double stepX;

  // Calculate derivative for all three dimensions successively.
  // The algorithm follows the one in TF1::Derivative() :
  //   df(x) = (4 D(h/2) - D(h)) / 3
  // with D(h) = (f(x + h) - f(x - h)) / (2 h).
  //
  // Could perhaps do better by also using f(x) which would be stB.
  TVectorD rightShort(getDim()), rightFull(getDim());
  TVectorD leftShort(getDim()), leftFull(getDim());
  for (size_t i = 0; i < getDim(); ++i) {
    {
      genfit::StateOnPlane stateCopy(origState);
      double temp = stateCopy.getState()(i) + defaultStepX / 2;
      // Find the actual size of the step, which will differ from
      // defaultStepX due to roundoff.  This is the step-size we will
      // use for this direction.  Idea taken from Numerical Recipes,
      // 3rd ed., section 5.7.
      //
      // Note that if a number is exactly representable, it's double
      // will also be exact.  Outside denormals, this also holds for
      // halving.  Unless the exponent changes (which it only will in
      // the vicinity of zero) adding or subtracing doesn't make a
      // difference.
      //
      // We determine the roundoff error for the half-step.  If this
      // is exactly representable, the full step will also be.
      stepX = 2 * (temp - stateCopy.getState()(i));
      (stateCopy.getState())(i) = temp;
      extrapolateToPlane(stateCopy, destPlane);
      rightShort = stateCopy.getState();
    }
    {
      genfit::StateOnPlane stateCopy(origState);
      (stateCopy.getState())(i) -= stepX / 2;
      extrapolateToPlane(stateCopy, destPlane);
      leftShort = stateCopy.getState();
    }
    {
      genfit::StateOnPlane stateCopy(origState);
      (stateCopy.getState())(i) += stepX;
      extrapolateToPlane(stateCopy, destPlane);
      rightFull = stateCopy.getState();
    }
    {
      genfit::StateOnPlane stateCopy(origState);
      (stateCopy.getState())(i) -= stepX;
      extrapolateToPlane(stateCopy, destPlane);
      leftFull = stateCopy.getState();
    }

    // Calculate the derivatives for the individual components of
    // the track parameters.
    for (size_t j = 0; j < getDim(); ++j) {
      double derivFull = (rightFull(j) - leftFull(j)) / 2 / stepX;
      double derivShort = (rightShort(j) - leftShort(j)) / stepX;

      jacobian(j, i) = 1./3.*(4*derivShort - derivFull);
    }
  }
}


bool AbsTrackRep::switchPDGSign() {
  TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(-pdgCode_);
  if(particle != NULL) {
    pdgCode_ *= -1;
    return true;
  }
  return false;
}



void AbsTrackRep::Print(const Option_t*) const {
  printOut << "genfit::AbsTrackRep, pdgCode = " << pdgCode_ << ". PropDir = " << (int)propDir_ << "\n";
}


} /* End of namespace genfit */
