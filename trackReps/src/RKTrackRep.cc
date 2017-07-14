/* Copyright 2008-2013, Technische Universitaet Muenchen, Ludwig-Maximilians-Universität München
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch & Tobias Schlüter

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

#include "RKTrackRep.h"
#include "IO.h"
#include "RootEigenTransformations.h"

#include <Exception.h>
#include <FieldManager.h>
#include <MaterialEffects.h>
#include <MeasuredStateOnPlane.h>
#include <MeasurementOnPlane.h>

#include <TBuffer.h>
#include <TDecompLU.h>
#include <TMath.h>

#include <algorithm>

#define MINSTEP 0.001   // minimum step [cm] for Runge Kutta and iteration to POCA

namespace {
  // Use fast inversion instead of LU decomposition?
  const bool useInvertFast = false;
}

namespace genfit {


RKTrackRep::RKTrackRep() :
  AbsTrackRep(),
  lastStartState_(this),
  lastEndState_(this),
  RKStepsFXStart_(0),
  RKStepsFXStop_(0),
  fJacobian_(Matrix5x5::Zero()),
  fNoise_(Matrix5x5Sym::Zero()),
  useCache_(false),
  cachePos_(0)
{
  initArrays();
}


RKTrackRep::RKTrackRep(int pdgCode, char propDir) :
  AbsTrackRep(pdgCode, propDir),
  lastStartState_(this),
  lastEndState_(this),
  RKStepsFXStart_(0),
  RKStepsFXStop_(0),
  fJacobian_(Matrix5x5::Zero()),
  fNoise_(Matrix5x5Sym::Zero()),
  useCache_(false),
  cachePos_(0)
{
  initArrays();
}


RKTrackRep::~RKTrackRep() {
  ;
}


double RKTrackRep::extrapolateToPlane(StateOnPlane& state,
    const SharedPlanePtr& plane,
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    debugOut << "RKTrackRep::extrapolateToPlane()\n";
  }


  if (state.getPlane() == plane) {
    if (debugLvl_ > 0) {
      debugOut << "state is already defined at plane. Do nothing! \n";
    }
    return 0;
  }

  checkCache(state, &plane);

  Vector7 state7(getState7(state));

  TMatrixDSym* covPtr(nullptr);
  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != nullptr) {
    covPtr = &(static_cast<MeasuredStateOnPlane*>(&state)->getCov());
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  // actual extrapolation
  bool isAtBoundary(false);
  double flightTime( 0. );
  double coveredDistance( Extrap(*(state.getPlane()), *plane, getCharge(state), getMass(state), isAtBoundary, state7, flightTime, fillExtrapSteps, covPtr, false, stopAtBoundary) );

  if (stopAtBoundary && isAtBoundary) {
    state.setPlane(SharedPlanePtr(new DetPlane(TVector3(state7[0], state7[1], state7[2]),
                                                TVector3(state7[3], state7[4], state7[5]))));
  }
  else {
    state.setPlane(plane);
  }

  // back to 5D
  getState5(state, state7);
  setTime(state, getTime(state) + flightTime);

  lastEndState_ = state;

  return coveredDistance;
}


double RKTrackRep::extrapolateToLine(StateOnPlane& state,
    const TVector3& linePoint,
    const TVector3& lineDirection,
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    debugOut << "RKTrackRep::extrapolateToLine()\n";
  }

  checkCache(state, nullptr);

  static const unsigned int maxIt(1000);

  Vector7 state7(getState7(state));

  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != nullptr) {
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  double step(0.), lastStep(0.), maxStep(1.E99), angle(0), distToPoca(0), tracklength(0);
  double charge = getCharge(state);
  double mass = getMass(state);
  double flightTime = 0;
  TVector3 dir(state7[3], state7[4], state7[5]);
  TVector3 lastDir(0,0,0);
  TVector3 poca, poca_onwire;
  bool isAtBoundary(false);

  DetPlane startPlane(*(state.getPlane()));
  SharedPlanePtr plane(new DetPlane(linePoint, dir.Cross(lineDirection), lineDirection));
  unsigned int iterations(0);

  while(true){
    if(++iterations == maxIt) {
      Exception exc("RKTrackRep::extrapolateToLine ==> extrapolation to line failed, maximum number of iterations reached",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    lastStep = step;
    lastDir = dir;

    step = this->Extrap(startPlane, *plane, charge, mass, isAtBoundary, state7, flightTime, false, nullptr, true, stopAtBoundary, maxStep);
    tracklength += step;

    dir.SetXYZ(state7[3], state7[4], state7[5]);
    poca.SetXYZ(state7[0], state7[1], state7[2]);
    poca_onwire = pocaOnLine(linePoint, lineDirection, poca);

    // check break conditions
    if (stopAtBoundary && isAtBoundary) {
      plane->setON(dir, poca);
      break;
    }

    angle = fabs(dir.Angle((poca_onwire-poca))-TMath::PiOver2()); // angle between direction and connection to point - 90 deg
    distToPoca = (poca_onwire-poca).Mag();
    if (angle*distToPoca < 0.1*MINSTEP) break;

    // if lastStep and step have opposite sign, the real normal vector lies somewhere between the last two normal vectors (i.e. the directions)
    // -> try mean value of the two (normalization not needed)
    if (lastStep*step < 0){
      dir += lastDir;
      maxStep = 0.5*fabs(lastStep); // make it converge!
    }

    startPlane = *plane;
    plane->setU(dir.Cross(lineDirection));
  }

  if (fillExtrapSteps) { // now do the full extrapolation with covariance matrix
    // make use of the cache
    lastEndState_.setPlane(plane);
    getState5(lastEndState_, state7);

    tracklength = extrapolateToPlane(state, plane, false, true);
    lastEndState_.getAuxInfo()(1) = state.getAuxInfo()(1); // Flight time
  }
  else {
    state.setPlane(plane);
    getState5(state, state7);
    state.getAuxInfo()(1) += flightTime;
  }

  if (debugLvl_ > 0) {
    debugOut << "RKTrackRep::extrapolateToLine(): Reached POCA after " << iterations+1 << " iterations. Distance: " << (poca_onwire-poca).Mag() << " cm. Angle deviation: " << dir.Angle((poca_onwire-poca))-TMath::PiOver2() << " rad \n";
  }

  lastEndState_ = state;

  return tracklength;
}


double RKTrackRep::extrapToPoint(StateOnPlane& state,
    const TVector3& point,
    const TMatrixDSym* G,
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    debugOut << "RKTrackRep::extrapolateToPoint()\n";
  }

  checkCache(state, nullptr);

  static const unsigned int maxIt(1000);

  Vector7 state7(getState7(state));

  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != nullptr) {
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  double step(0.), lastStep(0.), maxStep(1.E99), angle(0), distToPoca(0), tracklength(0);
  TVector3 dir(state7[3], state7[4], state7[5]);
  if (G != nullptr) {
    if(G->GetNrows() != 3) {
      Exception exc("RKTrackRep::extrapolateToLine ==> G is not 3x3",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }
    dir = TMatrix(*G) * dir;
  }
  TVector3 lastDir(0,0,0);

  TVector3 poca;
  bool isAtBoundary(false);

  DetPlane startPlane(*(state.getPlane()));
  SharedPlanePtr plane(new DetPlane(point, dir));
  unsigned int iterations(0);
  double charge = getCharge(state);
  double mass = getMass(state);
  double flightTime = 0;

  while(true){
    if(++iterations == maxIt) {
      Exception exc("RKTrackRep::extrapolateToPoint ==> extrapolation to point failed, maximum number of iterations reached",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    lastStep = step;
    lastDir = dir;

    step = this->Extrap(startPlane, *plane, charge, mass, isAtBoundary, state7, flightTime, false, nullptr, true, stopAtBoundary, maxStep);
    tracklength += step;

    dir.SetXYZ(state7[3], state7[4], state7[5]);
    if (G != nullptr) {
      dir = TMatrix(*G) * dir;
    }
    poca.SetXYZ(state7[0], state7[1], state7[2]);

    // check break conditions
    if (stopAtBoundary && isAtBoundary) {
      plane->setON(dir, poca);
      break;
    }

    angle = fabs(dir.Angle((point-poca))-TMath::PiOver2()); // angle between direction and connection to point - 90 deg
    distToPoca = (point-poca).Mag();
    if (angle*distToPoca < 0.1*MINSTEP) break;

    // if lastStep and step have opposite sign, the real normal vector lies somewhere between the last two normal vectors (i.e. the directions)
    // -> try mean value of the two
    if (lastStep*step < 0){
      if (G != nullptr) { // after multiplication with G, dir has not length 1 anymore in general
        dir.SetMag(1.);
        lastDir.SetMag(1.);
      }
      dir += lastDir;
      maxStep = 0.5*fabs(lastStep); // make it converge!
    }

    startPlane = *plane;
    plane->setNormal(dir);
  }

  if (fillExtrapSteps) { // now do the full extrapolation with covariance matrix
    // make use of the cache
    lastEndState_.setPlane(plane);
    getState5(lastEndState_, state7);

    tracklength = extrapolateToPlane(state, plane, false, true);
    lastEndState_.getAuxInfo()(1) = state.getAuxInfo()(1); // Flight time
  }
  else {
    state.setPlane(plane);
    getState5(state, state7);
    state.getAuxInfo()(1) += flightTime;
  }


  if (debugLvl_ > 0) {
    debugOut << "RKTrackRep::extrapolateToPoint(): Reached POCA after " << iterations+1 << " iterations. Distance: " << (point-poca).Mag() << " cm. Angle deviation: " << dir.Angle((point-poca))-TMath::PiOver2() << " rad \n";
  }

  lastEndState_ = state;

  return tracklength;
}


double RKTrackRep::extrapolateToCylinder(StateOnPlane& state,
    double radius,
    const TVector3& linePoint,
    const TVector3& lineDirection,
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    debugOut << "RKTrackRep::extrapolateToCylinder()\n";
  }

  checkCache(state, nullptr);

  static const unsigned int maxIt(1000);

  Vector7 state7(getState7(state));

  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != nullptr) {
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  double tracklength(0.), maxStep(1.E99);

  TVector3 dest, pos, dir;

  bool isAtBoundary(false);

  DetPlane startPlane(*(state.getPlane()));
  SharedPlanePtr plane(new DetPlane());
  unsigned int iterations(0);
  double charge = getCharge(state);
  double mass = getMass(state);
  double flightTime = 0;

  while(true){
    if(++iterations == maxIt) {
      Exception exc("RKTrackRep::extrapolateToCylinder ==> maximum number of iterations reached",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    pos.SetXYZ(state7[0], state7[1], state7[2]);
    dir.SetXYZ(state7[3], state7[4], state7[5]);

    // solve quadratic equation
    TVector3 AO = (pos - linePoint);
    TVector3 AOxAB = (AO.Cross(lineDirection));
    TVector3 VxAB  = (dir.Cross(lineDirection));
    float ab2    = (lineDirection * lineDirection);
    float a      = (VxAB * VxAB);
    float b      = 2 * (VxAB * AOxAB);
    float c      = (AOxAB * AOxAB) - (radius*radius * ab2);
    double arg = b*b - 4.*a*c;
    if(arg < 0) {
      Exception exc("RKTrackRep::extrapolateToCylinder ==> cannot solve",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }
    double term = sqrt(arg);
    double k1, k2;
    if (b<0) {
      k1 = (-b + term)/(2.*a);
      k2 = 2.*c/(-b + term);
    }
    else {
      k1 = 2.*c/(-b - term);
      k2 = (-b - term)/(2.*a);
    }

    // select smallest absolute solution -> closest cylinder surface
    double k = k1;
    if (fabs(k2)<fabs(k))
    k = k2;

    if (debugLvl_ > 0) {
      debugOut << "RKTrackRep::extrapolateToCylinder(); k = " << k << "\n";
    }

    dest = pos + k * dir;

    plane->setO(dest);
    plane->setUV((dest-linePoint).Cross(lineDirection), lineDirection);

    tracklength += this->Extrap(startPlane, *plane, charge, mass, isAtBoundary, state7, flightTime, false, nullptr, true, stopAtBoundary, maxStep);

    // check break conditions
    if (stopAtBoundary && isAtBoundary) {
      pos.SetXYZ(state7[0], state7[1], state7[2]);
      dir.SetXYZ(state7[3], state7[4], state7[5]);
      plane->setO(pos);
      plane->setUV((pos-linePoint).Cross(lineDirection), lineDirection);
      break;
    }

    if(fabs(k)<MINSTEP) break;

    startPlane = *plane;

  }

  if (fillExtrapSteps) { // now do the full extrapolation with covariance matrix
    // make use of the cache
    lastEndState_.setPlane(plane);
    getState5(lastEndState_, state7);

    tracklength = extrapolateToPlane(state, plane, false, true);
    lastEndState_.getAuxInfo()(1) = state.getAuxInfo()(1); // Flight time
  }
  else {
    state.setPlane(plane);
    getState5(state, state7);
    state.getAuxInfo()(1) += flightTime;
  }

  lastEndState_ = state;

  return tracklength;
}

  
double RKTrackRep::extrapolateToCone(StateOnPlane& state,
    double openingAngle,
    const TVector3& conePoint,
    const TVector3& coneDirection,
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    debugOut << "RKTrackRep::extrapolateToCone()\n";
  }

  checkCache(state, nullptr);

  static const unsigned int maxIt(1000);

  Vector7 state7(getState7(state));

  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != nullptr) {
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  double tracklength(0.), maxStep(1.E99);

  TVector3 dest, pos, dir;

  bool isAtBoundary(false);

  DetPlane startPlane(*(state.getPlane()));
  SharedPlanePtr plane(new DetPlane());
  unsigned int iterations(0);
  double charge = getCharge(state);
  double mass = getMass(state);
  double flightTime = 0;

  while(true){
    if(++iterations == maxIt) {
      Exception exc("RKTrackRep::extrapolateToCone ==> maximum number of iterations reached",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    pos.SetXYZ(state7[0], state7[1], state7[2]);
    dir.SetXYZ(state7[3], state7[4], state7[5]);

    // solve quadratic equation a k^2 + 2 b k + c = 0
    // a = (U . D)^2 - cos^2 alpha * U^2
    // b = (Delta . D) * (U . D) - cos^2 alpha * (U . Delta)
    // c = (Delta . D)^2 - cos^2 alpha * Delta^2
    // Delta = P - V, P track point, U track direction, V cone point, D cone direction, alpha opening angle of cone
    TVector3 cDirection = coneDirection.Unit();
    TVector3 Delta = (pos - conePoint);
    double DirDelta = cDirection * Delta;
    double Delta2 = Delta*Delta;
    double UDir = dir * cDirection;
    double UDelta = dir * Delta;
    double U2 = dir * dir;
    double cosAngle2 = cos(openingAngle)*cos(openingAngle);
    double a = UDir*UDir - cosAngle2*U2;
    double b = UDir*DirDelta - cosAngle2*UDelta;
    double c = DirDelta*DirDelta - cosAngle2*Delta2;
    
    double arg = b*b - a*c;
    if(arg < -1e-9) {
      Exception exc("RKTrackRep::extrapolateToCone ==> cannot solve",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    } else if(arg < 0) {
      arg = 0;
    }

    double term = sqrt(arg);
    double k1, k2;
    k1 = (-b + term) / a;
    k2 = (-b - term) / a;

    // select smallest absolute solution -> closest cone surface
    double k = k1;
    if(fabs(k2) < fabs(k)) {
      k = k2;
    }

    if (debugLvl_ > 0) {
      debugOut << "RKTrackRep::extrapolateToCone(); k = " << k << "\n";
    }

    dest = pos + k * dir;
    // debugOut << "In cone extrapolation ";
    // dest.Print();

    plane->setO(dest);
    plane->setUV((dest-conePoint).Cross(coneDirection), dest-conePoint);

    tracklength += this->Extrap(startPlane, *plane, charge, mass, isAtBoundary, state7, flightTime, false, nullptr, true, stopAtBoundary, maxStep);

    // check break conditions
    if (stopAtBoundary && isAtBoundary) {
      pos.SetXYZ(state7[0], state7[1], state7[2]);
      dir.SetXYZ(state7[3], state7[4], state7[5]);
      plane->setO(pos);
      plane->setUV((pos-conePoint).Cross(coneDirection), pos-conePoint);
      break;
    }

    if(fabs(k)<MINSTEP) break;

    startPlane = *plane;

  }

  if (fillExtrapSteps) { // now do the full extrapolation with covariance matrix
    // make use of the cache
    lastEndState_.setPlane(plane);
    getState5(lastEndState_, state7);

    tracklength = extrapolateToPlane(state, plane, false, true);
    lastEndState_.getAuxInfo()(1) = state.getAuxInfo()(1); // Flight time
  }
  else {
    state.setPlane(plane);
    getState5(state, state7);
    state.getAuxInfo()(1) += flightTime;
  }

  lastEndState_ = state;

  return tracklength;
}


double RKTrackRep::extrapolateToSphere(StateOnPlane& state,
    double radius,
    const TVector3& point, // center
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    debugOut << "RKTrackRep::extrapolateToSphere()\n";
  }

  checkCache(state, nullptr);

  static const unsigned int maxIt(1000);

  Vector7 state7(getState7(state));

  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != nullptr) {
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  double tracklength(0.), maxStep(1.E99);

  TVector3 dest, pos, dir;

  bool isAtBoundary(false);

  DetPlane startPlane(*(state.getPlane()));
  SharedPlanePtr plane(new DetPlane());
  unsigned int iterations(0);
  double charge = getCharge(state);
  double mass = getMass(state);
  double flightTime = 0;

  while(true){
    if(++iterations == maxIt) {
      Exception exc("RKTrackRep::extrapolateToSphere ==> maximum number of iterations reached",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    pos.SetXYZ(state7[0], state7[1], state7[2]);
    dir.SetXYZ(state7[3], state7[4], state7[5]);

    // solve quadratic equation
    TVector3 AO = (pos - point);
    double dirAO = dir * AO;
    double arg = dirAO*dirAO - AO*AO + radius*radius;
    if(arg < 0) {
      Exception exc("RKTrackRep::extrapolateToSphere ==> cannot solve",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }
    double term = sqrt(arg);
    double k1, k2;
    k1 = -dirAO + term;
    k2 = -dirAO - term;

    // select smallest absolute solution -> closest cylinder surface
    double k = k1;
    if (fabs(k2)<fabs(k))
    k = k2;

    if (debugLvl_ > 0) {
      debugOut << "RKTrackRep::extrapolateToSphere(); k = " << k << "\n";
    }

    dest = pos + k * dir;

    plane->setON(dest, dest-point);

    tracklength += this->Extrap(startPlane, *plane, charge, mass, isAtBoundary, state7, flightTime, false, nullptr, true, stopAtBoundary, maxStep);

    // check break conditions
    if (stopAtBoundary && isAtBoundary) {
      pos.SetXYZ(state7[0], state7[1], state7[2]);
      dir.SetXYZ(state7[3], state7[4], state7[5]);
      plane->setON(pos, pos-point);
      break;
    }

    if(fabs(k)<MINSTEP) break;

    startPlane = *plane;

  }

  if (fillExtrapSteps) { // now do the full extrapolation with covariance matrix
    // make use of the cache
    lastEndState_.setPlane(plane);
    getState5(lastEndState_, state7);

    tracklength = extrapolateToPlane(state, plane, false, true);
    lastEndState_.getAuxInfo()(1) = state.getAuxInfo()(1); // Flight time
  }
  else {
    state.setPlane(plane);
    getState5(state, state7);
    state.getAuxInfo()(1) += flightTime;
  }

  lastEndState_ = state;

  return tracklength;
}


double RKTrackRep::extrapolateBy(StateOnPlane& state,
    double step,
    bool stopAtBoundary,
    bool calcJacobianNoise) const {

  if (debugLvl_ > 0) {
    debugOut << "RKTrackRep::extrapolateBy()\n";
  }

  checkCache(state, nullptr);

  static const unsigned int maxIt(1000);

  Vector7 state7(getState7(state));

  bool fillExtrapSteps(false);
  if (dynamic_cast<MeasuredStateOnPlane*>(&state) != nullptr) {
    fillExtrapSteps = true;
  }
  else if (calcJacobianNoise)
    fillExtrapSteps = true;

  double tracklength(0.);

  TVector3 dest, pos, dir;

  bool isAtBoundary(false);

  DetPlane startPlane(*(state.getPlane()));
  SharedPlanePtr plane(new DetPlane());
  unsigned int iterations(0);
  double charge = getCharge(state);
  double mass = getMass(state);
  double flightTime = 0;

  while(true){
    if(++iterations == maxIt) {
      Exception exc("RKTrackRep::extrapolateBy ==> maximum number of iterations reached",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    pos.SetXYZ(state7[0], state7[1], state7[2]);
    dir.SetXYZ(state7[3], state7[4], state7[5]);

    dest = pos + 1.5*(step-tracklength) * dir;

    plane->setON(dest, dir);

    tracklength += this->Extrap(startPlane, *plane, charge, mass, isAtBoundary, state7, flightTime, false, nullptr, true, stopAtBoundary, (step-tracklength));

    // check break conditions
    if (stopAtBoundary && isAtBoundary) {
      pos.SetXYZ(state7[0], state7[1], state7[2]);
      dir.SetXYZ(state7[3], state7[4], state7[5]);
      plane->setON(pos, dir);
      break;
    }

    if (fabs(tracklength-step) < MINSTEP) {
      if (debugLvl_ > 0) {
        debugOut << "RKTrackRep::extrapolateBy(): reached after " << iterations << " iterations. \n";
      }
      pos.SetXYZ(state7[0], state7[1], state7[2]);
      dir.SetXYZ(state7[3], state7[4], state7[5]);
      plane->setON(pos, dir);
      break;
    }

    startPlane = *plane;

  }

  if (fillExtrapSteps) { // now do the full extrapolation with covariance matrix
    // make use of the cache
    lastEndState_.setPlane(plane);
    getState5(lastEndState_, state7);

    tracklength = extrapolateToPlane(state, plane, false, true);
    lastEndState_.getAuxInfo()(1) = state.getAuxInfo()(1); // Flight time
  }
  else {
    state.setPlane(plane);
    getState5(state, state7);
    state.getAuxInfo()(1) += flightTime;
  }

  lastEndState_ = state;

  return tracklength;
}


TVector3 RKTrackRep::getPos(const StateOnPlane& state) const {
  const Vector7 state7(getState7(state));

  return TVector3(state7[0], state7[1], state7[2]);
}


TVector3 RKTrackRep::getMom(const StateOnPlane& state) const {
  const Vector7 state7(getState7(state));

  TVector3 mom(state7[3], state7[4], state7[5]);
  mom.SetMag(getCharge(state)/state7[6]);
  return mom;
}


void RKTrackRep::getPosMom(const StateOnPlane& state, TVector3& pos, TVector3& mom) const {
  const Vector7 state7(getState7(state));

  pos.SetXYZ(state7[0], state7[1], state7[2]);
  mom.SetXYZ(state7[3], state7[4], state7[5]);
  mom.SetMag(getCharge(state)/state7[6]);
}


void RKTrackRep::getPosMomCov(const MeasuredStateOnPlane& state, TVector3& pos, TVector3& mom, TMatrixDSym& cov) const {
  getPosMom(state, pos, mom);
  cov.ResizeTo(6, 6);
  cov = get6DCov(state);
}


TMatrixDSym RKTrackRep::get6DCov(const MeasuredStateOnPlane& state) const {
  return eigenMatrixToRootMatrixSym<6>(transformPM6(state));
}


double RKTrackRep::getCharge(const StateOnPlane& state) const {

  if (dynamic_cast<const MeasurementOnPlane*>(&state) != nullptr) {
    Exception exc("RKTrackRep::getCharge - cannot get charge from MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  double pdgCharge( this->getPDGCharge() );

  // return pdgCharge with sign of q/p
  if (state.getState()(0) * pdgCharge < 0)
    return -pdgCharge;
  else
    return pdgCharge;
}


double RKTrackRep::getMomMag(const StateOnPlane& state) const {
  // p = q / qop
  double p = getCharge(state)/state.getState()(0);
  assert (p>=0);
  return p;
}


double RKTrackRep::getMomVar(const MeasuredStateOnPlane& state) const {

  if (dynamic_cast<const MeasurementOnPlane*>(&state) != nullptr) {
    Exception exc("RKTrackRep::getMomVar - cannot get momVar from MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  // p(qop) = q/qop
  // dp/d(qop) = - q / (qop^2)
  // (delta p) = (delta qop) * |dp/d(qop)| = delta qop * |q / (qop^2)|
  // (var p) = (var qop) * q^2 / (qop^4)

  // delta means sigma
  // cov(0,0) is sigma^2

  return state.getCov()(0,0) * pow(getCharge(state), 2)  / pow(state.getState()(0), 4);
}


double RKTrackRep::getSpu(const StateOnPlane& state) const {

  if (dynamic_cast<const MeasurementOnPlane*>(&state) != nullptr) {
    Exception exc("RKTrackRep::getSpu - cannot get spu from MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  const TVectorD& auxInfo = state.getAuxInfo();
  if (auxInfo.GetNrows() == 2
      || auxInfo.GetNrows() == 1) // backwards compatibility with old RKTrackRep
    return state.getAuxInfo()(0);
  else
    return 1.;
}

double RKTrackRep::getTime(const StateOnPlane& state) const {

  const TVectorD& auxInfo = state.getAuxInfo();
  if (auxInfo.GetNrows() == 2)
    return state.getAuxInfo()(1);
  else
    return 0.;
}


void RKTrackRep::calcForwardJacobianAndNoise(const Vector7& startState7, const DetPlane& startPlane,
                                             const Vector7& destState7, const DetPlane& destPlane) const {

  if (debugLvl_ > 0) {
    debugOut << "RKTrackRep::calcForwardJacobianAndNoise " << std::endl;
  }

  if (ExtrapSteps_.size() == 0) {
    Exception exc("RKTrackRep::calcForwardJacobianAndNoise ==> cache is empty. Extrapolation must run with a MeasuredStateOnPlane.",__LINE__,__FILE__);
    throw exc;
  }

  // The Jacobians returned from RKutta are transposed.
  Matrix7x7 jac(ExtrapSteps_.back().jac7_.transpose());
  Matrix7x7Sym noise(ExtrapSteps_.back().noise7_);
  for (int i = ExtrapSteps_.size() - 2; i >= 0; --i) {
    noise += jac * ExtrapSteps_[i].noise7_ * jac.transpose();
    jac *= ExtrapSteps_[i].jac7_.transpose();
  }

  Matrix5x7 J_pM(calcJ_pM_5x7(startState7, startPlane));
  Matrix7x5 J_Mp(calcJ_Mp_7x5(destState7, destPlane));

  fJacobian_ = J_Mp.transpose() * jac * J_pM.transpose();
  fNoise_ = J_Mp.transpose() * noise * J_Mp;

  if (debugLvl_ > 0) {
    debugOut << "total jacobian : " << std::endl << fJacobian_ << std::endl;
    debugOut << "total noise : " << std::endl << fNoise_ << std::endl;
  }

}


void RKTrackRep::getForwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const {

  jacobian.ResizeTo(5,5);
  jacobian = eigenMatrixToRootMatrix<5, 5>(fJacobian_);

  noise.ResizeTo(5,5);
  noise = eigenMatrixToRootMatrixSym<5>(fNoise_);

  // lastEndState_ = jacobian * lastStartState_  + deltaState
  deltaState.ResizeTo(5);
  // Calculate this without temporaries:
  //deltaState = lastEndState_.getState() - jacobian * lastStartState_.getState()
  deltaState = lastStartState_.getState();
  deltaState *= jacobian;
  deltaState -= lastEndState_.getState();
  deltaState *= -1;


  if (debugLvl_ > 0) {
    debugOut << "delta state : "; deltaState.Print();
  }
}


void RKTrackRep::getBackwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const {

  if (debugLvl_ > 0) {
    debugOut << "RKTrackRep::getBackwardJacobianAndNoise " << std::endl;
  }

  if (ExtrapSteps_.size() == 0) {
    Exception exc("RKTrackRep::getBackwardJacobianAndNoise ==> cache is empty. Extrapolation must run with a MeasuredStateOnPlane.",__LINE__,__FILE__);
    throw exc;
  }

  jacobian.ResizeTo(5,5);
  jacobian = eigenMatrixToRootMatrix<5, 5>(fJacobian_);
  if (!useInvertFast) {
    bool status = TDecompLU::InvertLU(jacobian, 0.0);
    if(status == 0){
      Exception e("cannot invert matrix, status = 0", __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }
  } else {
    double det;
    jacobian.InvertFast(&det);
    if(det < 1e-80){
      Exception e("cannot invert matrix, status = 0", __LINE__,__FILE__);
      e.setFatal();
      throw e;
    }
  }

  noise.ResizeTo(5,5);
  noise = eigenMatrixToRootMatrixSym<5>(fNoise_);
  noise.Similarity(jacobian);

  // lastStartState_ = jacobian * lastEndState_  + deltaState
  deltaState.ResizeTo(5);
  deltaState = lastStartState_.getState() - jacobian * lastEndState_.getState();
}


std::vector<genfit::MatStep> RKTrackRep::getSteps() const {

  // Todo: test

  if (RKSteps_.size() == 0) {
    Exception exc("RKTrackRep::getSteps ==> cache is empty.",__LINE__,__FILE__);
    throw exc;
  }

  std::vector<MatStep> retVal;
  retVal.reserve(RKSteps_.size());

  for (unsigned int i = 0; i<RKSteps_.size(); ++i) {
    retVal.push_back(RKSteps_[i].matStep_);
  }

  return retVal;
}


double RKTrackRep::getRadiationLenght() const {

  // Todo: test

  if (RKSteps_.size() == 0) {
    Exception exc("RKTrackRep::getRadiationLenght ==> cache is empty.",__LINE__,__FILE__);
    throw exc;
  }

  double radLen(0);

  for (unsigned int i = 0; i<RKSteps_.size(); ++i) {
    radLen += RKSteps_.at(i).matStep_.stepSize_ / RKSteps_.at(i).matStep_.material_.radiationLength;
  }

  return radLen;
}



void RKTrackRep::setPosMom(StateOnPlane& state, const TVector3& pos, const TVector3& mom) const {

  if (state.getRep() != this){
    Exception exc("RKTrackRep::setPosMom ==> state is defined wrt. another TrackRep",__LINE__,__FILE__);
    throw exc;
  }

  if (dynamic_cast<MeasurementOnPlane*>(&state) != nullptr) {
    Exception exc("RKTrackRep::setPosMom - cannot set pos/mom of a MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  if (mom.Mag2() == 0) {
    Exception exc("RKTrackRep::setPosMom - momentum is 0",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  // init auxInfo if that has not yet happened
  TVectorD& auxInfo = state.getAuxInfo();
  if (auxInfo.GetNrows() != 2) {
    bool alreadySet = auxInfo.GetNrows() == 1;  // backwards compatibility: don't overwrite old setting
    auxInfo.ResizeTo(2);
    if (!alreadySet)
      setSpu(state, 1.);
  }

  if (state.getPlane() != nullptr && state.getPlane()->distance(pos) < MINSTEP) { // pos is on plane -> do not change plane!

    Vector7 state7;

    state7[0] = pos.X();
    state7[1] = pos.Y();
    state7[2] = pos.Z();

    state7[3] = mom.X();
    state7[4] = mom.Y();
    state7[5] = mom.Z();

    // normalize dir
    double norm = 1. / sqrt(state7[3]*state7[3] + state7[4]*state7[4] + state7[5]*state7[5]);
    for (unsigned int i=3; i<6; ++i)
      state7[i] *= norm;

    state7[6] = getCharge(state) * norm;

    getState5(state, state7);

  }
  else { // pos is not on plane -> create new plane!

    // TODO: Raise Warning that a new plane has been created!
    SharedPlanePtr plane(new DetPlane(pos, mom));
    state.setPlane(plane);

    TVectorD& state5(state.getState());

    state5(0) = getCharge(state)/mom.Mag(); // q/p
    state5(1) = 0.; // u'
    state5(2) = 0.; // v'
    state5(3) = 0.; // u
    state5(4) = 0.; // v

    setSpu(state, 1.);
  }

}


void RKTrackRep::setPosMom(StateOnPlane& state, const TVectorD& state6) const {
  if (state6.GetNrows()!=6){
    Exception exc("RKTrackRep::setPosMom ==> state has to be 6d (x, y, z, px, py, pz)",__LINE__,__FILE__);
    throw exc;
  }
  setPosMom(state, TVector3(state6(0), state6(1), state6(2)), TVector3(state6(3), state6(4), state6(5)));
}


void RKTrackRep::setPosMomErr(MeasuredStateOnPlane& state, const TVector3& pos, const TVector3& mom, const TVector3& posErr, const TVector3& momErr) const {

  // TODO: test!

  setPosMom(state, pos, mom);

  const TVector3& U(state.getPlane()->getU());
  const TVector3& V(state.getPlane()->getV());
  TVector3 W(state.getPlane()->getNormal());

  double pw = mom * W;
  double pu = mom * U;
  double pv = mom * V;

  TMatrixDSym& cov(state.getCov());

  cov(0,0) = pow(getCharge(state), 2) / pow(mom.Mag(), 6) *
             (mom.X()*mom.X() * momErr.X()*momErr.X()+
              mom.Y()*mom.Y() * momErr.Y()*momErr.Y()+
              mom.Z()*mom.Z() * momErr.Z()*momErr.Z());

  cov(1,1) = pow((U.X()/pw - W.X()*pu/(pw*pw)),2.) * momErr.X()*momErr.X() +
             pow((U.Y()/pw - W.Y()*pu/(pw*pw)),2.) * momErr.Y()*momErr.Y() +
             pow((U.Z()/pw - W.Z()*pu/(pw*pw)),2.) * momErr.Z()*momErr.Z();

  cov(2,2) = pow((V.X()/pw - W.X()*pv/(pw*pw)),2.) * momErr.X()*momErr.X() +
             pow((V.Y()/pw - W.Y()*pv/(pw*pw)),2.) * momErr.Y()*momErr.Y() +
             pow((V.Z()/pw - W.Z()*pv/(pw*pw)),2.) * momErr.Z()*momErr.Z();

  cov(3,3) = posErr.X()*posErr.X() * U.X()*U.X() +
             posErr.Y()*posErr.Y() * U.Y()*U.Y() +
             posErr.Z()*posErr.Z() * U.Z()*U.Z();

  cov(4,4) = posErr.X()*posErr.X() * V.X()*V.X() +
             posErr.Y()*posErr.Y() * V.Y()*V.Y() +
             posErr.Z()*posErr.Z() * V.Z()*V.Z();

}




void RKTrackRep::setPosMomCov(MeasuredStateOnPlane& state, const TVector3& pos, const TVector3& mom, const TMatrixDSym& cov6x6) const {

  if (cov6x6.GetNcols()!=6 || cov6x6.GetNrows()!=6){
    Exception exc("RKTrackRep::setPosMomCov ==> cov has to be 6x6 (x, y, z, px, py, pz)",__LINE__,__FILE__);
    throw exc;
  }

  setPosMom(state, pos, mom); // charge does not change!
  transformM6P(rootMatrixSymToEigenMatrix<6>(cov6x6), getState7(state), state);
}

void RKTrackRep::setPosMomCov(MeasuredStateOnPlane& state, const TVectorD& state6, const TMatrixDSym& cov6x6) const {

  if (state6.GetNrows()!=6){
    Exception exc("RKTrackRep::setPosMomCov ==> state has to be 6d (x, y, z, px, py, pz)",__LINE__,__FILE__);
    throw exc;
  }

  if (cov6x6.GetNcols()!=6 || cov6x6.GetNrows()!=6){
    Exception exc("RKTrackRep::setPosMomCov ==> cov has to be 6x6 (x, y, z, px, py, pz)",__LINE__,__FILE__);
    throw exc;
  }

  TVector3 pos(state6(0), state6(1), state6(2));
  TVector3 mom(state6(3), state6(4), state6(5));
  setPosMom(state, pos, mom); // charge does not change!
  transformM6P(rootMatrixSymToEigenMatrix<6>(cov6x6), getState7(state), state);
}


void RKTrackRep::setChargeSign(StateOnPlane& state, double charge) const {

  if (dynamic_cast<MeasurementOnPlane*>(&state) != nullptr) {
    Exception exc("RKTrackRep::setChargeSign - cannot set charge of a MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  if (state.getState()(0) * charge < 0) {
    state.getState()(0) *= -1.;
  }
}


void RKTrackRep::setSpu(StateOnPlane& state, double spu) const {
  state.getAuxInfo().ResizeTo(2);
  (state.getAuxInfo())(0) = spu;
}

void RKTrackRep::setTime(StateOnPlane& state, double time) const {
  state.getAuxInfo().ResizeTo(2);
  (state.getAuxInfo())(1) = time;
}


double RKTrackRep::RKPropagate(Vector7& state7,
                        Matrix7x7* jacobianT,
                        Vector3& SA,
                        double S,
                        bool /*varField*/,
                        bool calcOnlyLastRowOfJ) const {
  // The algorithm is
  //  E Lund et al 2009 JINST 4 P04001 doi:10.1088/1748-0221/4/04/P04001
  //  "Track parameter propagation through the application of a new adaptive Runge-Kutta-Nyström method in the ATLAS experiment"
  //  http://inspirehep.net/search?ln=en&ln=en&p=10.1088/1748-0221/4/04/P04001&of=hb&action_search=Search&sf=earliestdate&so=d&rm=&rg=25&sc=0
  // where the transport of the Jacobian is described in
  //   L. Bugge, J. Myrheim  Nucl.Instrum.Meth. 160 (1979) 43-48
  //   "A Fast Runge-kutta Method For Fitting Tracks In A Magnetic Field"
  //   http://inspirehep.net/record/145692
  // and
  //   L. Bugge, J. Myrheim  Nucl.Instrum.Meth. 179 (1981) 365-381
  //   "Tracking And Track Fitting"
  //   http://inspirehep.net/record/160548

  // important fixed numbers
  static const Scalar EC = 0.000149896229;  // c/(2*10^12) resp. c/2Tera
  static const Scalar P3 = 1./3.;           // 1/3
  static const Scalar DLT = .0002;          // max. deviation for approximation-quality test
  // Aux parameters
  const Vector3& R = state7.block<3, 1>(0, 0);  // Start coordinates  [cm]  (x,  y,  z)
  const Vector3& A = state7.block<3, 1>(3, 0);  // Start directions         (ax, ay, az);   ax^2+ay^2+az^2=1

  Vector3 H0(Vector3::Zero());
  Vector3 H1(Vector3::Zero());
  Vector3 H2(Vector3::Zero());

  // Variables for Runge Kutta solver
  Vector3 X0, X1, X2, X3, X4, X5, X6;

  //
  // Runge Kutta Extrapolation
  //
  const Scalar S3 = P3*S;
  const Scalar S4 = 0.25*S;
  const Scalar PS2 = state7[6]*EC * S;

  // First point
  Vector3 r(R);  // starting point
  FieldManager::getInstance()->getFieldVal(r[0], r[1], r[2], H0[0], H0[1], H0[2]);       // magnetic field in 10^-1 T = kGauss
  H0 *= PS2;  // H0 is PS2*(Hx, Hy, Hz) @ R0

  X0 = A.cross(H0);  // X0 = (ax, ay, az) x H0
  X2 = X0 + A;       // X2 = (A0, B0, C0) + (ax, ay, az)
  X1 = X2 + A;       // X1 = (A0, B0, C0) + 2*(ax, ay, az)

  // Second point
  r += S4 * X1;
  FieldManager::getInstance()->getFieldVal(r[0], r[1], r[2], H1[0], H1[1], H1[2]);
  H1 *= PS2;  // H1 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * [(A0, B0, C0) + 2*(ax, ay, az)]

  X3 = X2.cross(H1) + A;  // X3 = (A2, B2, C2) x H1 + (ax, ay, az)
  X4 = X3.cross(H1) + A;  // X4 = (A3, B3, C3) x H1 + (ax, ay, az)
  X5 = 2 * X4 - A;        // X5 = 2*(A4, B4, C4) - (ax, ay, az)

  // Last point
  r = S * X4 + R;
  FieldManager::getInstance()->getFieldVal(r[0], r[1], r[2], H2[0], H2[1], H2[2]);
  H2 *= PS2;  // H2 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * (A4, B4, C4)

  X6 = X5.cross(H2);  // (A5, B5, C5) x H2


  //
  // Derivatives of track parameters
  //
  if(jacobianT != nullptr){
    Matrix7x7& J_T = *jacobianT;
    // jacobianT
    // 1 0 0 0 0 0 0  x
    // 0 1 0 0 0 0 0  y
    // 0 0 1 0 0 0 0  z
    // x x x x x x 0  a_x
    // x x x x x x 0  a_y
    // x x x x x x 0  a_z
    // x x x x x x 1  q/p
    Vector3 dX0, dX2, dX3, dX4, dX5, dX6;

    if (!calcOnlyLastRowOfJ) {
      for(unsigned int i = 0; i < 6; ++i) {
        //first point
        dX0 = J_T.block<1, 3>(i, 3).transpose().cross(H0);  // dX0 = dA x H0
        dX2 = dX0 + J_T.block<1, 3>(i, 3).transpose();  // dX2 = dX0 + dA

        //second point
        dX3 = J_T.block<1, 3>(i, 3).transpose() + dX2.cross(H1);  // dX3 = dA + dX2 x H1
        dX4 = J_T.block<1, 3>(i, 3).transpose() + dX3.cross(H1);  // dX4 = dA + dX3 x H1

        //last point
        dX5 = dX4 + dX4 - J_T.block<1, 3>(i, 3).transpose();  // dX5 = 2*dX4 - dA
        dX6 = dX5.cross(H2);  // dX6 = dX5 x H2

        // this gives the same results as multiplying the old with the new Jacobian
        // dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
        // dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
        J_T.block<1, 3>(i, 0) += S3 * (dX2 + dX3 + dX4).transpose();
        J_T.block<1, 3>(i, 3) = P3 * ((dX0 + 2*dX3 + dX5 + dX6)).transpose();
      }
    }

    J_T.block<1, 3>(6, 3) *= state7[6];

    //first point
    dX0 = J_T.block<1, 3>(6, 3).transpose().cross(H0) + X0;  // dX0 = dA x H0 + X0
    dX2 = dX0 + J_T.block<1, 3>(6, 3).transpose();  // dX2 = dX0 + dA

    //second point
    dX3 = J_T.block<1, 3>(6, 3).transpose() + dX2.cross(H1) + X3 - A;  // dX3 = dA + dX2 x H1 (+ X3 - A)
    dX4 = J_T.block<1, 3>(6, 3).transpose() + dX3.cross(H1) + X4 - A;  // dX4 = dA + dX3 x H1 (+ X4 - A)

    //last point
    dX5 = dX4 + dX4 - J_T.block<1, 3>(6, 3).transpose();  // dX5 = 2*dX4 - dA
    dX6 = dX5.cross(H2) + X6;  // dX6 = dX5 x H2 + X6

    // this gives the same results as multiplying the old with the new Jacobian
    // dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
    // dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
    J_T.block<1, 3>(6, 0) += S3/state7[6] * (dX2 + dX3 + dX4).transpose();
    J_T.block<1, 3>(6, 3) = P3/state7[6] * (dX0 + 2*dX3 + dX5 + dX6).transpose();
  }

  //
  // Track parameters in last point
  //
  state7.block<3, 1>(0, 0) += S3 * (X2 + X3 + X4);  // R  = R0 + S3*[(A2, B2, C2) +   (A3, B3, C3) + (A4, B4, C4)]
  Vector3 A_new = P3 * (X0 + 2 * X3 +X5 + X6);  // A  = 1/3*[(A0, B0, C0) + 2*(A3, B3, C3) + (A5, B5, C5) + (A6, B6, C6)]
  // FIXME: A_new should already be normalized here, or not?
  SA = A_new - A;  // SA = A_new - A_ol
  state7.block<3, 1>(3, 0) = A_new;
  state7.block<3, 1>(3, 0).normalize();

  // Test approximation quality on given step
  const Scalar EST = ((X1 + X6) - (X3 + X4)).cwiseAbs().rowwise().sum().maxCoeff();
  // EST = ||(ABC1+ABC6)-(ABC3+ABC4)||_1  =  ||(axzy x H0 + ABC5 x H2) - (ABC2 x H1 + ABC3 x H1)||_1
  if (debugLvl_ > 0) {
    debugOut << "    RKTrackRep::RKPropagate. Step = "<< S << "; quality EST = " << EST  << " \n";
  }

  // Prevent the step length increase from getting too large, this is
  // just the point where it becomes 10.
  if (EST < DLT*1e-5)
    return 10;

  // Step length increase for a fifth order Runge-Kutta, see e.g. 17.2
  // in Numerical Recipes.  FIXME: move to caller.
  return pow(DLT/EST, 1./5.);
}



void RKTrackRep::initArrays() const {
  noiseArray_ = Matrix7x7Sym::Zero();
  noiseProjection_ = Matrix7x7Sym::Identity();
  J_MMT_ = Matrix7x7::Zero();

  fJacobian_ = Matrix5x5::Identity();
  fNoise_ = Matrix5x5::Zero();
  limits_.reset();

  RKSteps_.reserve(100);
  ExtrapSteps_.reserve(100);

  lastStartState_.getAuxInfo().ResizeTo(2);
  lastEndState_.getAuxInfo().ResizeTo(2);
}

Vector7 RKTrackRep::getState7(const StateOnPlane& state) const {

    if (dynamic_cast<const MeasurementOnPlane*>(&state) != nullptr) {
        Exception exc("RKTrackRep::getState7 - cannot get pos or mom from a MeasurementOnPlane",__LINE__,__FILE__);
        exc.setFatal();
        throw exc;
    }

    const Vector3 U(TVector3ToEigenVector(state.getPlane()->getU()));
    const Vector3 V(TVector3ToEigenVector(state.getPlane()->getV()));
    const Vector3 O(TVector3ToEigenVector(state.getPlane()->getO()));
    const Vector3 W(TVector3ToEigenVector(state.getPlane()->getNormal()));

    assert(state.getState().GetNrows() == 5);
    const double* state5 = state.getState().GetMatrixArray();
    const Scalar qop = state5[0];
    const Scalar du = state5[1];
    const Scalar dv = state5[2];
    const Scalar u = state5[3];
    const Scalar v = state5[4];

    const Scalar spu = getSpu(state);

    Vector7 state7;
    state7.block<3, 1>(0, 0) = O + u * U + v * V;
    state7.block<3, 1>(3, 0) = (spu * (W + du * U + dv * V)).normalized();
    state7[6] = qop;
    return state7;
}


void RKTrackRep::getState5(StateOnPlane& state, const Vector7& state7) const {

    // state5: (q/p, u', v'. u, v)

    const Vector3 U(TVector3ToEigenVector(state.getPlane()->getU()));
    const Vector3 V(TVector3ToEigenVector(state.getPlane()->getV()));
    const Vector3 O(TVector3ToEigenVector(state.getPlane()->getO()));
    const Vector3 W(TVector3ToEigenVector(state.getPlane()->getNormal()));

    // force A to be in normal direction and set spu accordingly
    const Scalar AtW = state7.block<3, 1>(3, 0).dot(W);
    const Scalar spu = AtW < 0. ? -1 : 1;

    Vector5 state5;
    const Vector3 r = state7.block<3, 1>(0, 0);
    const Vector3 a = state7.block<3, 1>(3, 0);
    const Scalar qop = state7[6];

    state5[0] = qop;  // q/p
    state5[1] = a.dot(U) / AtW;  // u' = (A * U) / (A * W)
    state5[2] = a.dot(V) / AtW;  // v' = (A * V) / (A * W)
    state5[3] = (r - O).dot(U);  // u = (pos - O) * U
    state5[4] = (r - O).dot(V);  // v = (pos - O) * V

    state.setState(eigenVectorToRootVector<5>(state5));
    setSpu(state, spu);

}


Matrix5x7 RKTrackRep::calcJ_pM_5x7(const Vector7& state7, const DetPlane& plane) const {
    Matrix5x7 J_pM(Matrix5x7::Zero());

    const Vector3 normal(TVector3ToEigenVector(plane.getNormal()));
    const Vector3 U(TVector3ToEigenVector(plane.getU()));
    const Vector3 V(TVector3ToEigenVector(plane.getV()));

    Vector3 pTilde(state7.block<3, 1>(3, 0));

    const Scalar pTildeW = pTilde.dot(normal);
    const Scalar spu = pTildeW > 0 ? 1 : -1;
    pTilde *= spu / pTildeW;  // | pTilde * W | has to be 1 (definition of pTilde)

    const Scalar pTildeMag = pTilde.norm();
    const Scalar pTildeMag2 = pTildeMag * pTildeMag;

    const Scalar utpTildeOverpTildeMag2 = U.dot(pTilde) / pTildeMag2;
    const Scalar vtpTildeOverpTildeMag2 = V.dot(pTilde) / pTildeMag2;

    // J_pM matrix is d(x,y,z,ax,ay,az,q/p) / d(q/p,u',v',u,v)
    const Scalar fact = spu / pTildeMag;

    J_pM.block<1, 3>(3, 0) = U;  // d(x,y,z)/d(u)
    J_pM.block<1, 3>(4, 0) = V;  // d(x,y,z)/d(v)
    J_pM(0, 6) = 1.;             // d(q/p)/d(q/p)
    J_pM.block<1, 3>(1, 3) = fact * (U - pTilde * utpTildeOverpTildeMag2);  // d(ax,ay,az)/d(u')
    J_pM.block<1, 3>(2, 3) = fact * (V - pTilde * vtpTildeOverpTildeMag2);  // d(ax,ay,az)/d(v')
    return J_pM;
}


Matrix6x6Sym RKTrackRep::transformPM6(const MeasuredStateOnPlane& state) const {

  // get vectors and aux variables
  const Vector3 U = TVector3ToEigenVector(state.getPlane()->getU());
  const Vector3 V = TVector3ToEigenVector(state.getPlane()->getV());
  const Vector3 W = TVector3ToEigenVector(state.getPlane()->getNormal());

  const TVectorD& state5(state.getState());
  const Scalar spu = getSpu(state);

  const Vector3 pTilde = spu * (W + state5(1) * U + state5(2) * V);

  const Scalar pTildeMag = pTilde.norm();
  const Scalar pTildeMag2 = pTildeMag * pTildeMag;

  const Scalar utpTildeOverpTildeMag2 = U.dot(pTilde) / pTildeMag2;
  const Scalar vtpTildeOverpTildeMag2 = V.dot(pTilde) / pTildeMag2;

  //J_pM matrix is d(x,y,z,px,py,pz) / d(q/p,u',v',u,v)

  const Scalar qop = state5(0);
  const Scalar p = getCharge(state)/qop; // momentum

  Matrix5x6 J_pM_5x6(Matrix5x6::Zero());

  J_pM_5x6.block<1, 3>(0, 3) = -1. * p / (pTildeMag * qop) * pTilde;  // d(px,py,pz)/d(q/p)
  const Scalar fact = p * spu / pTildeMag;
  J_pM_5x6.block<1, 3>(1, 3) = fact * (U - pTilde * utpTildeOverpTildeMag2);  // d(px,py,pz)/d(u')
  J_pM_5x6.block<1, 3>(2, 3) = fact * (V - pTilde * vtpTildeOverpTildeMag2);  // d(px,py,pz)/d(v')
  J_pM_5x6.block<1, 3>(3, 0) = U;  // d(x,y,z)/d(u)
  J_pM_5x6.block<1, 3>(4, 0) = V;  // d(x,y,z)/d(v)

  // cov6x6 = J_pM^T * cov5x5 * J_pM
  return J_pM_5x6.transpose() * rootMatrixToEigenMatrix<5, 5>(state.getCov()) * J_pM_5x6;
}


Matrix7x5 RKTrackRep::calcJ_Mp_7x5(const Vector7& state7, const DetPlane& plane) const {
    Matrix7x5 J_Mp(Matrix7x5::Zero());

    const Vector3& A = state7.block<3, 1>(3, 0);

    const Vector3 U(TVector3ToEigenVector(plane.getU()));
    const Vector3 V(TVector3ToEigenVector(plane.getV()));
    const Vector3 W(TVector3ToEigenVector(plane.getNormal()));

    const Scalar AtU = A.dot(U);
    const Scalar AtV = A.dot(V);
    const Scalar AtW = A.dot(W);

    // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,ax,ay,az,q/p)
    const Scalar fact = 1. / (AtW*AtW);

    J_Mp.block<3, 1>(3, 1) = fact * (AtW * U - AtU * W);  // d(u')/d(ax,ay,az)
    J_Mp.block<3, 1>(3, 2) = fact * (AtW * V - AtV * W);  // d(v')/d(ax,ay,az)
    J_Mp(6, 0) = 1.;                                      // d(q/p)/d(q/p)
    J_Mp.block<3, 1>(0, 3) = U;                           //d(u)/d(x,y,z)
    J_Mp.block<3, 1>(0, 4) = V;                           //d(v)/d(x,y,z)
    return J_Mp;
}


void RKTrackRep::transformM6P(const Matrix6x6Sym& cov, const Vector7& state7, MeasuredStateOnPlane& state) const {

    // get vectors and aux variables
    const Vector3 U = TVector3ToEigenVector(state.getPlane()->getU());
    const Vector3 V = TVector3ToEigenVector(state.getPlane()->getV());
    const Vector3 W = TVector3ToEigenVector(state.getPlane()->getNormal());

    const Scalar AtU = state7.block<3, 1>(3, 0).dot(U);
    const Scalar AtV = state7.block<3, 1>(3, 0).dot(V);
    const Scalar AtW = state7.block<3, 1>(3 ,0).dot(W);

    // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,px,py,pz)

    const Scalar qop = state7[6];
    const Scalar p = getCharge(state)/qop; // momentum

    Matrix6x5 J_Mp_6x5(Matrix6x5::Zero());

    J_Mp_6x5.block<3, 1>(0, 3) = U;  //d(u)/d(x,y,z)
    J_Mp_6x5.block<3, 1>(0, 4) = V;  //d(v)/d(x,y,z)
    J_Mp_6x5.block<3, 1>(3, 0) = (-1.) * qop / p * state7.block<3, 1>(3, 0);  // d(q/p)/d(px,py,pz)

    // d(u')/d(px,py,pz)
    const Scalar fact = 1./(p*AtW*AtW);
    J_Mp_6x5.block<3, 1>(3, 1) = fact * (AtW * U - AtU * W);
    J_Mp_6x5.block<3, 1>(3, 2) = fact * (AtW * V - AtV * W);

    // cov5x5 = J_Mp^T * cov6x6 * J_Mp
    state.setCov(eigenMatrixToRootMatrixSym<5>(J_Mp_6x5.transpose() * cov * J_Mp_6x5));
}


//
// Runge-Kutta method for tracking a particles through a magnetic field.
// Uses Nystroem algorithm (See Handbook Nat. Bur. of Standards, procedure 25.5.20)
// in the way described in
//  E Lund et al 2009 JINST 4 P04001 doi:10.1088/1748-0221/4/04/P04001
//  "Track parameter propagation through the application of a new adaptive Runge-Kutta-Nyström method in the ATLAS experiment"
//  http://inspirehep.net/search?ln=en&ln=en&p=10.1088/1748-0221/4/04/P04001&of=hb&action_search=Search&sf=earliestdate&so=d&rm=&rg=25&sc=0
//
// Input parameters:
//    SU     - plane parameters
//    SU[0]  - direction cosines normal to surface Ex
//    SU[1]  -          -------                    Ey
//    SU[2]  -          -------                    Ez; Ex*Ex+Ey*Ey+Ez*Ez=1
//    SU[3]  - distance to surface from (0,0,0) > 0 cm
//
//    state7 - initial parameters (coordinates(cm), direction,
//             charge/momentum (Gev-1)
//    cov      and derivatives this parameters  (7x7)
//
//    X         Y         Z         Ax        Ay        Az        q/P
//    state7[0] state7[1] state7[2] state7[3] state7[4] state7[5] state7[6]
//
//    dX/dp     dY/dp     dZ/dp     dAx/dp    dAy/dp    dAz/dp    d(q/P)/dp
//    cov[ 0]   cov[ 1]   cov[ 2]   cov[ 3]   cov[ 4]   cov[ 5]   cov[ 6]               d()/dp1
//
//    cov[ 7]   cov[ 8]   cov[ 9]   cov[10]   cov[11]   cov[12]   cov[13]               d()/dp2
//    ............................................................................    d()/dpND
//
// Authors: R.Brun, M.Hansroul, V.Perevoztchikov (Geant3)
//
bool RKTrackRep::RKutta(const Vector4& SU,
                        const DetPlane& plane,
                        double charge,
                        double mass,
                        Vector7& state7,
                        Matrix7x7* jacobianT,
                        Vector7* J_MMT_unprojected_lastRow,
                        double& coveredDistance,
                        double& flightTime,
                        bool& checkJacProj,
                        Matrix7x7Sym& noiseProjection,
                        StepLimits& limits,
                        bool onlyOneStep,
                        bool calcOnlyLastRowOfJ) const {

  // limits, check-values, etc. Can be tuned!
  static const Scalar Wmax           ( 3000. );           // max. way allowed [cm]
  static const Scalar AngleMax       ( 6.3 );           // max. total angle change of momentum. Prevents extrapolating a curler round and round if no active plane is found.
  static const Scalar Pmin           ( 4.E-3 );           // minimum momentum for propagation [GeV]
  static const unsigned int maxNumIt ( 1000 );    // maximum number of iterations in main loop
  // Aux parameters
  Vector3 R = state7.block<3, 1>(0, 0);  // Start coordinates  [cm]  (x,  y,  z)
  Vector3 A = state7.block<3, 1>(3, 0);  // Start directions         (ax, ay, az);   ax^2+ay^2+az^2=1
  Vector3 SA(Vector3::Zero());  // Start directions derivatives dA/S

  Scalar Way        ( 0. );                     // Sum of absolute values of all extrapolation steps [cm]
  Scalar momentum   ( fabs(charge/state7[6]) ); // momentum [GeV]
  Scalar relMomLoss ( 0 );                      // relative momentum loss in RKutta
  Scalar deltaAngle ( 0. );                     // total angle by which the momentum has changed during extrapolation
  Scalar An(0), S(0), Sl(0);


  if (debugLvl_ > 0) {
    debugOut << "RKTrackRep::RKutta \n";
    debugOut << "position: "; TVector3(R[0], R[1], R[2]).Print();
    debugOut << "direction: "; TVector3(A[0], A[1], A[2]).Print();
    debugOut << "momentum: " << momentum << " GeV\n";
    debugOut << "destination: "; plane.Print();
  }

  checkJacProj = false;

  // check momentum
  if(momentum < Pmin){
    std::ostringstream sstream;
    sstream << "RKTrackRep::RKutta ==> momentum too low: " << momentum*1000. << " MeV";
    Exception exc(sstream.str(),__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  unsigned int counter(0);

  // Step estimation (signed)
  S = estimateStep(state7, SU, plane, charge, relMomLoss, limits);

  //
  // Main loop of Runge-Kutta method
  //
  while (fabs(S) >= MINSTEP || counter == 0) {

    if(++counter > maxNumIt){
      Exception exc("RKTrackRep::RKutta ==> maximum number of iterations exceeded",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    if (debugLvl_ > 0) {
      debugOut << "------ RKutta main loop nr. " << counter-1 << " ------\n";
    }

    Vector3 ABefore(A);
    RKPropagate(state7, jacobianT, SA, S, true, calcOnlyLastRowOfJ); // the actual Runge Kutta propagation

    // update paths
    coveredDistance += S;       // add stepsize to way (signed)
    Way  += fabs(S);

    double beta = 1/hypot(1, mass*state7[6]/charge);
    flightTime += S / beta / 29.9792458; // in ns

    // check way limit
    if(Way > Wmax){
      std::ostringstream sstream;
      sstream << "RKTrackRep::RKutta ==> Total extrapolation length is longer than length limit : " << Way << " cm !";
      Exception exc(sstream.str(),__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    if (onlyOneStep) return(true);

    // if stepsize has been limited by material, break the loop and return. No linear extrapolation!
    if (limits.getLowestLimit().first == stp_momLoss) {
      if (debugLvl_ > 0) {
        debugOut<<" momLossExceeded -> return(true); \n";
      }
      return(true);
    }

    // if stepsize has been limited by material boundary, break the loop and return. No linear extrapolation!
    if (limits.getLowestLimit().first == stp_boundary) {
      if (debugLvl_ > 0) {
        debugOut<<" at boundary -> return(true); \n";
      }
      return(true);
    }


    // estimate Step for next loop or linear extrapolation
    Sl = S; // last S used
    limits.removeLimit(stp_fieldCurv);
    limits.removeLimit(stp_momLoss);
    limits.removeLimit(stp_boundary);
    limits.removeLimit(stp_plane);
    S = estimateStep(state7, SU, plane, charge, relMomLoss, limits);

    if (limits.getLowestLimit().first == stp_plane &&
        fabs(S) < MINSTEP) {
      if (debugLvl_ > 0) {
        debugOut<<" (at Plane && fabs(S) < MINSTEP) -> break and do linear extrapolation \n";
      }
      break;
    }
    if (limits.getLowestLimit().first == stp_momLoss &&
        fabs(S) < MINSTEP) {
      if (debugLvl_ > 0) {
        debugOut<<" (momLossExceeded && fabs(S) < MINSTEP) -> return(true), no linear extrapolation; \n";
      }
      RKSteps_.erase(RKSteps_.end()-1);
      --RKStepsFXStop_;
      return(true); // no linear extrapolation!
    }

    // check if total angle is bigger than AngleMax. Can happen if a curler should be fitted and it does not hit the active area of the next plane.
    Scalar arg = ABefore.dot(A);
    arg = arg > 1 ? 1 : arg;
    arg = arg < -1 ? -1 : arg;
    deltaAngle += acos(arg);
    if (fabs(deltaAngle) > AngleMax){
      std::ostringstream sstream;
      sstream << "RKTrackRep::RKutta ==> Do not get to an active plane! Already extrapolated " << deltaAngle * 180 / TMath::Pi() << "°.";
      Exception exc(sstream.str(),__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    // check if we went back and forth multiple times -> we don't come closer to the plane!
    if (counter > 3){
      if (S                            *RKSteps_.at(counter-1).matStep_.stepSize_ < 0 &&
          RKSteps_.at(counter-1).matStep_.stepSize_*RKSteps_.at(counter-2).matStep_.stepSize_ < 0 &&
          RKSteps_.at(counter-2).matStep_.stepSize_*RKSteps_.at(counter-3).matStep_.stepSize_ < 0){
        Exception exc("RKTrackRep::RKutta ==> Do not get closer to plane!",__LINE__,__FILE__);
        exc.setFatal();
        throw exc;
      }
    }

  } //end of main loop


  //
  // linear extrapolation to plane
  //
  if (limits.getLowestLimit().first == stp_plane) {

    if (fabs(Sl) > 0.001*MINSTEP){
      if (debugLvl_ > 0) {
        debugOut << " RKutta - linear extrapolation to surface\n";
      }
      Sl = 1./Sl;        // Sl = inverted last Stepsize Sl

      // normalize SA
      // SA/Sl = delta A / delta way; local derivative of A with respect to the length of the way
      SA *= Sl;
      // calculate A
      // S  = distance to surface
      // A = A + S * SA*Sl
      // normalize A -> 1/|A|
      state7.block<3, 1>(3, 0) += S * SA;
      state7.block<3, 1>(3, 0).normalize();
      A = state7.block<3, 1>(3, 0);
      // R = R + S*(A - 0.5*S*SA); approximation for final point on surface
      state7.block<3, 1>(0, 0) += S * (A - 0.5 * S * SA);
      R = state7.block<3, 1>(0, 0);
      coveredDistance += S;
      Way  += fabs(S);

      double beta = 1/hypot(1, mass*state7[6]/charge);
      flightTime += S / beta / 29.9792458; // in ns;
    }
    else if (debugLvl_ > 0)  {
      debugOut << " RKutta - last stepsize too small -> can't do linear extrapolation! \n";
    }

    //
    // Project Jacobian of extrapolation onto destination plane
    //
    if (jacobianT != nullptr) {

      // projected jacobianT
      // x x x x x x 0
      // x x x x x x 0
      // x x x x x x 0
      // x x x x x x 0
      // x x x x x x 0
      // x x x x x x 0
      // x x x x x x 1

      if (checkJacProj && RKSteps_.size()>0){
        Exception exc("RKTrackRep::Extrap ==> covariance is projected onto destination plane again",__LINE__,__FILE__);
        throw exc;
      }

      if (debugLvl_ > 0) {
        //debugOut << "  Jacobian^T of extrapolation before Projection:\n";
        //RKTools::printDim(*jacobianT, 7,7);
        debugOut << "  Project Jacobian of extrapolation onto destination plane\n";
      }
      An = A[0]*SU[0] + A[1]*SU[1] + A[2]*SU[2];
      An = (fabs(An) > 1.E-7 ? 1./An : 0); // 1/A_normal

      *J_MMT_unprojected_lastRow = (*jacobianT).block<1, 7>(6, 0).transpose();

      for (int row = calcOnlyLastRowOfJ ? 6 : 0; row < 7; ++row) {
        const Scalar norm = (*jacobianT).block<1, 3>(row, 0).transpose().dot(SU.block<3, 1>(0, 0)) * An;  // dR_normal / A_normal
        (*jacobianT).block<1, 3>(row, 0).transpose() -= norm * A;
        (*jacobianT).block<1, 3>(row, 3).transpose() -= norm * SA;
      }
      checkJacProj = true;


      if (debugLvl_ > 0) {
        //debugOut << "  Jacobian^T of extrapolation after Projection:\n";
        //RKTools::printDim(*jacobianT, 7,7);
      }

      if (!calcOnlyLastRowOfJ) {
        for (int iRow = 0; iRow < 3; ++iRow) {
          for (int iCol = 0; iCol < 3; ++iCol) {
            noiseProjection(iRow, iCol)       = (iRow == iCol) - An * SU[iCol] * A[iRow];
            noiseProjection((iRow + 3), iCol) =                - An * SU[iCol] * SA[iRow];
          }
        }

        // noiseProjection will look like this:
        // x x x 0 0 0 0
        // x x x 0 0 0 0
        // x x x 0 0 0 0
        // x x x 1 0 0 0
        // x x x 0 1 0 0
        // x x x 0 0 1 0
        // 0 0 0 0 0 0 1
      }

    }
  } // end of linear extrapolation to surface

  return(true);

}


double RKTrackRep::estimateStep(const Vector7& state7,
                                const Vector4& SU,
                                const DetPlane& plane,
                                const Scalar& charge,
                                Scalar& relMomLoss,
                                StepLimits& limits) const {

  if (useCache_) {
    if (cachePos_ >= RKSteps_.size()) {
      useCache_ = false;
    }
    else {
      if (RKSteps_.at(cachePos_).limits_.getLowestLimit().first == stp_plane) {
        // we need to step exactly to the plane, so don't use the cache!
        useCache_ = false;
        RKSteps_.erase(RKSteps_.begin() + cachePos_, RKSteps_.end());
      }
      else {
        if (debugLvl_ > 0) {
          debugOut << " RKTrackRep::estimateStep: use stepSize " << cachePos_ << " from cache: " << RKSteps_.at(cachePos_).matStep_.stepSize_ << "\n";
        }
        //for(int n = 0; n < 1*7; ++n) RKSteps_[cachePos_].state7_[n] = state7[n];
        ++RKStepsFXStop_;
        limits = RKSteps_.at(cachePos_).limits_;
        return RKSteps_.at(cachePos_++).matStep_.stepSize_;
      }
    }
  }

  limits.setLimit(stp_sMax, 25.); // max. step allowed [cm]

  const auto& R = state7.block<3, 1>(0, 0);
  const auto& A = state7.block<3, 1>(3, 0);

  if (debugLvl_ > 0) {
    debugOut << " RKTrackRep::estimateStep \n";
    debugOut << "  position:  " << R << std::endl;
    debugOut << "  direction: " << A << std::endl;
  }

  // calculate SL distance to surface
  Scalar Dist = SU[3] - R.dot(SU.block<3, 1>(0, 0));  // Distance between start coordinates and surface
  Scalar An = A.dot(SU.block<3, 1>(0, 0));            // An = dir * N;  component of dir normal to surface

  const Scalar SLDist = fabs(An) > 1.E-10 ? Dist/An : fabs(Dist * 1.E10); // signed

  limits.setLimit(stp_plane, SLDist);
  limits.setStepSign(SLDist);

  if (debugLvl_ > 0) {
    debugOut << "  Distance to plane: " << Dist << "\n";
    debugOut << "  SL distance to plane: " << SLDist << "\n";
    if (limits.getStepSign()>0) 
      debugOut << "  Direction is  pointing towards surface.\n";
    else  
      debugOut << "  Direction is pointing away from surface.\n";
  }
  // DONE calculate SL distance to surface

  //
  // Limit according to curvature and magnetic field inhomogenities
  // and improve stepsize estimation to reach plane
  //
  Scalar fieldCurvLimit( limits.getLowestLimitSignedVal() ); // signed
  std::pair<Scalar, Scalar> distVsStep (9.E99, 9.E99); // first: smallest straight line distances to plane; second: RK steps

  static constexpr unsigned int maxNumIt = 10;
  unsigned int counter(0);

  while (fabs(fieldCurvLimit) > MINSTEP) {

    if(++counter > maxNumIt){
      // if max iterations are reached, take a safe value
      // (in previous iteration, fieldCurvLimit has been not more than doubled)
      // and break.
      fieldCurvLimit *= 0.5;
      break;
    }

    Vector7 state7_temp(state7);
    Vector3 SA(Vector3::Zero());

    const auto& R_temp = state7_temp.block<3, 1>(0, 0);
    const auto& A_temp = state7_temp.block<3, 1>(3, 0);

    const Scalar q = RKPropagate(state7_temp, nullptr, SA, fieldCurvLimit, true);
    if (debugLvl_ > 0) {
      debugOut << "  maxStepArg = " << fieldCurvLimit << "; q = " << q  << " \n";
    }

    // remember steps and resulting SL distances to plane for stepsize improvement
    // calculate distance to surface
    Dist = SU[3] - R_temp.dot(SU.block<3, 1>(0, 0));  // Distance between position and surface
    An = A_temp.dot(SU.block<3, 1>(0, 0));            // An = dir * N;  component of dir normal to surface

    if (fabs(Dist/An) < fabs(distVsStep.first)) {
      distVsStep.first = Dist/An;
      distVsStep.second = fieldCurvLimit;
    }

    // resize limit according to q never grow step size more than
    // two-fold to avoid infinite grow-shrink loops with strongly
    // inhomogeneous fields.
    if (q>2) {
      fieldCurvLimit *= 2;
      break;
    }

    fieldCurvLimit *= q * 0.95;

    if (fabs(q-1) < 0.25 || // good enough!
        fabs(fieldCurvLimit) > limits.getLowestLimitVal()) // other limits are lower!
      break;
  }
  if (fabs(fieldCurvLimit) < MINSTEP)
    limits.setLimit(stp_fieldCurv, MINSTEP);
  else
    limits.setLimit(stp_fieldCurv, fieldCurvLimit);

  double stepToPlane(limits.getLimitSigned(stp_plane));
  if (fabs(distVsStep.first) < 8.E99) {
    stepToPlane = distVsStep.first + distVsStep.second;
  }
  limits.setLimit(stp_plane, stepToPlane);


  //
  // Select direction
  //
  // auto select
  if (propDir_ == 0 || !plane.isFinite()){
    if (debugLvl_ > 0) {
      debugOut << "  auto select direction";
      if (!plane.isFinite()) debugOut << ", plane is not finite";
      debugOut << ".\n";
    }
  }
  // see if straight line approximation is ok
  else if ( limits.getLimit(stp_plane) < 0.2*limits.getLimit(stp_fieldCurv) ){
    if (debugLvl_ > 0) {
      debugOut << "  straight line approximation is fine.\n";
    }

    // if direction is pointing to active part of surface
    if( plane.isInActive(state7[0], state7[1], state7[2],  state7[3], state7[4], state7[5]) ) {
      if (debugLvl_ > 0) {
        debugOut << "  direction is pointing to active part of surface. \n";
      }
    }
    // if we are near the plane, but not pointing to the active area, make a big step!
    else {
      limits.removeLimit(stp_plane);
      limits.setStepSign(propDir_);
      if (debugLvl_ > 0) {
        debugOut << "  we are near the plane, but not pointing to the active area. make a big step! \n";
      }
    }
  }
  // propDir_ is set and we are not pointing to an active part of a plane -> propDir_ decides!
  else {
    if (limits.getStepSign() * propDir_ < 0){
      limits.removeLimit(stp_plane);
      limits.setStepSign(propDir_);
      if (debugLvl_ > 0) {
        debugOut << "  invert Step according to propDir_ and make a big step. \n";
      }
    }
  }


  // call stepper and reduce stepsize if step not too small
  static const RKStep defaultRKStep;
  RKSteps_.push_back( defaultRKStep );
  std::vector<RKStep>::iterator lastStep = RKSteps_.end() - 1;
  lastStep->state7_ = state7;
  ++RKStepsFXStop_;

  if(limits.getLowestLimitVal() > MINSTEP){ // only call stepper if step estimation big enough
    Vector7 state7_temp(state7);
    MaterialEffects::getInstance()->stepper(this,
                                            state7_temp,
                                            charge/state7[6], // |p|
                                            relMomLoss,
                                            pdgCode_,
                                            lastStep->matStep_.material_,
                                            limits,
                                            true);
  } else { //assume material has not changed
    if  (RKSteps_.size()>1) {
      lastStep->matStep_.material_ = (lastStep - 1)->matStep_.material_;
    }
  }

  if (debugLvl_ > 0) {
    debugOut << "   final limits:\n";
    limits.Print();
  }

  double finalStep = limits.getLowestLimitSignedVal();

  lastStep->matStep_.stepSize_ = finalStep;
  lastStep->limits_ = limits;

  if (debugLvl_ > 0) {
    debugOut << "  --> Step to be used: " << finalStep << "\n";
  }

  return finalStep;

}


TVector3 RKTrackRep::pocaOnLine(const TVector3& linePoint, const TVector3& lineDirection, const TVector3& point) const {

  TVector3 retVal(lineDirection);

  double t = 1./(retVal.Mag2()) * ((point*retVal) - (linePoint*retVal));
  retVal *= t;
  retVal += linePoint;
  return retVal; // = linePoint + t*lineDirection

}


double RKTrackRep::Extrap(const DetPlane& startPlane,
                          const DetPlane& destPlane,
                          double charge,
                          double mass,
                          bool& isAtBoundary,
                          Vector7& state7,
                          double& flightTime,
                          bool fillExtrapSteps,
                          TMatrixDSym* cov, // 5D
                          bool onlyOneStep,
                          bool stopAtBoundary,
                          double maxStep) const
{

  static const unsigned int maxNumIt(500);
  unsigned int numIt(0);

  double coveredDistance(0.);

  const TVector3 W(destPlane.getNormal());
  Vector4 SU;
  SU << W.X(), W.Y(), W.Z(), destPlane.distance(0., 0., 0.);

  // make SU vector point away from origin
  if (W*destPlane.getO() < 0) {
    SU[0] *= -1;
    SU[1] *= -1;
    SU[2] *= -1;
  }


  Vector7 startState7(state7);

  while(true){

    if (debugLvl_ > 0) {
      debugOut << "\n============ RKTrackRep::Extrap loop nr. " << numIt << " ============\n";
      debugOut << "Start plane: "; startPlane.Print();
      debugOut << "fillExtrapSteps " << fillExtrapSteps << "\n";
    }

    if(++numIt > maxNumIt){
      Exception exc("RKTrackRep::Extrap ==> maximum number of iterations exceeded",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    // initialize jacobianT with unit matrix
    J_MMT_ = Matrix7x7::Identity();

    isAtBoundary = false;

    // propagation
    bool checkJacProj = false;
    limits_.reset();
    limits_.setLimit(stp_sMaxArg, maxStep-fabs(coveredDistance));

    Vector7 J_MMT_unprojected_lastRow;
    J_MMT_unprojected_lastRow << 0, 0, 0, 0, 0, 0, 1;
    bool success = RKutta(SU, destPlane, charge, mass, state7, &J_MMT_, &J_MMT_unprojected_lastRow,
                          coveredDistance, flightTime, checkJacProj, noiseProjection_,
                          limits_, onlyOneStep, !fillExtrapSteps);
    if(not success) {
      Exception exc("RKTrackRep::Extrap ==> Runge Kutta propagation failed",__LINE__,__FILE__);
      exc.setFatal();
      throw exc;
    }

    bool atPlane(limits_.getLowestLimit().first == stp_plane);
    if (limits_.getLowestLimit().first == stp_boundary)
      isAtBoundary = true;


    if (debugLvl_ > 0) {
      debugOut<<"RKSteps \n";
      for (std::vector<RKStep>::iterator it = RKSteps_.begin(); it != RKSteps_.end(); ++it){
        debugOut << "stepSize = " << it->matStep_.stepSize_ << "\t";
        it->matStep_.material_.Print();
      }
      debugOut<<"\n";
    }



    // call MatFX
    Matrix7x7Sym* noise = fillExtrapSteps ? &(noiseArray_ = Matrix7x7Sym::Zero()) : nullptr;

    unsigned int nPoints(RKStepsFXStop_ - RKStepsFXStart_);
    if (/*!fNoMaterial &&*/ nPoints>0){
      // momLoss has a sign - negative loss means momentum gain
      double momLoss = MaterialEffects::getInstance()->effects(RKSteps_,
                                                               RKStepsFXStart_,
                                                               RKStepsFXStop_,
                                                               fabs(charge/state7[6]), // momentum
                                                               pdgCode_,
                                                               noise);

      RKStepsFXStart_ = RKStepsFXStop_;

      if (debugLvl_ > 0) {
        debugOut << "momLoss: " << momLoss << " GeV; relative: " << momLoss/fabs(charge/state7[6])
            << "; coveredDistance = " << coveredDistance << "\n";
        if (debugLvl_ > 1 && noise != nullptr) {
          debugOut << "7D noise: " << std::endl << *noise << std::endl;
        }
      }

      // do momLoss only for defined 1/momentum .ne.0
      if(fabs(state7[6])>1.E-10) {

        if (debugLvl_ > 0) {
          debugOut << "correct state7 with dx/dqop, dy/dqop ...\n";
        }

        // Correct coveredDistance and flightTime and momLoss if checkJacProj == true
        // The idea is to calculate the state correction (based on the mometum loss) twice:
        // Once with the unprojected Jacobian (which preserves coveredDistance),
        // and once with the projected Jacobian (which is constrained to the plane and does NOT preserve coveredDistance).
        // The difference of these two corrections can then be used to calculate a correction factor.
        if (checkJacProj && fabs(coveredDistance) > MINSTEP) {
          const Scalar dqop = charge/(fabs(charge/state7[6])-momLoss) - state7[6];
          const Vector3 state7_correction_unprojected(0.5 * dqop * J_MMT_unprojected_lastRow.block<3, 1>(0, 0));
          const Vector3 state7_correction_projected(0.5 * dqop * J_MMT_.block<1, 3>(6, 0).transpose());
          const Vector3 delta_state(state7_correction_unprojected - state7_correction_projected);
          const Scalar Dist = (delta_state.block<3, 1>(0, 0).dot(state7.block<3, 1>(3, 0))) > 0
                              ? -1 * delta_state.norm() : delta_state.norm();  // sign: delta * a

          const Scalar correctionFactor = 1. + Dist / coveredDistance;
          flightTime *= correctionFactor;
          momLoss *= correctionFactor;
          coveredDistance = coveredDistance + Dist;

          if (debugLvl_ > 0) {
            debugOut << "correctionFactor-1 = " << correctionFactor-1. << "; Dist = " << Dist << "\n";
            debugOut << "corrected momLoss: " << momLoss << " GeV; relative: " << momLoss/fabs(charge/state7[6])
                << "; corrected coveredDistance = " << coveredDistance << "\n";
          }
        }

        // correct state7 with dx/dqop, dy/dqop ... Greatly improves extrapolation accuracy!
        const Scalar qop = charge/(fabs(charge/state7[6])-momLoss);
        const Scalar dqop = qop - state7[6];
        state7[6] = qop;
        state7.block<6, 1>(0, 0) += 0.5 * dqop * J_MMT_.block<1, 6>(6, 0).transpose();
        state7.block<3, 1>(3, 0).normalize();  // normalize direction, just to make sure
      }
    } // finished MatFX


    // fill ExtrapSteps_
    if (fillExtrapSteps) {
      static const ExtrapStep defaultExtrapStep;
      ExtrapSteps_.push_back(defaultExtrapStep);
      std::vector<ExtrapStep>::iterator lastStep = ExtrapSteps_.end() - 1;

      // Store Jacobian of this step for final calculation.
      lastStep->jac7_ = J_MMT_;

      if( checkJacProj == true ){
        //project the noise onto the destPlane
        noiseArray_ = noiseProjection_ * noiseArray_ * noiseProjection_.transpose();

        if (debugLvl_ > 1) {
          debugOut << "7D noise projected onto plane: " << std::endl << noiseArray_ << std::endl;
        }
      }

      // Store this step's noise for final calculation.
      lastStep->noise7_ = noiseArray_;

      if (debugLvl_ > 2) {
        debugOut<<"ExtrapSteps \n";
        for (std::vector<ExtrapStep>::iterator it = ExtrapSteps_.begin(); it != ExtrapSteps_.end(); ++it){
          debugOut << "7D Jacobian: " << std::endl << it->jac7_ << std::endl;
          debugOut << "7D noise:    " << std::endl << it->noise7_ << std::endl;
        }
        debugOut<<"\n";
      }
    }



    // check if at boundary
    if (stopAtBoundary and isAtBoundary) {
      if (debugLvl_ > 0) {
        debugOut << "stopAtBoundary -> break; \n ";
      }
      break;
    }

    if (onlyOneStep) {
      if (debugLvl_ > 0) {
        debugOut << "onlyOneStep -> break; \n ";
      }
      break;
    }

    //break if we arrived at destPlane
    if(atPlane) {
      if (debugLvl_ > 0) {
        debugOut << "arrived at destPlane with a distance of  " << destPlane.distance(state7[0], state7[1], state7[2]) << " cm left. ";
        if (destPlane.isInActive(state7[0], state7[1], state7[2],  state7[3], state7[4], state7[5]))
          debugOut << "In active area of destPlane. \n";
        else
          debugOut << "NOT in active area of plane. \n";

        debugOut << "  position:  "; TVector3(state7[0], state7[1], state7[2]).Print();
        debugOut << "  direction: "; TVector3(state7[3], state7[4], state7[5]).Print();
      }
      break;
    }

  }

  if (fillExtrapSteps) {
    // propagate cov and add noise
    calcForwardJacobianAndNoise(startState7, startPlane, state7, destPlane);

    if (cov != nullptr) {
      cov->Similarity(eigenMatrixToRootMatrix<5, 5>(fJacobian_));
      *cov += eigenMatrixToRootMatrixSym<5>(fNoise_);
    }

    if (debugLvl_ > 0) {
      if (cov != nullptr) {
        debugOut << "final covariance matrix after Extrap: "; cov->Print();
      }
    }
  }

  return coveredDistance;
}


void RKTrackRep::checkCache(const StateOnPlane& state, const SharedPlanePtr* plane) const {

  if (state.getRep() != this){
    Exception exc("RKTrackRep::checkCache ==> state is defined wrt. another TrackRep",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  if (dynamic_cast<const MeasurementOnPlane*>(&state) != nullptr) {
    Exception exc("RKTrackRep::checkCache - cannot extrapolate MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  cachePos_ = 0;
  RKStepsFXStart_ = 0;
  RKStepsFXStop_ = 0;
  ExtrapSteps_.clear();
  initArrays();


  if (plane &&
      lastStartState_.getPlane() &&
      lastEndState_.getPlane() &&
      state.getPlane() == lastStartState_.getPlane() &&
      state.getState() == lastStartState_.getState() &&
      (*plane)->distance(getPos(lastEndState_)) <= MINSTEP) {
    useCache_ = true;

    // clean up cache. Only use steps with same sign.
    double firstStep(0);
    for (unsigned int i=0; i<RKSteps_.size(); ++i) {
      if (i == 0) {
        firstStep = RKSteps_.at(0).matStep_.stepSize_;
        continue;
      }
      if (RKSteps_.at(i).matStep_.stepSize_ * firstStep < 0) {
        if (RKSteps_.at(i-1).matStep_.material_ == RKSteps_.at(i).matStep_.material_) {
          RKSteps_.at(i-1).matStep_.stepSize_ += RKSteps_.at(i).matStep_.stepSize_;
        }
        RKSteps_.erase(RKSteps_.begin()+i, RKSteps_.end());
      }
    }

    if (debugLvl_ > 0) {
        debugOut << "RKTrackRep::checkCache: use cached material and step values.\n";
    }
  }
  else {

    if (debugLvl_ > 0) {
      debugOut << "RKTrackRep::checkCache: can NOT use cached material and step values.\n";

      if (plane != nullptr) {
        if (state.getPlane() != lastStartState_.getPlane()) {
          debugOut << "state.getPlane() != lastStartState_.getPlane()\n";
        }
        else {
          if (! (state.getState() == lastStartState_.getState())) {
            debugOut << "state.getState() != lastStartState_.getState()\n";
          }
          else if (lastEndState_.getPlane().get() != nullptr) {
            debugOut << "distance " << (*plane)->distance(getPos(lastEndState_)) << "\n";
          }
        }
      }
    }

    useCache_ = false;
    RKSteps_.clear();

    lastStartState_.setStatePlane(state.getState(), state.getPlane());
  }
}


bool RKTrackRep::isSameType(const AbsTrackRep* other) {
  if (dynamic_cast<const RKTrackRep*>(other) == nullptr)
    return false;

  return true;
}


bool RKTrackRep::isSame(const AbsTrackRep* other) {
  if (getPDG() != other->getPDG())
    return false;

  return isSameType(other);
}


void RKTrackRep::Streamer(TBuffer &R__b)
{
   // Stream an object of class genfit::RKTrackRep.

   //This works around a msvc bug and should be harmless on other platforms
   typedef ::genfit::RKTrackRep thisClass;
   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      ::genfit::AbsTrackRep::Streamer(R__b);
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      R__b.CheckByteCount(R__s, R__c, thisClass::IsA());
      lastStartState_.setRep(this);
      lastEndState_.setRep(this);
   } else {
      ::genfit::AbsTrackRep::Streamer(R__b);
      R__c = R__b.WriteVersion(thisClass::IsA(), kTRUE);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

} /* End of namespace genfit */
