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
  fJacobian_(5,5),
  fNoise_(5),
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
  fJacobian_(5,5),
  fNoise_(5),
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

  // to 7D
  M1x7 state7 = {{0, 0, 0, 0, 0, 0, 0}};
  getState7(state, state7);

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

  // to 7D
  M1x7 state7;
  getState7(state, state7);

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

  // to 7D
  M1x7 state7;
  getState7(state, state7);

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

  // to 7D
  M1x7 state7;
  getState7(state, state7);

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

  // to 7D
  M1x7 state7;
  getState7(state, state7);

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

  // to 7D
  M1x7 state7;
  getState7(state, state7);

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

  // to 7D
  M1x7 state7;
  getState7(state, state7);

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
  M1x7 state7;
  getState7(state, state7);

  return TVector3(state7[0], state7[1], state7[2]);
}


TVector3 RKTrackRep::getMom(const StateOnPlane& state) const {
  M1x7 state7 = {{0, 0, 0, 0, 0, 0, 0}};
  getState7(state, state7);

  TVector3 mom(state7[3], state7[4], state7[5]);
  mom.SetMag(getCharge(state)/state7[6]);
  return mom;
}


void RKTrackRep::getPosMom(const StateOnPlane& state, TVector3& pos, TVector3& mom) const {
  M1x7 state7 = {{0, 0, 0, 0, 0, 0, 0}};
  getState7(state, state7);

  pos.SetXYZ(state7[0], state7[1], state7[2]);
  mom.SetXYZ(state7[3], state7[4], state7[5]);
  mom.SetMag(getCharge(state)/state7[6]);
}


void RKTrackRep::getPosMomCov(const MeasuredStateOnPlane& state, TVector3& pos, TVector3& mom, TMatrixDSym& cov) const {
  getPosMom(state, pos, mom);
  cov.ResizeTo(6,6);
  transformPM6(state, *((M6x6*) cov.GetMatrixArray()));
}


TMatrixDSym RKTrackRep::get6DCov(const MeasuredStateOnPlane& state) const {
  TMatrixDSym cov(6);
  transformPM6(state, *((M6x6*) cov.GetMatrixArray()));

  return cov;
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


void RKTrackRep::calcForwardJacobianAndNoise(const M1x7& startState7, const DetPlane& startPlane,
					     const M1x7& destState7, const DetPlane& destPlane) const {

  if (debugLvl_ > 0) {
    debugOut << "RKTrackRep::calcForwardJacobianAndNoise " << std::endl;
  }

  if (ExtrapSteps_.size() == 0) {
    Exception exc("RKTrackRep::calcForwardJacobianAndNoise ==> cache is empty. Extrapolation must run with a MeasuredStateOnPlane.",__LINE__,__FILE__);
    throw exc;
  }

  // The Jacobians returned from RKutta are transposed.
  TMatrixD jac(TMatrixD::kTransposed, TMatrixD(7, 7, ExtrapSteps_.back().jac7_.begin()));
  TMatrixDSym noise(7, ExtrapSteps_.back().noise7_.begin());
  for (int i = ExtrapSteps_.size() - 2; i >= 0; --i) {
    noise += TMatrixDSym(7, ExtrapSteps_[i].noise7_.begin()).Similarity(jac);
    jac *= TMatrixD(TMatrixD::kTransposed, TMatrixD(7, 7, ExtrapSteps_[i].jac7_.begin()));
  }

  // Project into 5x5 space.
  M1x3 pTilde = {{startState7[3], startState7[4], startState7[5]}};
  const TVector3& normal = startPlane.getNormal();
  double pTildeW = pTilde[0] * normal.X() + pTilde[1] * normal.Y() + pTilde[2] * normal.Z();
  double spu = pTildeW > 0 ? 1 : -1;
  for (unsigned int i=0; i<3; ++i) {
    pTilde[i] *= spu/pTildeW; // | pTilde * W | has to be 1 (definition of pTilde)
  }
  M5x7 J_pM;
  calcJ_pM_5x7(J_pM, startPlane.getU(), startPlane.getV(), pTilde, spu);
  M7x5 J_Mp;
  calcJ_Mp_7x5(J_Mp, destPlane.getU(), destPlane.getV(), destPlane.getNormal(), *((M1x3*) &destState7[3]));
  jac.Transpose(jac); // Because the helper function wants transposed input.
  RKTools::J_pMTTxJ_MMTTxJ_MpTT(J_Mp, *(M7x7 *)jac.GetMatrixArray(),
				J_pM, *(M5x5 *)fJacobian_.GetMatrixArray());
  RKTools::J_MpTxcov7xJ_Mp(J_Mp, *(M7x7 *)noise.GetMatrixArray(),
			   *(M5x5 *)fNoise_.GetMatrixArray());

  if (debugLvl_ > 0) {
    debugOut << "total jacobian : "; fJacobian_.Print();
    debugOut << "total noise : "; fNoise_.Print();
  }

}


void RKTrackRep::getForwardJacobianAndNoise(TMatrixD& jacobian, TMatrixDSym& noise, TVectorD& deltaState) const {

  jacobian.ResizeTo(5,5);
  jacobian = fJacobian_;

  noise.ResizeTo(5,5);
  noise = fNoise_;

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
  jacobian = fJacobian_;
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
  noise = fNoise_;
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

    M1x7 state7;

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

  M1x7 state7;
  getState7(state, state7);

  const M6x6& cov6x6_( *((M6x6*) cov6x6.GetMatrixArray()) );

  transformM6P(cov6x6_, state7, state);

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

  M1x7 state7;
  getState7(state, state7);

  const M6x6& cov6x6_( *((M6x6*) cov6x6.GetMatrixArray()) );

  transformM6P(cov6x6_, state7, state);

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



double RKTrackRep::RKPropagate(M1x7& state7,
                        M7x7* jacobianT,
                        M1x3& SA,
                        double S,
                        bool varField,
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
  static const double EC  ( 0.000149896229 );  // c/(2*10^12) resp. c/2Tera
  static const double P3  ( 1./3. );           // 1/3
  static const double DLT ( .0002 );           // max. deviation for approximation-quality test
  // Aux parameters
  M1x3&   R           = *((M1x3*) &state7[0]);       // Start coordinates  [cm]  (x,  y,  z)
  M1x3&   A           = *((M1x3*) &state7[3]);       // Start directions         (ax, ay, az);   ax^2+ay^2+az^2=1
  double  S3(0), S4(0), PS2(0);
  M1x3     H0 = {{0.,0.,0.}}, H1 = {{0.,0.,0.}}, H2 = {{0.,0.,0.}};
  M1x3     r = {{0.,0.,0.}};
  // Variables for Runge Kutta solver
  double   A0(0), A1(0), A2(0), A3(0), A4(0), A5(0), A6(0);
  double   B0(0), B1(0), B2(0), B3(0), B4(0), B5(0), B6(0);
  double   C0(0), C1(0), C2(0), C3(0), C4(0), C5(0), C6(0);

  //
  // Runge Kutta Extrapolation
  //
  S3 = P3*S;
  S4 = 0.25*S;
  PS2 = state7[6]*EC * S;

  // First point
  r[0] = R[0];           r[1] = R[1];           r[2]=R[2];
  FieldManager::getInstance()->getFieldVal(r[0], r[1], r[2], H0[0], H0[1], H0[2]);       // magnetic field in 10^-1 T = kGauss
  H0[0] *= PS2; H0[1] *= PS2; H0[2] *= PS2;     // H0 is PS2*(Hx, Hy, Hz) @ R0
  A0 = A[1]*H0[2]-A[2]*H0[1]; B0 = A[2]*H0[0]-A[0]*H0[2]; C0 = A[0]*H0[1]-A[1]*H0[0]; // (ax, ay, az) x H0
  A2 = A[0]+A0              ; B2 = A[1]+B0              ; C2 = A[2]+C0              ; // (A0, B0, C0) + (ax, ay, az)
  A1 = A2+A[0]              ; B1 = B2+A[1]              ; C1 = C2+A[2]              ; // (A0, B0, C0) + 2*(ax, ay, az)

  // Second point
  if (varField) {
    r[0] += A1*S4;         r[1] += B1*S4;         r[2] += C1*S4;
    FieldManager::getInstance()->getFieldVal(r[0], r[1], r[2], H1[0], H1[1], H1[2]);
    H1[0] *= PS2; H1[1] *= PS2; H1[2] *= PS2; // H1 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * [(A0, B0, C0) + 2*(ax, ay, az)]
  }
  else { H1 = H0; };
  A3 = B2*H1[2]-C2*H1[1]+A[0]; B3 = C2*H1[0]-A2*H1[2]+A[1]; C3 = A2*H1[1]-B2*H1[0]+A[2]; // (A2, B2, C2) x H1 + (ax, ay, az)
  A4 = B3*H1[2]-C3*H1[1]+A[0]; B4 = C3*H1[0]-A3*H1[2]+A[1]; C4 = A3*H1[1]-B3*H1[0]+A[2]; // (A3, B3, C3) x H1 + (ax, ay, az)
  A5 = A4-A[0]+A4            ; B5 = B4-A[1]+B4            ; C5 = C4-A[2]+C4            ; //    2*(A4, B4, C4) - (ax, ay, az)

  // Last point
  if (varField) {
    r[0]=R[0]+S*A4;         r[1]=R[1]+S*B4;         r[2]=R[2]+S*C4;  //setup.Field(r,H2);
    FieldManager::getInstance()->getFieldVal(r[0], r[1], r[2], H2[0], H2[1], H2[2]);
    H2[0] *= PS2; H2[1] *= PS2; H2[2] *= PS2; // H2 is PS2*(Hx, Hy, Hz) @ (x, y, z) + 0.25*S * (A4, B4, C4)
  }
  else { H2 = H0; };
  A6 = B5*H2[2]-C5*H2[1]; B6 = C5*H2[0]-A5*H2[2]; C6 = A5*H2[1]-B5*H2[0]; // (A5, B5, C5) x H2


  //
  // Derivatives of track parameters
  //
  if(jacobianT != nullptr){

    // jacobianT
    // 1 0 0 0 0 0 0  x
    // 0 1 0 0 0 0 0  y
    // 0 0 1 0 0 0 0  z
    // x x x x x x 0  a_x
    // x x x x x x 0  a_y
    // x x x x x x 0  a_z
    // x x x x x x 1  q/p
    M7x7& J = *jacobianT;

    double   dA0(0), dA2(0), dA3(0), dA4(0), dA5(0), dA6(0);
    double   dB0(0), dB2(0), dB3(0), dB4(0), dB5(0), dB6(0);
    double   dC0(0), dC2(0), dC3(0), dC4(0), dC5(0), dC6(0);

    int start(0);

    if (!calcOnlyLastRowOfJ) {

      if (!varField) {
        // d(x, y, z)/d(x, y, z) submatrix is unit matrix
        J(0, 0) = 1;  J(1, 1) = 1;  J(2, 2) = 1;
        // d(ax, ay, az)/d(ax, ay, az) submatrix is 0
        // start with d(x, y, z)/d(ax, ay, az)
        start = 3;
      }

      for(int i=start; i<6; ++i) {

        //first point
        dA0 = H0[2]*J(i, 4)-H0[1]*J(i, 5);    // dA0/dp }
        dB0 = H0[0]*J(i, 5)-H0[2]*J(i, 3);    // dB0/dp  } = dA x H0
        dC0 = H0[1]*J(i, 3)-H0[0]*J(i, 4);    // dC0/dp }

        dA2 = dA0+J(i, 3);        // }
        dB2 = dB0+J(i, 4);        //  } = (dA0, dB0, dC0) + dA
        dC2 = dC0+J(i, 5);        // }

        //second point
        dA3 = J(i, 3)+dB2*H1[2]-dC2*H1[1];    // dA3/dp }
        dB3 = J(i, 4)+dC2*H1[0]-dA2*H1[2];    // dB3/dp  } = dA + (dA2, dB2, dC2) x H1
        dC3 = J(i, 5)+dA2*H1[1]-dB2*H1[0];    // dC3/dp }

        dA4 = J(i, 3)+dB3*H1[2]-dC3*H1[1];    // dA4/dp }
        dB4 = J(i, 4)+dC3*H1[0]-dA3*H1[2];    // dB4/dp  } = dA + (dA3, dB3, dC3) x H1
        dC4 = J(i, 5)+dA3*H1[1]-dB3*H1[0];    // dC4/dp }

        //last point
        dA5 = dA4+dA4-J(i, 3);      // }
        dB5 = dB4+dB4-J(i, 4);      //  } =  2*(dA4, dB4, dC4) - dA
        dC5 = dC4+dC4-J(i, 5);      // }

        dA6 = dB5*H2[2]-dC5*H2[1];      // dA6/dp }
        dB6 = dC5*H2[0]-dA5*H2[2];      // dB6/dp  } = (dA5, dB5, dC5) x H2
        dC6 = dA5*H2[1]-dB5*H2[0];      // dC6/dp }

        // this gives the same results as multiplying the old with the new Jacobian
        J(i, 0) += (dA2+dA3+dA4)*S3;  J(i, 3) = ((dA0+2.*dA3)+(dA5+dA6))*P3; // dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
        J(i, 1) += (dB2+dB3+dB4)*S3;  J(i, 4) = ((dB0+2.*dB3)+(dB5+dB6))*P3; // dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
        J(i, 2) += (dC2+dC3+dC4)*S3;  J(i, 5) = ((dC0+2.*dC3)+(dC5+dC6))*P3;
      }

    } // end if (!calcOnlyLastRowOfJ)

    J(6, 3) *= state7[6]; J(6, 4) *= state7[6]; J(6, 5) *= state7[6];

    //first point
    dA0 = H0[2]*J(6, 4)-H0[1]*J(6, 5) + A0;    // dA0/dp }
    dB0 = H0[0]*J(6, 5)-H0[2]*J(6, 3) + B0;    // dB0/dp  } = dA x H0 + (A0, B0, C0)
    dC0 = H0[1]*J(6, 3)-H0[0]*J(6, 4) + C0;    // dC0/dp }

    dA2 = dA0+J(6, 3);        // }
    dB2 = dB0+J(6, 4);        //  } = (dA0, dB0, dC0) + dA
    dC2 = dC0+J(6, 5);        // }

    //second point
    dA3 = J(6, 3)+dB2*H1[2]-dC2*H1[1] + (A3-A[0]);    // dA3/dp }
    dB3 = J(6, 4)+dC2*H1[0]-dA2*H1[2] + (B3-A[1]);    // dB3/dp  } = dA + (dA2, dB2, dC2) x H1
    dC3 = J(6, 5)+dA2*H1[1]-dB2*H1[0] + (C3-A[2]);    // dC3/dp }

    dA4 = J(6, 3)+dB3*H1[2]-dC3*H1[1] + (A4-A[0]);    // dA4/dp }
    dB4 = J(6, 4)+dC3*H1[0]-dA3*H1[2] + (B4-A[1]);    // dB4/dp  } = dA + (dA3, dB3, dC3) x H1
    dC4 = J(6, 5)+dA3*H1[1]-dB3*H1[0] + (C4-A[2]);    // dC4/dp }

    //last point
    dA5 = dA4+dA4-J(6, 3);      // }
    dB5 = dB4+dB4-J(6, 4);      //  } =  2*(dA4, dB4, dC4) - dA
    dC5 = dC4+dC4-J(6, 5);      // }

    dA6 = dB5*H2[2]-dC5*H2[1] + A6;      // dA6/dp }
    dB6 = dC5*H2[0]-dA5*H2[2] + B6;      // dB6/dp  } = (dA5, dB5, dC5) x H2 + (A6, B6, C6)
    dC6 = dA5*H2[1]-dB5*H2[0] + C6;      // dC6/dp }

    // this gives the same results as multiplying the old with the new Jacobian
    J(6, 0) += (dA2+dA3+dA4)*S3/state7[6];  J(6, 3) = ((dA0+2.*dA3)+(dA5+dA6))*P3/state7[6]; // dR := dR + S3*[(dA2, dB2, dC2) +   (dA3, dB3, dC3) + (dA4, dB4, dC4)]
    J(6, 1) += (dB2+dB3+dB4)*S3/state7[6];  J(6, 4) = ((dB0+2.*dB3)+(dB5+dB6))*P3/state7[6]; // dA :=     1/3*[(dA0, dB0, dC0) + 2*(dA3, dB3, dC3) + (dA5, dB5, dC5) + (dA6, dB6, dC6)]
    J(6, 2) += (dC2+dC3+dC4)*S3/state7[6];  J(6, 5) = ((dC0+2.*dC3)+(dC5+dC6))*P3/state7[6];

  }

  //
  // Track parameters in last point
  //
  R[0] += (A2+A3+A4)*S3;   A[0] += (SA[0]=((A0+2.*A3)+(A5+A6))*P3-A[0]);  // R  = R0 + S3*[(A2, B2, C2) +   (A3, B3, C3) + (A4, B4, C4)]
  R[1] += (B2+B3+B4)*S3;   A[1] += (SA[1]=((B0+2.*B3)+(B5+B6))*P3-A[1]);  // A  =     1/3*[(A0, B0, C0) + 2*(A3, B3, C3) + (A5, B5, C5) + (A6, B6, C6)]
  R[2] += (C2+C3+C4)*S3;   A[2] += (SA[2]=((C0+2.*C3)+(C5+C6))*P3-A[2]);  // SA = A_new - A_old

  // normalize A
  double CBA ( 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]) ); // 1/|A|
  A[0] *= CBA; A[1] *= CBA; A[2] *= CBA;


  // Test approximation quality on given step
  double EST ( fabs((A1+A6)-(A3+A4)) +
               fabs((B1+B6)-(B3+B4)) +
               fabs((C1+C6)-(C3+C4))  );  // EST = ||(ABC1+ABC6)-(ABC3+ABC4)||_1  =  ||(axzy x H0 + ABC5 x H2) - (ABC2 x H1 + ABC3 x H1)||_1
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
  std::fill(noiseArray_.begin(), noiseArray_.end(), 0);
  std::fill(noiseProjection_.begin(), noiseProjection_.end(), 0);
  for (unsigned int i=0; i<7; ++i) // initialize as diagonal matrix
    noiseProjection_[i*8] = 1;
  std::fill(J_MMT_.begin(), J_MMT_.end(), 0);

  fJacobian_.UnitMatrix();
  fNoise_.Zero();
  limits_.reset();

  RKSteps_.reserve(100);
  ExtrapSteps_.reserve(100);

  lastStartState_.getAuxInfo().ResizeTo(2);
  lastEndState_.getAuxInfo().ResizeTo(2);
}


void RKTrackRep::getState7(const StateOnPlane& state, M1x7& state7) const {

  if (dynamic_cast<const MeasurementOnPlane*>(&state) != nullptr) {
    Exception exc("RKTrackRep::getState7 - cannot get pos or mom from a MeasurementOnPlane",__LINE__,__FILE__);
    exc.setFatal();
    throw exc;
  }

  const TVector3& U(state.getPlane()->getU());
  const TVector3& V(state.getPlane()->getV());
  const TVector3& O(state.getPlane()->getO());
  const TVector3& W(state.getPlane()->getNormal());

  assert(state.getState().GetNrows() == 5);
  const double* state5 = state.getState().GetMatrixArray();

  double spu = getSpu(state);

  state7[0] = O.X() + state5[3]*U.X() + state5[4]*V.X(); // x
  state7[1] = O.Y() + state5[3]*U.Y() + state5[4]*V.Y(); // y
  state7[2] = O.Z() + state5[3]*U.Z() + state5[4]*V.Z(); // z

  state7[3] = spu * (W.X() + state5[1]*U.X() + state5[2]*V.X()); // a_x
  state7[4] = spu * (W.Y() + state5[1]*U.Y() + state5[2]*V.Y()); // a_y
  state7[5] = spu * (W.Z() + state5[1]*U.Z() + state5[2]*V.Z()); // a_z

  // normalize dir
  double norm = 1. / sqrt(state7[3]*state7[3] + state7[4]*state7[4] + state7[5]*state7[5]);
  for (unsigned int i=3; i<6; ++i) state7[i] *= norm;

  state7[6] = state5[0]; // q/p
}


void RKTrackRep::getState5(StateOnPlane& state, const M1x7& state7) const {

  // state5: (q/p, u', v'. u, v)

  double spu(1.);

  const TVector3& O(state.getPlane()->getO());
  const TVector3& U(state.getPlane()->getU());
  const TVector3& V(state.getPlane()->getV());
  const TVector3& W(state.getPlane()->getNormal());

  // force A to be in normal direction and set spu accordingly
  double AtW( state7[3]*W.X() + state7[4]*W.Y() + state7[5]*W.Z() );
  if (AtW < 0.) {
    //fDir *= -1.;
    //AtW *= -1.;
    spu = -1.;
  }

  double* state5 = state.getState().GetMatrixArray();

  state5[0] = state7[6]; // q/p
  state5[1] = (state7[3]*U.X() + state7[4]*U.Y() + state7[5]*U.Z()) / AtW; // u' = (A * U) / (A * W)
  state5[2] = (state7[3]*V.X() + state7[4]*V.Y() + state7[5]*V.Z()) / AtW; // v' = (A * V) / (A * W)
  state5[3] = ((state7[0]-O.X())*U.X() +
               (state7[1]-O.Y())*U.Y() +
               (state7[2]-O.Z())*U.Z()); // u = (pos - O) * U
  state5[4] = ((state7[0]-O.X())*V.X() +
               (state7[1]-O.Y())*V.Y() +
               (state7[2]-O.Z())*V.Z()); // v = (pos - O) * V

  setSpu(state, spu);

}

void RKTrackRep::calcJ_pM_5x7(M5x7& J_pM, const TVector3& U, const TVector3& V, const M1x3& pTilde, double spu) const {
  /*if (debugLvl_ > 1) {
    debugOut << "RKTrackRep::calcJ_pM_5x7 \n";
    debugOut << "  U = "; U.Print();
    debugOut << "  V = "; V.Print();
    debugOut << "  pTilde = "; RKTools::printDim(pTilde, 3,1);
    debugOut << "  spu = " << spu << "\n";
  }*/

  std::fill(J_pM.begin(), J_pM.end(), 0);

  const double pTildeMag = sqrt(pTilde[0]*pTilde[0] + pTilde[1]*pTilde[1] + pTilde[2]*pTilde[2]);
  const double pTildeMag2 = pTildeMag*pTildeMag;

  const double utpTildeOverpTildeMag2 = (U.X()*pTilde[0] + U.Y()*pTilde[1] + U.Z()*pTilde[2]) / pTildeMag2;
  const double vtpTildeOverpTildeMag2 = (V.X()*pTilde[0] + V.Y()*pTilde[1] + V.Z()*pTilde[2]) / pTildeMag2;

  //J_pM matrix is d(x,y,z,ax,ay,az,q/p) / d(q/p,u',v',u,v)   (out is 7x7)

   // d(x,y,z)/d(u)
  J_pM[21] = U.X(); // [3][0]
  J_pM[22] = U.Y(); // [3][1]
  J_pM[23] = U.Z(); // [3][2]
  // d(x,y,z)/d(v)
  J_pM[28] = V.X(); // [4][0]
  J_pM[29] = V.Y(); // [4][1]
  J_pM[30] = V.Z(); // [4][2]
  // d(q/p)/d(q/p)
  J_pM[6] = 1.; // not needed for array matrix multiplication
  // d(ax,ay,az)/d(u')
  double fact = spu / pTildeMag;
  J_pM[10] = fact * ( U.X() - pTilde[0]*utpTildeOverpTildeMag2 ); // [1][3]
  J_pM[11] = fact * ( U.Y() - pTilde[1]*utpTildeOverpTildeMag2 ); // [1][4]
  J_pM[12] = fact * ( U.Z() - pTilde[2]*utpTildeOverpTildeMag2 ); // [1][5]
  // d(ax,ay,az)/d(v')
  J_pM[17] = fact * ( V.X() - pTilde[0]*vtpTildeOverpTildeMag2 ); // [2][3]
  J_pM[18] = fact * ( V.Y() - pTilde[1]*vtpTildeOverpTildeMag2 ); // [2][4]
  J_pM[19] = fact * ( V.Z() - pTilde[2]*vtpTildeOverpTildeMag2 ); // [2][5]

  /*if (debugLvl_ > 1) {
    debugOut << "  J_pM_5x7_ = "; RKTools::printDim(J_pM_5x7_, 5,7);
  }*/
}


void RKTrackRep::transformPM6(const MeasuredStateOnPlane& state,
                              M6x6& out6x6) const {

  // get vectors and aux variables
  const TVector3& U(state.getPlane()->getU());
  const TVector3& V(state.getPlane()->getV());
  const TVector3& W(state.getPlane()->getNormal());

  const TVectorD& state5(state.getState());
  double spu(getSpu(state));

  M1x3 pTilde;
  pTilde[0] = spu * (W.X() + state5(1)*U.X() + state5(2)*V.X()); // a_x
  pTilde[1] = spu * (W.Y() + state5(1)*U.Y() + state5(2)*V.Y()); // a_y
  pTilde[2] = spu * (W.Z() + state5(1)*U.Z() + state5(2)*V.Z()); // a_z

  const double pTildeMag = sqrt(pTilde[0]*pTilde[0] + pTilde[1]*pTilde[1] + pTilde[2]*pTilde[2]);
  const double pTildeMag2 = pTildeMag*pTildeMag;

  const double utpTildeOverpTildeMag2 = (U.X()*pTilde[0] + U.Y()*pTilde[1] + U.Z()*pTilde[2]) / pTildeMag2;
  const double vtpTildeOverpTildeMag2 = (V.X()*pTilde[0] + V.Y()*pTilde[1] + V.Z()*pTilde[2]) / pTildeMag2;

  //J_pM matrix is d(x,y,z,px,py,pz) / d(q/p,u',v',u,v)       (out is 6x6)

  const double qop = state5(0);
  const double p = getCharge(state)/qop; // momentum

  M5x6 J_pM_5x6;
  std::fill(J_pM_5x6.begin(), J_pM_5x6.end(), 0);

  // d(px,py,pz)/d(q/p)
  double fact = -1. * p / (pTildeMag * qop);
  J_pM_5x6[3] = fact * pTilde[0]; // [0][3]
  J_pM_5x6[4] = fact * pTilde[1]; // [0][4]
  J_pM_5x6[5] = fact * pTilde[2]; // [0][5]
  // d(px,py,pz)/d(u')
  fact = p * spu / pTildeMag;
  J_pM_5x6[9]  = fact * ( U.X() - pTilde[0]*utpTildeOverpTildeMag2 ); // [1][3]
  J_pM_5x6[10] = fact * ( U.Y() - pTilde[1]*utpTildeOverpTildeMag2 ); // [1][4]
  J_pM_5x6[11] = fact * ( U.Z() - pTilde[2]*utpTildeOverpTildeMag2 ); // [1][5]
  // d(px,py,pz)/d(v')
  J_pM_5x6[15] = fact * ( V.X() - pTilde[0]*vtpTildeOverpTildeMag2 ); // [2][3]
  J_pM_5x6[16] = fact * ( V.Y() - pTilde[1]*vtpTildeOverpTildeMag2 ); // [2][4]
  J_pM_5x6[17] = fact * ( V.Z() - pTilde[2]*vtpTildeOverpTildeMag2 ); // [2][5]
  // d(x,y,z)/d(u)
  J_pM_5x6[18] = U.X(); // [3][0]
  J_pM_5x6[19] = U.Y(); // [3][1]
  J_pM_5x6[20] = U.Z(); // [3][2]
  // d(x,y,z)/d(v)
  J_pM_5x6[24] = V.X(); // [4][0]
  J_pM_5x6[25] = V.Y(); // [4][1]
  J_pM_5x6[26] = V.Z(); // [4][2]


  // do the transformation
  // out = J_pM^T * in5x5 * J_pM
  const M5x5& in5x5_ = *((M5x5*) state.getCov().GetMatrixArray());
  RKTools::J_pMTxcov5xJ_pM(J_pM_5x6, in5x5_, out6x6);

}

void RKTrackRep::calcJ_Mp_7x5(M7x5& J_Mp, const TVector3& U, const TVector3& V, const TVector3& W, const M1x3& A) const {

  /*if (debugLvl_ > 1) {
    debugOut << "RKTrackRep::calcJ_Mp_7x5 \n";
    debugOut << "  U = "; U.Print();
    debugOut << "  V = "; V.Print();
    debugOut << "  W = "; W.Print();
    debugOut << "  A = "; RKTools::printDim(A, 3,1);
  }*/

  std::fill(J_Mp.begin(), J_Mp.end(), 0);

  const double AtU = A[0]*U.X() + A[1]*U.Y() + A[2]*U.Z();
  const double AtV = A[0]*V.X() + A[1]*V.Y() + A[2]*V.Z();
  const double AtW = A[0]*W.X() + A[1]*W.Y() + A[2]*W.Z();

  // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,ax,ay,az,q/p)   (in is 7x7)

  // d(u')/d(ax,ay,az)
  double fact = 1./(AtW*AtW);
  J_Mp[16] = fact * (U.X()*AtW - W.X()*AtU); // [3][1]
  J_Mp[21] = fact * (U.Y()*AtW - W.Y()*AtU); // [4][1]
  J_Mp[26] = fact * (U.Z()*AtW - W.Z()*AtU); // [5][1]
  // d(v')/d(ax,ay,az)
  J_Mp[17] = fact * (V.X()*AtW - W.X()*AtV); // [3][2]
  J_Mp[22] = fact * (V.Y()*AtW - W.Y()*AtV); // [4][2]
  J_Mp[27] = fact * (V.Z()*AtW - W.Z()*AtV); // [5][2]
  // d(q/p)/d(q/p)
  J_Mp[30] = 1.; // [6][0]  - not needed for array matrix multiplication
  //d(u)/d(x,y,z)
  J_Mp[3]  = U.X(); // [0][3]
  J_Mp[8]  = U.Y(); // [1][3]
  J_Mp[13] = U.Z(); // [2][3]
  //d(v)/d(x,y,z)
  J_Mp[4]  = V.X(); // [0][4]
  J_Mp[9]  = V.Y(); // [1][4]
  J_Mp[14] = V.Z(); // [2][4]

  /*if (debugLvl_ > 1) {
    debugOut << "  J_Mp_7x5_ = "; RKTools::printDim(J_Mp, 7,5);
  }*/

}


void RKTrackRep::transformM6P(const M6x6& in6x6,
                              const M1x7& state7,
                              MeasuredStateOnPlane& state) const { // plane and charge must already be set!

  // get vectors and aux variables
  const TVector3& U(state.getPlane()->getU());
  const TVector3& V(state.getPlane()->getV());
  const TVector3& W(state.getPlane()->getNormal());

  const double AtU = state7[3]*U.X() + state7[4]*U.Y() + state7[5]*U.Z();
  const double AtV = state7[3]*V.X() + state7[4]*V.Y() + state7[5]*V.Z();
  const double AtW = state7[3]*W.X() + state7[4]*W.Y() + state7[5]*W.Z();

  // J_Mp matrix is d(q/p,u',v',u,v) / d(x,y,z,px,py,pz)       (in is 6x6)

  const double qop = state7[6];
  const double p = getCharge(state)/qop; // momentum

  M6x5 J_Mp_6x5;
  std::fill(J_Mp_6x5.begin(), J_Mp_6x5.end(), 0);

  //d(u)/d(x,y,z)
  J_Mp_6x5[3]  = U.X(); // [0][3]
  J_Mp_6x5[8]  = U.Y(); // [1][3]
  J_Mp_6x5[13] = U.Z(); // [2][3]
  //d(v)/d(x,y,z)
  J_Mp_6x5[4]  = V.X(); // [0][4]
  J_Mp_6x5[9]  = V.Y(); // [1][4]
  J_Mp_6x5[14] = V.Z(); // [2][4]
  // d(q/p)/d(px,py,pz)
  double fact = (-1.) * qop / p;
  J_Mp_6x5[15] = fact * state7[3]; // [3][0]
  J_Mp_6x5[20] = fact * state7[4]; // [4][0]
  J_Mp_6x5[25] = fact * state7[5]; // [5][0]
  // d(u')/d(px,py,pz)
  fact = 1./(p*AtW*AtW);
  J_Mp_6x5[16] = fact * (U.X()*AtW - W.X()*AtU); // [3][1]
  J_Mp_6x5[21] = fact * (U.Y()*AtW - W.Y()*AtU); // [4][1]
  J_Mp_6x5[26] = fact * (U.Z()*AtW - W.Z()*AtU); // [5][1]
  // d(v')/d(px,py,pz)
  J_Mp_6x5[17] = fact * (V.X()*AtW - W.X()*AtV); // [3][2]
  J_Mp_6x5[22] = fact * (V.Y()*AtW - W.Y()*AtV); // [4][2]
  J_Mp_6x5[27] = fact * (V.Z()*AtW - W.Z()*AtV); // [5][2]

  // do the transformation
  // out5x5 = J_Mp^T * in * J_Mp
  M5x5& out5x5_ = *((M5x5*) state.getCov().GetMatrixArray());
  RKTools::J_MpTxcov6xJ_Mp(J_Mp_6x5, in6x6, out5x5_);

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
bool RKTrackRep::RKutta(const M1x4& SU,
                        const DetPlane& plane,
                        double charge,
                        double mass,
                        M1x7& state7,
                        M7x7* jacobianT,
                        M1x7* J_MMT_unprojected_lastRow,
                        double& coveredDistance,
                        double& flightTime,
                        bool& checkJacProj,
                        M7x7& noiseProjection,
                        StepLimits& limits,
                        bool onlyOneStep,
                        bool calcOnlyLastRowOfJ) const {

  // limits, check-values, etc. Can be tuned!
  static const double Wmax           ( 3000. );           // max. way allowed [cm]
  static const double AngleMax       ( 6.3 );           // max. total angle change of momentum. Prevents extrapolating a curler round and round if no active plane is found.
  static const double Pmin           ( 4.E-3 );           // minimum momentum for propagation [GeV]
  static const unsigned int maxNumIt ( 1000 );    // maximum number of iterations in main loop
  // Aux parameters
  M1x3&   R          ( *((M1x3*) &state7[0]) );  // Start coordinates  [cm]  (x,  y,  z)
  M1x3&   A          ( *((M1x3*) &state7[3]) );  // Start directions         (ax, ay, az);   ax^2+ay^2+az^2=1
  M1x3    SA         = {{0.,0.,0.}};             // Start directions derivatives dA/S
  double  Way        ( 0. );                     // Sum of absolute values of all extrapolation steps [cm]
  double  momentum   ( fabs(charge/state7[6]) ); // momentum [GeV]
  double  relMomLoss ( 0 );                      // relative momentum loss in RKutta
  double  deltaAngle ( 0. );                     // total angle by which the momentum has changed during extrapolation
  double  An(0), S(0), Sl(0), CBA(0);


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

    M1x3 ABefore = {{ A[0], A[1], A[2] }};
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
    double arg = ABefore[0]*A[0] + ABefore[1]*A[1] + ABefore[2]*A[2];
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
      SA[0]*=Sl;  SA[1]*=Sl;  SA[2]*=Sl; // SA/Sl = delta A / delta way; local derivative of A with respect to the length of the way

      // calculate A
      A[0] += SA[0]*S;    // S  = distance to surface
      A[1] += SA[1]*S;    // A = A + S * SA*Sl
      A[2] += SA[2]*S;

      // normalize A
      CBA = 1./sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);  // 1/|A|
      A[0] *= CBA; A[1] *= CBA; A[2] *= CBA;

      R[0] += S*(A[0]-0.5*S*SA[0]);    // R = R + S*(A - 0.5*S*SA); approximation for final point on surface
      R[1] += S*(A[1]-0.5*S*SA[1]);
      R[2] += S*(A[2]-0.5*S*SA[2]);


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
      double norm;
      int i=0;
      if (calcOnlyLastRowOfJ)
        i = 42;

      double* jacPtr = jacobianT->begin();

      for(unsigned int j=42; j<49; j+=7) {
        (*J_MMT_unprojected_lastRow)[j-42] = jacPtr[j];
      }

      for(; i<49; i+=7) {
        norm = (jacPtr[i]*SU[0] + jacPtr[i+1]*SU[1] + jacPtr[i+2]*SU[2]) * An;  // dR_normal / A_normal
        jacPtr[i]   -= norm*A [0];   jacPtr[i+1] -= norm*A [1];   jacPtr[i+2] -= norm*A [2];
        jacPtr[i+3] -= norm*SA[0];   jacPtr[i+4] -= norm*SA[1];   jacPtr[i+5] -= norm*SA[2];
      }
      checkJacProj = true;


      if (debugLvl_ > 0) {
        //debugOut << "  Jacobian^T of extrapolation after Projection:\n";
        //RKTools::printDim(*jacobianT, 7,7);
      }

      if (!calcOnlyLastRowOfJ) {
        for (int iRow = 0; iRow < 3; ++iRow) {
          for (int iCol = 0; iCol < 3; ++iCol) {
            noiseProjection[iRow*7 + iCol]       = (iRow == iCol) - An * SU[iCol] * A[iRow];
            noiseProjection[(iRow + 3)*7 + iCol] =                - An * SU[iCol] * SA[iRow];
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


double RKTrackRep::estimateStep(const M1x7& state7,
                                const M1x4& SU,
                                const DetPlane& plane,
                                const double& charge,
                                double& relMomLoss,
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

  if (debugLvl_ > 0) {
    debugOut << " RKTrackRep::estimateStep \n";
    debugOut << "  position:  "; TVector3(state7[0], state7[1], state7[2]).Print();
    debugOut << "  direction: "; TVector3(state7[3], state7[4], state7[5]).Print();
  }

  // calculate SL distance to surface
  double Dist ( SU[3] - (state7[0]*SU[0] +
                         state7[1]*SU[1] +
                         state7[2]*SU[2])  );  // Distance between start coordinates and surface
  double An ( state7[3]*SU[0] +
              state7[4]*SU[1] +
              state7[5]*SU[2]   );              // An = dir * N;  component of dir normal to surface

  double SLDist; // signed
  if (fabs(An) > 1.E-10)
    SLDist = Dist/An;
  else {
    SLDist = Dist*1.E10;
    if (An<0) SLDist *= -1.;
  }

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
  double fieldCurvLimit( limits.getLowestLimitSignedVal() ); // signed
  std::pair<double, double> distVsStep (9.E99, 9.E99); // first: smallest straight line distances to plane; second: RK steps

  static const unsigned int maxNumIt = 10;
  unsigned int counter(0);

  while (fabs(fieldCurvLimit) > MINSTEP) {

    if(++counter > maxNumIt){
      // if max iterations are reached, take a safe value
      // (in previous iteration, fieldCurvLimit has been not more than doubled)
      // and break.
      fieldCurvLimit *= 0.5;
      break;
    }

    M1x7 state7_temp = state7;
    M1x3 SA = {{0, 0, 0}};

    double q ( RKPropagate(state7_temp, nullptr, SA, fieldCurvLimit, true) );
    if (debugLvl_ > 0) {
      debugOut << "  maxStepArg = " << fieldCurvLimit << "; q = " << q  << " \n";
    }

    // remember steps and resulting SL distances to plane for stepsize improvement
    // calculate distance to surface
    Dist = SU[3] - (state7_temp[0] * SU[0] +
                    state7_temp[1] * SU[1] +
                    state7_temp[2] * SU[2]); // Distance between position and surface

    An = state7_temp[3] * SU[0] +
         state7_temp[4] * SU[1] +
         state7_temp[5] * SU[2];    // An = dir * N;  component of dir normal to surface

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
    M1x7 state7_temp = {{ state7[0], state7[1], state7[2], state7[3], state7[4], state7[5], state7[6] }};

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
                          M1x7& state7,
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
  double dqop(0.);

  const TVector3 W(destPlane.getNormal());
  M1x4 SU = {{W.X(), W.Y(), W.Z(), destPlane.distance(0., 0., 0.)}};

  // make SU vector point away from origin
  if (W*destPlane.getO() < 0) {
    SU[0] *= -1;
    SU[1] *= -1;
    SU[2] *= -1;
  }


  M1x7 startState7 = state7;

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
    for(int i = 0; i < 7*7; ++i) J_MMT_[i] = 0;
    for(int i=0; i<7; ++i) J_MMT_[8*i] = 1.;

    M7x7* noise = nullptr;
    isAtBoundary = false;

    // propagation
    bool checkJacProj = false;
    limits_.reset();
    limits_.setLimit(stp_sMaxArg, maxStep-fabs(coveredDistance));

    M1x7 J_MMT_unprojected_lastRow = {{0, 0, 0, 0, 0, 0, 1}};

    if( ! RKutta(SU, destPlane, charge, mass, state7, &J_MMT_, &J_MMT_unprojected_lastRow,
		 coveredDistance, flightTime, checkJacProj, noiseProjection_,
		 limits_, onlyOneStep, !fillExtrapSteps) ) {
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
    if(fillExtrapSteps) {
      noise = &noiseArray_;
      for(int i = 0; i < 7*7; ++i) noiseArray_[i] = 0; // set noiseArray_ to 0
    }

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
          debugOut << "7D noise: \n";
          RKTools::printDim(noise->begin(), 7, 7);
        }
      }

      // do momLoss only for defined 1/momentum .ne.0
      if(fabs(state7[6])>1.E-10) {

        if (debugLvl_ > 0) {
          debugOut << "correct state7 with dx/dqop, dy/dqop ...\n";
        }

        dqop = charge/(fabs(charge/state7[6])-momLoss) - state7[6];

        // Correct coveredDistance and flightTime and momLoss if checkJacProj == true
        // The idea is to calculate the state correction (based on the mometum loss) twice:
        // Once with the unprojected Jacobian (which preserves coveredDistance),
        // and once with the projected Jacobian (which is constrained to the plane and does NOT preserve coveredDistance).
        // The difference of these two corrections can then be used to calculate a correction factor.
        if (checkJacProj && fabs(coveredDistance) > MINSTEP) {
          M1x3 state7_correction_unprojected = {{0, 0, 0}};
          for (unsigned int i=0; i<3; ++i) {
            state7_correction_unprojected[i] = 0.5 * dqop * J_MMT_unprojected_lastRow[i];
            //debugOut << "J_MMT_unprojected_lastRow[i] " << J_MMT_unprojected_lastRow[i] << "\n";
            //debugOut << "state7_correction_unprojected[i] " << state7_correction_unprojected[i] << "\n";
          }

          M1x3 state7_correction_projected = {{0, 0, 0}};
          for (unsigned int i=0; i<3; ++i) {
            state7_correction_projected[i] = 0.5 * dqop * J_MMT_[6*7 + i];
            //debugOut << "J_MMT_[6*7 + i] " << J_MMT_[6*7 + i] << "\n";
            //debugOut << "state7_correction_projected[i] " << state7_correction_projected[i] << "\n";
          }

          // delta distance
          M1x3 delta_state = {{0, 0, 0}};
          for (unsigned int i=0; i<3; ++i) {
            delta_state[i] = state7_correction_unprojected[i] - state7_correction_projected[i];
          }

          double Dist( sqrt(delta_state[0]*delta_state[0]
              + delta_state[1]*delta_state[1]
              + delta_state[2]*delta_state[2] )  );

          // sign: delta * a
          if (delta_state[0]*state7[3] + delta_state[1]*state7[4] + delta_state[2]*state7[5] > 0)
            Dist *= -1.;

          double correctionFactor( 1. + Dist / coveredDistance );
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
        double qop( charge/(fabs(charge/state7[6])-momLoss) );
        dqop = qop - state7[6];
        state7[6] = qop;

        for (unsigned int i=0; i<6; ++i) {
          state7[i] += 0.5 * dqop * J_MMT_[6*7 + i];
        }
        // normalize direction, just to make sure
        double norm( 1. / sqrt(state7[3]*state7[3] + state7[4]*state7[4] + state7[5]*state7[5]) );
        for (unsigned int i=3; i<6; ++i)
          state7[i] *= norm;

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
        RKTools::Np_N_NpT(noiseProjection_, noiseArray_);

        if (debugLvl_ > 1) {
          debugOut << "7D noise projected onto plane: \n";
          RKTools::printDim(noiseArray_.begin(), 7, 7);
        }
      }

      // Store this step's noise for final calculation.
      lastStep->noise7_ = noiseArray_;

      if (debugLvl_ > 2) {
        debugOut<<"ExtrapSteps \n";
        for (std::vector<ExtrapStep>::iterator it = ExtrapSteps_.begin(); it != ExtrapSteps_.end(); ++it){
          debugOut << "7D Jacobian: "; RKTools::printDim((it->jac7_.begin()), 5,5);
          debugOut << "7D noise:    "; RKTools::printDim((it->noise7_.begin()), 5,5);
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
      cov->Similarity(fJacobian_);
      *cov += fNoise_;
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


double RKTrackRep::momMag(const M1x7& state7) const {
  // FIXME given this interface this function cannot work for charge =/= +-1
  return fabs(1/state7[6]);
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
