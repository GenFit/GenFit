/* Copyright 2008-2026, Technische Universitaet Muenchen, DESY
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch & Christian Wessel

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

#include "ShellMaterialInterface.h"
#include "Exception.h"
#include "MathHelpers.h"

#include <algorithm>
#include <math.h>


namespace genfit {


bool ShellMaterialInterface::initTrack(double posX, double posY, double posZ,
                                       double dirX, double dirY, double dirZ) {
  State& s = state();
  const bool shell = inShell(posX, posY, posZ);
  const bool changed = s.valid && (shell != s.inShell);
  s.pos[0] = posX; s.pos[1] = posY; s.pos[2] = posZ;
  s.dir[0] = dirX; s.dir[1] = dirY; s.dir[2] = dirZ;
  s.inShell = shell;
  s.valid = true;
  return changed;
}


Material ShellMaterialInterface::getMaterialParameters() {
  return state().inShell ? shellMaterial_ : bulkMaterial_;
}


double
ShellMaterialInterface::findNextBoundary(const RKTrackRep* rep,
                                         const M1x7& stateOrig,
                                         double sMax, // signed
                                         bool varField)
{
  const double delta(1.E-2); // cm, distance limit beneath which straight-line steps are taken.
  const double epsilon(1.E-1); // cm, allowed upper bound on arch deviation

  M1x3 SA;
  M1x7 state7, oldState7;
  oldState7 = stateOrig;

  int stepSign(sMax < 0 ? -1 : 1);

  double s = 0;  // trajectory length to boundary

  const unsigned maxIt = 300;
  unsigned it = 0;

  // Initialize to the current location (set by the caller via initTrack).
  double safety, slDist;
  safetyStep(fabs(sMax) - s, safety, slDist);
  if (slDist < safety)
    slDist = safety;
  double step = slDist;

  while (1) {
    if (++it > maxIt) {
      Exception exc("ShellMaterialInterface::findNextBoundary ==> maximum number of iterations exceeded", __LINE__, __FILE__);
      exc.setFatal();
      throw exc;
    }

    // No boundary in sight?
    if (s + safety > fabs(sMax)) {
      return stepSign * (s + safety);
    }

    // We are at the boundary?
    if (slDist < delta) {
      return stepSign * (s + slDist);
    }

    // Follow curved arch, then see if we may have missed a boundary.
    // Always propagate complete way from original start to avoid
    // inconsistent extrapolations.
    state7 = stateOrig;
    rep->RKPropagate(state7, nullptr, SA, stepSign * (s + step), varField);

    // Straight line distance^2 between extrapolation finish and
    // the end of the previously determined safe segment.
    double dist2 = (square(state7[0] - oldState7[0])
        + square(state7[1] - oldState7[1])
        + square(state7[2] - oldState7[2]));
    // Maximal lateral deviation^2.
    double maxDeviation2 = 0.25 * (step * step - dist2);

    if (step > safety && maxDeviation2 > epsilon * epsilon) {
      // Take a shorter step to reliably estimate material, never shorter
      // than safety.
      step = std::max(step / 2, safety);
    } else {
      const State saved = state();
      bool volChanged = initTrack(state7[0], state7[1], state7[2],
          stepSign * state7[3], stepSign * state7[4], stepSign * state7[5]);

      if (volChanged) {
        state() = saved; // move back to start
        if (step <= safety) {
          return stepSign * (s + step);
        }
        step = std::max(step / 2, safety);
      } else {
        // we are in the new place, the step was safe, advance
        s += step;
        oldState7 = state7;

        safetyStep(fabs(sMax) - s, safety, slDist);
        if (slDist < safety)
          slDist = safety;
        step = slDist;
      }
    }
  }
}


bool ShellMaterialInterface::inShell(double x, double y, double z) const {
  if (fabs(z) > zHalfLength_) return false;
  const double r2 = x * x + y * y;
  for (unsigned i = 0; i + 1 < shellRadii_.size(); i += 2)
    if (r2 >= shellRadii_[i] * shellRadii_[i] && r2 <= shellRadii_[i + 1] * shellRadii_[i + 1]) return true;
  return false;
}


void ShellMaterialInterface::safetyStep(double maxDist, double& safety, double& step) const {
  const State& s = state();
  const double x = s.pos[0], y = s.pos[1], z = s.pos[2];
  const double dx = s.dir[0], dy = s.dir[1], dz = s.dir[2];
  const double r = sqrt(x * x + y * y);

  // safety: distance to the nearest surface (never overestimates the true
  // distance to the boundary of the current region)
  safety = 1.E30;
  for (double R : shellRadii_) safety = std::min(safety, fabs(r - R));
  safety = std::min(safety, fabs(zHalfLength_ - fabs(z)));

  // step: nearest positive surface crossing along the direction
  double nearest = 1.E30;
  const double a = dx * dx + dy * dy;
  const double b = 2. * (x * dx + y * dy);
  for (double R : shellRadii_) {
    if (a < 1.E-30) continue;
    const double c = x * x + y * y - R * R;
    const double disc = b * b - 4. * a * c;
    if (disc <= 0.) continue;
    const double sq = sqrt(disc);
    for (double t : {(-b - sq) / (2. * a), (-b + sq) / (2. * a)}) {
      if (t <= 1.E-10) continue;
      if (fabs(z + t * dz) <= zHalfLength_) nearest = std::min(nearest, t);
    }
  }
  if (fabs(dz) > 1.E-30) {
    for (double zPlane : {-zHalfLength_, zHalfLength_}) {
      const double t = (zPlane - z) / dz;
      if (t <= 1.E-10) continue;
      const double rCross = sqrt(square(x + t * dx) + square(y + t * dy));
      for (unsigned i = 0; i + 1 < shellRadii_.size(); i += 2) {
        if (rCross >= shellRadii_[i] && rCross <= shellRadii_[i + 1]) {
          nearest = std::min(nearest, t);
          break;
        }
      }
    }
  }
  step = std::min(nearest, maxDist);
}


} /* End of namespace genfit */
