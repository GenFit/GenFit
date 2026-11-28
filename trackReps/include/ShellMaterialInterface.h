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

/** @addtogroup RKTrackRep
 * @{
 */

#pragma once

#include "AbsMaterialInterface.h"

#include <vector>


namespace genfit {

/**
 * @brief Thread-safe analytic AbsMaterialInterface for coaxial cylindrical shells.
 *
 * Thread-safe alternative to TGeoMaterialInterface: ROOT's TGeo navigation is not
 * safe for concurrent use, so material lookups through it must be serialized,
 * while this implementation computes all boundaries in closed form and keeps its
 * per-call state in a thread_local block, so several threads may share one
 * instance and nothing is serialized.
 *
 * Geometry: coaxial cylindrical shells (pairs of inner/outer radii, common half
 * length) embedded in a bulk material. The materials are constant per region:
 * shellMaterial inside a shell, bulkMaterial everywhere else. Intended for
 * detectors well described by coaxial shells.
 *
 * The step search mirrors TGeoMaterialInterface::findNextBoundary() (same
 * delta/epsilon/RKPropagate arch control); only the geometry queries are
 * replaced by closed-form distances to the cylinder and endcap surfaces.
 */
class ShellMaterialInterface : public AbsMaterialInterface {

 public:

  /**
   * @param bulkMaterial  material outside the shells
   * @param shellMaterial material inside the shells
   * @param shellRadii    [rmin1, rmax1, rmin2, rmax2, ...]: inner/outer radii of the shells
   * @param zHalfLength   common half length of the shells; endcaps at z = +-zHalfLength
   */
  ShellMaterialInterface(const Material& bulkMaterial,
                         const Material& shellMaterial,
                         const std::vector<double>& shellRadii,
                         double zHalfLength)
      : bulkMaterial_(bulkMaterial), shellMaterial_(shellMaterial),
        shellRadii_(shellRadii), zHalfLength_(zHalfLength) {}

  /** @brief Initialize the navigator at given position and with given
      direction.  Returns true if the region (shell / bulk) changed.
   */
  bool initTrack(double posX, double posY, double posZ,
                 double dirX, double dirY, double dirZ) override;

  Material getMaterialParameters() override;

  /** @brief Make a step (following the curvature) until step length
   * sMax or the next boundary is reached.  After making a step to a
   * boundary, the position has to be beyond the boundary, i.e. the
   * current material has to be that beyond the boundary.  The actual
   * step made is returned.
   */
  double findNextBoundary(const RKTrackRep* rep,
                          const M1x7& state7,
                          double sMax,
                          bool varField = true) override;

 private:

  struct State {
    double pos[3] = {0., 0., 0.};
    double dir[3] = {0., 0., 1.};
    bool inShell = false;
    bool valid = false;
  };

  // Per-thread navigator state (initTrack/findNextBoundary/getMaterialParameters
  // form one usage sequence per thread, as documented in AbsMaterialInterface).
  static State& state() {
    static thread_local State s;
    return s;
  }

  //! Position inside one of the cylindrical shells (within +-zHalfLength)?
  bool inShell(double x, double y, double z) const;

  //! Lower bound of the distance to any boundary (any surface) and the distance
  //! to the nearest boundary surface along the current direction, capped at maxDist.
  void safetyStep(double maxDist, double& safety, double& step) const;

  Material bulkMaterial_;
  Material shellMaterial_;
  std::vector<double> shellRadii_;
  double zHalfLength_;

};

} /* End of namespace genfit */
/** @} */
