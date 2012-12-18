/* Copyright 2008-2010, Technische Universitaet Muenchen,
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
/** @addtogroup genfit
 * @{
 */

#ifndef GFABSWIREHIT_H
#define GFABSWIREHIT_H


#include "GFAbsRecoHit.h"

/** @brief Abstract hit class for hits in wire detectors (Straw tubes and drift chambers)
 *  which do not measure the coordinate along the wire.
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Lia Lavezzi (INFN Pavia, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * This hit class is not valid for any kind of plane orientation
 * choice: to use it you MUST choose a plane described by u 
 * and v axes with v coincident with the wire (and u orthogonal
 * to it, obviously).
 * The hit will be described by 7 coordinates:
 * w_x1, w_y1, w_z1, w_x2, w_y2, w_z2, rdrift
 * where w_ji (with j = x, y, z and i = 1, 2) are the wire
 * extremities coordinates; rdrift = distance from the wire (u 
 * coordinate in the plane)
 *
 */
class GFAbsWireHit : public GFAbsRecoHit {
 public:

  // Constructors/Destructors ---------
  /** @brief Constructor for inheriting hits with a higher dimensionality (e.g. position along wire, energy loss) */
  GFAbsWireHit(unsigned int dim = 7);

  virtual ~GFAbsWireHit(){;}

  // Operations ----------------------
  virtual void getMeasurement(const GFAbsTrackRep* rep,
                              const GFDetPlane& pl,
                              const TVectorD& statePred,
                              const TMatrixDSym& covPred,
                              TVectorD& m,
                              TMatrixDSym& V);

   /** @brief Get detector plane 
    * Calls GFAbsTrackRep::extrapolateToLine for POCA.
    * The detector plane will contain the wire as plane vector v.
    * The origin of the plane lies on the wire (first 3 hit coordinates).
    * The plane vector u will be perpendicular to the track direction at the POCA.
    * u = +-1 * (track direction) x (wire direction)
    * The direction of u will be selected according to fLeftRight.
    */
  virtual const GFDetPlane& getDetPlane(GFAbsTrackRep* rep);

  double getMaxDistance(){return fMaxdistance;}
  void setMaxDistance(double d){fMaxdistance = d;}
  
  /**
   * select how to resolve the left/right ambiguity:
   * -1: negative (left) side on vector (track direction) x (wire direction)
   * 0: auto select (take side with smallest distance to track)
   * 1: positive (right) side on vector (track direction) x (wire direction)
   */
  void setLeftRightResolution(int lr);
  int getLeftRightResolution() const {return fLeftRight;}

 protected:

  // Protected Data Members ------------
  double fMaxdistance;
  int fLeftRight;

  // Protected Methods -----------------
  /** @brief Check if the detector plane is valid
   */
  void checkPlane(const GFDetPlane&);

 public:
  ClassDef(GFAbsWireHit,1);

};

#endif

/** @} */
