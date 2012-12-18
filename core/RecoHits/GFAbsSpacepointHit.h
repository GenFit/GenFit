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

#ifndef GFABSSPACEPOINTHIT_H
#define GFABSSPACEPOINTHIT_H


#include "GFAbsRecoHit.h"

/** @brief Abstract hit class implementing a space point hit geometry.
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 * RecoHits for detectors measuring 3D space points should inherit 
 * from GFAbsSpacepointHit.
 *
 * For a space point the detector plane has to be defined with respect to
 * a track representation. GFAbsSpacepointHit implements a scheme where the
 * detectorplane is chosen perpendicular to the track.
 * In a track fit, only two of the three coordinates of a space point are
 * independent (the track is a one-dimensional object). Therefore the 3D
 * data of the hit is used to define a proper detector plane into which the
 * hit coordinates are then projected.
 */

class GFAbsSpacepointHit : public GFAbsRecoHit {
 public:

  // Constructors/Destructors ---------
  /** @brief Constructor for inheriting hits with a higher dimensionality (e.g. energy loss) */
  GFAbsSpacepointHit(unsigned int dim = 3) : GFAbsRecoHit(dim) {assert(dim >= 3);}

  virtual ~GFAbsSpacepointHit(){;}

  
  // Operations ----------------------
  virtual void getMeasurement(const GFAbsTrackRep* rep,
                              const GFDetPlane& pl,
                              const TVectorD& statePred,
                              const TMatrixDSym& covPred,
                              TVectorD& m,
                              TMatrixDSym& V);

  /** @brief Get detector plane perpendicular to track.
   *
   * The detector plane is contructed from the position of the hit and
   * the track representation. For this the track is extrapolated to the
   * point of closest approach to the hit.
   */
  virtual const GFDetPlane& getDetPlane(GFAbsTrackRep* rep);

 private:

 public:
  ClassDef(GFAbsSpacepointHit,1);
};

#endif

/** @} */
