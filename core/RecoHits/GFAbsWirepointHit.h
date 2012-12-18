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
 * @{ */



#ifndef GFABSWIREPOINTHIT_H
#define GFABSWIREPOINTHIT_H

#include "GFAbsWireHit.h"


/** @brief Abstract hit class for hits in wire detectors (Straw tubes and drift chambers)
 *  which can measure the coordinate along the wire
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Lia Lavezzi (INFN Pavia, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, original author)
 *
 * This  is not valid for any kind of plane orientation
 * choice: to use it you MUST choose a plane described by u 
 * and v axes with v coincident with the wire (and u orthogonal
 * to it, obviously).
 * The hit will be described by 8 coordinates:
 * w_x1, w_y1, w_z1, w_x2, w_y2, w_z2, rdrift, zreco
 * where w_ji (with j = x, y, z and i = 1, 2) are the wire
 * extremities coordinates; rdrift = distance from the wire (u 
 * coordinate in the plane) and zreco = coordinate along the
 * wire (w2 -w1) (in the plane reference frame, v coordinate).
 *
 */
class GFAbsWirepointHit : public GFAbsWireHit {
 public:

  // Constructors/Destructors ---------
  /** @brief Constructor for inheriting hits with a higher dimensionality (e.g. position along wire, energy loss) */
  GFAbsWirepointHit(unsigned int dim = 8) : GFAbsWireHit(dim) {assert(dim >= 8);}

  virtual ~GFAbsWirepointHit(){;}

  // Operations ----------------------
  virtual void getMeasurement(const GFAbsTrackRep* rep,
                              const GFDetPlane& pl,
                              const TVectorD& statePred,
                              const TMatrixDSym& covPred,
                              TVectorD& m,
                              TMatrixDSym& V);


 public:
  ClassDef(GFAbsWirepointHit,1);

};

#endif

/** @} */
