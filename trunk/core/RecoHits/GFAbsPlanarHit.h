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

#ifndef GFABSPLANARHIT_H
#define GFABSPLANARHIT_H


#include "GFAbsRecoHit.h"

/** @brief Policy class implementing a planar hit geometry. 
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Johannes Rauch  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 * RecoHits for planar detectors should inherit 
 * from GFAbsPlanarHit.
 *
 * The main feature of this type of hit is, that the detector plane
 * is completely defined by the detector hardware. Derived RecoHits need only
 * to supply the physical detector plane from their geometry database.
 */

class GFAbsPlanarHit : public GFAbsRecoHit {
public:

  // Constructors/Destructors ---------
  /** @brief Dimensionality is usually 1 (strip hit, default) or 2 (pixel hit)
   *  If the dimesionality is higher (e.g. energy loss), getMeasurement() has to be
   *  adapted in the inheriting hit.
   *  */
  GFAbsPlanarHit(unsigned int dim = 1) : GFAbsRecoHit(dim) {assert(dim >=1);}
  
  virtual ~GFAbsPlanarHit() {}

  // Accessors -----------------------
  
  /** @brief Returns the physical detector plane.
   */
  const GFDetPlane& getDetPlane(GFAbsTrackRep*) {return fDetPlane;}
  
  virtual void getMeasurement(const GFAbsTrackRep* rep,
                              const GFDetPlane& pl,
                              const TVectorD& statePred,
                              const TMatrixDSym& covPred,
                              TVectorD& m,
                              TMatrixDSym& V);

  // Modifiers -----------------------

  /** @brief Set physical detector plane. 
   * Needs to be called before hit can be used.
   *
   * For the planar detector the detector plane is fixed by the geometry of the
   * detector hardware. This method should be called in the constructor of
   * any derived RecoHit in order to setup the geometry of this hit.
   */
  void setDetPlane(const GFDetPlane& p){fDetPlane = p;}

 public:
  ClassDef(GFAbsPlanarHit,1)

};

#endif

/** @} */
