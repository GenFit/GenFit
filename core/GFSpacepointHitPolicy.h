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

#ifndef GFSPACEPOINTHITPOLICY_H
#define GFSPACEPOINTHITPOLICY_H


#include "TMatrixT.h"
#include "TObject.h"

#include "GFDetPlane.h"

class GFAbsRecoHit;
class GFAbsTrackRep;

/** @brief Policy class implementing a space point hit geometry. 
 *
 *  @author Christian H&ouml;ppner (Technische Universit&auml;t M&uuml;nchen, original author)
 *  @author Sebastian Neubert  (Technische Universit&auml;t M&uuml;nchen, original author)
 * 
 * RecoHits for detectors measuring 3D space points should inherit 
 * from RecoHitIfc<GFSpacepointHitPolicy>.
 *
 * For a space point the detector plane has to be defined with respect to
 * a track representation. GFSpacepointHitPolicy implements a scheme where the
 * detectorplane is chosen perpendicular to the track.
 * In a track fit only 2 of the three coordinates of a space point are 
 * independent (the track is a one-dimensional object). Therefore the 3D
 * data of the hit is used to define a proper detector plane into which the
 * hit coordinates are then projected.
 */

class GFSpacepointHitPolicy : public TObject {
public:


  GFSpacepointHitPolicy(){;}
  
  // Operations ----------------------
   /** @brief Get detector plane perpendicular to track.
    *
    * The detector plane is contructed from the position of the hit and
    * the track representation. For this the track is extrapolated to the
    * point of closest approach to the hit.
    */
  const GFDetPlane& detPlane(GFAbsRecoHit*, GFAbsTrackRep*);

  /** @brief Hit coordinates in detector plane.
   */
  TMatrixT<double> hitCoord(GFAbsRecoHit*,const GFDetPlane&);

  /** @brief Hit covariances in detector plane.
   */
  TMatrixT<double> hitCov(GFAbsRecoHit*,const GFDetPlane&);

  virtual ~GFSpacepointHitPolicy(){;}

  const std::string& getName(){return fPolicyName;}
 private:
  static const std::string fPolicyName;

  // Private Data Members ------------
  GFDetPlane fPlane;

  // Private Methods -----------------

 public:
  ClassDef(GFSpacepointHitPolicy,1);
};

#endif

/** @} */
